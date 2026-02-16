#!/bin/bash

# human 2020 genome (this code is a minor modification of code in https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps)


# Set up source directory
source="reference_sources"
mkdir -p "$source"

# URLs for downloading the source files
fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"

# Download FASTA if it doesn't exist
if [ ! -f "$fasta_in" ]; then
    wget -qO- "$fasta_url" | zcat > "$fasta_in"
fi

# Download GTF if it doesn't exist
if [ ! -f "$gtf_in" ]; then
    wget -qO- "$gtf_url" | zcat > "$gtf_in"
fi

# Modify sequence headers in the Ensembl FASTA to match the file "GRCh38.primary_assembly.genome.fa" from GENCODE
# Input FASTA: >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
# Output FASTA: >chr1 1
fasta_modified="genome.fa"  # Set output file to the current directory
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffix from transcript, gene, and exon IDs in order to match previous Cell Ranger reference packages
# Input GTF: ... gene_id "ENSG00000223972.5"; ...
# Output GTF: ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="genes.gtf.tmp"  # Temporary file for GTF modification
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

# Define string patterns for GTF tags
BIOTYPE_PATTERN="(protein_coding|lncRNA|IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene|rRNA|Mt_rRNA)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""

# Construct the gene ID allowlist
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "gene_allowlist"

# Filter the GTF file based on the gene allowlist
gtf_filtered="genes.gtf"  # Final output file
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"

# Clean up intermediate files (but don't delete the final files)
rm -f "$gtf_modified" "gene_allowlist"

# Remove the entire reference_sources directory
rm -rf "$source"
