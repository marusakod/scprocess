# snakemake rule for setting up scprocess_data directory:
# downloading files from scprocess_data repo
# building an index for alevin
# converting genes.gtf file to txt.gz file

rule make index:
  output:
    index_json = SCPROCESS_DATA_DIR + '/alevin_fry_home/' + INDEX_DIR_NAME + '/index/simpleaf_index.json
  threads: 16
  retries: 5
  resources:
    mem_mb      =  MB_RUN_SETUP
  conda:
    alevin_fry.yml
  shell:
    """



    """
