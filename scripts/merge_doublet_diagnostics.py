#!/usr/bin/env python3

"""Combine readable per-run scDblFinder CSV.gz diagnostics."""

import argparse
import csv
import gzip
from pathlib import Path


def merge_csv_gz(input_files, output_file):
  """Merge gzip CSVs using the ordered union of sample-specific columns."""
  input_paths = [Path(input_file) for input_file in input_files]
  merged_header = []
  seen_columns = set()
  n_rows = 0
  output_path = Path(output_file)
  output_path.parent.mkdir(parents=True, exist_ok=True)

  for input_path in input_paths:
    if not input_path.exists():
      raise FileNotFoundError(input_path)
    if input_path.stat().st_size == 0:
      continue
    with gzip.open(input_path, "rt", newline="") as input_handle:
      header = next(csv.reader(input_handle), None)
    if header is None:
      continue
    for column in header:
      if column not in seen_columns:
        merged_header.append(column)
        seen_columns.add(column)

  if not merged_header:
    raise ValueError(f"No diagnostic columns were available for {output_path}")

  with gzip.open(output_path, "wt", newline="") as output_handle:
    writer = csv.DictWriter(output_handle, fieldnames=merged_header)
    writer.writeheader()
    for input_path in input_paths:
      if input_path.stat().st_size == 0:
        continue
      with gzip.open(input_path, "rt", newline="") as input_handle:
        reader = csv.DictReader(input_handle)
        for row in reader:
          writer.writerow(row)
          n_rows += 1

  if n_rows == 0:
    raise ValueError(f"No diagnostic rows were available for {output_path}")


def parse_args(argv=None):
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--real-output", required=True)
  parser.add_argument("--artificial-output", required=True)
  parser.add_argument("--summary-output", required=True)
  parser.add_argument("--real-inputs", nargs="+", required=True)
  parser.add_argument("--artificial-inputs", nargs="+", required=True)
  parser.add_argument("--summary-inputs", nargs="+", required=True)
  return parser.parse_args(argv)


def main(argv=None):
  args = parse_args(argv)
  merge_csv_gz(args.real_inputs, args.real_output)
  merge_csv_gz(args.artificial_inputs, args.artificial_output)
  merge_csv_gz(args.summary_inputs, args.summary_output)


if __name__ == "__main__":
  main()
