import gzip
import pickle
import re
from pathlib import Path
from typing import List, Union, Optional

import numpy as np
import pandas as pd
import vcf

# %%
from search_your_dna.util import load_vcf_to_df

pgs_df = pd.read_csv("/home/s/src/search_your_dna/data/PGS000318.txt", sep="\t", skiprows=9)

# %%

cache_file_name = "/home/s/src/search_your_dna/data/tmp_vcf_records.pkl"
vcf_file_paths = [Path(f) for f in [
    "/home/s/Dropbox/Siim/health/genetest_2020/GFX0237425.cnv.vcf.gz",
    "/home/s/Dropbox/Siim/health/genetest_2020/GFX0237425.filtered.indel.vcf.gz",
    "/home/s/Dropbox/Siim/health/genetest_2020/GFX0237425.filtered.snp.vcf.gz",
    "/home/s/Dropbox/Siim/health/genetest_2020/GFX0237425.sv.vcf.gz"
]]


def cache_vcf_records(records: List, overwrite: bool = False) -> None:
    if overwrite and Path(cache_file_name).exists():
        raise Exception("Won't overwrite existing cache file. Please set flag specifically to do that")
    with open(cache_file_name, "wb") as f:
        print(f"Caching records to {cache_file_name}")

        pickle.dump(records, f)


def read_vcf_files(vcf_file_paths):
    if Path(cache_file_name).exists():
        print(f"Reading cached files from {cache_file_name}")
        with open(cache_file_name, 'rb') as f:
            records = pickle.load(f)
    else:
        records = []
        for vcf_file_path in vcf_file_paths:
            print(f"Reading in source vcf file {vcf_file_path}")
            vcf_reader = vcf.Reader(filename=str(vcf_file_path))
            for record in vcf_reader:
                records.append(record)
        cache_vcf_records(records=records)
    return records


# records = read_vcf_files(vcf_file_paths=vcf_file_paths)


def vcf_record_to_df(records):
    result = pd.DataFrame(
        columns=["chr", "pos", "ref", "alt", "affected_start", "affected_end", "heterozygosity", "quality"])
    for record in records:
        result.loc[len(result.index)] = [
            record.CHROM,
            record.POS,
            record.REF,
            record.ALT,
            record.affected_start,
            record.affected_end,
            record.heterozygosity,
            record.QUAL,
        ]
    return result


def get_vcf_file_header_line_number(file_name: Union[str, Path]) -> int:
    header_pattern = "#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO\s+FORMAT"  # has one more column with sample number

    with gzip.open(str(file_name), "r") as f:
        line_number = 0
        for line in f:
            if re.search(header_pattern, line.decode("utf-8")):
                return line_number
            line_number += 1
    raise Exception(f"Couldn't find header in file {file_name}. Expected header: {header_pattern}")


def read_raw_zipped_vcf_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_vcf_file_header_line_number(file_name=file_name)
    result = pd.read_csv(file_name, sep="\s+", skiprows=header_row_number, dtype=str)
    result["POS"] = result["POS"].astype(np.int64)
    return result


if __name__ == '__main__':
    zipped_vcf_file = "data/vcf_records.parquet.gz"
    vcf_df = load_vcf_to_df(vcf_file_paths, cache_file_name=zipped_vcf_file)
    # do something
