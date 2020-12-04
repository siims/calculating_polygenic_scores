import gzip
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, List, Union, Dict, Set

import numpy as np
import pandas as pd

chrom_list = [
    '1',
    '2',
    '3',
    '4',
    '5',
    '6',
    '7',
    '8',
    '9',
    '10',
    '11',
    '12',
    '13',
    '14',
    '15',
    '16',
    '17',
    '18',
    '19',
    '20',
    '21',
    '22',
    'MT',
    'X',
    'Y'
]


def get_file_header_line_number(file_name: Union[str, Path], header_pattern: str) -> int:
    with gzip.open(str(file_name), "r") as f:
        line_number = 0
        for line in f:
            if re.search(header_pattern, line.decode("utf-8")):
                return line_number
            line_number += 1
    raise Exception(f"Couldn't find header in file {file_name}. Expected header: {header_pattern}")


def get_vcf_file_header_line_number(file_name: Union[str, Path]) -> int:
    return get_file_header_line_number(
        file_name=file_name,
        header_pattern="#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO"
    )


def get_polygenic_score_file_header_line_number(file_name: Union[str, Path]) -> int:
    return get_file_header_line_number(
        file_name=file_name,
        header_pattern="rsID\s+chr_name\s+chr_position\s+effect_allele"
    )


def read_raw_zipped_vcf_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_vcf_file_header_line_number(file_name=file_name)
    result = pd.read_csv(file_name, sep="\s+", skiprows=header_row_number, dtype=str)
    result["POS"] = result["POS"].astype(np.int64)
    return result


def read_raw_zipped_polygenic_score_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_polygenic_score_file_header_line_number(file_name=file_name)
    result = pd.read_csv(file_name, sep="\s+", skiprows=header_row_number, dtype=str)
    result["effect_weight"] = result["effect_weight"].astype(np.float)
    result["chr_name"] = result["chr_name"].astype(np.int64)
    result["chr_position"] = result["chr_position"].astype(np.int64)
    return result


def load_vcf_to_df(vcf_files: List[Union[str, Path]], cache_file_name: str = "data/vcf_records.parquet.gz"):
    if Path(cache_file_name).exists():
        return pd.read_parquet(cache_file_name)

    dfs = []
    for vcf_file_path in vcf_files:
        print(f"Reading in source vcf file {vcf_file_path}")
        dfs.append(read_raw_zipped_vcf_file(vcf_file_path))
    raw_vcf_data = pd.concat(dfs, ignore_index=True)
    raw_vcf_data.to_parquet(cache_file_name)
    return raw_vcf_data


def load_polygenic_score_file_to_df(file_name: Union[str, Path]) -> pd.DataFrame:
    return read_raw_zipped_polygenic_score_file(file_name=file_name)


# section 2

def _get_reads_in_current_position(pileupcolumn):
    if len(pileupcolumn.pileups) == 0:
        raise RuntimeWarning("No reads found")
    reads_at_current_position = []
    for pileupread in pileupcolumn.pileups:
        if pileupread.is_del:
            reads_at_current_position.append("D")
        else:
            # print(pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position])
            # print('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
            reads_at_current_position.append(pileupread.alignment.query_sequence[pileupread.query_position])
    return reads_at_current_position


def get_chrom_reads_in_pos(alignment_data, chrom: Union[int, str], positions: Set[int]) -> Dict[int, List[str]]:
    sequence = defaultdict()
    for pileupcolumn in alignment_data.pileup(str(chrom), 0):
        # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n), "pileups", len(pileupcolumn.pileups))
        pos = pileupcolumn.pos + 1  # NOTE: pileup is 0 based, thus +1 is needed
        if pos in positions:
            try:
                sequence[pos] = _get_reads_in_current_position(pileupcolumn=pileupcolumn)
            except RuntimeWarning:
                print(f"Chromosome {chrom} position {pos} does not have any READS")
    return sequence


def get_read_values_for_allele(alignment_data, chrom: Union[int, str], pos: int) -> Dict[int, List[str]]:
    sequence = defaultdict()
    for pileupcolumn in alignment_data.pileup(str(chrom), pos - 1, pos + 1):
        # print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n), "pileups", len(pileupcolumn.pileups))
        if pos == pileupcolumn.pos + 1:  # NOTE: pileup is 0 based, thus +1 is needed
            try:
                sequence[pos] = _get_reads_in_current_position(pileupcolumn=pileupcolumn)
            except RuntimeWarning:
                print(f"Chromosome {chrom} position {pos} does not have any READS")
    return sequence


def calc_genotype_for_chrom_snp_reads(chrom_snp_reads: Dict[int, List[str]]) -> pd.DataFrame:
    chrom_snp_genotypes = defaultdict()
    for pos, reads in chrom_snp_reads.items():
        chrom_snp_genotypes[pos] = genotype_from_reads(reads)

    chrom_snp_genotypes_df = pd.DataFrame.from_dict(chrom_snp_genotypes, orient="index", columns=["genotype"])
    chrom_snp_genotypes_df.index.set_names("pos", inplace=True)
    chrom_snp_genotypes_df.reset_index(inplace=True)
    return chrom_snp_genotypes_df


def genotype_from_reads(reads):
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "D": 0}
    for read in reads:
        counts[read] += 1
    sorted_count_keys = sorted(counts, key=counts.__getitem__, reverse=True)
    sorted_count_values = [counts[k] for k in sorted_count_keys]
    if sorted_count_values[0] / sum(sorted_count_values) > 0.9:
        return f"{sorted_count_keys[0]}{sorted_count_keys[0]}"
    else:
        return f"{sorted_count_keys[0]}{sorted_count_keys[1]}"


def calculate_chromosome_read_values(alignment_data, loci_df: pd.DataFrame) -> Dict[str, Any]:
    chromosome_read_values = defaultdict()
    for entry in loci_df.to_dict(orient="records"):
        chrom = entry["chr_name"]
        pos = entry["chr_position"]
        if chrom not in chromosome_read_values:
            chromosome_read_values[chrom] = {}
        allele_read_values = get_read_values_for_allele(alignment_data, chrom, int(pos))

        chromosome_read_values[chrom] = {**chromosome_read_values[chrom], **allele_read_values}
    return chromosome_read_values


def calc_genotypes(alignment_data, loci_df: pd.DataFrame) -> pd.DataFrame:
    chromosome_read_values = calculate_chromosome_read_values(alignment_data=alignment_data, loci_df=loci_df)

    seq = pd.DataFrame(columns=["chr", "pos", "genotype"])
    for chrom, pos_reads in chromosome_read_values.items():
        for pos, reads in pos_reads.items():
            allele = genotype_from_reads(reads)
            seq = seq.append({"chr": chrom, "pos": pos, "genotype": allele}, ignore_index=True)
    seq["genotype"] = seq["genotype"].astype(str)
    return seq


def get_my_genotypes_for_pgs(
        alignment_data,
        pgs_df: pd.DataFrame,
        cache_file_name: str,
        use_filter: bool = False
) -> pd.DataFrame:
    cache_file = f"data/{cache_file_name}"
    if not Path(cache_file).exists():
        if use_filter:
            pgs_df_abs_weight = np.abs(pgs_df["effect_weight"])
            pgs_df = pgs_df[pgs_df_abs_weight > pgs_df_abs_weight.mean()]
        my_genotypes = calc_genotypes(alignment_data=alignment_data, loci_df=pgs_df)
        my_genotypes.to_csv(cache_file, index=False)
    else:
        my_genotypes = pd.read_csv(cache_file, index_col=None)
    return my_genotypes


def merge_pgs_with_my_genotype(pgs_df: pd.DataFrame, my_genome_df: pd.DataFrame) -> pd.DataFrame:
    merged_df = my_genome_df.merge(pgs_df, left_on=["chr", "pos"], right_on=["chr_name", "chr_position"])
    return merged_df[["chr", "pos", "genotype", "effect_allele", "reference_allele", "effect_weight"]]


def filter_out_none_effect_alleles(merged_pgs_with_my_genotype):
    return merged_pgs_with_my_genotype[
        (merged_pgs_with_my_genotype["genotype"].map(lambda x: x[0]) == merged_pgs_with_my_genotype["effect_allele"])
        | (merged_pgs_with_my_genotype["genotype"].map(lambda x: x[1]) == merged_pgs_with_my_genotype["effect_allele"])
        ]


def filter_out_effect_alleles(merged_pgs_with_my_genotype):
    return merged_pgs_with_my_genotype[
        (merged_pgs_with_my_genotype["genotype"].map(lambda x: x[0]) != merged_pgs_with_my_genotype["effect_allele"])
        & (merged_pgs_with_my_genotype["genotype"].map(lambda x: x[1]) != merged_pgs_with_my_genotype["effect_allele"])
        ]


def get_genotype_for_chrom_pos(alignment_data, chrom: str, pos: int) -> str:
    reads = get_read_values_for_allele(alignment_data=alignment_data, chrom=chrom, pos=pos)
    if len(reads) != 0:
        return genotype_from_reads(reads[pos])
    else:
        raise Exception(f"no reads found for chr{chrom}:{pos}")


def get_my_snps_for_chromosome(snp_db_file: str, chrom: str) -> pd.DataFrame:
    cache_file = Path(f"data/my_chrom_{chrom}_snp.csv")
    if cache_file.exists():
        res_df = pd.read_csv(cache_file, index_col=None)
    else:
        _conn = sqlite3.connect(snp_db_file)
        snp_for_chrom = pd.read_sql_query(f"SELECT * FROM all_snp_pos WHERE chrom = '{chrom}'", con=_conn)
        my_chrom_snp_reads = get_chrom_reads_in_pos(chrom, set(snp_for_chrom["pos"].to_list()))
        print(f"in database #{snp_for_chrom.shape[0]} SNPs; in my genome found #{len(my_chrom_snp_reads)}")
        res_df = calc_genotype_for_chrom_snp_reads(my_chrom_snp_reads)
        res_df.to_csv(cache_file, index=False)
    res_df["chrom"] = chrom
    return res_df
