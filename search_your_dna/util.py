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

KNOWN_CONTIG_CHROM_MAP = {
    "CM000663.2": "1",
    "CM000664.2": "2",
    "CM000665.2": "3",
    "CM000666.2": "4",
    "CM000667.2": "5",
    "CM000668.2": "6",
    "CM000669.2": "7",
    "CM000670.2": "8",
    "CM000671.2": "9",
    "CM000672.2": "10",
    "CM000673.2": "11",
    "CM000674.2": "12",
    "CM000675.2": "13",
    "CM000676.2": "14",
    "CM000677.2": "15",
    "CM000678.2": "16",
    "CM000679.2": "17",
    "CM000680.2": "18",
    "CM000681.2": "19",
    "CM000682.2": "20",
    "CM000683.2": "21",
    "CM000684.2": "22",
    "CM000685.2": "X",
    "CM000686.2": "Y",
    "J01415.2": "MT",
    "1": "1",
    "2": "2",
    "3": "3",
    "4": "4",
    "5": "5",
    "6": "6",
    "7": "7",
    "8": "8",
    "9": "9",
    "10": "10",
    "11": "11",
    "12": "12",
    "13": "13",
    "14": "14",
    "15": "15",
    "16": "16",
    "17": "17",
    "18": "18",
    "19": "19",
    "20": "20",
    "21": "21",
    "22": "22",
    "X": "X",
    "Y": "Y",
    "MT": "MT"
}


def is_alignment_supported(alignment_data):
    supported_reference_genomes = [
        "grch38.p7"
    ]
    for supported_reference_genome in supported_reference_genomes:
        if supported_reference_genome in str(alignment_data.header.to_dict()).lower():
            return True
    return False


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
        header_pattern="rsID\t"
    )


def read_raw_zipped_vcf_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_vcf_file_header_line_number(file_name=file_name)
    result = pd.read_csv(file_name, sep="\s+", skiprows=header_row_number, dtype=str)
    result["POS"] = result["POS"].astype(np.int64)
    return result


def read_raw_zipped_polygenic_score_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_polygenic_score_file_header_line_number(file_name=file_name)
    result = pd.read_csv(file_name, sep="\t", skiprows=header_row_number, dtype=str)
    result["effect_weight"] = result["effect_weight"].astype(np.float)
    if "chr_position" in result.columns:
        result["chr_position"] = result["chr_position"].astype('float').astype("Int64")  # cast to float before because of known bug https://github.com/pandas-dev/pandas/issues/25472
    result.rename(columns={"rsID": "rsid"}, inplace=True)
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


def _get_contig(alignment_data, chrom):
    for contig in alignment_data.header.references:
        if KNOWN_CONTIG_CHROM_MAP.get(contig) == chrom:
            return contig
    raise AssertionError(f"Chomosome {chrom} not found for known contigs {KNOWN_CONTIG_CHROM_MAP.keys()} "
                         f"from all alignment data contigs {alignment_data.header.references}")


def get_chrom_reads_in_pos(alignment_data, chrom: Union[int, str], positions: Set[int]) -> Dict[int, List[str]]:
    chrom = str(chrom)
    contig = _get_contig(alignment_data=alignment_data, chrom=chrom)
    sequence = defaultdict()
    for pileup_column in alignment_data.pileup(contig, 0):
        # print ("\ncoverage at base %s = %s" % (pileup_column.pos, pileup_column.n), "pileups", len(pileup_column.pileups))
        pos = pileup_column.pos + 1  # NOTE: pileup is 0 based, thus +1 is needed
        if pos in positions:
            try:
                sequence[pos] = _get_reads_in_current_position(pileupcolumn=pileup_column)
            except RuntimeWarning:
                ...
                # print(f"Chromosome {chrom} position {pos} does not have any READS")
    return sequence


def get_read_values_for_allele(alignment_data, chrom: Union[int, str], pos: int) -> Dict[int, List[str]]:
    chrom = str(chrom)
    contig = _get_contig(alignment_data=alignment_data, chrom=chrom)
    sequence = defaultdict()
    for pileup_column in alignment_data.pileup(contig, pos - 1, pos + 1):
        # print ("\ncoverage at base %s = %s" % (pileup_column.pos, pileup_column.n), "pileups", len(pileup_column.pileups))
        if pos == pileup_column.pos + 1:  # NOTE: pileup is 0 based, thus +1 is needed
            try:
                sequence[pos] = _get_reads_in_current_position(pileupcolumn=pileup_column)
            except RuntimeWarning:
                ...
                # print(f"Chromosome {chrom} position {pos} does not have any READS")
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
    if sorted_count_values[0] / sum(sorted_count_values) > 0.9:  # 0.9 is threshold for excluding wrong reads
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
        cache_file_name: str
) -> pd.DataFrame:
    assert is_alignment_supported(alignment_data)
    cache_file = f"data/my_genotype_in_pos_hg38/{cache_file_name}"
    if not Path(cache_file).exists():
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


def get_my_snps_for_chromosome(alignment_data, snp_db_file: str, chrom: str) -> pd.DataFrame:
    assert is_alignment_supported(alignment_data=alignment_data) # db has only hg38
    cache_file = Path(f"data/my_genotype_in_pos_hg38/my_chrom_{chrom}_snp.csv")
    if cache_file.exists():
        res_df = pd.read_csv(cache_file, index_col=None)
    else:
        _conn = sqlite3.connect(snp_db_file)
        snp_for_chrom = pd.read_sql_query(f"SELECT * FROM all_snp_pos WHERE chrom = '{chrom}'", con=_conn)
        my_chrom_snp_reads = get_chrom_reads_in_pos(
            alignment_data=alignment_data, chrom=chrom, positions=set(snp_for_chrom["pos"].to_list())
        )
        print(f"in database #{snp_for_chrom.shape[0]} SNPs; in my genome found #{len(my_chrom_snp_reads)}")
        res_df = calc_genotype_for_chrom_snp_reads(my_chrom_snp_reads)
        res_df.to_csv(cache_file, index=False)
    res_df["chrom"] = chrom
    return res_df


