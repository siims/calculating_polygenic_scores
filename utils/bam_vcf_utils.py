import gzip
import re
import sqlite3
from collections import defaultdict
from pathlib import Path
from typing import Any, List, Union, Dict, Set, Optional

import numpy as np
import pandas as pd
import rsidx

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
    "chr1": "1",
    "chr2": "2",
    "chr3": "3",
    "chr4": "4",
    "chr5": "5",
    "chr6": "6",
    "chr7": "7",
    "chr8": "8",
    "chr9": "9",
    "chr10": "10",
    "chr11": "11",
    "chr12": "12",
    "chr13": "13",
    "chr14": "14",
    "chr15": "15",
    "chr16": "16",
    "chr17": "17",
    "chr18": "18",
    "chr19": "19",
    "chr20": "20",
    "chr21": "21",
    "chr22": "22",
    "chrX": "X",
    "chrY": "Y",
    "chrM": "MT",
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
    "MT": "MT",
}


def is_alignment_supported(alignment_data):
    supported_reference_genomes = ["grch38.p7"]
    for supported_reference_genome in supported_reference_genomes:
        if supported_reference_genome in str(alignment_data.header.to_dict()).lower():
            return True
    return False


def get_file_header_line_number(file_name: Union[str, Path], header_pattern: str) -> int:
    is_compressed = Path(file_name).suffix == ".gz"
    if is_compressed:
        csv_open = gzip.open
        line_parser = lambda text: line.decode("utf-8")
    else:
        csv_open = open
        line_parser = lambda text: text

    with csv_open(file_name, "r") as f:
        line_number = 0
        for line in f:
            if re.search(header_pattern, line_parser(line)):
                return line_number
            line_number += 1
    raise Exception(f"Couldn't find header in file {file_name}. Expected header: {header_pattern}")


def get_vcf_file_header_line_number(file_name: Union[str, Path]) -> int:
    return get_file_header_line_number(
        file_name=file_name, header_pattern="#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO"
    )


def read_raw_zipped_vcf_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_vcf_file_header_line_number(file_name=file_name)
    result = pd.read_csv(file_name, sep="\s+", skiprows=header_row_number, dtype=str)
    result["POS"] = result["POS"].astype(np.int64)
    return result


def read_raw_zipped_polygenic_score_file(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_file_header_line_number(file_name=file_name, header_pattern="rsID\t")
    result = pd.read_csv(file_name, sep="\t", skiprows=header_row_number, dtype=str)
    result["effect_weight"] = result["effect_weight"].astype(np.float)
    result.rename(columns={"rsID": "rsid"}, inplace=True)
    return result[["rsid", "effect_weight"]]


def read_raw_zipped_polygenic_score_file_with_chrom_pos(file_name: Union[str, Path]) -> pd.DataFrame:
    header_row_number = get_file_header_line_number(file_name=file_name, header_pattern="chr_name\t")
    result = pd.read_csv(file_name, sep="\t", skiprows=header_row_number, dtype=str)
    result["effect_weight"] = result["effect_weight"].astype(np.float)
    # cast to float before because of known bug https://github.com/pandas-dev/pandas/issues/25472
    result["chr_position"] = result["chr_position"].astype("float").astype("Int64")
    result.rename(columns={"chr_name": "chrom", "chr_position": "pos"}, inplace=True)

    hg_build = _get_pgs_file_human_genome_build(file_name)

    result = result[["chrom", "pos", "effect_weight", "reference_allele", "effect_allele"]]
    result.attrs["metadata"] = {"hg_build": hg_build}
    return result


def _get_pgs_file_human_genome_build(file_name):
    is_compressed = Path(file_name).suffix == ".gz"
    if is_compressed:
        csv_open = gzip.open
        line_parser = lambda text: line.decode("utf-8")
    else:
        csv_open = open
        line_parser = lambda text: text
    is_hg19, is_hg38 = False, False
    with csv_open(file_name, "r") as pgs_file:
        for line in pgs_file:
            line = line_parser(line)
            if "grch37" in line.lower():
                is_hg19 = True
            if "grch38" in line.lower():
                is_hg38 = True
    assert (
        not is_hg19 or not is_hg38
    ), f"need to know what human genome build is used, but not known for pgs {file_name}"
    return "hg19" if is_hg19 else "hg38"


def load_vcf_to_df(vcf_files: List[Union[str, Path]], cache_file_name: str = "data/vcf_records.parquet"):
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
    raise AssertionError(
        f"Chomosome {chrom} not found for known contigs {KNOWN_CONTIG_CHROM_MAP.keys()} "
        f"from all alignment data contigs {alignment_data.header.references}"
    )


def get_chrom_reads_in_pos(
    alignment_data, positions: Set[int], contig: Optional[str] = None, chrom: Optional[Union[int, str]] = None
) -> Dict[int, List[str]]:
    if contig is None:
        contig = _get_contig(alignment_data=alignment_data, chrom=str(chrom))
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


def get_chrom_reads_in_range(
    alignment_data, start: int, stop: int, contig: Optional[str] = None, chrom: Optional[Union[int, str]] = None
) -> Dict[int, List[str]]:
    if contig is None:
        contig = _get_contig(alignment_data=alignment_data, chrom=str(chrom))
    sequence = defaultdict()
    for pileup_column in alignment_data.pileup(contig, start, stop):
        # print ("\ncoverage at base %s = %s" % (pileup_column.pos, pileup_column.n), "pileups", len(pileup_column.pileups))
        pos = pileup_column.pos + 1  # NOTE: pileup is 0 based, thus +1 is needed
        try:
            sequence[pos] = _get_reads_in_current_position(pileupcolumn=pileup_column)
        except RuntimeWarning:
            ...
            # print(f"Chromosome {chrom} position {pos} does not have any READS")
    return sequence


def get_read_values_for_allele(
    alignment_data, pos: int, chrom: Optional[Union[int, str]] = None, contig: Optional[str] = None
) -> Dict[int, List[str]]:
    """
    Need either chrom or contig. If both provided uses contig by default.
    :return: dictionary with pos as keys and list of nucleotide values as values
    """
    chrom = str(chrom)
    if contig is None:
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


def calc_genotype_for_chrom_snp_reads(
    chrom_snp_reads: Dict[int, List[str]], chrom: str, sex: str = "male"
) -> pd.DataFrame:
    chrom_snp_genotypes = defaultdict()
    for pos, reads in chrom_snp_reads.items():
        chrom_snp_genotypes[pos] = genotype_from_reads(reads, chrom=chrom, sex=sex)

    chrom_snp_genotypes_df = pd.DataFrame.from_dict(chrom_snp_genotypes, orient="index", columns=["genotype"])
    chrom_snp_genotypes_df.index.set_names("pos", inplace=True)
    chrom_snp_genotypes_df.reset_index(inplace=True)
    return chrom_snp_genotypes_df


def genotype_from_reads(reads, chrom: str, sex: str = "male"):
    counts = {"A": 0, "C": 0, "G": 0, "T": 0, "D": 0}
    for read in reads:
        counts[read] += 1
    sorted_count_keys = sorted(counts, key=counts.__getitem__, reverse=True)
    sorted_count_values = [counts[k] for k in sorted_count_keys]
    if sex == "male" and (chrom in ["X", "Y"]):
        return f"{sorted_count_keys[0]}"
    else:
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
        allele_read_values = get_read_values_for_allele(alignment_data, chrom=chrom, pos=int(pos))

        chromosome_read_values[chrom] = {**chromosome_read_values[chrom], **allele_read_values}
    return chromosome_read_values


def calc_genotypes(alignment_data, loci_df: pd.DataFrame) -> pd.DataFrame:
    chromosome_read_values = calculate_chromosome_read_values(alignment_data=alignment_data, loci_df=loci_df)

    seq = pd.DataFrame(columns=["chr", "pos", "genotype"])
    for chrom, pos_reads in chromosome_read_values.items():
        for pos, reads in pos_reads.items():
            allele = genotype_from_reads(reads, chrom=chrom)
            seq = seq.append({"chr": chrom, "pos": pos, "genotype": allele}, ignore_index=True)
    seq["genotype"] = seq["genotype"].astype(str)
    return seq


def get_my_genotypes_for_pgs(alignment_data, pgs_df: pd.DataFrame, cache_file_name: str) -> pd.DataFrame:
    assert is_alignment_supported(alignment_data)
    cache_file = f"data/my_genotype_in_pos_hg38/{cache_file_name}"
    if not Path(cache_file).exists():
        my_genotypes = calc_genotypes(alignment_data=alignment_data, loci_df=pgs_df)
        my_genotypes.to_csv(cache_file, index=False)
    else:
        my_genotypes = pd.read_csv(cache_file, index_col=None)
    return my_genotypes


def merge_pgs_with_my_genotype(pgs_df: pd.DataFrame, genome_df: pd.DataFrame) -> pd.DataFrame:
    merged_df = genome_df.merge(pgs_df, left_on=["chr", "pos"], right_on=["chr_name", "chr_position"])
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
        return genotype_from_reads(reads[pos], chrom=chrom)
    else:
        raise Exception(f"no reads found for chr{chrom}:{pos}")


def search_for_rsids(rsids: List[str], my_vcf_file: str) -> List[str]:
    if Path(my_vcf_file).suffix == ".gz":
        my_vcf_file = my_vcf_file[:-3]

    file_my_vcf_tabix_indexed = my_vcf_file + ".gz"
    file_my_vcf_rsidx_indexed = my_vcf_file + ".rsidx"
    with sqlite3.connect(file_my_vcf_rsidx_indexed) as db:
        return list(rsidx.search.search(rsids, db, file_my_vcf_tabix_indexed))


def calc_alt_contigs_to_use(
    region_contig_read_counts: Dict[str, Dict[str, int]], assembly_metadata_df: pd.DataFrame
) -> pd.DataFrame:

    alt_contigs_to_use = pd.DataFrame(columns=["chrom", "start", "stop", "contig", "region"])

    for region, contig_read_count in region_contig_read_counts.items():
        region_metadata_df = assembly_metadata_df[assembly_metadata_df["region_name"] == region]
        chrom = region_metadata_df["chromosome"].iloc[0]
        chrom_start = region_metadata_df["chromosome_start"].iloc[0]
        chrom_stop = region_metadata_df["chromosome_stop"].iloc[0]

        regions_contig_with_highest_coverage = sorted(contig_read_count.items(), key=lambda item: item[1])[-1]
        current_contig = regions_contig_with_highest_coverage[0]
        if current_contig != "main":
            alt_contigs_to_use = alt_contigs_to_use.append(
                {"chrom": chrom, "start": chrom_start, "stop": chrom_stop, "contig": current_contig, "region": region},
                ignore_index=True,
            )
    return alt_contigs_to_use
