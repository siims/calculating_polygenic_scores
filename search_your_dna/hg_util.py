import pandas as pd
from search_your_dna.util import get_file_header_line_number


def get_assembly_metadata_df(assembly_report_file: str, assembly_regions_file: str) -> pd.DataFrame:
    assembly_report_file_header_line_number = get_file_header_line_number(
        assembly_report_file, header_pattern="# Sequence-Name\t"  # header starts with `# Sequence-Name`
    )
    assembly_report_df = pd.read_csv(
        assembly_report_file, sep="\t", skiprows=assembly_report_file_header_line_number, dtype=str
    )
    assembly_report_df.columns = assembly_report_df.columns.str.lower()
    assembly_report_df.columns = assembly_report_df.columns.str.replace("-", "_")
    assembly_report_df.columns = assembly_report_df.columns.str.replace("# ", "")
    assembly_regions_file_header_line_number = get_file_header_line_number(
        assembly_regions_file, header_pattern="# Region-Name\t"  # header starts with `# Region-Name`
    )
    assembly_regions_df = pd.read_csv(
        assembly_regions_file, sep="\t", skiprows=assembly_regions_file_header_line_number, dtype=str
    )
    assembly_regions_df.columns = assembly_regions_df.columns.str.lower()
    assembly_regions_df.columns = assembly_regions_df.columns.str.replace("-", "_")
    assembly_regions_df.columns = assembly_regions_df.columns.str.replace("# ", "")
    assembly_regions_df["chromosome_start"] = assembly_regions_df["chromosome_start"].astype("float").astype("int64")
    assembly_regions_df["chromosome_stop"] = assembly_regions_df["chromosome_stop"].astype("float").astype("int64")
    assembly_metadata_df = pd.merge(
        assembly_regions_df,
        assembly_report_df,
        how="inner",  # inner join to exclude unlocalized contigs
        left_on="scaffold_genbank_accn",
        right_on="genbank_accn",
    )
    return assembly_metadata_df[assembly_metadata_df["ucsc_style_name"] != "na"]
