import re
import sqlite3
import time

import pandas as pd
from pathlib import Path
from typing import Union, Iterable, Tuple, List


def create_snp_db_schema(snp_db_file):
    conn = sqlite3.connect(snp_db_file)
    cursor = conn.cursor()
    cursor.execute('''create table if not exists
    	all_snp_pos
    (
    	chrom text,
    	pos int,
    	rsid var(25)
    )''')
    conn.commit()
    conn.close()


def create_snp_db_chrom_pos_index(snp_db_file):
    conn = sqlite3.connect(snp_db_file)
    cursor = conn.cursor()
    cursor.execute('''create unique index if not exists
    	all_snp_pos_chrom_pos_index
    on
    	all_snp_pos
    (
    	chrom,
    	pos
    )''')
    conn.commit()
    conn.close()


def create_snp_db_rsid_index(snp_db_file):
    conn = sqlite3.connect(snp_db_file)
    cursor = conn.cursor()
    cursor.execute('''create unique index if not exists
    	all_snp_pos_rsid_index
    on
    	all_snp_pos
    (
    	rsid
    )''')
    conn.commit()
    conn.close()


def chunked_insert(conn, all_values: Iterable[Tuple[str, int, str]], chunk_size=10000):
    cursor = conn.cursor()
    all_values = list(all_values)
    print(f"Inserting chrom {all_values[0][0]} values, totalling {len(all_values)}")
    for i in range(0, len(all_values), chunk_size):
        if i + chunk_size > len(all_values):
            cursor.execute(f"INSERT INTO all_snp_pos (chrom,pos,rsid) VALUES {all_values[i:].__repr__()[1:-1]};")
        else:
            cursor.execute(
                f"INSERT INTO all_snp_pos (chrom,pos,rsid) VALUES {all_values[i:i + chunk_size].__repr__()[1:-1]};")
        conn.commit()


def database_is_already_populated(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM all_snp_pos limit 1")
    res = cursor.fetchall()
    return len(res) != 0


def persist_all_snps_to_db(conn, file_name: Union[str, Path]) -> None:
    if database_is_already_populated(conn):
        print("Database already populated, not populating again.")
        return
    header_pattern = "#CHROM\s+POS\s+ID\s+REF\s+ALT\s+QUAL\s+FILTER\s+INFO"
    snps = set()
    with open(str(file_name), "r") as f:
        passed_header = False
        last_chrom = "1"
        for line_text in f:
            if not passed_header:
                if re.search(header_pattern, line_text):
                    passed_header = True
            else:
                line_parts = line_text.split("\t")
                current_chrom = line_parts[0]
                if last_chrom != current_chrom:
                    chunked_insert(conn=conn, all_values=snps)
                    snps = set()
                    last_chrom = current_chrom
                snps.add((current_chrom, int(line_parts[1]), line_parts[2]))

    # persist also final snps
    chunked_insert(conn=conn, all_values=snps)


def insert_genotype_to_db(snp_db_file: str, genotype_df: pd.DataFrame, genotype_col_name: str = "genotype"):
    _conn = sqlite3.connect(snp_db_file)
    cur = _conn.cursor()
    for entry in genotype_df.to_dict(orient="records"):
        pos = entry["pos"]
        genotype = entry[genotype_col_name]
        chrom = entry["chrom"]
        query = f"UPDATE all_snp_pos SET {genotype_col_name} = '{genotype}' WHERE chrom = '{chrom}' and pos = {pos}"
        cur.execute(query)
    _conn.commit()


def query_my_genotypes_for_rsids(
        snp_db_file: str,
        rsids: List[str],
        genotype_col_name: str = "genotype"
) -> pd.DataFrame:
    _conn = sqlite3.connect(snp_db_file)
    rsids_str = "(" + rsids.__repr__()[1:-1] + ")"
    query = f"SELECT rsid, {genotype_col_name} FROM all_snp_pos WHERE rsid in {rsids_str}"
    print(f"Executing query {query[:100]}.. With number of rsids: {len(rsids)}")
    start = time.time()
    res = pd.read_sql_query(sql=query, con=_conn)
    end = time.time()
    print("\tTook totally:", end - start)
    return res
