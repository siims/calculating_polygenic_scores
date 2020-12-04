import re
import sqlite3
from pathlib import Path
from typing import Union, Iterable


def create_snp_db_schema(snp_db_file):
    conn = sqlite3.connect(snp_db_file)
    cursor = conn.cursor()
    cursor.execute('''create table if not exists
    	all_snp_pos
    (
    	chrom text,
    	pos int
    )''')
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


def chunked_insert(conn, all_values: Iterable, chunk_size=10000):
    cursor = conn.cursor()
    all_values = list(all_values)
    print(f"Inserting chrom {all_values[0][0]} values, totalling {len(all_values)}")
    for i in range(0, len(all_values), chunk_size):
        if i + chunk_size > len(all_values):
            cursor.execute(f"INSERT INTO all_snp_pos (chrom,pos) VALUES {all_values[i:].__repr__()[1:-1]};")
        else:
            cursor.execute(
                f"INSERT INTO all_snp_pos (chrom,pos) VALUES {all_values[i:i + chunk_size].__repr__()[1:-1]};")
        conn.commit()


def database_is_already_populated(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM all_snp_pos limit 1")
    res = cursor.fetchall()
    return len(res) != 0


def persist_all_snps_to_db(conn, file_name: Union[str, Path]) -> None:
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
                snps.add((current_chrom, line_parts[1]))

    # persist also final snps
    chunked_insert(conn=conn, all_values=snps)
