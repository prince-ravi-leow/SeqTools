#!/usr/bin/env python3

"""
Tools for converting sequence and diamond blast data to and from pd.DataFrame
"""

import gzip
import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO, SeqRecord


def SeqIO_load(input_file: str) -> SeqRecord:
    if Path(input_file).suffix == ".gz":
        base_file, gz_ext = os.path.splitext(str(input_file))
        file_type = base_file.split(".")[1]
        with gzip.open(input_file, "rt") as input_handle:
            records = list(SeqIO.parse(input_handle, file_type))
    else:
        file_type = str(input_file).split(".")[-1]
        records = list(SeqIO.parse(input_file, file_type))
    return records


def fasta_to_df(input_file: str) -> pd.DataFrame:
    records = SeqIO_load(input_file)
    df = pd.DataFrame(
        {
            "protein": [record.id for record in records],
            "seq": [str(record.seq) for record in records],
        }
    )
    return df.astype("string")


def df_to_fasta(
    df: pd.DataFrame, output_file: str, embed_metadata: bool = False
) -> None:
    if not embed_metadata:
        fasta_str = "\n".join(">" + df["protein"] + "\n" + df["seq"])
    else:
        fasta_str = "\n".join(
            ">"
            + df["response"]
            + ":"
            + df["protein"]
            + ":"
            + df["pann"]
            + ":"
            + df["pcat"]
            + "\n"
            + df["seq"]
        )
    with open(output_file, "w") as f:
        f.write(fasta_str)


def dblast_to_df(dblast_align_file: str) -> pd.DataFrame:
    df = pd.read_table(
        dblast_align_file,
        header=None,
        names=[
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )
    return df[df["qseqid"] != df["sseqid"]]
