#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from Bio import SeqIO

def setup_args():
    """
    设置命令行参数
    """
    parser = argparse.ArgumentParser(description='Reorder and summarize Flye assembly contigs.')
    parser.add_argument('--fasta_file', type=str, required=True, help='Path to the Flye assembly.fasta file.')
    parser.add_argument('--output_fasta', type=str, default="reordered_flye_assembly.fasta", help='Path to output the reordered fasta file.')
    return parser.parse_args()

def process_fasta(args):


    # 读取FASTA文件，并过滤掉小于2000bp的序列，除非是成环序列
    records = [rec for rec in SeqIO.parse(args.fasta_file, "fasta")]
    
    # 按照长度对记录进行排序（从长到短）
    records.sort(key=lambda x: len(x.seq), reverse=True)

    # 将排序后的fasta写出
    SeqIO.write(records, args.output_fasta, "fasta")


def main():
    args = setup_args()
    process_fasta(args)

if __name__ == '__main__':
    main()
