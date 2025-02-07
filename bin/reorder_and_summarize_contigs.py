#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reorder_and_summarize_contigs.py

Description:
    This script processes assembly outputs from Flye to reorder contigs by length,
    rename them sequentially, and update statistics file accordingly. It also calculates
    the total number of sequence chunks based on a given chunk size and identifies the
    longest contig.

Usage:
    python reorder_and_summarize_contigs.py <fasta_file> <stats_file> <chunk_size>

Arguments:
    fasta_file: Path to the Flye assembly.fasta file.
    stats_file: Path to the Flye assembly statistics TSV file (flye_out_stat.tsv).
    chunk_size: Integer specifying the size of chunks for contig division.

Outputs:
    sorted_<fasta_file>: Reordered and renamed FASTA file.
    updated_<stats_file>: Updated statistics file with new contig names.
    assembly_summary.txt: File containing the name of the longest contig and total chunk count.
"""

import argparse
import pandas as pd
from Bio import SeqIO
import math

def get_long_reads_circ(stats_file):
    """
    从统计文件中读取并返回所有成环的contig名称。
    
    参数:
    - stats_file: 统计文件的路径 (包含'circ.'列的tsv文件)
    
    返回:
    - long_reads_circ: 成环contig的名称列表
    """
    # 读取统计文件
    stats_df = pd.read_csv(stats_file, sep='\t')
    
    # 获取成环(contig中circ.列值为'Y')的contig名称
    long_reads_circ = stats_df[stats_df['circ.'] == 'Y']['seq_name'].tolist()
    
    return long_reads_circ

def process_fasta_and_stats(fasta_file, stats_file, chunk_size, out_prefix):
    # 读取统计文件
    stats_df = pd.read_csv(stats_file, sep='\t')

    # 获取成环(contig中circ.列值为'Y')的contig名称
    long_reads_circ = stats_df[stats_df['circ.'] == 'Y']['#seq_name'].tolist()

    # 读取FASTA文件,并过滤掉小于2000bp的序列，除非是成环序列
    records = [rec for rec in SeqIO.parse(fasta_file, "fasta") if len(rec.seq) >= 2000 or rec.id in long_reads_circ]
    
    # 获取序列的长度并排序
    records.sort(key=lambda x: len(x.seq), reverse=True)

    # 重命名并创建新旧名称映射
    name_mapping = {}
    new_names = []
    for i, record in enumerate(records):
        new_name = f"contig_{i+1}"
        new_names.append(new_name)
        name_mapping[record.id] = new_name
        record.id = new_name
        record.name = new_name
        record.description = ''

    # 将排序后的fasta写出
    SeqIO.write(records, out_prefix+"_sorted_assembly.fa", "fasta")

    # 读取统计文件
    # stats_df = pd.read_csv(stats_file, sep='\t')

    # 根据name_mapping过滤stats_df，只保留过滤后的记录
    stats_df = stats_df[stats_df['#seq_name'].isin(name_mapping.keys())]

    # 插入新的列
    stats_df['new_seq_name'] = stats_df['#seq_name'].map(name_mapping)

    # 保存新的统计文件
    stats_df.to_csv(out_prefix+"_updated_flye_stat.tsv", sep='\t', index=False)

    # 输出最长的contig名称和总块数到文件
    longest_contig_name = new_names[0]
    total_chunks = sum(math.ceil(len(record.seq) / chunk_size) for record in records)
    
    with open(out_prefix+"_assembly_summary.txt", "w") as f:
        f.write(f"longest_contig_name:{longest_contig_name}\n")
        f.write(f"chunk_num:{total_chunks}\n")

def setup_args():
    parser = argparse.ArgumentParser(description='Reorder and summarize contigs from Flye assembly outputs.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fasta_file', type=str, help='Path to the Flye assembly.fasta file.', required=True)
    parser.add_argument('--stats_file', type=str, help='Path to the Flye assembly statistics TSV file.', required=True)
    parser.add_argument('--chunk_size', type=int, help='Chunk size for contig division.', required=True)
    parser.add_argument('--prefix', type=str, help='Out file prefix.', default="sampleID")
    # parser.add_argument('--output_fasta', type=str, default="reordered_flye_assembly.fasta", help='Path to output the reordered fasta file.')
    # parser.add_argument('--output_stats', type=str, default="updated_flye_stat.tsv", help='Path to output the updated statistics file.')
    # parser.add_argument('--summary_file', type=str, default="assembly_summary.txt", help='Path to output the summary file.')
    return parser.parse_args()


if __name__ == '__main__':
    args = setup_args()
    # 使用这个函数
    # long_reads_circ = get_long_reads_circ(args.stats_file) 
    process_fasta_and_stats(args.fasta_file, args.stats_file, args.chunk_size, args.prefix)

# process_fasta_and_stats(args.fasta_file, args., args.chunk_size)
