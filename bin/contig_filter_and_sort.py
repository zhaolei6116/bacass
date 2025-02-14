#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
from Bio import SeqIO
import csv


def select_qualified_contigs(stats_path):
    """
    从统计文件中筛选出合格的 contig，
    筛选条件是，1、第一条contig 是最长的contig ，直接保存，不做判断。 2、其余的contig dep 是需要高于第一条contig dep 的 5分之1。且需要成环。才能保存。
    并返回 contig 名称列表。
    
    Args:
        stats_path (str): 统计文件路径。     
    
    Returns:
        list: 合格的 contig 名称列表。
    """
    
    with open(stats_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        reader = list(reader)
    
    threshold_dep = float(reader[0][4]) / 5   # 这里确定第一条就是最长的，最好的contig，最长的contig是可以保证的，前期通过排序，但这不一定是深度最大的，这可能是一个漏洞。
    qualified_contigs = [line[0] for line in reader[1:] if float(line[4]) >= threshold_dep and line[5] == 'Y']
    qualified_contigs.append(reader[0][0])
    print(qualified_contigs)
    
    return qualified_contigs, header, reader

def read_fna_and_record_order(fna_path, qualified_contigs, out_fna):
    """
    读取 FNA 文件，过滤不合格的contig， 按长度排序， 并重新命名contig fasta序列名字。

    Args:
        fna_path (str): FNA 文件路径。
        qualified_contigs (list): 合格的 contig 名称列表。
        
    Returns:
        tuple: 包含 ()，
            - old_names: 合格的 SeqRecord 列表, 与fna 顺序一致。
            - name_mapping (dict): 记录每个合格 contig 名字对应关系（旧名称 -> 新序号）。
        生成新的fna文件。过滤并排序后的。
    """
    
    # 读取FASTA文件,并过滤掉低dep contig。
    records = [rec for rec in SeqIO.parse(fna_path, "fasta") if  rec.id in qualified_contigs]

    # 获取序列的长度并排序
    records.sort(key=lambda x: len(x.seq), reverse=True)

    # # 重命名并创建新旧名称映射
    # name_mapping = {}
    old_names = []
    for i, record in enumerate(records):
    #     new_name = f"contig_{i+1}"
    #     name_mapping[record.id] = new_name
        old_names.append(record.id)
    #     record.id = new_name
    #     record.name = new_name
    #     record.description = ''
    
    # 将排序后的fasta写出
    SeqIO.write(records, out_fna, "fasta")
    
    return  old_names


def update_stats(old_names, header, reader, output_stats):
    """
    更新统计文件。
    
    Args:
        old_names (list): 原始统计文件contig name。
        name_mapping (dict): contig name old:new。
        header (list): 输出统计文件标题行。
        reader (双重list): [[原始文件行信息][]].
        stats_path: 重新生成的统计文件， 输出文件；
    """

    # 过滤
    filtered_contig = [entry for entry in reader if entry[0] in old_names]

    # 按照 old_names 顺序排序
    order_dict = {name: index for index, name in enumerate(old_names)}
    filtered_contig.sort(key=lambda x: order_dict[x[0]])

    # 替换 old_name
    # for entry in filtered_contig:
    #     old_name = entry[0]
    #     if old_name in name_mapping:
    #         entry[0] = name_mapping[old_name]
    
    
    # 写入新的统计文件
    with open(output_stats, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        filtered_contig.insert(0, header)
        writer.writerows(filtered_contig)

def setup_args():
    parser = argparse.ArgumentParser(description='Reorder and summarize contigs from Flye assembly outputs.',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--contig_stat', type=str, help='Path to the contig stat file.', required=True)
    parser.add_argument('--in_fna', type=str, help='Path to the consense fna file.', required=True)
    parser.add_argument('--out_fna', type=str, help='Out new fna file.', required=True)
    parser.add_argument('--out_stat', type=str, help='Out new stat file.', required=True)
    # parser.add_argument('--output_fasta', type=str, default="reordered_flye_assembly.fasta", help='Path to output the reordered fasta file.')
    # parser.add_argument('--output_stats', type=str, default="updated_flye_stat.tsv", help='Path to output the updated statistics file.')
    # parser.add_argument('--summary_file', type=str, default="assembly_summary.txt", help='Path to output the summary file.')
    return parser.parse_args()

def main():
    args = setup_args()
    qualified_contigs, header, reader = select_qualified_contigs(args.contig_stat)
    old_names = read_fna_and_record_order(args.in_fna, qualified_contigs, args.out_fna)
    # read_fna_and_record_order(args.in_fna, qualified_contigs, args.out_fna)
    update_stats(old_names, header, reader, args.out_stat)


if __name__ == '__main__':
    main()
