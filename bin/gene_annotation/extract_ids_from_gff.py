#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
from collections import defaultdict
import argparse

def extract_ids_from_gff(gff_file):
    """从 GFF3 文件中提取 GO ID、KEGG ID 和 COG ID 并统计出现次数"""
    go_id_counter = defaultdict(int)
    kegg_id_counter = defaultdict(int)
    cog_id_counter = defaultdict(int)

    # 正则表达式匹配 GO ID、KEGG ID 和 COG ID
    go_id_pattern = re.compile(r'GO:\d{7}')
    kegg_id_pattern = re.compile(r'KEGG:K\d{5}')
    cog_id_pattern = re.compile(r'COG\d{4}')

    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):  # 跳过注释行
                continue
            # 查找所有匹配的 ID
            go_ids = go_id_pattern.findall(line)
            kegg_ids = kegg_id_pattern.findall(line)
            cog_ids = cog_id_pattern.findall(line)

            # 统计 GO ID
            for go_id in go_ids:
                go_id_counter[go_id] += 1

            # 统计 KEGG ID
            for kegg_id in kegg_ids:
                kegg_id_counter[kegg_id] += 1

            # 统计 COG ID
            for cog_id in cog_ids:
                cog_id_counter[cog_id] += 1

    return go_id_counter, kegg_id_counter, cog_id_counter

def save_ids_with_counts(counter, output_file):
    """将 ID 及其出现次数保存到文件中"""
    with open(output_file, 'w') as file:
        for id, count in sorted(counter.items()):
            file.write(f"{id}\t{count}\n")

def save_ids_with_counts2(counter, output_file):
    """将 ID 及其出现次数保存到文件中"""
    with open(output_file, 'w') as file:
        for id, count in sorted(counter.items()):
            idk = id.split(":")[1]
            file.write(f"{idk}\t{count}\n")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="从 GFF3 文件中提取 GO ID、KEGG ID 和 COG ID 并统计出现次数")
    parser.add_argument('-i', '--input', required=True, help="输入的 GFF3 文件路径")
    parser.add_argument('--go_output', required=True, help="输出的 GO ID 文件路径")
    parser.add_argument('--kegg_output', required=True, help="输出的 KEGG ID 文件路径")
    parser.add_argument('--cog_output', required=True, help="输出的 COG ID 文件路径")

    # 解析命令行参数
    args = parser.parse_args()

    # 提取 ID 并统计出现次数
    go_id_counter, kegg_id_counter, cog_id_counter = extract_ids_from_gff(args.input)

    # 将结果保存到文件
    save_ids_with_counts(go_id_counter, args.go_output)
    save_ids_with_counts2(kegg_id_counter, args.kegg_output)
    save_ids_with_counts(cog_id_counter, args.cog_output)

    print(f"提取完成！结果已保存到以下文件：")
    print(f"- GO ID: {args.go_output}")
    print(f"- KEGG ID: {args.kegg_output}")
    print(f"- COG ID: {args.cog_output}")

if __name__ == "__main__":
    main()
