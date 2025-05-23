#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import json
import argparse
from collections import defaultdict, Counter

def load_phi_json(phi_json_file):
    """加载预处理好的 PHI ID 到表型列表的映射"""
    with open(phi_json_file, 'r') as f:
        phi_data = json.load(f)
    return phi_data

def process_blast_results(blast_file, phi_data, output_file1, output_file2, output_file3, sep="\t"):
    """处理 BLAST 结果，结合 PHI 表型信息生成最终结果"""

    # 动态识别 BLAST 结果的分隔符（支持 CSV 或 TSV）
    blast_results = pd.read_csv(
        blast_file,
        sep=sep,
        header=None,
        names=[
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
    )

    # 确保 sseqid 是字符串类型
    blast_results["sseqid"] = blast_results["sseqid"].astype(str)

    # 存储输出文件1的数据
    file1_data = []
    phenotype_counter = Counter()  # 统计每种表型总出现次数（用于文件2）
    phi_id_counter = Counter()     # 统计每个 PHI ID 的出现次数（用于文件3）
    
    for _, row in blast_results.iterrows():
        gene_id = row['qseqid']
        phi_id = row['sseqid']

        phenotypes = phi_data.get(phi_id, [])  # 获取对应的所有表型

        # 如果没有匹配到 PHI 表型，至少也要保留这一行（可选）
        if not phenotypes:
            file1_data.append({
                "qseqid": gene_id,
                "PHI_ID": phi_id,
                "Mutant Phenotype": ""
            })
        else:
            # 收集当前 PHI ID 所有表型用于统计
            for p in phenotypes:
                phenotype_counter[p] += 1
            
            # 合并表型为字符串
            phenotype_str = ";".join(sorted(set(phenotypes)))  # 去重 + 排序
            file1_data.append({
                "qseqid": gene_id,
                "PHI_ID": phi_id,
                "Mutant Phenotype": phenotype_str
            })

            # 更新 PHI ID 计数器
            phi_id_counter[phi_id] += 1

    # 写入文件1
    file1_df = pd.DataFrame(file1_data)
    file1_df.to_csv(output_file1, sep=sep, index=False)

    # 写入文件2：Mutant Phenotype 统计
    stats_data = [{"Mutant Phenotype": k, "Count": v} for k, v in phenotype_counter.items()]
    stats_df = pd.DataFrame(stats_data)
    stats_df.sort_values(by="Count", ascending=False, inplace=True)
    stats_df.to_csv(output_file2, sep=sep, index=False)

    # 写入文件3：每个 PHI ID 及其对应的 Mutant Phenotype 和出现次数
    file3_data = []
    for phi_id, count in phi_id_counter.items():
        phenotypes = phi_data.get(phi_id, [])
        phenotype_str = ";".join(sorted(set(phenotypes)))  # 去重 + 排序
        file3_data.append({
            "PHI_ID": phi_id,
            "Mutant Phenotype": phenotype_str,
            "Counts": count
        })

    file3_df = pd.DataFrame(file3_data)
    file3_df.sort_values(by="PHI_ID", inplace=True)
    file3_df.to_csv(output_file3, sep=sep, index=False)

    print(f"文件1已保存到 {output_file1}")
    print(f"文件2已保存到 {output_file2}")
    print(f"文件3已保存到 {output_file3}")

def main():
    parser = argparse.ArgumentParser(description="结合 BLAST 结果和 PHI 注释生成注释结果")
    parser.add_argument(
        "-b", "--blast", required=True, help="BLAST 结果文件路径（支持 CSV 或 TSV）"
    )
    parser.add_argument(
        "-p", "--phi_json", required=True, help="预处理好的 PHI 映射 JSON 文件路径"
    )
    parser.add_argument(
        "-o1", "--output1", required=True, help="输出文件1路径（qseqid, PHI_ID, Mutant Phenotype）"
    )
    parser.add_argument(
        "-o2", "--output2", required=True, help="输出文件2路径（Mutant Phenotype 统计）"
    )
    parser.add_argument(
        "-o3", "--output3", required=True, help="输出文件3路径（每个 PHI ID 的 Mutant Phenotype 及其出现次数）"
    )
    parser.add_argument(
        "-s", "--separator", default="\t", help="输入/输出文件分隔符（默认制表符）"
    )
    args = parser.parse_args()

    phi_data = load_phi_json(args.phi_json)
    process_blast_results(
        args.blast,
        phi_data,
        args.output1,
        args.output2,
        args.output3,
        sep=args.separator
    )

if __name__ == "__main__":
    main()
