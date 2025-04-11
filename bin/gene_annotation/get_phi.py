#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse

def load_phi_base(phi_base_file):
    """加载并处理 phi-base_current.csv 文件"""
    # 跳过前两行，使用第二行作为表头
    phi_base_df = pd.read_csv(phi_base_file, skiprows=1, dtype=str, na_values=None)
    
    # 删除前两列（可能是空列或无关列）
    phi_base_df = phi_base_df.iloc[:, 2:]
    
    return phi_base_df

def process_blast_results(blast_file, phi_base_df, output_file1, output_file2, sep="\t"):
    """处理 BLAST 结果，结合 phi-base 库生成最终结果"""
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
    
    # 确保 sseqid 和 PHI_MolConn_ID 是字符串类型以避免合并错误
    blast_results["sseqid"] = blast_results["sseqid"].astype(str)
    phi_base_df["PHI_MolConn_ID"] = phi_base_df["PHI_MolConn_ID"].astype(str)
    
    # 左连接合并数据
    merged_df = pd.merge(
        blast_results,
        phi_base_df,
        left_on="sseqid",
        right_on="PHI_MolConn_ID",
        how="left"
    )
    
    merged_df.rename(columns={
        "sseqid": "PHI_ID",      # sseqid → PHI_ID
        "qseqid": "GeneID"       # qseqid → GeneID
    }, inplace=True)
    
    
    # 生成文件1：包含 qseqid, PHI_MolConn_ID 及 phi-base 所有列
    file1_columns = ["GeneID", "PHI_ID"] + [col for col in phi_base_df if col != "PHI_MolConn_ID"]
    file1_df = merged_df[file1_columns]
    file1_df.to_csv(output_file1, sep=sep, index=False)
    
    # 生成文件2：包含 qseqid, PHI_ID, Mutant Phenotype
    # 检查必要列是否存在（避免 KeyError）
    required_columns = ["GeneID", "PHI_ID", "Mutant Phenotype"]
    file2_df = merged_df[required_columns]
    file2_df.to_csv(output_file2, sep=sep, index=False)

    print(f"文件1已保存到 {output_file1}")
    print(f"文件2已保存到 {output_file2}")

def main():
    parser = argparse.ArgumentParser(description="结合 BLAST 结果和 phi-base 库生成注释结果")
    parser.add_argument(
        "-b", "--blast", required=True, help="BLAST 结果文件路径（支持 CSV 或 TSV）"
    )
    parser.add_argument(
        "-p", "--phi_base", required=True, help="phi-base_current.csv 文件路径"
    )
    parser.add_argument(
        "-o1", "--output1", required=True, help="输出文件1路径（保留所有注释信息）"
    )
    parser.add_argument(
        "-o2", "--output2", required=True, help="输出文件2路径（仅 PHI_ID 和表型）"
    )
    parser.add_argument(
        "-s", "--separator", default="\t", help="输入/输出文件分隔符（默认制表符）"
    )
    args = parser.parse_args()

    phi_base_df = load_phi_base(args.phi_base)
    process_blast_results(
        args.blast,
        phi_base_df,
        args.output1,
        args.output2,
        sep=args.separator
    )

if __name__ == "__main__":
    main()
