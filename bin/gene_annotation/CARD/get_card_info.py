#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import pandas as pd

def load_card(card_file):
    """加载 CARD 库的 JSON 文件"""
    with open(card_file, "r") as f:
        card_data = json.load(f)
    
    # 将 CARD 数据转换为字典，key 是 accession，value 是 name
    card_dict = {item["accession"]: item["name"] for item in card_data}
    return card_dict

def process_blast_results(blast_file, card_dict, output_file):
    """处理 BLAST 结果，结合 CARD 库生成最终结果"""
    # 加载 BLAST 结果
    blast_results = pd.read_csv(blast_file, sep="\t", header=None, names=[
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ])

    # 生成结果
    results = []
    for _, row in blast_results.iterrows():
        qseqid = row["qseqid"]
        sseqid = row["sseqid"]
        name = card_dict.get(sseqid, "Unknown")  # 查询 name
        results.append((qseqid, sseqid, name))

    # 保存结果
    result_df = pd.DataFrame(results, columns=["ProtID", "ARO ID", "Name"])
    result_df.to_csv(output_file, sep="\t", index=False)

    print(f"结果已保存到 {output_file}")

def main():
    import argparse

    # 设置命令行参数
    parser = argparse.ArgumentParser(description="结合 BLAST 结果和 CARD 库生成三列信息")
    parser.add_argument("-b", "--blast", required=True, help="BLAST 结果文件路径")
    parser.add_argument("-c", "--card", required=True, help="CARD 库文件路径（JSON 格式）")
    parser.add_argument("-o", "--output", required=True, help="输出文件路径")
    args = parser.parse_args()

    # 加载 CARD 库
    card_dict = load_card(args.card)

    # 处理 BLAST 结果
    process_blast_results(args.blast, card_dict, args.output)

if __name__ == "__main__":
    main()
