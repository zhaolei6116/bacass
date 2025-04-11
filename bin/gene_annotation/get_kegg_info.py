#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import json
from bioservices import KEGG
import time

# 初始化 KEGG 服务
kegg = KEGG()

def load_kegg_ids(input_file):
    """加载 KEGG ID 和 counts 文件"""
    return pd.read_csv(input_file, sep="\t", header=None, names=["KEGG_ID", "Count"])

def fetch_kegg_data(kegg_id, max_retries=3):
    """获取 KEGG ID 的 pathway 信息，支持重试机制"""
    retries = 0
    while retries < max_retries:
        try:
            pathway_info = kegg.get(kegg_id)
            if pathway_info:
                return kegg.parse(pathway_info)
            else:
                return None
        except Exception as e:
            print(f"Error fetching data for {kegg_id}: {e}. Retrying ({retries + 1}/{max_retries})...")
            retries += 1
            time.sleep(2)  # 等待 2 秒后重试
    print(f"Failed to fetch data for {kegg_id} after {max_retries} retries.")
    return None

def get_pathway_class(pathway_id, max_retries=3):
    """获取 pathway ID 的一级类群和二级类群信息，支持重试机制"""
    retries = 0
    while retries < max_retries:
        try:
            pathway_detail = kegg.get(pathway_id)
            if pathway_detail:
                pathway_detail_data = kegg.parse(pathway_detail)
                class_info = pathway_detail_data.get("CLASS", "")
                if class_info:
                    class_parts = class_info.split(";")
                    level1 = class_parts[0].strip() if len(class_parts) > 0 else "Unknown"
                    level2 = class_parts[1].strip() if len(class_parts) > 1 else "Unknown"
                    return level1, level2
                else:
                    return "Unknown", "Unknown"
            else:
                return "Unknown", "Unknown"
        except Exception as e:
            print(f"Error fetching pathway detail for {pathway_id}: {e}. Retrying ({retries + 1}/{max_retries})...")
            retries += 1
            time.sleep(2)  # 等待 2 秒后重试
    print(f"Failed to fetch pathway detail for {pathway_id} after {max_retries} retries.")
    return "Unknown", "Unknown"

def process_kegg_data(kegg_data):
    """处理 KEGG 数据，获取 pathway 信息和类群信息"""
    pathway_info_list = []
    all_pathway_data = {}

    for kegg_id in kegg_data["KEGG_ID"]:
        pathway_data = fetch_kegg_data(kegg_id)
        if pathway_data:
            all_pathway_data[kegg_id] = pathway_data
            pathway_ids = pathway_data.get("PATHWAY", {})
            for pathway_id, pathway_name in pathway_ids.items():
                level1, level2 = get_pathway_class(pathway_id)
                pathway_info_list.append({
                    "KEGG_ID": kegg_id,
                    "Pathway_ID": pathway_id,
                    "Pathway_Name": pathway_name,
                    "Level1": level1,
                    "Level2": level2
                })
        else:
            print(f"No pathway information found for {kegg_id}")

    return pathway_info_list, all_pathway_data

def save_results(kegg_pathway_mapping, level2_counts, pathway_data, output_mapping, output_counts, output_json):
    """保存结果到文件"""
    kegg_pathway_mapping.to_csv(output_mapping, index=False, sep="\t")
    level2_counts.to_csv(output_counts, index=True, sep="\t")
    with open(output_json, "w") as f:
        json.dump(pathway_data, f, indent=4)

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="KEGG 数据分析脚本")
    parser.add_argument('-i', '--input', required=True, help="输入的 KEGG ID 文件路径")
    parser.add_argument('-o1', '--output_mapping', required=True, help="输出的 KEGG ID 与 pathway 映射文件路径")
    parser.add_argument('-o2', '--output_counts', required=True, help="输出的二级类群统计文件路径")
    parser.add_argument('-o3', '--output_json', required=True, help="输出的 pathway_data JSON 文件路径")
    args = parser.parse_args()

    # 加载 KEGG ID 文件
    kegg_data = load_kegg_ids(args.input)

    # 处理 KEGG 数据
    pathway_info_list, all_pathway_data = process_kegg_data(kegg_data)

    # 将 pathway 信息转换为 DataFrame
    pathway_info_df = pd.DataFrame(pathway_info_list)

    # 合并 pathway 信息与 counts 信息
    kegg_pathway_mapping = pd.merge(
        kegg_data,
        pathway_info_df,
        on="KEGG_ID",
        how="left"
    )

    # 统计二级类群的数量
    level2_counts = kegg_pathway_mapping.groupby(["Level1", "Level2"]).size().reset_index(name="Counts")

    # 按 Level1 和 Level2 排序
    level2_counts = level2_counts.sort_values(by=["Level1", "Level2"])

    # 保存结果
    save_results(kegg_pathway_mapping, level2_counts, all_pathway_data, args.output_mapping, args.output_counts, args.output_json)

    print("KEGG ID 与 pathway 的映射关系已保存到", args.output_mapping)
    print("二级类群统计结果已保存到", args.output_counts)
    print("pathway_data 信息已保存到", args.output_json)

if __name__ == "__main__":
    main()
