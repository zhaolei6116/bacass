#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json

def load_terms_from_json(json_file):
    """从 JSON 文件中加载 GO 术语信息"""
    with open(json_file, 'r') as file:
        return json.load(file)

def load_go_ids_from_file(go_id_file):
    """从文件中加载 GO ID 列表（每行一个 ID 和 Number）"""
    go_ids = []
    with open(go_id_file, 'r') as file:
        for line in file:
            if line.strip():  # 跳过空行
                parts = line.strip().split("\t")
                if len(parts) == 2:  # 确保每行有两列
                    go_ids.append(parts)
    return go_ids

def lookup_go_info(go_ids, terms):
    """根据 GO ID 列表查询对应的 name 和 namespace，并保留 Number 信息"""
    results = {}
    for go_id, number in go_ids:
        # 查询 GO ID 对应的 name 和 namespace
        term_info = terms.get(go_id, {"name": "Not Found", "namespace": "Not Found"})
        # 添加 Number 信息
        results[go_id] = {
            "name": term_info["name"],
            "namespace": term_info["namespace"],
            "Number": number
        }
    return results

def save_results_to_file(results, output_file):
    """将查询结果保存到文件中"""
    with open(output_file, 'w') as file:
        # 写入表头
        file.write("GO ID\tName\tNamespace\tNumber\n")
        # 写入每个 GO ID 的信息
        for go_id, info in results.items():
            file.write(f"{go_id}\t{info['name']}\t{info['namespace']}\t{info['Number']}\n")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="查询 GO ID 对应的 name 和 namespace")
    parser.add_argument('-i', '--input', required=True, help="输入的 JSON 文件路径（terms.json）")
    parser.add_argument('-o', '--output', required=True, help="输出的结果文件路径")
    parser.add_argument('-g', '--go_ids_file', required=True, help="包含 GO ID 列表的文件路径（每行一个 ID 和 Number）")

    # 解析命令行参数
    args = parser.parse_args()

    # 加载 JSON 数据
    terms = load_terms_from_json(args.input)

    # 从文件中加载 GO ID 列表
    go_ids = load_go_ids_from_file(args.go_ids_file)

    # 查询 GO ID 对应的 name 和 namespace
    results = lookup_go_info(go_ids, terms)

    # 保存查询结果到文件
    save_results_to_file(results, args.output)

    print(f"查询完成！结果已保存到 {args.output}")

if __name__ == "__main__":
    main()
