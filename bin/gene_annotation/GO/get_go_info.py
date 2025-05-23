#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def load_terms_from_json(json_file):
    """从 JSON 文件中加载 GO 术语信息"""
    with open(json_file, 'r') as file:
        return json.load(file)

def load_go_ids_from_file(go_id_file):
    """
    使用 pandas 加载包含 GO ID 的文件（两列：GO ID 和 Count）
    """
    df = pd.read_csv(go_id_file, sep='\t', header=None, names=['GO ID', 'Count'], dtype=str)
    df.dropna(inplace=True)  # 去除空行
    # 将 Count 转换为数值类型
    df['Count'] = pd.to_numeric(df['Count'], errors='coerce')
    df.dropna(subset=['Count'], inplace=True)  # 移除转换失败的行
    return df

def lookup_go_info(go_df, terms):
    """
    根据 GO ID 列表查询对应的 name 和 namespace，并保留 Count 信息
    返回 DataFrame
    """
    go_data = []
    for _, row in go_df.iterrows():
        go_id = row['GO ID']
        count = row['Count']
        term_info = terms.get(go_id, {"name": "Not Found", "namespace": "Not Found"})
        go_data.append({
            "GO ID": go_id,
            "Name": term_info["name"],
            "Namespace": term_info["namespace"],
            "Count": count
        })
    return pd.DataFrame(go_data)

def save_results_to_file(df, output_file):
    """将查询结果保存到 TSV 文件中"""
    df.to_csv(output_file, sep='\t', index=False)
    print(f"查询完成！结果已保存到 {output_file}")

def generate_namespace_statistics(df, stats_output_file, plot_output_file="namespace_distribution.png"):
    """生成 namespace 统计结果并绘制更现代风格的柱状图"""

    # 计算每个 Namespace 的 Count 总和
    namespace_counts = df.groupby('Namespace')['Count'].sum()

    # 保存统计结果
    namespace_counts.to_csv(stats_output_file, sep='\t', header=["Total Count"])
    print(f"Namespace 统计结果已保存到 {stats_output_file}")

    # 设置 seaborn 样式
    sns.set_style("whitegrid")
    sns.set_context("talk")

    # 绘制柱状图
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        x=namespace_counts.index,
        y=namespace_counts.values,
        hue=namespace_counts.index,  # 避免警告
        palette="viridis",
        legend=False
    )

    # 显示柱子上的数值标签，并调整位置避免重叠
    for p in ax.patches:
        height = p.get_height()
        ax.annotate(f'{int(height)}',
                    (p.get_x() + p.get_width() / 2., height),
                    ha='center', va='bottom', fontsize=10, color='black')

    # 调整Y轴范围，给顶部留出空间以容纳数值标签
    max_count = namespace_counts.max()
    ax.set_ylim(top=max_count + max_count * 0.1)  # 增加10%的空间


    # ax.bar_label(ax.containers[0])  # 显示数字标签
    plt.title('Sum of Counts per GO Term Namespace', fontsize=16)
    plt.xlabel('Namespace', fontsize=14)
    plt.ylabel('Total Count', fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.tight_layout()

    # 保存图像
    plt.savefig(plot_output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Namespace 分布柱状图已保存为 {plot_output_file}")

def main():
    parser = argparse.ArgumentParser(description="查询 GO ID 对应的 name 和 namespace")
    parser.add_argument('-i', '--input', required=True, help="输入的 JSON 文件路径（terms.json）")
    parser.add_argument('-o', '--output', required=True, help="输出的结果文件路径（TSV 格式）")
    parser.add_argument('-g', '--go_ids_file', required=True, help="包含 GO ID 列表的文件路径（每行一个 ID 和 Count，TSV）")
    parser.add_argument('-s', '--stats_output', required=True, help="namespace 统计结果文件路径（TSV）")
    parser.add_argument('--plot_output', default="namespace_distribution.png", help="柱状图输出文件名（默认：namespace_distribution.png）")

    args = parser.parse_args()

    # 加载 JSON 数据
    terms = load_terms_from_json(args.input)

    # 加载 GO ID 数据
    go_df = load_go_ids_from_file(args.go_ids_file)

    # 查询 GO 信息
    results_df = lookup_go_info(go_df, terms)

    # 保存完整查询结果
    save_results_to_file(results_df, args.output)

    # 生成 namespace 统计及图表
    generate_namespace_statistics(results_df, args.stats_output, args.plot_output)

if __name__ == "__main__":
    main()