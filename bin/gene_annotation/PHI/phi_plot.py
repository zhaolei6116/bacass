#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def plot_pie_chart(input_file, output_file):
    # 读取数据
    df = pd.read_csv(input_file, sep='\t')

    # 检查必须列是否存在
    if 'Mutant Phenotype' not in df.columns or 'Count' not in df.columns:
        raise ValueError("输入文件必须包含 'Mutant Phenotype' 和 'Count' 列")

    labels = df['Mutant Phenotype'].tolist()
    sizes = df['Count'].tolist()

    # 设置绘图风格
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 8))

    # 计算百分比
    total = sum(sizes)
    percentages = [f"{(s / total * 100):.2f}%" for s in sizes]

    # 绘制饼图（不显示标签）
    wedges, texts = plt.pie(
        sizes,
        startangle=90,
        colors=sns.color_palette("Set3"),
        wedgeprops=dict(width=0.4)
    )

    # 图例设置：显示表型名 + 百分比
    legend_labels = [f"{label} ({pct})" for label, pct in zip(labels, percentages)]
    plt.legend(wedges, legend_labels,
               title="Mutant Phenotype",
               loc="center left",
               bbox_to_anchor=(1, 0, 0.5, 1),
               fontsize=10)

    plt.setp(texts, visible=False)  # 隐藏内部文本
    plt.title("Distribution of Mutant Phenotypes", fontsize=14, pad=20)

    # 保存图片
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"饼图已保存到 {output_file}")

def main():
    parser = argparse.ArgumentParser(description="根据 TSV 文件生成 Mutant Phenotype 的饼图")
    parser.add_argument("-i", "--input", required=True, help="输入文件路径（TSV 格式）")
    parser.add_argument("-o", "--output", default="pie_chart.png", help="输出图片路径（默认 pie_chart.png）")
    args = parser.parse_args()

    plot_pie_chart(args.input, args.output)

if __name__ == "__main__":
    main()
