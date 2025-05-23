#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def generate_pie_chart(input_file, output_file):
    """
    从给定的输入文件读取物种名称和对应的计数，并生成饼图保存到输出文件。
    
    参数:
        input_file: 包含物种和计数信息的TSV文件路径。
        output_file: 图表保存的路径。
    """
    # 使用pandas读取数据
    df = pd.read_csv(input_file, delimiter='\t')
    
    # 检查列数量是否为2
    if len(df.columns) < 2:
        raise ValueError("输入文件必须至少包含两列（类别和计数）")

    labels = df.iloc[:,0].tolist()  # 第一列是分类名称
    sizes = df.iloc[:,1].tolist()   # 第二列是对应的计数

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

    # 图例设置：显示分类名 + 百分比
    legend_labels = [f"{label} ({pct})" for label, pct in zip(labels, percentages)]
    plt.legend(wedges, legend_labels,
               title="Species",
               loc="center left",
               bbox_to_anchor=(1, 0, 0.5, 1),
               fontsize=10)

    plt.setp(texts, visible=False)  # 隐藏内部文本
    plt.title("Distribution of Species", fontsize=14, pad=20)

    # 保存图片
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"饼图已保存到 {output_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="根据 TSV 文件生成饼图")
    parser.add_argument('-i', '--input', type=str, required=True, help='输入 TSV 文件路径')
    parser.add_argument('-o', '--output', type=str, required=True, help='输出图片文件路径')

    args = parser.parse_args()

    generate_pie_chart(args.input, args.output)
