#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

def extract_english_description(desc):
    """
    从描述中提取括号内的英文文本。
    
    参数:
        description: 包含中英文描述的字符串。
    
    返回:
        str: 括号内的英文文本。
    """
    
    if '(' in desc and ')' in desc:
        start = desc.find('(') + 1
        end = desc.find(')')
        return desc[start:end].strip()
    else:
        return "Unknown"


def generate_pie_chart(input_file, output_file):
    """
    从给定的输入文件读取类别名称和对应的计数，并生成空心饼图保存到输出文件。
    
    参数:
        input_file: 包含类别、描述和计数信息的TSV文件路径。
        output_file: 图表保存的路径。
    """
    # 使用pandas读取数据
    df = pd.read_csv(input_file, delimiter='\t')
    
    # 确保列名正确
    if len(df.columns) < 3 or not all(col in df.columns for col in ['CAZy Module', 'Description', 'Count']):
        raise ValueError("输入文件必须包含三列：'CAZy Module', 'Description', 'Count'")
    
    labels = df['CAZy Module'].tolist()      # 第一列用于内部楔形块
    descriptions = df['Description'].tolist()  # 第二列用于提取英文图例
    sizes = df['Count'].tolist()             # 第三列是计数
    
    # 打印原始描述
    # print("Original Descriptions:")
    # for desc in descriptions:
    #     print(desc)

    # 提取括号中的英文描述作为图例标签
    legend_labels = [extract_english_description(desc) for desc in descriptions]

    # 打印提取后的英文描述
    # print("Extracted English Descriptions:")
    # for label in legend_labels:
    #     print(label)

    # 设置绘图风格
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))

    # 计算百分比
    total = sum(sizes)
    percentages = [f"{s / total * 100:.2f}%" for s in sizes]

    # 绘制空心饼图（甜甜圈）
    wedges, texts = ax.pie(
        sizes,
        startangle=90,
        colors=sns.color_palette("Set3", n_colors=len(sizes)),
        wedgeprops=dict(width=0.4),
        textprops={'fontsize': 10}
    )

    # 添加百分比到图例
    legend_labels_with_pct = [f"{label} ({pct})" for label, pct in zip(legend_labels, percentages)]
    
    # 图例设置：显示英文描述 + 百分比
    ax.legend(wedges, legend_labels_with_pct,
              title="Categories",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=10)

    # 隐藏内部标签
    plt.setp(texts, visible=False)

    # 在中心添加一个白色圆圈形成“甜甜圈”效果
    centre_circle = plt.Circle((0, 0), 0.3, color='white')
    ax.add_artist(centre_circle)

    plt.title("Distribution of CAZy Modules", fontsize=16, pad=20)

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
