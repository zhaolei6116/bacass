import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


def load_data(input_file):
    """加载 Level1/Level2/Count 格式的 TSV 文件"""
    df = pd.read_csv(input_file, sep='\t')
    
    # 确保列名正确（去除可能的空格）
    df.columns = df.columns.str.strip()
    
    return df


def prepare_plot_data(df):
    """
    准备用于绘图的数据，并生成 Level1 到颜色的映射
    """
    # 获取所有唯一的一级分类并排序
    unique_level1 = df['Level1'].unique()
    palette = sns.color_palette("husl", len(unique_level1))
    level1_to_color = dict(zip(unique_level1, palette))

    # 构造 y_label 为 Level2 名称
    df['y_label'] = df['Level2']

    return df, level1_to_color, unique_level1


def plot_horizontal_bar(df, level1_to_color, output_file):
    """绘制横向条形图"""

    fig, ax = plt.subplots(figsize=(15, len(df) * 0.45))  # 动态设置高度

    # 获取当前轴的最大宽度值
    max_count = df['Count'].max()

    label_buffer = 50
    plt.xlim(0, max_count + label_buffer)

    # 绘制每个条形
    for idx, row in df.iterrows():
        bar = ax.barh(
            y=row['y_label'],
            width=row['Count'],
            color=level1_to_color[row['Level1']],
            edgecolor='black'
        )

        # 根据最大宽度动态调整文本位置，确保不超出边界
        # text_x_position = min(row['Count'] + max(1, row['Count'] * 0.02), max_count * 1.1)

        # 在柱子右侧添加数值标签
        ax.text(
            x=row['Count']+3,
            # row['Count'] + max(1, row['Count'] * 0.02),  # 右边留点空间
            y=bar[0].get_y() + bar[0].get_height() / 2,
            s=str(row['Count']),
            va='center',
            ha='left',
            fontsize=9,
            color='black'
        )
    # ax.spines['right'].set_visible(False)  # 去掉右边框线
    # === 给 Y 轴标签上色 ===
    yticks = ax.get_yticklabels()
    for tick, (_, row) in zip(yticks, df.iterrows()):
        tick.set_color(level1_to_color[row['Level1']])

    # 设置图例
    legend_patches = [
        mpatches.Patch(color=color, label=level1)
        for level1, color in level1_to_color.items()
    ]
    ax.legend(handles=legend_patches, title="Level1", loc='center left', bbox_to_anchor=(1.05, 0.5))

    # 坐标轴反转，使第一个类别在最上方
    ax.invert_yaxis()

    # 设置图表标题和标签
    ax.set_xlabel('Count')
    ax.set_ylabel('Functional Categories (Level2)')
    ax.set_title('KEGG Functional Classification by Level2')
    ax.grid(axis='x', linestyle='--', alpha=0.7)

    # 自动调整布局
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"图表已保存至: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="基于 Level1/Level2/Count 格式的 TSV 文件生成横向柱状图")
    parser.add_argument('-i', '--input', required=True, help="输入文件路径 (TSV 格式)")
    parser.add_argument('-o', '--output', required=True, help="输出图像路径 (例如: kegg_level2_barplot.png)")

    args = parser.parse_args()

    # 加载数据
    df = load_data(args.input)

    # 准备绘图数据
    df_plot, level1_to_color, _ = prepare_plot_data(df)

    # 绘图
    plot_horizontal_bar(df_plot, level1_to_color, args.output)


if __name__ == "__main__":
    main()