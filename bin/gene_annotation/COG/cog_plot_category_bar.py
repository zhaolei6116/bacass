import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns


def load_data(input_file):
    """加载 category_counts.tsv 文件"""
    return pd.read_csv(input_file, sep='\t')


def prepare_plot_data(df):
    """准备用于绘图的数据，并生成 Group 到颜色的映射"""
    # 为每个 Group 分配颜色
    group_list = df[['Group_ID', 'Group_Description']].drop_duplicates().sort_values('Group_ID')
    group2color = dict(zip(group_list['Group_ID'], sns.color_palette("husl", len(group_list))))

    # 构造 y_labels 和对应的颜色
    df['y_label'] = df['Category_ID'].apply(lambda x: f'[{x}]') + ' ' + df['Category_Description'] 

    return df, group2color, group_list


def plot_horizontal_bar(df, group2color, group_list, output_file):
    """绘制横向条形图"""
    fig, ax = plt.subplots(figsize=(10, len(df) * 0.5))

    # 绘制每个条形
    for idx, row in df.iterrows():
        bar = ax.barh(
            y=row['y_label'],
            width=row['Total_Counts'],
            color=group2color[row['Group_ID']],
            edgecolor='black'
        )

        ax.text(
            row['Total_Counts'] + 3,  # X位置稍微向右移动一点，使其不直接位于条形的边缘上
            bar[0].get_y() + bar[0].get_height() / 2.0,  # Y位置为条形的中心
            str(row['Total_Counts']),  # 显示的文本内容为'Total_Counts'的值
            ha='left', va='center'
        )

    # === 给 Y 轴标签上色 ===
    yticks = ax.get_yticklabels()
    for tick, (_, row) in zip(yticks, df.iterrows()):
        tick.set_color(group2color[row['Group_ID']])        

    # 设置图例
    legend_patches = [
        mpatches.Patch(color=group2color[grp], label=desc)
        for grp, desc in zip(group_list['Group_ID'], group_list['Group_Description'])
    ]
    ax.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1.05, 0.5)) # 调整图例位置到右侧中间
    ax.invert_yaxis()

    # 图表美化
    ax.set_xlabel("Count")
    ax.set_ylabel("Category")
    ax.set_title("Category Counts by Group")
    ax.grid(axis='x', linestyle='--', alpha=0.7)

    # 自动调整布局
    plt.tight_layout(rect=[0,0,0.85,1])
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"图表已保存至: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="基于 category_counts.tsv 生成横向柱状图")
    parser.add_argument('-i', '--input', required=True, help="输入文件路径 (category_counts.tsv)")
    parser.add_argument('-o', '--output', required=True, help="输出图像路径 (例如: category_barplot.png)")

    args = parser.parse_args()

    # 加载数据
    df = load_data(args.input)

    # 准备绘图数据
    df_plot, group2color, group_list = prepare_plot_data(df)

    # 绘图
    plot_horizontal_bar(df_plot, group2color, group_list, args.output)


if __name__ == "__main__":
    main()
