import argparse
import pandas as pd
import json
import time

def load_kegg_ids(input_file):
    """加载 KEGG ID 和 counts 文件"""
    return pd.read_csv(input_file, sep="\t", header=None, names=["KEGG_ID", "Count"])


def load_local_databases(kegg_db_path, pathway_db_path):
    """加载本地数据库"""
    with open(kegg_db_path, "r") as f:
        kegg_to_pathway = json.load(f)
    with open(pathway_db_path, "r") as f:
        pathway_to_class = json.load(f)
    return kegg_to_pathway, pathway_to_class


def process_kegg_data(kegg_data, kegg_to_pathway, pathway_to_class):
    """处理 KEGG 数据，使用本地数据库"""
    pathway_info_list = []

    for _, row in kegg_data.iterrows():
        kegg_id = row["KEGG_ID"]
        count = row["Count"]

        # 从本地数据库获取 pathway_id 列表
        pathway_ids = kegg_to_pathway.get(kegg_id, [])
        for pathway_id in pathway_ids:
            # 从本地数据库获取分类信息
            pathway_info = pathway_to_class.get(pathway_id)
            if pathway_info:
                pathway_info_list.append({
                    "KEGG_ID": kegg_id,
                    "Count": count,
                    "Pathway_ID": pathway_id,
                    "Pathway_Name": pathway_info["Pathway_Name"],
                    "Level1": pathway_info["Level1"],
                    "Level2": pathway_info["Level2"]
                })

    return pd.DataFrame(pathway_info_list)


def save_results(df, output_mapping, output_counts):
    """保存结果到文件（不再输出 all_pathway_data）"""
    df.to_csv(output_mapping, index=False, sep="\t")

    # 统计二级类群数量
    level2_counts = df.groupby(["Level1", "Level2"])["Count"].sum().reset_index()
    level2_counts = level2_counts.sort_values(by=["Level1", "Level2"])
    level2_counts.to_csv(output_counts, index=False, sep="\t")


def main():
    parser = argparse.ArgumentParser(description="KEGG 数据分析脚本（离线版）")
    parser.add_argument('-i', '--input', required=True, help="输入的 KEGG ID 文件路径")
    parser.add_argument('-o1', '--output_mapping', required=True, help="输出的 KEGG ID 与 pathway 映射文件路径")
    parser.add_argument('-o2', '--output_counts', required=True, help="输出的二级类群统计文件路径")
    parser.add_argument('-d1', '--kegg_db', required=True, help="KEGG ID 到 Pathway 的本地数据库")
    parser.add_argument('-d2', '--pathway_db', required=True, help="Pathway ID 到类别的本地数据库")
    args = parser.parse_args()

    # 计时开始
    start_time = time.time()

    # 加载 KEGG ID 文件
    kegg_data = load_kegg_ids(args.input)

    # 加载本地数据库
    kegg_to_pathway, pathway_to_class = load_local_databases(args.kegg_db, args.pathway_db)

    # 处理 KEGG 数据
    pathway_info_df = process_kegg_data(kegg_data, kegg_to_pathway, pathway_to_class)

    # 保存结果（不再合并、不再输出 JSON）
    save_results(pathway_info_df, args.output_mapping, args.output_counts)

    # 输出完成信息
    print("KEGG ID 与 pathway 的映射关系已保存到", args.output_mapping)
    print("二级类群统计结果已保存到", args.output_counts)
    print(f"处理完成，耗时 {time.time() - start_time:.2f} 秒")


if __name__ == "__main__":
    main()