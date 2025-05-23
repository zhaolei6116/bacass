import argparse
import pandas as pd
import json

def load_cog_ids(cog_id_file):
    """从文件中加载 COG ID 列表及其 counts"""
    cog_ids = []
    with open(cog_id_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) == 2:
                    cog_ids.append((parts[0], int(parts[1])))
    return cog_ids

def load_cog2cat_mapping(def_tab_file):
    """从 def.tab 文件中加载 COG -> Category 的映射"""
    cog2cat = {}
    with open(def_tab_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                cog2cat[parts[0]] = parts[1]
    return cog2cat

def load_fun_json(fun_json_file):
    """加载 fun.json 文件，返回一个列表，其中每个元素是一个完整的 category 信息字典"""
    with open(fun_json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    return data

def generate_mapping_table(cog_ids, cog2cat, fun_json_data):
    """生成 COG ID 到 Category 和 Description 的映射表，并保留 counts"""
    result = []
    # 构建 Category_ID -> {Group_ID, Group_Description, Category_Description} 映射
    cat_info_map = {
        item["Category_ID"]: {
            "Group_ID": item["Group_ID"],
            "Group_Description": item["Group_Description"],
            "Category_Description": item["Category_Description"]
        }
        for item in fun_json_data
    }

    for cog, counts in cog_ids:
        category = cog2cat.get(cog, "NA")
        if category == "NA":
            result.append({
                "COG_ID": cog,
                "Group_ID": "Unknown",
                "Group_Description": "Unknown group",
                "Category_ID": "NA",
                "Category_Description": "Unknown category",
                "Counts": counts
            })
        else:
            # 如果 category 长度 > 1，拆分处理
            if len(category) > 1:
                for c in list(category):
                    info = cat_info_map.get(c, {
                        "Group_ID": "Unknown",
                        "Group_Description": "Unknown group",
                        "Category_Description": "Unknown category"
                    })
                    result.append({
                        "COG_ID": cog,
                        "Group_ID": info["Group_ID"],
                        "Group_Description": info["Group_Description"],
                        "Category_ID": c,
                        "Category_Description": info["Category_Description"],
                        "Counts": counts
                    })
            else:
                info = cat_info_map.get(category, {
                    "Group_ID": "Unknown",
                    "Group_Description": "Unknown group",
                    "Category_Description": "Unknown category"
                })
                result.append({
                    "COG_ID": cog,
                    "Group_ID": info["Group_ID"],
                    "Group_Description": info["Group_Description"],
                    "Category_ID": category,
                    "Category_Description": info["Category_Description"],
                    "Counts": counts
                })

    return pd.DataFrame(result)

def count_categories_by_json_order(df, fun_json_data):
    """
    按照 JSON 中 Category 出现的顺序统计每个 Category 的总 counts。
    返回 DataFrame 包含：
        Group_ID, Group_Description,
        Category_ID, Category_Description,
        Total_Counts
    """
    result = []

    # 先统计每个 Category 的总 counts
    category_counts = df.groupby("Category_ID")["Counts"].sum().to_dict()

    # 然后按照 JSON 中的顺序遍历每个 Category
    for item in fun_json_data:
        cat_id = item["Category_ID"]
        total = category_counts.get(cat_id, 0)
        result.append({
            "Group_ID": item["Group_ID"],
            "Group_Description": item["Group_Description"],
            "Category_ID": cat_id,
            "Category_Description": item["Category_Description"],
            "Total_Counts": total
        })

    return pd.DataFrame(result)

def main():
    parser = argparse.ArgumentParser(description="从 COG ID 列表中提取 Category 和 Description 信息，并统计每个 Category 的总 counts")
    parser.add_argument('-i', '--input', required=True, help="输入的 COG ID 文件路径（含 counts）")
    parser.add_argument('-d', '--def_tab', required=True, help="COG 定义文件路径（cog.def.tab）")
    parser.add_argument('-f', '--fun_tab', required=True, help="COG 功能描述文件路径（JSON 格式）")
    parser.add_argument('-o1', '--output_mapping', required=True, help="输出的 COG 映射表文件路径（TSV 格式）")
    parser.add_argument('-o2', '--output_counts', required=True, help="输出的 Category 统计文件路径（TSV 格式）")

    args = parser.parse_args()

    # 加载数据
    cog_ids = load_cog_ids(args.input)
    cog2cat = load_cog2cat_mapping(args.def_tab)
    fun_json_data = load_fun_json(args.fun_tab)

    # 生成 mapping 表
    mapping_df = generate_mapping_table(cog_ids, cog2cat, fun_json_data)

    # 按 JSON 顺序统计 Category 数量
    count_df = count_categories_by_json_order(mapping_df, fun_json_data)

    # 保存结果
    mapping_df.to_csv(args.output_mapping, index=False, sep='\t')
    count_df.to_csv(args.output_counts, index=False, sep='\t')

    print(f"COG 映射表已保存到 {args.output_mapping}")
    print(f"Category 统计结果已保存到 {args.output_counts}")

if __name__ == "__main__":
    main()
