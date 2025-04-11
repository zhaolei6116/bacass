

import argparse
import pandas as pd

def load_cog_ids(cog_id_file):
    """从文件中加载 COG ID 列表及其 counts"""
    cog_ids = []
    with open(cog_id_file, 'r') as f:
        for line in f:
            if line.strip():  # 跳过空行
                parts = line.strip().split("\t")
                if len(parts) == 2:  # 确保每行有两列
                    cog_ids.append((parts[0], int(parts[1])))  # (COG ID, counts)
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

def load_cat2desc_mapping(fun_tab_file):
    """从 fun.tab 文件中加载 Category -> Description 的映射"""
    cat2desc = {}
    with open(fun_tab_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                cat2desc[parts[0]] = parts[3]
    return cat2desc

def generate_mapping_table(cog_ids, cog2cat, cat2desc):
    """生成 COG ID 到 Category 和 Description 的映射表，并保留 counts"""
    result = []
    for cog, counts in cog_ids:
        category = cog2cat.get(cog, "NA")
        description = cat2desc.get(category, "Unknown category")
        result.append({
            "COG_ID": cog,
            "Category_ID": category,
            "Category_Description": description,
            "Counts": counts
        })
    return pd.DataFrame(result)


def generate_mapping_table2(cog_ids, cog2cat, cat2desc):
    """生成 COG ID 到 Category 和 Description 的映射表，并保留 counts"""
    result = []
    for cog, counts in cog_ids:
        category = cog2cat.get(cog, "NA")
        if category == "NA":
            description = "Unknown category"
        else:
            # 如果 category 的长度大于1，拆分成单个字符
            if len(category) > 1:
                categories = list(category)
                descriptions = [cat2desc.get(cat, "Unknown category") for cat in categories]
                description = "; ".join(descriptions)
            else:
                description = cat2desc.get(category, "Unknown category")
        result.append({
            "COG_ID": cog,
            "Category_ID": category,
            "Category_Description": description,
            "Counts": counts
        })
    return pd.DataFrame(result)


def count_categories(df):
    """统计每个 Category 的数量（基于 counts 列求和）"""
    return df.groupby(["Category_ID", "Category_Description"])["Counts"].sum().reset_index(name="Total_Counts")


def count_categories2(df, cat2desc):
    """统计每个 Category 的数量（基于 counts 列求和）"""
    # 创建一个空的列表来存储拆分后的行
    expanded_rows = []
    
    # 遍历每一行
    for _, row in df.iterrows():
        category_id = row['Category_ID']
        counts = row['Counts']
        
        # 如果 Category_ID 是多个字符，拆分成单个字符
        if len(category_id) > 1:
            for char in category_id:
                # 获取每个字符对应的描述信息
                description = cat2desc.get(char, "Unknown category")
                expanded_rows.append({
                    'Category_ID': char,
                    'Category_Description': description,
                    'Counts': counts
                })
        else:
            expanded_rows.append({
                'Category_ID': category_id,
                'Category_Description': row['Category_Description'],
                'Counts': counts
            })
    
    # 将扩展后的行转换为 DataFrame
    expanded_df = pd.DataFrame(expanded_rows)
    
    # 按 Category_ID 和 Category_Description 分组，求和
    result = expanded_df.groupby(["Category_ID", "Category_Description"])["Counts"].sum().reset_index(name="Total_Counts")
    
    return result

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="从 COG ID 列表中提取 Category 和 Description 信息，并统计每个 Category 的总 counts")
    parser.add_argument('-i', '--input', required=True, help="输入的 COG ID 文件路径（cog_ids_with_counts.txt）")
    parser.add_argument('-d', '--def_tab', required=True, help="COG 定义文件路径（cog-24.def.tab）")
    parser.add_argument('-f', '--fun_tab', required=True, help="COG 功能描述文件路径（cog-24.fun.tab）")
    parser.add_argument('-o1', '--output_mapping', required=True, help="输出的 COG 映射表文件路径（cog_category_mapping.csv）")
    parser.add_argument('-o2', '--output_counts', required=True, help="输出的 Category 统计文件路径（category_counts.csv）")

    # 解析命令行参数
    args = parser.parse_args()

    # 加载 COG ID 列表及其 counts
    cog_ids = load_cog_ids(args.input)

    # 加载 COG -> Category 映射
    cog2cat = load_cog2cat_mapping(args.def_tab)

    # 加载 Category -> Description 映射
    cat2desc = load_cat2desc_mapping(args.fun_tab)

    # 生成映射表
    mapping_df = generate_mapping_table2(cog_ids, cog2cat, cat2desc)

    # 统计每个 Category 的总 counts
    # count_df = count_categories(mapping_df)
    count_df = count_categories2(mapping_df, cat2desc)

    # 保存结果
    mapping_df.to_csv(args.output_mapping, index=False, sep="\t")
    count_df.to_csv(args.output_counts, index=False, sep="\t")

    print(f"COG 映射表已保存到 {args.output_mapping}")
    print(f"Category 统计结果已保存到 {args.output_counts}")

if __name__ == "__main__":
    main()
