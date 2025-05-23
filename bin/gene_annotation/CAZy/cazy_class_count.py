import argparse
import csv
import json
import re

def load_json(file_path):
    """加载 JSON 文件并返回字典"""
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def extract_class_keys(cazy_id_value):
    """从 CAZy ID 列值中提取所有类群键"""
    # 假设格式为 "AVW75449.1|CBM35|GH5_41" 或 "QJK91651.1|GH13_16"
    parts = cazy_id_value.split('|')[1:]  # 忽略第一个部分，即蛋白质ID
    
    parts = cazy_id_value.split('|')[1:]
    
    class_keys = []
    for part in parts:
        match = re.match(r'^([A-Z]+)', part)
        if match:
            # 使用 match.group(1) 获取匹配到的分类键
            class_keys.append(match.group(1))
    return class_keys

    

def count_classes(tsv_path, cazy_data):
    """统计每个类群的出现次数"""
    class_counts = {key: 0 for key in cazy_data.keys()}
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        # 添加表头
        # headers = ["Gene ID", "CAZy ID", "% Identical", "Length", "Mismatches", "Gap Open", "Gene Start", "Gene End", "CAZy Start", "CAZy End", "E Value", "Bit Score"]
        # next(reader)  # 如果已经有表头，则跳过此行；若无表头，这行可以注释掉
        for row in reader:
            cazy_id_value = row[1]  # 第二列 (CAZy ID)
            class_keys = extract_class_keys(cazy_id_value)
            for class_key in class_keys:
                if class_key in class_counts:
                    class_counts[class_key] += 1
    return class_counts

# def generate_classify_name(cazy_info):
#     """生成 Classify 列格式：chinese_name(classify, abbr)"""
#     return f"{cazy_info['chinese_name']}({cazy_info['classify']}, {cazy_info['abbr']})"

def write_output(output_path, class_counts, cazy_data):
    """写入输出 TSV 文件"""
    with open(output_path, 'w', encoding='utf-8', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['CAZy Module','Description', 'Count'])
        for class_key in cazy_data:
            # classify_name = generate_classify_name(cazy_data[class_key])
            modules = cazy_data[class_key]["abbr"]
            description = f"{cazy_data[class_key]["chinese_name"]}({cazy_data[class_key]["classify"]})"
            count = class_counts.get(class_key, 0)
            writer.writerow([modules,description,count])

def main(tsv_path, json_path, output_path):
    """主函数：统计类群并生成输出文件"""
    # 加载 JSON 数据
    cazy_data = load_json(json_path)
    # 统计类群出现次数
    class_counts = count_classes(tsv_path, cazy_data)
    # 写入输出文件
    write_output(output_path, class_counts, cazy_data)

if __name__ == '__main__':
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Generate CAZy class count from TSV and JSON files.')
    parser.add_argument('--tsv', required=True, help='Path to input TSV file without header')
    parser.add_argument('--json', required=True, help='Path to cazy_level1.json file')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    args = parser.parse_args()

    # 执行主函数
    main(args.tsv, args.json, args.output)