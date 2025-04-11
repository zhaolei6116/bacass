import argparse
import csv
import json
import re

def load_json(file_path):
    """加载 JSON 文件并返回字典"""
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def extract_class_key(diamond_value):
    """从 DIAMOND 列值中提取类群键（如 GT, GH）"""
    match = re.match(r'^([A-Za-z]+)', diamond_value)
    return match.group(1) if match else None

def count_classes(tsv_path, cazy_data):
    """统计每个类群的出现次数"""
    class_counts = {key: 0 for key in cazy_data.keys()}
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)  # 跳过表头
        for row in reader:
            diamond_value = row[4]  # 第 5 列 (DIAMOND)
            class_key = extract_class_key(diamond_value)
            if class_key in class_counts:
                class_counts[class_key] += 1
    return class_counts

def generate_classify_name(cazy_info):
    """生成 Classify 列格式：chinese_name(classify, abbr)"""
    return f"{cazy_info['chinese_name']}({cazy_info['classify']}, {cazy_info['abbr']})"

def write_output(output_path, class_counts, cazy_data):
    """写入输出 TSV 文件"""
    with open(output_path, 'w', encoding='utf-8', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['Classify', 'Count'])
        for class_key in cazy_data:
            classify_name = generate_classify_name(cazy_data[class_key])
            count = class_counts[class_key]
            writer.writerow([classify_name, count])

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
    parser.add_argument('--tsv', required=True, help='Path to out1.tsv file')
    parser.add_argument('--json', required=True, help='Path to cazy_level1.json file')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    args = parser.parse_args()

    # 执行主函数
    main(args.tsv, args.json, args.output)
