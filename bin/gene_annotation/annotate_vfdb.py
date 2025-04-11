#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
import json

def load_json(file_path):
    """加载 JSON 文件并返回字典"""
    with open(file_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def parse_blast(blast_path):
    """解析 BLAST 文件，逐行返回数据"""
    with open(blast_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        # next(reader)  # 跳过表头
        for row in reader:
            yield row

def get_vf_info(setb_data, vfs_data, sseqid):
    """根据 sseqid 获取 VF 信息"""
    setb_info = setb_data.get(sseqid)
    if not setb_info:
        return None

    vf_id = setb_info.get('VF_ID')
    if vf_id in vfs_data:
        vf_info = vfs_data[vf_id]
        return {
            'VFID': vf_id,
            'VF_Name': vf_info.get('VF_Name', ''),
            'Bacteria': vf_info.get('Bacteria', ''),
            'VFCID': vf_info.get('VFCID', ''),
            'VFcategory': vf_info.get('VFcategory', ''),
            'Function': vf_info.get('Function', ''),
            'Mechanism': vf_info.get('Mechanism', '')
        }
    else:
        # 从 setB.json 获取替代信息
        return {
            'VFID': vf_id,
            'VF_Name': setb_info.get('VF_name', ''),
            'Bacteria': setb_info.get('Taxonomy', ''),
            'VFCID': setb_info.get('VF_category_id', ''),
            'VFcategory': setb_info.get('VF_category_level1', ''),
            'Function': '',
            'Mechanism': ''
        }

def main(blast_path, setb_path, vfs_path, output_path):
    """主函数：关联数据并生成 TSV 文件"""
    # 加载 JSON 数据
    setb_data = load_json(setb_path)
    vfs_data = load_json(vfs_path)

    # 写入 TSV 文件
    with open(output_path, 'w', encoding='utf-8', newline='') as out_file:
        writer = csv.writer(out_file, delimiter='\t')
        writer.writerow(['GeneID', 'VF_gene_id', 'VFID', 'VF_Name', 'Bacteria', 'VFCID', 'VFcategory', 'Function', 'Mechanism'])

        # 解析 BLAST 文件并关联信息
        for row in parse_blast(blast_path):
            qseqid, sseqid = row[0], row[1]
            vf_info = get_vf_info(setb_data, vfs_data, sseqid)
            if vf_info:
                writer.writerow([
                    qseqid,
                    sseqid,
                    vf_info['VFID'],
                    vf_info['VF_Name'],
                    vf_info['Bacteria'],
                    vf_info['VFCID'],
                    vf_info['VFcategory'],
                    vf_info['Function'],
                    vf_info['Mechanism']
                ])

if __name__ == '__main__':
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Annotate BLAST results with VFDB database.')
    parser.add_argument('--blast', required=True, help='Path to BLAST output file (blast.out)')
    parser.add_argument('--setb', required=True, help='Path to setB.json file')
    parser.add_argument('--vfs', required=True, help='Path to VFs.json file')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    args = parser.parse_args()

    # 执行主函数
    main(args.blast, args.setb, args.vfs, args.output)