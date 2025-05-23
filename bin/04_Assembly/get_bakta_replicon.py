#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse

def parse_arguments():
    """
    使用argparse模块解析命令行参数
    """
    parser = argparse.ArgumentParser(description="Process contig information and generate output file.")
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    args = parser.parse_args()
    return args

def read_input_file(input_file):
    """
    读取输入文件并返回DataFrame
    """
    df = pd.read_csv(input_file, sep='\t')
    return df

def process_contig_data(df):
    """
    处理contig数据并生成输出DataFrame
    """
    output_data = {
        'original sequence id': df['Contig'],
        'new sequence id': df['Contig'],
        'type': '-',
        'topology': df['Circular'].apply(lambda x: 'circular' if x == 'Y' else 'linear'),
        'name': '-'
    }
    output_df = pd.DataFrame(output_data)
    return output_df

def write_output_file(output_df, output_file):
    """
    将处理后的数据写入输出文件
    """
    output_df.to_csv(output_file, sep='\t', index=False)

def main():
    # 解析命令行参数
    args = parse_arguments()

    # 读取输入文件
    df = read_input_file(args.input)

    # 处理contig数据
    output_df = process_contig_data(df)

    # 写入输出文件
    write_output_file(output_df, args.output)

if __name__ == "__main__":
    main()