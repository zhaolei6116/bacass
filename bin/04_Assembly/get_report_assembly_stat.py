#!/usr/bin/env python3

import argparse

def read_file_to_dict(filename):
    """
    读取文件，将每行拆分存到列表，然后生成一个字典。
    
    参数:
    filename -- 文件路径
    
    返回:
    一个字典，其中键和值分别来自文件的两列。
    """
    
    # 读取文件并拆分每行到两个列表
    with open(filename, 'r', encoding='utf-8') as fileHandler:
        lines = fileHandler.readlines()
        list1 = lines[0].strip().split('\t')
        list2 = lines[1].strip().split('\t')
        
    result_dict = dict(zip(list1, list2))
    
    return result_dict



def main():
    parser = argparse.ArgumentParser(description='stat data file from Nanostat ')
    parser.add_argument('filename', type=str, help='The path to the TSV file.')
    parser.add_argument('outfile', type=str, help='The path to the out TSV file.')
    args = parser.parse_args()

    result_dict = read_file_to_dict(args.filename)

    head_line = ["Genome Size", "Contig Number", "N50", "N50n", "Longest Contig", "Shortest Contig", "Mean Length"]
    key_line = ["total_length", "number", "N50", "N50n", "longest", "shortest", "mean_length"]
    
    # 修改输出格式
    result_dict["total_length"] = f'{int(result_dict["total_length"]):,}' 
    result_dict["N50"] = f'{int(result_dict["N50"]):,}'
    result_dict["longest"] = f'{int(result_dict["longest"]):,}'
    result_dict["shortest"] = f'{int(result_dict["shortest"]):,}'
    result_dict["mean_length"] = f'{float(result_dict["mean_length"]):,}'

    value_line = [result_dict[k] for k in key_line]



    with open(args.outfile, "w") as outF:
        for line in [head_line, value_line]:
            outF.write("\t".join(line)+"\n")


if __name__ == "__main__":
    main()
