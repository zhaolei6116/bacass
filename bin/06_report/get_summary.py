#!/usr/bin/env python3

import csv
import argparse

__version__ = "v0.0.1"
__date__    = "2024.12.10"
__author__  = "zhaolei"


def get_info_dict(infile):
    '''
    获取项目信息，ID， 物种名称， 测序平台 信息字典
    '''
    with open(infile, 'r', encoding='utf-8') as inF:
       reader = csv.reader(inF, delimiter='\t')
       rows = list(reader)
       out_dict1 = dict(zip(rows[0], rows[1]))

    return out_dict1


def get_seq_info():
    '''
    获取测序信息
    '''

    return None

def get_ass_info():
    '''
    获取组装信息
    '''

    return None
  
def get_result(pjid, spname, dict_seq_info, dict_assembly_info, file_o):
    '''
    获取输出结果
    '''
    out_dict = {}
    out_list = ["项目编号", "物种名称", "测序平台", "测序数据量(bp)", "Reads条数", "最长Reads(bp)", "基因组大小(bp)", "contig 个数", "原始测序深度(X)"]

    out_dict["项目编号"] = pjid
    out_dict["物种名称"] = spname if spname != "_" else "bacteria"
    out_dict["测序平台"] = "ONT"
    out_dict["测序数据量(bp)"] = dict_seq_info["Number of base"]
    out_dict["Reads条数"] = dict_seq_info["Number of Reads"]
    out_dict["最长Reads(bp)"] = dict_seq_info["Longest Reads(Q)"].split(" ")[0]
    out_dict["基因组大小(bp)"] = dict_assembly_info["Genome Size"]
    out_dict["contig 个数"] = dict_assembly_info["Contig Number"]
    out_dict["原始测序深度(X)"] = f'{float(dict_seq_info["Number of base"].replace("," , ""))/float(dict_assembly_info["Genome Size"].replace("," , "")):.2f}' 

    with open(file_o, "w") as outF:
        for info in out_list:
            outF.write(info + "\t" + out_dict[info] + "\n")


  

def setup_args():
    '''
    '''
    parser = argparse.ArgumentParser(description = 'merge raw and clean state info into qc table',
            formatter_class=argparse.RawTextHelpFormatter,
            epilog = f"Version: {__version__}\nDate: {__date__}\nAuthor: {__author__}")
    parser.add_argument('-v','--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-p','--pj_id', type=str, help='The project id ',  required = True)
    parser.add_argument('-sp','--species_name', type=str, help='The species name',  required = True)
    parser.add_argument('-s','--seq_info', type=str, help='The path to seq info', required = True)
    parser.add_argument('-a','--assembly_info',type=str, help='The path to assembly info', required = True)
    parser.add_argument('-o','--outfile', type=str, help='The path to output file', default="pj_summary.xls")
    args = parser.parse_args()

    return args

def main():
    args = setup_args()
    # 获取输入
    pjid = args.pj_id
    spname = args.species_name
    file_s = args.seq_info
    file_a = args.assembly_info
    file_o = args.outfile
    # 解析文件到字典
    # dict_pj_info = get_info_dict(file_p)
    dict_seq_info = get_info_dict(file_s)
    dict_assembly_info = get_info_dict(file_a)
    # 生成结果文件
    get_result(pjid, spname, dict_seq_info, dict_assembly_info, file_o)


if __name__ == "__main__":
    main()
    




