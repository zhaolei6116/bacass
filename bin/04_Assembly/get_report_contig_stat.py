#!/usr/bin/env python3

import argparse
import csv
from collections import namedtuple

__version__ = "v0.0.1"
__date__    = "2024.11.25"
__author__  = "zhaolei"


def gc_parse(infile):
    '''
    infile format : chr, length, #A, #C, #G, #T
    '''
    gc_ratio = {}
    contig_list = []
    with open(infile, 'r', encoding='utf-8') as inF:
        reader = csv.reader(inF, delimiter='\t')
        for row in reader:
            gc_ratio[row[0]] = (int(row[3]) + int(row[4]))/int(row[1]) * 100
            gc_ratio[row[0]] = round(gc_ratio[row[0]], 2)
            contig_list.append(row[0])
    
    return gc_ratio, contig_list


def cov_dep_parse(infile):

    cov_dep = {}
    with open(infile, 'r', encoding='utf-8') as inF:
        reader = csv.reader(inF, delimiter='\t')
        next(reader)
        for row in reader:
            cov_dep[row[0]] = {}
            cov_dep[row[0]]['length'] = f'{int(row[2]):,}'
            cov_dep[row[0]]['cov'] = format_number_str(row[5])
            cov_dep[row[0]]['dep'] = f'{float(row[6]):.2f}' 
    
    return cov_dep

def get_circ_info(infile):
    circs = {} 
    with open(infile,  newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        print(reader)
        print(reader.fieldnames)
        # 修改表头字段名，去掉开头的 '#' 和结尾的 '.'
        new_fieldnames = [field.lstrip('#').rstrip('.') for field in reader.fieldnames]
        print(new_fieldnames)

        for row in reader:
             circs[row['new_seq_name']] = row['circ.']
        
    return circs





def format_number_str(num_str):
    num = float(num_str)  # 将字符串转换为浮点数
    # 检查是否有小数部分，并根据小数位数进行格式化
    if num % 1 == 0:
        return str(int(num))  # 没有小数则不显示小数点
    elif num * 10 % 1 == 0:
        return f'{num:.1f}'  # 有一位小数则显示一位
    else:
        return f'{num:.2f}'  # 有两位小数则显示两位



def get_result(gc_ratio, contig_list, cov_dep, circs, outfile):

    head_line = ["Contig", "Length", "GC content(%)", "Coverage", "Sequence depth(X)", "Circular"]

    with open(outfile, "w") as outF:
        outF.write("\t".join(head_line) + "\n")
        for contig in contig_list:
            temp_list = [contig, cov_dep[contig]['length'], str(gc_ratio[contig]), cov_dep[contig]['cov'], cov_dep[contig]['dep'], circs[contig]]
            outF.write("\t".join(temp_list) + "\n")
    

def setup_args():
    parser = argparse.ArgumentParser(description = 'get contig stat',
            formatter_class=argparse.RawTextHelpFormatter,
            epilog = f"Version: {__version__}\nDate: {__date__}\nAuthor: {__author__}")
    parser.add_argument('-v','--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-gc', type=str, help='The path to gc stat', required = True)
    parser.add_argument('-cov', type=str, help='The path to cov and depth stat', required = True)
    parser.add_argument('-circ', type=str, help='The path to flye stat', required = True)
    parser.add_argument('-o', type=str, help='output file', default="contig_stat.xls")
    args = parser.parse_args()

    return args


def main():
    args = setup_args()
    gc_ratio, contig_list = gc_parse(args.gc)
    cov_dep = cov_dep_parse(args.cov)
    circs = get_circ_info(args.circ)
    get_result(gc_ratio, contig_list, cov_dep, circs, args.o)



if __name__ == "__main__":
    main()

