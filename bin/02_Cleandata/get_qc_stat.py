#!/usr/bin/env python3

import csv
import argparse

__version__ = "v0.0.1"
__date__    = "2024.11.21"
__author__  = "zhaolei"



class MetricsReader:
    def __init__(self, filename):
        self._filename = filename
        self._read_file_and_set_attributes()
        self.num_format_set()

    def _read_file_and_set_attributes(self):
        with open(self._filename, 'r', encoding='utf-8') as file:
            reader = csv.reader(file, delimiter='\t')
            next(reader)  # 
            for row in reader:
                if row:  #
                    clean_key = self.clean_key(row[0])
                    try:
                        value = float(row[1].strip())
                    except ValueError:
                        value = row[1].strip()
                    setattr(self, clean_key, value)

    def clean_key(self, key):
        # 替换掉非法字符
        return ''.join(char if char.isalnum() or char == '_' else '' for char in key)

    def num_format_set(self):
        self.number_of_bases = f'{int(self.number_of_bases):,}'
        self.number_of_reads = f'{int(self.number_of_reads):,}'
        self.n50 = f'{int(self.n50):,}'
        self.longest_read_with_Q1 = self.convert_num(self.longest_read_with_Q1)
        self.mean_read_length = f'{self.mean_read_length:,}'
        self.median_read_length = f'{self.median_read_length:,}'
        self.ReadsQ5 = self.convert_num_2(self.ReadsQ5)
        self.ReadsQ10 = self.convert_num_2(self.ReadsQ10)
        self.ReadsQ15 = self.convert_num_2(self.ReadsQ15)




    def convert_num(self, innum):
        innum = f'{int(innum.split(" ")[0]):,}' + " " +innum.split(" ")[1]
        return innum

    def convert_num_2(self, innum):
        arr = innum.split(" ")
        innum = f'{int(arr[0]):,}' + " " + arr[1] + " " + arr[2]
        return innum


def setup_args():
    parser = argparse.ArgumentParser(description = 'merge raw and clean state info into qc table',
            formatter_class=argparse.RawTextHelpFormatter,
            epilog = f"Version: {__version__}\nDate: {__date__}\nAuthor: {__author__}")
    parser.add_argument('-v','--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-raw', type=str, help='The path to raw_stat.tsv', required = True)
    parser.add_argument('-clean', type=str, help='The path to clean_stat.tsv', required = True)
    parser.add_argument('-o', type=str, help='output file', default="qc_stat.xls")
    args = parser.parse_args()

    return args    
    

def merge_stat(raw_stat_sub, clean_stat_sub, qc_stat):
    
    with open(qc_stat, "w") as qc:
        head_line = ["Type","Number of base", "Number of Reads", "Reads N50", "Longest Reads(Q)", "Mean Reads Length", "Median Reads Length", "Mean Quality", "Reads >Q10", "Reads >Q15"]
        attr_list = ["number_of_bases", 'number_of_reads', 'n50', 'longest_read_with_Q1', 'mean_read_length', 'median_read_length', 'mean_qual', 'ReadsQ10', 'ReadsQ15']

        raw_stat_line = ["Raw Data"]
        clean_stat_line = ["Clean Data"]

        
        for attr in attr_list:
            raw_stat_line.append(str(getattr(raw_stat_sub, attr, 'None')))
            clean_stat_line.append(str(getattr(clean_stat_sub, attr, 'None')))
        
        for line in [head_line, raw_stat_line, clean_stat_line]:
            qc.write("\t".join(line)+"\n")
        



if __name__ == "__main__":
    # 获取输入
    args = setup_args()
    raw_stat = args.raw
    clean_stat = args.clean
    qc_stat = args.o
    
    # 生成nanostat对象
    raw_stat_sub =  MetricsReader(raw_stat)
    clean_stat_sub = MetricsReader(clean_stat)
    
    # 合并数据，生成输出文件
    merge_stat(raw_stat_sub, clean_stat_sub, qc_stat)

    
    # 访问属性
    # print(dir(metrics_reader))  # 打印属性名

    # ['ReadsQ10', 'ReadsQ12', 'ReadsQ15', 'ReadsQ5', 'ReadsQ7', 'clean_key', 'filename', 
    #  'highest_Q_read_with_length1', 'highest_Q_read_with_length2', 'highest_Q_read_with_length3', 
    #  'highest_Q_read_with_length4', 'highest_Q_read_with_length5', 'longest_read_with_Q1', 
    #  'longest_read_with_Q2', 'longest_read_with_Q3', 'longest_read_with_Q4', 'longest_read_with_Q5', 
    #  'mean_qual', 'mean_read_length', 'median_qual', 'median_read_length', 'n50', 'number_of_bases', 'number_of_reads', 'read_length_stdev']
