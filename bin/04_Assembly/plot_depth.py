#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from concurrent.futures import ProcessPoolExecutor
import time

def xaxis_major_formatter(x, pos=None):
    if x % 1000 == 0:
        return f"{int(x / 1000)}k"
    elif x % 100 == 0:
        return f"{x / 1000:.1f}k"
    return ''

def calculate_coverage(bin_depth, thresholds):
    length = len(bin_depth)
    coverage = [(i, (bin_depth >= i).sum() / length) for i in thresholds]
    return coverage

def plot_single_contig(args):
    sub_df, figure_dir, sample_name, group = args
    length = len(sub_df)
    
    # 计算覆盖度
    thresholds = [10, 20, 50, 100, 200, 500]
    coverage = calculate_coverage(sub_df['depth'], thresholds)
    title = [f">={i}X: {cov:.2%}" for i, cov in coverage]
    plot_title = f"{sample_name}\n" + \
                 f"{title[0]:^20}{title[1]:^20}{title[2]:^20}\n" + \
                 f"{title[3]:^20}{title[4]:^20}{title[5]:^20}"
    
    # 计算每个bin的中间位置
    positions = (sub_df['start'] + sub_df['end']) / 2
    depth_avg = sub_df['depth'].values
    bin_widths = sub_df['end'] - sub_df['start']
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(positions, depth_avg, width=bin_widths, color='blue', align='center', edgecolor='none')
    ax.set_xlabel('Position')
    ax.set_ylabel("Mean Depth")
    ax.set_title(plot_title)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=10))
    ax.xaxis.set_major_formatter(xaxis_major_formatter)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis="x", labelsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(figure_dir, f"{sample_name}-{group}.coverage.png"))
    plt.close(fig)  # 释放内存

def plot_coverage(pre_binned_dp, figure_dir, sample_name):
    # 假设预下采样的数据文件有四列：contig, start, end, mean_depth
    dtype = {
        'contig': 'category',
        'start': np.int32,
        'end': np.int32,
        'depth': np.float32
    }
    df = pd.read_csv(pre_binned_dp, sep="\t", header=None, names=['contig', 'start', 'end', 'depth'], dtype=dtype)
    
    grouped = df.groupby("contig", observed=True)
    args_list = []
    for group in grouped.groups:
        sub_df = grouped.get_group(group)
        args_list.append((sub_df, figure_dir, sample_name, group))
    
    # 并行处理多个contig
    with ProcessPoolExecutor() as executor:
        executor.map(plot_single_contig, args_list)

if __name__ == '__main__':
    usage = f"usage: {os.path.basename(__file__)} <pre_binned_dp.bed> <fig_dir> <sample_name>"
    if len(sys.argv) != 4:
        raise Exception(usage + "\n")
    pre_binned_dp = sys.argv[1]
    figure_dir = sys.argv[2]
    sample = sys.argv[3]
    
    # 确保输出目录存在
    os.makedirs(figure_dir, exist_ok=True)
    start_time = time.time()
    plot_coverage(pre_binned_dp, figure_dir, sample)
    end_time = time.time()
    print(f"Total runtime: {end_time - start_time:.2f} seconds")

