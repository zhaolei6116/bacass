
import os
import time
import pandas as pd
from datetime import datetime, timedelta
import subprocess

# 配置参数
START_OFFSET = 10  # 前10小时
CHECK_INTERVAL = 3600  # 每小时执行
DATA_DIR = "LimsData"
ALL_TSV = "all.tsv"
java_config = ""
java_jar = ""

# 初始化 all.tsv 文件
if not os.path.exists(ALL_TSV):
    pd.DataFrame(columns=["SampleID", "Status", "StartTime", "EndTime"]).to_csv(ALL_TSV, sep="\t", index=False)

def fetch_data(start_time, end_time):
    """调用Java接口获取样本信息"""
    cmd = f'java -jar {java_jar} --startTime="{start_time}" --endTime="{end_time}" --config "{java_config}"'
    subprocess.run(cmd, shell=True)

def parse_sample_files():
    """解析样本信息"""
    all_samples = []
    for root, dirs, files in os.walk(DATA_DIR):
        for file in files:
            if file.endswith(".xls"):
                sample_path = os.path.join(root, file)
                sample_data = pd.read_excel(sample_path, skiprows=1)
                all_samples.append(sample_data)
    return pd.concat(all_samples, ignore_index=True)

def update_all_tsv(new_samples):
    """更新 all.tsv 文件"""
    all_df = pd.read_csv(ALL_TSV, sep="\t")
    for _, sample in new_samples.iterrows():
        if sample["SampleID"] not in all_df["SampleID"].values:
            all_df = all_df.append({
                "SampleID": sample["SampleID"],
                "Status": "未分析",
                "StartTime": sample["StartTime"],
                "EndTime": sample["EndTime"]
            }, ignore_index=True)
    all_df.to_csv(ALL_TSV, sep="\t", index=False)

def start_analysis():
    """启动未分析样本的分析流程"""
    all_df = pd.read_csv(ALL_TSV, sep="\t")
    pending_samples = all_df[all_df["Status"] == "未分析"]
    for _, sample in pending_samples.iterrows():
        sample_id = sample["SampleID"]
        analysis_path = os.path.join("analysis", sample_id)
        os.makedirs(analysis_path, exist_ok=True)

        # 创建输入文件
        input_tsv = os.path.join(analysis_path, "input.tsv")
        with open(input_tsv, "w") as f:
            f.write("SampleID\tPath\n")
            f.write(f"{sample_id}\t/sample/path/{sample_id}\n")

        # 启动分析
        cmd = f"source /Miniconda3/miniconda3/bin/activate /conda/env/nextflow/v24.04.4 && "
        cmd += f"nextflow run /path/to/main.nf --input {input_tsv} --id {sample_id} -resume"
        subprocess.run(cmd, shell=True)

        # 更新状态
        all_df.loc[all_df["SampleID"] == sample_id, "Status"] = "分析中"
    all_df.to_csv(ALL_TSV, sep="\t", index=False)

# 主流程
while True:
    end_time = datetime.now()
    start_time = end_time - timedelta(hours=START_OFFSET)
    fetch_data(start_time.strftime("%Y-%m-%d %H:%M:%S"), end_time.strftime("%Y-%m-%d %H:%M:%S"))

    new_samples = parse_sample_files()
    update_all_tsv(new_samples)
    start_analysis()

    time.sleep(CHECK_INTERVAL)
