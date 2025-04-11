#!/bin/bash

# 获取 run.sh 的绝对路径
RUN_SCRIPT_PATH=$(readlink -f "$0")

# 获取 run.sh 所在的目录
SCRIPT_DIR=$(dirname "$RUN_SCRIPT_PATH")

# 激活 Conda 环境并运行 get_info.py
source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate

# 定义所有区域
# REGIONS=("WH" "SH" "BJ" "GZ" "TZ")

python $SCRIPT_DIR/get_and_management_info.py  

