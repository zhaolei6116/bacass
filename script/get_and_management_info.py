
import logging
import os
from datetime import datetime, timedelta
import json
import argparse
import subprocess



def load_config(config_file):
    """
    加载配置文件
    :param config_file: 配置文件路径
    :return: 配置字典
    """
    
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"配置文件 {config_file} 不存在！")
    
    with open(config_file, "r") as f:
        config = json.load(f)

    # 获取配置文件所在目录    
    config_dir = os.path.dirname(os.path.abspath(config_file))

    # 将相对路径转换为绝对路径
    for key in ["backup_dir", "java_jar"]:
        if key in config and not os.path.isabs(config[key]):
            config[key] = os.path.join(config_dir, config[key])
    
    # 转换 java_configs 中的所有路径为绝对路径
    if "java_configs" in config:
        for region, path in config["java_configs"].items():
            if not os.path.isabs(path):
                config["java_configs"][region] = os.path.join(config_dir, path)

    return config


def setup_logging(log_dir="logs"):
    """
    设置日志路径和文件名，并初始化日志系统
    :param log_dir: 日志目录
    :return: None
    """
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    log_file = os.path.join(log_dir, f"log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
    
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    logging.getLogger().addHandler(console)
    
    logging.info(f"日志系统初始化完成，日志文件：{log_file}")

def check_sample_info_db(sample_info_path="sample_info_db.tsv"):
    """
    检查 sample_info_db.tsv 是否存在
    :param sample_info_path: 样本信息文件路径
    :return: bool
    """
    logging.info(f"检测样本信息库文件")
    if os.path.exists(sample_info_path):
        logging.info(f"{sample_info_path} 文件存在，准备加载数据")
        return True
    else:
        logging.info(f"{sample_info_path} 文件不存在，将在首次加载数据时创建")
        return False


def fetch_sample_data(java, start_time, end_time, java_jar, java_config, region):
    """
    调用 Java 接口获取样本信息
    :param java: Java 可执行文件路径
    :param start_time: 起始时间，格式 "YYYY-MM-DD HH:MM:SS"
    :param end_time: 结束时间，格式 "YYYY-MM-DD HH:MM:SS"
    :param java_jar: Java 程序路径
    :param java_config: Java 配置文件路径
    :param region: 当前处理的区域
    """
    logging.info(f"调用 Java 接口获取样本信息 (区域: {region}, 配置文件: {java_config})")
    cmd = (
        f'{java} -jar {java_jar} '
        f'--startTime="{start_time}" --endTime="{end_time}" '
        f'--config="{java_config}"'
    )

    logging.info(f"Executing command:{cmd}")
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            if "\"data\":" in result.stdout:
                logging.info("样本信息获取成功")
                logging.info(f"Java 接口返回信息:\n{result.stdout}\n{result.stderr}")
                if "\"data\":[]" in result.stdout:
                    logging.warning(f"区域 {region} 返回结果为空，无满足条件的数据。\n\n"+"*"*50)
                else:
                    logging.info(f"区域 {region} 返回结果包含有效数据。\n\n"+"*"*50)
            else:
                logging.error(f"区域 {region} 样本信息获取失败")
                logging.error(f"Java 接口错误信息:\n{result.stdout}\n{result.stderr}")
        else:
            logging.error(f"区域 {region} 样本信息获取失败")
            logging.error(f"Java 接口错误信息:\n{result.stdout}\n{result.stderr}")
    except Exception as e:
        logging.critical(f"调用 Java 接口时发生错误 (区域: {region}): {e}")

def setup_args():
    """
    设置外部参数输入。
    支持的参数：
    -c 或 --config: 配置文件路径，默认为 "config.json"。
    -s 或 --sample_info: 样本信息文件路径，默认为 "sample_info_db.tsv"。
    """
    script_absolute_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_absolute_path)

    parser = argparse.ArgumentParser(description="样本信息处理脚本")

    # 配置文件路径参数
    
    parser.add_argument(
        "-c", "--config",
        type=str,
        required=False,
        default=os.path.join(script_dir, "config.json"),
        help=f"配置文件路径，默认为 {os.path.join(script_dir, "config.json")}"
    )

    # 添加区域参数
    # parser.add_argument(
    #     "-r", "--region",
    #     required=True,
    #     choices=["WH", "SH", "BJ", "GZ", "TZ"],
    #     help="Region to use"
    # )

    # 样本信息文件路径参数
    # parser.add_argument(
    #     "-s", "--sample_info",
    #     type=str,
    #     default="sample_info_db.tsv",
    #     help="样本信息文件路径，默认为 'sample_info_db.tsv'"
    # )

    return parser.parse_args()

def main():
    # 获取命令行参数
    args = setup_args()
    

    # 加载配置文件
    config = load_config(args.config)

    # 初始化日志系统
    setup_logging(log_dir=config["log_dir"])

    # 样本信息库文件路径
    # sample_info_path = args.sample_info

    # 检查样本信息文件
    # check_out = check_sample_info_db(sample_info_path)

    # 定义所有区域
    regions = ["WH", "TZ", "BJ", "GZ", "SH"]

    # 定义时间范围
    end_time = datetime.now()
    start_time = end_time - timedelta(hours=config["START_OFFSET"])

    # 格式化时间
    start_time_str = start_time.strftime("%Y-%m-%d %H:%M:%S")
    end_time_str = end_time.strftime("%Y-%m-%d %H:%M:%S")

    # 循环处理每个区域
    for region in regions:
        logging.info(f"Processing region: {region}")
        java_config = config["java_configs"].get(region)
        if not java_config:
            logging.error(f"No java_config found for region: {region}")
            continue

        # 调用 Java 接口
        fetch_sample_data(
            config["java"],
            start_time_str,
            end_time_str,
            config["java_jar"],
            java_config,
            region
        )



if __name__ == "__main__":
    main()
    



    # 调用Java接口
    # fetch_out = fetch_sample_data(config["java"],
    #                               start_time.strftime("%Y-%m-%d %H:%M:%S"), 
    #                               end_time.strftime("%Y-%m-%d %H:%M:%S"), 
    #                               config["java_jar"], config["java_configs"][args.region])
    
    
    

    

    



