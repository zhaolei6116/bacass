import os
import logging
import argparse
from kaleido.scopes.plotly import PlotlyScope
import plotly.io as pio

# 配置日志记录
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def convert_html_to_static_image(html_path, output_format):
    """
    将 HTML 文件转换为静态图片。
    
    参数:
        html_path (str): 输入的 HTML 文件路径。
        output_format (str): 输出的图片格式（如 "png", "jpg", "pdf"）。
    """
    try:
        # 确保输入文件存在
        if not os.path.exists(html_path):
            logging.error(f"HTML 文件不存在: {html_path}")
            return
        
        # 读取 HTML 文件并加载为 Plotly 图表对象
        fig = pio.read_html(html_path)
        
        # 创建 PlotlyScope 对象
        scope = PlotlyScope()
        
        # 构造输出文件路径（保持文件名不变，仅更改后缀）
        output_path = os.path.splitext(html_path)[0] + f".{output_format}"
        
        # 使用 Kaleido 将图表转换为静态图片
        image_bytes = scope.transform(fig, format=output_format)
        
        # 将生成的图片写入文件
        with open(output_path, "wb") as f:
            f.write(image_bytes)
        
        logging.info(f"成功保存静态图片: {output_path}")
    
    except Exception as e:
        logging.error(f"转换失败: {e}")

def batch_convert_html_to_images(input_dir, output_format):
    """
    批量将目录中的所有 HTML 文件转换为静态图片。
    
    参数:
        input_dir (str): 包含 HTML 文件的目录路径。
        output_format (str): 输出的图片格式（如 "png", "jpg", "pdf"）。
    """
    # 获取目录中的所有 HTML 文件
    html_files = [f for f in os.listdir(input_dir) if f.endswith(".html")]
    
    if not html_files:
        logging.warning(f"目录中未找到 HTML 文件: {input_dir}")
        return
    
    logging.info(f"开始批量转换，共找到 {len(html_files)} 个 HTML 文件")
    
    for html_file in html_files:
        html_path = os.path.join(input_dir, html_file)
        convert_html_to_static_image(html_path, output_format)

def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="将 HTML 文件转换为静态图片")
    
    # 添加命令行参数
    parser.add_argument(
        "-m", "--mode", 
        choices=["single", "batch"], 
        required=True, 
        help="转换模式：'single' 表示单独转换，'batch' 表示批量转换"
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="输入路径：如果是单独转换，指定单个 HTML 文件路径；如果是批量转换，指定包含 HTML 文件的目录路径"
    )
    parser.add_argument(
        "-f", "--format", 
        default="png", 
        choices=["png", "jpg", "pdf"], 
        help="输出图片格式，默认为 'png'"
    )
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 根据模式执行转换
    if args.mode == "single":
        # 单独转换
        if not args.input.endswith(".html"):
            logging.error("单独转换模式下，输入路径必须是 .html 文件")
            return
        convert_html_to_static_image(args.input, args.format)
    elif args.mode == "batch":
        # 批量转换
        if not os.path.isdir(args.input):
            logging.error("批量转换模式下，输入路径必须是目录")
            return
        batch_convert_html_to_images(args.input, args.format)

if __name__ == "__main__":
    main()