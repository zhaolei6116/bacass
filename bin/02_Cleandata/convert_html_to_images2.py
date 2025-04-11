import os
import logging
import argparse
from kaleido.scopes.plotly import PlotlyScope
import plotly.io as pio

# 配置日志记录
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def convert_plotly_to_static_image(input_path, output_format):
    """
    将 Plotly 图表（JSON 或 HTML）转换为静态图片。
    
    参数:
        input_path (str): 输入文件路径（JSON 或 HTML）。
        output_format (str): 输出的图片格式（如 "png", "jpg", "pdf"）。
    """
    try:
        # 确保输入文件存在
        if not os.path.exists(input_path):
            logging.error(f"文件不存在: {input_path}")
            return
        
        # 加载 Plotly 图表对象
        if input_path.endswith(".json"):
            fig = pio.read_json(input_path)
        elif input_path.endswith(".html"):
            # 如果是 HTML 文件，尝试提取 JSON 数据
            try:
                with open(input_path, "r", encoding="utf-8") as f:
                    html_content = f.read()
                # 提取 JSON 数据（假设 HTML 中嵌入了 Plotly 的 JSON 数据）
                start_marker = 'Plotly.newPlot('
                end_marker = ');'
                json_str = html_content[html_content.find(start_marker) + len(start_marker):html_content.find(end_marker)]
                fig = pio.from_json(json_str)
            except Exception as e:
                logging.error(f"无法从 HTML 文件中提取 Plotly 数据: {e}")
                return
        else:
            logging.error("不支持的文件格式：仅支持 .json 和 .html 文件")
            return
        
        # 创建 PlotlyScope 对象
        # scope = PlotlyScope()
        
        # 构造输出文件路径（保持文件名不变，仅更改后缀）
        output_path = os.path.splitext(input_path)[0] + f".{output_format}"
        fig.write_image(output_path, format=output_format)
        # 使用 Kaleido 将图表转换为静态图片
        # image_bytes = scope.transform(fig, format=output_format)
        
        # 将生成的图片写入文件
        # with open(output_path, "wb") as f:
        #     f.write(image_bytes)
        
        logging.info(f"成功保存静态图片: {output_path}")
    
    except Exception as e:
        logging.error(f"转换失败: {e}")

def batch_convert_plots_to_images(input_dir, output_format):
    """
    批量将目录中的所有 Plotly 文件（JSON 或 HTML）转换为静态图片。
    
    参数:
        input_dir (str): 包含 Plotly 文件的目录路径。
        output_format (str): 输出的图片格式（如 "png", "jpg", "pdf"）。
    """
    # 获取目录中的所有 JSON 和 HTML 文件
    plot_files = [f for f in os.listdir(input_dir) if f.endswith((".json"))]
    
    if not plot_files:
        logging.warning(f"目录中未找到 Plotly Json 文件: {input_dir}")
        return
    
    logging.info(f"开始批量转换，共找到 {len(plot_files)} 个 Plotly Json 文件")
    
    for plot_file in plot_files:
        plot_path = os.path.join(input_dir, plot_file)
        convert_plotly_to_static_image(plot_path, output_format)

def main():
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="将 Plotly 图表文件（JSON 或 HTML）转换为静态图片")
    
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
        help="输入路径：如果是单独转换，指定单个文件路径；如果是批量转换，指定包含文件的目录路径"
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
        if not args.input.endswith((".json", ".html")):
            logging.error("单独转换模式下，输入路径必须是 .json 或 .html 文件")
            return
        convert_plotly_to_static_image(args.input, args.format)
    elif args.mode == "batch":
        # 批量转换
        if not os.path.isdir(args.input):
            logging.error("批量转换模式下，输入路径必须是目录")
            return
        batch_convert_plots_to_images(args.input, args.format)

if __name__ == "__main__":
    main()