import base64
import configparser
import pandas as pd
from jinja2 import Environment, FileSystemLoader
import os
import re

__version__ = "v0.0.2"
__date__    = "2024.12.09"
__author__  = "zhaolei"


class ReportGenerator:
    def __init__(self, config_path, template_path):
        self.config = self._parse_config(config_path)
        self.env = Environment(loader=FileSystemLoader(os.path.dirname(template_path)))
        self.template = self.env.get_template(os.path.basename(template_path))
        self.data = self._process_config_attributes()
        
    def _parse_config(self, config_path):
        """解析配置文件并返回配置对象"""
        config = configparser.ConfigParser()
        config.read(config_path, encoding='utf-8')
        return config['ReportConfig']

    def _process_config_attributes(self):
        """根据属性名规则处理表格和图片，生成渲染所需的数据字典"""
        data = {}
        
        for key, value in self.config.items():
            if key.endswith("table"):
                if key == "pj_summary_table":
                    data[key] = self._table_to_html2(value)
                else:    
                    data[key] = self._table_to_html(value)
            elif key.endswith("png"):
                data[key] = self._image_to_base64(value)
            else:
                data[key] = value
        print(data.keys())
        return data

    def _table_to_html(self, table_path):
        """将 CSV 表格文件转换为 HTML 表格字符串"""
        try:
            df = pd.read_csv(table_path, sep='\t')
            return df.to_html(index=False, classes='table table-bordered table-striped table-hover')
        except Exception as e:
            print(f"Error reading table {table_path}: {e}")
            return "<p>表格加载失败</p>"
    
    def _table_to_html2(self, table_path):
        """将 CSV 表格文件转换为 HTML 表格字符串"""
        try:
            df = pd.read_csv(table_path, sep='\t',  header=None)
            temp_html = df.to_html(header=False, index=False, classes='table table-bordered table-striped table-hover')
            temp_html = re.sub(r'<td>物种名称</td>\n\s*<td>', r'<td>物种名称</td>\n\s*<td style="font-style: italic;">')
            return temp_html
        except Exception as e:
            print(f"Error reading table {table_path}: {e}")
            return "<p>表格加载失败</p>"


    def _image_to_base64(self, image_path):
        """将图片文件转换为 Base64 编码字符串"""
        try:
            with open(image_path, "rb") as img_file:
                return base64.b64encode(img_file.read()).decode('utf-8')
        except Exception as e:
            print(f"Error reading image {image_path}: {e}")
            return ""

    def render_to_html(self, output_path):
        """将数据渲染到模板并生成 HTML 报告文件"""
        html_content = self.template.render(self.data)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"报告已生成：{output_path}")


def setup_args():
    script_absolute_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_absolute_path)

    parser = argparse.ArgumentParser(description = 'generate bacteria genome assembly report',
            formatter_class=argparse.RawTextHelpFormatter,
            epilog = f"Version: {__version__}\nDate: {__date__}\nAuthor: {__author__}")
    parser.add_argument('-v','--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-c', '--config', type=str, help='The path to config', required = True, default=script_dir+"/config.ini") 
    parser.add_argument('-t', '--template', type=str, help='The path to template report', required = True, default=script_dir+"/report.html")
    parser.add_argument('-p', '--prefix', type=str, help='output file prefix, 一般为项目号', default=None)
    args = parser.parse_args()

    return args  


def main():
    args = setup_args()
    config_path = args.config
    template_path = args.template
    out_prefix = args.prefix

    if out_prefix:
        output_path = os.getcwd()+ out_prefix + "report.html"
    else:
        output_path = os.getcwd() + "report.html"
    

    report_generator = ReportGenerator(config_path, template_path)
    report_generator.render_to_html(output_path)

# 使用示例
if __name__ == "__main__":
    config_path = "config.ini"
    template_path = "demo_template.html"
    output_path = "report.html"

    report_generator = ReportGenerator(config_path, template_path)
    report_generator.render_to_html(output_path)
