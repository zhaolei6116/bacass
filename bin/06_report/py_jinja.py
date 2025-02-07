import pandas as pd
from jinja2 import Environment, FileSystemLoader


'''
传递的参数
project_name
species
report_date 
project_sn
species
platform
total_pass_reads_bases
pipeline_png
sequennce_state_table
ont_reads_length_png
assemble_table
checkm_png
pla_state_table
genome_circos_png
plasmid_circos_png
software_table
'''




# 读取数据文件并生成 HTML 表格
# sequencing_df = pd.read_csv('table/seq_stat.tsv', sep='\t')   
software_df = pd.read_csv('table/software.tsv', sep='\t')
software_table_html = software_df.to_html(classes='table table-bordered table-striped table-hover', index=False)
# assembly_table_html = assembly_df.to_html(classes='table', index=False)

# 配置 Jinja2 模板环境
env = Environment(loader=FileSystemLoader('.'))
template = env.get_template('demo_template.html')

# 渲染模板的数据
# template_data = {
#     "report_title": "细菌完成图报告",
#     "report_subtitle": "样本分析报告",
#     "logo_path": "细菌完成图-plasmidsaurus/images/logo.png",
#     "sequencing_description": "测序数据统计的描述内容。",
#     "sequencing_table_title": "测序数据统计表",
#     "sequencing_table": sequencing_table_html,
#     "contamination_description": "数据污染评估的描述内容。",
#     "contamination_image": "细菌完成图-plasmidsaurus/images/top10_sp_hist.png",
#     "contamination_caption": "图1. 数据污染评估图",
#     "assembly_description": "基因组组装结果的描述内容。",
#     "assembly_table": assembly_table_html,
#     "circle_description": "基因组圈图的描述内容。",
#     "circle_image": "细菌完成图-plasmidsaurus/images/circle.png",
#     "circle_caption": "图2. 基因组圈图",
#     "methods_description": "分析方法的描述内容。"
# }
template_data = {
    
    "software_table": software_table_html 
}


# 渲染 HTML 内容
html_content = template.render(template_data)

# 输出到 HTML 文件
with open('report.html', 'w', encoding='utf-8') as f:
    f.write(html_content)
