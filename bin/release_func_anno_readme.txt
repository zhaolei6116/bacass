Release/
├── Assembly
│   ├── sample_name_contig_stat.xls -> 组装结果contig统计 
│   ├── sample_name_genome_stat.tsv -> 组装结果基因组统计
│   ├── sample_name.fasta ->  基因组组装结果fasta格式
│   ├── sample_name.fq ->  基因组组装fastq格式
│   └── output_figures -> 组装基因组测序深度图，按contig 分开展示
├── Decontamination
│   └── top10.tsv -> 物种检测，top10物种统计 
├── Genome_Annotion
│   ├── sample_name.embl -> EMBL格式的注释和序列  
│   ├── sample_name.faa -> CDS/sORF氨基酸FASTA序列
│   ├── sample_name.ffn -> CDS/sORF核苷酸FASTA序列
│   ├── sample_name.fna -> 基因组FASTA序列
│   ├── sample_name.gbff -> GenBank格式的注释和序列
│   ├── sample_name.gff3 -> GFF3格式的注释和序列
│   ├── sample_name.hypotheticals.faa -> 假设蛋白质CDS氨基酸FASTA序列
│   ├── sample_name.hypotheticals.tsv -> 假设蛋白质CDS的详细信息表
│   ├── sample_name.inference.tsv -> 注释结果的推断指标信息
│   ├── sample_name.json -> 注释和序列信息的json格式
│   ├── sample_name.png -> 圆形基因组注释图png格式 
│   ├── sample_name.svg -> 圆形基因组注释图svg格式
│   ├── sample_name.tsv -> 注释结果tsv格式
│   └── sample_name.txt -> 基因组注释总结
├── Gene_function_annotation
│   ├── CARD
│   │   └── sample_name_card_out.tsv ->  CARD 库注释结果统计表
│   ├── CAZy
│   │   ├── sample_name_cazy.out.tsv -> CAZy 库注释结果统计表
│   │   └── sample_name_cazy.png -> CAZy 库注释结果统计饼图
│   ├── COG
│   │   ├── sample_name_category_barplot.png -> COG数据库注释按类群展示柱状图
│   │   ├── sample_name_cog_category_counts.tsv -> COG数据库注释按类群统计结果
│   │   └── sample_name_cog_category.tsv -> COG数据库结果统计表
│   ├── GO
│   │   ├── sample_name_go_anno.tsv -> GO 数据库注释结果统计表
│   │   ├── sample_name_go.png -> GO 数据库按类群统计结果展示柱状图
│   │   └── sample_name_go_stats.tsv -> GO 数据库按类群统计结果
│   ├── KEGG
│   │   ├── sample_name_kegg_level2_counts.tsv -> 按pathway 类型 level2级别统计结果
│   │   ├── sample_name_kegg_pathway_mapping.tsv -> KEGG 数据库注释结果，KEGG ID 与 pathway 对应关系表。
│   │   └── sample_name_kegg_plot.png -> 按pathway 类型 level2级别统计柱状图
│   ├── NR
│   │   ├── sample_name_nr_species_count.tsv -> NR 数据库注释结果表
│   │   ├── sample_name_nr_top5.png -> 注释结果统计中前五的物种饼图
│   │   └── sample_name_nr_top5.tsv -> 注释结果统计中前五的物种统计表
│   ├── PHI
│   │   ├── sample_name_phi_info.tsv -> PHI 数据库注释结果表 
│   │   ├── sample_name_phi_Mutant_Phenotype_stat.png -> PHI 数据库注释结果按 Mutant_Phenotype 统计图片展示
│   │   ├── sample_name_phi_Mutant_Phenotype_stat.tsv -> PHI 数据库注释结果按 Mutant_Phenotype 统计
│   │   └── sample_name_phi_stat.tsv -> PHI 数据库注释结果按 PHI ID 统计
│   ├── Swissprot
│   │   └── sample_name_swissprot.out.tsv -> Swissprot 数据库注释结果统计表
│   └── VFDB
│   └── sample_name_vfdb_out.tsv -> VFDB 数据库注释结果统计表
└── Report
    └── sample_name_report.html -> 组装结果报告

