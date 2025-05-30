Release/
├── Assembly
│   ├── sample_name_contig_stat.xls -> 组装结果contig统计 
│   ├── sample_name_genome_stat.tsv -> 组装结果基因组统计
│   ├── sample_name.fasta ->  基因组组装结果
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
├── Rawdata
│   └── sample_name.raw.fastq.gz -> 测序数据 
└── Report
    └── sample_name_report.html -> 组装结果报告

