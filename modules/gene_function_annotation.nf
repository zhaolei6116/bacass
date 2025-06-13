nextflow.enable.dsl=2

process extract_gene_id {
  // publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("go.list"), path("kegg.list"), path("cog.list"), emit:sample_info_tuple
    
  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  """
  python3 ${projectDir}/bin/gene_annotation/extract_ids_from_gff.py  -i ${sample_info_map.cds_gff3}  --go_output go.list --kegg_output kegg.list --cog_output cog.list
  """
}

process nr_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/NR/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_nr_species_count.tsv"), path("${sample_info_map.sample_label}_nr_top5.tsv"), path("${sample_info_map.sample_label}_nr_top5.png"), emit:sample_info_tuple
    
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  """
  # 比对
  ${params.software.diamond}  blastp \
      --db ${params.database.nr} \
      -q ${sample_info_map.cds_faa}  \
      -e 1e-5  --query-cover 50 --subject-cover 30  \
      --outfmt  6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
      --max-target-seqs 3 \
      -o  ${sample_info_map.sample_label}_nr_blast.xls
  
  # 获取taxid
  cut -f 13 ${sample_info_map.sample_label}_nr_blast.xls| tr ';' '\\n' | sort | uniq -c |  awk '\$2 ~ /^[0-9]+\$/ {print \$2 "\\t" \$1}' | sort -k2,2 -nr |${params.software.taxonkit}   lineage  -n  --data-dir ${params.database.taxonomy} |sed 1i'#Taxid\\tCount\\tLineage\\tSpecies' > ${sample_info_map.sample_label}_nr_seq2taxid.list
  
  # 得到物种和计数
  sh  ${projectDir}/bin/gene_annotation/NR/tax_sort_count.sh  ${sample_info_map.sample_label}_nr_seq2taxid.list  ${sample_info_map.sample_label}_nr_species_count.tsv  ${sample_info_map.sample_label}_nr_top5.tsv
  
  # 作图
  ${params.software.python} ${projectDir}/bin/gene_annotation/NR/generate_nr_pie_chart.py -i ${sample_info_map.sample_label}_nr_top5.tsv -o ${sample_info_map.sample_label}_nr_top5.png

  """

}

process go_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/GO/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val(sample_info_map)

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_go_anno.tsv"), path("${sample_info_map.sample_label}_go_stats.tsv"), path("${sample_info_map.sample_label}_go.png"),  emit:sample_info_tuple

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  ${params.software.python}  ${projectDir}/bin/gene_annotation/GO/get_go_info.py   -i ${params.database.go_term}  -g  ${sample_info_map.go_list} -o ${sample_info_map.sample_label}_go_anno.tsv   -s ${sample_info_map.sample_label}_go_stats.tsv  --plot_output ${sample_info_map.sample_label}_go.png
  """

}

process cog_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/COG/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_cog_category.tsv"), path("${sample_info_map.sample_label}_cog_category_counts.tsv"),  path("${sample_info_map.sample_label}_category_barplot.png"), emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  """
  ${params.software.python} ${projectDir}/bin/gene_annotation/COG/cog_info_extractor.py  -d  ${params.database.cog_def}  -f ${params.database.cog_fun}  -i  ${sample_info_map.cog_list}   -o1 ${sample_info_map.sample_label}_cog_category.tsv  -o2 ${sample_info_map.sample_label}_cog_category_counts.tsv

  ${params.software.python} ${projectDir}/bin/gene_annotation/COG/cog_plot_category_bar.py -i ${sample_info_map.sample_label}_cog_category_counts.tsv -o ${sample_info_map.sample_label}_category_barplot.png
  """

}

process kegg_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/KEGG/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 5
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_kegg_pathway_mapping.tsv"), path("${sample_info_map.sample_label}_kegg_level2_counts.tsv"), path("${sample_info_map.sample_label}_kegg_plot.png"),  emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  ${params.software.python} ${projectDir}/bin/gene_annotation/KEGG/get_kegg_infor.py  -d1 ${params.database.kegg2pathway_json} -d2 ${params.database.full_kegg_pathway_json} -i ${sample_info_map.kegg_list}  -o1 ${sample_info_map.sample_label}_kegg_pathway_mapping.tsv -o2 ${sample_info_map.sample_label}_kegg_level2_counts.tsv 
  ${params.software.python} ${projectDir}/bin/gene_annotation/KEGG/kegg_plot.py  -i ${sample_info_map.sample_label}_kegg_level2_counts.tsv -o ${sample_info_map.sample_label}_kegg_plot.png
  """

}

process swissprot_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/Swissprot/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_swissprot.out.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  """
  # 比对
  ${params.software.diamond}  blastp --db ${params.database.swiss_prot}  -q ${sample_info_map.cds_faa}  -e 1e-5  --query-cover 50 --subject-cover 30  -o  ${sample_info_map.sample_label}_swiss-prot_blast.xls --max-target-seqs 2
  # 信息获取
  ${params.software.python} ${projectDir}/bin/gene_annotation/Swissprot/parse_blast_swissprot.py  -b  ${sample_info_map.sample_label}_swiss-prot_blast.xls  -s  ${params.database.swiss_prot_id_description} -o  ${sample_info_map.sample_label}_swissprot.out.tsv
  """

}

process card_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/CARD/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_card_out.tsv"),   emit:sample_info_tuple
    
  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  # 比对
  ${params.software.diamond}  blastp --db ${params.database.card}  -q ${sample_info_map.cds_faa}  -e 1e-5  --query-cover 50 --subject-cover 30  -o  ${sample_info_map.sample_label}_card_blast.xls --max-target-seqs 2
  # 信息获取
  ${params.software.python}  ${projectDir}/bin/gene_annotation/CARD/get_card_info.py  -b ${sample_info_map.sample_label}_card_blast.xls  -c ${params.database.card_aro_json}  -o ${sample_info_map.sample_label}_card_out.tsv

  """

}

process phi_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/PHI/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_phi_info.tsv"), path("${sample_info_map.sample_label}_phi_Mutant_Phenotype_stat.tsv"),  path("${sample_info_map.sample_label}_phi_stat.tsv"), path("${sample_info_map.sample_label}_phi_Mutant_Phenotype_stat.png"), emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  # 比对
  ${params.software.diamond}  blastp --db ${params.database.phi}  -q ${sample_info_map.cds_faa}  -e 1e-5  --query-cover 50 --subject-cover 30  -o  ${sample_info_map.sample_label}_phi_blast.xls --max-target-seqs 2
  
  # 信息获取
  ${params.software.python} ${projectDir}/bin/gene_annotation/PHI/phi_anno.py  -p  ${params.database.phi_base_info} -b ${sample_info_map.sample_label}_phi_blast.xls  -o1 ${sample_info_map.sample_label}_phi_info.tsv  -o2 ${sample_info_map.sample_label}_phi_Mutant_Phenotype_stat.tsv -o3 ${sample_info_map.sample_label}_phi_stat.tsv

  ${params.software.python} ${projectDir}/bin/gene_annotation/PHI/phi_plot.py -i ${sample_info_map.sample_label}_phi_Mutant_Phenotype_stat.tsv  -o ${sample_info_map.sample_label}_phi_Mutant_Phenotype_stat.png

  """
}

process vfdb_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/VFDB/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_vfdb_out.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  
  script:
  """
  # 比对
  ${params.software.diamond}  blastp --db ${params.database.vfdb}  -q ${sample_info_map.cds_faa}  -e 1e-20  --query-cover 50 --subject-cover 30  -o  ${sample_info_map.sample_label}_vfdb_blast.xls --max-target-seqs 1
  
  # 信息获取
  ${params.software.python} ${projectDir}/bin/gene_annotation/VFDB/annotate_vfdb.py  --setb  ${params.database.vfdb_setb} --vfs ${params.database.vfdb_vfs} --blast ${sample_info_map.sample_label}_vfdb_blast.xls  --output ${sample_info_map.sample_label}_vfdb_out.tsv
  
  """
}

process cazy_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/CAZy/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map),  path("${sample_info_map.sample_label}_cazy.out.tsv"), path("${sample_info_map.sample_label}_cazy.png"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  
  script:
  """
  # 比对
  # run_dbcan ${sample_info_map.cds_faa}  protein --out_dir ./ --dia_cpu 32  --tools diamond  --out_pre ${sample_info_map.sample_label}_  --db_dir ${params.database.cazy_dbcan}

  ${params.software.diamond} blastp --db ${params.database.cazy_diamond} -q ${sample_info_map.cds_faa} -e 1e-102   --max-target-seqs 1 -o  ${sample_info_map.sample_label}_diamond.tsv

  # 信息获取
  ${params.software.python} ${projectDir}/bin/gene_annotation/CAZy/cazy_class_count.py  --tsv ${sample_info_map.sample_label}_diamond.tsv  --json ${params.database.cazy_level1_json}  --output ${sample_info_map.sample_label}_cazy.out.tsv

  ${params.software.python} ${projectDir}/bin/gene_annotation/CAZy/cazy_generate_pie_pic.py  -i ${sample_info_map.sample_label}_cazy.out.tsv -o ${sample_info_map.sample_label}_cazy.png
  
  """
}

