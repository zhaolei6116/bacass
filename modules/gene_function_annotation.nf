nextflow.enable.dsl=2

process extract_gene_id {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/",  mode: "rellink", overwrite: true
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
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_seq2taxid.list"), path("${sample_info_map.sample_label}_sorted_species_count.tsv"), path("${sample_info_map.sample_label}_top5.tsv"),  emit:sample_info_tuple
    path("${sample_info_map.sample_label}_nr_blast.xls")
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  """
  # 比对
  ${params.software.diamond}  blastp --db ${params.database.nr} -q ${sample_info_map.cds_faa}  -e 1e-5  --query-cover 50 --subject-cover 30  --outfmt  6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids   -o  ${sample_info_map.sample_label}_nr_blast.xls
  # 获取taxid
  cut -f 13 ${sample_info_map.sample_label}_nr_blast.xls| tr ';' '\n' | sort | uniq -c |  awk '$2 ~ /^[0-9]+$/ {print $2 "\t" $1}' | sort -k2,2 -nr | taxonkit lineage  -n  -i 2 |sed 1i'#Taxid\tCount\tLineage\tSpecies' > ${sample_info_map.sample_label}_seq2taxid.list
  # 得到物种和计数
  sh  ${projectDir}/bin/gene_annotation/tax_sort_count.sh  ${sample_info_map.sample_label}_seq2taxid.list  ${sample_info_map.sample_label}_sorted_species_count.tsv  ${sample_info_map.sample_label}_top5.tsv

  # cut -f2 seq2taxid.list |sort|uniq -c|sed 's#^\s*##g'|awk '{print \$2"\t"\$1}'|taxonkit lineage -i 1 -L -n| taxonkit reformat -i 3 -f "k__{k};K__{K};p__{p};c__{c};o__{o};f__{f};g__{g};s__{s}" -F |sed 1i'Taxid\tCount\tSpecies\tLineage' > ${sample_info_map.sample_label}_lineage.info.tsv
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
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_go_anno.out.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  python3  ${projectDir}/bin/gene_annotation/get_go_info.py   -i ${params.database.go_term}  -o ${sample_info_map.sample_label}_go_anno.out.tsv  -g ${sample_info_map.go_list}
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
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_cog_category.tsv"), path("${sample_info_map.sample_label}_category_counts.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  """
  python ${projectDir}/bin/gene_annotation/cog_info_extractor.py  -d  ${params.database.cog_def}  -f ${params.database.cog_fun}  -i  ${sample_info_map.cog_list}   -o1 ${sample_info_map.sample_label}_cog_category.tsv  -o2 ${sample_info_map.sample_label}_category_counts.tsv
  """

}

process kegg_anno {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Gene_annotation/KEGG/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_kegg_pathway_mapping.tsv"), path("${sample_info_map.sample_label}_kegg_level2_counts.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  python ${projectDir}/bin/gene_annotation/get_kegg_info.py   -i ${sample_info_map.kegg_list}  -o1 ${sample_info_map.sample_label}_kegg_pathway_mapping.tsv -o2 ${sample_info_map.sample_label}_kegg_level2_counts.tsv -o3 ${sample_info_map.sample_label}_pathway_data.json
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
  python ${projectDir}/bin/gene_annotation/parse_blast_swissprot.py  -b  ${sample_info_map.sample_label}_swiss-prot_blast.xls  -s  ${params.database.swiss_prot_id_description} -o  ${sample_info_map.sample_label}_swissprot.out.tsv
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
  python  ${projectDir}/bin/gene_annotation/get_card_info.py  -b ${sample_info_map.sample_label}_card_blast.xls  -c ${params.database.card_aro_json}  -o ${sample_info_map.sample_label}_card_out.tsv

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
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_phi_info.tsv"), path("${sample_info_map.sample_label}_phi_out.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
   
  script:
  """
  # 比对
  ${params.software.diamond}  blastp --db ${params.database.phi}  -q ${sample_info_map.cds_faa}  -e 1e-5  --query-cover 50 --subject-cover 30  -o  ${sample_info_map.sample_label}_phi_blast.xls --max-target-seqs 2
  
  # 信息获取
  python ${projectDir}/bin/gene_annotation/get_phi.py  -p  ${params.database.phi_base_info} -b ${sample_info_map.sample_label}_phi_blast.xls  -o1 ${sample_info_map.sample_label}_phi_info.tsv  -o2 ${sample_info_map.sample_label}_phi_out.tsv

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
  ${params.software.diamond}  blastp --db ${params.database.vfdb}  -q ${sample_info_map.cds_faa}  -e 1e-5  --query-cover 50 --subject-cover 30  -o  ${sample_info_map.sample_label}_vfdb_blast.xls --max-target-seqs 2
  
  # 信息获取
  python ${projectDir}/bin/gene_annotation/annotate_vfdb.py  --setb  ${params.database.vfdb_setb} --vfs ${params.database.vfdb_vfs} --blast ${sample_info_map.sample_label}_vfdb_blast.xls  --output ${sample_info_map.sample_label}_vfdb_out.tsv
  
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
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_overview.txt"), path("${sample_info_map.sample_label}_cazy.out.tsv"),   emit:sample_info_tuple
    

  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  
  script:
  """
  # 比对
  run_dbcan ${sample_info_map.cds_faa}  protein --out_dir ./ --dia_cpu 32  --tools diamond  --out_pre ${sample_info_map.sample_label}_  --db_dir ${params.database.cazy_dbcan}

  # 信息获取
  python ${projectDir}/bin/gene_annotation/cazy_class_count.py  --tsv ${sample_info_map.sample_label}_overview.txt  --json ${params.database.cazy_level1_json}  --output ${sample_info_map.sample_label}_cazy.out.tsv
  
  """
}

