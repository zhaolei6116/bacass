nextflow.enable.dsl=2

process kraken2 {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_name}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.kraken.report"),   emit:sample_info_tuple
    path("*")

  script:  
  """
  ${params.software.kraken2} \
    -db ${params.database.kraken2_db} \
    --threads 16 \
    --confidence 0.1 \
    --classified-out  ${sample_info_map.sample_name}_classified.fastq  \
    --unclassified-out  ${sample_info_map.sample_name}_unclassified.fastq \
    --output ${sample_info_map.sample_name}.kraken.out  \
	  --report ${sample_info_map.sample_name}.kraken.report \
    ${sample_info_map.clean_filt_data}
	  
  """

}


process braken {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_name}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.bracken.out"),   emit:sample_info_tuple
    path("*")
  
  script:
  """
  ${params.software.braken} \
    -d ${params.database.kraken2_db} \
    -i ${sample_info_map.kraken_report} \
    -o ${sample_info_map.sample_name}.bracken.out \
    -w ${sample_info_map.sample_name}.bracken.new.report \
    -r 300 \
    -l S \
    -t 16
  """
}

process top_10 {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_name}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("top10.png"),   emit:sample_info_tuple
    path("*")
  
  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  script:
  """
  ${params.software.csvtk}  sort -t -k 7:nr  ${sample_info_map.braken_out}|${params.software.csvtk} head|${params.software.csvtk} cut -t -f 1,6 >  top10.tsv
  python  top10_plot.py  top10.tsv
  """
  
}

process report {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_name}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("top10.png"),   emit:sample_info_tuple
    path("*")

  
  script:
  """
  source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate
  mkdir images
  mkdir tables

  ln -fs ${sample_info_map.clean_qc}  tables/qc_stat.xls
  ln -fs ${sample_info_map.clean_data_png}  images/ont_reads_length.png
  ln -fs ${sample_info_map.top10_png} images/top10.png
  ln -fs ${sample_info_map.genome_stat} tables/assemble_stat_table.tsv
  ln -fs ${sample_info_map.contig_stat}  tables/contig_stat.xls
  ln -fs ${sample_info_map.circos_png}  images/genome_circos.png
  ln -fs ${sample_info_map.coverage_png}  images/depth_coverage.png

  get_summary.py  -p  ${params.csv}  -s  tables/qc_stat.xls   -a tables/assemble_stat_table.tsv  -o tables/pj_summary.xls
  generate_report.py  -p ${sample_info_map.sample_name}_
  """

}


process check_and_modify {
    input:
    val input_value

    output:
    val input_value

    exec:
    // 使用 Groovy 判断和处理输入
    def modified_value = input_value
    if (input_value =~ /[Mm]$/) {
        // 如果以 M/m 结尾，保留原值
        modified_value = input_value.toUpperCase()
    } else {
        // 提取数字部分并除以 1000000，加上 M
        def number = input_value.replaceAll(/[^\d]/, '') as BigDecimal
        modified_value = (number / 1000000).setScale(6, BigDecimal.ROUND_HALF_UP) + 'M'
    }
    println(modified_value)
    input_value = modified_value

    // 将结果传递到 Bash 脚本
    """
    echo "${modified_value}" > result.txt
    """
}

process splitRegions {
    // split the bam reference sequences into overlapping sub-regions
    // label "medaka"
    memory "4 GB"
    errorStrategy  { return 'retry'}
    maxRetries 2
    tag "${sample_info_map.sample_name}.Try${task.attempt}"
    
    input:
        val(sample_info_map)
    
    output:
        val(sample_info_map),  emit:sample_info
        path("output.txt"),    emit:regions_txt  
    
    beforeScript "source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka"
    """
    #!/usr/bin/env python

    import itertools
    import medaka.common

    regions = itertools.chain.from_iterable(
        x.split(${params.chunk_size}, overlap=1000, fixed_size=False)
        for x in medaka.common.get_bam_regions("${sample_info_map.reads2ref_bam}"))
    
    with open("output.txt", "w") as outfile:
        for reg in regions:
            # don't ask...just grep &split!
            outfile.write("${sample_info_map.sample_id}" + '&split!' + str(reg) + "\\n")
    """
}



workflow {
    input_ch = Channel.from("100", "150M", "300m", "250", "5000000", "123456789")

    input_ch
        | check_and_modify
        | view
}
