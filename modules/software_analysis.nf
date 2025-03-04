nextflow.enable.dsl=2

process nanostat1 {
  
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/01_rawdata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}_raw_stat.tsv"),  emit: sample_info_tuple
      

  script:
  def fastq_path = sample_info_map.barcode ? "${sample_info_map.fastq_path}/${sample_info_map.barcode}/" : "${sample_info_map.fastq_path}/"
    
    """
    ${params.software.nanostat} --fastq  $fastq_path/*fastq.gz  --tsv -t 5 -n ${sample_info_map.sample_name}_raw_stat.tsv
    """

}

process fastcat {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/01_rawdata/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("fastq_chunks/${sample_info_map.sample_name}.raw.fastq.gz"),   emit:sample_info_tuple
    path("*")

  script:
  def fastq_path = sample_info_map.barcode ? "${sample_info_map.fastq_path}/${sample_info_map.barcode}/" : "${sample_info_map.fastq_path}/"


  """
  source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/software/fastcat/v0.18.6/
  # 测试1
  mkdir fastcat_stats
  mkdir fastq_chunks

  fastcat \
  -a 300 -q 10 \
	-s ${sample_info_map.sample_name} \
	-H  \
	-f fastcat_stats/per-file-stats.tsv \
	-r fastcat_stats/per-reads-stats.tsv \
	$fastq_path | bgzip > fastq_chunks/${sample_info_map.sample_name}.raw.fastq.gz
  """

}

process filtlong1 {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/01_rawdata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    tuple val(sample_info_map), path(rawdata)
  
  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.cut_raw.fastq.gz"), emit:sample_info_tuple

  script:
    """
    ${params.software.filtlong} \
      --min_length 3000 \
      --keep_percent 85 \
      --min_mean_q 15 \
      ${rawdata} |  gzip > ${sample_info_map.sample_name}.cut_raw.fastq.gz
    """

}

process cut_data {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/01_rawdata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"
  
  input:
    tuple val(sample_info_map), path(rawdata)
  
  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.cut_raw.fastq.gz"), emit:sample_info_tuple

  script:
  """
  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/rasusa reads \
       -o ${sample_info_map.sample_name}.cut_raw.fastq.gz  \
       -b 1.5G  \
       -O g  \
       ${rawdata}
  """
}

process nanostat2 {
  
  publishDir path:"${params.out_dir}/${sample_info.sample_id}/01_rawdata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info.sample_name}.Try${task.attempt}" 

  input:
    val sample_info

  output:
    tuple val(sample_info), path("${sample_info.sample_name}_cut_raw_stat.tsv"),  emit: sample_info_tuple
      

  script:
    
    """
    ${params.software.nanostat} --fastq  ${sample_info.raw_data}  --tsv -t 5 -n ${sample_info.sample_name}_cut_raw_stat.tsv
    """

}


  
process ont_barcoder {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.clean.fastq.gz"),   emit:sample_info_tuple
    
  
  
  script:

    def file_path = "${sample_info_map.raw_data}"
    def folder_path = new File(file_path).getParent()

  """
  ${params.software.ont_barcoder} \
    --input_path ${folder_path} \
	  --save_path  clean_fastq \
	  --barcode_kits SQK-NBD114-96 \
	  --detect_mid_strand_barcodes \
	  --compress_fastq \
	  --enable_trim_barcodes


	${params.software.zcat}   \$(find ./clean_fastq/barcode*  -name "*.fastq.gz")| bgzip > ${sample_info_map.sample_name}.clean.fastq.gz
  """

}

process fastplong {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.clean.fastq.gz"),  emit:sample_info_tuple
    
  script:
  """
  ${params.software.fastplong}  -i  ${sample_info_map.raw_data}   -o  ${sample_info_map.sample_name}.clean.fastq.gz  --failed_out  fastplong.failed.fastq.gz   -h  fastplong.html  -j  fastplong.json  -w 4
  """
    
}



process filtlong {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    tuple val(sample_info_map), path(cut_clean_data)
  
  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}.clean.filt.fastq.gz"), emit:sample_info_tuple

  script:
    """
    ${params.software.filtlong} \
      --min_length 1000 \
      --keep_percent 95 \
      ${cut_clean_data} |  gzip > ${sample_info_map.sample_name}.clean.filt.fastq.gz
    """

}

process nanoplot {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
  tuple val(sample_info_map), path("${sample_info_map.sample_name}_clean_stat.tsv"), path("nanoplot/${sample_info_map.sample_name}_LengthvsQualityScatterPlot_dot.png"), emit: sample_info_tuple
  path("nanoplot") 
  

  beforeScript 'source /data/software/miniconda/bin/activate'
  
  script:
    
  """
  ${params.software.nanoplot} \
  -t 10 \
	--fastq ${sample_info_map.clean_filt_data} \
	--maxlength 9000000000 \
	--plots dot \
	-c pink \
	-f png  \
	--tsv_stats  \
	-o nanoplot \
	-p ${sample_info_map.sample_name}_ \
	--N50

  mv nanoplot/${sample_info_map.sample_name}_NanoStats.txt  ./${sample_info_map.sample_name}_clean_stat.tsv

  """
}

process qc_stat {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_id}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}" 
  
  input:
    val sample_info_map

  output:
  tuple val(sample_info_map), path("${sample_info_map.sample_name}_qc_stat.xls"), emit: sample_info_tuple
  
  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'


  script:
    
  """  
  get_qc_stat.py -raw ${sample_info_map.raw_qc}  -clean  ${sample_info_map.clean_qc}   -o ${sample_info_map.sample_name}_qc_stat.xls
  """
}

process kraken2 {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/03_Decontamination/",  mode: "rellink", overwrite: true
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

process bracken {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/03_Decontamination/",  mode: "rellink", overwrite: true
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
  ${params.software.bracken} \
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
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/03_Decontamination/",  mode: "rellink", overwrite: true
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
  # source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate
  ${params.software.csvtk}  sort -t -k 7:nr  ${sample_info_map.bracken_out}|${params.software.csvtk} head|${params.software.csvtk} cut -t -f 1,6 >  top10.tsv
  top10_plot.py  top10.tsv
  """
  
}

process flye {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("flye_out/assembly.fasta"), path("flye_out/assembly_info.txt"), emit:sample_info_tuple
  
  
  beforeScript "source /data/software/miniconda/bin/activate"
  

  script:

  // def genome_size = sample_info_map.genome_size
  // if (sample_info_map.genome_size =~ [Mm]$/) {
  //     genome_size = sample_info_map.genome_size.toUpperCase()
  // } else {
  //     def number = sample_info_map.replaceAll(/[^\d]/, '') as BigDecimal
  //     genome_size = (number / 1000000).setScale(3, BigDecimal.ROUND_HALF_UP) + "M"
  // }
  // println(genome_size)


  // 使用 Groovy 判断和处理输入
  def genome_size = sample_info_map.genome_size
  def gz = "--genome-size ${genome_size}"
  def rd = "-g ${genome_size}"

  if (sample_info_map.genome_size =~ /-/) {
      gz = ""
      rd = "-b 1G "
  } else if (sample_info_map.genome_size =~ /[Mm]$/) {
      // 如果以 M/m 结尾，保留原值
      genome_size = sample_info_map.genome_size.toUpperCase()
      gz = "--genome-size ${genome_size}"
      rd = "-g ${genome_size} -c 250 "
  } else {
      // 提取数字部分并除以 1000000，加上 M
      def number = sample_info_map.genome_size.replaceAll(/[^\d]/, '') as BigDecimal
      genome_size = (number / 1000000).setScale(6, BigDecimal.ROUND_HALF_UP) + 'M'
      gz = "--genome-size ${genome_size}"
      rd = "-g ${genome_size} -c 250 "
  }
  

  """
  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/rasusa reads \
       -o downsample.fastq.gz  \
       ${rd}  \
       -O g  \
       ${sample_info_map.clean_filt_data}

  flye --nano-hq  downsample.fastq.gz  \
      ${gz} \
      --out-dir flye_out \
      --threads 10 \
      --min-overlap 1000
      
  # if [[ \$? -eq 0 ]]; then
  # mv flye_out/assembly.fasta     ./${sample_info_map.sample_name}.draft_assembly.fasta
  # mv flye_out/assembly_info.txt  ./${sample_info_map.sample_name}_flye_stats.tsv
  # ${params.software.bgzip}  ${sample_info_map.sample_name}.draft_assembly.fasta
  """
}

process reorder_and_summarize {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
        tuple val(sample_info_map), path(flye_fa), path(flye_out_stat)
  output:
        tuple val(sample_info_map), path("${sample_info_map.sample_name}_sorted_assembly.fa"), path("${sample_info_map.sample_name}_updated_flye_stat.tsv"), path("${sample_info_map.sample_name}_assembly_summary.txt"), emit:sample_info_tuple
  
  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'

  script:
  """
  reorder_and_summarize_contigs.py   --fasta_file  ${flye_fa}  --stats_file  ${flye_out_stat}  --chunk_size ${params.chunk_size} --prefix ${sample_info_map.sample_name}
  """
    
}

process alignReads {
    // label "wfbacterialgenomes"
    publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/",  mode: "rellink", overwrite: true
    memory "8 GB"
    errorStrategy  { return 'retry'}
    maxRetries 2
    tag "${sample_info_map.sample_name}.Try${task.attempt}"


    input:
        val(sample_info_map)
    output:
        tuple val(sample_info_map), path("${sample_info_map.sample_name}.reads2ref.bam"), path("${sample_info_map.sample_name}.reads2ref.bam.bai"), emit:sample_info_tuple
    
    
    """
    source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
    mini_align -i ${sample_info_map.clean_filt_data}  -r ${sample_info_map.flye_fa}   -p "${sample_info_map.sample_name}.reads2ref" -t 16 -m
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
        tuple val(sample_info_map), path("output.txt"),  emit:sample_info
         
    
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
            outfile.write(str(reg) + "\\n")
    """
}

process medakaHdf {
    // run medaka inference for an individual region
    // label "medaka"
    // cpus 2
    // medaka rarely uses more than 8 GB, but sometimes it does happen
    // memory { task.attempt == 1 ? "8 GB" : "15 GB" }

    errorStrategy { task.exitStatus == 137 ? "retry" : "terminate" }
    maxRetries 2
    tag "${sample_info_map.sample_name}.Try${task.attempt}"

    input:
        tuple val(sample_info_map), val(region)            
        val type

    output:
        tuple val(sample_info_map), path("*consensus_probs.hdf"),  emit: info_hdf_tuple

    script:
    assert type in ["consensus", "variant"]

    """
    #test1
    source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka    
    medaka --version
    
    medaka inference  ${sample_info_map.reads2ref_bam}   ${sample_info_map.sample_name}.consensus_probs.hdf  \
        --threads 4 --regions ${region} --model dna_r10.4.1_e8.2_400bps_sup@v4.3.0:${type}
    """

}


process medakaConsensus {
    publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/",  mode: "rellink", overwrite: true
    label "medaka"
    errorStrategy { task.exitStatus == 137 ? "retry" : "terminate" }
    maxRetries 2
    tag "${sample_info_map.sample_name}.Try${task.attempt}"
    
    input:
        tuple val(sample_info_map),path("consensus_probs*.hdf")
            
    output:
        tuple val(sample_info_map), path("${sample_info_map.sample_name}.medaka.order.fasta"),  emit: sample_info_tuple
    
    shell:
    """
    source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
    medaka sequence --threads 4  consensus_probs*.hdf   ${sample_info_map.flye_fa}    ${sample_info_map.sample_name}.medaka.fasta
    add_model_to_fasta.sh dna_r10.4.1_e8.2_400bps_sup@v4.3.0  "${sample_info_map.sample_name}.medaka.fasta"
    /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/python /nas02/project/zhaolei/pipeline/bacteria_genome_assembly/bin/fa_order.py --fasta_file ${sample_info_map.sample_name}.medaka.fasta --output_fasta ${sample_info_map.sample_name}.medaka.order.fasta
    """
}

process medakaVariant {
    label "medaka"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta),
            path("consensus_probs*.hdf"),
            path("align.bam"),
            path("align.bam.bai"),
            path("ref.fasta.gz")
    output:
        tuple val(meta), path("${meta.alias}.medaka.vcf.gz"), emit: variants
        tuple val(meta), path("${meta.alias}.variants.stats"), emit: variant_stats
    // note: extension on ref.fasta.gz might not be accurate but shouldn't (?) cause
    //       issues. Also the first step may create an index if not already existing so
    //       the alternative reference.* will break.
    """
    medaka variant ref.fasta.gz consensus_probs*.hdf vanilla.vcf
    bcftools sort vanilla.vcf > vanilla.sorted.vcf

    medaka tools annotate \
        vanilla.sorted.vcf ref.fasta.gz align.bam "${meta.alias}.medaka.vcf"

    bgzip -i "${meta.alias}.medaka.vcf"
    bcftools stats  "${meta.alias}.medaka.vcf.gz" > "${meta.alias}.variants.stats"
    """
}


process re_align_bac {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/Realign",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
      tuple val(sample_info_map),path(medaka_consence)

  output:
      tuple val(sample_info_map),path(medaka_consence), path("${sample_info_map.sample_name}_calls_to_consence.bam"), path("${sample_info_map.sample_name}_genome_stat.tsv"), path("${sample_info_map.sample_name}_contig_stat.xls"), path("output_figures/${sample_info_map.sample_name}-contig_1.coverage.png"),  emit:sample_info_tuple
        
      path("*")

  //  beforeScript = "source /nas02/software/conda/Miniconda3/miniconda3/bin/activate" 

  script:
  """
  export JAVA_HOME=/nas02/software/jdk/jdk-11.0.1/
  export JAVA_LD_LIBRARY_PATH=${JAVA_HOME}/lib/
  export PATH=${JAVA_HOME}/bin:$PATH
  # 反比得到bam
  source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
  mini_align -i ${sample_info_map.clean_filt_data}  -r ${medaka_consence}   -p  ${sample_info_map.sample_name}_calls_to_consence  -t 16 -m
  

  source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate base
  # 统计基因组信息
  ${params.software.assembly_stats} -t ${medaka_consence} > assembly_stats.tsv
  get_report_assembly_stat.py assembly_stats.tsv  ${sample_info_map.sample_name}_genome_stat.tsv

  #gc
  ${params.software.seqtk} comp ${medaka_consence}  | cut -f1,2,3,4,5,6  > assembly_gc.list
  #cov and dep, for contig 整体的
  ${params.software.samtools}  coverage  ${sample_info_map.sample_name}_calls_to_consence.bam  -o ${sample_info_map.sample_name}_cov_dep.stat
  # contig info summary
  get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_name}_cov_dep.stat -o ${sample_info_map.sample_name}_contig_stat.xls

  ${params.software.faidx}  ${medaka_consence} -i   chromsizes > genome.size
  ${params.software.bedtools} makewindows -g genome.size  -w ${params.dep_png_bin_size}  > dep_bin.bed
  ${params.software.bedtools} coverage  -a dep_bin.bed  -b ${sample_info_map.sample_name}_calls_to_consence.bam -mean  > bed_mean_cov.txt
  plot_depth.py bed_mean_cov.txt output_figures ${sample_info_map.sample_name}
  """
}

process re_align {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/Realign",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
      tuple val(sample_info_map),path(medaka_consence)

  output:
      tuple val(sample_info_map),path(medaka_consence), path("${sample_info_map.sample_name}_calls_to_consence.bam"),   emit:sample_info_tuple
      path("*")
  //  beforeScript = "source /nas02/software/conda/Miniconda3/miniconda3/bin/activate" 

  script:
  """
  # 反比得到bam
  source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
  mini_align -i ${sample_info_map.clean_filt_data}  -r ${medaka_consence}   -p  ${sample_info_map.sample_name}_calls_to_consence  -t 16 -m
  """
}

process assembly_stat {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/Realign/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    tuple val(sample_info_map),path(medaka_consence),path(consence_bam)

  output:
    tuple val(sample_info_map),path("${sample_info_map.sample_name}_assembly.fna"), path(consence_bam), path("${sample_info_map.sample_name}_genome_stat.tsv"), path("${sample_info_map.sample_name}_contig_stat.xls"), path("output_figures/${sample_info_map.sample_name}-${sample_info_map.longest_contig}.coverage.png"),path("${sample_info_map.sample_name}_replicon.tsv"), emit:sample_info_tuple
    path("*")


  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'

  script:
  """
  #gc
  ${params.software.seqtk} comp ${medaka_consence}  | cut -f1,2,3,4,5,6  > assembly_gc.list
  
  #cov and dep, for contig 整体的
  ${params.software.samtools}  coverage  ${consence_bam}  -o ${sample_info_map.sample_name}_cov_dep.stat
  
  # contig info summary
  get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_name}_cov_dep.stat -circ ${sample_info_map.flye_stat}   -o ${sample_info_map.sample_name}_contig_stat1.xls
  contig_filter_and_sort.py --contig_stat  ${sample_info_map.sample_name}_contig_stat1.xls  --in_fna  ${medaka_consence} --out_fna ${sample_info_map.sample_name}_assembly.fna  --out_stat  ${sample_info_map.sample_name}_contig_stat.xls
  get_bakta_replicon.py -i ${sample_info_map.sample_name}_contig_stat.xls -o ${sample_info_map.sample_name}_replicon.tsv


  # 统计基因组信息
  ${params.software.assembly_stats} -t ${sample_info_map.sample_name}_assembly.fna > assembly_stats.tsv
  get_report_assembly_stat.py assembly_stats.tsv  ${sample_info_map.sample_name}_genome_stat.tsv


  ${params.software.faidx}  ${sample_info_map.sample_name}_assembly.fna  -i   chromsizes > genome.size
  ${params.software.bedtools} makewindows -g genome.size  -w ${params.dep_png_bin_size}  > dep_bin.bed
  ${params.software.bedtools} coverage  -a dep_bin.bed  -b ${consence_bam} -mean  > bed_mean_cov.txt
  plot_depth.py bed_mean_cov.txt output_figures ${sample_info_map.sample_name}
  """
}

process contig_stat {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/Realign/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    tuple val(sample_info_map),path(medaka_consence),path(consence_bam)

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}_contig_stat.xls")
    
  
  script:
  """
  #gc
  ${params.software.seqtk} comp ${medaka_consence}  | cut -f1,2,3,4,5,6  > assembly_gc.list
  #cov and dep, for contig 整体的
  ${params.software.samtools}  coverage  ${consence_bam}  -o ${sample_info_map.sample_name}_cov_dep.stat
  # contig info summary
  get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_name}_cov_dep.stat -o ${sample_info_map.sample_name}_contig_stat.xls
  """

}

process depth_stat_png {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/04_Assembly/Realign/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    tuple val(sample_info_map), path(medaka_consence),path(consence_bam)

  output:
    tuple val(sample_info_map), path("output_figures/${sample_info_map.sample_name}-contig_1.coverage.png")
    path("output_figures")
    
    
  
  script:
  """
  ${params.software.faidx}  ${medaka_consence} -i   chromsizes > genome.size
  ${params.software.bedtools} makewindows -g genome.size  -w ${params.dep_png_bin_size}  > dep_bin.bed
  ${params.software.bedtools} coverage  -a dep_bin.bed  -b ${consence_bam} -mean  > bed_mean_cov.txt
  plot_depth.py bed_mean_cov.txt output_figures ${sample_info_map.sample_name}
  """

}

process bakata {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/05_Genome_Annotion/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
      val(sample_info_map)

  output:
      tuple val(sample_info_map), path("${sample_info_map.sample_name}.png"), path("${sample_info_map.sample_name}.tsv"), emit: sample_info_tuple
      path("*")


  // beforeScript = "source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/project/zhaolei/software/conda/conda_env/bakta/"

  script:
  def species_name = sample_info_map.species_name.replace(' ', '_')

  """
  source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/project/zhaolei/software/conda/conda_env/bakta/
  /nas02/project/zhaolei/software/conda/conda_env/bakta/bin/bakta \
    -d ${params.database.bakata_db} \
    -f \
    -p ${sample_info_map.sample_name} \
    -o ./ \
    --species ${species_name} \
    --replicons ${sample_info_map.replicon} \
    --complete \
    --threads 10 \
    ${sample_info_map.medaka_consence}
  """
}

process get_report {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/06_Report/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_name}_report.html"),   emit:sample_info_tuple
    
  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  def species_name = sample_info_map.species_name.replaceAll(/^"(.*)"$/, '$1')

  """
  mkdir images
  mkdir tables

  ln -fs ${sample_info_map.qc_stat}  tables/qc_stat.xls
  ln -fs ${sample_info_map.clean_data_png}  images/ont_reads_length.png
  ln -fs ${sample_info_map.top10_png} images/top10.png
  ln -fs ${sample_info_map.genome_stat} tables/assemble_stat_table.tsv
  ln -fs ${sample_info_map.contig_stat}  tables/contig_stat.xls
  ln -fs ${sample_info_map.circos_png}  images/genome_circos.png
  ln -fs ${sample_info_map.depth_png}  images/depth_coverage.png

  get_summary.py  -p ${sample_info_map.project_id} -sp "${species_name}"  -s  tables/qc_stat.xls   -a tables/assemble_stat_table.tsv -cds ${sample_info_map.cds_tsv}   -o tables/pj_summary.xls
  generate_report.py  -p ${sample_info_map.sample_name}_
  """

}

process get_release {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_id}/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_name}.Try${task.attempt}"

  input:
    val(sample_info_map)

  output:
    path("release.sh")

  script:
  """
cat > release.sh << 'EOF'
dir=\$(pwd)
mkdir -p Release/Rawdata
# mkdir -p Release/Cleandata
mkdir -p Release/Decontamination
mkdir -p Release/Assembly
mkdir -p Release/Genome_Annotion
mkdir -p Release/Report
ln -fs  \${dir}/01_rawdata/${sample_info_map.sample_name}.cut_raw.fastq.gz               Release/Rawdata/${sample_info_map.sample_name}.raw.fastq.gz
# ln -fs  ${sample_info_map.raw_qc}                 Release/Rawdata/
# ln -fs  ${sample_info_map.clean_filt_data}        Release/Cleandata
# ln -fs  ${sample_info_map.clean_qc}               Release/Cleandata
# ln -fs  ${sample_info_map.qc_stat}                Release/Cleandata
ln -fs  \${dir}/03_Decontamination/top10.tsv        Release/Decontamination
# ln -fs    \${dir}/04_Assembly/${sample_info_map.sample_name}_updated_flye_stat.tsv             Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_name}_assembly.fna        Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_name}_genome_stat.tsv           Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_name}_contig_stat.xls           Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/output_figures                           Release/Assembly/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.tsv                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.gff3                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.gbff                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.embl                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.fna                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.ffn                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.faa                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.inference.tsv                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.hypotheticals.tsv                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.hypotheticals.faa                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.txt                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.png                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.svg                          Release/Genome_Annotion/
ln -fs    \${dir}/05_Genome_Annotion/${sample_info_map.sample_name}.json                          Release/Genome_Annotion/
ln -fs    ${sample_info_map.report_html}            Release/Report
ln -fs    ${projectDir}/bin/release_readme.txt      Release/Readme.txt

${params.software.zip}  -r  ${sample_info_map.sample_id}.zip   Release/
/usr/bin/ossutil64  cp ${sample_info_map.sample_id}.zip   oss://cwbiobj${sample_info_map.report_path}  -c ${params.database.ossconfig_bj}
echo -e "${sample_info_map.sample_id}\tseqconfirm\t${sample_info_map.report_path}\t-\t-" > ${sample_info_map.sample_id}.judge.txt
/nas02/software/jdk-22.0.1/bin/java -jar /nas02/project/fengdong/04.pipe/pipe.xml/Bin/FTGSpipeV3/10.Judge/CwbioPutDataLims.jar --config /nas02/project/fengdong/04.pipe/pipe.xml/Bin/FTGSpipeV3/10.Judge/config/config.properties.BJ  --path ${sample_info_map.sample_id}.judge.txt

EOF
  """
}



workflow  {
    in_ch = Channel.fromList([[project_id:"SD241222201239", sample_name:"1222-3", fastq_path:"/BSequenator03/24122301/no_sample/20241223_0225_P2S-02106-A_PBC02152_a58c68fb/fastq_pass/CWBar0215", barcode:"CWBar1055", species_name:"Escherichia coli",genome_size:"4M", data_volume:"500000000", barcode_type:"CW-Bar16bp-2208-1", sample_id:"B22412220003"],[project_id:"SD241222201239", sample_name:"1222-4", fastq_path:"/BSequenator03/24122301/no_sample/20241223_0225_P2S-02106-A_PBC02152_a58c68fb/fastq_pass/CWBar0216", barcode:"CWBar1056", species_name:"Escherichia coli", genome_size:"4M", data_volume:"500000000", barcode_type:"CW-Bar16bp-2208-1", sample_id:"B22412220004"]])
      
    
    // in_ch.subscribe {maplist -> 
        // println "$maplist"}
    
    // in_ch.view()
    // ch_input.view()
                
    nanostat(in_ch)
    nanostat.out.raw_qc.view()
    nanostat.out.sample_info_map.view()
    

    // nanostat.out.view()

    // fastcat(ch_input)
    // fastcat.out.outp.view()
    // nanoplot(fastcat.out.list)
    // nanoplot.out.view()

}
