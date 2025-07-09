nextflow.enable.dsl=2

process fastcat {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/01_rawdata/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("fastq_chunks/${sample_info_map.sample_label}.raw.fastq.gz"),   emit:sample_info_tuple
    path("*")

  script:
  def fastq_path = sample_info_map.barcode ? "${sample_info_map.fastq_path}/${sample_info_map.barcode}/" : "${sample_info_map.fastq_path}"


  """
  source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/software/fastcat/v0.18.6/
  # 测试1
  mkdir fastcat_stats
  mkdir fastq_chunks

  fastcat \
  -a 300 -q 10 \
	-s ${sample_info_map.sample_label} \
	-H  \
	-f fastcat_stats/per-file-stats.tsv \
	-r fastcat_stats/per-reads-stats.tsv \
	$fastq_path | bgzip > fastq_chunks/${sample_info_map.sample_label}.raw.fastq.gz
  """

}

process cut_data {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/01_rawdata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"
  
  input:
    tuple val(sample_info_map), path(rawdata)
  
  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}.cut_raw.fastq.gz"), emit:sample_info_tuple

  script:
  """
  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/rasusa reads \
       -o ${sample_info_map.sample_label}.cut_raw.fastq.gz  \
       -b 1.5G  \
       -O g  \
       ${rawdata}
  """
}

process nanostat {
  
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/01_rawdata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_raw_stat.tsv"),  emit: sample_info_tuple
      

  script:
    
    """
    ${params.software.nanostat} --fastq  ${sample_info_map.raw_data}  --tsv -t 5 -n ${sample_info_map.sample_label}_raw_stat.tsv
    """

}


  
process ont_barcoder {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}.clean.fastq.gz"),   emit:sample_info_tuple
    
  
  
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


	${params.software.zcat}   \$(find ./clean_fastq/barcode*  -name "*.fastq.gz")| bgzip > ${sample_info_map.sample_label}.clean.fastq.gz
  """

}

process fastplong {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}.clean.fastq.gz"),  emit:sample_info_tuple
    
  script:

   def adapter = "AAGGTTAA${sample_info_map.barcode_seq}CAGCACCT"
  """
  ${params.software.fastplong}  \
      -i  ${sample_info_map.raw_data}  \
      -o  ${sample_info_map.sample_label}.clean.fastq.gz   \
      --trimming_extension 0  \
      --start_adapter  ${adapter}  \
      --distance_threshold 0.25  \
      --length_required 50  \
      --disable_quality_filtering 

  """
    
}



process filtlong {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    tuple val(sample_info_map), path(cut_clean_data)
  
  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}.clean.filt.fastq.gz"), emit:sample_info_tuple

  script:
    """
    ${params.software.filtlong} \
      --min_length 1000 \
      --keep_percent 95 \
      ${cut_clean_data} |  gzip > ${sample_info_map.sample_label}.clean.filt.fastq.gz
    """

}

process nanoplot {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
  tuple val(sample_info_map), path("${sample_info_map.sample_label}_clean_stat.tsv"), path("nanoplot/${sample_info_map.sample_label}_LengthvsQualityScatterPlot_dot.png"), emit: sample_info_tuple
  path("nanoplot") 
  
  
  script:
    
  """
  source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate /nas02/project/zhaolei/software/conda/nanoplot
  NanoPlot \
  -t 10 \
	--fastq ${sample_info_map.clean_filt_data} \
	--maxlength 9000000000 \
	--plots dot \
	-c pink \
	-f png  \
	--tsv_stats  \
	-o nanoplot \
	-p ${sample_info_map.sample_label}_ \
	--N50

  mv nanoplot/${sample_info_map.sample_label}_NanoStats.txt  ./${sample_info_map.sample_label}_clean_stat.tsv

  """
}

process qc_stat {
  publishDir path:"${params.out_dir}/${sample_info_map.sample_label}/02_cleandata/", mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}" 
  
  input:
    val sample_info_map

  output:
  tuple val(sample_info_map), path("${sample_info_map.sample_label}_qc_stat.xls"), emit: sample_info_tuple
  
  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'


  script:
    
  """  
  ${params.software.python} ${projectDir}/bin/02_Cleandata/get_qc_stat.py -raw ${sample_info_map.raw_qc}  -clean  ${sample_info_map.clean_qc}   -o ${sample_info_map.sample_label}_qc_stat.xls
  """
}

process kraken2 {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}" 

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}.kraken.report"),   emit:sample_info_tuple
    path("*")

  script:  
  """
  ${params.software.kraken2} \
    -db ${params.database.kraken2_db} \
    --threads 16 \
    --confidence 0.1 \
    --classified-out  ${sample_info_map.sample_label}_classified.fastq  \
    --unclassified-out  ${sample_info_map.sample_label}_unclassified.fastq \
    --output ${sample_info_map.sample_label}.kraken.out  \
	  --report ${sample_info_map.sample_label}.kraken.report \
    ${sample_info_map.clean_filt_data}
	  
  """
}

process bracken {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}.bracken.out"),   emit:sample_info_tuple
    path("*")
  
  script:
  """
  ${params.software.bracken} \
    -d ${params.database.kraken2_db} \
    -i ${sample_info_map.kraken_report} \
    -o ${sample_info_map.sample_label}.bracken.out \
    -w ${sample_info_map.sample_label}.bracken.new.report \
    -r 300 \
    -l S \
    -t 16
  """
}

process top_10 {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/03_Decontamination/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

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
  ${params.software.python} ${projectDir}/bin/03_Decontamination/top10_plot.py  top10.tsv
  """
  
}

process flye {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

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
  def genome_size = ""
  def gz = ""
  def rd = ""

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
      // def number = sample_info_map.genome_size.replaceAll(/[^\d]/, '') as BigDecimal
      // genome_size = (number / 1000000).setScale(6, BigDecimal.ROUND_HALF_UP) + 'M'
      // gz = "--genome-size ${genome_size}"
      // rd = "-g ${genome_size} -c 250 "
      BigDecimal number = sample_info_map.genome_size as BigDecimal
      genome_size = number.setScale(2, BigDecimal.ROUND_HALF_UP) + 'M'
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
      --meta
      
  # if [[ \$? -eq 0 ]]; then
  # mv flye_out/assembly.fasta     ./${sample_info_map.sample_label}.draft_assembly.fasta
  # mv flye_out/assembly_info.txt  ./${sample_info_map.sample_label}_flye_stats.tsv
  # ${params.software.bgzip}  ${sample_info_map.sample_label}.draft_assembly.fasta
  """
}

process reorder_and_summarize {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
        tuple val(sample_info_map), path(flye_fa), path(flye_out_stat)
  output:
        tuple val(sample_info_map), path("${sample_info_map.sample_label}_sorted_assembly.fa"), path("${sample_info_map.sample_label}_updated_flye_stat.tsv"), path("${sample_info_map.sample_label}_assembly_summary.txt"), emit:sample_info_tuple
  
  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'

  script:
  """
  ${params.software.python} ${projectDir}/bin/04_Assembly/reorder_and_summarize_contigs.py   --fasta_file  ${flye_fa}  --stats_file  ${flye_out_stat}  --chunk_size ${params.chunk_size} --prefix ${sample_info_map.sample_label}
  """
    
}

process alignReads {
    // label "wfbacterialgenomes"
    publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/",  mode: "rellink", overwrite: true
    memory "8 GB"
    errorStrategy  { return 'retry'}
    maxRetries 2
    tag "${sample_info_map.sample_label}.Try${task.attempt}"


    input:
        val(sample_info_map)
    output:
        tuple val(sample_info_map), path("${sample_info_map.sample_label}.reads2ref.bam"), path("${sample_info_map.sample_label}.reads2ref.bam.bai"), emit:sample_info_tuple
    
    
    """
    source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
    mini_align -i ${sample_info_map.clean_filt_data}  -r ${sample_info_map.flye_fa}   -p "${sample_info_map.sample_label}.reads2ref" -t 16 -m
    """
}

process splitRegions {
    // split the bam reference sequences into overlapping sub-regions
    // label "medaka"
    memory "4 GB"
    errorStrategy  { return 'retry'}
    maxRetries 2
    tag "${sample_info_map.sample_label}.Try${task.attempt}"
    
    input:
        val(sample_info_map)
    
    output:
        tuple val(sample_info_map), path("output.txt"),  emit:sample_info_tuple
         
    
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
    tag "${sample_info_map.sample_label}.Try${task.attempt}"

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
    
    medaka inference  ${sample_info_map.reads2ref_bam}   ${sample_info_map.sample_label}.consensus_probs.hdf  \
        --threads 4 --regions ${region} --model dna_r10.4.1_e8.2_400bps_sup@v4.3.0:${type}
    """

}


process medakaConsensus {
    publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/",  mode: "rellink", overwrite: true
    label "medaka"
    errorStrategy { task.exitStatus == 137 ? "retry" : "terminate" }
    maxRetries 2
    tag "${sample_info_map.sample_label}.Try${task.attempt}"
    
    input:
        tuple val(sample_info_map),path("consensus_probs*.hdf")
            
    output:
        tuple val(sample_info_map), path("${sample_info_map.sample_label}.medaka.order.fasta"), path("${sample_info_map.sample_label}.medaka.fastq"), emit: sample_info_tuple
        
    
    shell:
    """
    source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
    medaka sequence --threads 4   --qualities  consensus_probs*.hdf   ${sample_info_map.flye_fa}    ${sample_info_map.sample_label}.medaka.fastq
    ${params.software.seqkit} fq2fa  ${sample_info_map.sample_label}.medaka.fastq  -o  ${sample_info_map.sample_label}.medaka.fasta
    ${projectDir}/bin/04_Assembly/add_model_to_fasta.sh dna_r10.4.1_e8.2_400bps_sup@v4.3.0  "${sample_info_map.sample_label}.medaka.fasta"

    ${params.software.python}  ${projectDir}/bin/04_Assembly/fa_order.py --fasta_file ${sample_info_map.sample_label}.medaka.fasta --output_fasta ${sample_info_map.sample_label}.medaka.order.fasta
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
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/Realign",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
      tuple val(sample_info_map),path(medaka_consence)

  output:
      tuple val(sample_info_map),path(medaka_consence), path("${sample_info_map.sample_label}_calls_to_consence.bam"), path("${sample_info_map.sample_label}_genome_stat.tsv"), path("${sample_info_map.sample_label}_contig_stat.xls"), path("output_figures/${sample_info_map.sample_label}-contig_1.coverage.png"),  emit:sample_info_tuple
        
      path("*")

  //  beforeScript = "source /nas02/software/conda/Miniconda3/miniconda3/bin/activate" 

  script:
  """
  export JAVA_HOME=/nas02/software/jdk/jdk-11.0.1/
  export JAVA_LD_LIBRARY_PATH=${JAVA_HOME}/lib/
  export PATH=${JAVA_HOME}/bin:$PATH
  # 反比得到bam
  source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
  mini_align -i ${sample_info_map.clean_filt_data}  -r ${medaka_consence}   -p  ${sample_info_map.sample_label}_calls_to_consence  -t 16 -m
  

  source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate base
  # 统计基因组信息
  ${params.software.assembly_stats} -t ${medaka_consence} > assembly_stats.tsv
  ${params.software.python}  ${projectDir}/bin/04_Assembly/get_report_assembly_stat.py assembly_stats.tsv  ${sample_info_map.sample_label}_genome_stat.tsv

  #gc
  ${params.software.seqtk} comp ${medaka_consence}  | cut -f1,2,3,4,5,6  > assembly_gc.list
  #cov and dep, for contig 整体的
  ${params.software.samtools}  coverage  ${sample_info_map.sample_label}_calls_to_consence.bam  -o ${sample_info_map.sample_label}_cov_dep.stat
  # contig info summary
  ${params.software.python}  ${projectDir}/bin/04_Assembly/get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_label}_cov_dep.stat -o ${sample_info_map.sample_label}_contig_stat.xls

  ${params.software.faidx}  ${medaka_consence} -i   chromsizes > genome.size
  ${params.software.bedtools} makewindows -g genome.size  -w ${params.dep_png_bin_size}  > dep_bin.bed
  ${params.software.bedtools} coverage  -a dep_bin.bed  -b ${sample_info_map.sample_label}_calls_to_consence.bam -mean  > bed_mean_cov.txt
  ${params.software.python}  ${projectDir}/bin/04_Assembly/plot_depth.py bed_mean_cov.txt output_figures ${sample_info_map.sample_label}
  """
}

process re_align {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/Realign",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
      tuple val(sample_info_map), path(medaka_consence), path(medaka_consence_fq)

  output:
      tuple val(sample_info_map),path(medaka_consence), path(medaka_consence_fq), path("${sample_info_map.sample_label}_calls_to_consence.bam"),   emit:sample_info_tuple
      path("*")
  //  beforeScript = "source /nas02/software/conda/Miniconda3/miniconda3/bin/activate" 

  script:
  """
  # 反比得到bam
  source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate  /nas02/project/zhaolei/software/conda/conda_env/medaka
  mini_align -i ${sample_info_map.clean_filt_data}  -r ${medaka_consence}   -p  ${sample_info_map.sample_label}_calls_to_consence  -t 16 -m
  """
}

process assembly_stat {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/Realign/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    tuple val(sample_info_map),path(medaka_consence), path(medaka_consence_fq), path(consence_bam)

  output:
    tuple val(sample_info_map),path("${sample_info_map.sample_label}_assembly.fna"), path(consence_bam), path("${sample_info_map.sample_label}_genome_stat.tsv"), path("${sample_info_map.sample_label}_contig_stat.xls"), path("output_figures/${sample_info_map.sample_label}-${sample_info_map.longest_contig}.coverage.png"),path("${sample_info_map.sample_label}_replicon.tsv"), emit:sample_info_tuple
    path("${sample_info_map.sample_label}_assembly.fq")


  beforeScript 'source /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'

  script:
  
  def species_name = sample_info_map.species_name.replaceAll(/^"(.*)"$/, '$1')

  """
  #gc
  ${params.software.seqtk} comp ${medaka_consence}  | cut -f1,2,3,4,5,6  > assembly_gc.list
  
  #cov and dep, for contig 整体的
  ${params.software.samtools}  coverage  ${consence_bam}  -o ${sample_info_map.sample_label}_cov_dep.stat
  
  # 获取 phylum 信息
  species_name="${species_name}" 
  phylum=\$(echo "\${species_name}" | cut -f1 -d ' ' | ${params.software.taxonkit} name2taxid  --data-dir ${params.database.taxonomy} | ${params.software.taxonkit} lineage -i 2 --data-dir ${params.database.taxonomy} | ${params.software.taxonkit} reformat -i 3 -f "p__{p}" -F  --data-dir ${params.database.taxonomy} |cut -f4)
  echo \${phylum}
  # 判断是否属于 Actinomycetota
  if [[ "\${phylum}" == "p__Actinomycetota" ]]; then
    echo "Species belongs to Actinomycetota. Skipping contig_filter_and_sort.py..."
    ${params.software.python}  ${projectDir}/bin/04_Assembly/get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_label}_cov_dep.stat -circ ${sample_info_map.flye_stat}   -o ${sample_info_map.sample_label}_contig_stat1.xls
    ${params.software.python}  ${projectDir}/bin/04_Assembly/contig_filter_and_sort_2.py --contig_stat  ${sample_info_map.sample_label}_contig_stat1.xls  --in_fna  ${medaka_consence} --out_fna ${sample_info_map.sample_label}_assembly.fna  --out_stat  ${sample_info_map.sample_label}_contig_stat.xls
    cut -f1 ${sample_info_map.sample_label}_contig_stat.xls | sed '1d'  > contig.list
    ${params.software.seqtk}  subseq  ${medaka_consence_fq}  contig.list > ${sample_info_map.sample_label}_assembly.fq
    
  else
    # 执行 contig_filter_and_sort.py
    # contig info summary
    ${params.software.python}  ${projectDir}/bin/04_Assembly/get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_label}_cov_dep.stat -circ ${sample_info_map.flye_stat}   -o ${sample_info_map.sample_label}_contig_stat1.xls
    ${params.software.python}  ${projectDir}/bin/04_Assembly/contig_filter_and_sort.py --contig_stat  ${sample_info_map.sample_label}_contig_stat1.xls  --in_fna  ${medaka_consence} --out_fna ${sample_info_map.sample_label}_assembly.fna  --out_stat  ${sample_info_map.sample_label}_contig_stat.xls
    cut -f1 ${sample_info_map.sample_label}_contig_stat.xls | sed '1d'  > contig.list
    ${params.software.seqtk}  subseq  ${medaka_consence_fq}  contig.list > ${sample_info_map.sample_label}_assembly.fq 
  fi

  ${params.software.python}  ${projectDir}/bin/04_Assembly/get_bakta_replicon.py -i ${sample_info_map.sample_label}_contig_stat.xls -o ${sample_info_map.sample_label}_replicon.tsv

  # 统计基因组信息
  ${params.software.assembly_stats} -t ${sample_info_map.sample_label}_assembly.fna > assembly_stats.tsv
  ${params.software.python}  ${projectDir}/bin/04_Assembly/get_report_assembly_stat.py assembly_stats.tsv  ${sample_info_map.sample_label}_genome_stat.tsv


  ${params.software.faidx}  ${sample_info_map.sample_label}_assembly.fna  -i   chromsizes > genome.size
  ${params.software.bedtools} makewindows -g genome.size  -w ${params.dep_png_bin_size}  > dep_bin.bed
  ${params.software.bedtools} coverage  -a dep_bin.bed  -b ${consence_bam} -mean  > bed_mean_cov.txt
  ${params.software.python}  ${projectDir}/bin/04_Assembly/plot_depth.py bed_mean_cov.txt output_figures ${sample_info_map.sample_label}
  """
}

process contig_stat {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/Realign/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    tuple val(sample_info_map),path(medaka_consence),path(consence_bam)

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_contig_stat.xls")
    
  
  script:
  """
  #gc
  ${params.software.seqtk} comp ${medaka_consence}  | cut -f1,2,3,4,5,6  > assembly_gc.list
  #cov and dep, for contig 整体的
  ${params.software.samtools}  coverage  ${consence_bam}  -o ${sample_info_map.sample_label}_cov_dep.stat
  # contig info summary
  ${params.software.python}  ${projectDir}/bin/04_Assembly/get_report_contig_stat.py -gc assembly_gc.list -cov  ${sample_info_map.sample_label}_cov_dep.stat -o ${sample_info_map.sample_label}_contig_stat.xls
  """

}

process depth_stat_png {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/04_Assembly/Realign/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 0
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    tuple val(sample_info_map), path(medaka_consence),path(consence_bam)

  output:
    tuple val(sample_info_map), path("output_figures/${sample_info_map.sample_label}-contig_1.coverage.png")
    path("output_figures")
    
    
  
  script:
  """
  ${params.software.faidx}  ${medaka_consence} -i   chromsizes > genome.size
  ${params.software.bedtools} makewindows -g genome.size  -w ${params.dep_png_bin_size}  > dep_bin.bed
  ${params.software.bedtools} coverage  -a dep_bin.bed  -b ${consence_bam} -mean  > bed_mean_cov.txt
  ${params.software.python}  ${projectDir}/bin/04_Assembly/plot_depth.py bed_mean_cov.txt output_figures ${sample_info_map.sample_label}
  """

}

process bakata {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/05_Genome_Annotation/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
      val(sample_info_map)

  output:
      tuple val(sample_info_map), path("${sample_info_map.sample_label}.png"), path("${sample_info_map.sample_label}.faa"), path("${sample_info_map.sample_label}.gff3"),path("${sample_info_map.sample_label}.tsv"), emit: sample_info_tuple
      path("*")


  // beforeScript = "source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/project/zhaolei/software/conda/conda_env/bakta/"

  script:
  def species_name = sample_info_map.species_name.replace(' ', '_')

  """
  source /nas02/software/conda/Miniconda3/miniconda3/bin/activate /nas02/project/zhaolei/software/conda/conda_env/bakta/
  /nas02/project/zhaolei/software/conda/conda_env/bakta/bin/bakta \
    -d ${params.database.bakata_db} \
    -f \
    -p ${sample_info_map.sample_label} \
    -o ./ \
    --species ${species_name} \
    --replicons ${sample_info_map.replicon} \
    --complete \
    --threads 10 \
    ${sample_info_map.medaka_consence}
  """
}

process get_report {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/06_Report/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val sample_info_map

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_report.html"),   emit:sample_info_tuple
    
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

  ${params.software.python}  ${projectDir}/bin/06_report/get_summary.py  -p ${sample_info_map.project_id} -sp "${species_name}"  -s  tables/qc_stat.xls   -a tables/assemble_stat_table.tsv -cds ${sample_info_map.cds_tsv}   -o tables/pj_summary.xls
  ${params.software.python}  ${projectDir}/bin/06_report/generate_report.py  -p ${sample_info_map.sample_label}
  """

}

process get_anno_report {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/07_Report/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    tuple( val(sample_info_map),
          path(nr_species_count_tsv),
          path(nr_top5_tsv),
          path(nr_top5_png),
          path(go_anno_tsv),
          path(go_stats_tsv),
          path(go_png),
          path(cog_category_tsv),
          path(cog_category_counts_tsv),
          path(category_barplot_png),
          path(kegg_pathway_mapping_tsv),
          path(kegg_level2_counts_tsv),
          path(kegg_plot_png),
          path(swissprot_out_tsv),
          path(card_out_tsv),
          path(phi_info_tsv),
          path(phi_Mutant_Phenotype_stat_tsv),
          path(phi_stat_tsv),
          path(phi_Mutant_Phenotype_stat_png),
          path(vfdb_out_tsv),
          path(cazy_out_tsv),
          path(cazy_png),
    )

  output:
    tuple val(sample_info_map), path("${sample_info_map.sample_label}_report.html"),   emit:sample_info_tuple
    
  beforeScript 'source  /nas02/project/zhaolei/software/conda/conda_env/bininfo/bin/activate'
  
  script:
  def species_name = sample_info_map.species_name.replaceAll(/^"(.*)"$/, '$1')

  """
  mkdir images
  mkdir tables

  ln -fs  ${sample_info_map.qc_stat}  tables/qc_stat.xls
  ln -fs  ${sample_info_map.clean_data_png}  images/ont_reads_length.png
  ln -fs  ${sample_info_map.top10_png} images/top10.png
  ln -fs  ${sample_info_map.genome_stat} tables/assemble_stat_table.tsv
  ln -fs  ${sample_info_map.contig_stat}  tables/contig_stat.xls
  ln -fs  ${sample_info_map.circos_png}  images/genome_circos.png
  ln -fs  ${sample_info_map.depth_png}  images/depth_coverage.png
  mv  ${nr_top5_tsv} tables/nr_table.tsv  
  mv  ${nr_top5_png} images/nr.png
  mv  ${go_png}  images/go.png
  mv  ${category_barplot_png} images/cog.png
  mv  ${kegg_plot_png} images/kegg.png
  head -n 11 ${swissprot_out_tsv} > tables/SwissProt_table.tsv 
  head -n 11 ${card_out_tsv} > tables/card_table.tsv
  head -n 11  ${phi_stat_tsv} > tables/phi_table.tsv
  mv ${phi_Mutant_Phenotype_stat_png} images/phi.png
  head -n 11 ${vfdb_out_tsv} > tables/vfdb_table.tsv
  head -n 11 ${cazy_out_tsv} > tables/cazy_table.tsv
  mv ${cazy_png}  images/cazy.png
 
  ${params.software.python}  ${projectDir}/bin/06_report/get_summary.py  -p ${sample_info_map.project_id} -sp "${species_name}"  -s  tables/qc_stat.xls   -a tables/assemble_stat_table.tsv -cds ${sample_info_map.cds_tsv}   -o tables/pj_summary.xls
  ${params.software.python}  ${projectDir}/bin/06_report/generate_report.py  -p ${sample_info_map.sample_label} -t ${projectDir}/bin/06_report/demo_anno_template.html  -c ${projectDir}/bin/06_report/anno_gene_report_config.ini
  """

}

process get_release {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val(sample_info_map)

  output:
    path("release.sh")

  script:
  def labMapping = ['B': 'BJ', 'W': 'WH', 'T': 'TZ', 'S': 'SH', 'G': 'GZ']
  def firstLetter = sample_info_map.sample_id?.charAt(0)
  def labCode = labMapping["${firstLetter}"]


  if (labCode == null) {
      println "警告：样本来源实验室不明确，使用默认BJ"
      labCode = "BJ"  // 或者其他默认值
  }



  """
cat > release.sh << 'EOF'
dir=\$(pwd)
mkdir -p Rawdata
# mkdir -p Release/Cleandata
mkdir -p Release/Decontamination
mkdir -p Release/Assembly
mkdir -p Release/Genome_Annotation
mkdir -p Release/Report
ln -fs  \${dir}/01_rawdata/${sample_info_map.sample_label}.cut_raw.fastq.gz               Rawdata/${sample_info_map.sample_label}.raw.fastq.gz
# ln -fs  ${sample_info_map.raw_qc}                 Release/Rawdata/
# ln -fs  ${sample_info_map.clean_filt_data}        Release/Cleandata
# ln -fs  ${sample_info_map.clean_qc}               Release/Cleandata
# ln -fs  ${sample_info_map.qc_stat}                Release/Cleandata
ln -fs  \${dir}/03_Decontamination/top10.tsv        Release/Decontamination
# ln -fs    \${dir}/04_Assembly/${sample_info_map.sample_label}_updated_flye_stat.tsv             Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_assembly.fna        Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_assembly.fq        Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_genome_stat.tsv           Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_contig_stat.xls           Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/output_figures                           Release/Assembly/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.tsv                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.gff3                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.gbff                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.embl                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.fna                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.ffn                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.faa                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.inference.tsv                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.hypotheticals.tsv                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.hypotheticals.faa                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.txt                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.png                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.svg                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.json                          Release/Genome_Annotation/
ln -fs    ${sample_info_map.report_html}            Release/Report
ln -fs    ${projectDir}/bin/release_readme.txt      Release/Readme.txt

${params.software.zip}  -r  ${sample_info_map.sample_id}.zip   Release/
${params.software.zip}  -r  ${sample_info_map.sample_id}_rawdata.zip  Rawdata/
${params.software.ossutil64}  cp -f  ${sample_info_map.sample_id}.zip   oss://cwbiobj${sample_info_map.report_path}  -c ${params.database.ossconfig_bj}
${params.software.ossutil64}  cp -f ${sample_info_map.sample_id}_rawdata.zip   oss://cwbiobj${sample_info_map.report_raw_path}  -c ${params.database.ossconfig_bj}
echo -e "${sample_info_map.sample_id}\\tseqconfirm\\t${sample_info_map.report_path}\\t-\\t-" > ${sample_info_map.sample_id}.judge.txt
${params.software.java} -jar ${projectDir}/script/resources/CwbioPutDataLims.V2.jar --config ${projectDir}/script/resources/config/config.properties.${labCode}  --path ${sample_info_map.sample_id}.judge.txt

EOF
  """
}


process get_release_2 {
  publishDir path: "${params.out_dir}/${sample_info_map.sample_label}/",  mode: "rellink", overwrite: true
  errorStrategy  { return 'retry'}
  maxRetries 2
  tag "${sample_info_map.sample_label}.Try${task.attempt}"

  input:
    val(sample_info_map)

  output:
    path("release.sh")

  script:
  def labMapping = ['B': 'BJ', 'W': 'WH', 'T': 'TZ', 'S': 'SH', 'G': 'GZ']
  def firstLetter = sample_info_map.sample_id?.charAt(0)
  def labCode = labMapping["${firstLetter}"]


  if (labCode == null) {
      println "警告：样本来源实验室不明确，使用默认BJ"
      labCode = "BJ"  // 或者其他默认值
  }



  """
cat > release.sh << 'EOF'
dir=\$(pwd)
mkdir -p Rawdata
# mkdir -p Release/Cleandata
mkdir -p Release/Decontamination
mkdir -p Release/Assembly
mkdir -p Release/Genome_Annotation
mkdir -p Release/Gene_function_annotation
mkdir -p Release/Report
ln -fs  \${dir}/01_rawdata/${sample_info_map.sample_label}.cut_raw.fastq.gz                         Rawdata/${sample_info_map.sample_label}.raw.fastq.gz
# ln -fs  ${sample_info_map.raw_qc}                                                                 Release/Rawdata/
# ln -fs  ${sample_info_map.clean_filt_data}                                                        Release/Cleandata
# ln -fs  ${sample_info_map.clean_qc}                                                               Release/Cleandata
# ln -fs  ${sample_info_map.qc_stat}                                                                Release/Cleandata
ln -fs  \${dir}/03_Decontamination/top10.tsv                                                        Release/Decontamination
# ln -fs    \${dir}/04_Assembly/${sample_info_map.sample_label}_updated_flye_stat.tsv               Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_assembly.fna                  Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_assembly.fq                   Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_genome_stat.tsv               Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/${sample_info_map.sample_label}_contig_stat.xls               Release/Assembly/
ln -fs    \${dir}/04_Assembly/Realign/output_figures                                                Release/Assembly/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.tsv                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.gff3                         Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.gbff                         Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.embl                         Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.fna                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.ffn                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.faa                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.inference.tsv                Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.hypotheticals.tsv            Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.hypotheticals.faa            Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.txt                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.png                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.svg                          Release/Genome_Annotation/
ln -fs    \${dir}/05_Genome_Annotation/${sample_info_map.sample_label}.json                         Release/Genome_Annotation/
ln -fs    \${dir}/06_Gene_annotation/*                                                              Release/Gene_function_annotation
ln -fs    ${sample_info_map.report_html}            Release/Report
ln -fs    ${projectDir}/bin/release_func_anno_readme.txt      Release/Readme.txt

${params.software.zip}  -r  ${sample_info_map.sample_id}.zip   Release/
${params.software.zip}  -r  ${sample_info_map.sample_id}_rawdata.zip  Rawdata/
${params.software.ossutil64}  cp -f  ${sample_info_map.sample_id}.zip   oss://cwbiobj${sample_info_map.report_path}  -c ${params.database.ossconfig_bj}
${params.software.ossutil64}  cp -f ${sample_info_map.sample_id}_rawdata.zip   oss://cwbiobj${sample_info_map.report_raw_path}  -c ${params.database.ossconfig_bj}
echo -e "${sample_info_map.sample_id}\\tseqconfirm\\t${sample_info_map.report_path}\\t-\\t-" > ${sample_info_map.sample_id}.judge.txt
${params.software.java} -jar ${projectDir}/script/resources/CwbioPutDataLims.V2.jar --config ${projectDir}/script/resources/config/config.properties.${labCode}  --path ${sample_info_map.sample_id}.judge.txt

EOF
  """
}
