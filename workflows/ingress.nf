nextflow.enable.dsl=2

include { getBarcode } from "../lib/kit.nf"

/*
 * ingress.nf
 *
 * 用来判断 Nextflow 参数 (params) ，并生成样本信息的 Channel
 * 再用 DSL2 的 workflow + emit 语法导出
 */

// 这里假设你在 nextflow.config 或命令行定义了以下参数
// params.manual      = params.manual      ?: false
// params.csv         = params.csv         ?: null
// params.project_id  = params.project_id  ?: null
// params.sample_name = params.sample_name ?: null
// params.fastq_path  = params.fastq_path  ?: null
// params.barcode     = params.barcode     ?: null
// params.species_name= params.species_name?: null
// params.genome_size = params.genome_size ?: null

// 定义一个子工作流 "ingress_flow" (名字可自定义)
workflow ingress_flow {
    /*
     * 这里我们创建一个 Channel 用来发射样本信息 (samples)
     * 如果是 --manual 模式，则只创建单个样本
     * 如果给了 --csv ，则从 CSV 文件中批量读取
     * 如果两者都没有，则报错
     */

    // 声明一个变量先
    def samplesCh

    main:
    // 根据不同逻辑去生成 samplesCh
    if (params.manual) {
        // 检查必要参数
        if (!params.project_id || !params.sample_name || !params.fastq_path || !params.barcode || !params.species_name || !params.genome_size) {
            error "在手动模式下，必须提供全部单样本参数 (project_id, sample_name, fastq_path, barcode, species_name, genome_size)！"
        }

        samplesCh = Channel.of([
            project_id   : params.project_id,
            sample_name  : params.sample_name,
            fastq_path   : params.fastq_path,
            barcode      : params.barcode,
            species_name : params.species_name,
            genome_size  : params.genome_size,
            sample_id    : params.sample_id,
        ])
    }
    else if (params.csv) {
        // 批量模式：从 CSV 读取
        if (!file(params.csv).exists()) {
            error "指定的 CSV 文件不存在：${params.csv}"
        }

        samplesCh = Channel
            .fromPath(params.csv)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                 def barcode_version = row.Barcode_type
                 def barcode_id = row.Barcode
                 def barcode_seq = getBarcode(barcode_version, barcode_id)
                
                [
                    project_id   : row.Client,
                    sample_id    : row.Detect_no,
                    sample_name  : row.Sample_name,
                    sample_label : "${row.Detect_no}_${row.Sample_name}",
                    barcode_type : row.Barcode_type,
                    barcode      : row.Barcode,
                    barcode_seq  : "${barcode_seq}",
                    fastq_path   : row.Path,
                    species_name : row.Species_name,
                    genome_size  : row.Genome_size,
                    data_volume  : row.Data_volume,
                    report_path  : row.Report_path,                    
                ]
            }
    }
    else {
        error "请使用 --manual 或 --csv 参数来指定输入模式"
    }

    // 最后，用 emit: 将这个 Channel 命名为 "samples"
    // 这样在外部就可以通过 "ingress_flow.out.samples" 拿到它
    emit:
        samples = samplesCh
}
