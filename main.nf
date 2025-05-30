#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// println params


// import subworkflows
include { ingress_flow } from "./workflows/ingress.nf"

include { raw;
          clean;
          decontamination;
          denovo;
          report;
          report2;
          release;
          release2;
} from "./workflows/analysis"

include { gene_anno } from "./workflows/gene_func_anno.nf"




// workflow {
//        def ch_samples = ingress_flow()
//        def ch_nanostat_result = nanostat(ch_samples.samples)
//        def ch_fastcat_result = fastcat(ch_nanostat_result.sample_info_tuple1)
//        def ch_ont_barcoder_result = ont_barcoder(ch_fastcat_result.sample_info_tuple2)

//       //  ch_fastcat_result.sample_info_map2.subscribe {out -> println "$out"}
      
// }

workflow {
   def ch_samples = ingress_flow().view()
   def ch_raw_result = raw(ch_samples)
   def ch_clean_result = clean(ch_raw_result)
   def ch_decontamination_result = decontamination(ch_clean_result)
   def ch_denovo_result = denovo(ch_decontamination_result)

   // 使用 branch 根据 project_name 分流
   def result_channels = ch_denovo_result.branch { it -> 
      anno_analysis: it.project_name == "细菌完成图（标准分析）"
      basic_analysis: it.project_name == "细菌完成图（基础分析）"
      unknown: true
   }
   
   result_channels.basic_analysis.view { v -> "基础分析:  $v"}
   result_channels.anno_analysis.view {v -> "标准分析: $v"}

   // 如果有进入 anno_analysis 分支的数据，则执行对应的流程
   if (result_channels.anno_analysis) {
       def ch_gene_anno = gene_anno(result_channels.anno_analysis)
       def ch_report2 = report2(ch_gene_anno)
       def ch_release2 = release2(ch_report2)
       ch_release2.view()
   }

   // 如果有进入 report 分支的数据，则执行对应的流程
   if (result_channels.basic_analysis) {
       def ch_report_result = report(result_channels.basic_analysis)
       def ch_release = release(ch_report_result)
       ch_release.view()
   }

   if (result_channels.unknown) {
       result_channels.unknown.view { "未知项目名: ${it.project_name}" }
   }
   
   
   // ch_report_result.subscribe { println "Updated Map => $it" }
}




