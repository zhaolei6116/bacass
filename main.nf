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
         //  report2;
          release;
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

   if (params.anno_analysys){
      gene_anno(ch_denovo_result)
      def ch_report_in = gene_anno.out.map {it -> it[0]}.view()
      report(ch_report_in)|release
      
   }else {
      def ch_report_result = report(ch_denovo_result)
      release(ch_report_result)
      release.out.view()
   }
   
   
   // ch_report_result.subscribe { println "Updated Map => $it" }
}




