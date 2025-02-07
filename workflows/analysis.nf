#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { nanostat; 
          fastcat; 
          ont_barcoder;
          fastplong;
          nanoplot;
          filtlong;
          cut_data;
          qc_stat;
          kraken2;
          bracken;
          top_10;
          flye;
          reorder_and_summarize;
          alignReads;
          splitRegions;
          medakaHdf;
          medakaConsensus;
          bakata;
          re_align;
          assembly_stat;
          get_report;
          get_release;
} from "../modules/software_analysis.nf"


workflow raw {
  take:
    samples_map
  
  main:
    fastcat(samples_map)
    fastcat_out_ch = fastcat.out.sample_info_tuple

    cut_data(fastcat_out_ch)
    cut_data_out_ch = cut_data.out.sample_info_tuple
                    .map {samples_info, raw_data ->
                          samples_info + ["raw_data" : raw_data.toString()]
                        }

    nanostat(cut_data_out_ch)
    nanostat_out_ch = nanostat.out.sample_info_tuple
                        .map {samples_info, raw_qc  ->
                          samples_info + ["raw_qc":raw_qc.toString()]
                        }
                        

    // fastcat(nanostat_out_ch)
    // fastcat_out_ch = fastcat.out.sample_info_tuple
    //                     .map {samples_info, raw_data ->
    //                       samples_info + ["raw_data" : raw_data.toString()]
    //                     }
  emit:
    raw_out_ch = nanostat_out_ch

}


workflow clean {
  take:
    sample_info_map
  
  main:
    // ont_barcoder(sample_info_map)
    // ont_barcoder_out_ch = ont_barcoder.out.sample_info_tuple
    //                         .map {samples_info, clean_data ->
    //                       samples_info + ["clean_data":clean_data.toString()]
    //                     }
    
    fastplong(sample_info_map)
    fastplong_out_ch = fastplong.out.sample_info_tuple
                           
    
    // filtlong(ont_barcoder_out_ch)
    filtlong(fastplong_out_ch)
    filtlong_out_ch = filtlong.out.sample_info_tuple
                          .map { samples_info, clean_filt_data ->
                            samples_info + ["clean_filt_data": clean_filt_data.toString()]
                          }
    
    nanoplot(filtlong_out_ch)
    nanoplot_out_ch = nanoplot.out.sample_info_tuple
                          .map { samples_info, clean_qc, clean_data_png ->
                            samples_info + ["clean_qc" : clean_qc, "clean_data_png":clean_data_png]
                          }
    
    qc_stat(nanoplot_out_ch)
    qc_stat_out_ch = qc_stat.out.sample_info_tuple
                        .map { samples_info, qc_stat  ->
                          samples_info + ["qc_stat" : qc_stat.toString()]
                        }
  
  emit:
    clean_out_ch = qc_stat_out_ch

}


workflow decontamination {
  take:
    sample_info_map
  
  main:
    kraken2(sample_info_map)
    kraken_out_ch = kraken2.out.sample_info_tuple
                      .map { samples_info, kraken_report  ->
                        samples_info + ["kraken_report" : kraken_report]
                      }
    
    bracken(kraken_out_ch)
    bracken_out_ch = bracken.out.sample_info_tuple
                      .map { samples_info, bracken_out ->
                        samples_info + ["bracken_out" : bracken_out]
                      }
    
    top_10(bracken_out_ch)
    top10_out_ch = top_10.out.sample_info_tuple
                      .map { samples_info, top10_png ->
                        samples_info + ["top10_png" : top10_png]
                      }


  emit:
    decontamination_out_ch = top10_out_ch
}


workflow denovo {
  take:
    sample_info_map
  
  main:
    flye(sample_info_map)
    // flye_out_ch = flye.out.sample_info_tuple
    //                 .map { samples_info, draft_assembly_fa -> 
    //                   samples_info + ["draft_assembly_fa" : draft_assembly_fa]
                    // }
    flye_out_ch = flye.out.sample_info_tuple
    
    reorder_and_summarize(flye_out_ch)
    reorder_and_summarize_out_ch = reorder_and_summarize.out.sample_info_tuple
                         .map { samples_info, flye_fa, flye_stat, summaryFile ->
                         def lines = summaryFile.text.readLines()
                         def longest = lines[0].split(':')[1].trim()
                         def chunk   = lines[1].split(':')[1].trim().toInteger()

                         return samples_info + ["flye_fa":flye_fa, "flye_stat":flye_stat, "longest_contig":longest, "chunk":chunk]

                         }


    
    
    alignReads(reorder_and_summarize_out_ch)
    alignReads_out_ch = alignReads.out.sample_info_tuple
                        .map { samples_info, reads2ref_bam, reads2ref_bam_bai ->
                          samples_info + ["reads2ref_bam": reads2ref_bam, "reads2ref_bam_bai": reads2ref_bam_bai]
                        }
    
    splitRegions(alignReads_out_ch)

    // sample_info_ch = splitRegions.out.sample_info
    //                   .map { samples_info ->
    //                   return tuple(samples_info.sample_id, samples_info)
    // }
    sample_info_and_regin_ch = splitRegions.out.sample_info
                              .map { samples_info, regin_txt ->
                                      def lines = regin_txt.text.readLines()
                                      return tuple(samples_info, lines)
                              }
                              .flatMap { samples_info, lines ->
                                      lines.collect { line -> 
                                      tuple(groupKey(samples_info, lines.size()), line)
                                      }
                               }

    // regions_ch = splitRegions.out.regions_txt.splitText()
    //               .map { 
    //                 it -> return tuple(it.split(/&split!/)[0], it.split(/&split!/)[1].trim())
    //               }
    
    // splitRegions_out_ch = sample_info_ch.combine(regions_ch, by:0).map{it[1..-1]}  // [1:-1], 去掉第一个值，即sample name，目前channel = [sample_info_map, region]

    // medakaHdf(splitRegions_out_ch, "consensus")
    medakaHdf(sample_info_and_regin_ch, "consensus")
    medakaHdf_out_ch = medakaHdf.out.info_hdf_tuple.groupTuple().view()   //队列中每个channel [sample_info_map, region1, region2 ……]

    
    medakaConsensus_out_ch = medakaConsensus(medakaHdf_out_ch)
    
    re_align(medakaConsensus_out_ch.sample_info_tuple)
    assembly_stat(re_align.out.sample_info_tuple)
    assembly_stat_out_ch =  assembly_stat.out.sample_info_tuple
                         .map {samples_info, medaka_consence, consence_bam, genome_stat, contig_stat, depth_png -> 
                         samples_info + ["medaka_consence":medaka_consence , "consence_bam":consence_bam, "genome_stat":genome_stat, "contig_stat":contig_stat, "depth_png":depth_png]
                         }

    
    bakata(assembly_stat_out_ch)
    bakata_out_ch = bakata.out.sample_info_tuple
                      .map { samples_info, circle_png ->
                      samples_info + ["circos_png":circle_png]
                      }
    
  emit:
    denovo_out_ch = bakata_out_ch
    // denovo_out_ch = splitRegions_out_ch

}

workflow report {
  take:
    sample_info_map
  
  main:
    get_report(sample_info_map)
    get_report_out_ch = get_report.out.sample_info_tuple
                          .map {samples_info, report_html ->
                          samples_info + ["report_html": report_html]
                          }

  emit:
    report_out = get_report_out_ch
}

workflow release {
  take:
    sample_info_map

  main:
    get_release(sample_info_map)

  emit:
    get_release.out
}
