#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {extract_gene_id;
         nr_anno;
         go_anno; 
         cog_anno; 
         kegg_anno; 
         swissprot_anno; 
         card_anno;
         phi_anno;
         vfdb_anno; 
         cazy_anno;
         } from "../modules/gene_function_annotation.nf"


workflow gene_anno {
  take:
    sample_info_map
  main:
    extract_gene_id(sample_info_map)
    extract_gene_id_out_ch = extract_gene_id.out.sample_info_tuple.map { samples_info, go_list, kegg_list, cog_list ->
                            samples_info + ["go_list": go_list, "kegg_list": kegg_list, "cog_list": cog_list]
                          }
    
    nr_anno(extract_gene_id_out_ch)
    go_anno(extract_gene_id_out_ch)
    cog_anno(extract_gene_id_out_ch)
    kegg_anno(extract_gene_id_out_ch)
    swissprot_anno(extract_gene_id_out_ch)
    card_anno(extract_gene_id_out_ch)
    phi_anno(extract_gene_id_out_ch)
    vfdb_anno(extract_gene_id_out_ch)
    cazy_anno(extract_gene_id_out_ch)

    // nr_anno.out.sample_info_tuple.view()
    
    gene_anno_out = nr_anno.out.sample_info_tuple
         .join(go_anno.out.sample_info_tuple)
         .join(cog_anno.out.sample_info_tuple)
         .join(kegg_anno.out.sample_info_tuple)
         .join(swissprot_anno.out.sample_info_tuple)
         .join(card_anno.out.sample_info_tuple)
         .join(phi_anno.out.sample_info_tuple)
         .join(vfdb_anno.out.sample_info_tuple)
         .join(cazy_anno.out.sample_info_tuple)
         
  emit:
    gene_anno_out 
}