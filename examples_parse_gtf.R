source("functions_parse_gtf.R")

#load a gtf, convert to GRanges, isolate TSS, extend symetrically
gene_df = parse_gtf("examples_data/tiny.gtf")
gene_gr = gtf.df2gr(ref_df = gene_df)
gene_gr_tss = gtf.to_tss(gene_gr)
gene_gr_tss2kb = gtf.extend(ref_gr = gene_gr_tss, ext = 1000)

#the equivalent with magrittr
require(magrittr)
gene_gr_tss2kb.m = parse_gtf("examples_data/tiny.gtf") %>% 
  gtf.df2gr() %>% 
  gtf.to_tss() %>% 
  gtf.extend(ext = 1000)

#load exons from gtf, convert to GRanges, isolate TSS, extend assymetrically
exon_df = parse_gtf("examples_data/tiny.gtf", feature_type = "exon", rownames_attrib = "exon_id", additional_attrib = c("transcript_id"))

#load transcripts from gtf, convert to GRanges, extend assymetrically
trans_df = parse_gtf("examples_data/tiny.gtf")
trans_gr = gtf.df2gr(ref_df = trans_df)
trans_gr_10kup = gtf.extend_direction(trans_gr, ext = 10000, apply_upstream = T)
trans_gr_10kup_5kdown = gtf.extend_direction(trans_gr_10kup, ext = 5000, apply_upstream = F)

#the equivalent with magrittr
require(magrittr)
trans_gr_10kup_5kdown.m = parse_gtf("examples_data/tiny.gtf") %>% 
  gtf.df2gr() %>% 
  gtf.extend_direction(ext = 10^4, apply_upstream = T) %>% 
  gtf.extend_direction(ext = .5*10^4, apply_upstream = F)
