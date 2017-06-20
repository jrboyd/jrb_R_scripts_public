#This family of funcitons is intended to allow easy loading and manipulation of gtf files.
#Specifically, gtf files downloaded from gencode: http://www.gencodegenes.org/releases/current.html
#Many of the manipulation functions are applicable generally to data.frame or GenomicRanges objects.

require(GenomicRanges)

#load data from a gtf_file and handle parsing and filtering to create and return a data.frame
#gtf_file - the file to be read from
#rownames_attrib - the attribute to use as the rownames of output data.frame.  "gene_id" by default.
#feature_type - only load entries with the specified feature_type (gene, transcript, exon ...).  "gene" by default.
#default_attrib - exposes standard attributes if you want to override to ignore them, default default_attrib = c("gene_id", "gene_name", "gene_type")
#additional_attrib - use an array to parse other attributes. none by default, c(). 
#  a usefule example: additional_attrib = c("tag", "transcript_support_level", "level")
#  "gene_id" and "gene_name" are hardcoded to be parsed.
#  see http://www.gencodegenes.org/gencodeformat.html for attribute names and descriptions
parse_gtf = function(gtf_file, rownames_attrib = 'gene_id', feature_type = 'gene', default_attrib = c("gene_id", "gene_name", "gene_type"), additional_attrib = c()) {
  print('loading gtf contents...')
  raw_lines = read.table(gtf_file, sep = "\n", stringsAsFactors = F)[, 1]
  get_col = function(i) {
    unlist(lapply(strsplit(raw_lines, "\t"), function(x) x[[i]]))
  }
  print('filtering gtf contents...')
  types = get_col(3)
  keep = types == feature_type
  raw_lines = raw_lines[keep]
  
  attribs = get_col(9)
  all_attribs = strsplit(attribs, "; +")
  get_attrib = function(key) {
    out = lapply(all_attribs, function(x) {
      keep = grepl(key, x)
      str = "NA"
      if (sum(keep) > 0) 
        str = strsplit(x[keep], " ")[[1]][2]
      return(str)
    })
    return(unlist(out))
  }
  print('parsing gtf attributes')
  # gene_id = get_attrib("gene_id")
  # gene_name = get_attrib("gene_name")
  chrm = get_col(1)
  strand = get_col(7)
  start = as.numeric(get_col(4))
  end = as.numeric(get_col(5))
  rnames = get_attrib(rownames_attrib)
  if(sum(duplicated(rnames)) > 0){
    warnings("the rowname_attrib was not unique, using arbitraty number instead. rowname_attrib will be included as column.")
    additional_attrib = c(additional_attrib, rownames_attrib)
    rnames = 1:length(all_attribs)
  }
  ref_dict = data.frame(seqnames = chrm, start, end, strand, row.names = rnames, 
                        stringsAsFactors = F)
  for(attrib in default_attrib){
    ref_dict = cbind(ref_dict, get_attrib(attrib))
    colnames(ref_dict)[ncol(ref_dict)] = attrib
  }
  for(attrib in additional_attrib){
    ref_dict = cbind(ref_dict, get_attrib(attrib))
    colnames(ref_dict)[ncol(ref_dict)] = attrib
  }
  
  return(ref_dict)
}

#convert the supplied data.frame, as output by parse_gtf(), to a GenomicRanges object
#if this doesn't work, verify your data.frame includes the column names seqnames, start, and end.
#if it still doesn't work, try updating GenomicRanges:
# source("https://bioconductor.org/biocLite.R"); biocLite("GenomicRanges")
#ref_df - data.frame with seqnames, start, and end. output from parse_gtf() works.
gtf.df2gr = function(ref_df){
  #verify required colnames are present
  req_cn = c("seqnames", "start", "end")
  if(length(intersect(req_cn, colnames(ref_df))) != length(req_cn)){
    stop(paste("missing colnames in ref_df:", setdiff(req_cn, colnames(ref_df))))
  }
  
  ref_gr = GRanges(ref_df)
  return(ref_gr)
}

#for input GRanges object, reduce to tss using strand info
#follow gtf format 1-based indexing
#ref_gr - GRanges data with strand +/-
#tes_instead - returns downstread instead of upstream
#bed_indexing - treats ranges as being 0-based index like a bed-file
gtf.to_tss = function(ref_gr, tes_instead = F, bed_indexing = F){
  req_strands = c("+", "-")
  if(length(union(unique(as.character(strand(ref_gr))), req_strands)) != length(req_strands)){
    stop("invalid strand detected, must be limited to + and/or -.")
  }
  is_pos = as.character(strand(ref_gr)) == "+"
  if(tes_instead) is_pos = !is_pos
  index_shift = 0
  if(bed_indexing) index_shift = 1
  end(ref_gr[is_pos]) = start(ref_gr[is_pos]) + index_shift
  start(ref_gr[!is_pos]) = end(ref_gr[!is_pos]) - index_shift
  return(ref_gr)
}

#extends all regions bidirectionally by ext.
#an ext of 100 will increase total width by 200 bp
#ref_gr - a GenomicRanges object to be extended
#ext - the number of bp to extend by in either direction
gtf.extend = function(ref_gr, ext){
  start(ref_gr) = start(ref_gr) - ext
  end(ref_gr) = end(ref_gr) + ext
  return(ref_gr)
}

#extends regions sensitive to feature direction using strand info.
#an ext of 100 will increase total width by 100 bp either up (default) or downstream
#ref_gr - a GenomicRanges object to be extended
#ext - the number of bp to extend by in either direction
#apply_upstream - if not upstream, then will apply downstream
gtf.extend_direction = function(ref_gr, ext, apply_upstream = T){
  req_strands = c("+", "-")
  if(length(union(unique(as.character(strand(ref_gr))), req_strands)) != length(req_strands)){
    stop("invalid strand detected, must be limited to + and/or -.")
  }
  is_pos = as.character(strand(ref_gr)) == "+"
  if(!apply_upstream) is_pos = !is_pos
  index_shift = 0
  start(ref_gr[is_pos]) = start(ref_gr[is_pos]) - ext
  end(ref_gr[!is_pos]) = end(ref_gr[!is_pos]) + ext
  return(ref_gr)
}
  
#write a loaded gtf to bed
gtf.write_bed = function(ref_gr, bed_file, convert_index_1_to_0 = T){
  #convert to GRange if necessary
  if(is.data.frame(ref_gr)) ref_gr = gtf.df2gr(ref_gr)
  mat = cbind(as.character(seqnames(ref_gr)), 
              start(ref_gr), 
              end(ref_gr), 
              names(ref_gr), 
              rep(0, length(ref_gr)), 
              as.character(strand(ref_gr)))
  write.table(mat, file = bed_file, sep = "\t", col.names = F, row.names = F, quote = F)
}
