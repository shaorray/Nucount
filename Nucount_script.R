# R script for nucleosomal reads assignment
#
# Source code: https://github.com/shaorray/Nucount
#
# MIT License
# Copyright (c) 2022 Rui Shao <shaorray@hotmail.com>

# --------------------------------------------- initialization ------------------------------------------------- #
args <- commandArgs(TRUE)

suppressMessages(require(optparse, quietly = T))

suppressMessages(require(dplyr, quietly = T))
suppressMessages(require(S4Vectors, quietly = T))

suppressMessages(require(Rsamtools, quietly = T))
suppressMessages(require(rtracklayer, quietly = T))

suppressMessages(require(GenomeInfoDb, quietly = T))
suppressMessages(require(BiocGenerics, quietly = T))
suppressMessages(require(IRanges, quietly = T))

suppressMessages(require(foreach, quietly = T))
suppressMessages(require(doParallel, quietly = T))
suppressMessages(require(ggplot2, quietly = T))

options("scipen" = 999, "digits" = 4)
options(warn = -1)
registerDoParallel(cores = 4)

# --------------------------------------------- specifications ------------------------------------------------- #

option_list = list(
  # bam input
  make_option(c("-b", "--bam"), action = "store", default = NULL, type = 'character',
              help = "Bam file input group 1"),
  make_option(c("-p", "--paired"), action = "store", default = TRUE, type = 'logical',
              help = "Paired-end reads [default %default]"),
  
  # annotations
  make_option(c("-n", "--nuc_fltd"), action = "store", default = NULL, type = 'character',
              help = "Filtered TSS flanking nucleosome GRanges in RData format."),
  make_option(c("-t", "--tss_fltd"), action = "store", default = NULL, type = 'character',
              help = "Filtered TSS GRanges in RData format."),
  
  # tuning algorithm
  make_option(c("-e", "--is_EM"), action = "store_false", default = FALSE,
              help = "Perform EM optimization of density function [default %default]"),
  make_option("--iter_num", default = 2,
              help = "Number of EM iteration [default %default]"),
  
  # output
  make_option("--nuc_num", action = "store", default = 10, type = 'integer',
              help = "Number of nucleosomes away from TSS to keep [default %default]"),
  make_option("--type", action = "store", default = "table", type = 'character',
              help = "Output file format, BW or table [default %default]"),
  make_option(c("-o", "--out"), action = "store", default = "est_Nucount", type = 'character',
              help = "Output file name"),
  
  make_option(c("-q", "--quietly"), action = "store_false", default = FALSE,
              help = "Mute verbos and messages")
)
arguments <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)
opt <- arguments$options


# parse arguments
bam_file <- opt$bam
is_paired_end <- opt$paired

nuc_file <- opt$nuc_fltd
tss_file <- opt$tss_fltd

is_EM <- opt$is_EM
iter_num <- opt$iter_num

nuc_num <- opt$nuc_num
out_name <- opt$out
out_type <- opt$type

is_queitly <- opt$quietly

if (!is_queitly) {
  message("\n----------------------------------------------\n")
  message(paste("BAM:\n", bam_file, "\n"))
  message("----------------------------------------------\n")
  message(paste("Is read paired:\n", is_paired_end, "\n"))
  message("----------------------------------------------\n")
  message(paste("Nucleosome file:\n", nuc_file, "\n"))
  message("----------------------------------------------\n")
  message(paste("TSS file:\n", tss_file, "\n"))
  message("----------------------------------------------\n")
  message(paste("Perform EM:\n", is_EM, "\n"))
  message("----------------------------------------------\n")
  message(paste("EM iteration time:\n", iter_num, "\n"))
  message("----------------------------------------------\n")
  message(paste("Output name:\n", out_name, "\n"))
  message("----------------------------------------------\n")
  message(paste("Output type:\n", out_type, "\n"))
  message("----------------------------------------------\n")
}

if (!file.exists("Nucount_output"))
  dir.create("Nucount_output")

# ----------------------------------------------- functions ---------------------------------------------------- #
prob_den_fun <-  function(x) exp(-(x / 20)^2 / 280) # initialize with a normal distribution with mean = 0, sd = 20

bam_nucelsome_reads <- function(bam_file, is_paired_end = T, nuc_pos_gr, tss_gr) {
  # function (1)
  
  if (!file.exists(paste0(bam_file,'.bai'))) Rsamtools::indexBam(bam_file)
  if (!is_queitly) message("[1] BAM OK.\n")
  
  if (is_paired_end) {
    sbw = c('pos', 'qwidth','strand','rname', 'mrnm', 'mpos', 'isize')
    flg = Rsamtools::scanBamFlag(isUnmappedQuery = F,
                                 isProperPair = T,
                                 isNotPassingQualityControls = F,
                                 isFirstMateRead = T,
                                 isSecondaryAlignment = F)
  } else {
    sbw = c('pos', 'qwidth','strand','rname')
    flg = Rsamtools::scanBamFlag(isUnmappedQuery = F,
                                 isNotPassingQualityControls = F,
                                 isSecondaryAlignment = F)
  }
  param = Rsamtools::ScanBamParam(what = sbw,
                                  flag = flg,
                                  which = tss_gr)
  
  srg = Rsamtools::scanBam(bam_file, param = param)
  
  if (is_paired_end) { # discard non- or  > triple nucleosomal reads
    all_index = with(srg, !is.na(mpos) & rname == mrnm & abs(isize) > 10 & abs(isize) < 1000)
    srg <- lapply(srg, function(x) x = x[all_index])
  }
  
  if (!is_queitly)
    message(paste("Number of reads:", sum(unlist(lapply(srg, function(x) length(x[[1]])))), "\n"))
  srg # short read ranges
}

attribute_read_to_nucleosome <- function(tss_tmp, srg_tmp, tmp_nuc, is_paired_end, prob_den_fun) {
  # function (2)
  
  .strand = as.character(strand(tss_tmp))
  .start = start(tss_tmp)
  
  tmp_nuc_ranges = sort(ranges(shift(tmp_nuc, shift = -.start)))
  
  if (is_paired_end) {
    read_ranges = GRanges(seqnames = "1",
                          IRanges(start = with(srg_tmp, ifelse(.strand == "+", pos, mpos)) - .start, 
                          width = abs(srg_tmp$isize)), 
                          strand = srg_tmp$strand)
  } else {
    srg_tmp_start = with(srg_tmp, ifelse(strand == "+", pos, pos + qwidth)) - .start
    read_ranges = IRanges(start = srg_tmp_start, width = 1)
  }
  
  out_attribution = rep(NA, length(read_ranges)) # initiate read-nucleosome attribution
  
  # Step1: allocate unique read-nucleosome overlaps
  rov_cnt = countQueryHits(findOverlaps(query = read_ranges, subject = tmp_nuc_ranges))
  uni_idx = rov_cnt == 1
  uni_mtch = findOverlaps(query = read_ranges[uni_idx], subject = tmp_nuc_ranges)
  uni_mtch_ls = split(queryHits(uni_mtch), subjectHits(uni_mtch))
  
  tmp_count = countSubjectHits(uni_mtch)
  
  out_attribution[uni_idx][queryHits(uni_mtch)] = subjectHits(uni_mtch)
  
  # Step2: EM on multi-overlapped nucleosomal reads
  if (sum(!uni_idx) > 0) {
    read_ranges_multi = read_ranges[!uni_idx]
    multi_mat = cbind(start = start(read_ranges_multi), end = end(read_ranges_multi))
    
    # keep reads within 1 kb
    keep_idx = apply(multi_mat, 1, function(x) 
      pmin(min(abs(x[1] - start(tmp_nuc_ranges))),
           min(abs(x[2] - start(tmp_nuc_ranges)))) < 1000)
    
    # discard reads without any overlap of nucleosome centers
    if (is_paired_end) keep_idx = keep_idx & (rov_cnt > 1)[!uni_idx]  
    
    multi_mat = multi_mat[keep_idx, ]
    
    # sum the estimated probability from all qualified reads 
    if (is_paired_end) {
      tmp_nuc_prob = 
        sapply(seq_along(tmp_nuc_ranges), 
               function(x) {
                 x_nuc_start = start(tmp_nuc_ranges[x])
                 x_nuc_prob = sum(pmax(prob_den_fun(x - multi_mat[, 1]), 
                                       prob_den_fun(x - multi_mat[, 2])))
                 # specify for unique overlapped reads
                 if (x %in% names(uni_mtch_ls)) {
                   tmp_read_ranges = read_ranges[uni_idx][uni_mtch_ls[[as.character(x)]]]
                   x_nuc_prob = x_nuc_prob +
                     sum(pmax(prob_den_fun(start(tmp_nuc_ranges)[x] - start(tmp_read_ranges)), 
                              prob_den_fun(start(tmp_nuc_ranges)[x] - end(tmp_read_ranges))))
                 }
                 x_nuc_prob
               }
        )
    } else { # single-end reads must towards their dedicated nucleosomes
      tmp_nuc_prob = 
        sapply(seq_along(tmp_nuc_ranges), 
               function(x) {
                 x_nuc_start = start(tmp_nuc_ranges[x])
                 x_nuc_prob = 0
                 if (x %in% names(uni_mtch_ls)) {
                   tmp_read_ranges = read_ranges[uni_idx][uni_mtch_ls[[as.character(x)]]]
                   x_nuc_prob = sum(prob_den_fun(x_nuc_start - start(tmp_read_ranges)))
                 }
                 x_nuc_gaps = x_nuc_start - multi_mat[, 1]
                 x_nuc_prob = x_nuc_prob +
                   sum(prob_den_fun(x_nuc_gaps) * # filter positioning
                         xor(x_nuc_gaps > 0, srg_tmp$strand[!uni_idx][keep_idx] == "-") ) 
                 x_nuc_prob
               }
        )
    }
    
    tmp_attribution = apply(multi_mat, 1,
                            function(x)
                              which.max(tmp_nuc_prob * prob_den_fun(x - start(tmp_nuc_ranges))))
    
    out_attribution[!uni_idx][keep_idx] = tmp_attribution
  } # end of EM
  
  out_attribution
}


iterate_pdf <- function(gene_ids, tss_gr, nuc_pos_gr, srg, is_paired_end, nuc_num) {
  # function (3)
  
  read_nucleosome_coverage = foreach(gene_id = gene_ids, .combine = "c") %dopar% {
    n = which(tss_gr$gene_id == gene_id)
    srg_tmp = srg[[n]]
    if (length(srg_tmp[[1]]) == 0) return(NULL)
    tss_tmp = tss_gr[n]
    
    tmp_nuc = nuc_pos_gr[nuc_pos_gr$gene_id == tss_tmp$gene_id]
    
    .strand = as.character(strand(tss_tmp))
    .start = start(tss_tmp)
    
    tmp_attribution = attribute_read_to_nucleosome(tss_tmp, srg_tmp, tmp_nuc, is_paired_end, prob_den_fun)
    
    # relative position between reads and nucleosome
    delta_pos = with(srg_tmp, pos + ifelse(strand == "+", 0, qwidth)) - start(sort(tmp_nuc))[tmp_attribution]
    
    if (is_paired_end) 
      delta_pos = delta_pos +
      with(srg_tmp, ifelse(strand == "+", mpos + qwidth, mpos)) - start(sort(tmp_nuc))[tmp_attribution]
    
    if (.strand == "-") delta_pos = delta_pos * (-1)
    
    delta_pos[which(gsub("-", "", tmp_nuc$nuc_pos[tmp_attribution]) %in% seq_len(nuc_num))]
  } # end of foreach
  approxfun(density(read_nucleosome_coverage), yleft = 0, yright = 0)
}


# main function
get_nucelsome_count <- function(bam_file,            # bam file input
                                is_paired_end = T,   # if reads paried end
                                nuc_pos_gr,          # saved nucleosome GRanges
                                nuc_num,             # the number of nucleosomes to use, so symmetric outputs = nuc_num * 2  
                                tss_gr,              # saved TSS GRanges
                                is_EM = F,           # if update probability density function
                                iter_num = 2      # number of EM iteration
                                ) {

  tss_gr = tss_gr[tss_gr$gene_id %in% nuc_pos_gr$gene_id]
  nuc_pos_gr = nuc_pos_gr[nuc_pos_gr$gene_id %in% tss_gr$gene_id]
  
  # process bam file
  srg = bam_nucelsome_reads(bam_file, is_paired_end, nuc_pos_gr, tss_gr)
  
  # initiate counts
  nuc_pos_gr$nuc_read_count = 0
  
  # Optional step: update probability density with the iterations of EM algorithm
  if (is_EM) {
    
    if (!is_queitly)
      message("\n[2] Updating probability density function.")
    
    # record changes
    iter_records = data.frame(position = (-1000) : 1000,
                              probability = prob_den_fun((-1000) : 1000) / sum(prob_den_fun((-1000) : 1000)),
                              iteration = 0)
    for (i in seq_len(iter_num)) {
      gene_ids = sample(tss_gr$gene_id, 2000) # use 2000 genes for PDF update
      prob_den_fun = iterate_pdf(gene_ids, tss_gr, nuc_pos_gr, srg, is_paired_end, nuc_num)
      iter_records = rbind(iter_records,
                           data.frame(position = (-1000) : 1000,
                                      probability = prob_den_fun((-1000) : 1000) / sum(prob_den_fun((-1000) : 1000)),
                                      iteration = rep(i, 2001)))
      if (!is_queitly)
        message(paste("\nIteration loop:", i))
    }
    iter_records$iteration = factor(as.character(iter_records$iteration))
    
    ggsave(filename = paste0("Nucount_output/PDF_iteration_", out_name, ".png"),
           plot = ggplot(data = iter_records, aes(x = position, y = probability, color = iteration)) +
             geom_line(lwd = 1) +
             theme_minimal() +
             xlab("Position [bp]") + ylab("Probability") + ggtitle(out_name),
           device = "png", 
           width = 5, height = 4)
    
    
  } else {
    if (!is_queitly)
      message("\n[2] Default probability density function is applied.")
  }
  
  # loop through targets
  nuc_pos_gr$nuc_read_count <- foreach(gene_id = unique(nuc_pos_gr$gene_id),
                                       .combine = "c") %dopar% 
    {
      n = which(tss_gr$gene_id == gene_id)
      tss_tmp = tss_gr[n]
      .strand = as.character(strand(tss_tmp))
      .start = start(tss_tmp)
      
      # center to TSS interval start
      tmp_nuc = nuc_pos_gr[nuc_pos_gr$gene_id == tss_tmp$gene_id]
      
      
      srg_tmp = srg[[n]]
      if (length(srg_tmp[[1]]) == 0) return(rep(0, length(tmp_nuc)))
      
      
      # estimate reads to nucleosome, combine unique and multiple reads
      tmp_attribution = attribute_read_to_nucleosome(tss_tmp, srg_tmp, tmp_nuc, is_paired_end, prob_den_fun)
      
      tmp_attr_counts = `names<-`(rep(0, length(tmp_nuc)), seq_along(tmp_nuc))
      tmp_read_tbl = table(tmp_attribution)
      tmp_attr_counts[names(tmp_read_tbl)] = tmp_read_tbl
      
      # append reads count
      if (.strand == "-") tmp_attr_counts = rev(tmp_attr_counts)
      
      tmp_attr_counts
    } # the end of "foreach" 
  
  if (!is_queitly)
    message("\n[3] Reads allocation done.")
  
  nuc_pos_gr[as.numeric(gsub("-", "", nuc_pos_gr$nuc_pos)) %in% seq_len(nuc_num)]
}

# -------------------------------------------------- run ------------------------------------------------------- #
# load pre-processed nucleosome and TSS GRanges 

nuc_fltd_gr <- readRDS(nuc_file)
nuc_fltd_gr <- GenomicRanges::resize(nuc_fltd_gr, width = 1, fix = "center")
tss_gr <- readRDS(tss_file)

count_gr <- get_nucelsome_count(bam_file      = bam_file,
                                is_paired_end = is_paired_end,
                                nuc_pos_gr    = nuc_fltd_gr,
                                nuc_num       = nuc_num,
                                tss_gr        = tss_gr,
                                is_EM         = T,
                                iter_num      = iter_num)


# ------------------------------------------------- output ----------------------------------------------------- #
# Nucount results
if (out_type == "table") {
  write.table(mcols(count_gr), file = paste0("Nucount_output/", out_name, ".txt"),
              quote = F, sep = "\t", row.names = F, col.names = T)
} else if (out_type == "bw"){
  rtracklayer::export.bw(IRanges::coverage(count_gr + 25, weight = as.numeric(count_gr$nuc_read_count)),
                         con = paste0("Nucount_output/", out_name, ".bw"))
}


