# R script for nucleosomal reads assignment
#
# Source code: https://github.com/shaorray/Nucount
#
# MIT License
# Copyright (c) 2022 Rui Shao <shaorray@hotmail.com>

# --------------------------------------------- initialization ------------------------------------------------- #
args <- commandArgs(TRUE)

pkgs <- c("optparse", "dplyr", "foreach", "doParallel", "ggplot2", "BiocManager")
bioc_pkgs <- c("S4Vectors", "Rsamtools", "rtracklayer", "GenomeInfoDb", "IRanges",
               "GenomicRanges", "BiocGenerics")

suppressMessages(for (pkg in pkgs) {
  if (!require(pkg, character.only = T)) install.packages(pkg)
  require(pkg, quietly = T, character.only = T)
})

suppressMessages(for (pkg in bioc_pkgs) {
  if (!require(pkg, character.only = T)) BiocManager::install(pkg)
  require(pkg, quietly = T, character.only = T)
})

options("scipen" = 999, "digits" = 4)
options(warn = -1)

# --------------------------------------------- specifications ------------------------------------------------- #

option_list = list(
  # bam input
  make_option(c("-b", "--bam"), action = "store", default = NULL, type = 'character',
              help = "Bam file, single or paired end."),
  make_option(c("-p", "--paired"), action = "store", default = TRUE, type = 'logical',
              help = "Paired-end reads [default %default]"),
  
  # annotations
  make_option(c("-n", "--nucleosome"), action = "store", default = NULL, type = 'character',
              help = "Filtered TSS flanking nucleosome GRanges in RData format."),
  make_option(c("-g", "--gene_target"), action = "store", default = NULL, type = 'character',
              help = "Filtered gene TSS GRanges in RData format."),
  
  # tuning algorithm
  make_option(c("-e", "--is_EM"), action = "store", default = FALSE, type = 'logical',
              help = "Perform EM optimization of density function [default %default]"),
  make_option("--iter_num", default = 2,
              help = "Number of EM iteration [default %default]"),
  
  make_option(c("-t", "--thread"), action = "store", default = 1, type = 'numeric',
              help = "Number of CPU threads [default %default]"),
  make_option(c("-s", "--scaling"), action = "store", default = 1, type = 'numeric',
              help = "Sample size scaling number [default %default]"),
  
  
  # output
  make_option("--nuc_num", action = "store", default = 10, type = 'integer',
              help = "Number of nucleosomes away from TSS to keep [default %default]"),
  make_option("--type", action = "store", default = "table", type = 'character',
              help = "Output file format (BW, table, or matrix) [default %default]"),
  make_option(c("-o", "--out"), action = "store", default = "est_Nucount", type = 'character',
              help = "Output file name"),
  
  make_option(c("-q", "--quietly"), action = "store", default = FALSE, type = 'logical',
              help = "Mute verbos and messages [default %default]")
)
arguments <- parse_args(OptionParser(option_list=option_list), positional_arguments=TRUE)
opt <- arguments$options


# parse arguments
bam_file <- opt$bam
is_paired_end <- opt$paired

nuc_file <- opt$nucleosome
gene_target_file <- opt$gene_target

is_EM <- opt$is_EM
iter_num <- opt$iter_num

sample_scale <- opt$scaling

nuc_num <- opt$nuc_num
out_name <- opt$out
out_type <- opt$type

is_queitly <- opt$quietly

if (!is_queitly) {
  message("\n----------------------------------------------\n")
  message(paste("Nucount v1.0: ", date(), "\n"))
  message("----------------------------------------------\n")
  message(paste("BAM:\n\n", bam_file, "\n"))
  message("----------------------------------------------\n")
  message(paste("Is read paired:\n\n", is_paired_end, "\n"))
  message("----------------------------------------------\n")
  message(paste("Sample size factor:\n\n", sample_scale, "\n"))
  message("----------------------------------------------\n")
  message(paste("Nucleosome file:\n\n", nuc_file, "\n"))
  message("----------------------------------------------\n")
  message(paste("Gene target file:\n\n", gene_target_file, "\n"))
  message("----------------------------------------------\n")
  message(paste("Perform EM:\n\n", is_EM, "\n"))
  message("----------------------------------------------\n")
  message(paste("EM iteration time:\n\n", iter_num, "\n"))
  message("----------------------------------------------\n")
  message(paste("Output name:\n\n", out_name, "\n"))
  message("----------------------------------------------\n")
  message(paste("Output type:\n\n", out_type, "\n"))
  message("----------------------------------------------\n")
}

# create output folder if not exists
if (!file.exists("Nucount_output"))
  dir.create("Nucount_output")


# indicate the number of multi-thread cores
doParallel::registerDoParallel(cores = opt$thread)

# ----------------------------------------------- functions ---------------------------------------------------- #

# initialize with a normal distribution with mean = 0, sd = 20
prob_den_fun = function(x) exp(-(x / 20)^2 / 280) # function (0)

bam_nucelsome_reads = function(bam_file, is_paired_end = T, tss_gr) {
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
    srg <- lapply(srg, 
                  function(x) {
                    idx = with(x, !is.na(mpos) & rname == mrnm & abs(isize) > 10 & abs(isize) < 1000)
                    lapply(x, function(y) y[idx])})
  }
  
  if (!is_queitly)
    message(paste("Number of reads:", sum(unlist(lapply(srg, function(x) length(x[[1]])))), "\n"))
  
  srg # short read ranges
}


attribute_read_to_nucleosome = function(tss_tmp, srg_tmp, nuc_tmp, is_paired_end, prob_den_fun) {
  # function (2)
  
  .start = start(tss_tmp)
  
  nuc_tmp_ranges = sort(ranges(shift(nuc_tmp, shift = -.start)))
  
  if (is_paired_end) {
    read_ranges = IRanges(start = with(srg_tmp, ifelse(strand == "+", pos, mpos)) - .start, 
                          width = abs(srg_tmp$isize))
  } else {
    srg_tmp_start = with(srg_tmp, ifelse(strand == "+", pos, pos + qwidth)) - .start
    read_ranges = IRanges(start = srg_tmp_start, width = 1)
  }
  
  out_attribution = rep(NA, length(read_ranges)) # initiate read-nucleosome attribution
  
  # Step1: allocate unique read-nucleosome overlaps
  rov_cnt = countQueryHits(findOverlaps(query = read_ranges, subject = nuc_tmp_ranges))
  uni_idx = rov_cnt == 1
  uni_mtch = findOverlaps(query = read_ranges[uni_idx], subject = nuc_tmp_ranges)
  uni_mtch_ls = split(queryHits(uni_mtch), subjectHits(uni_mtch))
  
  tmp_count = countSubjectHits(uni_mtch)
  
  out_attribution[uni_idx][queryHits(uni_mtch)] = subjectHits(uni_mtch)
  
  # Step2: EM on multi-overlapped nucleosomal reads
  if (sum(!uni_idx) > 0) {
    read_ranges_multi = read_ranges[!uni_idx]
    multi_mat = cbind(start = start(read_ranges_multi), end = end(read_ranges_multi))
    
    # keep reads within 1 kb
    keep_idx = apply(multi_mat, 1, function(x) 
      pmin(min(abs(x[1] - start(nuc_tmp_ranges))),
           min(abs(x[2] - start(nuc_tmp_ranges)))) < 1000)
    
    # discard reads without any overlap of nucleosome centers
    if (is_paired_end) keep_idx = keep_idx & (rov_cnt > 1)[!uni_idx]  
    if (sum(!keep_idx) == nrow(multi_mat)) {
      return(out_attribution)
    } else {
      multi_mat = matrix(multi_mat[keep_idx, ], ncol = 2)
    }
    
    # sum the estimated probability from all qualified reads 
    if (is_paired_end) {
      nuc_tmp_prob = 
        sapply(seq_along(nuc_tmp_ranges), 
               function(x) {
                 x_nuc_start = start(nuc_tmp_ranges[x])
                 x_nuc_prob = sum(pmax(prob_den_fun(x - multi_mat[, 1]), 
                                       prob_den_fun(x - multi_mat[, 2])))
                 # specify for unique overlapped reads
                 if (x %in% names(uni_mtch_ls)) {
                   tmp_read_ranges = read_ranges[uni_idx][uni_mtch_ls[[as.character(x)]]]
                   x_nuc_prob = x_nuc_prob +
                     sum(pmax(prob_den_fun(start(nuc_tmp_ranges)[x] - start(tmp_read_ranges)), 
                              prob_den_fun(start(nuc_tmp_ranges)[x] - end(tmp_read_ranges))))
                 }
                 x_nuc_prob
               }
        )
    } else {
      # single-end reads must towards their dedicated nucleosomes
      nuc_tmp_prob = 
        sapply(seq_along(nuc_tmp_ranges), 
               function(x) {
                 x_nuc_start = start(nuc_tmp_ranges[x])
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
                              which.max(nuc_tmp_prob * prob_den_fun(x - start(nuc_tmp_ranges))))
    
    out_attribution[!uni_idx][keep_idx] = tmp_attribution
  } # end of EM
  
  out_attribution
}


iterate_pdf = function(gene_ids, tss_gr, nuc_pos_gr, srg, is_paired_end, nuc_num) {
  # function (3)
  
  read_nucleosome_coverage = foreach(gene_id = gene_ids, .combine = "c") %dopar% {
    n = which(tss_gr$gene_id == gene_id)
    srg_tmp = srg[[n]]
    if (length(srg_tmp[[1]]) == 0) return(NULL)
    tss_tmp = tss_gr[n]
    
    nuc_tmp = nuc_pos_gr[nuc_pos_gr$gene_id == tss_tmp$gene_id]
    
    tmp_attribution = attribute_read_to_nucleosome(tss_tmp, 
                                                   srg_tmp,
                                                   nuc_tmp, 
                                                   is_paired_end,
                                                   prob_den_fun)
    
    if (sum(is.na(tmp_attribution)) == 0) return(NULL)
    
    # relative position between reads and nucleosome
    delta_pos = with(srg_tmp, pos + ifelse(strand == "+", 0, qwidth)) - start(sort(nuc_tmp))[tmp_attribution]
    
    if (is_paired_end) { # include mate read end
      delta_pos = delta_pos +
        with(srg_tmp, ifelse(strand == "+", mpos + qwidth, mpos)) - start(sort(nuc_tmp))[tmp_attribution]
    }
      
    if (as.character(strand(tss_tmp)) == "-") delta_pos = delta_pos * (-1)
    
    delta_pos[which(gsub("-", "", nuc_tmp$nuc_pos[tmp_attribution]) %in% seq_len(nuc_num))]
  } # end of foreach
  
  approxfun(density(read_nucleosome_coverage), yleft = 0, yright = 0)
}


# main function
get_nucelsome_count = function(bam_file,            # bam file input
                               is_paired_end = T,   # if reads paired end
                               nuc_pos_gr,          # saved nucleosome GRanges
                               nuc_num,             # the number of nucleosomes to use, so symmetric outputs = nuc_num * 2  
                               tss_gr,              # saved TSS GRanges
                               is_EM = T,           # if update probability density function
                               iter_num = 2,        # number of EM iterations, pdf could fluctuate depending on ChIP type, manually pick a small number
                               out_name,            # user defined sample name
                               is_queitly = FALSE   # if show message
) {

  tss_gr = tss_gr[tss_gr$gene_id %in% nuc_pos_gr$gene_id]
  nuc_pos_gr = nuc_pos_gr[nuc_pos_gr$gene_id %in% tss_gr$gene_id]
  gene_ids = tss_gr$gene_id
  
  # process bam file
  srg = bam_nucelsome_reads(bam_file, is_paired_end, tss_gr)
  
  # initiate counts
  nuc_pos_gr$nuc_read_count = 0
  
  # Optional step: update probability density with the iterations of EM algorithm
  if (is_EM) {
    
    if (!is_queitly)
      message("\n[2] Updating probability density function.")
    
    # record changes
    positions = (-1000) : 1000
    iter_records = data.frame(position = positions,
                              probability = prob_den_fun(positions) / sum(prob_den_fun(positions)),
                              iteration = 0)
    for (i in seq_len(iter_num)) {
      
      if (length(gene_ids) > 4e3) { # use only 2000 genes for PDF update
        use_gene_ids = sample(gene_ids, 2000, replace = T) 
      } else {
        use_gene_ids = gene_ids
      }
      
      prob_den_fun <<- iterate_pdf(gene_ids      = use_gene_ids,
                                   tss_gr        = tss_gr,
                                   nuc_pos_gr    = nuc_pos_gr,
                                   srg           = srg, 
                                   is_paired_end = is_paired_end, 
                                   nuc_num       = nuc_num)
      
      iter_records = rbind(iter_records,
                           data.frame(position = positions,
                                      probability = prob_den_fun(positions) / sum(prob_den_fun(positions)),
                                      iteration = rep(i, 2001))
                           )
      if (!is_queitly)
        message(paste("\nIteration loop:", i))
    }
    iter_records$iteration = factor(as.character(iter_records$iteration))
    
    # print PDF iteration records 
    if (!file.exists("Nucount_output")) dir.create("Nucount_output")
    
    ggsave(filename = paste0("Nucount_output/PDF_iteration_", out_name, ".png"),
           plot = ggplot(data = iter_records,
                         aes(x = position, y = probability, color = iteration)) +
             geom_line(lwd = 1) +
             ggplot2::theme_minimal() +
             xlab("Position [bp]") + ylab("Probability") + ggtitle(out_name),
           device = "png",
           width = 5, height = 4)
    
  } else {
    
    if (!is_queitly)
      message("\n[2] Default probability density function is applied.")
    
  }
  
  # loop through targets
  nuc_pos_gr$nuc_read_count = foreach(gene_id = gene_ids, .combine = "c") %dopar% 
    {
      n = which(tss_gr$gene_id == gene_id)
      tss_tmp = tss_gr[n]
      
      # center to TSS interval start
      nuc_tmp = nuc_pos_gr[nuc_pos_gr$gene_id == tss_tmp$gene_id]
      
      srg_tmp = srg[[n]]
      if (length(srg_tmp[[1]]) == 0) return(rep(0, length(nuc_tmp)))
      
      
      # estimate reads to nucleosome, combine unique and multiple reads
      tmp_attribution = attribute_read_to_nucleosome(tss_tmp,
                                                     srg_tmp, 
                                                     nuc_tmp, 
                                                     is_paired_end,
                                                     prob_den_fun)
      
      tmp_attr_counts = `names<-`(rep(0, length(nuc_tmp)), seq_along(nuc_tmp))
      tmp_read_tbl = table(tmp_attribution)
      tmp_attr_counts[names(tmp_read_tbl)] = tmp_read_tbl
      
      # append reads count
      if (as.character(strand(tss_tmp)) == "-") tmp_attr_counts = rev(tmp_attr_counts)
      
      tmp_attr_counts
    } # the end of "foreach" 
  
  if (!is_queitly)
    message("\n[3] Reads allocation done.")
  
  nuc_pos_gr[as.numeric(gsub("-", "", nuc_pos_gr$nuc_pos)) %in% seq_len(nuc_num)]
}

assign_tss_nucleosome = function(nucleosome_gr, tss_gr, n_nucleosomes, is_queitly) {
  # this function fetches the indicated number of nucleosome around TSSs
  
  if (!is_queitly)
    message("\n[-] Nucleosome-TSS assigment begins.")
  
  nucleosome_gr = sort(nucleosome_gr)
  
  if (is.null(tss_gr$gene_id)) {
    gene_ids = seq_along(tss_gr)
  } else {
    gene_ids = tss_gr$gene_id
  } 
  
  mtch = findOverlaps(nucleosome_gr, tss_gr)
  mtch = split(queryHits(mtch), subjectHits(mtch))
  
  nuc_pos_gr = foreach(tss = names(mtch), .combine = "c") %dopar% {
    
    tss_i = as.numeric(tss)
    nuc_tmp = nucleosome_gr[mtch[[tss]]]
    tmp_tss = tss_gr[tss_i]
    
    dist_tss = start(nuc_tmp) - start(resize(tmp_tss, 0, "center"))
    
    if (as.character(strand(tmp_tss)) == "-") {
      dist_tss = rev(dist_tss * -1) 
      nuc_tmp = rev(nuc_tmp)
    } 
    
    if (sum(dist_tss >= 0) < n_nucleosomes | sum(dist_tss < 0) < n_nucleosomes) {
      nuc_tmp = GRanges()
    } else {
      nuc_tmp = c(tail(nuc_tmp[dist_tss < 0], n_nucleosomes), 
                  head(nuc_tmp[dist_tss >= 0], n_nucleosomes))
      nuc_tmp$gene_id = gene_ids[tss_i]
      nuc_tmp$nuc_pos = c(rev(seq_len(n_nucleosomes) * -1), seq_len(n_nucleosomes))
    }
    nuc_tmp
  }
  
  if (!is_queitly)
    message("\n[0] Nucleosome-TSS assigment done.")
  
  return(nuc_pos_gr)
}

column_to_matrix = function(count_gr, tss_gr, nuc_num) {
  # this function converts count_gr "nuc_read_count" to a nucleosome position matrix by TSS gene_ids
  
  mat_count = t(matrix(count_gr$nuc_read_count, nrow = nuc_num *2))
  colnames(mat_count) = unique(count_gr$nuc_pos)
  rownames(mat_count) = unique(count_gr$gene_id)
  
  mat_out = matrix(NA, ncol = nuc_num *2, nrow = length(tss_gr))
  rownames(mat_out) = tss_gr$gene_id
  colnames(mat_out) = colnames(mat_count)
  
  mat_out[rownames(mat_count), ] = mat_count
  
  mat_out
}


# ------------------------------------------------- input ------------------------------------------------------ #
# load TSS ranges, from formats {rds, BED, GTF, etc.} 

gene_target_file_format <- toupper(sapply(strsplit(gene_target_file, '\\.'), tail, 1))
if (gene_target_file_format %in% c("RDS", "RDATA")) {
  tss_gr <- readRDS(gene_target_file)
} else {
  gene_gr <- rtracklayer::import(gene_target_file, format = gene_target_file_format)
  seqlevelsStyle(gene_gr) <- "UCSC"
  gene_gr <- gene_gr[seqnames(gene_gr) %in% paste0("chr", c(1:100, "X", "Y"))]
  seqlevels(gene_gr) <- as.character(unique(seqnames(gene_gr)))
  tss_gr <- GenomicFeatures::promoters(gene_gr, upstream = 5e3, downstream = 5e3)
  if (is.null(tss_gr$gene_id)) tss_gr$gene_id = paste0("tss_", seq_along(tss_gr))
}

# load pre-processed nucleosomes ranges
nuc_gr <- readRDS(nuc_file)

if (is.null(nuc_gr$nuc_pos)) {
  # add 5 flanking nucleosomes to buffer boundary bias
  with_progress(
    nuc_gr <- assign_tss_nucleosome(nucleosome_gr = nuc_gr, 
                                    tss_gr        = tss_gr, 
                                    n_nucleosomes = nuc_num + 5, 
                                    is_queitly    = is_queitly)
    )
}

nuc_fltd_gr <- GenomicRanges::resize(nuc_gr, width = 1, fix = "center")

# --------------------------------------------------- run ------------------------------------------------------ #

count_gr <- get_nucelsome_count(bam_file      = bam_file,
                                is_paired_end = is_paired_end,
                                nuc_pos_gr    = nuc_fltd_gr,
                                nuc_num       = nuc_num,
                                tss_gr        = tss_gr,
                                is_EM         = is_EM,
                                iter_num      = iter_num,
                                out_name      = out_name,
                                is_queitly    = is_queitly)

count_gr$nuc_read_count <- count_gr$nuc_read_count / sample_scale

# ------------------------------------------------- output ----------------------------------------------------- #
# save Nucount results
dir.create(paste0("Nucount_output/", out_name))

if (out_type == "table") {
  write.table(mcols(count_gr),
              file = paste0("Nucount_output/", out_name, "/", out_name, "_tab.txt"),
              quote = F, sep = "\t",
              row.names = F,
              col.names = T)
} else if (out_type == "bw") {
  rtracklayer::export.bw(
    IRanges::coverage(count_gr + 25, # 51 bp in coverage width
                      weight = as.numeric(count_gr$nuc_read_count)),
    con = paste0("Nucount_output/", out_name, "/", out_name, ".bw")
    )
} else {
  write.table(column_to_matrix(count_gr, tss_gr, nuc_num),
              file = paste0("Nucount_output/", out_name, "/", out_name, "_mat.txt"),
              quote = F, sep = "\t",
              row.names = T,
              col.names = T)
}

