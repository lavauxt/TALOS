.mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
`%||%` <- function(x, y) if (is.null(x)) y else x
source("R/utils.R")
e <- new.env(); sys.source("R/utils.R", e)
.exons_df <- e$.exons_df
source("R/engine_bam_extractor.R")
source("R/engine_graph_ptd.R")

mk <- function(qname, flag, pos, cigar, mpos = NA, sa = NA, rname = "chr11", mrnm = "=", mapq = 60, qwidth = 150)
  data.frame(qname=qname, flag=flag, rname=rname, pos=pos, qwidth=qwidth, mapq=mapq, cigar=cigar,
             mrnm=mrnm, mpos=mpos, isize=0L, sa=sa, stringsAsFactors=FALSE)
S <- 2600L; E <- 6700L
dt <- rbind(
  # 3 split reads, both primary (right-clipped at E) and supplementary (left-clipped at S)
  do.call(rbind, lapply(1:3, function(i) rbind(
    mk(paste0("sr",i), 1+2+64, E-79, "80M70S", sa="chr11,2600,+,80S70M,60,0;"),
    mk(paste0("sr",i), 1+2+64+2048, S, "80H70M", sa=paste0("chr11,", E-79, ",+,80M70S,60,0;"))))),
  # 2 everted discordant pairs (left reverse at 2650, right forward at 6500)
  mk("dp1", 1+16+64, 2650, "150M", mpos=6500), mk("dp1", 1+32+128, 6500, "150M", mpos=2650),
  mk("dp2", 1+16+64, 2700, "150M", mpos=6450), mk("dp2", 1+32+128, 6450, "150M", mpos=2700),
  # normal pair (FR, 300bp)  -> must be ignored
  mk("np", 1+2+32+64, 4000, "150M", mpos=4200), mk("np", 1+2+16+128, 4200, "150M", mpos=4000),
  # deletion split read: right-clipped piece upstream, left-clipped downstream -> ignored
  mk("del", 65, 3000, "80M70S", sa="chr11,5000,+,80S70M,60,0;"),
  # translocation SA on another chromosome -> ignored
  mk("tra", 65, 3500, "80M70S", sa="chr9,1000,+,80S70M,60,0;"),
  # inversion (opposite strand) -> ignored
  mk("inv", 65, 6000, "80M70S", sa="chr11,3000,-,80S70M,60,0;")
)
sr <- .parse_sa_tags(dt); dp <- .extract_discordant_and_shadow(dt)
print(sr); print(dp)
stopifnot(nrow(sr) == 3, all(sr$start == S), all(sr$end == E), nrow(dp) == 2)

ex_plus <- data.frame(start = 1000*(1:10), end = 1000*(1:10)+100, exon_idx = 1:10)
cyc <- .graph_cycles_from_evidence(rbind(sr, dp), ex_plus)
str(cyc[[1]][c("ptd_signature","dup_start_exon","dup_end_exon","support_reads","dup_start","dup_end","dup_length","evidence_score","n_split","n_pairs")])
stopifnot(length(cyc) == 1, cyc[[1]]$dup_start_exon == "Exon_3", cyc[[1]]$dup_end_exon == "Exon_6",
          cyc[[1]]$dup_start == S, cyc[[1]]$dup_end == E, identical(cyc[[1]]$internal_path, paste0("Exon_", 3:6)))

# minus strand: transcript exon 1 is the rightmost
ex_minus <- data.frame(start = 1000*(1:10), end = 1000*(1:10)+100, exon_idx = 10:1)
cm <- .graph_cycles_from_evidence(rbind(sr, dp), ex_minus)
cat("minus:", cm[[1]]$ptd_signature, paste(cm[[1]]$internal_path, collapse=">"), "\n")
stopifnot(cm[[1]]$dup_start_exon == "Exon_5", cm[[1]]$dup_end_exon == "Exon_8")

# single exon duplication (exon 4 only)
ev1 <- data.frame(qname=c("a","b"), start=3900L, end=4300L, prob=0.95, evidence_type="split_read", precise=TRUE)
c1 <- .graph_cycles_from_evidence(ev1, ex_plus)
cat("single exon:", c1[[1]]$ptd_signature, "\n"); stopifnot(c1[[1]]$dup_start_exon == "Exon_4")
# single read is gated out
stopifnot(length(.graph_cycles_from_evidence(ev1[1,], ex_plus)) == 0)
cat("GRAPH TESTS OK\n")
