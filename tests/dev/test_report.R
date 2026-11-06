`%||%` <- function(x, y) if (is.null(x)) y else x
for (f in c("utils.R","report_assets.R","report.R")) sys.source(file.path("R", f), globalenv())

set.seed(3)
ex <- data.frame(start = 10000 + (0:9) * 3000, end = 10000 + (0:9) * 3000 + 150, exon_idx = 1:10)
cfg <- list(gene = "KMT2A", chrom = "chr11", strand = "+", transcript = "NM_001197104", build = "hg19",
            genomic_start = 9000, genomic_end = 40000, cds_offset = 100L,
            all_exons = ex, target_exons = ex[2:9, ], gene_settings = list(min_support = 5))
cov <- c(rpois(15899, 120), rpois(4101, 200), rpois(23000, 120))
mkrow <- function(pos, len, vaf, hgvs, types, rle = NA) data.frame(
  Sample = "S1", Gene = "KMT2A", Genome = "hg19", GenomicPosition = pos, Length = len,
  DupExons = if (len > 1000) "3-4" else NA, AlleleFrequency = vaf, SupportingReads = 25L, WildtypeReads = 200,
  DepthAtBreakpoint = 225L, HGVS_cDNA = hgvs, HGVS_Protein = NA, Region = "intronic", ExonNumber = NA,
  Orientation = "tandem", Hotspot = FALSE, HotspotName = NA, ITD_Sequence = paste(rep("ACGT", 20), collapse=""),
  ITDCoverageRLE = rle, MergedFromN = 1L, MergeFragmentWarning = FALSE, CandidateTypes = types,
  LeftSoftclipCount = 10L, RightSoftclipCount = 12L, stringsAsFactors = FALSE)
df <- rbind(mkrow(15900, 4101, 0.12, "c.4000-50_4500+20dup", "ptd_clip+ptd_graph"),
            mkrow(25010, 60, 0.31, "c.6000_6059dup", "bwd_jump", "3:20,5:40"),
            mkrow(30500, 0, 0.05, NA, "ptd_clip"))
out <- talos_html_report(df, gene_configs = list(KMT2A = cfg), output_file = "tests/dev/out/test_report.html",
                         title = "TALOS Report - KMT2A S1", show_config = TRUE,
                         coverages = list(S1_KMT2A = cov), plotly_js = "auto")
cat("size", file.size(out), "\n")
# empty result + ALU mode
talos_html_report(data.frame(), output_file = "tests/dev/out/empty.html")
alu <- data.frame(Sample="S1", Gene="KMT2A", Genome="hg19", InsertionSite=20000L, ALUSubtype="AluY",
  ALUOrientation="sense", TSD_Length=12L, TSD_Sequence="ACGTACGTACGT", PolyALength=20L, ALU5pTruncation=0L,
  ALU_Sequence="GGCCGGGCGCGG", ALU_RawScore=0.9, SupportingReads=8L, WildtypeReads=50L, DepthAtBreakpoint=58L,
  AlleleFrequency=0.14, MeanSupportMAPQ=58, StrandBias=0.5, Hotspot=FALSE, HotspotName=NA, Region="intronic", ExonNumber=NA)
talos_html_report(alu, gene_configs = list(KMT2A = cfg), output_file = "tests/dev/out/alu.html", mode = "alu")
