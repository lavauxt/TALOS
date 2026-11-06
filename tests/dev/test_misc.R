e <- new.env(); for (f in list.files("R", full.names = TRUE)) sys.source(f, e)
row <- data.frame(Length = 0L, GenomicPosition = 100L, HGVS_cDNA = "c.1_2dup", Gene = "FLT3", AlleleFrequency = 0.1,
  DepthAtBreakpoint = 50L, SupportingReads = 5L, WildtypeReads = 45L, HGVS_Protein = NA, StrandBias = 0.5,
  Hotspot = TRUE, HotspotName = "TKD; JMD", Orientation = "tandem", MedianMicrohomology = 3)
cat(e$.format_info(row), "\n")
df <- data.frame(GenomicPosition = c(1000L, 1002L, 5100L, 3000L), Length = c(4101L, 0L, 0L, 0L),
                 CandidateTypes = c("ptd_graph", "ptd_clip", "ptd_clip", "ptd_clip"))
print(e$.fold_clip_calls_into_graph(df, 10L, message))
