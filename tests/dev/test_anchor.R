source("R/engine_candidates.R")
set.seed(1)
ref <- paste(sample(c("A","C","G","T"), 3000, TRUE), collapse="")
S <- 1001L; L <- 300L; E <- S + L - 1L
sub <- function(a,b) substr(ref,a,b)
for (c in c(25L, 35L, 60L)) {
  lead <- sub(E - c + 1L, E)                 # read aligned from S, leading clip = end of copy 1
  r1 <- .estimate_itd_from_anchors(lead, NA, ref, orig_bp_local = S, kmer_len = L)
  trail <- sub(S, S + c - 1L)                # read aligned up to E, trailing clip = start of copy 2
  r2 <- .estimate_itd_from_anchors(NA, trail, ref, orig_bp_local = S, kmer_len = L, trail_bp_local = E)
  cat(sprintf("clip=%d lead=%s trail=%s\n", c, r1, r2))
  stopifnot(r1 == L, r2 == L)
}
# wrong kmer_len -> banded pass empty -> unbanded pass recovers the true length
stopifnot(.estimate_itd_from_anchors(sub(E-29L,E), NA, ref, S, kmer_len = 150L) == L)
cat("ANCHOR TESTS OK\n")
