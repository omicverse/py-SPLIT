args <- commandArgs(trailingOnly=TRUE)
upstream <- normalizePath(args[[1]], winslash='/')
out <- args[[2]]
suppressPackageStartupMessages(library(Matrix))
source(file.path(upstream, 'R/method_agnostic_purification.R'))
source(file.path(upstream, 'R/residual_counts_reassignment.R'))
read_matrix <- function(name) as.matrix(read.csv(file.path(out,name), row.names=1, check.names=FALSE))
x <- read_matrix('split_counts.csv')
w <- read_matrix('split_weights.csv')
ref <- read_matrix('split_reference.csv')
obs <- read.csv(file.path(out, 'split_obs.csv'), row.names=1, check.names=FALSE)
primary <- setNames(obs$first_type, rownames(obs))
purified <- rctd_free_purify(x, w, ref, primary, DO_run_in_chunks=TRUE)$purified_counts
write.csv(as.matrix(t(purified)),file.path(out,'split_purified_r.csv'),row.names=FALSE)
# Shared balance decision isolates the redistribution kernel from scoring.
balanced <- x
selected <- obs$split_balance_status == 'purified'
balanced[,selected] <- as.matrix(purified[,selected])
nn <- as.matrix(read.csv(file.path(out,'split_neighbors.csv'),header=FALSE))
network <- list(cell_id=colnames(x), nn_idx=nn)
network$second_type_neighbors_no_reject <- lapply(seq_len(nrow(nn)), function(i) {
  which(obs$first_type[nn[i,]] == obs$second_type[i])
})
status <- setNames(obs$split_balance_status, rownames(obs))
reassigned <- reassign_residual_counts(x, balanced, network, status, mode='count_proportional')
write.csv(as.matrix(t(reassigned)),file.path(out,'split_reassigned_r.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R-sessionInfo.txt'))
