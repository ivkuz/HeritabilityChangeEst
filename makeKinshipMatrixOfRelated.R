library(data.table)

threshold_grm <- function(prefix_in, prefix_out, 
                          threshold = 0.05, 
                          chunk_size = 1e7) {
  
  # Read sample IDs
  id <- read.table(paste0(prefix_in, ".grm.id"))
  n  <- nrow(id)
  total_vals <- n * (n + 1) / 2
  
  # Open input and output binary connections
  con_in  <- file(paste0(prefix_in,  ".grm.bin"), "rb")
  con_out <- file(paste0(prefix_out, ".grm.bin"), "wb")
  
  vals_processed <- 0
  print("started")
  
  while (vals_processed < total_vals) {
    
    # Read next chunk
    to_read <- min(chunk_size, total_vals - vals_processed)
    chunk   <- readBin(con_in, what = "numeric", 
                       n = to_read, size = 4)
    
    # Threshold - keep only relatives
    chunk[chunk < threshold] <- 0
    
    # Write chunk
    writeBin(as.numeric(chunk), con_out, size = 4)
    
    vals_processed <- vals_processed + to_read
    
    print(paste0(vals_processed, " processed out of ", format(total_vals, scientific = T)))
    
  }
  
  close(con_in)
  close(con_out)
  
  # # Copy N file unchanged - no modification needed
  # file.copy(paste0(prefix_in,  ".grm.N.bin"),
  #           paste0(prefix_out, ".grm.N.bin"))
  
  # # Write ID file
  # write.table(id, paste0(prefix_out, ".grm.id"),
  #             quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Copy id file unchanged - no modification needed
  file.copy(paste0(prefix_in,  ".grm.id"),
            paste0(prefix_out, ".grm.id"))
  
  cat("Done. Processed", total_vals, "values in",
      ceiling(total_vals / chunk_size), "chunks.\n")
}

threshold_grm("~/EA_heritability/h2/grm/ldak/kins/kinships.all", "~/EA_heritability/h2/grm/ldak/kins/kinships.all.threshold_0.05", threshold = 0.05, chunk_size = 1e8)







gxe_grm <- function(prefix_in, prefix_out,
                          group1_ids, group2_ids,
                          chunk_size = 1e7) {
  
  # Read sample IDs
  id <- read.table(paste0(prefix_in, ".grm.id"))
  n  <- nrow(id)
  total_vals <- n * (n + 1) / 2
  
  # Build a group lookup vector (1 or 2) indexed by row position in id file
  sample_names <- id[[2]]
  group_vec <- integer(n)
  group_vec[sample_names %in% group1_ids] <- 1L
  group_vec[sample_names %in% group2_ids] <- 2L
  
  # Pre-compute the group for each sample index (1-based)
  # GRM lower triangle is stored row by row:
  # (1,1), (2,1),(2,2), (3,1),(3,2),(3,3), ...
  # For value at flat index k, we need the (row, col) pair.
  # We precompute row/col for all values upfront for efficiency.
  # Row i contains i values (cols 1..i), cumulative count = i*(i+1)/2
  # So for flat index k (1-based): row i = ceil((sqrt(8k+1)-1)/2)
  
  ### To avoid computing row/col on-the-fly per chunk, we'll track position
  ### with a running (row, col) cursor instead.

  # Open connections
  con_in  <- file(paste0(prefix_in,  ".grm.bin"), "rb")
  con_out <- file(paste0(prefix_out, ".grm.bin"), "wb")
  
  vals_processed <- 0L
  cur_row <- 1L  # current row in lower triangle (1-based)
  cur_col <- 1L  # current col in lower triangle (1-based, col <= row)
  
  message("Started")
  
  while (vals_processed < total_vals) {
    
    to_read <- min(chunk_size, total_vals - vals_processed)
    chunk   <- readBin(con_in, what = "numeric", n = to_read, size = 4)
    
    # Build cross-group mask for this chunk
    # Walk through each value's (row, col) position
    mask <- logical(to_read)
    # r <- cur_row
    # c <- cur_col
    # 
    # for (k in seq_len(to_read)) {
    #   mask[k] <- (group_vec[r] != 0L & group_vec[c] != 0L &   # both assigned
    #                 group_vec[r] != group_vec[c])                 # different groups
    #   # Advance cursor
    #   if (c == r) { r <- r + 1L; c <- 1L } else { c <- c + 1L }
    # }
    # 
    # # Update cursor for next chunk
    # cur_row <- r
    # cur_col <- c

    
    # Vectorised (row, col) computation for the chunk
    # flat indices (1-based) for this chunk
    
    # idx <- vals_processed + seq_len(to_read)          # global 1-based positions
    idx <- as.numeric(vals_processed) + seq_len(to_read)  # double arithmetic
    
    # row_idx <- as.integer(ceiling((sqrt(8 * idx + 1) - 1) / 2))
    # col_idx <- as.integer(idx - row_idx * (row_idx - 1L) / 2)
    row_idx <- ceiling((sqrt(8 * idx + 1) - 1) / 2)  # keep as double, drop as.integer()
    col_idx <- idx - row_idx * (row_idx - 1) / 2      # double arithmetic throughout
    
    # mask <- (group_vec[row_idx] != 0L &
    #            group_vec[col_idx] != 0L &
    #            group_vec[row_idx] != group_vec[col_idx])    
    mask <- (group_vec[as.integer(row_idx)] != 0L &
               group_vec[as.integer(col_idx)] != 0L &
               group_vec[as.integer(row_idx)] != group_vec[as.integer(col_idx)])
    
    # Zero out cross-group pairs
    chunk[mask] <- 0
    
    writeBin(as.numeric(chunk), con_out, size = 4)
    
    vals_processed <- vals_processed + to_read
    message(vals_processed, " processed out of ", 
            format(total_vals, scientific = TRUE))
  }
  
  close(con_in)
  close(con_out)
  
  file.copy(paste0(prefix_in,  ".grm.id"),
            paste0(prefix_out, ".grm.id"))
  
  cat("Done. Processed", total_vals, "values in",
      ceiling(total_vals / chunk_size), "chunks.\n")
}

gxe_grm(prefix_in = "~/EA_heritability/h2/grm/ldak/kins/kinships.all", 
        prefix_out = "~/EA_heritability/h2/grm/ldak/kins/kinships.all.byera", 
        group1_ids = fread("~/EA_heritability/gcta/data/s_ldak_list_all_15.tsv")$V2, 
        group2_ids = fread("~/EA_heritability/gcta/data/ps_ldak_list_all_15.tsv")$V2, 
        chunk_size = 1e7)

groups12 <- rbind(fread("~/EA_heritability/gcta/data/s_ldak_list_all_15.tsv"), 
                  fread("~/EA_heritability/gcta/data/ps_ldak_list_all_15.tsv"))
write.table(groups12, "~/EA_heritability/gcta/data/all_ldak_list_all.tsv", 
            row.names = F, col.names = F, quote = F, sep = "\t")
