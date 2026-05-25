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
