count_file <- file.path(
  "data", "raw", "GSE171213_AllSample.counts.tsv.gz"
)

con <- gzfile(count_file, open = "rt")
preview <- readLines(con, n = 3)
close(con)

for (i in seq_along(preview)) {
  fields <- strsplit(preview[i], "\t", fixed = TRUE)[[1]]

  cat("\nLine", i, "has", length(fields), "fields\n")
  cat("First 12 fields:\n")
  print(fields[1:min(12, length(fields))])

  cat("Last 8 fields:\n")
  print(tail(fields, 8))
}