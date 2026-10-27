# Project-GumDetective: inspect downloaded GEO files

project_dir <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(project_dir, "data", "raw")

files <- list.files(raw_dir, recursive = TRUE, full.names = TRUE)

if (length(files) == 0) {
  stop("No files found in data/raw. Download the GEO supplementary data there first.")
}

file_info <- file.info(files)

print(data.frame(
  file = files,
  size_MB = round(file_info$size / 1024^2, 2)
))