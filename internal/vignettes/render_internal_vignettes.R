args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (!length(file_arg)) {
  stop("Could not determine the script path from commandArgs().")
}

script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
docs_dir <- dirname(script_path)
output_dir <- file.path(docs_dir, "rendered")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- sort(list.files(
  docs_dir,
  pattern = "^[0-9]{2}-.*\\.Rmd$",
  full.names = TRUE
))

if (!length(inputs)) {
  stop("No vignette sources were found.")
}

for (input in inputs) {
  message("Rendering ", basename(input), " ...")
  rmarkdown::render(
    input = input,
    output_dir = output_dir,
    clean = TRUE,
    envir = new.env(parent = globalenv())
  )
}
