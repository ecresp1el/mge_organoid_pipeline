# Safely inspect a serialized RDS object without triggering computation

## ---- user input -------------------------------------------------------------
input_rds_path <- "/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds"

## ---- setup output -----------------------------------------------------------
input_dir <- dirname(input_rds_path)
output_dir <- file.path(input_dir, "r_preliminary_output_folder")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

text_output_path <- file.path(output_dir, "object_structure_overview.txt")
provenance_output_path <- file.path(output_dir, "object_structure_provenance.rds")

## ---- helpers ----------------------------------------------------------------
get_dimensions <- function(x) {
  dims <- NULL
  dims <- tryCatch(dim(x), error = function(e) NULL)
  if (is.null(dims) && is.data.frame(x)) {
    dims <- c(nrow(x), ncol(x))
  }
  dims
}

get_empty_status <- function(x) {
  if (is.null(x)) {
    return("NULL")
  }
  len <- tryCatch(length(x), error = function(e) NA_integer_)
  if (is.na(len)) {
    return("length_unknown")
  }
  if (len == 0L) {
    return("empty")
  }
  "non_empty"
}

capture_structure <- function(x) {
  capture.output(
    utils::str(
      x,
      max.level = 2,
      vec.len = 5,
      list.len = 5,
      give.attr = FALSE
    )
  )
}

## ---- load object ------------------------------------------------------------
read_error <- NULL
serialized_object <- tryCatch(
  readRDS(input_rds_path),
  error = function(e) {
    read_error <<- conditionMessage(e)
    NULL
  }
)

if (!is.null(read_error)) {
  writeLines(
    c(
      paste("Input path:", input_rds_path),
      "Status: FAILED to read RDS file.",
      paste("Error:", read_error),
      "",
      "Next steps:",
      "- Ensure any required packages for the serialized object are installed.",
      "- For Seurat objects, install the SeuratObject package in your R session."
    ),
    con = text_output_path
  )

  provenance_payload <- list(
    input_path = input_rds_path,
    class = NA_character_,
    size_bytes = NA_real_,
    read_error = read_error,
    session_info = utils::sessionInfo(),
    timestamp = Sys.time()
  )

  saveRDS(provenance_payload, provenance_output_path)

  stop("Failed to read RDS file. See output files for details.")
}

## ---- collect metadata -------------------------------------------------------
object_class <- class(serialized_object)
object_type <- typeof(serialized_object)
object_size <- utils::object.size(serialized_object)
object_dimensions <- get_dimensions(serialized_object)
object_empty_status <- get_empty_status(serialized_object)
structure_preview <- capture_structure(serialized_object)

## ---- write outputs ----------------------------------------------------------
writeLines(
  c(
    paste("Input path:", input_rds_path),
    paste("Class:", paste(object_class, collapse = ", ")),
    paste("Type:", object_type),
    paste("Size (bytes):", as.numeric(object_size)),
    paste("Size (pretty):", format(object_size, units = "auto")),
    paste("Empty status:", object_empty_status),
    if (!is.null(object_dimensions)) {
      paste("Dimensions:", paste(object_dimensions, collapse = " x "))
    } else {
      "Dimensions: NA"
    },
    "",
    "Structure preview (depth-limited):",
    structure_preview
  ),
  con = text_output_path
)

provenance_payload <- list(
  input_path = input_rds_path,
  class = object_class,
  size_bytes = as.numeric(object_size),
  session_info = utils::sessionInfo(),
  timestamp = Sys.time()
)

saveRDS(provenance_payload, provenance_output_path)

message("Inspection complete. Outputs written to: ", output_dir)
