#' Read an ANNOVAR frequency column
#'
#' Reads a single frequency column from a tab-delimited ANNOVAR output and
#' converts missing-value markers to numeric missing values.
#'
#' @param annovar_file Path to an ANNOVAR tab-delimited output file.
#' @param column Name of the frequency column to read.
#'
#' @return A numeric vector preserving the original ANNOVAR row order.
#'
#' @noRd
#'
#' @keywords internal
#' @importFrom data.table fread
.read_annovar_frequency_column <- function(annovar_file, column) {
  valid_file <- is.character(annovar_file) &&
    length(annovar_file) == 1L &&
    !is.na(annovar_file) &&
    nzchar(annovar_file)

  if (!valid_file) {
    stop(
      "`annovar_file` must be one non-empty character value.",
      call. = FALSE
    )
  }

  valid_column <- is.character(column) &&
    length(column) == 1L &&
    !is.na(column) &&
    nzchar(column)

  if (!valid_column) {
    stop(
      "`column` must be one non-empty character value.",
      call. = FALSE
    )
  }

  if (!file.exists(annovar_file)) {
    stop(
      sprintf("ANNOVAR file does not exist: %s", annovar_file),
      call. = FALSE
    )
  }

  header <- fread(
    file = annovar_file,
    sep = "\t",
    nrows = 0L,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE,
    nThread = 1L
  )

  column_matches <- which(names(header) == column)

  if (length(column_matches) == 0L) {
    stop(
      sprintf(
        "Column `%s` was not found in the ANNOVAR file.",
        column
      ),
      call. = FALSE
    )
  }

  if (length(column_matches) > 1L) {
    stop(
      sprintf(
        "Column `%s` occurs more than once in the ANNOVAR file.",
        column
      ),
      call. = FALSE
    )
  }

  frequency_table <- fread(
    file = annovar_file,
    sep = "\t",
    select = column,
    colClasses = list(character = column),
    na.strings = NULL,
    data.table = FALSE,
    check.names = FALSE,
    showProgress = FALSE,
    nThread = 1L
  )

  raw_values <- trimws(
    as.character(frequency_table[[column]])
  )

  missing <- is.na(raw_values) |
    !nzchar(raw_values) |
    raw_values == "." |
    toupper(raw_values) == "NA"

  frequency <- rep(NA_real_, length(raw_values))

  frequency[!missing] <- suppressWarnings(
    as.numeric(raw_values[!missing])
  )

  invalid <- !missing & (
    is.na(frequency) |
      !is.finite(frequency)
  )

  if (any(invalid)) {
    invalid_values <- unique(raw_values[invalid])
    preview <- invalid_values[
      seq_len(min(5L, length(invalid_values)))
    ]

    stop(
      sprintf(
        paste0(
          "Column `%s` contains non-numeric or non-finite ",
          "value(s): %s."
        ),
        column,
        paste(shQuote(preview), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  frequency
}

#' Add ABRaOM frequencies by ANNOVAR row order
#'
#' Reads an ABRaOM frequency column from an ANNOVAR output and assigns its
#' values directly to the allele records stored in an `ICAAlleles` object.
#'
#' @details
#' The ANNOVAR file and `x` must originate from the same VCF and retain the
#' same split-allele order. This function deliberately performs no coordinate
#' matching, normalization, validation, or join.
#'
#' The only correspondence check is that the ANNOVAR frequency column and
#' the allele-level `DataFrame` have the same number of rows. This prevents
#' unintended vector recycling but does not verify genomic alignment.
#'
#' Missing values represented by `"."`, an empty field, or `"NA"` are
#' converted to `NA_real_`.
#'
#' The ICA header and object metadata are preserved.
#'
#' @param x An `ICAAlleles` object returned by [split_ica_alleles()].
#' @param annovar_file Path to the corresponding tab-delimited ANNOVAR
#'   output.
#' @param column Name of the ABRaOM frequency column.
#'
#' @return The modified and validated `ICAAlleles` object with
#'   `af_abraom` populated.
#'
#' @examples
#' annovar_file <- tempfile(fileext = ".txt")
#'
#' writeLines(
#'     c(
#'         "Chr\tStart\tabraom_freq",
#'         "chr1\t100\t0.125",
#'         "chr2\t200\t."
#'     ),
#'     annovar_file
#' )
#'
#' allele_table <- S4Vectors::DataFrame(
#'     site_index = c(1L, 2L),
#'     allele_index = c(1L, 1L),
#'     chromosome = c("chr1", "chr2"),
#'     position = c(100L, 200L),
#'     refAllele = c("A", "G"),
#'     altAllele = c("C", "T"),
#'     variant_id = c(
#'         "chr1-100-A-C",
#'         "chr2-200-G-T"
#'     ),
#'     af_gnomad_genome = rep(NA_real_, 2L),
#'     af_gnomad_exome = rep(NA_real_, 2L),
#'     af_gnomad = rep(NA_real_, 2L),
#'     af_1000g = rep(NA_real_, 2L),
#'     af_abraom = rep(NA_real_, 2L),
#'     dbsnp_id = rep(NA_character_, 2L)
#' )
#'
#' x <- ICAAlleles(alleles = allele_table)
#'
#' x <- add_abraom_frequency(
#'     x,
#'     annovar_file
#' )
#'
#' allele_data(x)$af_abraom
#' unlink(annovar_file)
#'
#' @export
add_abraom_frequency <- function(
    x,
    annovar_file,
    column = "abraom_freq") {
  if (!methods::is(x, "ICAAlleles")) {
    stop(
      "`x` must be an ICAAlleles object.",
      call. = FALSE
    )
  }

  methods::validObject(x)

  alleles <- allele_data(x)

  abraom_frequency <- .read_annovar_frequency_column(
    annovar_file = annovar_file,
    column = column
  )

  if (length(abraom_frequency) != nrow(alleles)) {
    stop(
      sprintf(
        paste0(
          "Positional assignment requires identical row counts: ",
          "`x` contains %d allele records, but column `%s` in ",
          "the ANNOVAR file contains %d rows."
        ),
        nrow(alleles),
        column,
        length(abraom_frequency)
      ),
      call. = FALSE
    )
  }

  alleles$af_abraom <- abraom_frequency
  allele_data(x) <- alleles

  x
}
