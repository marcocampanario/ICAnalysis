#' Validate an ICAPositions object
#'
#' @param x An object to validate.
#' @param argument Character scalar naming the argument being validated.
#'
#' @return `x`, invisibly.
#'
#' @noRd
#' @keywords internal
.assert_ica_positions <- function(x, argument = "x") {
  if (!methods::is(x, "ICAPositions")) {
    stop(
      sprintf(
        "`%s` must be an `ICAPositions` object.",
        argument
      ),
      call. = FALSE
    )
  }

  invisible(x)
}


#' Check whether a VCF INFO value is missing
#'
#' Considers `NULL`, zero-length values, missing values, empty strings, and
#' `"."` to represent missing VCF information.
#'
#' @param value A potentially nested VCF INFO value.
#'
#' @return A logical scalar.
#'
#' @noRd
#' @keywords internal
.is_missing_vcf_value <- function(value) {
  if (is.null(value) || length(value) == 0L) {
    return(TRUE)
  }

  values <- unlist(
    value,
    recursive = TRUE,
    use.names = FALSE
  )

  if (length(values) == 0L) {
    return(TRUE)
  }

  text <- trimws(as.character(values))

  all(
    is.na(values) |
      !nzchar(text) |
      text == "."
  )
}


#' Extract one field from a VCF INFO record
#'
#' @param vcf_info A named list or data-frame-like VCF INFO record.
#' @param field Character scalar naming the field to extract.
#'
#' @return The field value, or `NULL` when the field is unavailable.
#'
#' @noRd
#' @keywords internal
.extract_vcf_info_field <- function(vcf_info, field) {
  if (.is_missing_vcf_value(vcf_info)) {
    return(NULL)
  }

  if (
    is.list(vcf_info) &&
    is.null(names(vcf_info)) &&
    length(vcf_info) == 1L &&
    (
      is.list(vcf_info[[1L]]) ||
      is.data.frame(vcf_info[[1L]])
    )
  ) {
    vcf_info <- vcf_info[[1L]]
  }

  if (
    is.null(names(vcf_info)) ||
    !(field %in% names(vcf_info))
  ) {
    return(NULL)
  }

  value <- vcf_info[[field]]

  if (.is_missing_vcf_value(value)) {
    NULL
  } else {
    value
  }
}


#' Parse a numeric VCF INFO vector
#'
#' Parses numeric values represented as atomic vectors or comma-separated
#' character strings. Input order and duplicated values are preserved.
#'
#' @param value A potentially nested VCF INFO value.
#'
#' @return An unnamed numeric vector. Missing or empty inputs return
#'   `NA_real_`.
#'
#' @noRd
#' @keywords internal
.parse_vcf_numeric_vector <- function(value) {
  if (.is_missing_vcf_value(value)) {
    return(NA_real_)
  }

  values <- unlist(
    value,
    recursive = TRUE,
    use.names = FALSE
  )

  values <- unlist(
    strsplit(
      as.character(values),
      split = ",",
      fixed = TRUE
    ),
    recursive = FALSE,
    use.names = FALSE
  )

  values <- trimws(values)
  values[
    is.na(values) |
      !nzchar(values) |
      values == "."
  ] <- NA_character_

  parsed <- suppressWarnings(as.numeric(values))

  if (length(parsed) == 0L) {
    NA_real_
  } else {
    unname(parsed)
  }
}


#' Parse a scalar numeric VCF INFO value
#'
#' @param value A potentially nested VCF INFO value.
#'
#' @return A numeric scalar. Only the first parsed value is returned.
#'
#' @noRd
#' @keywords internal
.parse_vcf_numeric_scalar <- function(value) {
  parsed <- .parse_vcf_numeric_vector(value)

  if (length(parsed) == 0L || is.na(parsed[[1L]])) {
    NA_real_
  } else {
    parsed[[1L]]
  }
}


#' Construct site identifiers
#'
#' Constructs identifiers in the form
#' `chromosome-position-refAllele-alt1,alt2`.
#'
#' @param chromosome Character vector of chromosome labels.
#' @param position Integer vector of genomic positions.
#' @param ref_allele Character vector of reference alleles.
#' @param alt_alleles List-like object containing alternative alleles.
#'
#' @return A character vector with one identifier per site.
#'
#' @noRd
#' @keywords internal
.make_ica_site_id <- function(
    chromosome,
    position,
    ref_allele,
    alt_alleles
) {
  alt_text <- vapply(
    alt_alleles,
    function(value) {
      value <- as.character(
        unlist(
          value,
          recursive = TRUE,
          use.names = FALSE
        )
      )

      paste(value, collapse = ",")
    },
    character(1L)
  )

  paste(
    chromosome,
    position,
    ref_allele,
    alt_text,
    sep = "-"
  )
}
