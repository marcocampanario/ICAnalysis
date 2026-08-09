#' Collapse unique nested character values
#'
#' @param value A potentially nested object.
#' @param separator Character scalar used to join unique values.
#'
#' @return A character scalar or `NA_character_`.
#'
#' @noRd
#'
#' @keywords internal
.collapse_unique_character <- function(value, separator = ",") {
  if (.is_missing_vcf_value(value)) {
    return(NA_character_)
  }

  values <- trimws(
    as.character(
      unlist(
        value,
        recursive = TRUE,
        use.names = FALSE
      )
    )
  )

  values <- unique(
    values[
      !is.na(values) &
        nzchar(values) &
        values != "."
    ]
  )

  if (length(values) == 0L) {
    NA_character_
  } else {
    paste(values, collapse = separator)
  }
}


#' Extract a named field from a nested annotation
#'
#' @param value A potentially nested annotation object.
#' @param field Character scalar naming the field.
#'
#' @return The extracted field or `NULL`.
#'
#' @noRd
#'
#' @keywords internal
.extract_nested_field <- function(value, field) {
  if (is.null(value)) {
    return(NULL)
  }

  if (
    (is.list(value) || is.data.frame(value)) &&
    !is.null(names(value)) &&
    field %in% names(value)
  ) {
    return(value[[field]])
  }

  if (!is.list(value)) {
    return(NULL)
  }

  nested <- lapply(
    value,
    .extract_nested_field,
    field = field
  )

  keep <- !vapply(
    nested,
    .is_missing_vcf_value,
    logical(1L)
  )

  nested <- nested[keep]

  if (length(nested) == 0L) {
    NULL
  } else {
    unlist(
      nested,
      recursive = TRUE,
      use.names = FALSE
    )
  }
}


#' Extract a named annotation block
#'
#' @param variant A parsed ICA variant annotation.
#' @param annotation_name Character scalar naming the annotation.
#'
#' @return The annotation block or `NULL`.
#'
#' @noRd
#'
#' @keywords internal
.extract_annotation_block <- function(
    variant,
    annotation_name
) {
  if (is.null(variant)) {
    return(NULL)
  }

  if (
    (is.list(variant) || is.data.frame(variant)) &&
    !is.null(names(variant)) &&
    annotation_name %in% names(variant)
  ) {
    return(variant[[annotation_name]])
  }

  if (!is.list(variant)) {
    return(NULL)
  }

  for (item in variant) {
    block <- .extract_annotation_block(
      item,
      annotation_name
    )

    if (!is.null(block)) {
      return(block)
    }
  }

  NULL
}


#' Extract the overall allele frequency from a population annotation
#'
#' @param variant A parsed ICA variant annotation.
#' @param annotation_name Character scalar naming an ICA population
#'   annotation.
#'
#' @return A numeric scalar containing `allAf`, or `NA_real_`.
#'
#' @noRd
#'
#' @keywords internal
.extract_population_all_af <- function(
    variant,
    annotation_name
) {
  block <- .extract_annotation_block(
    variant,
    annotation_name
  )

  values <- .extract_nested_field(
    block,
    "allAf"
  )

  if (.is_missing_vcf_value(values)) {
    return(NA_real_)
  }

  values <- suppressWarnings(
    as.numeric(as.character(values))
  )
  values <- values[is.finite(values)]

  if (length(values) == 0L) {
    NA_real_
  } else {
    mean(values)
  }
}


#' Select the preferred gnomAD frequency
#'
#' Selects gnomAD Genomes when it contains a finite value, including zero.
#' gnomAD Exomes is used only when the Genomes value is unavailable.
#'
#' @param genome_af Numeric vector of gnomAD Genomes frequencies.
#' @param exome_af Numeric vector of gnomAD Exomes frequencies.
#'
#' @return A numeric vector with the same length as the inputs.
#'
#' @noRd
#'
#' @keywords internal
.select_preferred_gnomad_af <- function(
    genome_af,
    exome_af
) {
  if (length(genome_af) != length(exome_af)) {
    stop(
      "`genome_af` and `exome_af` must have equal lengths.",
      call. = FALSE
    )
  }

  genome_af <- suppressWarnings(
    as.numeric(as.character(genome_af))
  )
  exome_af <- suppressWarnings(
    as.numeric(as.character(exome_af))
  )

  result <- rep(NA_real_, length(genome_af))

  genome_available <- is.finite(genome_af)
  result[genome_available] <- genome_af[genome_available]

  use_exome <- !genome_available & is.finite(exome_af)
  result[use_exome] <- exome_af[use_exome]

  result
}


#' Extract dbSNP identifiers
#'
#' @param variant A parsed ICA variant annotation.
#'
#' @return A comma-separated character scalar or `NA_character_`.
#'
#' @noRd
#'
#' @keywords internal
.extract_dbsnp_id <- function(variant) {
  .collapse_unique_character(
    .extract_nested_field(variant, "dbsnp"),
    separator = ","
  )
}


#' Normalize ICA variant annotations to individual records
#'
#' @param value A parsed ICA `variants` element.
#'
#' @return A list containing one annotation record per element.
#'
#' @noRd
#'
#' @keywords internal
.normalize_ica_variant_records <- function(value) {
  if (is.null(value) || length(value) == 0L) {
    return(list())
  }

  if (is.data.frame(value)) {
    if (nrow(value) == 0L) {
      return(list())
    }

    return(
      lapply(
        seq_len(nrow(value)),
        function(row_index) {
          value[row_index, , drop = FALSE]
        }
      )
    )
  }

  if (is.list(value)) {
    # A named list represents the fields of one variant record.
    if (!is.null(names(value))) {
      return(list(value))
    }

    records <- list()

    for (item in value) {
      records <- c(
        records,
        .normalize_ica_variant_records(item)
      )
    }

    return(records)
  }

  list(value)
}


#' Match ICA variant annotations to alternative alleles
#'
#' Matching first uses the ICA `vid` field, followed by an exact
#' `altAllele` match. Remaining records are matched positionally when their
#' number equals the number of unmatched alternative alleles.
#'
#' @param variants A parsed ICA `variants` element.
#' @param alt_alleles Character vector of VCF alternative alleles.
#' @param chromosome Character scalar containing the chromosome.
#' @param position Integer scalar containing the genomic position.
#' @param ref_allele Character scalar containing the reference allele.
#'
#' @return A list with one annotation record per alternative allele.
#'
#' @noRd
#'
#' @keywords internal
.match_ica_variants_to_alts <- function(
    variants,
    alt_alleles,
    chromosome,
    position,
    ref_allele
) {
  number_of_alts <- length(alt_alleles)

  if (number_of_alts == 0L) {
    return(list())
  }

  records <- .normalize_ica_variant_records(variants)
  matched <- rep(list(NULL), number_of_alts)

  if (length(records) == 0L) {
    return(matched)
  }

  first_character_field <- function(record, field) {
    value <- .extract_nested_field(record, field)

    if (.is_missing_vcf_value(value)) {
      return(NA_character_)
    }

    value <- trimws(
      as.character(
        unlist(
          value,
          recursive = TRUE,
          use.names = FALSE
        )
      )
    )

    value <- value[
      !is.na(value) &
        nzchar(value)
    ]

    if (length(value) == 0L) {
      NA_character_
    } else {
      value[[1L]]
    }
  }

  normalize_vid <- function(value) {
    if (is.na(value) || !nzchar(value)) {
      return(NA_character_)
    }

    value <- tolower(trimws(value))
    sub("^chr", "", value)
  }

  expected_vids <- paste(
    sub(
      "^chr",
      "",
      as.character(chromosome),
      ignore.case = TRUE
    ),
    position,
    ref_allele,
    alt_alleles,
    sep = "-"
  )
  expected_vids <- vapply(
    expected_vids,
    normalize_vid,
    character(1L)
  )

  record_vids <- vapply(
    records,
    first_character_field,
    character(1L),
    field = "vid"
  )
  record_vids <- vapply(
    record_vids,
    normalize_vid,
    character(1L)
  )

  record_alts <- vapply(
    records,
    first_character_field,
    character(1L),
    field = "altAllele"
  )

  used_records <- rep(FALSE, length(records))

  # Prefer the original ICA variant identifier.
  for (allele_index in seq_len(number_of_alts)) {
    hits <- which(
      !used_records &
        !is.na(record_vids) &
        record_vids == expected_vids[[allele_index]]
    )

    if (length(hits) == 1L) {
      matched[[allele_index]] <- records[[hits]]
      used_records[[hits]] <- TRUE
    }
  }

  # Use an exact ALT match for records not matched by vid.
  for (allele_index in seq_len(number_of_alts)) {
    if (!is.null(matched[[allele_index]])) {
      next
    }

    hits <- which(
      !used_records &
        !is.na(record_alts) &
        toupper(record_alts) ==
        toupper(alt_alleles[[allele_index]])
    )

    if (length(hits) == 1L) {
      matched[[allele_index]] <- records[[hits]]
      used_records[[hits]] <- TRUE
    }
  }

  # ICA normally preserves ALT order. Use it as the final fallback.
  unmatched_alts <- which(
    vapply(matched, is.null, logical(1L))
  )
  unused_records <- which(!used_records)

  if (
    length(unmatched_alts) > 0L &&
    length(unmatched_alts) == length(unused_records)
  ) {
    for (index in seq_along(unmatched_alts)) {
      matched[[unmatched_alts[[index]]]] <-
        records[[unused_records[[index]]]]
    }
  }

  matched
}
