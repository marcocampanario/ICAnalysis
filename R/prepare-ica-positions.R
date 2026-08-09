#' Prepare ICA position records for analysis
#'
#' Extracts site-level allele-frequency and missingness fields from the
#' `vcfInfo` column of an `ICAPositions` object. The original ICA position
#' columns are retained, including `vcfInfo`.
#'
#' @param x An `ICAPositions` object created by
#'   [AnnoJSON::ICAPositions()].
#'
#' @details
#' `AF`, `AF_EXCL`, and `AF_used` are stored as numeric list-columns because
#' multiallelic sites can contain one frequency value per alternative allele.
#'
#' `AF_used` contains the complete `AF_EXCL` vector when `AF_EXCL` has at
#' least one non-missing value. Otherwise, it contains the original `AF`
#' vector. Selection is performed at the site level rather than separately
#' for each alternative allele.
#'
#' `site_id` has the format
#' `chromosome-position-refAllele-alt1,alt2`. Because site identifiers may
#' not be unique, `site_index` should be treated as the primary site key.
#'
#' Calling this function more than once recalculates and replaces the eight
#' preparation columns without removing other columns.
#'
#' @return An `ICAPositions` object containing the original position columns
#'   followed by `site_index`, `site_id`, `F_MISSING`, `AF`, `AF_EXCL`,
#'   `AF_used`, `uses_adjusted_af`, and `is_biallelic`.
#'
#' @examples
#' example_file <- system.file(
#'     "extdata",
#'     "subset_100_af.json.gz",
#'     package = "AnnoJSON"
#' )
#'
#' parsed <- AnnoJSON::read_ica_json(
#'     example_file,
#'     verbose = FALSE
#' )
#'
#' x <- AnnoJSON::ICAPositions(
#'     positions = S4Vectors::DataFrame(
#'         parsed$tbl_positions,
#'         check.names = FALSE
#'     ),
#'     header = parsed$header,
#'     metadata = list(
#'         source_file = basename(example_file)
#'     )
#' )
#'
#' prepared <- prepare_ica_positions(x)
#' names(AnnoJSON::position_data(prepared))
#'
#' @importFrom AnnoJSON position_data position_data<-
#' @importFrom S4Vectors I
#'
#' @export
prepare_ica_positions <- function(x) {
  .assert_ica_positions(x)

  positions <- position_data(x)

  prepared_columns <- c(
    "site_index",
    "site_id",
    "F_MISSING",
    "AF",
    "AF_EXCL",
    "AF_used",
    "uses_adjusted_af",
    "is_biallelic"
  )

  # Remove previously calculated versions to make the function idempotent.
  positions <- positions[
    ,
    !(names(positions) %in% prepared_columns),
    drop = FALSE
  ]

  vcf_info <- positions$vcfInfo

  f_missing <- vapply(
    vcf_info,
    function(info) {
      .parse_vcf_numeric_scalar(
        .extract_vcf_info_field(
          info,
          "F_MISSING"
        )
      )
    },
    numeric(1L)
  )

  af <- lapply(
    vcf_info,
    function(info) {
      .parse_vcf_numeric_vector(
        .extract_vcf_info_field(
          info,
          "AF"
        )
      )
    }
  )

  af_excl <- lapply(
    vcf_info,
    function(info) {
      .parse_vcf_numeric_vector(
        .extract_vcf_info_field(
          info,
          "AF_EXCL"
        )
      )
    }
  )

  uses_adjusted_af <- vapply(
    af_excl,
    function(value) {
      any(!is.na(value))
    },
    logical(1L)
  )

  af_used <- Map(
    function(af_value, af_excl_value, use_adjusted) {
      if (use_adjusted) {
        af_excl_value
      } else {
        af_value
      }
    },
    af,
    af_excl,
    uses_adjusted_af
  )

  is_biallelic <- vapply(
    positions$altAlleles,
    function(value) {
      length(
        unlist(
          value,
          recursive = TRUE,
          use.names = FALSE
        )
      ) == 1L
    },
    logical(1L)
  )

  positions$site_index <- seq_len(nrow(positions))

  positions$site_id <- .make_ica_site_id(
    chromosome = positions$chromosome,
    position = positions$position,
    ref_allele = positions$refAllele,
    alt_alleles = positions$altAlleles
  )

  positions$F_MISSING <- f_missing
  positions$AF <- I(unname(af))
  positions$AF_EXCL <- I(unname(af_excl))
  positions$AF_used <- I(unname(af_used))
  positions$uses_adjusted_af <- uses_adjusted_af
  positions$is_biallelic <- is_biallelic

  position_data(x) <- positions
  x
}
