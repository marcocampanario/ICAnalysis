#' Construct an ICAAlleles object
#'
#' Creates an `ICAAlleles` object from allele-level records while preserving
#' the corresponding ICA header and analysis metadata.
#'
#' @details
#' The input must follow the allele-level schema defined by the
#' `ICAAlleles` class. A base `data.frame` is converted to an
#' `S4Vectors::DataFrame` before validation.
#'
#' Row names are removed because allele records are identified explicitly by
#' `site_index` and `allele_index`.
#'
#' This constructor does not split ICA position records. Use
#' [split_ica_alleles()] to construct allele records from an
#' `ICAPositions` object.
#'
#' @param alleles An `S4Vectors::DataFrame` or base `data.frame` containing
#'   one row per alternative allele. When `NULL`, an empty object with the
#'   complete allele schema is created.
#' @param header A list containing the original ICA header.
#' @param metadata A list containing provenance or analysis metadata.
#'
#' @return A valid `ICAAlleles` object.
#'
#' @examples
#' empty_object <- ICAAlleles()
#'
#' allele_data <- S4Vectors::DataFrame(
#'     site_index = 1L,
#'     allele_index = 1L,
#'     chromosome = "chr1",
#'     position = 100L,
#'     refAllele = "A",
#'     altAllele = "C",
#'     variant_id = "chr1-100-A-C",
#'     af_gnomad_genome = 0.01,
#'     af_gnomad_exome = 0.02,
#'     af_gnomad = 0.01,
#'     af_1000g = NA_real_,
#'     af_abraom = NA_real_,
#'     dbsnp_id = "rs123"
#' )
#'
#' object <- ICAAlleles(
#'     alleles = allele_data,
#'     header = list(genomeAssembly = "GRCh38"),
#'     metadata = list(source = "ICAPositions")
#' )
#'
#' object
#'
#' @importFrom S4Vectors DataFrame
#' @export
ICAAlleles <- function(
    alleles = NULL,
    header = list(),
    metadata = list()) {
  if (!is.list(header)) {
    stop(
      "`header` must be a list.",
      call. = FALSE
    )
  }

  if (!is.list(metadata)) {
    stop(
      "`metadata` must be a list.",
      call. = FALSE
    )
  }

  if (is.null(alleles)) {
    return(
      methods::new(
        "ICAAlleles",
        header = header,
        metadata = metadata
      )
    )
  }

  if (!methods::is(alleles, "DataFrame")) {
    if (!is.data.frame(alleles)) {
      stop(
        paste0(
          "`alleles` must be an S4Vectors::DataFrame, ",
          "a base data.frame, or NULL."
        ),
        call. = FALSE
      )
    }

    alleles <- DataFrame(
      alleles,
      check.names = FALSE
    )
  }

  rownames(alleles) <- NULL

  methods::new(
    "ICAAlleles",
    alleles = alleles,
    header = header,
    metadata = metadata
  )
}
