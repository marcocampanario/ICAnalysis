#' Split ICA positions into alternative-allele records
#'
#' Converts a prepared `ICAPositions` object into an `ICAAlleles` object
#' containing one row per alternative allele.
#'
#' @param x An `ICAPositions` object processed with
#'   [prepare_ica_positions()].
#'
#' @details
#' Input site order and ALT-allele order are preserved. `site_index` links
#' each allele record to its original ICA position, while `allele_index`
#' identifies the ALT allele within that site.
#'
#' Both gnomAD Genomes and gnomAD Exomes frequencies are retained.
#' `af_gnomad` prefers the Genomes value whenever it is available and uses
#' Exomes only as a fallback. Consequently, `af_gnomad` must not be used
#' alone to determine population-database presence.
#'
#' `af_abraom` is initialized with `NA_real_`. It can subsequently be
#' populated by [add_abraom_frequency()].
#'
#' The original ICA header and metadata are retained. Transformation
#' provenance, including the number of source sites and resulting allele
#' records, is added to the output metadata.
#'
#' @return An `ICAAlleles` object. Its allele-level `DataFrame`, accessible
#'   with [allele_data()], contains the columns `site_index`, `allele_index`,
#'   `chromosome`, `position`, `refAllele`, `altAllele`, `variant_id`,
#'   `af_gnomad_genome`, `af_gnomad_exome`, `af_gnomad`, `af_1000g`,
#'   `af_abraom`, and `dbsnp_id`.
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
#'     header = parsed$header
#' )
#'
#' x <- prepare_ica_positions(x)
#' alleles <- split_ica_alleles(x)
#'
#' alleles
#' allele_data(alleles)
#'
#' @importFrom AnnoJSON position_data ica_header ica_metadata
#' @importFrom S4Vectors DataFrame
#' @export
split_ica_alleles <- function(x) {
  .assert_ica_positions(x)

  positions <- position_data(x)

  if (!"site_index" %in% names(positions)) {
    stop(
      paste0(
        "`x` has not been prepared. Run ",
        "`prepare_ica_positions(x)` before splitting ALT alleles."
      ),
      call. = FALSE
    )
  }

  alt_alleles <- lapply(
    positions$altAlleles,
    function(value) {
      as.character(
        unlist(
          value,
          recursive = TRUE,
          use.names = FALSE
        )
      )
    }
  )

  alt_counts <- lengths(alt_alleles)
  number_of_alleles <- sum(alt_counts)

  site_rows <- rep.int(
    seq_len(nrow(positions)),
    alt_counts
  )

  allele_index <- unname(sequence(alt_counts))

  matched_variants <- vector(
    mode = "list",
    length = number_of_alleles
  )

  output_index <- 0L

  for (site_row in seq_len(nrow(positions))) {
    number_of_alts <- alt_counts[[site_row]]

    if (number_of_alts == 0L) {
      next
    }

    indexes <- output_index + seq_len(number_of_alts)

    matched_variants[indexes] <-
      .match_ica_variants_to_alts(
        variants = positions$variants[[site_row]],
        alt_alleles = alt_alleles[[site_row]],
        chromosome = positions$chromosome[[site_row]],
        position = positions$position[[site_row]],
        ref_allele = positions$refAllele[[site_row]]
      )

    output_index <- output_index + number_of_alts
  }

  chromosome <- as.character(
    positions$chromosome[site_rows]
  )
  position <- as.integer(
    positions$position[site_rows]
  )
  ref_allele <- as.character(
    positions$refAllele[site_rows]
  )
  alt_allele <- as.character(
    unlist(
      alt_alleles,
      recursive = FALSE,
      use.names = FALSE
    )
  )

  af_gnomad_genome <- vapply(
    matched_variants,
    .extract_population_all_af,
    numeric(1L),
    annotation_name = "gnomad"
  )

  af_gnomad_exome <- vapply(
    matched_variants,
    .extract_population_all_af,
    numeric(1L),
    annotation_name = "gnomad-exome"
  )

  af_1000g <- vapply(
    matched_variants,
    .extract_population_all_af,
    numeric(1L),
    annotation_name = "oneKg"
  )

  dbsnp_id <- vapply(
    matched_variants,
    .extract_dbsnp_id,
    character(1L)
  )

  af_gnomad <- .select_preferred_gnomad_af(
    genome_af = af_gnomad_genome,
    exome_af = af_gnomad_exome
  )

  variant_id <- paste(
    chromosome,
    position,
    ref_allele,
    alt_allele,
    sep = "-"
  )

  allele_table <- DataFrame(
    site_index = as.integer(
      positions$site_index[site_rows]
    ),
    allele_index = as.integer(allele_index),
    chromosome = chromosome,
    position = position,
    refAllele = ref_allele,
    altAllele = alt_allele,
    variant_id = variant_id,
    af_gnomad_genome = af_gnomad_genome,
    af_gnomad_exome = af_gnomad_exome,
    af_gnomad = af_gnomad,
    af_1000g = af_1000g,
    af_abraom = rep(
      NA_real_,
      number_of_alleles
    ),
    dbsnp_id = dbsnp_id,
    check.names = FALSE
  )

  output_metadata <- ica_metadata(x)

  output_metadata$split_ica_alleles <- list(
    source_class = "ICAPositions",
    n_sites = nrow(positions),
    n_alleles = number_of_alleles
  )

  ICAAlleles(
    alleles = allele_table,
    header = ica_header(x),
    metadata = output_metadata
  )
}
