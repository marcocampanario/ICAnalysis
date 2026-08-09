make_test_ica_positions <- function(prepared = TRUE) {
  example_file <- system.file(
    "extdata",
    "subset_100_af.json.gz",
    package = "AnnoJSON"
  )

  if (!nzchar(example_file)) {
    stop(
      "The AnnoJSON example file was not found.",
      call. = FALSE
    )
  }

  parsed <- AnnoJSON::read_ica_json(
    example_file,
    workers = 1L,
    progress = FALSE,
    verbose = FALSE
  )

  object <- AnnoJSON::ICAPositions(
    positions = S4Vectors::DataFrame(
      parsed$tbl_positions,
      check.names = FALSE
    ),
    header = parsed$header,
    metadata = list(
      test_marker = "retained"
    )
  )

  if (prepared) {
    object <- prepare_ica_positions(object)
  }

  object
}


make_multiallelic_test_positions <- function() {
  variant_records <- data.frame(
    vid = c(
      "1-100-A-G",
      "1-100-A-T",
      "1-100-A-C"
    ),
    chromosome = rep("chr1", 3L),
    begin = rep(100L, 3L),
    end = rep(100L, 3L),
    refAllele = rep("A", 3L),
    altAllele = c("G", "T", "C"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  variant_records$gnomad <- I(
    list(
      NULL,
      data.frame(allAf = 0),
      data.frame(allAf = 0.10)
    )
  )

  variant_records[["gnomad-exome"]] <- I(
    list(
      data.frame(allAf = 0.40),
      data.frame(allAf = 0.60),
      data.frame(allAf = 0.20)
    )
  )

  variant_records$oneKg <- I(
    list(
      data.frame(allAf = 0.03),
      NULL,
      data.frame(allAf = 0.05)
    )
  )

  variant_records$dbsnp <- I(
    list(
      "rsG",
      "rsT",
      "rsC"
    )
  )

  positions <- S4Vectors::DataFrame(
    chromosome = "chr1",
    position = 100L,
    refAllele = "A",
    altAlleles = S4Vectors::I(
      list(c("C", "G", "T"))
    ),
    vcfInfo = S4Vectors::I(
      list(list())
    ),
    quality = 99,
    filters = "PASS",
    cytogeneticBand = NA_character_,
    samples = S4Vectors::I(
      list(list())
    ),
    variants = S4Vectors::I(
      list(variant_records)
    ),
    site_index = 1L,
    site_id = "chr1-100-A-C,G,T",
    F_MISSING = 0,
    AF = S4Vectors::I(
      list(c(0.10, 0.20, 0.30))
    ),
    AF_EXCL = S4Vectors::I(
      list(c(NA_real_, NA_real_, NA_real_))
    ),
    AF_used = S4Vectors::I(
      list(c(0.10, 0.20, 0.30))
    ),
    uses_adjusted_af = FALSE,
    is_biallelic = FALSE,
    check.names = FALSE
  )

  AnnoJSON::ICAPositions(
    positions = positions,
    header = list(
      annotator = "Illumina Connected Annotations 3.23.0",
      genomeAssembly = "GRCh38",
      samples = c("sample1", "sample2")
    ),
    metadata = list(
      test_marker = "multiallelic"
    )
  )
}


test_that("split_ica_alleles returns an ICAAlleles object", {
  object <- make_test_ica_positions()
  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_true(
    methods::is(observed, "ICAAlleles")
  )

  expect_true(
    methods::is(observed_data, "DataFrame")
  )

  expect_true(
    methods::validObject(observed, test = TRUE)
  )
})


test_that("split ICA allele data has the expected schema", {
  object <- make_test_ica_positions()
  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_identical(
    names(observed_data),
    c(
      "site_index",
      "allele_index",
      "chromosome",
      "position",
      "refAllele",
      "altAllele",
      "variant_id",
      "af_gnomad_genome",
      "af_gnomad_exome",
      "af_gnomad",
      "af_1000g",
      "af_abraom",
      "dbsnp_id"
    )
  )
})


test_that("split_ica_alleles creates one row per ALT allele", {
  object <- make_test_ica_positions()
  positions <- AnnoJSON::position_data(object)

  expected_alleles <- sum(
    lengths(positions$altAlleles)
  )

  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_equal(
    nrow(observed_data),
    expected_alleles
  )
})


test_that("site and ALT allele order are preserved", {
  object <- make_test_ica_positions()
  positions <- AnnoJSON::position_data(object)

  alt_counts <- lengths(positions$altAlleles)

  expected_site_index <- as.integer(
    rep.int(
      positions$site_index,
      alt_counts
    )
  )

  expected_allele_index <- as.integer(
    unname(sequence(alt_counts))
  )

  expected_alt_alleles <- as.character(
    unlist(
      positions$altAlleles,
      recursive = FALSE,
      use.names = FALSE
    )
  )

  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_identical(
    observed_data$site_index,
    expected_site_index
  )

  expect_identical(
    observed_data$allele_index,
    expected_allele_index
  )

  expect_identical(
    observed_data$altAllele,
    expected_alt_alleles
  )
})


test_that("variant identifiers preserve VCF allele representation", {
  object <- make_test_ica_positions()
  positions <- AnnoJSON::position_data(object)

  alt_counts <- lengths(positions$altAlleles)

  site_rows <- rep.int(
    seq_len(nrow(positions)),
    alt_counts
  )

  expected_variant_id <- paste(
    as.character(positions$chromosome[site_rows]),
    as.integer(positions$position[site_rows]),
    as.character(positions$refAllele[site_rows]),
    as.character(
      unlist(
        positions$altAlleles,
        recursive = FALSE,
        use.names = FALSE
      )
    ),
    sep = "-"
  )

  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_identical(
    observed_data$variant_id,
    expected_variant_id
  )
})


test_that("ABRaOM frequencies are initialized as numeric missing values", {
  object <- make_test_ica_positions()
  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_true(
    is.numeric(observed_data$af_abraom)
  )

  expect_true(
    all(is.na(observed_data$af_abraom))
  )
})


test_that("ICA header is preserved after allele splitting", {
  object <- make_test_ica_positions()
  observed <- split_ica_alleles(object)

  expect_identical(
    ica_header(observed),
    AnnoJSON::ica_header(object)
  )
})


test_that("original metadata and splitting provenance are retained", {
  object <- make_test_ica_positions()
  observed <- split_ica_alleles(object)

  observed_metadata <- ica_metadata(observed)
  observed_data <- allele_data(observed)

  expect_identical(
    observed_metadata$test_marker,
    "retained"
  )

  expect_identical(
    observed_metadata$split_ica_alleles$source_class,
    "ICAPositions"
  )

  expect_equal(
    observed_metadata$split_ica_alleles$n_sites,
    nrow(AnnoJSON::position_data(object))
  )

  expect_equal(
    observed_metadata$split_ica_alleles$n_alleles,
    nrow(observed_data)
  )
})


test_that("variant annotations are matched to the original ALT order", {
  object <- make_multiallelic_test_positions()

  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_identical(
    observed_data$altAllele,
    c("C", "G", "T")
  )

  expect_identical(
    observed_data$variant_id,
    c(
      "chr1-100-A-C",
      "chr1-100-A-G",
      "chr1-100-A-T"
    )
  )

  expect_equal(
    observed_data$af_gnomad_genome,
    c(0.10, NA_real_, 0)
  )

  expect_equal(
    observed_data$af_gnomad_exome,
    c(0.20, 0.40, 0.60)
  )

  expect_equal(
    observed_data$af_1000g,
    c(0.05, 0.03, NA_real_)
  )

  expect_identical(
    observed_data$dbsnp_id,
    c("rsC", "rsG", "rsT")
  )
})


test_that("gnomAD Genomes is preferred and Exomes is the fallback", {
  object <- make_multiallelic_test_positions()

  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  expect_equal(
    observed_data$af_gnomad,
    c(
      0.10,
      0.40,
      0
    )
  )
})


test_that("a zero gnomAD Genomes AF is retained over positive Exomes AF", {
  object <- make_multiallelic_test_positions()

  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  t_allele <- observed_data$altAllele == "T"

  expect_equal(
    observed_data$af_gnomad_genome[t_allele],
    0
  )

  expect_equal(
    observed_data$af_gnomad_exome[t_allele],
    0.60
  )

  expect_equal(
    observed_data$af_gnomad[t_allele],
    0
  )
})


test_that("allele indexes restart within each source site", {
  object <- make_test_ica_positions()
  observed <- split_ica_alleles(object)
  observed_data <- allele_data(observed)

  indexes_by_site <- split(
    observed_data$allele_index,
    observed_data$site_index
  )

  expected_indexes <- lapply(
    indexes_by_site,
    function(indexes) {
      seq_along(indexes)
    }
  )

  expect_identical(
    unname(indexes_by_site),
    unname(expected_indexes)
  )
})


test_that("split_ica_alleles rejects unprepared ICAPositions objects", {
  object <- make_test_ica_positions(
    prepared = FALSE
  )

  expect_error(
    split_ica_alleles(object),
    "has not been prepared"
  )
})


test_that("split_ica_alleles rejects unsupported inputs", {
  object <- S4Vectors::DataFrame(
    chromosome = "chr1",
    position = 100L
  )

  expect_error(
    split_ica_alleles(object),
    "ICAPositions"
  )
})
