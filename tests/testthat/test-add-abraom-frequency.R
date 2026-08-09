make_abraom_test_object <- function(n = 3L) {
  indexes <- seq_len(n)

  alleles <- S4Vectors::DataFrame(
    site_index = as.integer(indexes),
    allele_index = rep(1L, n),
    chromosome = paste0("chr", indexes),
    position = as.integer(indexes * 100L),
    refAllele = rep("A", n),
    altAllele = rep("C", n),
    variant_id = paste0(
      "chr",
      indexes,
      "-",
      indexes * 100L,
      "-A-C"
    ),
    af_gnomad_genome = rep(NA_real_, n),
    af_gnomad_exome = rep(NA_real_, n),
    af_gnomad = rep(NA_real_, n),
    af_1000g = rep(NA_real_, n),
    af_abraom = rep(NA_real_, n),
    dbsnp_id = rep(NA_character_, n)
  )

  ICAAlleles(
    alleles = alleles,
    header = list(
      annotator = "Illumina Connected Annotations 3.23.0",
      genomeAssembly = "GRCh38"
    ),
    metadata = list(
      test_marker = "retained"
    )
  )
}


test_that("ABRaOM frequencies are assigned by row order", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "Chr\tStart\tRef\tAlt\tabraom_freq",
      "chr22\t999\tG\tT\t0.25",
      "chr18\t888\tC\tA\t.",
      "chr7\t777\tT\tG\t0"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()
  original_variant_ids <- allele_data(object)$variant_id

  observed <- add_abraom_frequency(
    object,
    annovar_file
  )

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

  expect_identical(
    observed_data$variant_id,
    original_variant_ids
  )

  expect_equal(
    observed_data$af_abraom,
    c(0.25, NA_real_, 0)
  )
})


test_that("ANNOVAR coordinates are not used for assignment", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "Chr\tStart\tabraom_freq",
      "chr99\t999999\t0.30",
      "chr88\t888888\t0.20",
      "chr77\t777777\t0.10"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()

  observed <- add_abraom_frequency(
    object,
    annovar_file
  )

  expect_equal(
    allele_data(observed)$af_abraom,
    c(0.30, 0.20, 0.10)
  )
})


test_that("ABRaOM assignment preserves header and metadata", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "abraom_freq",
      "0.1",
      "0.2",
      "0.3"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()

  original_header <- ica_header(object)
  original_metadata <- ica_metadata(object)

  observed <- add_abraom_frequency(
    object,
    annovar_file
  )

  expect_identical(
    ica_header(observed),
    original_header
  )

  expect_identical(
    ica_metadata(observed),
    original_metadata
  )
})


test_that("existing ABRaOM frequencies are replaced", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "abraom_freq",
      "0.01",
      "0.02",
      "0.03"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()
  object_data <- allele_data(object)

  object_data$af_abraom <- c(0.9, 0.8, 0.7)
  allele_data(object) <- object_data

  observed <- add_abraom_frequency(
    object,
    annovar_file
  )

  expect_equal(
    allele_data(observed)$af_abraom,
    c(0.01, 0.02, 0.03)
  )
})


test_that("a custom ANNOVAR column can be selected", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "ABRAOM_AF\tanother_column",
      "0.11\t1",
      ".\t2",
      "0.33\t3"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()

  observed <- add_abraom_frequency(
    object,
    annovar_file,
    column = "ABRAOM_AF"
  )

  expect_equal(
    allele_data(observed)$af_abraom,
    c(0.11, NA_real_, 0.33)
  )
})


test_that("ABRaOM assignment requires identical row counts", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "abraom_freq",
      "0.1",
      "0.2"
    ),
    annovar_file
  )

  object <- make_abraom_test_object(
    n = 3L
  )

  expect_error(
    add_abraom_frequency(
      object,
      annovar_file
    ),
    "identical row counts"
  )
})


test_that("missing ANNOVAR columns are reported", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "another_frequency",
      "0.1",
      "0.2",
      "0.3"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()

  expect_error(
    add_abraom_frequency(
      object,
      annovar_file
    ),
    "was not found"
  )
})


test_that("duplicated ANNOVAR frequency columns are rejected", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "abraom_freq\tabraom_freq",
      "0.1\t0.2",
      "0.2\t0.3",
      "0.3\t0.4"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()

  expect_error(
    add_abraom_frequency(
      object,
      annovar_file
    ),
    "occurs more than once"
  )
})


test_that("invalid ABRaOM values are not silently converted", {
  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "abraom_freq",
      "0.1",
      "not-a-frequency",
      "0.3"
    ),
    annovar_file
  )

  object <- make_abraom_test_object()

  expect_error(
    add_abraom_frequency(
      object,
      annovar_file
    ),
    "non-numeric or non-finite"
  )
})


test_that("missing ANNOVAR files are reported", {
  object <- make_abraom_test_object()

  missing_file <- tempfile(
    fileext = ".txt"
  )

  expect_error(
    add_abraom_frequency(
      object,
      missing_file
    ),
    "does not exist"
  )
})


test_that("add_abraom_frequency requires ICAAlleles", {
  object <- S4Vectors::DataFrame(
    af_abraom = NA_real_
  )

  annovar_file <- tempfile(fileext = ".txt")
  on.exit(unlink(annovar_file), add = TRUE)

  writeLines(
    c(
      "abraom_freq",
      "0.1"
    ),
    annovar_file
  )

  expect_error(
    add_abraom_frequency(
      object,
      annovar_file
    ),
    "must be an ICAAlleles object"
  )
})
