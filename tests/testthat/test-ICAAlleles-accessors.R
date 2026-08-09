make_accessor_test_object <- function() {
  alleles <- S4Vectors::DataFrame(
    site_index = c(1L, 1L),
    allele_index = c(1L, 2L),
    chromosome = c("chr1", "chr1"),
    position = c(100L, 100L),
    refAllele = c("A", "A"),
    altAllele = c("C", "G"),
    variant_id = c(
      "chr1-100-A-C",
      "chr1-100-A-G"
    ),
    af_gnomad_genome = c(0.01, NA_real_),
    af_gnomad_exome = c(0.02, NA_real_),
    af_gnomad = c(0.01, NA_real_),
    af_1000g = c(NA_real_, NA_real_),
    af_abraom = c(NA_real_, NA_real_),
    dbsnp_id = c("rs123", NA_character_)
  )

  ICAAlleles(
    alleles = alleles,
    header = list(
      annotator = "Illumina Connected Annotations 3.23.0",
      genomeAssembly = "GRCh38"
    ),
    metadata = list(
      source = "ICAPositions"
    )
  )
}


test_that("allele_data returns the allele DataFrame", {
  object <- make_accessor_test_object()
  observed <- allele_data(object)

  expect_true(methods::is(observed, "DataFrame"))
  expect_identical(nrow(observed), 2L)
  expect_identical(
    observed$variant_id,
    c("chr1-100-A-C", "chr1-100-A-G")
  )
})


test_that("ica_header returns the retained ICA header", {
  object <- make_accessor_test_object()
  observed <- ica_header(object)

  expect_identical(
    observed$genomeAssembly,
    "GRCh38"
  )

  expect_identical(
    observed$annotator,
    "Illumina Connected Annotations 3.23.0"
  )
})


test_that("ica_metadata returns ICAAlleles metadata", {
  object <- make_accessor_test_object()
  observed <- ica_metadata(object)

  expect_identical(
    observed$source,
    "ICAPositions"
  )
})


test_that("allele_data replacement modifies and validates the object", {
  object <- make_accessor_test_object()
  replacement <- allele_data(object)

  replacement$af_abraom <- c(0.05, 0.10)
  allele_data(object) <- replacement

  expect_equal(
    allele_data(object)$af_abraom,
    c(0.05, 0.10)
  )

  expect_identical(
    ica_header(object)$genomeAssembly,
    "GRCh38"
  )

  expect_identical(
    ica_metadata(object)$source,
    "ICAPositions"
  )

  expect_true(
    methods::validObject(object, test = TRUE)
  )
})


test_that("allele_data replacement removes row names", {
  object <- make_accessor_test_object()
  replacement <- allele_data(object)

  rownames(replacement) <- c("first", "second")
  allele_data(object) <- replacement

  expect_null(
    rownames(allele_data(object))
  )
})


test_that("allele_data replacement rejects invalid schemas", {
  object <- make_accessor_test_object()
  invalid <- allele_data(object)

  invalid$af_gnomad_genome <- NULL

  expect_error(
    allele_data(object) <- invalid,
    "must contain exactly one column"
  )
})
