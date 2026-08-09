make_valid_allele_data <- function() {
  S4Vectors::DataFrame(
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
}


test_that("ICAAlleles constructs an empty valid object", {
  object <- ICAAlleles()

  expect_true(methods::is(object, "ICAAlleles"))
  expect_true(methods::validObject(object, test = TRUE))
  expect_identical(nrow(object@alleles), 0L)
})


test_that("ICAAlleles constructs an object from a DataFrame", {
  alleles <- make_valid_allele_data()

  object <- ICAAlleles(
    alleles = alleles,
    header = list(genomeAssembly = "GRCh38"),
    metadata = list(source = "ICAPositions")
  )

  expect_true(methods::is(object, "ICAAlleles"))
  expect_true(methods::validObject(object, test = TRUE))
  expect_identical(nrow(object@alleles), 2L)

  expect_identical(
    object@header$genomeAssembly,
    "GRCh38"
  )

  expect_identical(
    object@metadata$source,
    "ICAPositions"
  )
})


test_that("ICAAlleles converts base data frames", {
  alleles <- as.data.frame(
    make_valid_allele_data(),
    stringsAsFactors = FALSE
  )

  object <- ICAAlleles(alleles)

  expect_true(
    methods::is(object@alleles, "DataFrame")
  )
  expect_true(methods::validObject(object, test = TRUE))
})


test_that("ICAAlleles removes row names", {
  alleles <- make_valid_allele_data()
  rownames(alleles) <- c("first", "second")

  object <- ICAAlleles(alleles)

  expect_null(rownames(object@alleles))
})


test_that("ICAAlleles validates header and metadata", {
  expect_error(
    ICAAlleles(header = "GRCh38"),
    "`header` must be a list"
  )

  expect_error(
    ICAAlleles(metadata = "split"),
    "`metadata` must be a list"
  )
})


test_that("ICAAlleles rejects unsupported allele inputs", {
  expect_error(
    ICAAlleles(alleles = list(site_index = 1L)),
    "DataFrame.*data.frame.*NULL"
  )
})


test_that("ICAAlleles invokes class validity", {
  alleles <- make_valid_allele_data()
  alleles$site_index <- as.numeric(alleles$site_index)

  expect_error(
    ICAAlleles(alleles),
    "site_index.*integer"
  )
})
