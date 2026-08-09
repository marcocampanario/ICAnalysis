test_that("the empty ICAAlleles prototype is valid", {
  object <- methods::new("ICAAlleles")

  expect_true(methods::is(object, "ICAAlleles"))
  expect_true(methods::is(object@alleles, "DataFrame"))
  expect_identical(nrow(object@alleles), 0L)
  expect_true(methods::validObject(object, test = TRUE))

  expect_identical(
    names(object@alleles),
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


test_that("ICAAlleles accepts multiple ALT alleles from one site", {
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
    dbsnp_id = c("rs1", NA_character_)
  )

  object <- methods::new(
    "ICAAlleles",
    alleles = alleles,
    header = list(genomeAssembly = "GRCh38"),
    metadata = list(source = "ICAPositions")
  )

  expect_true(methods::validObject(object, test = TRUE))
  expect_identical(nrow(object@alleles), 2L)
})


test_that("ICAAlleles rejects incomplete allele schemas", {
  alleles <- S4Vectors::DataFrame(
    site_index = 1L,
    allele_index = 1L
  )

  expect_error(
    methods::new(
      "ICAAlleles",
      alleles = alleles
    ),
    "must contain exactly one column"
  )
})


test_that("ICAAlleles rejects invalid indexes and positions", {
  alleles <- S4Vectors::DataFrame(
    site_index = 0L,
    allele_index = 1L,
    chromosome = "chr1",
    position = 100L,
    refAllele = "A",
    altAllele = "C",
    variant_id = "chr1-100-A-C",
    af_gnomad_genome = NA_real_,
    af_gnomad_exome = NA_real_,
    af_gnomad = NA_real_,
    af_1000g = NA_real_,
    af_abraom = NA_real_,
    dbsnp_id = NA_character_
  )

  expect_error(
    methods::new(
      "ICAAlleles",
      alleles = alleles
    ),
    "site_index.*positive"
  )
})
