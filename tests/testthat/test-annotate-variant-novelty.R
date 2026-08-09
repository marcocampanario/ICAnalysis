test_that("variant novelty follows population AF and dbSNP status", {
  alleles <- S4Vectors::DataFrame(
    variant_id = paste0("variant", seq_len(8L)),
    af_gnomad_genome = c(
      0.1, 0, NA, NA, 0, NA, 0, 0
    ),
    af_gnomad_exome = c(
      NA, 0.2, NA, NA, 0, NA, 0, 0
    ),
    af_1000g = c(
      NA, NA, 0.3, NA, 0, NA, 0, 0.2
    ),
    af_abraom = c(
      NA, NA, NA, 0.4, NA, NA, 0, NA
    ),
    dbsnp_id = c(
      NA, NA, NA, NA, "rs5", NA, ".", "rs8"
    )
  )

  observed <- annotate_variant_novelty(alleles)

  expect_true(methods::is(observed, "DataFrame"))

  expect_identical(
    observed$variant_novelty,
    c(
      "known",
      "known",
      "known",
      "known",
      "absent",
      "novel",
      "novel",
      "known"
    )
  )
})


test_that("gnomAD Exomes is evaluated independently from Genomes", {
  alleles <- S4Vectors::DataFrame(
    af_gnomad_genome = 0,
    af_gnomad_exome = 0.00001,
    af_1000g = NA_real_,
    af_abraom = NA_real_,
    dbsnp_id = NA_character_
  )

  observed <- annotate_variant_novelty(alleles)

  expect_identical(
    observed$variant_novelty,
    "known"
  )
})


test_that("missing ABRaOM frequencies do not change novelty", {
  alleles <- S4Vectors::DataFrame(
    af_gnomad_genome = c(NA_real_, NA_real_),
    af_gnomad_exome = c(NA_real_, NA_real_),
    af_1000g = c(NA_real_, NA_real_),
    af_abraom = c(NA_real_, NA_real_),
    dbsnp_id = c("rs123", NA_character_)
  )

  observed <- annotate_variant_novelty(alleles)

  expect_identical(
    observed$variant_novelty,
    c("absent", "novel")
  )
})


test_that("existing novelty annotations are recalculated", {
  alleles <- S4Vectors::DataFrame(
    af_gnomad_genome = 0.1,
    af_gnomad_exome = NA_real_,
    af_1000g = NA_real_,
    af_abraom = NA_real_,
    dbsnp_id = NA_character_,
    variant_novelty = "incorrect"
  )

  observed <- annotate_variant_novelty(alleles)

  expect_identical(
    observed$variant_novelty,
    "known"
  )

  expect_equal(
    sum(names(observed) == "variant_novelty"),
    1L
  )
})


test_that("required novelty columns are validated", {
  alleles <- S4Vectors::DataFrame(
    af_gnomad_genome = NA_real_
  )

  expect_error(
    annotate_variant_novelty(alleles),
    "missing required column"
  )
})


test_that("population frequencies must be numeric", {
  alleles <- S4Vectors::DataFrame(
    af_gnomad_genome = ".",
    af_gnomad_exome = NA_real_,
    af_1000g = NA_real_,
    af_abraom = NA_real_,
    dbsnp_id = NA_character_
  )

  expect_error(
    annotate_variant_novelty(alleles),
    "must be numeric"
  )
})


test_that("annotate_variant_novelty requires a DataFrame", {
  alleles <- data.frame(
    af_gnomad_genome = NA_real_,
    af_gnomad_exome = NA_real_,
    af_1000g = NA_real_,
    af_abraom = NA_real_,
    dbsnp_id = NA_character_
  )

  expect_error(
    annotate_variant_novelty(alleles),
    "S4Vectors::DataFrame"
  )
})
