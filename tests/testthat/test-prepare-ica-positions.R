testthat::test_that(
  "prepare_ica_positions retains original ICA columns",
  {
    example_file <- system.file(
      "extdata",
      "subset_100_af.json.gz",
      package = "AnnoJSON"
    )

    parsed <- AnnoJSON::read_ica_json(
      example_file,
      verbose = FALSE
    )

    object <- AnnoJSON::ICAPositions(
      positions = S4Vectors::DataFrame(
        parsed$tbl_positions,
        check.names = FALSE
      ),
      header = parsed$header,
      metadata = list(
        source_file = basename(example_file)
      )
    )

    before <- AnnoJSON::position_data(object)
    original_names <- names(before)

    prepared <- prepare_ica_positions(object)
    after <- AnnoJSON::position_data(prepared)

    expected_new_columns <- c(
      "site_index",
      "site_id",
      "F_MISSING",
      "AF",
      "AF_EXCL",
      "AF_used",
      "uses_adjusted_af",
      "is_biallelic"
    )

    testthat::expect_s4_class(
      prepared,
      "ICAPositions"
    )
    testthat::expect_identical(
      names(after),
      c(original_names, expected_new_columns)
    )
    testthat::expect_identical(
      after[, original_names, drop = FALSE],
      before
    )
    testthat::expect_identical(
      AnnoJSON::ica_header(prepared),
      AnnoJSON::ica_header(object)
    )
    testthat::expect_identical(
      AnnoJSON::ica_metadata(prepared),
      AnnoJSON::ica_metadata(object)
    )
    testthat::expect_identical(
      methods::validObject(prepared, test = TRUE),
      TRUE
    )
  }
)


testthat::test_that(
  "site indexes and biallelic status are calculated correctly",
  {
    example_file <- system.file(
      "extdata",
      "subset_100_af.json.gz",
      package = "AnnoJSON"
    )

    parsed <- AnnoJSON::read_ica_json(
      example_file,
      verbose = FALSE
    )

    object <- AnnoJSON::ICAPositions(
      positions = S4Vectors::DataFrame(
        parsed$tbl_positions,
        check.names = FALSE
      ),
      header = parsed$header
    )

    prepared <- prepare_ica_positions(object)
    positions <- AnnoJSON::position_data(prepared)

    expected_biallelic <- vapply(
      positions$altAlleles,
      function(value) {
        length(unlist(value, use.names = FALSE)) == 1L
      },
      logical(1L)
    )

    testthat::expect_identical(
      positions$site_index,
      seq_len(nrow(positions))
    )
    testthat::expect_identical(
      positions$is_biallelic,
      expected_biallelic
    )
    testthat::expect_length(
      positions$F_MISSING,
      nrow(positions)
    )
    testthat::expect_length(
      positions$AF,
      nrow(positions)
    )
  }
)


testthat::test_that(
  "AF_used follows the site-level AF_EXCL rule",
  {
    example_file <- system.file(
      "extdata",
      "subset_100_af.json.gz",
      package = "AnnoJSON"
    )

    parsed <- AnnoJSON::read_ica_json(
      example_file,
      verbose = FALSE
    )

    object <- AnnoJSON::ICAPositions(
      positions = S4Vectors::DataFrame(
        parsed$tbl_positions,
        check.names = FALSE
      ),
      header = parsed$header
    )

    positions <- AnnoJSON::position_data(object)

    positions$altAlleles[[1L]] <- c("A", "T")

    info <- positions$vcfInfo[[1L]]
    info$F_MISSING <- "0.25"
    info$AF <- "0.10,0.10"
    info$AF_EXCL <- "0.20,0.30"
    positions$vcfInfo[[1L]] <- info

    AnnoJSON::position_data(object) <- positions

    prepared <- prepare_ica_positions(object)
    observed <- AnnoJSON::position_data(prepared)

    testthat::expect_identical(
      observed$F_MISSING[[1L]],
      0.25
    )
    testthat::expect_identical(
      observed$AF[[1L]],
      c(0.10, 0.10)
    )
    testthat::expect_identical(
      observed$AF_EXCL[[1L]],
      c(0.20, 0.30)
    )
    testthat::expect_identical(
      observed$AF_used[[1L]],
      c(0.20, 0.30)
    )
    testthat::expect_true(
      observed$uses_adjusted_af[[1L]]
    )
    testthat::expect_false(
      observed$is_biallelic[[1L]]
    )
  }
)


testthat::test_that(
  "prepare_ica_positions is idempotent",
  {
    example_file <- system.file(
      "extdata",
      "subset_100_af.json.gz",
      package = "AnnoJSON"
    )

    parsed <- AnnoJSON::read_ica_json(
      example_file,
      verbose = FALSE
    )

    object <- AnnoJSON::ICAPositions(
      positions = S4Vectors::DataFrame(
        parsed$tbl_positions,
        check.names = FALSE
      ),
      header = parsed$header
    )

    once <- prepare_ica_positions(object)
    twice <- prepare_ica_positions(once)

    testthat::expect_identical(
      AnnoJSON::position_data(once),
      AnnoJSON::position_data(twice)
    )
  }
)


testthat::test_that(
  "prepare_ica_positions requires ICAPositions",
  {
    testthat::expect_error(
      prepare_ica_positions(list()),
      "`x` must be an `ICAPositions` object.",
      fixed = TRUE
    )
  }
)
