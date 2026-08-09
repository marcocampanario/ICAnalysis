example_file <- system.file(
    "extdata",
    "subset_100_af.json.gz",
    package = "AnnoJSON"
)

parsed <- AnnoJSON::read_ica_json(
    example_file,
    verbose = FALSE
)

x <- AnnoJSON::ICAPositions(
    positions = S4Vectors::DataFrame(
        parsed$tbl_positions,
        check.names = FALSE
    ),
    header = parsed$header,
    metadata = list(
        source_file = basename(example_file)
    )
)

prepared <- prepare_ica_positions(x)
names(AnnoJSON::position_data(prepared))

show(prepared)
positions <- position_data(prepared)

### IF ABRAOM NEEDED

split <- split_ica_alleles(prepared)
alleles <- allele_data(split)
alleles$af_abraom <- abraom_frequency
allele_data(x) <- alleles

###

prepared <- annotate_variant_novelty(prepared)
