#' Access ICA allele data
#'
#' Returns the allele-level `DataFrame` stored in an `ICAAlleles` object.
#'
#' @param x An `ICAAlleles` object.
#'
#' @return An `S4Vectors::DataFrame` containing one row per alternative
#'   allele.
#'
#' @examples
#' object <- ICAAlleles()
#' allele_data(object)
#'
#' @export
methods::setGeneric(
  "allele_data",
  function(x) {
    methods::standardGeneric("allele_data")
  }
)


#' @export
#' @noRd
methods::setMethod(
  "allele_data",
  signature(x = "ICAAlleles"),
  function(x) {
    x@alleles
  }
)

#' Replace ICA allele data
#'
#' Replaces the allele-level `DataFrame` stored in an `ICAAlleles` object.
#' The modified object is validated before it is returned.
#'
#' @param x An `ICAAlleles` object.
#' @param value An `S4Vectors::DataFrame` following the required
#'   `ICAAlleles` schema.
#'
#' @return The modified and validated `ICAAlleles` object.
#'
#' @examples
#' object <- ICAAlleles()
#' replacement <- allele_data(object)
#' allele_data(object) <- replacement
#'
#' @name allele_data-set
#' @export
methods::setGeneric(
  "allele_data<-",
  function(x, value) {
    methods::standardGeneric("allele_data<-")
  }
)


#' @export
#' @noRd
methods::setReplaceMethod(
  "allele_data",
  signature(
    x = "ICAAlleles",
    value = "DataFrame"
  ),
  function(x, value) {
    rownames(value) <- NULL
    x@alleles <- value

    methods::validObject(x)
    x
  }
)

#' @importFrom AnnoJSON ica_header
#' @export
#' @noRd
AnnoJSON::ica_header

#' Access the ICA header from an ICAAlleles object
#'
#' Returns the original Illumina Connected Annotations header retained when
#' the allele-level object was created.
#'
#' @param x An `ICAAlleles` object.
#'
#' @return A list containing the ICA header.
#'
#' @examples
#' object <- ICAAlleles(
#'     header = list(
#'         annotator = "Illumina Connected Annotations 3.23.0",
#'         genomeAssembly = "GRCh38"
#'     )
#' )
#'
#' ica_header(object)
#'
#' @name ica_header-ICAAlleles-method
#' @export
methods::setMethod(
  "ica_header",
  signature(x = "ICAAlleles"),
  function(x) {
    x@header
  }
)

#' @importFrom AnnoJSON ica_metadata
#' @export
#' @noRd
AnnoJSON::ica_metadata

#' Access ICAAlleles metadata
#'
#' Returns provenance and analysis metadata stored in an `ICAAlleles`
#' object.
#'
#' @param x An `ICAAlleles` object.
#'
#' @return A list containing the object metadata.
#'
#' @examples
#' object <- ICAAlleles(
#'     metadata = list(
#'         source = "ICAPositions",
#'         operation = "split alternative alleles"
#'     )
#' )
#'
#' ica_metadata(object)
#'
#' @name ica_metadata-ICAAlleles-method
#' @export
methods::setMethod(
  "ica_metadata",
  signature(x = "ICAAlleles"),
  function(x) {
    x@metadata
  }
)
