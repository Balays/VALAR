###
###

.pkgname <- "BSgenome.PRV.MdBio.1.0"

.seqnames <- NULL

.circ_seqs <- character(0)

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Suid alphaherpesvirus 1",
        common_name="PRV",
        genome="LT934125.1",
        provider="MDBIO",
        release_date="2017-05-22",
        source_url="-- information not available --",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "PRV"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

