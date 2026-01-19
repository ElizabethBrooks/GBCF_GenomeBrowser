###
###

.pkgname <- "BSgenome.Dmagna.LRV0"

.seqnames <- c("LG1", "LG10", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "mitochondrion", "scaffold_11", "scaffold_12", "scaffold_13", "scaffold_14", "scaffold_15", "scaffold_16", "scaffold_17", "scaffold_18", "scaffold_19", "scaffold_20", "scaffold_21", "scaffold_22", "scaffold_23", "scaffold_24", "scaffold_25", "scaffold_26", "scaffold_27", "scaffold_28", "scaffold_29", "scaffold_30", "scaffold_31", "scaffold_32", "scaffold_33", "scaffold_34", "scaffold_35", "scaffold_36", "scaffold_37", "scaffold_38",  "scaffold_39", "scaffold_40", "scaffold_41", "scaffold_42", "scaffold_44", "scaffold_45", "scaffold_46", "scaffold_47", "scaffold_48", "scaffold_49", "scaffold_50", "scaffold_51", "scaffold_52", "scaffold_53", "scaffold_54", "scaffold_55", "scaffold_56", "scaffold_57", "scaffold_58", "scaffold_59", "scaffold_60", "scaffold_61", "scaffold_62", "scaffold_63", "scaffold_64", "scaffold_65", "scaffold_66", "scaffold_67", "scaffold_68", "scaffold_69", "scaffold_70", "scaffold_71")

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
        organism="Daphnia magna",
        common_name="D. magna",
        genome="Dmagna.LRV0_1",
        provider="MEP",
        release_date="2023/04",
        source_url="n/a",
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

    old_objname <- "Dmagna"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

