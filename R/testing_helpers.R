## These functions are meant to facilitate testing only

# This data.frame the options for the `type` parameter in
# `makeExampleDataSet` in its `$type` column, as well as other values required
# to create and verify the object that should be created.
.ds.info <- read.table(textConnection("
type                    class          package
DESeq2                  DESeqDataSet   DESeq2
edgeR                   DGEList        edgeR
limma                   EList          edgeR
voom                    EList          edgeR
voomWithQualityWeights  EList          edgeR
ExpressionSet           eSet           Biobase
"), stringsAsFactors=FALSE, header=TRUE)

# Makes example datasets for different objects we tidy.
#
# This is simply a wrapper to the \code{\link[DESeq2]{makeExampleDESeqDataSet}}
# function, which then transforms the output to the appropriate type.

# @param type The type of object to output. So far we can build objects for
#   \itemize{
#   \item{\code{"DESeq2"}}{
#         A \code{\link[DESeq2]{DESeqDataSet}}
#     }
#     \item{\code{"edgeR"}}{
#         A \code{\link[edgeR]{DGEList}}
#     }
#     \item{\code{"limma"}}{
#         An \code{\link[limma]{EList}}
#     }
#     \item{\code{"voom"}}{
#         An \code{\link[limma]{EList}} which results from passing a
#         \code{\link[edgeR]{DGEList}} through \code{\link[limma]{voom}}
#     }
#     \item{\code{"voomWithQualityWeights"}}{
#         Same as voom, but run through
#         \code{\link[limma]{voomWithQualityWeights}}
#     }
#     \item{\code{"ExpressionSet"}}{
#       An \code{\link[Biobase]{ExpressionSet}}
#     }
#   }
# @return The desired output object
#' @importFrom stats model.matrix
makeExampleDataSet <-
    function(type = 'DESeq2',
             n = 1000, m = 12, betaSD = 0, interceptMean = 4,
             interceptSD = 2,
             dispMeanRel = function(x) 4/x + 0.1, sizeFactors = rep(1, m)) {
    type <- match.arg(type, .ds.info$type)
    ds.info <- .ds.info[.ds.info$type == type,]
    if (nrow(ds.info) != 1L) {
        stop("Unrecognized data type: ", type ,".\n  Available options are: ",
             paste(.ds.info$type, collapse=',', sep=' '))
    }
    if (!requireNamespace("DESeq2")) {
        stop("DESeq2 required to create dataset")
    }
    pkg <- ds.info$package
    if (type != 'DESeq2' && !requireNamespace(pkg, character.only=TRUE)) {
        stop(pkg, " package required to create data of this type")
    }

    dds <- DESeq2::makeExampleDESeqDataSet(n = n, m = m, betaSD = betaSD,
                                           interceptMean = interceptMean,
                                           interceptSD = interceptSD,
                                           dispMeanRel = dispMeanRel)
    dds <- DESeq2::estimateSizeFactors(dds)
    xprs <- log2(DESeq2::counts(dds, normalized=TRUE) + 1)
    pd <- as.data.frame(SummarizedExperiment::colData(dds))
    fd <- data.frame(
        gene=rownames(dds),
        seqnames=as.character(SummarizedExperiment::seqnames(SummarizedExperiment::rowRanges(dds))),
        start=SummarizedExperiment::start(SummarizedExperiment::rowRanges(dds)),
        end=SummarizedExperiment::end(SummarizedExperiment::rowRanges(dds)),
        stringsAsFactors=FALSE)

    if (type == 'DESeq2') {
        out <- dds
    } else if (type %in% c('edgeR', 'limma', 'voom', 'voomWithQualityWeights')) {
        out <- edgeR::DGEList(DESeq2::counts(dds))
        out <- edgeR::calcNormFactors(out)
        out$samples <- cbind(out$samples, pd)
        out$genes <- fd
        if (type != 'edgeR') {
            mm <- model.matrix(DESeq2::design(dds), out$samples)
            if (type == 'voomWithQualityWeights') {
                out <- limma::voomWithQualityWeights(out, mm, plot=FALSE)
            } else {
                out <- limma::voom(out, mm, plot=FALSE)
            }
            if (type == 'limma') {
                out$weights <- NULL
                out$design <- NULL
                out$targets$group <- NULL
            }
        }
    } else if (type == 'ExpressionSet') {
        out <- ExpressionSet(xprs)
        Biobase::pData(out) <- pd
        Biobase::fData(out) <- fd
    } else {
        stop(type, " not yet accounted for")
    }

    if (!inherits(out, ds.info$class)) {
        warning("The output object is a ", class(out), " but we expected ",
                ds.info$class)
    }

    out
}


