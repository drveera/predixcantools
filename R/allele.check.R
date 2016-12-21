#' allele.check
#'
#' allele.check
#'
#' @param summary summary data frame
#' @param dbweights db weights produced by function extract_weights
#' @param snp.col snp column in summary
#' @param effect.allele.col effect allele column in summary
#' @param alt.allele.col alternative allele column in summary
#' @param beta.col beta column in summary; if you have odds ration, change it beta by log(odds ratio)
#'
#' @import data.table
#' 
#' @export
allele.check <- function(summary,
                         dbweights,
                         snp.col="rsid",
                         effect.allele.col="a1",
                         alt.allele.col="a2",
                         beta.col = "beta"){

    ##subset summary file

    cat("preprocessing summary file:")
    summary <- data.table(summary)
    if(!snp.col == 'rsid'){
        summary$rsid <- summary[,snp.col,with=FALSE]
    }

    if(! effect.allele.col == 'a1'){
        summary$a1 <- summary[,effect.allele.col,with=FALSE]
    }

    if(! alt.allele.col == 'a2'){
        summary$a2 <- summary[,alt.allele.col,with=FALSE]
    }

    if(! beta.col == 'beta'){
        summary$beta <- summary[,beta.col,with=FALSE]
    }

    ##remove duplicated in summary
    summary <- summary[!duplicated(summary$rsid),]
    cat("ok \n")
    cat("subsetting only the SNPs in weight:")
    dfm <- merge(summary,dbweights, by= "rsid")
    dfm <<- data.table(dfm)
    cat("ok \n")
    cat("processing for possible strand flips:")
    ## seperate variants with strand flips
    dfm.flip <- dfm[!(dfm$ref_allele == a2 | dfm$ref_allele == a1)]
    dfm.noflip <- dfm[(dfm$ref_allele == a2 | dfm$ref_allele == a1)]
    ##flip the strand
    dfm.flip$a1 <- sapply(dfm.flip$a1,flipstrand)
    dfm.flip$a2 <- sapply(dfm.flip$a2,flipstrand)
    ##take only those flipped and drop others
    dfm.flip <- dfm.flip[(ref_allele == a2 | ref_allele == a1) & (eff_allele == a2 | eff_allele == a1)]
    ##merge back
    dfm <- rbind(dfm.flip,dfm.noflip)
    cat("ok \n")
    cat("matching the reference and alt alleles with weights:")
    ##split the dfm with respect to allele match
    dfm$amatch <- with(dfm, a1 == eff_allele & a2 == ref_allele)
    dfm$amatch <- ifelse(dfm$amatch,"match","nomatch")
    ##split
    dfm.split <- split(dfm,dfm$amatch)
    matchdfm <- dfm.split$match
    nomatchdfm <- dfm.split$nomatch
    nomatchdfm$beta <- nomatchdfm$beta * -1
    a1 <- nomatchdfm$a2
    a2 <- nomatchdfm$a1
    nomatchdfm$a1 <- a1
    nomatchdfm$a2 <- a2
    cat("ok \n")
    ##combine back
    dfm <- rbind(matchdfm,nomatchdfm)
    cat("done \n")
    return(dfm)
}
