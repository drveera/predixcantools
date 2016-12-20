#' allele.check
#'
#' allele.check
#'
#' @param summary summary data frame
#' @param db a single or a vector of dbfile names(with full absolute path) as character vector
#' @param snp.col snp column in summary
#' @param effect.allele.col effect allele column in summary
#' @param alt.allele.col alternative allele column in summary
#' @param beta.col beta column in summary; if you have odds ration, change it beta by log(odds ratio)
#' 
#' @export
allele.check <- function(summary,
                         db,
                         snp.col="rsid",
                         effect.allele.col="a1",
                         alt.allele.col="a2",
                         beta.col = "beta"){
    dblist <- list()
    for(i in 1:length(db)){
        ##create a connection
        dbcon <- dbConnect(dbDriver("SQLite"), dbname=db[i])
        ##read the weights
        weights <- dbReadTable(dbcon,"weights")
        dblist[[i]] <- data.table(weights)
        dbDisconnect(dbcon)
    }
    dbweights <- do.call(rbind,dblist)
    ##remove duplicates
    dbweights <- dbweights[!duplicated(rsid)]

    ##subset summary file
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
    summary <- summary[!duplicated(rsid)]
    
    dfm <- merge(summary,dbweights, by= "rsid")
    dfm <- data.table(dfm)

    ## seperate variants with strand flips
    dfm.flip <- dfm[!(ref_allelle == a2 | ref_allele == a1)]
    dfm.noflip <- dfm[(ref_allelle == a2 | ref_allele == a1)]
    ##flip the strand
    dfm.flip$a1 <- flipstrand(dfm.flip$a1)
    dfm.flip$a2 <- flipstrand(dfm.flip$a2)
    ##take only those flipped and drop others
    dfm.flip <- dfm.flip[(ref_allele == a2 | ref_allele == a1) & (eff_allele == a2 | eff_allele == a1)]
    ##merge back
    dfm <- rbind(dfm.flip,dfm.noflip)
    ##split the dfm with respect to allele match
    dfm$amatch <- with(dfm, a1 == eff_allele & a2 == ref_allele)
    dfm$amatch <- ifelse(dfm$amatch,"match","nomatch")
    ##split
    dfm.split <- split(dfm,dfm$amatch)
    match <- dfm.split$match
    nomatch <- dfm.split$nomatch
    nomatch$beta <- nomatch$beta * -1
    a1 <- nomatch$a2
    a2 <- nomatch$a1
    nomatch$a1 <- a1
    nomatch$a2 <- a2
    ##combine back
    dfm <- rbind(match,nomatch)
    return(dfm)
}
