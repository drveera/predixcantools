#' extract_weights
#'
#' extract_weights
#'
#' @param db a chracter vector of db filenames with absolute path to read
#'
#' @export
extract_weights <- function(db){
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
    dbweights <- dbweights[!duplicated(dbweights$rsid),]
    return(dbweights)
}
