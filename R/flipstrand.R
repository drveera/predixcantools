#' flipstrand
#'
#' flipstrand
#'
#' @param x a vector of alleles
#'
#' @export
flipstrand <- function(x){
    if(x == 'A'){return('T')}
    if(x == 'T'){return('A')}
    if(x == 'G'){return('C')}
    if(x == 'C'){return('G')}
}
