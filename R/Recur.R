######################################################################################################
#' @title Recur
#' @aliases Recur
#' @export
#' @description Recur: Preparação dos dados
#' @param time: tempo até ocorrência
#' @param event: identificador de ocorrência de evento ou censura
#' @param id: identificador do indivíduo
#' @param IndRec: Indicador se o indivíduo apresenta uma ou mais recorrência
#' @export
######################################################################################################
Recur <- function(time, event = NULL, id = NULL, IndRec=NULL){
  n <- length(time)
  if(is.null(event)){
    event <- rep(1, n)
  }
  if(is.null(id)){
    id <- rep(1, n)
  }
  if(is.null(event)){
    IndRec <- rep(2, n)
  }
  resp <- cbind(time, event, id, IndRec)
  colnames(resp) <- c("time", "event", "id","IndRec")
  return(resp)
}
