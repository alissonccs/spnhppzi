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
Recur1 <- function(time, event = NULL, id = NULL,SP_ID, IndRec=NULL){
  n <- length(time)
  if(is.null(event)){
    event <- rep(1, n)
  }
  if(is.null(id)){
    id <- rep(1, n)
  }
  if(is.null(SP_ID)){
    SP_ID <- rep(1, n)
  }
  if(is.null(event)){
    IndRec <- rep(2, n)
  }
  resp <- cbind(time, event, id,SP_ID, IndRec)
  colnames(resp) <- c("time", "event", "id","SP_ID","IndRec")
  return(resp)
}
