######################################################################################################
#' @title Recur
#' @aliases Recur
#' @export
#' @description
#' Prepares recurrent event data for analysis by structuring key variables,
#' including event times, censoring indicators, and recurrence status.
#'
#' @param time Numeric. Time until the occurrence of an event.
#' @param event Integer (0 or 1). Indicator of event occurrence (`1`) or censoring (`0`).
#' @param id Integer or Character. Unique identifier for each individual.
#' @param IndRec Integer (0 or 1). Indicator of whether the individual has at least one recurrence (`1`) or not (`0`).
#'
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
