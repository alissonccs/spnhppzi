######################################################################################################
#' @title func_sim
#' @aliases func_sim
#' @export
#' @description func_sim: Executa simulações
#' @param time: tempo até ocorrência
#' @param event: identificador de ocorrência de evento ou censura
#' @param id: identificador do indivíduo
#' @param IndRec: Indicador se o indivíduo apresenta uma ou mais recorrência
#' @param N: N
#' @param dist.x: Distribuição das covariáveis
#' @param par.x: Parâmetros das distribuições associadas às covariáveis
#' @param beta.x: Coeficientes associados às covariáveis
#' @param dist.z: Distribuição dos efeitos aleatórios
#' @param par.z: Parâmetros da distribuição dos efeitos aleatórios
#' @param dist.rec: Distribuição dos eventos recorrentes
#' @param par.rec: Parâmetros da distribuição dos eventos recorrentes
#' @param fun.min: Período de acompanhamento mínimo
#' @param fun.max: Período de acompanhamento máximo
#' @param cens.prob: Percentual de indivíduos censurados antes do fim do período de acompanhamento
#' @param dfree: Tempo fora de risco após experimentar um primeiro evento
#' @param pfree: Proporção de indivíduos fora de risco (sem utilizar covariáveis)
#' @param pfree1: Probabilidade de um indivíduo não recorrer (regressão logística)
#' @param logist: Uso da regressão logística
#' @param ent.dis.z: indicadora - entrada dos parâmetros da normal (1) ou da lognormal(0) para os efeitos aleatórios
#' @param mu.z: média da normal para os efeitos aleatórios (caso ent.dis.z=1)
#' @param sigma.z desvio da normal para os efeitos aleatórios (caso ent.dis.z=1)
#' @export
######################################################################################################

func_sim<-function(N,
                   dist.x,
                   par.x,
                   beta.x,
                   dist.z,
                   par.z,
                   dist.rec,
                   par.rec,
                   fun.min,
                   fun.max,
                   cens.prob,
                   dfree,
                   pfree,
                   pfree1,
                   logist,
                   ent.dis.z,
                   mu.z,
                   sigma.z){
  simdata <- simrecev(N, fun.min, fun.max, cens.prob, dist.x, par.x, beta.x,
                     dist.z, par.z, dist.rec, par.rec, pfree,pfree1, dfree,logist,ent.dis.z,
                     mu.z,sigma.z)
  print(prop.table(table(as.vector(table(simdata$id)))))
  simdata <- simdata %>%
    group_by(id)%>%
    mutate(ngrupo=n(),
           IndRec=case_when(ngrupo>2~1, TRUE~0),
           max_stop=max(stop)
           # stop1=case_when(max_stop==stop~5,TRUE~stop),
           # max_stop1=max(stop1)
    )

  # View(simdata)
  # saveRDS(simdata,file="/home/alisson/Documents/R/NHPP/Tabelas/simdata1.rds")
  return(simdata)
}
