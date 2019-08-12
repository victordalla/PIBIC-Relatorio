library(magrittr)
library(tibble)

BetaParams <- function(mu, sig) {
  sig <- ifelse(sig > 0.5^2, 0.5^2, sig)
  alpha <- ((1 - mu)/sig - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

ComputeTree <- function(Exame, 
                        REx, RArr, RCa, 
                        CEx, CArr, CCa) {
    
    if (Exame == 'Cariótipo com Array') {
      return(data.frame(
        Exame = Exame,
        Custo_esper  = CCa + (1-RCa)*CArr, 
        Efetiv_esper = 100*RCa + 100*(1-RCa)*RArr
        ))
    }
    # Cariótipo com Array e Exoma
    else if (Exame == 'Cariótipo') {
      return(data.frame(
        Exame = Exame,
        Custo_esper  = CCa + (1-RCa)*CArr + (1-RCa)*(1-RArr)*CEx, 
        Efetiv_esper = 100*RCa + 100*(1-RCa)*RArr + 100*(1-RCa)*(1-RArr)*REx
        ))
    }

    else if (Exame == 'Microarray') {
      return(data.frame(
        Exame = Exame,
        Custo_esper  = CArr + (1-RArr)*CEx, 
        Efetiv_esper = 100*RArr + 100*(1-RArr)*REx
        ))
    }
    
    else if (Exame == 'Exoma') {
      return(data.frame(
        Exame = Exame,
        Custo_esper  = CEx + (1-REx)*CArr, 
        Efetiv_esper = 100*REx + 100*(1-REx)*RArr
        ))
    }
  }

ComputeRCEI <- function(df, Exame, base) {
    df$Custo_incr[df$Exame == Exame]  <- df$Custo_esper[df$Exame == Exame] - df$Custo_esper[df$Exame == base]
    df$Efetiv_incr[df$Exame == Exame] <- df$Efetiv_esper[df$Exame == Exame] - df$Efetiv_esper[df$Exame == base]
    df$RCEI[df$Exame == Exame]        <- df$Custo_incr[df$Exame == Exame] / df$Efetiv_incr[df$Exame == Exame]
    return(df)
  }

param <- list(
  Exoma = list(
    Rend = list(base = 0.32), 
    Custo  = list(base = 6950)
    ), 
  Microarray = list(
    Rend = list(base = 0.14), 
    Custo  = list(base = 3490)
    ),
  Cariotipo = list(
    Rend = list(base = 0.03), 
    Custo  = list(base = 1377)
    )
  )
nrun <- 1000

run_simulation <- function() {

## Exoma
REx <- c(27, 32, 41, 32.5, 27.2) / 100
REx_simul <- BetaParams(mean(REx), var(REx)) %$%
  rbeta(nrun, alpha, beta)
CEx <- 6000
CV_Ex <- 0.02
CEx_simul <- rgamma(nrun, mean(CEx) * CV_Ex, rate = CV_Ex)

## Array
RArr <- c(
  10, 14, 10, 9.8, 13.7, 18.6, 5.1, 11, 16.8, 13.6, 20, 16.7, 16, 55, 
  5.6, 35, 5.6, 10, 9.6, 7, 13.8, 6.9, 7.6, 6.3, 16.4, 5.3, 15.6, 7.8, 
  6.4, 18, 11.8, 17.1
  ) / 100
RArr_simul <- BetaParams(mean(RArr), var(RArr)) %$%
  rbeta(nrun, alpha, beta)
CArr <- 3400
CV_Arr <- 0.1
CArr_simul <- rgamma(nrun, mean(CArr) * CV_Arr, rate = CV_Arr)

## Cariotipo
RCa_simul <- BetaParams(0.03, 0.00001) %$%
  rbeta(nrun, alpha, beta)
CCa <- 1400
CV_Ca <- 0.15
CCa_simul <- rgamma(nrun, mean(CCa) * CV_Ca, rate = CV_Ca)

## bind data
anal_sens <- rbind(
  ComputeTree(
    'Cariótipo', 
    REx_simul, RArr_simul, RCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Microarray', 
    REx_simul, RArr_simul, RCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Exoma', 
    REx_simul, RArr_simul, RCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Cariótipo com Array', 
    REx_simul, RArr_simul, RCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
  )
  )
anal_sens$CER <- anal_sens$Custo_esper / anal_sens$Efetiv_esper


## RCEI
anal_sens %<>% 
  ComputeRCEI('Microarray', 'Cariótipo') %>%  
  ComputeRCEI('Exoma', 'Cariótipo')


## tornado data
base_cariotipo <- ComputeTree(
  'Cariótipo', 
  param$Exoma$Rend$base, param$Microarray$Rend$base, param$Cariotipo$Rend$base, 
  param$Exoma$Custo$base, param$Microarray$Custo$base, param$Cariotipo$Custo$base
  )

list(
  anal_sens      = anal_sens, 
  anal_sens_det  = anal_sens_det, 
  base_cariotipo = base_cariotipo, 
  REx_simul  = REx_simul, CEx_simul   = CEx_simul, 
  RArr_simul = RArr_simul, CArr_simul = CArr_simul, 
  RCa_simul  = RCa_simul, CCa_simul   = CCa_simul
)
}