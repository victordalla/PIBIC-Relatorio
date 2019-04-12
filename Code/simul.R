library(magrittr)
library(tibble)

BetaParams <- function(mu, sig) {
  sig <- ifelse(sig > 0.5^2, 0.5^2, sig)
  alpha <- ((1 - mu)/sig - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

ComputeTree <- function(Exame, 
                        REx_simul, RArr_simul, RCa_simul, 
                        CEx_simul, CArr_simul, CCa_simul) {
    
    if (Exame == 'Cariótipo com Array') {
      return(data.frame(
        Exame = 'Cariótipo com Array',
        Custo_esper  = RCa_simul*CCa_simul + (CCa_simul+CArr_simul)*(1-RCa_simul)*RArr_simul, 
        Efetiv_esper = 100*RCa_simul + 100*(1-RCa_simul)*RArr_simul
        ))
    }

    else if (Exame == 'Cariótipo com Exoma') {
      return(data.frame(
        Exame = 'Cariótipo com Exoma',
        Custo_esper  = RCa_simul*CCa_simul + (CCa_simul+CArr_simul)*(1-RCa_simul)*RArr_simul + (CCa_simul+CArr_simul+CEx_simul)*(1-RCa_simul)*(1-RArr_simul), 
        Efetiv_esper = 100*RCa_simul + 100*(1-RCa_simul)*RArr_simul + 100*(1-RCa_simul)*(1-RArr_simul)*REx_simul
        ))
    }

    else if (Exame == 'Microarray') {
      return(data.frame(
        Exame = 'Microarray',
        Custo_esper  = RArr_simul*CArr_simul + (CArr_simul+CEx_simul)*(1-RArr_simul), 
        Efetiv_esper = 100*RArr_simul + 100*(1-RArr_simul)*REx_simul
        ))
    }
    
    else if (Exame == 'Exoma') {
      return(data.frame(
        Exame = 'Exoma',
        Custo_esper  = REx_simul*CEx_simul + (CEx_simul+CArr_simul)*(1-REx_simul), 
        Efetiv_esper = 100*REx_simul + 100*(1-REx_simul)*RArr_simul
        ))
    }
  }

ComputeRCEI <- function(df, Exame, base) {
    df$Custo_incr[df$Exame == Exame]  <- df$Custo_esper[df$Exame == Exame] - df$Custo_esper[df$Exame == base]
    df$Efetiv_incr[df$Exame == Exame]   <- df$Efetiv_esper[df$Exame == Exame] - df$Efetiv_esper[df$Exame == base]
    df$RCEI[df$Exame == Exame]        <- df$Custo_incr[df$Exame == Exame] / df$Efetiv_incr[df$Exame == Exame]
    return(df)
  }

  
param <- list(
  Exoma = list(
    Rend = list(base = 0.27), 
    custo  = list(base = 7500)
    ), 
  Microarray = list(
    Rend = list(base = 0.13), 
    custo  = list(base = 4054)
    ),
  Cariotipo = list(
    Rend = list(base = 0.03), 
    custo  = list(base = 1390)
    )
  )
nrun <- 1000

run_simulation <- function() {
  
## Exoma
REx <- c(27, 32, 41, 32.5, 27.2) / 100
set.seed(42)
REx_simul <- BetaParams(mean(REx), var(REx)) %$%
  rbeta(nrun, alpha, beta)
CEx <- c(5500, 7500, 9900)
CV_Ex <- 50
set.seed(42)
CEx_simul <- rgamma(nrun, mean(CEx) * CV_Ex, rate = CV_Ex)

## Array
RArr <- c(
  10, 14, 10, 9.8, 13.7, 18.6, 5.1, 11, 16.8, 13.6, 20, 16.7, 16, 55, 
  5.6, 35, 5.6, 10, 9.6, 7, 13.8, 6.9, 7.6, 6.3, 16.4, 5.3, 15.6, 7.8, 
  6.4, 18, 11.8, 17.1
  ) / 100
set.seed(42)
RArr_simul <- BetaParams(mean(RArr), var(RArr)) %$%
  rbeta(nrun, alpha, beta)
CArr <- c(4054, 3000)
CV_Arr <- 30
set.seed(42)
CArr_simul <- rgamma(nrun, mean(CArr) * CV_Arr, rate = CV_Arr)

## Cariotipo
RCa <- c(0.01, 0.02, 0.03)
set.seed(42)
RCa_simul <- BetaParams(mean(RCa), var(RCa)) %$%
  rbeta(nrun, alpha, beta)
CCa <- c(1390, 1500)
CV_Ca <- 20
set.seed(42)
CCa_simul <- rgamma(nrun, mean(CCa) * CV_Ca, rate = CV_Ca)

## bind data
anal_sens <- rbind(
  ComputeTree(
    'Cariótipo com Array', 
    REx_simul, RArr_simul, RCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Cariótipo com Exoma', 
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
    )
  )
anal_sens$CER <- anal_sens$Custo_esper / anal_sens$Efetiv_esper


## RCEI
anal_sens %<>% ComputeRCEI('Cariótipo com Exoma', 'Cariótipo com Array') %>% 
  ComputeRCEI('Microarray', 'Cariótipo com Array') %>%  
  ComputeRCEI('Exoma', 'Cariótipo com Array')

## tornado data
base_cariotipo <- ComputeTree(
  'Cariótipo com Array', 
  param$Exoma$Rend$base, param$Microarray$Rend$base, param$Cariotipo$Rend$base, 
  param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
  )

list(
  anal_sens      = anal_sens, 
  base_cariotipo = base_cariotipo, 
  REx_simul  = REx_simul, CEx_simul   = CEx_simul, 
  RArr_simul = RArr_simul, CArr_simul = CArr_simul, 
  RCa_simul  = RCa_simul, CCa_simul   = CCa_simul
)
}