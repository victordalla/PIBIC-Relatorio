BetaParams <- function(mu, sig) {
  sig <- ifelse(sig > 0.5^2, 0.5^2, sig)
  alpha <- ((1 - mu)/sig - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

ComputeTree <-
  function(Exame, 
           EfEx_simul, EfArr_simul, EfCa_simul, 
           CEx_simul, CArr_simul, CCa_simul) {
    
    if (Exame == 'Cariótipo') {
      return(data.frame(Exame = 'Cariótipo',
        Custo_esper = EfCa_simul*CCa_simul + (CCa_simul+CArr_simul)*(1-EfCa_simul)*EfArr_simul, 
        Efetiv_esper = 100*EfCa_simul + 100*(1-EfCa_simul)*EfArr_simul))
    }
    
    else if (Exame == 'Cariótipo com Exoma') {
      return(data.frame(Exame = 'Cariótipo com Exoma',
        Custo_esper = EfCa_simul*CCa_simul + (CCa_simul+CArr_simul)*(1-EfCa_simul)*EfArr_simul + (CCa_simul+CArr_simul+CEx_simul)*(1-EfCa_simul)*(1-EfArr_simul), 
        Efetiv_esper = 100*EfCa_simul + 100*(1-EfCa_simul)*EfArr_simul + 100*(1-EfCa_simul)*(1-EfArr_simul)*EfEx_simul))
    }

    else if (Exame == 'Microarray') {
      return(data.frame(Exame = 'Microarray',
        Custo_esper = EfArr_simul*CArr_simul + (CArr_simul+CEx_simul)*(1-EfArr_simul), 
        Efetiv_esper = 100*EfArr_simul + 100*(1-EfArr_simul)*EfEx_simul))
    }
    
    else if (Exame == 'Exoma') {
      return(data.frame(Exame = 'Exoma',
        Custo_esper = EfEx_simul*CEx_simul + (CEx_simul+CArr_simul)*(1-EfEx_simul), 
        Efetiv_esper = 100*EfEx_simul + 100*(1-EfEx_simul)*EfArr_simul))
    }
  }

ComputeICER <- function(df, Exame, base) {
    df$Custo_incr[df$Exame == Exame] <- 
      df$Custo_esper[df$Exame == Exame] - df$Custo_esper[df$Exame == base]
    df$Efetiv_incr[df$Exame == Exame] <- 
      df$Efetiv_esper[df$Exame == Exame] - df$Efetiv_esper[df$Exame == base]
    df$ICER[df$Exame == Exame] <- 
      df$Custo_incr[df$Exame == Exame] / df$Efetiv_incr[df$Exame == Exame]
    
    return(df)
  }

param <- list(
  Exoma      = list(efetiv = list(base = 0.27, mean = 0, var = 0), 
                    custo = list(base = 7500, mean = 0, var = 0)), 
  Microarray = list(efetiv = list(base = 0.13, mean = 0, var = 0), 
                    custo = list(base = 4054, mean = 0, var = 0)),
  Cariotipo  = list(efetiv = list(base = 0.03, mean = 0, var = 0), 
                    custo = list(base = 1390, mean = 0, var = 0))
  )
nrun <- 1000


## Exoma
EfEx <- c(27, 32, 41, 32.5, 27.2) / 100
param$Exoma$efetiv$mean <- mean(EfEx)
param$Exoma$efetiv$var <- var(EfEx)
EfEx_simul <- 
  rbeta(
    nrun, 
    BetaParams(param$Exoma$efetiv$mean, param$Exoma$efetiv$var)$alpha, 
    BetaParams(param$Exoma$efetiv$mean, param$Exoma$efetiv$var)$beta
    )

CEx <- c(7500, 8000, 9500)
param$Exoma$custo$mean <- mean(log(CEx))
param$Exoma$custo$var <- var(log(CEx))
CEx_simul <- 
  rlnorm(
    nrun, 
    param$Exoma$custo$mean, 
    param$Exoma$custo$var
    )

## Array
EfArr <- c(10, 14, 10, 9.8, 13.7, 18.6, 5.1, 11, 16.8, 13.6, 20, 16.7, 16, 55, 
           5.6, 35, 5.6, 10, 9.6, 7, 13.8, 6.9, 7.6, 6.3, 16.4, 5.3, 15.6, 7.8, 
           6.4, 18, 11.8, 17.1) / 100
param$Microarray$efetiv$mean <- mean(EfArr)
param$Microarray$efetiv$var <- var(EfArr)
EfArr_simul <- 
  rbeta(
    nrun, 
    BetaParams(param$Microarray$efetiv$mean, param$Microarray$efetiv$var)$alpha, 
    BetaParams(param$Microarray$efetiv$mean, param$Microarray$efetiv$var)$beta
    )

CArr <- c(4054, 3000)
param$Microarray$custo$mean <- mean(log(CArr))
param$Microarray$custo$var <- var(log(CArr))
CArr_simul <- 
  rlnorm(
    nrun, 
    param$Microarray$custo$mean, 
    param$Microarray$custo$var
    )

## Cariotipo
EfCa <- c(0.025, 0.03, 0.035)
param$Cariotipo$efetiv$mean <- mean(EfCa)
param$Cariotipo$efetiv$var <- var(EfCa)
EfCa_simul <- 
  rbeta(
    nrun, 
    BetaParams(param$Cariotipo$efetiv$mean, param$Cariotipo$efetiv$var)$alpha, 
    BetaParams(param$Cariotipo$efetiv$mean, param$Cariotipo$efetiv$var)$beta
    )

CCa <- c(1390, 1500)
param$Cariotipo$custo$mean <- mean(log(CCa))
param$Cariotipo$custo$var <- var(log(CCa))
CCa_simul <- 
  rlnorm(
    nrun, 
    param$Cariotipo$custo$mean, 
    param$Cariotipo$custo$var
    )

## bind data
anal_sens <- rbind(
  ComputeTree('Cariótipo', 
              EfEx_simul, EfArr_simul, EfCa_simul, 
              CEx_simul, CArr_simul, CCa_simul), 
  ComputeTree('Cariótipo com Exoma', 
              EfEx_simul, EfArr_simul, EfCa_simul, 
              CEx_simul, CArr_simul, CCa_simul), 
  ComputeTree('Microarray', 
              EfEx_simul, EfArr_simul, EfCa_simul, 
              CEx_simul, CArr_simul, CCa_simul), 
  ComputeTree('Exoma', 
              EfEx_simul, EfArr_simul, EfCa_simul, 
              CEx_simul, CArr_simul, CCa_simul)
  )
anal_sens$CER <- anal_sens$Custo_esper / anal_sens$Efetiv_esper


## ICER
library(magrittr)
anal_sens %<>% ComputeICER('Cariótipo com Exoma', 'Cariótipo') %>% 
  ComputeICER('Microarray', 'Cariótipo') %>%  
  ComputeICER('Exoma', 'Cariótipo')


## tornado data
library(tibble)
base_cariotipo <- 
  ComputeTree('Cariótipo', 
              param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
              param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
              )
