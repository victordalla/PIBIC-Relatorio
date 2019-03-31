library(magrittr)
library(tibble)

BetaParams <- function(mu, sig) {
  sig <- ifelse(sig > 0.5^2, 0.5^2, sig)
  alpha <- ((1 - mu)/sig - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

ComputeTree <- function(Exame, 
                        EfEx_simul, EfArr_simul, EfCa_simul, 
                        CEx_simul, CArr_simul, CCa_simul) {
    
    if (Exame == 'Cariótipo com Array') {
      return(data.frame(
        Exame = 'Cariótipo com Array',
        Custo_esper  = EfCa_simul*CCa_simul + (CCa_simul+CArr_simul)*(1-EfCa_simul)*EfArr_simul, 
        Efetiv_esper = 100*EfCa_simul + 100*(1-EfCa_simul)*EfArr_simul
        ))
    }

    else if (Exame == 'Cariótipo com Exoma') {
      return(data.frame(
        Exame = 'Cariótipo com Exoma',
        Custo_esper  = EfCa_simul*CCa_simul + (CCa_simul+CArr_simul)*(1-EfCa_simul)*EfArr_simul + (CCa_simul+CArr_simul+CEx_simul)*(1-EfCa_simul)*(1-EfArr_simul), 
        Efetiv_esper = 100*EfCa_simul + 100*(1-EfCa_simul)*EfArr_simul + 100*(1-EfCa_simul)*(1-EfArr_simul)*EfEx_simul
        ))
    }

    else if (Exame == 'Microarray') {
      return(data.frame(
        Exame = 'Microarray',
        Custo_esper  = EfArr_simul*CArr_simul + (CArr_simul+CEx_simul)*(1-EfArr_simul), 
        Efetiv_esper = 100*EfArr_simul + 100*(1-EfArr_simul)*EfEx_simul
        ))
    }
    
    else if (Exame == 'Exoma') {
      return(data.frame(
        Exame = 'Exoma',
        Custo_esper  = EfEx_simul*CEx_simul + (CEx_simul+CArr_simul)*(1-EfEx_simul), 
        Efetiv_esper = 100*EfEx_simul + 100*(1-EfEx_simul)*EfArr_simul
        ))
    }
  }

ComputeICER <- function(df, Exame, base) {
    df$Custo_incr[df$Exame == Exame]  <- df$Custo_esper[df$Exame == Exame] - df$Custo_esper[df$Exame == base]
    df$Efetiv_incr[df$Exame == Exame] <- df$Efetiv_esper[df$Exame == Exame] - df$Efetiv_esper[df$Exame == base]
    df$ICER[df$Exame == Exame]        <- df$Custo_incr[df$Exame == Exame] / df$Efetiv_incr[df$Exame == Exame]
    return(df)
  }

  
param <- list(
  Exoma = list(
    efetiv = list(base = 0.27), 
    custo  = list(base = 7500)
    ), 
  Microarray = list(
    efetiv = list(base = 0.13), 
    custo  = list(base = 4054)
    ),
  Cariotipo = list(
    efetiv = list(base = 0.03), 
    custo  = list(base = 1390)
    )
  )
nrun <- 1000

run_simulation <- function() {
  
## Exoma
EfEx <- c(27, 32, 41, 32.5, 27.2) / 100
set.seed(42)
EfEx_simul <- BetaParams(mean(EfEx), var(EfEx)) %$%
  rbeta(nrun, alpha, beta)
CEx <- c(5500, 7500, 9900)
CV_Ex <- 50
set.seed(42)
CEx_simul <- rgamma(nrun, mean(CEx) * CV_Ex, rate = CV_Ex)

## Array
EfArr <- c(
  10, 14, 10, 9.8, 13.7, 18.6, 5.1, 11, 16.8, 13.6, 20, 16.7, 16, 55, 
  5.6, 35, 5.6, 10, 9.6, 7, 13.8, 6.9, 7.6, 6.3, 16.4, 5.3, 15.6, 7.8, 
  6.4, 18, 11.8, 17.1
  ) / 100
set.seed(42)
EfArr_simul <- BetaParams(mean(EfArr), var(EfArr)) %$%
  rbeta(nrun, alpha, beta)
CArr <- c(4054, 3000)
CV_Arr <- 30
set.seed(42)
CArr_simul <- rgamma(nrun, mean(CArr) * CV_Arr, rate = CV_Arr)

## Cariotipo
EfCa <- c(0.01, 0.02, 0.03)
set.seed(42)
EfCa_simul <- BetaParams(mean(EfCa), var(EfCa)) %$%
  rbeta(nrun, alpha, beta)
CCa <- c(1390, 1500)
CV_Ca <- 20
set.seed(42)
CCa_simul <- rgamma(nrun, mean(CCa) * CV_Ca, rate = CV_Ca)

## bind data
anal_sens <- rbind(
  ComputeTree(
    'Cariótipo com Array', 
    EfEx_simul, EfArr_simul, EfCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Cariótipo com Exoma', 
    EfEx_simul, EfArr_simul, EfCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Microarray', 
    EfEx_simul, EfArr_simul, EfCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    ), 
  ComputeTree(
    'Exoma', 
    EfEx_simul, EfArr_simul, EfCa_simul, 
    CEx_simul, CArr_simul, CCa_simul
    )
  )
anal_sens$CER <- anal_sens$Custo_esper / anal_sens$Efetiv_esper


## ICER
anal_sens %<>% ComputeICER('Cariótipo com Exoma', 'Cariótipo com Array') %>% 
  ComputeICER('Microarray', 'Cariótipo com Array') %>%  
  ComputeICER('Exoma', 'Cariótipo com Array')

## tornado data
base_cariotipo <- ComputeTree(
  'Cariótipo com Array', 
  param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
  param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
  )

list(
  anal_sens      = anal_sens, 
  base_cariotipo = base_cariotipo, 
  EfEx_simul  = EfEx_simul, CEx_simul   = CEx_simul, 
  EfArr_simul = EfArr_simul, CArr_simul = CArr_simul, 
  EfCa_simul  = EfCa_simul, CCa_simul   = CCa_simul
)
}