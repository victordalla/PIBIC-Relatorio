library(ggplot2)
library(dplyr)

## ICER
ggplot(filter(anal_sens, Exame != "Cariótipo"), aes(Efetiv_incr, Custo_incr, col = Exame)) + 
  geom_point(alpha = 0.6) + 
  scale_y_continuous(labels = scales::dollar_format("R$")) +
  labs(x = "Efetividade incremental", y = "Custo incremental") + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  theme_classic() # + stat_ellipse()

ggplot(anal_sens, aes(Efetiv_esper, Custo_esper, col = Exame)) + 
  geom_point(alpha = 0.4) + 
  scale_y_continuous(labels = scales::dollar_format("R$")) +
  labs(x = "Efetividade esperada", y = "Custo esperado") + 
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  theme_classic()

## tornado
tor_exoma <- data_frame(
  fator = 'Efetividade do Exoma', 
  min = min(EfEx_simul), max = max(EfEx_simul), 
  
  ICER_base = 
    ComputeTree('Exoma', 
                param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
                param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base) %>% 
    rbind(base_cariotipo) %>% 
    ComputeICER('Exoma', 'Cariótipo') %>% 
    `[`(1, 'ICER'),
  
  ICER_min = 
    ComputeTree('Exoma', 
                min(EfEx_simul), param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
                param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base) %>% 
    rbind(base_cariotipo) %>% 
    ComputeICER('Exoma', 'Cariótipo') %>% 
    `[`(1, 'ICER'), 
  
  ICER_max = 
    ComputeTree('Exoma', 
                max(EfEx_simul), param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
                param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base) %>% 
    rbind(base_cariotipo) %>% 
    ComputeICER('Exoma', 'Cariótipo') %>% 
    `[`(1, 'ICER')
) %>% 
  rbind(data_frame(
    fator = 'Custo do Exoma', 
    min = min(CEx_simul), max = max(CEx_simul), 
    
    ICER_base = 
      ComputeTree('Exoma', 
                  param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
                  param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base) %>% 
      rbind(base_cariotipo) %>% 
      ComputeICER('Exoma', 'Cariótipo') %>% 
      `[`(1, 'ICER'),
    
    ICER_min = 
      ComputeTree('Exoma', 
                  param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
                  min(CEx_simul), param$Microarray$custo$base, param$Cariotipo$custo$base) %>% 
      rbind(base_cariotipo) %>% 
      ComputeICER('Exoma', 'Cariótipo') %>% 
      `[`(1, 'ICER'), 
    
    ICER_max = 
      ComputeTree('Exoma', 
                  param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
                  max(CEx_simul), param$Microarray$custo$base, param$Cariotipo$custo$base) %>% 
      rbind(base_cariotipo) %>% 
      ComputeICER('Exoma', 'Cariótipo') %>% 
      `[`(1, 'ICER')
  )
  )
ggplot(tor_exoma) +
  geom_segment(aes(x=fator, xend=fator, y=ICER_min, yend=ICER_max), color="grey", size = 4) +
  geom_point(aes(x=fator, y=ICER_min), color=rgb(0.2,0.7,0.1,0.5), size=3) +
  geom_point(aes(x=fator, y=ICER_max), color=rgb(0.7,0.2,0.1,0.5), size=3) +
  coord_flip() +
  geom_hline(yintercept = tor_exoma$ICER_base[1]) +
  theme_classic() +
  theme(legend.position = "none", panel.border = element_blank()) +
  xlab("") +
  ylab(paste('ICER da análise de sensitividade.', 
             'A Linha vertical é o ICER encontrado na árvore de decisão', sep = '\n'))

tor_exoma <- data_frame(
  fator = rep(tor_exoma$fator, 2), grupo = c('min', 'min', 'max', 'max'),
  valor = c(tor_exoma$min, tor_exoma$max),
  ICER_dev = c(tor_exoma$ICER_min, tor_exoma$ICER_max) - tor_exoma$ICER_base, 
  ICER_base = tor_exoma$ICER_base[1], 
  ajuste = c(30, -30, -30, 30)
)
ggplot(tor_exoma, aes(x = as.factor(fator), y = ICER_dev, fill = grupo)) + 
  geom_bar(stat = "identity", position = "identity", show.legend = FALSE) + 
  geom_text(aes(y = round(ICER_dev) + ajuste, label = round(valor, 2)), size = 3, vjust = 0.5) +
  xlab("") +
  coord_flip() + 
  ylab(paste('Desvios do ICER da análise de sensitividade', 
             'com relação ao ICER encontrado na árvore de decisão', sep = '\n'))
