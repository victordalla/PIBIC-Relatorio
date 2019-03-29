library(ggplot2)
library(dplyr)

ggplot(tibble(x = CEx_simul), aes(x)) + 
  ggdistribute::stat_density_ci(ci_width = 0.95) + 
  labs(x = "", y = "densidade estimada ")

## ICER
filter(anal_sens, Exame != "Cariótipo com Array") %>% 
  ggplot(aes(Efetiv_incr, Custo_incr, col = Exame, fill = Exame)) + 
  stat_ellipse(geom = "polygon", alpha = 0.3) + 
  geom_point(alpha = 0.7) + 
  scale_y_continuous(labels = scales::dollar_format("R$")) + 
  scale_x_continuous(labels = function(x) paste0(x, "%")) + 
  labs(x = "Efetividade incremental", y = "Custo incremental") + 
  theme_classic() 
  # geom_hline(yintercept = 0) + 
  # geom_vline(xintercept = 0) + 
  # theme(
  #   axis.line    = element_blank(), 
  #   axis.text.x  = element_blank(),
  #   axis.ticks.x = element_blank(), 
  #   axis.title.y = element_blank(),
  #   axis.text.y  = element_blank(),
  #   axis.ticks.y = element_blank()
  # )

ggplot(anal_sens, aes(Efetiv_esper, Custo_esper, col = Exame, fill = Exame)) + 
  stat_ellipse(geom = "polygon", alpha = 0.1) + 
  geom_point(alpha = 0.3) + 
  scale_y_continuous(labels = scales::dollar_format("R$")) +
  stat_ellipse() + 
  labs(x = "Efetividade esperada", y = "Custo esperado") + 
  # geom_hline(yintercept = 0, linetype = 'dashed') + 
  # geom_vline(xintercept = 0, linetype = 'dashed') + 
  theme_classic()

## tornado
tor_exoma <- tibble(
  fator = 'Efetividade do Exoma', 
  min = min(EfEx_simul), max = max(EfEx_simul), 
  
  ICER_base = ComputeTree(
    'Exoma', 
    param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
    param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
    ) %>% 
    rbind(base_cariotipo) %>% 
    ComputeICER('Exoma', 'Cariótipo') %$% 
    ICER[1],
  
  ICER_min = ComputeTree(
    'Exoma', 
    min(EfEx_simul), param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
    param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
    ) %>% 
    rbind(base_cariotipo) %>% 
    ComputeICER('Exoma', 'Cariótipo') %$% 
    ICER[1], 
  
  ICER_max = ComputeTree(
    'Exoma', 
    max(EfEx_simul), param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
    param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
    ) %>% 
    rbind(base_cariotipo) %>% 
    ComputeICER('Exoma', 'Cariótipo') %$% 
    ICER[1],
) %>% 
  rbind(data_frame(
    fator = 'Custo do Exoma', 
    min = min(CEx_simul), max = max(CEx_simul), 
    
    ICER_base = ComputeTree(
      'Exoma', 
      param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
      param$Exoma$custo$base, param$Microarray$custo$base, param$Cariotipo$custo$base
      ) %>% 
      rbind(base_cariotipo) %>% 
      ComputeICER('Exoma', 'Cariótipo') %$% 
      ICER[1],
    
    ICER_min = ComputeTree(
      'Exoma', 
      param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
      min(CEx_simul), param$Microarray$custo$base, param$Cariotipo$custo$base
      ) %>% 
      rbind(base_cariotipo) %>% 
      ComputeICER('Exoma', 'Cariótipo') %$% 
      ICER[1],
    
    ICER_max = ComputeTree(
      'Exoma', 
      param$Exoma$efetiv$base, param$Microarray$efetiv$base, param$Cariotipo$efetiv$base, 
      max(CEx_simul), param$Microarray$custo$base, param$Cariotipo$custo$base
      ) %>% 
      rbind(base_cariotipo) %>% 
      ComputeICER('Exoma', 'Cariótipo') %$% 
      ICER[1]
  )
  )
ggplot(tor_exoma) +
  geom_segment(aes(x=fator, xend=fator, y=ICER_min, yend=ICER_max), color="grey", size = 4) +
  geom_point(aes(x=fator, y=ICER_min), color=rgb(0.2,0.7,0.1,0.5), size=3) +
  geom_point(aes(x=fator, y=ICER_max), color=rgb(0.7,0.2,0.1,0.5), size=3) +
  geom_hline(yintercept = tor_exoma$ICER_base[1]) +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "none", panel.border = element_blank()) +
  labs(
    x = '', 
    y = 'ICER da análise de sensitividade', 
    subtitle = 'A Linha vertical é o ICER encontrado na árvore de decisão'
    )

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
  coord_flip() + 
  labs(
    x = '', 
    y = 'Desvios do ICER da análise de sensitividade\ncom relação ao ICER encontrado na árvore de decisão'
  )