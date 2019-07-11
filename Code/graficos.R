library(ggplot2)
library(dplyr)

## RCEI
filter(anal_sens, Exame != "Cariótipo") %>% 
  ggplot(aes(Efetiv_incr, Custo_incr, col = Exame, fill = Exame)) + 
  stat_ellipse(geom = "polygon", alpha = 0.3) + 
  geom_point(alpha = 0.7) + 
  scale_y_continuous(labels = scales::dollar_format("R$")) + 
  scale_x_continuous(labels = function(x) paste0(x, "%")) + 
  labs(x = "Efetividade incremental", y = "Custo incremental") + 
  theme_classic() 


ggplot(anal_sens, aes(Efetiv_esper, Custo_esper, col = Exame, fill = Exame)) + 
  stat_ellipse(geom = "polygon", alpha = 0.1) + 
  geom_point(alpha = 0.3) + 
  scale_y_continuous(labels = scales::dollar_format("R$")) +
  stat_ellipse() + 
  labs(x = "Efetividade esperada", y = "Custo esperado") + 
  # geom_hline(yintercept = 0, linetype = "dashed") + 
  # geom_vline(xintercept = 0, linetype = "dashed") + 
  theme_classic()

ggplot(tor_exoma, aes(x = fator, y = RCEI_dev, fill = grupo)) + 
  geom_bar(stat = "identity", position = "identity", show.legend = FALSE) + 
  geom_text(
    aes(y = round(RCEI_dev) + ajuste, label = round(valor, 2)), 
    size = 3, vjust = 0.5
  ) +
  coord_flip() + 
  labs(
    title = "Desvios do RCEI da análise de sensitividade\ncom relação ao RCEI encontrado na árvore de decisão", 
    x = ""
  )
