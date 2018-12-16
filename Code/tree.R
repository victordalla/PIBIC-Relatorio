library(data.tree)

start <- Node$new("Qual exame fazer primeiro?")
  exoma <- start$AddChild("Exoma - Diagnóstico esclarido?")
    exoma.nao <- exoma$AddChild("Não - Microarray - Diagnóstico esclarido?", prob = 0.77)
      exoma.nao.microarray.nao <- exoma.nao$AddChild("Não", prob = 0.9, custo = 11554)
      exoma.nao.microarray.sim <- exoma.nao$AddChild("Sim", prob = 0.1, custo = 11554)
    exoma.sim <- exoma$AddChild("Sim", prob = 0.23, custo = 7500)
    
  microarray <- start$AddChild("Microarray - Diagnóstico esclarido?")
  microarray.nao <- microarray$AddChild("Não - Exoma - Diagnóstico esclarido?", prob = 0.9)
  microarray.nao.exoma.nao <- microarray.nao$AddChild("Não", prob = 0.77, custo = 11554)
  microarray.nao.exoma.sim <- microarray.nao$AddChild("Sim", prob = 0.23, custo = 11554)
  microarray.sim <- microarray$AddChild("Sim", prob = 0.1, custo = 4054)
  

print(start, "prob", "custo")
      
plot(start)
