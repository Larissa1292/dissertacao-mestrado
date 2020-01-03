# Simulação de MC para o modelo 3 (Somente com erro de medida)
# Modelo somente com erro de medida, logo não tem erro de classificacao. Com isso, pi0 = 0 e pi1 = 0.

#Gerar M3 e aplicar em M1

require(fExtremes)
require(mvtnorm)
require(sn)
library(optimParallel)

argumentos = commandArgs(trailingOnly = TRUE) #comando para inserir os valores dos parametros via terminal

cl <- makeCluster(3)     # definindo o num de clusters no processador
setDefaultCluster(cl=cl) # definindo 'cl' como cluster padrao

#### Definindo os parâmetros iniciais (via terminal) ####

pi0 <- as.numeric(argumentos[1])
pi1 <- as.numeric(argumentos[2])
beta0 <- as.numeric(argumentos[3])
beta1 <- as.numeric(argumentos[4])
lambda <- as.numeric(argumentos[5])
sig <- as.numeric(argumentos[6])
R <- as.numeric(argumentos[7]) #num de replicas de Monte Carlo
n <- as.numeric(argumentos[8]) # tamanho da amostra

#### Definindo os parâmetros iniciais (para rodar direto pelo rstudio) ####

# pi0 <- 0
# pi1 <- 0
# beta0 <-0 
# beta1 <- 1
# lambda <- 2
# sig <- 0.2
# R <- 500 #num de replicas de Monte Carlo
# n <- 10000  # tamanho da amostra

#### Vetor de parâmetros ####

parametros <- c(beta0, beta1, lambda)

#### Definindo a função de Log-Verossimilhança ####

inicio <- Sys.time()

#fixando a semente
set.seed(1992)

m3_loglik <- function(theta, w, y, sig, n){
  # sig = 0.2
  # n = 10000
  
  #### Definindo expressões e valores para a esp.condicional ####
  
  #theta <- c(0, 1, 2) #vetor para testar sem precisar rodar a funcao m3
  gama <- c(theta[2], theta[3] / sig) #Definir como vetor linha
  mu.w <- cbind(theta[1] * rep(1, n),-theta[3] * w / sig)
  media <- rep(0, 2)
  covariancia <- diag(2) + (sig ^ 2 * (as.matrix(gama) %*% t(gama)))
  up <- cbind(theta[1] + theta[2] * w, rep(0, n))
  
  ### Com base nas contas temos: up_i = mu.w + (gama * wi) = [Beta0 + Beta1 * wi  0]
  
  #### Definindo o resultado da Esperança Condicional ####
  
  # prob <- vector() #inicializando um vetor para armazenar os valores da 'funcao prob'
  # 
  # for (k in 1:n) {
  #   esp <- 2 * mvtnorm::pmvnorm(mean = media, sigma = covariancia, lower = c(-Inf,-Inf), upper = up[k, ])
  #   prob[k] <- esp[1] # funcao p = pi0 + (1 - pi0 - pi1) * E_X|W, neste caso pi0 = pi1 = 0
  # }
  
  prob <- pnorm(theta[1] + theta[2] * w) # Aplicando M1 nos dados do M3
  
  ## Se prob = 1, assumir 0.999999999; Se prob = 0, assumir 0.000000001:
  p <- ifelse(prob == 1, 0.999999999, prob)
  p <- ifelse(prob == 0, 0.000000001, prob)
  loglike <- sum(y * log(p)) + sum((1 - y) * log(1 - p)) #Log-verossimilhanca
  return(loglike)
}  


#vetor de zeros para armazenar as estimativas de cada parâmetro
#emv.pi0 <- rep(0, R)
#emv.pi1 <- rep(0, R)
emv.beta0 <- rep(0, R)
emv.beta1 <- rep(0, R)
emv.lambda <- rep(0, R)

#### Processo de  Monte Carlo ####
# i <- 1 # p/ testar o programar sem rodar o 'for'

for(i in 1:R){
  print(i)
  print(Sys.time() - inicio)
  
  #### Passo 1: Gerar w_i amostras da U(-4, 4) ####
  
  w <- runif(n,-4, 4)
  print(Sys.time() - inicio)  
  
  #### Passo 2: Gerar x_i amostras da Skew Normal ####
  
  x <- rsn(n = n, xi = w, omega = sig ^ 2, alpha = parametros[3])
  print(Sys.time() - inicio)
  
  #### Passo 3: Gerar y_i da Bernoulli ####
  
  y <- rbinom(n = n, size = 1, prob = pnorm(parametros[1] + parametros[2] * x))
  print(Sys.time() - inicio)
  
  #### Passo 4: Gerar Ytil ####
  
  ytil <- y #no modelo 3 temos que Y_T = Y
  
  
  print(Sys.time() - inicio)
  
  #### Passo 5: otimizacao ####
  
  tryCatch(  {
    otimizacao <- optimParallel(
      #par = c(0.001, 0.99, 0.01),
      #par = c(0.001, 0.99, 0.0009), #chutes iniciais para esquema 1 (sig=0.1)
      #par = c(0.001, 0.99, 0.0009), #chutes iniciais para esquema 2 (sig=0.7)
      par = c(0.001, 0.99, 1.98), #chutes iniciais para esquema 3 (sig=0.7)
      fn = m3_loglik,
      method = "L-BFGS-B",
      control = list(fnscale = -1),
      lower = c(-Inf,-Inf,-Inf),
      upper = c(Inf, Inf, Inf),
      w = w,
      y = ytil,
      sig = sig,
      n = n
    )
    ## O R faz minimização por default, então para maximizar devo usar "control=list(fnscale=-1)"
    
    print(Sys.time() - inicio)
    if (otimizacao$convergence == 0) { #0: indica convergencia completa
      #emv.pi0[i] = otimizacao$par[1]
      #emv.pi1[i] = otimizacao$par[2]
      emv.beta0[i] = otimizacao$par[1]
      emv.beta1[i] = otimizacao$par[2]
      emv.lambda[i] = otimizacao$par[3]
    }
    #else{
    # falhas = falhas + 1
    #}
    
  }, error = function(cond)
    NA)
  
}
# fim Monte Carlo

# calculando as médias das estimativas de cada parâmetro
#pi0medio <- mean(emv.pi0)
#pi1medio <- mean(emv.pi1)
beta0medio <- mean(emv.beta0)
beta1medio <- mean(emv.beta1)
lambdamedio <- mean(emv.lambda)

# calculando o erro padrão das estimativas de cada parâmetro
#pi0sd <- sd(emv.pi0)
#pi1sd <- sd(emv.pi1)
beta0sd <- sd(emv.beta0)
beta1sd <- sd(emv.beta1)
lambdasd <- sd(emv.lambda)

#### Calculando viés e erro quadratico medio (eqm) ####
# calculando o viés (vies = media - valor verdadeiro do parametro)
#pi0vies <- pi0medio - pi0
#pi1vies <- pi1medio - pi1
beta0vies <- beta0medio - beta0
beta1vies <- beta1medio - beta1
lambdavies <- lambdamedio - lambda

#pi0viesrel <- pi0vies / pi0
#pi1viesrel <- pi1vies / pi1
beta0viesrel <- beta0vies / beta0
beta1viesrel <- beta1vies / beta1
lambdaviesrel <- lambdavies / lambda

##### EQM ####
# EQM(theta_chapeu) = Var(theta_chaeu) + b²(theta), b²(): viés
# b²(theta) = E(theta_chapeu) - theta

#pi0EQM <- var(emv.pi0) + (pi0vies) ^ 2
#pi1EQM <- var(emv.pi1) + (pi1vies) ^ 2
beta0EQM <- var(emv.beta0) + (beta0vies) ^ 2
beta1EQM <- var(emv.beta1) + (beta1vies) ^ 2
lambdaEQM <- var(emv.lambda) + (lambdavies) ^ 2

# Finalizando a contagem do tempo de execução do programa
fim <- Sys.time()
tempo <- fim - inicio
tempo
#### Lista com os resultados finais ####
resultado <- list(
  Num_obs = n,
  Replicas = R,
  # Pi0 = pi0,
  # Pi1 = pi1,
  Beta_0 = beta0,
  Beta_1 = beta1,
  Lambda = lambda,
  # Pi0_est_medio = pi0medio,
  # Pi1_est_medio = pi1medio,
  Beta0_est_medio = beta0medio,
  Beta1_est_medio = beta1medio,
  Lambda_est_medio = lambdamedio,
  # Pi0_sd = pi0sd,
  # P1_sd = pi1sd,
  Beta0_sd = beta0sd,
  Beta1_sd = beta1sd,
  Lambda_sd = lambdasd,
  # Pi0_est_vies = pi0vies,
  # Pi1_est_vies = pi1vies,
  Beta0_est_vies = beta0vies,
  Beta1_est_vies = beta1vies,
  Lambda_est_vies = lambdavies,
  # Pi0_est_vies_rel = pi0viesrel,
  # Pi1_est_vies_rel = pi1viesrel,
  Beta0_est_vies_rel = beta0viesrel,
  Beta1_est_vies_rel = beta1viesrel,
  Lambda_est_vies_rel = lambdaviesrel,
  # EQM_Pi0 = pi0EQM,
  # EQM_Pi1 = pi1EQM,
  EQM_Beta0 = beta0EQM,
  EQM_Beta1 = beta1EQM, 
  EQM_Lambda = lambdaEQM,
  # Num_Falhas = falhas,
  Tempo_Execução = tempo
)

# Imprimindo os resultados
resultado

# Salvando os resultados em um arquivo
write.csv(x = resultado, file = paste0("geraM3_aplicaM1_r", R,"_n", n, "_l", lambda,"_sig", sig,".csv"))
