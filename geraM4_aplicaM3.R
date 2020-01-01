# Simulação de MC para o modelo 4 (Com erro de classificacao e erro de medida)

# Gera M4 e aplica no M3

require(fExtremes)
require(mvtnorm)
require(sn)
require(optimParallel)

argumentos <- commandArgs(trailingOnly = TRUE) #comando para inserir os valores dos parametros via terminal


cl <- makeCluster(7)     # definindo o num de clusters no processador
setDefaultCluster(cl=cl) # definindo 'cl' como cluster padrao

#### Definindo os parâmetros iniciais (via terminal) ####

pi0 <- as.numeric(argumentos[1])
pi1 <- as.numeric(argumentos[2])
beta0 <- as.numeric(argumentos[3])
beta1 <- as.numeric(argumentos[4])
lambda <- as.numeric(argumentos[5])
sig <- as.numeric(argumentos[6]) # => sig² = 0.01
R <- as.numeric(argumentos[7]) #num de replicas de Monte Carlo
n <- as.numeric(argumentos[8]) # tamanho da amostra

#### Definindo os parâmetros iniciais (para rodar direto pelo rstudio) ####

# beta0 <- 0
# beta1 <- 1
# pi0 <- 0.05
# pi1 <- 0.05
# lambda <- 0.001
# sig <- 0.1 # => sig² = 0.01
# R <- 500 #num de replicas de Monte Carlo
# n <- 20000 # tamanho da amostra

#### Vetor de parâmetros ####

parametros <- c(pi0, pi1, beta0, beta1, lambda)

#### Definindo a função de Log-Verossimilhança ####

#Modelo 4: Há erro de medição e de classificação

inicio <- Sys.time()

#fixando a semente
set.seed(1992)

m4_loglik <- function(theta, w, y){
  sig = 0.7
  n = 10000
  #### Definindo expressões e valores para a esp.condicional ####
  
  #theta <- c(0.1, 0.2, 0, 1, 2) #vetor para testar sem precisar rodar a funcao m4
  gama <- c(theta[4], theta[5] / sig) #Definir como vetor linha
  mu.w <- cbind(theta[3] * rep(1, n),-theta[5] * w / sig)
  media <- rep(0, 2)
  covariancia <- diag(2) + (sig ^ 2 * (as.matrix(gama) %*% t(gama)))
  up <- cbind(theta[3] + theta[4] * w, rep(0, n))
  
  ### Com base nas contas temos: up_i = mu.w + (gama * wi) = [Beta0 + Beta1 * wi  0]
  
  #### Definindo o resultado da Esperança Condicional ####
  
  prob <- vector() #inicializando um vetor para armazenar os valores da 'funcao prob'
  
  for (k in 1:n) {
    esp <- 2 * mvtnorm::pmvnorm(mean = media, sigma = covariancia, lower = c(-Inf,-Inf), upper = up[k, ])
    prob[k] <- esp[1] # Aplicando M3 nos dados do M4
  }
  
  ## Se prob = 1, assumir 0.999999999; Se prob = 0, assumir 0.000000001:
  p <- ifelse(prob == 1, 0.999999999, prob)
  p <- ifelse(prob == 0, 0.000000001, prob)
  loglike <- sum(y * log(p)) + sum((1 - y) * log(1 - p)) #Log-verossimilhanca
  return(loglike)
}  


#vetor de zeros para armazenar as estimativas de cada parâmetro
emv.pi0 <- rep(0, R)
emv.pi1 <- rep(0, R)
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
  
  x <- rsn(n = n, xi = w, omega = sig ^ 2, alpha = parametros[5])
  print(Sys.time() - inicio)
  
  #### Passo 3: Gerar y_i da Bernoulli ####
  
  y <- rbinom(n = n, size = 1, prob = pnorm(parametros[3] + parametros[4] * x))
  print(Sys.time() - inicio)
  p.i <- ifelse(y == 0, pi0, pi1)
  
  uniformes <- runif(n, 0, 1)
  
  comparacao <- ifelse(uniformes < p.i, 1, 0)
  #Comparar cada elemento da Uniforme com o vetor y em (pi0, pi1)
  
  #### Passo 4: Gerar Ytil ####
  
  ytil <- abs(y - comparacao)
  
  print(Sys.time() - inicio)
  
  #### Passo 5: otimizacao ####
  
  tryCatch(  {
    otimizacao <- optimParallel(
      #par = c(0.09, 0.19, 0.001, 0.99, 1.98), #chutes iniciais
      #par = c(0.009, 0.009, 0.001, 0.99, 0.0009), #chutes iniciais para esquema 1 (sig=0.1)
      #par = c(0.03, 0.03, 0.001, 0.99, 0.0009), #chutes iniciais para esquema 2 (sig=0.7)
      par = c(0.03, 0.03, 0.001, 0.99, 1.98), #chutes iniciais para esquema 3 (sig=0.7)
      fn = m4_loglik,
      method = "L-BFGS-B",
      control = list(fnscale = -1),
      lower = c(0, 0,-Inf,-Inf,-Inf),
      upper = c(0.99999999999, 0.99999999999, Inf, Inf, Inf),
      w = w,
      y = ytil
    )
    ## O R faz minimização por default, então para maximizar devo usar "control=list(fnscale=-1)"
    
    print(Sys.time() - inicio)
    if (otimizacao$convergence == 0) { #0: indica convergencia completa
      emv.pi0[i] = otimizacao$par[1]
      emv.pi1[i] = otimizacao$par[2]
      emv.beta0[i] = otimizacao$par[3]
      emv.beta1[i] = otimizacao$par[4]
      emv.lambda[i] = otimizacao$par[5]
    }
    #else{
    # falhas = falhas + 1
    #}
    
  }, error = function(cond)
    NA)
  
}
# fim Monte Carlo

# calculando as médias das estimativas de cada parâmetro
pi0medio <- mean(emv.pi0)
pi1medio <- mean(emv.pi1)
beta0medio <- mean(emv.beta0)
beta1medio <- mean(emv.beta1)
lambdamedio <- mean(emv.lambda)

# calculando o erro padrão das estimativas de cada parâmetro
pi0sd <- sd(emv.pi0)
pi1sd <- sd(emv.pi1)
beta0sd <- sd(emv.beta0)
beta1sd <- sd(emv.beta1)
lambdasd <- sd(emv.lambda)

#### Calculando viés e erro quadratico medio (eqm) ####
# calculando o viés (vies = media - valor verdadeiro do parametro)
pi0vies <- pi0medio - pi0
pi1vies <- pi1medio - pi1
beta0vies <- beta0medio - beta0
beta1vies <- beta1medio - beta1
lambdavies <- lambdamedio - lambda

pi0viesrel <- pi0vies / pi0
pi1viesrel <- pi1vies / pi1
beta0viesrel <- beta0vies / beta0
beta1viesrel <- beta1vies / beta1
lambdaviesrel <- lambdavies / lambda

##### EQM ####
# EQM(theta_chapeu) = Var(theta_chaeu) + b²(theta), b²(): viés
# b²(theta) = E(theta_chapeu) - theta

pi0EQM <- var(emv.pi0) + (pi0vies) ^ 2
pi1EQM <- var(emv.pi1) + (pi1vies) ^ 2
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
  Pi0 = pi0,
  Pi1 = pi1,
  Beta_0 = beta0,
  Beta_1 = beta1,
  Lambda = lambda,
  Pi0_est_medio = pi0medio,
  Pi1_est_medio = pi1medio,
  Beta0_est_medio = beta0medio,
  Beta1_est_medio = beta1medio,
  Lambda_est_medio = lambdamedio,
  Pi0_sd = pi0sd,
  P1_sd = pi1sd,
  Beta0_sd = beta0sd,
  Beta1_sd = beta1sd,
  Lambda_sd = lambdasd,
  Pi0_est_vies = pi0vies,
  Pi1_est_vies = pi1vies,
  Beta0_est_vies = beta0vies,
  Beta1_est_vies = beta1vies,
  Lambda_est_vies = lambdavies,
  Pi0_est_vies_rel = pi0viesrel,
  Pi1_est_vies_rel = pi1viesrel,
  Beta0_est_vies_rel = beta0viesrel,
  Beta1_est_vies_rel = beta1viesrel,
  Lambda_est_vies_rel = lambdaviesrel,
  EQM_Pi0 = pi0EQM,
  EQM_Pi1 = pi1EQM,
  EQM_Beta0 = beta0EQM,
  EQM_Beta1 = beta1EQM, 
  EQM_Lambda = lambdaEQM,
  # Num_Falhas = falhas,
  Tempo_Execução = tempo
)

# Imprimindo os resultados
resultado

# Salvando os resultados em um arquivo
write.csv(x = resultado, file = paste0("geram4_aplicam3_r", R, "_n", n, "_l", lambda, "_sig", sig, ".csv"))
