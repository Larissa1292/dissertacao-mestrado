# Ajustes (Especificados no Trello) de acordo com a reuniao do dia 19/06


require(fExtremes)
require(mvtnorm)
require(sn)


#### Definindo os parâmetros iniciais ####

beta0 <- 0
beta1 <- 1
pi0 <- 0.1
pi1 <- 0.2
lambda <- 2
sig <- 0.2
R <- 100 #num de replicas de Monte Carlo
n <- 1000 # tamanho da amostra
#ytil <- matrix(0, nrow = n, ncol = R)
#opt1 <- matrix(0, R, n) #col = num de parametros a serem estimados

# pi0 <- 0.35       ## SAO OS VALORES DA SURUPA?
# pi1 <- 0.20
# beta1_true <- 1   ## ??? se ja definiu beta0 e beta1 !
# beta0_true <- -1   ###????

#### Vetor de parâmetros ####

parametros <- c(pi0, pi1, beta0, beta1, lambda)

#### Definindo a função de Log-Verossimilhança ####

#Modelo 4: Há erro de medição e de classificação

inicio <- Sys.time()

#fixando a semente
set.seed(1992)

m4_loglik <- function(theta, w, y){
  
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
    esp <- 2 * pmvnorm(mean = media, sigma = covariancia, lower = c(-Inf,-Inf), upper = up[k, ])
    prob[k] <- theta[1] + (1 - theta[1] - theta[2]) * esp[1] # funcao p = pi0 + (1 - pi0 - pi1) * E_X|W
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

#inicializando o contador de falhas
#falhas <- 0

#### Processo de  Monte Carlo ####
# i <- 1 # p/ testar o programar sem rodar o 'for'

for (i in 1:R) {
  print(i)
  
  #### Passo 1: Gerar w_i amostras da U(-4, 4) ####
  
  w <- runif(n,-4, 4)
  
  #### Passo 2: Gerar x_i amostras da Skew Normal ####
  
  x <- rsn(n = n, xi = w, omega = sig ^ 2, alpha = parametros[5])
  
  #### Passo 3: Gerar y_i da Bernoulli ####
  
  ##  ???????? ERRADO! REFAZER
  # p <-  pnorm(parametros[3] + parametros[4] * w) #probit
  # y <-  rbinom(n = length(x), size = 1, prob = p)   ### QUEM E ESTE Y ??????????
  # 
  #### Passo 4: Gerar Ytil ####
  
  # Generate the true binary response y_true, with covariate x
  
  lm <-  beta0 + beta1 * x    ###  AQUI TEM QUE SER beta0 e beta 1.... POR QUE DIFERENTES BETAS? ok
  pr.probit <- pnorm(lm)  
  y_true <- rbinom(n, 1, pr.probit)
  
  # # Generate the misclassifed variable 
  # 
  pr_pi0.probit <- pi0
  alpha0.probit <- rbinom(n, 1, pr_pi0.probit)  # alpha0=(Y=1|Y_T=0)
  
  pr_pi1.probit <- 1 - pi1
  alpha1.probit <- rbinom(n, 1, pr_pi1.probit)  # alpha1=(Y=1|Y_T=1)
  y <- vector()  ### ytil ? Y OBSERVADO!!!! 
  
  for(i in 1:n){
    y[i] <- ifelse(y_true[i]==1, alpha1.probit[i], alpha0.probit[i])
  }
  
  # QUAL A OUTRA MANEIRA DE GERAR Y ?
  # Considerando a prob de sucesso P(Y=1|W) =  pi0 + (1 - pi0 - pi1) * E_x|W{\Phi(beta0 + beta1*x)}
  
  
  #### Calculando a log-verossimilhanca para cada n ####
  m4_n <- function(theta) {
    m4_loglik(theta, w, y)   #### POR QUE UTILIZA X??  ESSE NAO SE OBSERVA!!! (ok)
  }
  
  
  #### Passo 5: otimizacao ####
  
  tryCatch({
    otimizacao <- optim(
      par = c(0.09, 0.1, 0.1, 0.9, 1.9),
      fn = m4_n,
      method = "L-BFGS-B",
      control = list(fnscale = -1),
      lower = c(0, 0,-Inf,-Inf,-Inf),
      upper = c(0.99999999999, 0.99999999999, Inf, Inf, Inf)
    )
    ## O R faz minimização por default, então para maximizar devo usar "control=list(fnscale=-1)"
    
    
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
  
  
} # fim Monte Carlo

# calculando as médias das estimativas de cada parâmetro
pi0medio <- mean(emv.pi0)
pi1medio <- mean(emv.pi1)
beta0medio <- mean(emv.beta0)
beta1medio <- mean(emv.beta1)
lambdamedio <- mean(emv.lambda)

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
