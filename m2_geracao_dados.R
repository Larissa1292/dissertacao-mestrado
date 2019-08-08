# M2 - Geração dos dados

parametros <- c(pi0, pi1, beta0, beta1)
#### Passo 1: Gerar w_i amostras da U(-4, 4) ####

w <- runif(n,-4, 4)

#### Passo 2: Gerar x_i amostras da Skew Normal ####

#x <- rsn(n = n, xi = w, omega = sig ^ 2, alpha = parametros[5])

x <- w
#### Passo 3: Gerar y_i da Bernoulli ####

y <- rbinom(n = n, size = 1, prob = pnorm(parametros[3] + parametros[4] * x))

p.i <- ifelse(y == 0, pi0, pi1)

uniformes <- runif(n, 0, 1)

comparacao <- ifelse(uniformes < p.i, 1, 0)
#Comparar cada elemento da Uniforme com o vetor y em (pi0, pi1)

#### Passo 4: Gerar Ytil ####

ytil <- abs(y - comparacao)
