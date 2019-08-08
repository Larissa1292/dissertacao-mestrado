#M3 - Gearação dos dados

parametros.m3 <- c( beta0, beta1, lambda)

n=10000

w <- runif(n,-4, 4)
#print(Sys.time() - inicio)  
#### Passo 2: Gerar x_i amostras da Skew Normal ####

x <- rsn(n = n, xi = w, omega = sig ^ 2, alpha = parametros.m3[3])
#print(Sys.time() - inicio)
#### Passo 3: Gerar y_i da Bernoulli ####

y <- rbinom(n = n, size = 1, prob = pnorm(parametros.m3[1] + parametros.m3[2] * x))
print(Sys.time() - inicio)

