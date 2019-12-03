## Aplicacao:

require(sn)

#Data of Sposto et al. "The effect os diagnostic misclassification on non cancer and 
# cancer mortality dose response in A-bomb survivors"

#### Abordagem 1 ####

x <- rep(0, 2)
y1 <- rep(1, 2784)
y2 <- rep(0, 10201)ww
y <- c(y1, y0)
head(y)
tail(y)


#### Abordagem 2 ####

x <- c(0, 0.018, 0.072, 0.137, 0.324, 0.693, 1.350, 2.350, 3.520, 4.430) #concentracao media de radiacao absorvida
y1 <- c(2784, 2105, 439, 523, 586, 339, 204, 57, 21, 13) # y = 1, morte por cancer
y0 <- c(10201, 7451, 1509, 1701, 1785, 826, 369, 86, 51, 23) # y = 0, nao morte por cancer
p <- c(0.2144, 0.2203, 0.2253, 0.2352, 0.2471, 0.2910, 0.3560, 0.3986, 0.2917, 0.3611) # proporcao de mortes por cancer

dados <- data.frame(x, y1, y0, p)

lm(y1 ~ x, data = dados)

#####

X <- dados[,1]
Y <- dados[,c(2,3)]
mod <- lm(Y~X)
mod <- lm(c(Y[,1],c[Y[,2]])~X)

## Exemplo do help do r
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)
model.matrix(lm.D9)

#############################

#Mean dose

w1 <- rep(0, 2784)
w2 <- rep(0, 10201)
w3 <- rep(1, 2105) * 0.018
w4 <- rep(1, 7451) * 0.018
w5 <- rep(1, 439) * 0.072
w6 <- rep(1, 1509) * 0.072
w7 <- rep(1, 523) * 0.137
w8 <- rep(1, 1701) * 0.137
w9 <- rep(1, 586) * 0.324
w10 <- rep(1, 1785) * 0.324
w11 <- rep(1, 339) * 0.693
w12 <- rep(1, 826) * 0.693
w13 <- rep(1, 204) * 1.350
w14 <- rep(1, 369) * 1.350
w15 <- rep(1, 57) * 2.350
w16 <- rep(1, 86) * 2.350
w17 <- rep(1, 21) * 3.520
w18 <- rep(1, 51) * 3.520
w19 <- rep(1, 13) * 4.430
w20 <- rep(1, 23) * 4.430

## Var. resposta (Y)

y1 <- rep(1, 2784)
y2 <- rep(0, 10201)
y3 <- rep(1, 2105)
y4 <- rep(0, 7451)
y5 <- rep(1, 439)
y6 <- rep(0, 1509)
y7 <- rep(1, 523)
y8 <- rep(0, 1701)
y9 <- rep(1, 586)
y10 <- rep(0, 1785)
y11 <- rep(1, 339)
y12 <- rep(0, 826)
y13 <- rep(1, 204)
y14 <- rep(0, 369)
y15 <- rep(1, 57)
y16 <- rep(0, 86)
y17 <- rep(1, 21)
y18 <- rep(0, 51)
y19 <- rep(1, 13)
y20 <- rep(0, 23)

## Juntando os y_i em um único Y, onde i = 20:

Y <- matrix(c(y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, y17, y18, 
              y19, y20), ncol = 1)

## Juntando os w_i em um único Y, onde i = 20:

W <- matrix(c(w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17, w18, 
             w19, w20), ncol = 1)

length(W) 
length(Y)

#Atribuindo nomes as colunas da matriz:

row_names <- c()
col_names <- c("W", "Y")

banco <- matrix(c(W,Y), nrow = length(W), ncol = 2, byrow = F, dimnames = list(row_names, col_names))

#Transformando o banco em data frame:
banco2 <- as.data.frame(banco)
class(banco2)

#### Modelo 1: Sem erro de classificação e sem erro de medida #### ok

set.seed(1992)

modelo1 <- glm(Y ~ W, family = binomial(link = "probit"), data = banco2)
summary(modelo1)

#### Modelo 2: Somente erro de classificação (X = W) ####

w <- banco2$W

x <- w

modelo2 <- glm(Y ~ W, family = binomial(link = "probit"), data = banco2)
summary(modelo2)

#### Modelo 3: Somente erro de medida (Y_t = Y) #### ok
# parametros <- c(beta0, beta1, lambda)

sig <- 0.1
lambda <- 2 #verificar como estimar

w <- banco2$W
x <- rsn(n = length(banco2$W), xi = w, omega = sig ^ 2, alpha = lambda)
banco2$X <- x

modelo3 <- glm(Y ~ X, family = binomial(link = "probit"), data = banco2)
summary(modelo3)

#### Modelo 4: Com erro de classificação e com erro de medida ####
sig <- 0.1
lambda <- 2 #verificar como estimar
pi0 <- 0.1
pi1 <- 0.2

w <- banco2$W
x <- rsn(n = length(banco2$W), xi = w, omega = sig ^ 2, alpha = lambda)
y.t <- banco2$Y

#gerando Y (observado)
p.i <- ifelse(y.t == 0, pi0, pi1)
uniformes <- runif(length(banco2$W), 0, 1)
comparacao <- ifelse(uniformes < p.i, 1, 0)
y <- abs(y.t - comparacao)
banco2$Ytil <- y 


modelo4 <- glm(Ytil ~ X, family = binomial(link = "probit"), data = banco2)
summary(modelo4)



##### Novas ideias sobre a origem dos dados para os diferentes modelos



#O parâmetro \beta_0 corresponde ao intercepto do plano com o eixo z. 
#Se x=(x_1, x_2)=(0,0) o parâmetro beta_0 fornece a resposta média nesse ponto.
#Caso contrário, não é possível interpretar o parâmetro beta_0.

#O parâmetro \beta_1 indica uma mudança na resposta média a cada unidade de mudança em 
# x_1, quando as demais variáveis são mantidas fixas.
