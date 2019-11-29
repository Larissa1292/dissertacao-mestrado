## Aplicacao:

#Data of Sposto et al. "The effect os diagnostic misclassification on non cancer and 
# cancer mortality dose response in A-bomb survivors"

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

## Juntando os w_i em um único W, onde i = 20:

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

modelo <- glm(Y ~ W, family = binomial(link = "probit"), data = banco2)
summary(modelo)
