# dissertacao-mestrado
Programa para simulação de Monte Carlo

Objetivo: Montar uma simulação de Monte Carlo para obter as estimativas dos Parâmetros.
Parâmetros a serem estimados: pi0, pi1, beta0, beta1, lambda.

Etapas:
1) Definir parâmetros iniciais;
2) Definir vetor de parâmetros;
3) Definir a função de Log-verossimilhança: para isso temos que definir antes E_x|w{ Phi(Beta0 + Beta1 * X)}, pois será necessário para definir a P(Y=1|W) = pi0 + (1 - pi0 - pi1) * E_x|w{ Phi(Beta0 + Beta1 * X)}.
    - Definir a acumulada da Normal Bivariada (definir: Gama, mu(W), média, covariância, limite superior da integral (up[]));
    - Definir P(Y=1|W) = pi0 + (1 - pi0 - pi1) * E_x|w{ Phi(Beta0 + Beta1 * X)};
    - Definir a função de Log-Verossimilhança.
4) Definir os vetores vazios para armazenar as estimativas dos parâmetros;
5) Iniciar processo de Monte Carlo (MC):
    - Passo 1: Gerar w_i amostras da U(-4, 4);
    - Passo 2: Gerar x_i amostras da Skew Normal;
    - Passo 3: Gerar y_i da Bernoulli
    - Passo 4: Gerar Ytil (y observado)
    - Passo 5: Iniciar a otimizacao
6) Calcular as médias das estimativas;
7) Calcular o desvio padrão, viés e EQM;
