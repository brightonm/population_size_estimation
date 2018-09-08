########## PARTIE 1 ##############

rm(list=ls()) #permet de vider la liste des variables

# Question 1 : simulation de l'échantillon (x1, ..., xn)
# On a décidé de créer un fonction pour faire varier
# plus facilement les valeurs de n, theta et n_0

simulation_echantillon <- function (theta, n_0, n) {
    prob <- n_0 / theta 
    echantillon <- rbinom(n, 1, prob) # simulation d'un échantillon de loi de bernoulli de paramètre prob
    return(echantillon)
}

# Test de la fonction avec, comme indiqué en question 1
# theta = 1000, n_0 = 50 et on choisit n = 100

theta <- 1000
n_0 <- 50
n <- 100

echantillon1 <- simulation_echantillon(theta, n_0, n)
echantillon1

cat("La moyenne empirique de l'échantillon est : ", mean(echantillon1))
cat("Alors que son espérance théorique est : ", n_0/theta)

cat("La variance empirique de l'échantillon est : ", ((n-1)/n)*var(echantillon1))
cat("Alors que sa variance théorique est : ", (n_0/theta)*(1-(n_0/theta)))


# Question 2 : Réalisation de T (variable aléatoire binomiale) sur l'exemple simulé 

t <- sum(echantillon1)
cat("La réalisation de la loi du nombre total de poissons bagués parmi les n poissons pêchés est : ", t)

# Question 3 : Estimation de theta^ et de theta~ qui sont les mêmes

estimation_emm_emv <- (n * n_0) / t
cat("L'estimation du maximum de vraisemblance et des moments du paramètre theta est donc : ", estimation_emm_emv)

# Question 4 : intervalles de confiance exacts et asymptotique de différents seuils pour theta
# d'après la formule du cours en utilisant les quantiles de la loi de Fischer-Snedecor

ic_exact <- function (alpha, echantillon, n_0) {
    n <- length(echantillon)
    t <- sum(echantillon)
    borne_inf <- n_0 * (1 + ((n-t+1)/t)*qf(alpha/2, 2*(n-t+1), 2*t)) # quantiles de la loi de fisher
    borne_sup <- n_0 * (1 + ((n-t)/(t+1))*qf(1-alpha/2, 2*(n-t), 2*(t+1)))
    return(c(borne_inf, borne_sup))
}

ic_exact_1 <- ic_exact(.01, echantillon1, n_0)
cat("L'intervalle de confiance exact pour theta de seuil de 1% est : ", ic_exact_1)

ic_exact_5 <- ic_exact(.05, echantillon1, n_0)
cat("L'intervalle de confiance exact pour theta de seuil de 5% est : ", ic_exact_5)

ic_exact_10 <- ic_exact(.10, echantillon1, n_0)
cat("L'intervalle de confiance exact pour theta de seuil de 10% est : ", ic_exact_10)

ic_exact_20 <- ic_exact(.20, echantillon1, n_0)
cat("L'intervalle de confiance exact pour theta de seuil de 20% est : ", ic_exact_20)

ic_asymptotique <- function(alpha, echantillon, n_0) {
    n <- length(echantillon)
    t <- sum(echantillon)
    p <- t / n 
    borne_inf <- n_0 / (p + qnorm(1-alpha/2)*sqrt(p*(1-p)/n)) 
    borne_sup <- n_0 / (p - qnorm(1-alpha/2)*sqrt(p*(1-p)/n))
    return(c(borne_inf, borne_sup))
}

ic_asymptotique_1 <- ic_asymptotique(.01, echantillon1, n_0)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 1% est : ", ic_asymptotique_1)

ic_asymptotique_5 <- ic_asymptotique(.05, echantillon1, n_0)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 5% est : ", ic_asymptotique_5)

ic_asymptotique_10 <- ic_asymptotique(.10, echantillon1, n_0)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 10% est : ", ic_asymptotique_10)

ic_asymptotique_20 <- ic_asymptotique(.20, echantillon1, n_0)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 20% est : ", ic_asymptotique_20)

# Question 5 : Calcul de P(estimateur = +inf)

proba_estimateur_inf <- (1 - n_0/theta)**n
cat("P(estimateur = +inf) =", proba_estimateur_inf)

# Question 6 : Simulations mettant en évidence le fait que P(estimateur = +inf) > 1/2
# pour n appartenant à [|0, 13|]
# Remarque : en R ln = log() et log en base 10 = log10()

sup_intervalle = floor(-log(2)/log(1-n_0/theta))

cat("On a P(estimateur = +inf) > 1/2 pour n appartenant à l'intervalle suivant : ", c(0, sup_intervalle))

cat("Testons en faisant des simulations sur notre exemple : ")

for (i in seq(from = 0, to = sup_intervalle + 2, by = 1))
    cat("Pour n = ", i , "on a P(estimateur = +inf) =", (1 - n_0/theta)**i, "\n")

############ PARTIE 2 ################

rm(list=ls()) #permet de vider la liste des variables

# Question 1 : Simulation de l'échantillon (y1, y2, ..., ym)

simulation_echantillon2 <- function (theta, n_0, m) {
  prob <- n_0 / theta 
  echantillon <- rgeom(m, prob) + 1 # simulation d'un échantillon m réalisations de v.a de loi géométrique
  return(echantillon)
}

n_0 <- 50 
theta <- 1000
m <- 100 #on choisit m = 100
echantillon2 = simulation_echantillon2(theta, n_0, m)
echantillon2

cat("La moyenne empirique de l'échantillon est : ", mean(echantillon2))
cat("Alors que son espérance théorique est : ", theta/n_0)

cat("La variance empirique de l'échantillon est : ", ((m-1)/m)*var(echantillon2))
cat("Alors que sa variance théorique est : ", (1-(n_0/theta))/(n_0/theta)**2)

# Question 2 : Réalisation de N (variable aléatoire binomiale négative) sur l'exemple simulé 

n2 <- sum(echantillon2)
cat("La réalisation de la loi du nombre total de poissons pêchés est : ", n2)

# Question 3 : Estimation de theta^ et de theta~ qui sont les mêmes

estimation_emm_emv2<- (n_0 * n2) / m
cat("L'estimation du maximum de vraisemblance et des moments du paramètre theta est donc : ", estimation_emm_emv2)

# Question 5 : Calculs d'intervalles de confiance

info_fisher <- m / ((estimation_emm_emv2 - n_0)*estimation_emm_emv2)



ic_asymptotique2 <- function(alpha) {
  borne_inf <- estimation_emm_emv2 - qnorm(1-alpha/2)/sqrt(info_fisher) 
  borne_sup <- estimation_emm_emv2 + qnorm(1-alpha/2)/sqrt(info_fisher)
  return(c(borne_inf, borne_sup))
}

ic_asymptotique2_1 <- ic_asymptotique2(.01)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 1% est : ", ic_asymptotique2_1)

ic_asymptotique2_5 <- ic_asymptotique2(.05)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 5% est : ", ic_asymptotique2_5)

ic_asymptotique2_10 <- ic_asymptotique2(.10)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 10% est : ", ic_asymptotique2_10)

ic_asymptotique2_20 <- ic_asymptotique2(.20)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 20% est : ", ic_asymptotique2_20)


############ PARTIE 3 ###############

rm(list=ls()) #permet de vider la liste des variables

# Question 1 

n_0 <- 100
n <- 1000
echantillon3 <- scan('Peche.txt')
t <- sum(echantillon3)
cat("La réalisation de la loi du nombre total de poissons bagués parmi les n poissons pêchés est : ", t)

estimation_emm_emv <- (n * n_0) / t
cat("L'estimation du maximum de vraisemblance et des moments du paramètre theta est donc : ", estimation_emm_emv)

ic_exact <- function (alpha, echantillon, n_0) {
  n <- length(echantillon)
  t <- sum(echantillon)
  borne_inf <- n_0 * (1 + ((n-t+1)/t)*qf(alpha/2, 2*(n-t+1), 2*t)) #quantiles de la loi de fisher
  borne_sup <- n_0 * (1 + ((n-t)/(t+1))*qf(1-alpha/2, 2*(n-t), 2*(t+1)))
  return(c(borne_inf, borne_sup))
}

ic_exact5 <- ic_exact(.05, echantillon3, n_0)
cat("L'intervalle de confiance exact pour theta de seuil de 5% selon la première stratégie est : ", ic_exact5)

ic_asymptotique <- function(alpha, echantillon, n_0) {
  n <- length(echantillon)
  t <- sum(echantillon)
  p <- t / n 
  borne_inf <- n_0 / (p + qnorm(1-alpha/2)*sqrt(p*(1-p)/n)) 
  borne_sup <- n_0 / (p - qnorm(1-alpha/2)*sqrt(p*(1-p)/n))
  return(c(borne_inf, borne_sup))
}

ic_asymptotique_5 <- ic_asymptotique(.05, echantillon3, n_0)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 5% est : ", ic_asymptotique_5)

# Question 2

transformation_x_vers_y <- function(vecteur) {
  n <- length(vecteur)
  echecs_en_cours <- 0 # variable comptant le nombre de poissons non bagués péchés en cours d'affilé 
  y <- c() # vecteur contenant les futurs réalisations de Y1, ..., Yn
  for (i in seq(1:n)) {
    echecs_en_cours <- echecs_en_cours + 1
    if (vecteur[i] == 1) {
      y <- c(y, echecs_en_cours)
      echecs_en_cours <- 0
    }
  }
  return(y)
}

echantillon4 <- transformation_x_vers_y(echantillon3)

m <- length(echantillon4)
n2 <- sum(echantillon4)
estimation_emm_emv2<- (n_0 * n2) / m
cat("L'estimation du maximum de vraisemblance et des moments du paramètre theta est donc : ", estimation_emm_emv2)

info_fisher <- m / ((estimation_emm_emv2 - n_0)*estimation_emm_emv2)



ic_asymptotique2 <- function(alpha) {
  borne_inf <- estimation_emm_emv2 - (qnorm(1-alpha/2)/sqrt(info_fisher))
  borne_sup <- estimation_emm_emv2 + (qnorm(1-alpha/2)/sqrt(info_fisher))
  return(c(borne_inf, borne_sup))
}

ic_asymptotique2_5 <- ic_asymptotique2(.05)
cat("L'intervalle de confiance asymptotique pour theta de seuil de 5% est : ", ic_asymptotique2_5)

# Question 3 : Graphe de probabilités pour la loi géométrique

# on enlève la dernière valeur pour éviter une erreur en calculant log(0)

echantillon4ord <- sort(echantillon4)
plot(echantillon4ord[1:(m-1)], log(1-seq(1:(m-1))/m), col='red', col.lab=rgb(0,0.5,0), ylab="Xi*", xlab="log(1-i/m)")
title("Graphe de probabilités loi géométrique")
# Superposition de la droite des moindres carres

abs<-echantillon4ord[1:(m-1)]
ord<-log(1 - seq(1:(m-1))/m)
reg<-lm(ord~abs)
lines(abs, fitted.values(reg))

############## PARTIE 4 #################

rm(list=ls()) #permet de vider la liste des variables

# Question 1 

ic_exact <- function (alpha, echantillon, n_0) {
  n1 <- length(echantillon)
  t1 <- sum(echantillon)
  borne_inf <- n_0 * (1 + ((n1-t1+1)/t1)*qf(alpha/2, 2*(n1-t1+1), 2*t1)) # quantile de la loi de Fisher
  borne_sup <- n_0 * (1 + ((n1-t1)/(t1+1))*qf(1-alpha/2, 2*(n1-t1), 2*(t1+1)))
  return(c(borne_inf, borne_sup))
}

simulation <- function (theta, n_0, n, m, alpha) {
  nb_bon_ic <- 0
  for (i in seq(1:m)) {
    echantillon <- rbinom(n, 1, n_0/theta)
    ic <- ic_exact(alpha, echantillon, n_0)
    # empecher l'arret du programme avec la production de NaN fréquente
    # de qf(alpha/2, 2 * (n1 - t1 + 1), 2 * t1) :
    if (!(is.nan(ic[1]) || is.nan(ic[2]))) {
      if ((theta >= ic[1]) && (theta <= ic[2])) {
        nb_bon_ic <- nb_bon_ic + 1
      }
    }
  }
  return((nb_bon_ic/m)*100) # on retourne le résultat en pourcentage
}

# Différents tests pour tirer des conclusions

accuracy <- simulation(1000, 50, 10000, 1000, .20)
cat("La proportion d'intervalles contenant la vraie valeur de theta à 95% est :", accuracy, "%")

print(simulation(1000,50,1000,1000,.10))
print(simulation(1000,50,100,1000,.01))

# Question 2 

loi_faible_grds_nb <- function (m, p, n, epsilon) {
  nb_sup_eps <- 0 # variable qui compte le nombre fois que la différence est supérieur à epsilon
  for (i in seq(1:m)) {
    echantillon <- rbinom(n, 1, p)
    if (abs(mean(echantillon) - p) >= epsilon) {
      nb_sup_eps <- nb_sup_eps + 1
    }
  }
  return(nb_sup_eps)
}

# Tests, variation du paramètre n

# On choisit epsilon = 0,001, p = 0,05, m = 100 et on trace en nb_sup_eps fonction de n

ord <- c()
abs <- c()
for (i in seq(1:6)) {
  abs <- c(abs, 10**i)
  ord <- c(ord, loi_faible_grds_nb(100, .05, 10**i, .001))
}

plot(abs, ord, log='x', type='o', col="red", xlab='n', col.lab=rgb(0,0.5,0), ylab="Nombre d'écarts supérieurs à epsilon (m=100)", col.lab=rgb(0, 0.5, 0))
title("Simulation loi faible des grands nombres")

# Question 3

theo_central_limite <- function (m, p, n) {
  vecteur_moyenne <- c()
  for (i in seq(1:m)) {
    echantillon <- rbinom(n, 1, p)
    vecteur_moyenne <- c(vecteur_moyenne, mean(echantillon))
  }
  return(vecteur_moyenne)
}

histolarg <- function(x, xlim=NULL, ...)
{
  # nombre de donnees
  n <- length(x) 
  # nombre de classes (regle de Sturges)
  if (n<12) k<-5 else k <- round(log2(n)+1) 
  # bornes des classes
  rangex <- max(x)-min(x)
  a0 <- min(x)-0.025*rangex
  ak <- max(x)+0.025*rangex
  bornes <- seq(a0, ak, length=k+1)
  # etendue du trace
  if (is.null(xlim))
    xlim<-c(bornes[1], bornes[k+1])
  # histogramme
  histx<-hist(x, prob=T, breaks=bornes, xlim=xlim, ...)
  # histx
}

# On fait varier n, histogramme qui ressemble à une loi normale

histolarg(theo_central_limite(10000, .05, 10000))

# Graphe de probabilités 

echantillon <- theo_central_limite(1000, .05, 1000)
plot(sort(echantillon)[1:999], qnorm(seq(1:999)/1000), col='red', xlab='Xi*', col.lab=rgb(0, 0.5, 0), ylab='inv_phi(i/999)', col.lab=rgb(0, 0.5, 0))
title("Graphe de probabilités de la loi normale")

    