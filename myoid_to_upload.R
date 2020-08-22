#library
library(rethinking)

#model for age
d <- list(A = Table_matrix$Age)
m1age <- ulam(
  alist(
    A ~ dnorm( mu, sigma),
    mu ~ dnorm (60, 20),
    sigma ~ dexp(1)
  ), data = d, chains = 4)

#model for size
d <- list(S = Table_matrix$Size)
m1size <-ulam(
  alist(
    S ~ normal( mu, sigma),
    mu ~ dnorm (5, 0.75),
    sigma ~ dexp(1)
  ), data = d, chains = 4)

#Model for mitosis
i <- !is.na(Table_matrix$Mitotic_count_10HPF)
d <- list(M = as.integer(Table_matrix$Mitotic_count_10HPF[i]))
m1mitotic <-ulam(
  alist(
    M ~ dzipois( p, lambda),
    logit(p) <- ap,
    log(lambda) <- al,
    ap ~ dnorm(-1.5, 1),
    al ~ dnorm (0.5, 0.5)
  ), data = d, chains = 4)