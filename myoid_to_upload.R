# Libraries
library(igraph) #to graph network structure
library(rethinking) #a package for doing Bayesian statistics
library(tidyverse) # we will use this for data wrangling

# Analysis

## Network analysis

### Loading the databases
Net <- read.csv("Net.csv")
nodes <- read.csv("nodes.csv")

### Constructicng the citation network object
net <- graph_from_data_frame(Net,vertices = nodes,directed = T)
net$layout <- layout_in_circle # Plotting citation network in circle

### Make the figure
plot(net, vertex.shape="none",  
     vertex.label.font=2, vertex.label.color="black",
     vertex.label.cex=.7, edge.color="gray85")


## Bayesian Analysis

### Loading the database
Table_matrix <- read.csv("Table_matrix.csv")

### Model for Age

#### Prior predictive simulation for age
new_title <- c("Normal( 60, 10)", "Normal( 60, 15)", "Normal( 60, 20)")
par(mfrow=c(1,3) ) # 3 plots in 1 row 
for ( d in 1:3) 
{ sigma <- c(10,15,20) 
sample_mu <- rnorm( 1e4 , 60 , sigma[d] )
sample_sigma <- rexp( 1e4 , 1 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h, xlim= c(0,120), main = new_title[d], ylim=c(0,0.041), xlab = "Age, years")
abline(v=c(0,120), col=c("blue", "red"), lty=c(1,2), lwd=c(1, 3))
}

#### Probabilistic Model
d <- list(A = Table_matrix$Age)
m1age <- ulam(
  alist(
    A ~ dnorm( mu, sigma),
    mu ~ dnorm (60, 20),
    sigma ~ dexp(1)
  ), data = d, chains = 4)

#### Model ceofficients
precis(m1age,2)

#### Posterior predicitive simulations
post <- extract.samples( m1age , n=1e4 )
sim.age <- rnorm(1000, post$mu, post$sigma)
p <- c(d$A[1:16],sim.age)
plo <- sample(p,1016)
sim.mean <- mean(post$mu)
sim.mean.sd <- c(sim.mean+sd(post$mu),sim.mean-sd(post$mu))
sim.sd.mean <- c(mean(post$mu)+2*mean(post$sigma),mean(post$mu)-2*mean(post$sigma))
idx <- rep(0, 1016) 
idx <- ifelse(plo %in% d$A[1:4], 1,
              ifelse(plo %in% d$A[5:16] , 2, 0))
idx <- as.factor(idx)

#### Plotting the posterior predictive check comparing predicted data with actual data 
par(mfrow=c(1,1))
plot(plo, col = c(col.alpha(rangi2,0.5) ,"black" ,"black")[idx],
     pch = c(1,19,21)[idx],  bg = "lightgray", cex = 1.2,
     xaxt="n", main= "Age of patients with Myoid tumors",
     ylab = "Age, years",ylim = c(0,90), xlab = "")
abline(h=c(sim.mean,sim.mean.sd), col=c("blue"), lty=c(1,2,2)) # mean and sd of mean
abline(h=c(sim.sd.mean), col=c("gray"), lty=c(2,2))       # 2*sigma in gray
#black solid: patient form our series
#gray dots: patients from the literatures
#blue shaded dots simulated patients

#### Prior predictive simulation for size
par(mfrow=c(1,3))
new_title <- c("Normal(5, 2)", "Normal(5, 1)", "Normal(5, 0.5)")
par(mfrow=c(1,3) )
for ( d in 1:3) 
{ sigma <- c(2,1,0.5) 
sample_mu <- rnorm( 1e4 ,5 , sigma[d] )
sample_sigma <- rexp( 1e4 , 1 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h, xlim= c(-3,10), main = new_title[d],  xlab = "Size, cm")
abline(v=c(0,120), col=c("blue", "red"), lty=c(1,2), lwd=c(1, 3))
}

#### Probabilistic Model
d <- list(S = Table_matrix$Size)
m1size <-ulam(
  alist(
    S ~ normal( mu, sigma),
    mu ~ dnorm (5, 0.75),
    sigma ~ dexp(1)
  ), data = d, chains = 4)

#### Model ceofficients
precis(m1size,2)

#### Posterior predicitive simulations
post <- extract.samples( m1size , n=1e4 )
p.1 <- c(d$S[1:16],rnorm(1000, post$mu, post$sigma))
p.2 <- c(rep(1, 4),rep(2,12), rep(0, 1000))
plo <- sample(1:1016)
idx <- as.factor(p.2[plo])
sim.mean <- mean(post$mu)
sim.mean.sd <- c(sim.mean+sd(post$mu),sim.mean-sd(post$mu))
sim.sd.mean <- c(mean(post$mu)+2*mean(post$sigma),mean(post$mu)-2*mean(post$sigma))

#### Plotting the posterior predictive check comparing predicted data with actual data 
par(mfrow=c(1,1))
plot(p.1[plo], col = c(col.alpha(rangi2,0.5) ,"black" ,"black")[idx],
     pch = c(1,19,21)[idx],  bg = "lightgray", cex = 1.2,
     xaxt="n", main= "Size of Myoid tumors",
     ylab = "Size, cm",ylim = c(0,10), xlab = "")
abline(h=c(sim.mean,sim.mean.sd), col=c("blue"), lty=c(1,2,2))            # mean and sd of mean
abline(h=c(sim.sd.mean), col=c("gray"), lty=c(2,2)) #2*sigma in gray
#black solid: patient form our series
#gray dots: patients from the literatures
#blue shaded dots simulated patients

### Model for Mitotic count

#### Prior predictive simulation for  mitotoic count
new_title <- c("ZIPoisson(0.5)", "ZIPoisson(1)", "ZIPoisson(1.5)")
par(mfrow=c(1,3) ) # 3 plots in 1 row 
for ( d in 1:3) 
{ lambda <- c(0.5,1,1.5) 
p <- runif( 1e4)
prior_m <- rzipois( 1e4 , p, lambda[d] )
dens( prior_m, xlim= c(0,10), ylim= c(0,8), main = new_title[d], xlab = "Mitosis/10HPF")
}

#### Probabilistic Model
i <- !is.na(Table_matrix$Mitotic_count_10HPF)
d <- list(M = as.integer(Table_matrix$Mitotic_count_10HPF[i]))
m1mitotic <-ulam(
  alist(
    M ~ dzipois( p, lambda),
    logit(p) <- ap,
    log(lambda) <- al,
    ap ~ dnorm(-1.5, 1), #more often than not we see mitosis
    al ~ dnorm (0.5, 0.5) #tumor have a low mitotic activity
  ), data = d, chains = 4)

#### Model ceofficients
precis(m1mitotic,2)
inv_logit(precis(m1mitotic, pars = "ap")) #. probabiligy of seeng mitosis
exp(precis(m1mitotic, pars = "al")) #true mitotic rate when all the mitosis are present on the slide 

#### Posterior predicitive simulations
post <- extract.samples( m1mitotic , n=1e4 )
t <- 1:1015
p <- c(d$M,rzipois(1000, inv_logit(post$ap), exp(post$al)))
idx <- c(rep(1,4),rep(2,11),rep(0,1000))
dat <- data.frame(p = p, idx = as.factor(idx))
t <- sample(t,1015)

#### Plotting the posterior predictive check comparing predicted data with actual data 
par(mfrow = c (1,1))
plot(dat$p[t], col = c(col.alpha(rangi2,0.5) ,"black" ,"black")[dat$idx[t]],
     pch = c(1,19,21)[dat$idx[t]],  bg = "lightgray", cex = 1.2,
     xaxt="n", main= "Mitotic count of Myoid tumors",
     ylab = "Number of mitosis/10HPF", xlab = "")
abline(h=mean(p), col=c("blue"))
#black solid: patient form our series
#gray dots: patients from the literatures
#blue shaded dots simulated patients

### Model for epithelial differentiation

#### Data wrangling
TM2 <- Table_matrix %>% gather(key = "feature", value = "score", 
                               "S100":"desmosomes")
TM2$score<- as.integer(TM2$score)
epit <- TM2 %>% filter(feature %in% c("CK5.6","CK8.18","CKAE1.AE3","desmosomes"))
epit$score <- ifelse(epit$score == 0, 0L, 1L)
feature.lelvels <- 1L:4L
epit$feature.new <- feature.lelvels[as.factor(epit$feature)]
i <- !is.na(epit$score)
dat <- list(
  S = epit$score[i],
  feature = epit$feature.new[i])

#### Prior predictive simulation for epithelial differentiation
new_title <- c("p = Normal(0, 0.15)", "p = Normal(0, 1.5)", "p = Normal(0, 10)")
adding <- c(FALSE, TRUE, TRUE)
par(mfrow=c(1,3) ) # 3 plots in 1 row 
for ( d in 1:3) 
{ sigma <- c(0.15,1.5 ,15 ) 
alpha <- rnorm(1e4, mean = 0, sd = sigma[d])
p <- inv_logit(alpha)
dens( p, main = new_title[d],xlim=c(0,1), xlab = "Probability of Epithelial Differentiation")}

#### Probabilistic Model
      # the following models are equivalent; the first (the centered version)
      # will result in serveral divergent transitions
mepit<- ulam( alist(         
  S ~ dbinom( 1 , p ),
  logit(p) <- a[feature] ,
  a[feature] ~ dnorm( a_bar , sigma ),
  a_bar ~ dnorm( 0 , 1.5 ) ,
  sigma ~ dexp( 1 )
), data=dat, chains=4 , cores=4 )

      # the non centered version you will have a smoother sampling see https://youtu.be/ZG3Oe35R5sY
i <- !is.na(epit$score)
dat <- list(
  S = epit$score[i],
  feature = epit$feature.new[i])

mepit<- ulam( alist(
  S ~ dbinom( 1 , p ),
  logit(p) <- a_bar + z[feature]*sigma ,
  z[feature] ~ dnorm( 0,1 ),
  a_bar ~ dnorm( 0 , 1.5 ) ,
  sigma ~ dexp( 1 )
), data=dat, chains=4 , cores=4 )

#### Model ceofficients
precis(mepit,2)

#### Posterior predicitive simulations
post <- extract.samples(mepit)
p_epit_bar <- inv_logit( post$a_bar )
p_epit <- inv_logit( post$z )
plot( precis( as.data.frame(p_epit) , xlim=c(0,1) ))
plot(precis(p_epit_bar))

#### Comparing probability densities 
par(mfrow=c(1,1))
dens(p_epit_bar, xlim= c(0,1), xlab ="Probability", main = "Epithelial Differentiation of Myoid tumors")
for(q in 1:4) dens(p_epit[,q], col=c("blue","red", "green", "purple")[q],add =TRUE)
legend("topleft", legend=c("Overall","CK5/6","CK8/18","CKAE1/AE3",
                           "Desmosomes"),  col=c("black","blue","red",
                                                 "green", "purple"), lty = 1, cex = 0.75)
