set_ulam_cmdstan(TRUE) #set cmdstan to true

m0 <- ulam(
  alist(
    nutcrackin ~ binomial(1,p),
    logit(p) <- ap ,
    ap ~ dnorm( 0 , 1 )
  ) , data=data_list , chains=4)

# to go from probability to log-odds use logit link 
logit(.99)
# to go from log-odds to proability/real scale use
logistic(.99)
inv_logit(.99)

post <- extract.samples(m0)
str(post)
dens(post$ap) #log odds scale
dens(logistic(post$ap)) #log odds scale
abline(v=mean(logistic(post$ap)))
precis( m0 ) # real probability scale
inv_logit(-1.27)
mean(focals_data_cropped$nutcrackin)


m1 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-a_id[id],
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)

plot(precis( m1 , depth=2 ))


m2 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-a_id[id],
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)

plot(precis( m1 , depth=2 ))

m3 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bPT*count_palm_std,
    # fixed priors
    c(ap,bPT) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)

plot(precis( m3 , depth=2 ))

m4 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bS*log_avg_stone_std ,
    # fixed priors
    c(ap,bS) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)
plot(precis( m4 , depth=2 ))
precis( m4 , depth=2 )

m5 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bPT*count_palm_std + bS*log_avg_stone_std ,
    # fixed priors
    c(ap,bS,bPT) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)
plot(precis( m5 , depth=2 ))
precis( m5 , depth=2 )

m6 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bPT*count_palm_std + bS*log_avg_stone_std 
    + bPTxS*count_palm_std*avg_stone_std,
    # fixed priors
    c(ap,bS,bPT,bPTxS) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)
plot(precis( m6 , depth=2 ))
precis( m6 , depth=2 )
compare(m1,m2,m3,m4,m5,m6)

######### varying effects of slopes per individual
set.seed(384)
m7 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bPT + a_id[id,2])*count_palm_std,
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bPT) ~ dnorm( 0 , 1 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16 , log_lik=FALSE, pars_omit = 'z_id')

precis(m7 , depth=3)

set.seed(482)
m8 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*log_avg_stone_std ,
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS) ~ dnorm( 0 , 1 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16, log_lik=FALSE, pars_omit = 'z_id')

precis(m8 , depth=3)

set.seed(183)
m9 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*log_avg_stone_std +  (bPT + a_id[id,3])*count_palm_std,
    # adaptive priors - non-centered
    transpars> matrix[id,3]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[3,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,bPT) ~ dnorm( 0 , 1 ),
    vector[3]:sigma_id ~ dexp(1),
    cholesky_factor_corr[3]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3,3]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16 , log_lik=FALSE, pars_omit = 'z_id')

precis(m9 , depth=3)

set.seed(345)
m10 <- ulam(
  alist(
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*log_avg_stone_std +
      (bPT + a_id[id,3])*count_palm_std + (bPTxS + a_id[id,4])*log_avg_stone_std*count_palm_std,
    # adaptive priors - non-centered
    transpars> matrix[id,4]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[4,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,bPT,bPTxS) ~ dnorm( 0 , 1 ),
    vector[4]:sigma_id ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16 , log_lik=FALSE, pars_omit = 'z_id')


precis(m10 , depth=3)
plot(precis(m10 , depth=2))
precis(m10 , depth=3 , pars="Rho_id")
plot(precis(m10 , depth=3 , pars="a_id"))

#compare(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10)

##plot preds stones
a_id_z <- matrix(0,2000,length(unique(data_list_daily$id))) #need to add zeros in VE to plot main effect
stone_seq <-  seq(from=min(data_list_daily$log_avg_stone_std) , 
                  to=max(data_list_daily$log_avg_stone_std) , length=30)
dpred <- list(
  id=rep(1,30), #list of ones we will replace
  log_avg_stone_std=stone_seq,
  count_palm_std=rep(0,30)
)

link2 <- link(m10, data=dpred , replace=list(id=a_id_z) )

str(link2)
plot(data_list_daily$nutcrackin ~ data_list_daily$log_avg_stone_std , pch="|" , 
     col=col.alpha(rangi2,0.01) , ylab="probability of nut cracking" , 
     xlab= "log stone density (standardized)")
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
  lines( link2[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue", alpha=0.1) , lty=1)
}
#axis( 1 , at= ( seq(from=0 , to=0.45 , by=0.05) - mean(dc$c70))/sd(dc$c70) , labels= seq(from=0 , to=0.45 , by=0.05) )

#plot daily id preds
plot(data_list_daily$nutcrackin ~ data_list_daily$log_avg_stone_std , pch="|" , 
     col=col.alpha(rangi2,0.01) , ylab="probability of nut cracking" , 
     xlab= "log stone density (standardized)")
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for ( i in 1:max(data_list_daily$id)){
  dpred <- list(
    id=rep(i,30), #list of ones we will replace
    log_avg_stone_std=stone_seq,
    count_palm_std=rep(0,30)
  )
  link2 <- link(m10, data=dpred )
  pred_median <- apply(link2 , 2 , median)
  lines(pred_median ~ stone_seq , lw=2, col=col.alpha("slateblue",0.5) , lty=1)
}
####nut cracking and palm trees

plot(data_list_daily$nutcrackin ~ data_list_daily$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
     xlab= "palm tree density/ 22 m^2 grid (standardized)")

palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                 to=max(data_list_daily$count_palm_std) , length=30)
dpred <- list(
  id=rep(1,30), #list of ones we will replace
  log_avg_stone_std=rep(0,30),
  count_palm_std=palm_seq
)

link2 <- link(m10, data=dpred , replace=list(id=a_id_z) )
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
  lines( link2[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
}


## per id plot
plot(data_list_daily$nutcrackin ~ data_list_daily$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
     xlab= "palm tree density/ 22 m^2 grid (standardized)")

lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for ( i in 1:max(data_list_daily$id)){
  dpred <- list(
    id=rep(i,30), #list of ones we will replace
    log_avg_stone_std=rep(0,30),
    count_palm_std=palm_seq
  )
  link2 <- link(m10, data=dpred )
  pred_median <- apply(link2 , 2 , median)
  lines(pred_median ~ palm_seq , lw=2, col=col.alpha("green4",0.5) , lty=1)
}





post <- extract.samples(m10)
str(post)


##### lets vctorize ST abundance



m8me <- ulam(
  alist(
    #stone model
    stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
    log(lambda) <-  a_bar + a_grid[stone_grid_index],
    #priors
    a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
    a_grid[stone_grid_index] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
    sigma_stone ~ dexp(1), #p
   #nutcracking model
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*a_grid[grid_id_follow] , #use the posterior of corresponding grid as predictor
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS) ~ dnorm( 0 , 1 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=4, log_lik=FALSE, pars_omit = 'z_id' , iter=4000 , control=list(adapt_delta=0.99))

precis(m8me , depth=3)
plot(precis(m8me , depth=2))

set.seed(183)
m9me <- ulam(
  alist(
    #stone model
    stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
    log(lambda) <-  a_bar + a_grid[stone_grid_index],
    #priors
    a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
    a_grid[stone_grid_index] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
    sigma_stone ~ dexp(1), #p
    #nutcracking model
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*a_grid[grid_id_follow] +  (bPT + a_id[id,3])*count_palm_std,
    # adaptive priors - non-centered
    transpars> matrix[id,3]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[3,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,bPT) ~ dnorm( 0 , 1 ),
    vector[3]:sigma_id ~ dexp(1),
    cholesky_factor_corr[3]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3,3]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=4, log_lik=FALSE, iter=4000,  pars_omit = 'z_id' , control=list(adapt_delta=0.99))

precis(m9me , depth=3)
plot(precis(m9me , depth=2))

set.seed(345)
m10me <- ulam(
  alist(
    #stone model
    stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
    log(lambda) <-  a_bar + a_grid[stone_grid_index],
    #priors
    a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
    a_grid[stone_grid_index] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
    sigma_stone ~ dexp(1), #p
    #nutcracking model
    nutcrackin ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*a_grid[grid_id_follow] +
      (bPT + a_id[id,3])*count_palm_std + (bPTxS + a_id[id,4])*a_grid[grid_id_follow]*count_palm_std,
    # adaptive priors - non-centered
    transpars> matrix[id,4]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[4,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,bPT,bPTxS) ~ dnorm( 0 , 1 ),
    vector[4]:sigma_id ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=4, log_lik=FALSE, iter=4000, pars_omit = 'z_id' , control=list(adapt_delta=0.99))


precis(m10me , depth=3)
plot(precis(m10me , depth=2))
precis(m10me , depth=3 , pars="Rho_id")
plot(precis(m10me , depth=3 , pars="a_id"))
