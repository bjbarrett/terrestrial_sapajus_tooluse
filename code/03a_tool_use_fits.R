set_ulam_cmdstan(TRUE) #set cmdstan to true

#below are fits looking at pure probability of tool use in a follow
#most basic logistic model estimating mean

m1 <- ulam(
  alist(
    nutcracking ~ binomial(1,p),
    logit(p) <- ap ,
    ap ~ dnorm( 0 , 1 )
  ) , data=data_list , chains=4)

# to go from probability to log-odds use logit link 
logit(.99)
# to go from log-odds to proability/real scale use
logistic(.99)
inv_logit(.99)

post <- extract.samples(m1)
str(post)
dens(post$ap) #log odds scale
dens(logistic(post$ap)) #log odds scale
abline(v=mean(logistic(post$ap)))
precis( m1 ) # real probability scale
inv_logit(-1.27)
mean(focals_data_cropped$nutcrackin)

# varying intercepts per individual
m2 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-a_id[id],
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)

plot(precis( m1 , depth=2 ))

#varying intercepts per individual and fixed effects palm
m3 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bPT*count_palm_std,
    # fixed priors
    c(ap,bPT) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)

plot(precis( m3 , depth=2 ))

#varying intercepts per individual and fixed effects stone
m4 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bS*log_avg_stone_std ,
    # fixed priors
    c(ap,bS) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)
plot(precis( m4 , depth=2 ))
precis( m4 , depth=2 )

#varying intercepts per individual and fixed effects stone and palm

m5 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bPT*count_palm_std + bS*log_avg_stone_std ,
    # fixed priors
    c(ap,bS,bPT) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)
plot(precis( m5 , depth=2 ))
precis( m5 , depth=2 )

#varying intercepts per individual and fixed effects stone and palm plus interaction

m6 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-a_id[id] + bPT*count_palm_std + bS*log_avg_stone_std 
    + bPTxS*count_palm_std*avg_stone_std,
    # fixed priors
    c(ap,bS,bPT,bPTxS) ~ dnorm( 0 , 1 ),
    a_id[id] ~ dnorm(ap,sigma_id),
    sigma_id ~ dexp(1)
    
  ) , data=data_list_daily , chains=4 , cores=4 , log_lik=TRUE)
plot(precis( m6 , depth=2 ))
precis( m6 , depth=2 )
compare(m2,m3,m4,m5,m6)

######### varying effects of slopes per individual
set.seed(384)
m7 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
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
    nutcracking ~ dbinom(1, p ),
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
    nutcracking ~ dbinom(1, p ),
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
    nutcracking ~ dbinom(1, p ),
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
plot(data_list_daily$nutcracking ~ data_list_daily$log_avg_stone_std , pch="|" , 
     col=col.alpha(rangi2,0.01) , ylab="probability of nut cracking" , 
     xlab= "log stone density (standardized)")
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
  lines( link2[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue", alpha=0.1) , lty=1)
}
#axis( 1 , at= ( seq(from=0 , to=0.45 , by=0.05) - mean(dc$c70))/sd(dc$c70) , labels= seq(from=0 , to=0.45 , by=0.05) )

#plot daily id preds
plot(data_list_daily$nutcracking ~ data_list_daily$log_avg_stone_std , pch="|" , 
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

plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
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
plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
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

##stone abundance
par( mfrow=c(1,3))
stone_samp <- c(-2,0,2)

for(j in 1:3){
  plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
       pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
       xlab= "palm tree density/ 22 m^2 grid (standardized)" , main=stone_samp[j] )
  
  palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                   to=max(data_list_daily$count_palm_std) , length=30)
  dpred <- list(
    id=rep(1,30), #list of ones we will replace
    log_avg_stone_std=rep(stone_samp[j],30),
    count_palm_std=palm_seq
  )
  
  link2 <- link(m10, data=dpred , replace=list(id=a_id_z) )
  pred_median <- apply(link2 , 2 , median)
  lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)
  
  for (j in sample( c(1:2000) , 100) ){
    lines( link2[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
  }
}
#dev.off()



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
    nutcracking ~ dbinom(1, p ),
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

m8gp <- ulam(
        alist(
            #stone model
            stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
            log(lambda) <-  a_bar + a_grid[stone_grid_index],
            #priors
            a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
            etasq ~ dexp( 2 ),
            rhosq ~ dexp( 0.5 ),
            # non-centered Gaussian Process prior
            transpars> vector[40]: a_grid <<- L_SIGMA * z,
            vector[40]: z ~ normal( 0 , 1 ),
            transpars> matrix[40,40]: L_SIGMA <<- cholesky_decompose( SIGMA ),
            transpars> matrix[40,40]: SIGMA <- cov_GPL2( Dmat110 , etasq , rhosq , 0.01 ),
            #nutcracking model
            nutcracking ~ dbinom(1, p ),
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
            
        ) , data=data_list_daily_2 , chains=4 , cores=4, log_lik=FALSE, pars_omit = 'z_id' , iter=4000 , control=list(adapt_delta=0.99))

plot(precis(m8gp, depth=2 , pars='a_grid'))
plot(precis(m8me, depth=2 , pars='a_grid'))
precis(m8gp)
precis(m8me)

##based on this, covariance of terrestriality rates drops off around 400m
#prior <- extract.prior(mzag_12)
post <- extract.samples(m8gp)
# plot the posterior median covariance function for p
plot( NULL , xlab="distance (hectometers)" , ylab="covariance" ,
      xlim=c(0,20) , ylim=c(0,3) )
# compute posterior mean covariance
x_seq <- seq( from=0 , to=10 , length.out=100 )
pmcov <- sapply( x_seq , function(x) post$etasq*exp(-post$rhosq*x^2) )
pmcov_mu <- apply( pmcov , 2 , mean )
lines( x_seq , pmcov_mu , lwd=2 )
# plot 60 functions sampled from posterior
for ( i in 1:50 )
    curve( post$etasq[i]*exp(-post$rhosq[i]*x^2) , add=TRUE ,
           col=col.alpha("black",0.14) )



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
    nutcracking ~ dbinom(1, p ),
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
a_grid_f_z <- matrix(0,2000,length(unique(data_list_daily$grid_id_follow))) #need to add zeros in VE to plot main effect
a_grid_s_z <- matrix(0,2000,length(unique(data_list_daily$stone_grid_index))) #need to add zeros in VE to plot main effect

post <- extract.samples(m9me)
median(post$a_bar)
str(post)
png(file = "plots/p_nc_m9me.png") 

    plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
         pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
         xlab= "palm tree density/ 22 m^2 grid (standardized)")
    
    palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                     to=max(data_list_daily$count_palm_std) , length=30)
    dpred <- list(
        id=rep(1,30), #list of ones we will replace
        count_palm_std=palm_seq,
        stone_grid_index=rep(1,30),
        grid_id_follow=rep(1,30)
    )
    
    link2 <- link(m9me, data=dpred , replace=list(id=a_id_z , stone_grid_index=a_grid_s_z , grid_id_follow=a_grid_f_z) )
    pred_median <- apply(link2$p , 2 , median)
    lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)
    
    for (j in sample( c(1:2000) , 100) ){
        lines( link2$p[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
    }
dev.off()

## per id plot
png(file = "plots/p_nc_all_inds_m9me.png") 
    
    plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
         pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
         xlab= "palm tree density/ 22 m^2 grid (standardized)")
    
    lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)
    
    for ( i in 1:max(data_list_daily$id)){
        dpred <- list(
            id=rep(i,30), #list of ones we will replace
            log_avg_stone_std=rep(median(post$a_bar),30),
            count_palm_std=palm_seq,
            stone_grid_index=rep(1,30),
            grid_id_follow=rep(1,30)
        )
        link2 <- link(m9me, data=dpred , replace=list(stone_grid_index=a_grid_s_z , grid_id_follow=a_grid_f_z) )
        pred_median <- apply(link2$p , 2 , median)
        lines(pred_median ~ palm_seq , lw=2, col=col.alpha("green4",0.5) , lty=1)
    }
dev.off()


png(file = "plots/stone_m9me.png") 
    p_link_stone <- function(x){
        id <- 1
        logit_scale <- with(post , ap + a_id[,id,1]*0 +  (bS + a_id[,id,2])*x +
                                (bPT + a_id[,id,3]*0)*0 )
        return(logistic(logit_scale))
    }
    
    stone_seq <-  seq(from=min(data_list_daily$log_avg_stone_std) , 
                      to=max(data_list_daily$log_avg_stone_std) , length=30)
    
    p_raw <- sapply( stone_seq, function(x) p_link_stone(x) )
    
    plot(data_list_daily$nutcracking ~ data_list_daily$log_avg_stone_std ,
         pch=19 , col=col.alpha("darkslateblue",0.004), ylab="probability of nut cracking",
         xlab= "log stone density/ 110 m^2 grid (standardized)")
    
    pred_median <- apply(p_raw , 2 , median)
    lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)
    
    for (j in sample( c(1:2000) , 100) ){
        lines( p_raw[j,] ~ stone_seq , lw=3, col=col.alpha("darkslateblue", alpha=0.1) , lty=1)
    }
dev.off()

m_tree_gp <- ulam(
    alist(
        #stone model
        count_palm_all ~ dpois( lambda_palm ), #define poisson likelihood
        log(lambda_palm) <-  a_bar_palm + a_cell[cell_id_all],
        #priors
        a_bar_palm ~ dnorm( 0 , 1.2 ), #mean effect on logscale
        etasq_palm ~ dexp( 1 ),
        rhosq_palm ~ dexp( 1 ),
        # non-centered Gaussian Process prior
        transpars> vector[1000]: a_cell <<- L_SIGMA_PALM * z_palm,
        vector[1000]: z_palm ~ normal( 0 , 1 ),
        transpars> matrix[1000,1000]: L_SIGMA_PALM <<- cholesky_decompose( SIGMA_PALM ),
        transpars> matrix[1000,1000]: SIGMA_PALM <- cov_GPL2( Dmat22 , etasq_palm , rhosq_palm , 0.01 )
) , data=data_list_daily_2 , chains=4 , cores=4, log_lik=FALSE, iter=1000,  pars_omit = 'z_id' , control=list(adapt_delta=0.9) , cmdstan=TRUE , threads=5 , refresh=10)
        
set.seed(183)
m9gp <- ulam(
    alist(
        #stone model
        stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
        log(lambda) <-  a_bar + a_grid[stone_grid_index],
        #priors
        a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 ),
        # non-centered Gaussian Process prior
        transpars> vector[40]: a_grid <<- L_SIGMA * z,
        vector[40]: z ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA <<- cholesky_decompose( SIGMA ),
        transpars> matrix[40,40]: SIGMA <- cov_GPL2( Dmat110 , etasq , rhosq , 0.01 ),
        #nutcracking model
        nutcracking ~ dbinom(1, p ),
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
        
    ) , data=data_list_daily_2 , chains=4 , cores=4, log_lik=FALSE, iter=4000,  pars_omit = 'z_id' , control=list(adapt_delta=0.99))

post <- extract.samples(m9gp)
p_link_stone <- function(x){
    id <- 1
    logit_scale <- with(post , ap + a_id[,id,1]*0 +  (bS + a_id[,id,2])*x +
                            (bPT + a_id[,id,3]*0)*0 )
    return(logistic(logit_scale))
}

stone_seq <-  seq(from=min(data_list_daily_2$log_avg_stone_std) , 
                  to=max(data_list_daily_2$log_avg_stone_std) , length=30)

p_raw <- sapply( stone_seq, function(x) p_link_stone(x) )

plot(data_list_daily_2$nutcracking ~ data_list_daily_2$log_avg_stone_std ,
     pch=19 , col=col.alpha("darkslateblue",0.004), ylab="probability of nut cracking",
     xlab= "log stone density/ 110 m^2 grid (standardized)")

pred_median <- apply(p_raw , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
    lines( p_raw[j,] ~ stone_seq , lw=3, col=col.alpha("darkslateblue", alpha=0.1) , lty=1)
}


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
    nutcracking ~ dbinom(1, p ),
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

plot(precis(m10me , depth=2))
plot(precis(m10me , depth=2))

precis(m10me , depth=3 , pars="Rho_id")
plot(precis(m10me , depth=3 , pars="a_id"))

#plot m10me
a_grid_f_z <- matrix(0,2000,length(unique(data_list_daily$grid_id_follow))) #need to add zeros in VE to plot main effect
a_grid_s_z <- matrix(0,2000,length(unique(data_list_daily$stone_grid_index))) #need to add zeros in VE to plot main effect

post <- extract.samples(m10me)
median(post$a_bar)
str(post)
plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
     xlab= "palm tree density/ 22 m^2 grid (standardized)")

palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                 to=max(data_list_daily$count_palm_std) , length=30)
dpred <- list(
  id=rep(1,30), #list of ones we will replace
  count_palm_std=palm_seq,
  stone_grid_index=rep(1,30),
  grid_id_follow=rep(1,30)
)

link2 <- link(m10me, data=dpred , replace=list(id=a_id_z , stone_grid_index=a_grid_s_z , grid_id_follow=a_grid_f_z) )
pred_median <- apply(link2$p , 2 , median)
lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
  lines( link2$p[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
}


## per id plot
plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
     xlab= "palm tree density/ 22 m^2 grid (standardized)")

lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for ( i in 1:max(data_list_daily$id)){
  dpred <- list(
    id=rep(i,30), #list of ones we will replace
    log_avg_stone_std=rep(median(post$a_bar),30),
    count_palm_std=palm_seq,
    stone_grid_index=rep(1,30),
    grid_id_follow=rep(1,30)
  )
  link2 <- link(m10me, data=dpred , replace=list(stone_grid_index=a_grid_s_z , grid_id_follow=a_grid_f_z) )
  pred_median <- apply(link2$p , 2 , median)
  lines(pred_median ~ palm_seq , lw=2, col=col.alpha("green4",0.5) , lty=1)
}

m11me <- ulam(
    alist(
        #stone model
        stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
        log(lambda) <-  a_bar + a_grid[stone_grid_index],
        #priors
        a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
        a_grid[stone_grid_index] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
        sigma_stone ~ dexp(1), #p
        #nutcracking model
        nutcracking ~ dbinom(1, p ),
        logit(p) <-ap + a_id[id,1] + a_g[grid_id_follow_dummy] + (bS + a_id[id,2])*a_grid[grid_id_follow] +
            (bPT + a_id[id,3])*count_palm_std + (bPTxS + a_id[id,4])*a_grid[grid_id_follow]*count_palm_std,
        # adaptive priors - non-centered
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        a_g[grid_id_follow_dummy] ~ dnorm( 0 , sigma_grid ), #priovarying effects centered on mean
        # fixed priors
        c(ap,bS,bPT,bPTxS) ~ dnorm( 0 , 1 ),
        vector[4]:sigma_id ~ dexp(1),
        sigma_grid ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
        
    ) , data=data_list_daily , chains=4 , cores=4, log_lik=FALSE, iter=4000, pars_omit = 'z_id' , control=list(adapt_delta=0.99))

m12me <- ulam(
    alist(
        #stone model
        stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
        log(lambda) <-  a_bar + a_grid[stone_grid_index],
        #priors
        a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
        a_grid[stone_grid_index] ~ dnorm( 0 , sigma_stone ), #priovarying effects centered on mean
        sigma_stone ~ dexp(1), #p
        #nutcracking model
        nutcracking ~ dbinom(1, p ),
        logit(p) <- a_id[id,1] + (bS + a_id[id,2])*a_grid[grid_id_follow] +
            (bPT + a_id[id,3])*count_palm_std + (bPTxS + a_id[id,4])*a_grid[grid_id_follow]*count_palm_std
        + bSEX[sex_index],
        # adaptive priors - non-centered
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        # fixed priors
        c(ap,bS,bPT,bPTxS) ~ dnorm( 0 , 1 ),
        bSEX[sex_index] ~ dnorm( ap , 1 ),
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
        
    ) , data=data_list_daily , chains=4 , cores=4, log_lik=FALSE, iter=4000, pars_omit = 'z_id' , control=list(adapt_delta=0.99))


precis(m12me , depth=2)
###add time

m6 <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
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
m7e <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bPT + a_id[id,2])*count_palm_std + b_exp*log_follow_length,
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bPT,b_exp) ~ dnorm( 0 , 1 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16 , log_lik=FALSE, pars_omit = 'z_id')

precis(m7e , depth=3)

set.seed(482)
m8e <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*log_avg_stone_std + b_exp*log_follow_length,
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,b_exp) ~ dnorm( 0 , 1 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16, log_lik=FALSE, pars_omit = 'z_id')

precis(m8 , depth=3)

set.seed(183)
m9e <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*log_avg_stone_std +  (bPT + a_id[id,3])*count_palm_std + b_exp*log_follow_length,
    # adaptive priors - non-centered
    transpars> matrix[id,3]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[3,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,bPT,b_exp) ~ dnorm( 0 , 1 ),
    vector[3]:sigma_id ~ dexp(1),
    cholesky_factor_corr[3]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3,3]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16 , log_lik=FALSE, pars_omit = 'z_id')

precis(m9e , depth=3)


set.seed(345)
m10e <- ulam(
  alist(
    nutcracking ~ dbinom(1, p ),
    logit(p) <-ap + a_id[id,1] +  (bS + a_id[id,2])*log_avg_stone_std +
      (bPT + a_id[id,3])*count_palm_std + (bPTxS + a_id[id,4])*log_avg_stone_std*count_palm_std
    + b_exp*log_follow_length,
    # adaptive priors - non-centered
    transpars> matrix[id,4]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[4,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bS,bPT,bPTxS,b_exp) ~ dnorm( 0 , 1 ),
    vector[4]:sigma_id ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    
  ) , data=data_list_daily , chains=4 , cores=16 , threads=16 , log_lik=FALSE, pars_omit = 'z_id')


precis(m10e , depth=3)
plot(precis(m10e , depth=2))
precis(m10e , depth=3 , pars="Rho_id")
plot(precis(m10e , depth=3 , pars="a_id"))


dpred <- list(
  id=rep(1,30), #list of ones we will replace
  log_avg_stone_std=stone_seq,
  count_palm_std=rep(0,30),
  log_follow_length=rep(log(5),30)
)

link2 <- link(m10e, data=dpred , replace=list(id=a_id_z) )

str(link2)
plot(data_list_daily$nutcracking ~ data_list_daily$log_avg_stone_std , pch="|" , 
     col=col.alpha(rangi2,0.01) , ylab="probability of nut cracking" , 
     xlab= "log stone density (standardized)")
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
  lines( link2[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue", alpha=0.1) , lty=1)
}
#axis( 1 , at= ( seq(from=0 , to=0.45 , by=0.05) - mean(dc$c70))/sd(dc$c70) , labels= seq(from=0 , to=0.45 , by=0.05) )

#plot daily id preds
plot(data_list_daily$nutcracking ~ data_list_daily$log_avg_stone_std , pch="|" , 
     col=col.alpha(rangi2,0.01) , ylab="probability of nut cracking" , 
     xlab= "log stone density (standardized)")
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for ( i in 1:max(data_list_daily$id)){
  dpred <- list(
    id=rep(i,30), #list of ones we will replace
    log_avg_stone_std=stone_seq,
    count_palm_std=rep(0,30),
    log_follow_length=rep(log(5),30)
  )
  link2 <- link(m10e, data=dpred )
  pred_median <- apply(link2 , 2 , median)
  lines(pred_median ~ stone_seq , lw=2, col=col.alpha("slateblue",0.5) , lty=1)
}

####nut cracking and palm trees

plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
     xlab= "palm tree density/ 22 m^2 grid (standardized)")

palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                 to=max(data_list_daily$count_palm_std) , length=30)
dpred <- list(
  id=rep(1,30), #list of ones we will replace
  log_avg_stone_std=rep(0,30),
  count_palm_std=palm_seq,
  log_follow_length=rep(log(5),30)
)

link2 <- link(m10e, data=dpred , replace=list(id=a_id_z) )
pred_median <- apply(link2 , 2 , median)
lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 100) ){
  lines( link2[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
}


## per id plot
plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
     xlab= "palm tree density/ 22 m^2 grid (standardized)")

lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for ( i in 1:max(data_list_daily$id)){
  dpred <- list(
    id=rep(i,30), #list of ones we will replace
    log_avg_stone_std=rep(0,30),
    count_palm_std=palm_seq,
    log_follow_length=rep(log(5),30)
  )
  link2 <- link(m10, data=dpred )
  pred_median <- apply(link2 , 2 , median)
  lines(pred_median ~ palm_seq , lw=2, col=col.alpha("green4",0.5) , lty=1)
}

##stone abundance
par( mfrow=c(1,3))
stone_samp <- c(-2,0,2)

for(j in 1:3){
  plot(data_list_daily$nutcracking ~ data_list_daily$count_palm_std ,
       pch=19 , col=col.alpha("green4",0.004), ylab="probability of nut cracking",
       xlab= "palm tree density/ 22 m^2 grid (standardized)" , main=stone_samp[j] )
  
  palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                   to=max(data_list_daily$count_palm_std) , length=30)
  dpred <- list(
    id=rep(1,30), #list of ones we will replace
    log_avg_stone_std=rep(stone_samp[j],30),
    count_palm_std=palm_seq,
    log_follow_length=rep(log(5),30)
  )
  
  link2 <- link(m10e, data=dpred , replace=list(id=a_id_z) )
  pred_median <- apply(link2 , 2 , median)
  lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)
  
  for (j in sample( c(1:2000) , 100) ){
    lines( link2[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
  }
}
#dev.off()



post <- extract.samples(m10)
str(post)
