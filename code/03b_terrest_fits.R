set_ulam_cmdstan(TRUE) #set cmdstan to true

#means only
m_h0 <- ulam(
  alist(
    height_dm ~ dzipois( p , lambda ),
    logit(p) <- ap,
    log(lambda) <- al,
    ap ~ dnorm( 0 , 1 ),
    al ~ dnorm( 3.5 , 2 )
  ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

precis( m_h0 )

hist(data_list_height$height_dm)

#vaying intercepts per individual of heights.
m_h1 <- ulam(
  alist(
    height_dm ~ dzipois( p , lambda ),
    logit(p) <- ap + a_id[id,1],
    log(lambda) <- al + a_id[id,2],
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    al ~ dnorm( 3.5 , 2 ),
    vector[2]:sigma_id ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
  ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

precis( m_h1 )
#ap is the probability of observing a zero (i.e. being terrestrial)
precis( m_h1 , depth=3)

#below has convergence issues, likely due to measuring scale ZAGamma might be better

# m_h1a <- ulam(
#   alist(
#     height_dm ~ dzipois( p , lambda ),
#     logit(p) <- ap + a_id[id,1] + a_f[follow_index,1],
#     log(lambda) <- al + a_id[id,2] + a_f[follow_index,2],
#     # adaptive priors - non-centered
#     transpars> matrix[id,2]:a_id <-
#       compose_noncentered( sigma_id , L_Rho_id , z_id ),
#     matrix[2,id]:z_id ~ normal( 0 , 1 ),
#     transpars> matrix[follow_index,2]:a_f <-
#       compose_noncentered( sigma_f , L_Rho_f , z_f ),
#     matrix[2,follow_index]:z_f ~ normal( 0 , 1 ),
#     # fixed priors
#     ap ~ dnorm( 0 , 1 ),
#     al ~ dnorm( 3.5 , 2 ),
#     vector[2]:sigma_id ~ dexp(1),
#     vector[2]:sigma_f ~ dexp(1),
#     cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
#     cholesky_factor_corr[2]:L_Rho_f ~ lkj_corr_cholesky( 3 ),
#     # compute ordinary correlation matrixes from Cholesky factors
#     gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id),
#     gq> matrix[2,2]:Rho_f <<- Chol_to_Corr(L_Rho_f)
#   ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

# precis(m_h1a)
# adding varying effects for day
m_h1a <- ulam(
  alist(
    height_dm ~ dzipois( p , lambda ),
    logit(p) <- ap + a_id[id,1] + a_d[day_index,1],
    log(lambda) <- al + a_id[id,2] + a_d[day_index,2],
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    transpars> matrix[day_index,2]:a_d <-
      compose_noncentered( sigma_d , L_Rho_d , z_d ),
    matrix[2,day_index]:z_d ~ normal( 0 , 1 ),
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    al ~ dnorm( 3.5 , 2 ),
    vector[2]:sigma_id ~ dexp(1),
    vector[2]:sigma_d ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id),
    gq> matrix[2,2]:Rho_d <<- Chol_to_Corr(L_Rho_d)
  ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

precis(m_h1a)
plot(precis(m_h1a , depth=3 , pars='a_d'))
precis(m_h1a, pars='Rho_d' , depth=3)
precis(m_h1a, pars='Rho_id' , depth=3)

#below does not work in ulam
# m_h1f <- ulam(
#   alist(
#     height ~ dzagamma2( p , mu , scale),
#     logit(p) <- ap + a_id[id,1],
#     log(mu) <- al + a_id[id,2],
#     # adaptive priors - non-centered
#     transpars> matrix[id,2]:a_id <-
#       compose_noncentered( sigma_id , L_Rho_id , z_id ),
#     matrix[2,id]:z_id ~ normal( 0 , 1 ),
#     # fixed priors
#     ap ~ dnorm( 0 , 1 ),
#     al ~ dnorm( 0 , 2 ),
#     vector[2]:sigma_id ~ dexp(1),
#     scale ~ dexp(1),
#     cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
#     # compute ordinary correlation matrixes from Cholesky factors
#     gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
#   ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

sort(unique(data_list_height$heights))

## add effects for nutcracking
m_h2 <- ulam(
  alist(
    height_dm ~ dzipois( p , lambda ),
    logit(p) <- ap + a_id[id,1] + (bNp + a_id[id,3])*nutcrackin ,
    log(lambda) <- al + a_id[id,2] + (bNl + a_id[id,4])*nutcrackin ,
    # adaptive priors - non-centered
    transpars> matrix[id,4]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[4,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
    al ~ dnorm( 3.5 , 2 ),
    vector[4]:sigma_id ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
  ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

precis( m_h2 , depth=3)
precis( m_h2 , depth=2)

post <- extract.samples(m_h2)

dens(logistic(post$ap) , xlim=c(0.05, 0.5 ) , xlab="probability of being terrestrial")
dens(logistic(post$ap + post$bNp) , add=TRUE , col="slateblue" )
legend(x = "topright", legend = c("tool use", "non tool use") , fill = c("slateblue","black") , bty='n' ) 
blog=1.5
for(i in 1:17){
  squanch <- logistic(post$ap + post$a_id[,i,1] )
  points( median(squanch) , (i-.1)/blog , pch=18)
  segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
  squanch <- logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2])
  points( median(squanch) , (i+.1)/blog , col="slateblue" , pch=16)
  segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i+.1)/blog , y1=(i+.1)/blog , col="slateblue" )
}


dens( (1-logistic(post$ap))*exp(post$al) , xlim=c(0, 200) , xlab="height in dm during follows" , show.HPDI=.99999 )
dens((1-logistic(post$ap + post$bNp))*exp(post$al + post$bNl) , add=TRUE , col="slateblue" , show.HPDI=.99999)
legend(x = "topright", legend = c("tool use", "non tool use") , fill = c("slateblue","black") , bty='n' ) 
blog=200
for(i in 1:17){
  squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3])
  points( median(squanch) , (i-.1)/blog , pch=18)
  segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
  squanch <- (1-logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2]))*exp(post$al + post$a_id[,i,3]+ post$bNl + post$a_id[,i,4])
  points( median(squanch) , (i+.1)/blog , col="slateblue" , pch=16)
  segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i+.1)/blog , y1=(i+.1)/blog , col="slateblue" )
}

for(i in 1:17){
  dens( (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3]) , xlim=c(0, 150) , add=TRUE)
}

for(i in 1:17){
  dens((1-logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2]))*exp(post$al + post$a_id[,i,3]+ post$bNl + post$a_id[,i,4]) , add=TRUE , col="slateblue" )
}

#if dot plot is desired
plot( 1,1, xlim=c(0,300) , ylim=c(0,18), col="white")
for(i in 1:17){
  squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3])
  points( median(squanch) , i-.1 , pch=18)
  segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=i-.1 , y1=i-.1  )
  squanch <- (1-logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2]))*exp(post$al + post$a_id[,i,3]+ post$bNl + post$a_id[,i,4])
  points( median(squanch) , i+.1 , col="slateblue" , pch=16)
  segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=i+.1 , y1=i+.1 , col="slateblue" )
}

## stop here and build this up by adding samples per follow as a predictor
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9042155/
# precis(m_h1a)
# adding varying effects for day
m_h1a <- ulam(
  alist(
    height_dm ~ dzipois( p , lambda ),
    logit(p) <- ap + a_id[id,1] + a_d[day_index,1],
    log(lambda) <- al + a_id[id,2] + a_d[day_index,2],
    # adaptive priors - non-centered
    transpars> matrix[id,2]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[2,id]:z_id ~ normal( 0 , 1 ),
    transpars> matrix[day_index,2]:a_d <-
      compose_noncentered( sigma_d , L_Rho_d , z_d ),
    matrix[2,day_index]:z_d ~ normal( 0 , 1 ),
    # fixed priors
    ap ~ dnorm( 0 , 1 ),
    al ~ dnorm( 3.5 , 2 ),
    vector[2]:sigma_id ~ dexp(1),
    vector[2]:sigma_d ~ dexp(1),
    cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    cholesky_factor_corr[2]:L_Rho_d ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id),
    gq> matrix[2,2]:Rho_d <<- Chol_to_Corr(L_Rho_d)
  ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

precis(m_h1a)
plot(precis(m_h1a , depth=3 , pars='a_d'))
precis(m_h1a, pars='Rho_d' , depth=3)
precis(m_h1a, pars='Rho_id' , depth=3)

#below does not work in ulam
# m_h1f <- ulam(
#   alist(
#     height ~ dzagamma2( p , mu , scale),
#     logit(p) <- ap + a_id[id,1],
#     log(mu) <- al + a_id[id,2],
#     # adaptive priors - non-centered
#     transpars> matrix[id,2]:a_id <-
#       compose_noncentered( sigma_id , L_Rho_id , z_id ),
#     matrix[2,id]:z_id ~ normal( 0 , 1 ),
#     # fixed priors
#     ap ~ dnorm( 0 , 1 ),
#     al ~ dnorm( 0 , 2 ),
#     vector[2]:sigma_id ~ dexp(1),
#     scale ~ dexp(1),
#     cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
#     # compute ordinary correlation matrixes from Cholesky factors
#     gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
#   ) , data=data_list_height , chains=4 , cores=4 , log_lik=TRUE)

sort(unique(data_list_height$heights))

## add effects for nutcracking
m_h2 <- ulam(
  alist(
    height_dm ~ dzipois( p , lambda ),
    logit(p) <- ap + a_id[id,1] + (bNp + a_id[id,3])*nutcrackin ,
    log(lambda) <- al + a_id[id,2] + (bNl + a_id[id,4])*nutcrackin ,
    # adaptive priors - non-centered
    transpars> matrix[id,4]:a_id <-
      compose_noncentered( sigma_id , L_Rho_id , z_id ),
    matrix[4,id]:z_id ~ normal( 0 , 1 ),
    # fixed priors
    c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
    al ~ dnorm( 3.5 , 2 ),
    vector[4]:sigma_id ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
  ) , data=data_list_height , chains=4 , cores=4)

precis( m_h2 , depth=3)
precis( m_h2 , depth=2)



##simulate zag
rzagamma2 <- function(n, p, mu, scale){
    ist_null <- rbinom(n, 1, prob=1-p)
    ist_null * rgamma2(n, mu=3, scale=2)
}

lala <- rzagamma2(n=500, p=0.3 , mu=3 , scale=2 )
dens(lala)
lala=list(y=lala)

mh2_alt <- ulam(
    alist(
        y|y==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        y|y>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( y | mu/scale , 1/scale) ) ,
        # y|y==0 ~ custom( log(p) ) ,
        # y|y>0 ~ custom( log1m(p) + gamma_lpdf( y | mu/scale , 1/scale) ) ,
        logit(p) <- ap,
        log(mu) <- al,
        ap ~ dnorm(0,1),
        al ~ dnorm( 0 , 2 ),
        scale~dexp(1)
    ) , 
    data=lala, chains=1   
    )

precis(mh2_alt)

### zag

mzag_0 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap,
        log(mu) <- al,
        ap ~ dnorm(0,1),
        al ~ dnorm( 1 , 2 ),
        scale~dexp(1)
    ) , 
    data=data_list_height, chains=4, cores=4   
)

precis(mzag_0 , depth=2)
logistic(-2.01)
exp(1.78)


mzag_1 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1],
        log(mu) <- al + a_id[id,2],
        # adaptive priors - non-centered
        transpars> matrix[id,2]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[2,id]:z_id ~ normal( 0 , 1 ),
        # fixed priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[2]:sigma_id ~ dexp(1),
        cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , 
    data=data_list_height, chains=4, cores=4   
)

precis(mzag_1)
mzag_2 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1] + (bNp + a_id[id,3])*nutcrackin ,
        log(mu) <- al + a_id[id,2] + (bNl + a_id[id,4])*nutcrackin ,
        # adaptive priors - non-centered
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        # fixed priors
        c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1) ,
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , 
    data=data_list_height, chains=4, cores=4   
)

trankplot(mzag_2, pars='a_id')
precis(mzag_2 , depth=3 , pars="Rho_id")
##pick up here we can add grid id
mzag_3 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1] + a_fz[follow_index] + (bNp + a_id[id,3])*nutcrackin ,
        log(mu) <- al + a_id[id,2] + a_fg[follow_index] + (bNl + a_id[id,4])*nutcrackin,
        # adaptive priors - non-centered
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        a_fz[follow_index] ~ dnorm(0,sigma_fz),
        a_fg[follow_index] ~ dnorm(0,sigma_fg),
        # fixed priors
        c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        c(scale,sigma_fz,sigma_fg) ~ dexp(1),
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        #cholesky_factor_corr[2]:L_Rho_f ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
        #gq> matrix[2,2]:Rho_f <<- Chol_to_Corr(L_Rho_f)
        
    ) , 
    data=data_list_height, chains=4, cores=4, control=list(adapt_delta=0.99)   
)
precis(mzag_3 , depth=1)

precis(mzag_3 , depth=3 , pars="a_fz")

## zero augmented gamma
post <- extract.samples(mzag_2)
str(post)

dens(logistic(post$ap) , xlim=c(0.05, 0.5 ) , xlab="probability of being terrestrial")
dens(logistic(post$ap + post$bNp) , add=TRUE , col="slateblue" )
legend(x = "topright", legend = c("tool use", "non tool use") , fill = c("slateblue","black") , bty='n' ) 
blog=1.5
for(i in 1:17){
    squanch <- logistic(post$ap + post$a_id[,i,1] )
    points( median(squanch) , (i-.1)/blog , pch=18)
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    squanch <- logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2])
    points( median(squanch) , (i+.1)/blog , col="slateblue" , pch=16)
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i+.1)/blog , y1=(i+.1)/blog , col="slateblue" )
}


dens( (1-logistic(post$ap))*exp(post$al) , xlim=c(0, 30) , xlab="height in m during follows" , show.HPDI=.99999 , ylim=c(0,3) )
dens((1-logistic(post$ap + post$bNp))*exp(post$al + post$bNl) , add=TRUE , col="slateblue" , show.HPDI=.99999)
legend(x = "topright", legend = c("tool use", "non tool use") , fill = c("slateblue","black") , bty='n' ) 
blog=10
for(i in 1:17){
    squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3])
    points( median(squanch) , (i-.1)/blog , pch=18)
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    squanch <- (1-logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2]))*exp(post$al + post$a_id[,i,3]+ post$bNl + post$a_id[,i,4])
    points( median(squanch) , (i+.1)/blog , col="slateblue" , pch=16)
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i+.1)/blog , y1=(i+.1)/blog , col="slateblue" )
}

for(i in 1:17){
    dens( (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3]) , xlim=c(0, 150) , add=TRUE)
}

for(i in 1:17){
    dens((1-logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2]))*exp(post$al + post$a_id[,i,3]+ post$bNl + post$a_id[,i,4]) , add=TRUE , col="slateblue" )
}

mzag_4 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_f[follow_index,1],
        log(mu) <- al + a_f[follow_index,2],
        # adaptive priors - non-centered
        transpars> matrix[follow_index,2]:a_f <-
            compose_noncentered( sigma_f , L_Rho_f , z_f ),
        matrix[2,follow_index]:z_f ~ normal( 0 , 1 ),
        # fixed priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        c(scale) ~ dexp(1),
        vector[2]:sigma_f ~ dexp(1),
        cholesky_factor_corr[2]:L_Rho_f ~ lkj_corr_cholesky( 3 ),
        #cholesky_factor_corr[2]:L_Rho_f ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[2,2]:Rho_f <<- Chol_to_Corr(L_Rho_f)
        #gq> matrix[2,2]:Rho_f <<- Chol_to_Corr(L_Rho_f)
        
    ) , 
    data=data_list_height, chains=4, cores=4, control=list(adapt_delta=0.99)   
)

mzag_5 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1] + a_g22[cell_id_index,1],
        log(mu) <- al + a_id[id,2] + a_g22[cell_id_index,2],
        # adaptive priors - non-centered
        transpars> matrix[cell_id_index,2]:a_g22 <-
            compose_noncentered( sigma_g22 , L_Rho_g22 , z_g22 ),
        matrix[2,cell_id_index]:z_g22 ~ normal( 0 , 1 ),
        transpars> matrix[id,2]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[2,cell_id_index]:z_g22 ~ normal( 0 , 1 ),
        matrix[2,id]:z_id ~ normal( 0 , 1 ),
        
        # fixed priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[2]:sigma_id ~ dexp(1),
        vector[2]:sigma_g22 ~ dexp(1),
        cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        cholesky_factor_corr[2]:L_Rho_g22 ~ lkj_corr_cholesky( 3 ),
        
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        gq> matrix[2,2]:Rho_g22 <<- Chol_to_Corr(L_Rho_g22)
        
    ) , 
    data=data_list_height2, chains=4, cores=4   
)


precis(mzag_5)

mzag_6 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1] + a_g22[cell_id_index,1] + (bNp + a_id[id,3] + a_g22[cell_id_index,3])*nutcrackin,
        log(mu) <- al + a_id[id,2] + a_g22[cell_id_index,2] + (bNl + a_id[id,4] + a_g22[cell_id_index,4])*nutcrackin,
        # adaptive priors - non-centered
        transpars> matrix[cell_id_index,4]:a_g22 <-
            compose_noncentered( sigma_g22 , L_Rho_g22 , z_g22 ),
        matrix[4,cell_id_index]:z_g22 ~ normal( 0 , 1 ),
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[4]:sigma_id ~ dexp(1),
        vector[4]:sigma_g22 ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        cholesky_factor_corr[4]:L_Rho_g22 ~ lkj_corr_cholesky( 3 ),
        
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        gq> matrix[4,4]:Rho_g22 <<- Chol_to_Corr(L_Rho_g22)
        
    ) , 
    data=data_list_height2, chains=4, cores=4 , iter=1500   
)


mzag_7 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1] + a_g22[cell_id_index,1] + a_g110[grid_id_follow_dummy,1] + (bNp + a_id[id,3] + a_g22[cell_id_index,3] + a_g110[grid_id_follow_dummy,3])*nutcrackin,
        log(mu) <- al + a_id[id,2] + a_g22[cell_id_index,2] + a_g110[grid_id_follow_dummy,2]+ (bNl + a_id[id,4] + a_g22[cell_id_index,4] + a_g110[grid_id_follow_dummy,4])*nutcrackin,
        # adaptive priors - non-centered
        transpars> matrix[cell_id_index,4]:a_g22 <-
            compose_noncentered( sigma_g22 , L_Rho_g22 , z_g22 ),
        matrix[4,cell_id_index]:z_g22 ~ normal( 0 , 1 ),
        
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),

        transpars> matrix[grid_id_follow_dummy,4]:a_g110 <-
            compose_noncentered( sigma_g110 , L_Rho_g110 , z_g110 ),
        matrix[4,grid_id_follow_dummy]:z_g110 ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[4]:sigma_id ~ dexp(1),
        vector[4]:sigma_g22 ~ dexp(1),
        vector[4]:sigma_g110 ~ dexp(1),
        
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        cholesky_factor_corr[4]:L_Rho_g22 ~ lkj_corr_cholesky( 3 ),
        cholesky_factor_corr[4]:L_Rho_g110 ~ lkj_corr_cholesky( 3 ),
        
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        gq> matrix[4,4]:Rho_g22 <<- Chol_to_Corr(L_Rho_g22),
        gq> matrix[4,4]:Rho_g110 <<- Chol_to_Corr(L_Rho_g110)
        
    ) , 
    data=data_list_height2, chains=4, cores=4 , iter=1500   
)

precis(mzag_7)
mzag_8 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1] + a_g22[cell_id_index,1] + a_g110[grid_id_follow_dummy,1] 
            + (bNp + a_id[id,3] + a_g22[cell_id_index,3] + a_g110[grid_id_follow_dummy,3])*nutcrackin
            + (bSp + a_id[id,5])*log_avg_stone_std +  (bPTp + a_id[id,7])*count_palm_std,
        log(mu) <- al + a_id[id,2] + a_g22[cell_id_index,2] + a_g110[grid_id_follow_dummy,2]
            + (bNl + a_id[id,4] + a_g22[cell_id_index,4] + a_g110[grid_id_follow_dummy,4])*nutcrackin
            + (bSl + a_id[id,6])*log_avg_stone_std +  (bPTl + a_id[id,8])*count_palm_std,
        # adaptive priors - non-centered
        transpars> matrix[cell_id_index,4]:a_g22 <-
            compose_noncentered( sigma_g22 , L_Rho_g22 , z_g22 ),
        matrix[4,cell_id_index]:z_g22 ~ normal( 0 , 1 ),
        
        transpars> matrix[id,8]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[8,id]:z_id ~ normal( 0 , 1 ),
        
        transpars> matrix[grid_id_follow_dummy,4]:a_g110 <-
            compose_noncentered( sigma_g110 , L_Rho_g110 , z_g110 ),
        matrix[4,grid_id_follow_dummy]:z_g110 ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bNp,bNl,bSp,bSl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[8]:sigma_id ~ dexp(1),
        vector[4]:sigma_g22 ~ dexp(1),
        vector[4]:sigma_g110 ~ dexp(1),
        
        cholesky_factor_corr[8]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        cholesky_factor_corr[4]:L_Rho_g22 ~ lkj_corr_cholesky( 3 ),
        cholesky_factor_corr[4]:L_Rho_g110 ~ lkj_corr_cholesky( 3 ),
        
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[8,8]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        gq> matrix[4,4]:Rho_g22 <<- Chol_to_Corr(L_Rho_g22),
        gq> matrix[4,4]:Rho_g110 <<- Chol_to_Corr(L_Rho_g110)
        
    ) , 
    data=data_list_height2, chains=4, cores=4 , iter=1500   
)

precis(mzag_8, depth=1)
precis(mzag_8, depth=3 , pars='sigma_id')
precis(mzag_8, depth=3 , pars='sigma_g110')
precis(mzag_8, depth=3 , pars='sigma_g22')
precis(mzag_8, depth=3 , pars='Rho_id')
precis(mzag_8, depth=3 , pars='Rho_g22')


post <- extract.samples(mzag_8)
str(post)
dens(rnorm(3000,0,2) , col="darkgreen" , lty=3 , ylim=c(0,15) , xlim=c(-3,3))
for(i in 4:9) dens(post[[i]] , add=TRUE)


# no spatial
mzag_9 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1]  + (bNp + a_id[id,3])*nutcrackin + (bSp + a_id[id,5])*log_avg_stone_std + (bPTp + a_id[id,7])*count_palm_std,
        log(mu) <- al + a_id[id,2] + (bNl + a_id[id,4])*nutcrackin + (bSl + a_id[id,6])*log_avg_stone_std + (bPTl + a_id[id,8])*count_palm_std,
        # adaptive priors - non-centered
      
        transpars> matrix[id,8]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[8,id]:z_id ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bNp,bNl,bSp,bSl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[8]:sigma_id ~ dexp(1),
        cholesky_factor_corr[8]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[8,8]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , 
    data=data_list_height2, chains=4, cores=4 , iter=1500   
)

post <- extract.samples(mzag_9)
precis(mzag_9 , depth=1)
str(post)
dens(rnorm(3000,0,2) , col="darkgreen" , lty=3 , ylim=c(0,15) , xlim=c(-3,3))
for(i in 2:6) dens(post[[i]] , add=TRUE)




mzag_10 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1]  + (bNp + a_id[id,3])*nutcrackin + (bSp + a_id[id,5])*log_avg_stone_std + (bPTp + a_id[id,7])*count_palm_std
        + bSxNCp*log_avg_stone_std*nutcrackin + bPTxNCp*count_palm_std*nutcrackin,
        log(mu) <- al + a_id[id,2] + (bNl + a_id[id,4])*nutcrackin + (bSl + a_id[id,6])*log_avg_stone_std + (bPTl + a_id[id,8])*count_palm_std
        + bSxNCl*log_avg_stone_std*nutcrackin + bPTxNCl*count_palm_std*nutcrackin,
        # adaptive priors - non-centered
        
        transpars> matrix[id,8]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[8,id]:z_id ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bNp,bNl,bSp,bSl,bPTp,bPTl,bSxNCp,bPTxNCp,bSxNCl,bPTxNCl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[8]:sigma_id ~ dexp(1),
        cholesky_factor_corr[8]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[8,8]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , 
    data=data_list_height2, chains=4, cores=4 , iter=1500   
)

mzag_11ntu <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1]  + (bSp + a_id[id,3])*log_avg_stone_std + (bPTp + a_id[id,5])*count_palm_std,
        log(mu) <- al + a_id[id,2]  + (bSl + a_id[id,4])*log_avg_stone_std + (bPTl + a_id[id,6])*count_palm_std,
        # adaptive priors - non-centered
        
        transpars> matrix[id,6]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[6,id]:z_id ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bSp,bSl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[6]:sigma_id ~ dexp(1),
        cholesky_factor_corr[6]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[6,6]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , 
    data=data_list_height_ntu, chains=4, cores=4 , iter=1500   
)

precis(mzag_11ntu)
mzag_11tu <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_id[id,1]  + (bSp + a_id[id,3])*log_avg_stone_std + (bPTp + a_id[id,5])*count_palm_std,
        log(mu) <- al + a_id[id,2]  + (bSl + a_id[id,4])*log_avg_stone_std + (bPTl + a_id[id,6])*count_palm_std,
        # adaptive priors - non-centered
        
        transpars> matrix[id,6]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[6,id]:z_id ~ normal( 0 , 1 ),
        
        # fixed priors
        c(ap,bSp,bSl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        vector[6]:sigma_id ~ dexp(1),
        cholesky_factor_corr[6]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[6,6]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , 
    data=data_list_height_tu, chains=4, cores=4 , iter=1500   
)

precis(mzag_11tu)
precis(mzag_11ntu)




mzag_12 <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) , #likelihood loop for heingt gamma
        logit(p) <- ap + a_gp[grid_id_follow_dummy], #bern
        log(mu) <- al + a_gl[grid_id_follow_dummy], #gamma
        # fixed priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        c(scale) ~ dexp(1),
        vector[39]: a_gp ~ multi_normal( 0 , SIGMA_p ),
        matrix[39,39]: SIGMA_p <- cov_GPL2( Dmat110 , etasq_p , rhosq_p , 0.01 ),
        etasq_p ~ exponential(1),
        rhosq_p ~ exponential(1),
        vector[39]: a_gl ~ multi_normal( 0 , SIGMA_l ),
        matrix[39,39]: SIGMA_l <- cov_GPL2( Dmat110 , etasq_l , rhosq_l , 0.01 ),
        etasq_l ~ exponential(1),
        rhosq_l ~ exponential(1)
    ) , 
    data=data_list_height2, chains=6, cores=6, iter=1500 , control=list(adapt_delta=0.99)   
)

mzag_12nc <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_gp[grid_id_follow_dummy],
        log(mu) <- al + a_gl[grid_id_follow_dummy],
        
        # priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        transpars> vector[39]: a_gp <<- L_SIGMA_p * z_gp,
        vector[39]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[39,39]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[39,39]: SIGMA_p <- cov_GPL2( Dmat110 , etasq_p , rhosq_p , 0.01 ),
        transpars> vector[39]: a_gl <<- L_SIGMA_l * z_gl,
        vector[39]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[39,39]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[39,39]: SIGMA_l <- cov_GPL2( Dmat110 , etasq_l , rhosq_l , 0.01 ),
        etasq_p ~ exponential(1),
        rhosq_p ~ exponential(1),
        etasq_l ~ exponential(1),
        rhosq_l ~ exponential(1)
    ) , 
    data=data_list_height2, chains=6, cores=6, iter=1000   
)


precis(mzag_12nc)
precis(mzag_12)

precis(mzag_12nc, depth=2, pars='a_gp')
##based on this, covariance of terrestriality rates drops off around 400m
#prior <- extract.prior(mzag_12)
post <- extract.samples(mzag_12)
# plot the posterior median covariance function for p
plot( NULL , xlab="distance (hectometers)" , ylab="covariance" ,
      xlim=c(0,20) , ylim=c(0,1) )
# compute posterior mean covariance
x_seq <- seq( from=0 , to=10 , length.out=100 )
pmcov <- sapply( x_seq , function(x) post$etasq_p*exp(-post$rhosq_p*x^2) )
pmcov_mu <- apply( pmcov , 2 , mean )
lines( x_seq , pmcov_mu , lwd=2 )
# plot 60 functions sampled from posterior
for ( i in 1:50 )
    curve( post$etasq_p[i]*exp(-post$rhosq_p[i]*x^2) , add=TRUE ,
           col=col.alpha("black",0.3) )


#no spatial covariance in height given they are not terrestrial
# plot the posterior median covariance function for lambda/gamma
plot( NULL , xlab="distance (hectometers)" , ylab="covariance" ,
      xlim=c(0,20) , ylim=c(0,1) )
# compute posterior mean covariance
x_seq <- seq( from=0 , to=10 , length.out=100 )
pmcov <- sapply( x_seq , function(x) post$etasq_l*exp(-post$rhosq_l*x^2) )
pmcov_mu <- apply( pmcov , 2 , mean )
lines( x_seq , pmcov_mu , lwd=2 )
# plot 60 functions sampled from posterior
for ( i in 1:50 )
    curve( post$etasq_l[i]*exp(-post$rhosq_l[i]*x^2) , add=TRUE ,
           col=col.alpha("black",0.3) )


# compute posterior median covariance among grid cells
K_p <- K_l <- matrix(0,nrow=39,ncol=39)
d <- data_list_height2
for ( i in 1:39 ){
    for ( j in 1:39 ){
        K_p[i,j] <- median(post$etasq_p) * exp( -median(post$rhosq_p) * d$Dmat110[i,j]^2 )
        K_l[i,j] <- median(post$etasq_l) * exp( -median(post$rhosq_l) * d$Dmat110[i,j]^2 )
    }
}
diag(K_p) <- median(post$etasq_p) + 0.01 # why is there 0.01 added?
diag(K_l) <- median(post$etasq_l) + 0.01
# ger corr matrix
Rho_p <- round(cov2cor(K_p) , 2)
Rho_l <- round(cov2cor(K_l) , 2)

mapview(mtx_distance_110)
# https://github.com/r-spatial/sf/issues/231
seal_coords <- do.call(rbind, st_geometry(mtx_distance_110)) %>% 
    as_tibble() %>% setNames(c("lon","lat"))
plot(mtx_distance_110)
mtx_distance_110$post_med  <- mtx_distance_110$post_med_p  <- mtx_distance_110$post_med_1<- NA
for (i in 1:39){
    mtx_distance_110$post_med_p[i] <- median(logistic(with(post , ap + a_gp[,i])))
    mtx_distance_110$post_med_l[i] <- median(exp(with(post , al + a_gl[,i])))
    mtx_distance_110$post_med[i] <-  median( with( post ,  exp(al + a_gl[,i])*(1- logistic(ap + a_gp[,i]) ) ) )
}

heatmap(Rho_p, Colv = NA, Rowv = NA, scale="column")
heatmap(Rho_l, Colv = NA, Rowv = NA, scale="column")

plot(seal_coords$lon,seal_coords$lat )

# overlay lines shaded by Rho
for( i in 1:39 )
    for ( j in 1:39 )
        if ( i < j )
            lines( c( seal_coords$lon[i],seal_coords$lon[j] ) , c(seal_coords$lat[i],seal_coords$lat[j]) , lwd=2 , col=col.alpha("black",Rho_p[i,j]^2) )


mtx_distance_110$post_med  <- mtx_distance_110$post_med_p  <- mtx_distance_110$post_med_1<- NA
for (i in 1:39){
    mtx_distance_110$post_med_p[i] <- median(logistic(with(post , ap + a_gp[,i])))
    mtx_distance_110$post_med_l[i] <- median(exp(with(post , al + a_gl[,i])))
    mtx_distance_110$post_med[i] <-  median( with( post ,  exp(al + a_gl[,i])*(1- logistic(ap + a_gp[,i]) ) ) )
}


plot(mtx_distance_110)

mzag_13nc <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_gp[grid_id_follow_dummy],
        log(mu) <- al + a_gl[grid_id_follow_dummy],
        # fixed priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        transpars> vector[574]: a_gp <<- L_SIGMA_p * z_gp,
        vector[574]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[574,574]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[574,574]: SIGMA_p <- cov_GPL2( Dmat110 , etasq_p , rhosq_p , 0.01 ),
        
        transpars> vector[574]: a_gl <<- L_SIGMA_l * z_gl,
        vector[574]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[574,574]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[574,574]: SIGMA_l <- cov_GPL2( Dmat110 , etasq_l , rhosq_l , 0.01 ),
        
        # vector[39]: a_gp ~ multi_normal( 0 , SIGMA_p ),
        # matrix[39,39]: SIGMA_p <- cov_GPL2( Dmat110 , etasq_p , rhosq_p , 0.01 ),
        etasq_p ~ exponential(1),
        rhosq_p ~ exponential(1),
        # vector[39]: a_gl ~ multi_normal( 0 , SIGMA_l ),
        # matrix[39,39]: SIGMA_l <- cov_GPL2( Dmat110 , etasq_l , rhosq_l , 0.01 ),
        etasq_l ~ exponential(1),
        rhosq_l ~ exponential(1)
    ) , 
    data=data_list_height3, chains=6, cores=6, iter=1200 , control=list(adapt_delta=0.9)   
)

### terrestriality predictied by s and palm
m1 <-ulam(
    alist(
        ### main model
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        logit(p) <- ap + a_gp[grid_id_follow],
        log(mu) <- al + a_gl[grid_id_follow],
        # priors
        ap ~ dnorm( 0 , 1 ),
        al ~ dnorm( 1 , 2 ),
        scale ~ dexp(1),
        transpars> vector[40]: a_gp <<- L_SIGMA_p * z_gp,
        vector[40]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[40,40]: SIGMA_p <- cov_GPL2( Dmat110_full , etasq_p , rhosq_p , 0.01 ),
        transpars> vector[40]: a_gl <<- L_SIGMA_l * z_gl,
        vector[40]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[40,40]: SIGMA_l <- cov_GPL2( Dmat110_full , etasq_l , rhosq_l , 0.01 ),
        c(etasq_p,rhosq_p,etasq_l,rhosq_l) ~ exponential(1)
        ) , 
    data=data_list_height2, chains=6, cores=6, iter=1000   
)

m1a <- ulam(
    alist(
        # model stone abundance as a GP
        #stone model
        stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
        log(lambda) <-  s_bar + s_grid[stone_grid_index],
        #priors
        s_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
        etasq_s ~ dexp( 2 ),
        rhosq_s ~ dexp( 0.5 ),
        # non-centered Gaussian Process prior
        transpars> vector[40]: s_grid <<- L_SIGMA_s * z_s,
        vector[40]: z_s ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_s <<- cholesky_decompose( SIGMA_s ),
        transpars> matrix[40,40]: SIGMA_s <- cov_GPL2( Dmat110_full , etasq_s , rhosq_s , 0.01 ),
        
        ### main model
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        # logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bSp + a_id[id,3])*log_avg_stone_std + (bPTp + a_id[id,5])*count_palm_std,
        # log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bSl + a_id[id,4])*log_avg_stone_std + (bPTl + a_id[id,6])*count_palm_std,
        logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bSp + a_id[id,3])*z_s[grid_id_follow],
        log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bSl + a_id[id,4])*z_s[grid_id_follow],
        
        # priors
        al ~ dnorm( 1 , 2 ),
        c(ap,bSp,bSl) ~ dnorm( 0 , 1 ),
        scale ~ dexp(1),
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        # bernouli component GP priors
        transpars> vector[40]: a_gp <<- L_SIGMA_p * z_gp,
        vector[40]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[40,40]: SIGMA_p <- cov_GPL2( Dmat110_full , etasq_p , rhosq_p , 0.01 ),
        # gamma component GP priors
        transpars> vector[40]: a_gl <<- L_SIGMA_l * z_gl,
        vector[40]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[40,40]: SIGMA_l <- cov_GPL2( Dmat110_full , etasq_l , rhosq_l , 0.01 ),
        c(etasq_p,rhosq_p,etasq_l,rhosq_l) ~ exponential(1)
    ) , 
    data=data_list_height2, chains=6, cores=6, iter=1000   
)

precis(m1a)
precis(m1a , depth=2 , pars="s_grid")


m1b <- ulam(
    alist(
        # model stone abundance as a GP
        #stone model
        stone_raw_grid ~ dpois( lambda ), #define poisson likelihood
        log(lambda) <-  s_bar + s_grid[stone_grid_index],
        #priors
        s_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
        etasq_s ~ dexp( 2 ),
        rhosq_s ~ dexp( 0.5 ),
        # non-centered Gaussian Process prior
        transpars> vector[40]: s_grid <<- L_SIGMA_s * z_s,
        vector[40]: z_s ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_s <<- cholesky_decompose( SIGMA_s ),
        transpars> matrix[40,40]: SIGMA_s <- cov_GPL2( Dmat110_full , etasq_s , rhosq_s , 0.01 ),
        
        ### main model
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        # logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bSp + a_id[id,3])*log_avg_stone_std + (bPTp + a_id[id,5])*count_palm_std,
        # log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bSl + a_id[id,4])*log_avg_stone_std + (bPTl + a_id[id,6])*count_palm_std,
        logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bSp + a_id[id,3])*z_s[grid_id_follow]+ (bPTp + a_id[id,5])*count_palm_std,
        log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bSl + a_id[id,4])*z_s[grid_id_follow]+ (bPTl + a_id[id,6])*count_palm_std,
        
        # priors
        al ~ dnorm( 1 , 2 ),
        c(ap,bSp,bSl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        scale ~ dexp(1),
        transpars> matrix[id,6]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[6,id]:z_id ~ normal( 0 , 1 ),
        vector[6]:sigma_id ~ dexp(1),
        cholesky_factor_corr[6]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[6,6]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        # bernouli component GP priors
        transpars> vector[40]: a_gp <<- L_SIGMA_p * z_gp,
        vector[40]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[40,40]: SIGMA_p <- cov_GPL2( Dmat110_full , etasq_p , rhosq_p , 0.01 ),
        # gamma component GP priors
        transpars> vector[40]: a_gl <<- L_SIGMA_l * z_gl,
        vector[40]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[40,40]: SIGMA_l <- cov_GPL2( Dmat110_full , etasq_l , rhosq_l , 0.01 ),
        c(etasq_p,rhosq_p,etasq_l,rhosq_l) ~ exponential(1)
    ) , 
    data=data_list_height2, chains=6, cores=6, iter=1000   
)

precis(m1b)
precis(m1b , depth=2 , pars="sigma_id")
precis(m1b , depth=3 , pars="z_s")

m1d <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,

        logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bNp + a_id[id,3])*nutcrackin,
        log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bNl + a_id[id,4])*nutcrackin,
        
        # priors
        al ~ dnorm( 1 , 2 ),
        c(ap,bNp,bNl) ~ dnorm( 0 , 1 ),
        scale ~ dexp(1),
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        # bernouli component GP priors
        transpars> vector[40]: a_gp <<- L_SIGMA_p * z_gp,
        vector[40]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[40,40]: SIGMA_p <- cov_GPL2( Dmat110_full , etasq_p , rhosq_p , 0.01 ),
        # gamma component GP priors
        transpars> vector[40]: a_gl <<- L_SIGMA_l * z_gl,
        vector[40]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[40,40]: SIGMA_l <- cov_GPL2( Dmat110_full , etasq_l , rhosq_l , 0.01 ),
        c(etasq_p,rhosq_p,etasq_l,rhosq_l) ~ exponential(1)
    ) , 
    data=data_list_height_zero, chains=6, cores=6, iter=1000 , control=list(adapt_delta=0.97)  
)

precis(m1d)

m1e <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        
        logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bNp + a_id[id,3])*nutcrackin + (bPTp + a_id[id,5])*count_palm_std,
        log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bNl + a_id[id,4])*nutcrackin + (bPTl + a_id[id,6])*count_palm_std,
        
        # priors
        al ~ dnorm( 1 , 2 ),
        c(ap,bNp,bNl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        scale ~ dexp(1),
        transpars> matrix[id,6]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[6,id]:z_id ~ normal( 0 , 1 ),
        vector[6]:sigma_id ~ dexp(1),
        cholesky_factor_corr[6]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[6,6]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        # bernouli component GP priors
        transpars> vector[40]: a_gp <<- L_SIGMA_p * z_gp,
        vector[40]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[40,40]: SIGMA_p <- cov_GPL2( Dmat110_full , etasq_p , rhosq_p , 0.01 ),
        # gamma component GP priors
        transpars> vector[40]: a_gl <<- L_SIGMA_l * z_gl,
        vector[40]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[40,40]: SIGMA_l <- cov_GPL2( Dmat110_full , etasq_l , rhosq_l , 0.01 ),
        c(etasq_p,rhosq_p,etasq_l,rhosq_l) ~ exponential(1)
    ) , 
    data=data_list_height_zero, chains=6, cores=6, iter=1000 , control=list(adapt_delta=0.97)  
)

precis(m1e)
m1f <- ulam(
    alist(
        height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) , #likelihood loop for terrestrial Bernoulli
        height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
        
        logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bNp + a_id[id,3])*nutcrackin 
        + (bPTp + a_id[id,5])*count_palm_std + (bNxPTp + a_id[id,7])*count_palm_std*nutcrackin,
        log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bNl + a_id[id,4])*nutcrackin 
        + (bPTl + a_id[id,6])*count_palm_std +  (bNxPTl + a_id[id,8])*count_palm_std*nutcrackin,
        
        # priors
        al ~ dnorm( 1 , 2 ),
        c(ap,bNp,bNl,bPTp,bPTl,bNxPTp,bNxPTl) ~ dnorm( 0 , 1 ),
        scale ~ dexp(1),
        transpars> matrix[id,8]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[8,id]:z_id ~ normal( 0 , 1 ),
        vector[8]:sigma_id ~ dexp(1),
        cholesky_factor_corr[8]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[8,8]:Rho_id <<- Chol_to_Corr(L_Rho_id),
        # bernouli component GP priors
        transpars> vector[40]: a_gp <<- L_SIGMA_p * z_gp,
        vector[40]: z_gp ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_p <<- cholesky_decompose( SIGMA_p ),
        transpars> matrix[40,40]: SIGMA_p <- cov_GPL2( Dmat110_full , etasq_p , rhosq_p , 0.01 ),
        # gamma component GP priors
        transpars> vector[40]: a_gl <<- L_SIGMA_l * z_gl,
        vector[40]: z_gl ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA_l <<- cholesky_decompose( SIGMA_l ),
        transpars> matrix[40,40]: SIGMA_l <- cov_GPL2( Dmat110_full , etasq_l , rhosq_l , 0.01 ),
        c(etasq_p,rhosq_p,etasq_l,rhosq_l) ~ exponential(1)
    ) , 
    data=data_list_height_zero, chains=6, cores=6, iter=1000 , control=list(adapt_delta=0.97)  
)

precis(m1f , depth=3 , pars="Rho_id")
