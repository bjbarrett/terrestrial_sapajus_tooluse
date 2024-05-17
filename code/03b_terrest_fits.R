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

dens(logistic(post$ap) , xlim=c(0.05, 0.5 ))
dens(logistic(post$ap + post$bNp) , add=TRUE , col="slateblue" )

dens((1-logistic(post$ap))*exp(post$al) , xlim=c(0, 150 ))
dens((1-logistic(post$ap + post$bNp))*exp(post$al + post$bNl) , add=TRUE , col="slateblue" )



## stop here and build this up