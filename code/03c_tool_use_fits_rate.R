#this does rates of tool use as a function of stuff, not just probabailty
#03c_tooluse fits 

m1zip <- ulam(
    alist(
        n_nuts ~ dzipois( p , lambda ),
        logit(p) <- ap + a_id[id,1],
        log(lambda) <- al + a_id[id,2] + log_follow_length,
        # adaptive priors - non-centered
        transpars> matrix[id,2]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[2,id]:z_id ~ normal( 0 , 1 ),
        # fixed priors
        c(ap,al) ~ dnorm( 0 , 1 ),
        vector[2]:sigma_id ~ dexp(1),
        cholesky_factor_corr[2]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[2,2]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , data=data_list_daily , iter=2000 , chains=4 , cores=4 , log_lik=FALSE, pars_omit = 'z_id')

wamp <- precis(m1zip , depth=2, digits=2)
write.csv ( wamp , file="tables/tu_rate_m1zip.csv")


m7zip <- ulam(
    alist(
        n_nuts ~ dzipois( p , lambda ),
        logit(p) <- ap + a_id[id,1] + (bPTp + a_id[id,2])*count_palm_std,
        log(lambda) <- al + a_id[id,3] + (bPTl + a_id[id,4])*count_palm_std + log_follow_length,
        # adaptive priors - non-centered
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        # fixed priors
        c(ap,al,bPTp,bPTl) ~ dnorm( 0 , 1 ),
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , data=data_list_daily , chains=4 , cores=4 , log_lik=FALSE, pars_omit = 'z_id')

precis(m7zip)

#################m7zip###########

wamp <- precis(m7zip , depth=2, digits=2)
write.csv ( wamp , file="tables/tu_stone_m7zip.csv")

post <- extract.samples(m7zip)
names(post)

png(file = "plots/precis_m7zip.png",  width = 6, height = 6,units = "in" , res=300) 
plot( precis(m7zipGP , depth=3 , 
             pars=c( "ap" , "bPTp" , "al" , "bPTl" , "sigma_id") ),
      labels=c("ap" , "bPTp" , "al" , "bPTl" , "sigma_id_ap" , "sigma_id_bSp" , "sigma_id_al" , "sigma_id_bSl" ) ) 
dev.off()

#plot predictions relative to raw data on rate scale

p_link_palm <- function(x){
    id <- 1
    logit_scale <- with(post , ap + a_id[,id,1]*0 + (bPTp + a_id[,id,2])*x)
    return(logistic(logit_scale))
}

l_link_palm <- function(x){
    id <- 1
    log_scale <- with(post , al + a_id[,id,3]*0 +  (bPTl + a_id[,id,4])*x)
    return(exp(log_scale))
}

palm_seq <-  seq(from=min(data_list_daily_2$count_palm_std) , 
                  to=max(data_list_daily_2$count_palm_std) , length=30)

p_raw <- sapply( palm_seq, function(x) (1-p_link_palm(x) )*l_link_palm(x) )

data_list_daily_2$nut_per_min <- data_list_daily_2$n_nuts/data_list_daily_2$follow_length

png(file = "plots/crack_rate_pop_mean_m7zip.png" ,  width = 6, height = 6,units = "in" , res=300) 
plot(data_list_daily_2$nut_per_min ~ data_list_daily_2$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.05), ylab="nut-cracking rate (per minute)",
     xlab= "palm tree density (standardized)" , ylim=c(0,1.6))

pred_median <- apply(p_raw , 2 , median)
lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 50) ){
    lines( p_raw[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
}
dev.off()

###plot predictions of prob of tool use
png(file = "plots/prob_crack_pop_mean_m7zip.png" ,  width = 6, height = 6,units = "in" , res=300) 
p_raw <- sapply( palm_seq, function(x) (1-p_link_palm(x) ) )

plot(data_list_daily_2$nutcracking ~ data_list_daily_2$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.05), ylab="probability of nut-cracking",
     xlab= "palm tree density (standardized)")

pred_median <- apply(p_raw , 2 , median)
lines(pred_median ~ palm_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 50) ){
    lines( p_raw[j,] ~ palm_seq , lw=3, col=col.alpha("green4", alpha=0.1) , lty=1)
}

dev.off()


###lets get per individual prediction
png(file = "plots/crack_rate_individuals_m7zip.png" ,  width = 6, height = 6,units = "in" , res=300) 

p_link_palm <- function(x,id){
    logit_scale <- with(post , ap + a_id[,id,1] +  (bPTp + a_id[,id,2])*x)
    return(logistic(logit_scale))
}

l_link_palm <- function(x,id){
    log_scale <- with(post , al + a_id[,id,3]*0 +  (bPTl + a_id[,id,4])*x)
    return(exp(log_scale))
}

palm_seq <-  seq(from=min(data_list_daily$count_palm_std) , 
                  to=max(data_list_daily$count_palm_std) , length=30)


plot(data_list_daily_2$nut_per_min ~ data_list_daily_2$count_palm_std ,
     pch=19 , col=col.alpha("green4",0.05), ylab="nut-cracking rate (per minute)",
     xlab= "palm tree density/ 22 m^2 grid (standardized)" , ylim=c(0,1.6) , cex.lab=1.5)

for(mono in 1:max(data_list_daily_2$id)){
    p_raw <- sapply( palm_seq, function(x) (1-p_link_palm(x , id=mono) )*l_link_palm(x,id=mono) )
    pred_median <- apply(p_raw , 2 , median)
    lines(pred_median ~ palm_seq , lw=2, col=col.alpha("green4", alpha=0.8) , lty=1)
}

dev.off()


# m8zip <- ulam(
#     alist(
#         n_nuts ~ dzipois( p , lambda ),
#         logit(p) <- ap + a_id[id,1] +  (bSp + a_id[id,2])*log_avg_stone_std,
#         log(lambda) <- al + a_id[id,3] +  (bSl + a_id[id,4])*log_avg_stone_std + log_follow_length,
#         # adaptive priors - non-centered
#         transpars> matrix[id,4]:a_id <-
#             compose_noncentered( sigma_id , L_Rho_id , z_id ),
#         matrix[4,id]:z_id ~ normal( 0 , 1 ),
#         # fixed priors
#         c(ap,al,bSp,bSl) ~ dnorm( 0 , 1 ),
#         vector[4]:sigma_id ~ dexp(1),
#         cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
#         # compute ordinary correlation matrixes from Cholesky factors
#         gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
#     ) , data=data_list_daily , chains=4 , cores=4 , log_lik=FALSE, pars_omit = 'z_id')
# 
# precis(m8zip , depth=3)
# 
# m9zip <- ulam(
#     alist(
#         n_nuts ~ dzipois( p , lambda ),
#         logit(p) <- ap + a_id[id,1] +  (bPTp + a_id[id,2])*count_palm_std +  (bSp + a_id[id,3])*log_avg_stone_std,
#         log(lambda) <- al + a_id[id,4] + (bPTl + a_id[id,5])*count_palm_std +  (bSl + a_id[id,6])*log_avg_stone_std + log_follow_length,
#         # adaptive priors - non-centered
#         transpars> matrix[id,6]:a_id <-
#             compose_noncentered( sigma_id , L_Rho_id , z_id ),
#         matrix[6,id]:z_id ~ normal( 0 , 1 ),
#         # fixed priors
#         c(ap,al,bSp,bSl,bPTp,bPTl) ~ dnorm( 0 , 1 ),
#         vector[6]:sigma_id ~ dexp(1),
#         cholesky_factor_corr[6]:L_Rho_id ~ lkj_corr_cholesky( 3 ),
#         # compute ordinary correlation matrixes from Cholesky factors
#         gq> matrix[6,6]:Rho_id <<- Chol_to_Corr(L_Rho_id)
#     ) , data=data_list_daily , chains=4 , cores=4 , log_lik=FALSE, pars_omit = 'z_id')

# 
# m7zipGP <- ulam(
#     alist(
#         #stone model
#         stone_raw_grid ~ dpois( lambda_stone ), #define poisson likelihood
#         log(lambda_stone) <-  a_bar + a_grid[stone_grid_index],
#         #priors
#         a_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
#         etasq ~ dexp( 2 ),
#         rhosq ~ dexp( 0.5 ),
#         # non-centered Gaussian Process prior
#         transpars> vector[40]: a_grid <<- L_SIGMA * z,
#         vector[40]: z ~ normal( 0 , 1 ),
#         transpars> matrix[40,40]: L_SIGMA <<- cholesky_decompose( SIGMA ),
#         transpars> matrix[40,40]: SIGMA <- cov_GPL2( Dmat110 , etasq , rhosq , 0.01 ),
#         ###nutcracking model
#         n_nuts ~ dzipois( p , lambda ),
#         logit(p) <- ap + a_id[id,1] +  (bPTp + a_id[id,2])*count_palm_std,
#         log(lambda) <- al + a_id[id,3] +  (bPTl + a_id[id,4])*count_palm_std + log_follow_length,
#         # adaptive priors - non-centered
#         transpars> matrix[id,4]:a_id <-
#             compose_noncentered( sigma_id , L_Rho_id , z_id ),
#         matrix[4,id]:z_id ~ normal( 0 , 1 ),
#         # fixed priors
#         c(ap,al,bPTp,bPTl) ~ dnorm( 0 , 1 ),
#         vector[4]:sigma_id ~ dexp(1),
#         cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
#         # compute ordinary correlation matrixes from Cholesky factors
#         gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
#     ) , data=data_list_daily_2 , chains=4 , cores=4 , log_lik=FALSE, pars_omit = 'z_id')
# 
# precis(m7zipGP)



### this is what we use to estimate the effects of stones
m8zipGP <- ulam(
    alist(
        n_nuts ~ dzipois( p , lambda ),
        logit(p) <- ap + a_id[id,1] +  (bSp + a_id[id,2])*s_grid[grid_id_follow],
        log(lambda) <- al + a_id[id,3] +  (bSl + a_id[id,4])*s_grid[grid_id_follow] + log_follow_length,
        # adaptive priors - non-centered
        transpars> matrix[id,4]:a_id <-
            compose_noncentered( sigma_id , L_Rho_id , z_id ),
        matrix[4,id]:z_id ~ normal( 0 , 1 ),
        # fixed priors
        c(bSl,al,bSp,ap) ~ dnorm( 0 , 1 ),
        vector[4]:sigma_id ~ dexp(1),
        cholesky_factor_corr[4]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
        
        #stone model
        stone_raw_grid ~ dpois( lambda_stone ), #define poisson likelihood
        log(lambda_stone) <-  s_bar + s_grid[stone_grid_index],
        #priors
        s_bar ~ dnorm( 3 , 2 ), #mean effect on logscale
        etasq ~ dexp( 2 ),
        rhosq ~ dexp( 0.5 ),
        # non-centered Gaussian Process prior
        transpars> vector[40]: s_grid <<- L_SIGMA * z,
        vector[40]: z ~ normal( 0 , 1 ),
        transpars> matrix[40,40]: L_SIGMA <<- cholesky_decompose( SIGMA ),
        transpars> matrix[40,40]: SIGMA <- cov_GPL2( Dmat110 , etasq , rhosq , 0.01 ),
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[4,4]:Rho_id <<- Chol_to_Corr(L_Rho_id)
    ) , data=data_list_daily_2 , chains=4 , cores=4 , log_lik=FALSE, pars_omit = 'z_id' , iter=2000)

wamp <- precis(m8zipGP , depth=2, digits=2)
write.csv ( wamp , file="tables/tu_stone_m8zipGP.csv")

post <- extract.samples(m8zipGP)
names(post)

png(file = "plots/precis_m8zipGP.png" ,  width = 6, height = 6,units = "in" , res=300 ) 
plot( precis(m8zipGP , depth=3 , 
             pars=c( "ap" , "bSp" , "al" , "bSl" , "sigma_id" , "s_bar" ,"etasq" ,"rhosq") ),
             labels=c("ap" , "bSp" , "al" , "bSl" , "sigma_id_ap" , "sigma_id_bSp" , "sigma_id_al" , "sigma_id_bSl" ,"s_bar" ,"etasq" ,"rhosq" ) ) 
dev.off()


##stone density covariance
png(file = "plots/covar_GP_stonedens_m8zipGP.png" ,  width = 6, height = 6,units = "in" , res=300) 
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
dev.off()

#plot predictions relative to raw data on rate scale

p_link_stone <- function(x){
    id <- 1
    logit_scale <- with(post , ap + a_id[,id,1]*0 +  (bSp + a_id[,id,2])*x)
    return(logistic(logit_scale))
}

l_link_stone <- function(x){
    id <- 1
    log_scale <- with(post , al + a_id[,id,3]*0 +  (bSl + a_id[,id,4])*x)
    return(exp(log_scale))
}

stone_seq <-  seq(from=min(data_list_daily$log_avg_stone_std) , 
                  to=max(data_list_daily$log_avg_stone_std) , length=30)

p_raw <- sapply( stone_seq, function(x) (1-p_link_stone(x) )*l_link_stone(x) )

data_list_daily_2$nut_per_min <- data_list_daily_2$n_nuts/data_list_daily_2$follow_length

png(file = "plots/crack_rate_pop_mean_m8zipGP.png"  ,  width = 6, height = 6,units = "in" , res=300) 
plot(data_list_daily_2$nut_per_min ~ data_list_daily_2$log_avg_stone_std ,
     pch=19 , col=col.alpha("darkslateblue",0.05), ylab="nut-cracking rate (per minute)",
     xlab= "log stone density (standardized)" , ylim=c(0,1.6))

pred_median <- apply(p_raw , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 50) ){
    lines( p_raw[j,] ~ stone_seq , lw=3, col=col.alpha("darkslateblue", alpha=0.1) , lty=1)
}
dev.off()

###plot predictions of prob of tool use
png(file = "plots/prob_crack_pop_mean_m8zipGP.png" ,  width = 6, height = 6,units = "in" , res=300) 
p_raw <- sapply( stone_seq, function(x) (1-p_link_stone(x) ) )

plot(data_list_daily_2$nutcracking ~ data_list_daily_2$log_avg_stone_std ,
     pch=19 , col=col.alpha("darkslateblue",0.05), ylab="probability of nut-cracking",
     xlab= "log stone density (standardized)")

pred_median <- apply(p_raw , 2 , median)
lines(pred_median ~ stone_seq , lw=2, col=1 , lty=1)

for (j in sample( c(1:2000) , 50) ){
    lines( p_raw[j,] ~ stone_seq , lw=3, col=col.alpha("darkslateblue", alpha=0.1) , lty=1)
}

dev.off()


###lets get per individual prediction
png(file = "plots/crack_rate_individuals_m8zipGP.png"  ,  width = 6, height = 6,units = "in" , res=300) 

p_link_stone <- function(x,id){
    logit_scale <- with(post , ap + a_id[,id,1] +  (bSp + a_id[,id,2])*x)
    return(logistic(logit_scale))
}

l_link_stone <- function(x,id){
    log_scale <- with(post , al + a_id[,id,3]*0 +  (bSl + a_id[,id,4])*x)
    return(exp(log_scale))
}

stone_seq <-  seq(from=min(data_list_daily$log_avg_stone_std) , 
                  to=max(data_list_daily$log_avg_stone_std) , length=30)


plot(data_list_daily_2$nut_per_min ~ data_list_daily_2$log_avg_stone_std ,
     pch=19 , col=col.alpha("darkslateblue",0.05), ylab="nut-cracking rate (per minute)",
     xlab= "log stone density/ 110 m^2 grid (standardized)" , ylim=c(0,1.6) , cex.lab=1.5)

for(mono in 1:max(data_list_daily_2$id)){
    p_raw <- sapply( stone_seq, function(x) (1-p_link_stone(x , id=mono) )*l_link_stone(x,id=mono) )
    pred_median <- apply(p_raw , 2 , median)
    lines(pred_median ~ stone_seq , lw=2, col=col.alpha("darkslateblue", alpha=0.8) , lty=1)
    }

dev.off()


#### lets do the palm abundnance me model

