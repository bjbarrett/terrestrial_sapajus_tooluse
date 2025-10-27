library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
data_whoa <- list(
    n_nuts = joined_df4$n_nuts,
    id = joined_df4$id,
    grid_id_follow = joined_df4$grid_id,# real grid id from follow, mising cell 15 cause no monkeys went there
    #grid_id_follow_dummy = as.integer(as.factor(joined_df4$grid_id)),# had to reassign ids so indexing works finwe with ulam
    log_follow_length=log(joined_df4$follow_length),
    
    height_m=blah$height,
    nutcracking=blah$nutcracking,
    id_h = blah$id,
    grid_id_follow_h = blah$grid_id# real grid id from follow, mising cell 15 cause no monkeys went there
    # grid_id_follow_dummy_h = as.integer(as.factor(blah$grid_id))# had to reassign ids so indexing works finwe with ulam
)

sort(unique(data_whoa$id))
sort(unique(data_whoa$gid))

sort(unique(data_whoa$grid_id_follow_h))
sort(unique(data_whoa$grid_id_follow))

str(data_whoa)

# m_cov <- ulam(
#     alist(
#             n_nuts ~ dzipois( p , lambda ),
#             logit(p) <- apn + a_id[id,1] + a_g110[grid_id_follow_dummy,1],
#             log(lambda) <- aln + a_id[id,2] +  a_g110[grid_id_follow_dummy,2] + log_follow_length,
#             
#             height_m|height_m==0 ~ custom( bernoulli_lpmf(1|p) ) ,
#             height_m|height_m>0 ~ custom( bernoulli_lpmf(0|p) + gamma_lpdf( height_m | mu/scale , 1/scale) ) ,
#             logit(p) <- aph + a_id[gid,3]  +  a_g110[grid_id_follow_dummy_h,3] + (bnp + a_id[gid,4]  +  a_g110[grid_id_follow_dummy_h,4] )*nutcracking,
#             log(mu) <- alh + a_id[gid,5] +  a_g110[grid_id_follow_dummy_h,5] + (bnl + a_id[gid,6]  +  a_g110[grid_id_follow_dummy_h,6])*nutcracking,
#             
#             
#             # adaptive priors - non-centered
#             transpars> matrix[id,6]:a_id <-
#                 compose_noncentered( sigma_id , L_Rho_id , z_id ),
#             matrix[6,id]:z_id ~ normal( 0 , 1 ),
#             
#             transpars> matrix[grid_id_follow,6]:a_g110 <-
#                 compose_noncentered( sigma_gid , L_Rho_gid , z_gid ),
#             
#             matrix[6,grid_id_follow]:z_gid ~ normal( 0 , 1 ),
#             # fixed priors
#             c(apn,aln,aph,alh,bnl,bnp) ~ dnorm( 0 , 1 ),
#             vector[6]:sigma_id ~ dexp(1),
#             vector[6]:sigma_gid ~ dexp(1),
#             cholesky_factor_corr[6]:L_Rho_id ~ lkj_corr_cholesky( 4 ),
#             cholesky_factor_corr[6]:L_Rho_gid ~ lkj_corr_cholesky( 4 ),
#             scale ~dexp(1),
#             # compute ordinary correlation matrixes from Cholesky factors
#             gq> matrix[6,6]:Rho_id <<- Chol_to_Corr(L_Rho_id),
#             gq> matrix[6,6]:Rho_gid <<- Chol_to_Corr(L_Rho_gid)
#             
#         
#     ) , data=data_whoa, chains=4 , cores=4, log_lik=FALSE, iter=2000 , control=list(adapt_delta=0.95))

    
scode <- "

"
mod <- cmdstan_model(scode)

file <- file.path("multivargrid.stan")
mod <- cmdstan_model(file , cpp_options = list(stan_threads = TRUE) )
fit_multi <- mod$sample(
    data = data_whoa,
    seed = 1303,
    adapt_delta = 0.96,
    chains = 4,
    parallel_chains = 4, 
    refresh = 100,
    iter_sampling = 1000,
    iter_warmup = 1000,
    threads=4,
    max_treedepth = 11
)

poster <- extract.samples(fit_multi)


precis(fit_multi , pars="Rho_id" , depth=3)


plot (precis(fit_multi , pars="Rho_grid" , depth=3) )

stanfit <- rstan::read_stan_csv(fit_multi$output_files())
post_multi <- extract(stanfit)
save(stanfit , file="fit_multi.rds")

str(poster)
wamp <- precis(fit_multi , pars="Rho_id" , depth=3)
wamp <- precis(fit_multi , pars=c("bnp","bnl","alh","aph",
        "aln","apn" ,"sigma_id" ,"sigma_grid" ,"scale" , "Rho_grid" , "Rho_id" ) , 
        depth=3 , digits=2)
row.names(wamp)
write.csv ( wamp[1:4] , file="tables/fit_mutli_rho_grid_gen_terr.csv")

precis(fit_multi , pars="a_g110" , depth=3)

precis(fit_ul , pars="a_g110" , depth=3)


# We predict that tool use is more frequent in regions where general terrestriality is higher. 
# So plot tool use rate in a cell  by height in cells where t, do varying effects of individuals and plot id look at covariance of 
# 
# We also predict that general terrestriality (ap_cell) is higher and spatially correlated with tool use terrestriality (bp_cell) 


###TU rate, give TU x general terrestriality


png(file = "plots/corr_mat_tu.png",  width = 7, height = 7,units = "in" , res=300) 
    par(mar=c(2,2,.2,.2) , mfrow=c(2,2) , oma=c(2 , 2 ,.25, .25))
    ###prob TU x general terrestriality
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,3] , poster$a_g110[1:50,i,1] , col="grey" , cex=0.5 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,3]) , median(poster$a_g110[,i,1]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    mtext("ap_t,grid (probability not using tools)" , side=2 , line=2.5 , outer=FALSE , cex=1)
    text(-2.8,2.8,"a.")
    ###prob TU x height, given terr height
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,5] , poster$a_g110[1:50,i,1] , col="grey" , cex=0.75 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,5]) , median(poster$a_g110[,i,1]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    text(-2.8,2.8,"b.")
    ###TU rate, give TU x general terrestriality
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,3] , poster$a_g110[1:50,i,2] , col="grey" , cex=0.75 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,3]) , median(poster$a_g110[,i,2]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    mtext("ap_h,grid (probability terrestriality)" , side=1 , line=2.5 , outer=FALSE , cex=1)
    mtext("al_t,grid (tool use rate | tool use)" , side=2 , line=2.5 , outer=FALSE , cex=1)
    text(-2.8,2.8,"c.")
    ###TU rate, give TU x general terrestriality
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,5] , poster$a_g110[1:50,i,2] , col="grey" , cex=0.75 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,5]) , median(poster$a_g110[,i,2]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    mtext("al_h,grid (height | not terrestrial)" , side=1 , line=2.5 , outer=FALSE , cex=1)
    text(-2.8,2.8,"d.")
dev.off()

# Rho_grid[1,3]  0.28 0.17 -0.02  0.53 1.00  1353.78
# Rho_grid[5,1]  0.10 0.18 -0.20  0.40 1.00  1628.03
# Rho_grid[3,2] -0.01 0.26 -0.44  0.40 1.01   475.64
# Rho_grid[5,2] -0.09 0.27 -0.49  0.34 1.00   479.04


png(file = "plots/corr_mat_genterr_toolterr.png",  width = 7, height = 7,units = "in" , res=300) 
    par(mar=c(2,2,.2,.2) , mfrow=c(2,2) , oma=c(2 , 2 ,.25, .25))
    ###prob TU x general terrestriality
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,3] , poster$a_g110[1:50,i,4] , col="grey" , cex=0.5 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,3]) , median(poster$a_g110[,i,4]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    mtext("bNp,grid (probability tool terrestriality)" , side=2 , line=2.5 , outer=FALSE , cex=1)
    text(-2.8,2.8,"a.")
    ###prob TU x height, given terr height
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,5] , poster$a_g110[1:50,i,4] , col="grey" , cex=0.75 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,5]) , median(poster$a_g110[,i,4]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    text(-2.8,2.8,"b.")
    ###TU rate, give TU x general terrestriality
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,3] , poster$a_g110[1:50,i,6] , col="grey" , cex=0.75 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,3]) , median(poster$a_g110[,i,6]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    mtext("ap_h,grid (probability general terrestriality)" , side=1 , line=2.5 , outer=FALSE , cex=1)
    mtext("bNl,grid (tool height | not terrestrial)" , side=2 , line=2.5 , outer=FALSE , cex=1)
    text(-2.8,2.8,"c.")
    ###TU rate, give TU x general terrestriality
    plot("NA" , xlim=c(-3,3) , ylim=c(-3,3))
    for(i in 1:40) points(poster$a_g110[1:50,i,5] , poster$a_g110[1:50,i,6] , col="grey" , cex=0.75 , pch="*")
    for(i in 1:40) points(median(poster$a_g110[,i,5]) , median(poster$a_g110[,i,6]) , col=rangi2 , pch=1 )
    abline(a=0 , b=1 , lty=3)
    mtext("al_h,grid (general height | not terrestrial)" , side=1 , line=2.5 , outer=FALSE , cex=1)
    text(-2.8,2.8,"d.")
dev.off()

# Rho_grid[4,3] -0.50 0.19 -0.77 -0.16 1.00  3652.85
# Rho_grid[4,5]  0.25 0.20 -0.09  0.56 1.00  1573.25
# Rho_grid[3,6]  0.06 0.22 -0.30  0.40 1.00  6022.47
# Rho_grid[5,6] -0.42 0.21 -0.73 -0.07 1.00  4910.12

