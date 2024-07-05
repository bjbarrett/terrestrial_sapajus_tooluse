# TERRESTRIALITY AND PREDICTORS
PropPositive <- function(x) length(which(x > 0) ) /length(x)
PropPositive(post$bSp) 
PropPositive(post$bSl) 
PropPositive(post$bPTp) 
PropPositive(post$bPTl) 
precis(m1b)
post <- extract.samples(m1b)

median(logistic(post$ap)) #overall terrest rate

#general prob terr
png(file = "plots/prob_general_terr.png") 
    dens(logistic(post$ap) , xlim=c(0, 0.4 ) , xlab="probability of being terrestrial")
    blog=1.5
    for(i in 1:17){
        squanch <- logistic(post$ap + post$a_id[,i,1] )
        points( median(squanch) , (i-.1)/blog , pch=18 , col="darkgreen")
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    }
dev.off()
## terr rate
median(exp(post$al)) 
png(file = "plots/height_overall.png") 
    dens( (1-logistic(post$ap))*exp(post$al) , xlim=c(0, 20) , xlab="height in m " , show.HPDI=.99999 )
    blog=20
    for(i in 1:17){
        squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3])
        points( median(squanch) , (i-.1)/blog , pch=18)
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    }
dev.off()

# covariance of tools and space
png(file = "plots/covary_prob_terr_space_gp.png") 
    plot( NULL , xlab="distance (hectometers)" , ylab="covariance" ,
          xlim=c(0,20) , ylim=c(0,1) , main="spatial covariance of terrestriality probability")
    # compute posterior mean covariance
    x_seq <- seq( from=0 , to=10 , length.out=100 )
    pmcov <- sapply( x_seq , function(x) post$etasq_p*exp(-post$rhosq_p*x^2) )
    pmcov_mu <- apply( pmcov , 2 , mean )
    lines( x_seq , pmcov_mu , lwd=2 )
    # plot 60 functions sampled from posterior
    for ( i in 1:50 )
        curve( post$etasq_p[i]*exp(-post$rhosq_p[i]*x^2) , add=TRUE ,
               col=col.alpha("black",0.3) )
dev.off()

#no spatial covariance in height given they are not terrestrial
# plot the posterior median covariance function for lambda/gamma
png(file = "plots/covary_prob_height|not_terr_space_gp.png") 
    plot( NULL , xlab="distance (hectometers)" , ylab="covariance" ,
          xlim=c(0,20) , ylim=c(0,1) , main="spatial covariance of height given in canopy" )
    # compute posterior mean covariance
    x_seq <- seq( from=0 , to=10 , length.out=100 )
    pmcov <- sapply( x_seq , function(x) post$etasq_l*exp(-post$rhosq_l*x^2) )
    pmcov_mu <- apply( pmcov , 2 , mean )
    lines( x_seq , pmcov_mu , lwd=2 )
    # plot 60 functions sampled from posterior
    for ( i in 1:50 )
        curve( post$etasq_l[i]*exp(-post$rhosq_l[i]*x^2) , add=TRUE ,
               col=col.alpha("black",0.3) )
dev.off()

##height against distance

##make preds of post mean

##plot preds stones as a function of height
falafel <- range(apply (post$z_s ,2, mean))
stone_seq <-  seq(from=min(falafel) , to=max(falafel) , length=30)
# generate preds to plot
p_preds <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bPTp*0) ) 
p_med <- apply( p_preds , 2 , median )
lambda_preds <- sapply( stone_seq, function(x) exp(  post$al + post$bSl*x + post$bPTl*0) ) 
lambda_med <- apply( lambda_preds , 2 , median )
joint_preds <- (1-p_preds)*(lambda_preds)
joint_med <- apply( joint_preds , 2 , median )

png(file = "plots/height_sd_maineff.png") 
    plot(data_list_height2$height_m~data_list_height2$log_avg_stone_std , col=col.alpha("black",0.5)) # flag above 30 m
    plot(data_list_height2$height_m~data_list_height2$log_avg_stone_std , 
         col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="log(stone density per 110 m square grid)" ,
         ylab="height (m)") 
    
    #plot post median
    lines(joint_med ~ stone_seq , lw=2, col="slateblue4" , lty=1)
    for (j in sample( c(1:2000) , 100) ){
        lines( joint_preds[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue4", alpha=0.05) , lty=1)
    }
dev.off()

data_list_height2$terr <- ifelse(data_list_height2$height_m==0,1,0) 


png(file = "plots/terr_sd_maineff.png") 
    plot(data_list_height2$terr~data_list_height2$log_avg_stone_std , 
         col=col.alpha("black",0.01) , ylim=c(0,1) , xlab="log(stone density per 110 m square grid)" ,
         ylab="probability of being terrestrial" , pch=19) 
    #plot post median
    lines(p_med ~ stone_seq , lw=2, col="slateblue4" , lty=1)
    for (j in sample( c(1:2000) , 100) ){
        lines( p_preds[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue4", alpha=0.05) , lty=1)
    }
dev.off()

##plot palm density as a function of height

##plot preds stones as a function of height
palm_seq <-  seq(from=min(data_list_height2$count_palm_std) , to=max(data_list_height2$count_palm_std) , length=30)
# generate preds to plot
p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bSp*0 + post$bPTp*x) ) 
p_med <- apply( p_preds , 2 , median )
lambda_preds <- sapply( palm_seq, function(x) exp(  post$al + post$bSl*0 + post$bPTl*x) ) 
lambda_med <- apply( lambda_preds , 2 , median )
joint_preds <- (1-p_preds)*(lambda_preds)
joint_med <- apply( joint_preds , 2 , median )

png(file = "plots/height_palm_maineff.png") 
    plot(data_list_height2$height_m~data_list_height2$count_palm_std, 
         col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm density per 22 m square grid  (standardized)" ,
         ylab="height (m)") 
    
    #plot post median
    lines(joint_med ~ palm_seq , lw=2, col="forestgreen" , lty=1)
    for (j in sample( c(1:2000) , 100) ){
        lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("forestgreen", alpha=0.05) , lty=1)
    }
dev.off()

png(file = "plots/terr_palm_maineff.png") 
    plot(data_list_height2$terr~data_list_height2$count_palm_std, 
         col=col.alpha("black",0.01) , ylim=c(0,1) , xlab="palm density per 22 m square grid  (standardized)" ,
         ylab="probability of being terrestrial" , pch=19) 
    
    #plot post median
    lines(p_med ~ palm_seq , lw=2, col="forestgreen" , lty=1)
    for (j in sample( c(1:2000) , 100) ){
        lines( p_preds[j,] ~ palm_seq , lw=3, col=col.alpha("forestgreen", alpha=0.05) , lty=1)
    }
dev.off()

#height in first point scan from tool use
plot(precis(m1d))
post <- extract.samples(m1d)
str(post)

png(file = "plots/prob_terr_first_scan.png") 
    dens(logistic(post$ap) , xlim=c(0.05, 0.9 ) , xlab="probability of being terrestrial")
    dens(logistic(post$ap + post$bNp) , add=TRUE , col="slateblue" )
    legend(x = "topright", legend = c("tool use", "non tool use") , fill = c("slateblue","black") , bty='n' ) 
    blog=5
    for(i in 1:17){
        squanch <- logistic(post$ap + post$a_id[,i,1] )
        points( median(squanch) , (i-.1)/blog , pch=18)
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
        squanch <- logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,3])
        points( median(squanch) , (i+.1)/blog , col="slateblue" , pch=16)
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i+.1)/blog , y1=(i+.1)/blog , col="slateblue" )
    }
dev.off()

png(file = "plots/height_first_scan.png") 
    dens( (1-logistic(post$ap))*exp(post$al) , xlim=c(0, 15) , xlab="height (m) at beginning of focal follow" , show.HPDI=.99999 , ylim=c(0,1.1))
    dens((1-logistic(post$ap + post$bNp))*exp(post$al + post$bNl) , add=TRUE , col="slateblue" , show.HPDI=.99999)
    legend(x = "topright", legend = c("tool use", "non tool use") , fill = c("slateblue","black") , bty='n'  , lw=2) 
    blog=30
    
    for(i in 1:17){
        squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,3])
        points( median(squanch) , (i-.1)/blog , pch=18)
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
        squanch <- (1-logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2]))*exp(post$al + post$a_id[,i,3]+ post$bNl + post$a_id[,i,4])
        points( median(squanch) , (i+.1)/blog , col="slateblue" , pch=16)
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i+.1)/blog , y1=(i+.1)/blog , col="slateblue" )
    }
dev.off()


