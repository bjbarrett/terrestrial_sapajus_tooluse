######TU RATES
post <- extract.samples(m1zip)

PropPositive <- function(x) length(which(x > 0) ) /length(x)
post <- extract.samples(m1zip)
median(post$ap)
png(file = "plots/tu_rate_m1zip.png",  width = 6, height = 6,units = "in" , res=300) 
    par(mar=c(3,0.1,0.1,0.1))
    
    trout <- (1-logistic(post$ap))*exp(post$al)
    dens(trout , xlim=c(0, 0.2) , col="NA" , xlab="nut-cracking rate (per minute)", yaxt="n" , ylab="" , cex.lab=1.5)
    shade( density(trout) ,
           PCI(trout , prob=0.9999) ,
           col=col.alpha("purple1" , alpha=0.75) )
    blog=0.7
    
    for(i in 1:17){
        squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,2])
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
        points( median(squanch) , (i-.1)/blog , pch=18 , cex=0.75 )
    }
    # points(data_list_height$height_m ,rep(-0.01, length(data_list_height$height_m)) ,
    #        pch=19 , cex=0.5 , col=col.alpha("purple1" , alpha=0.05))

dev.off()

median(trout)
HPDI(trout)
png(file = "plots/tu_prob_m1zip.png",  width = 6, height = 6,units = "in" , res=300) 
    par(mar=c(3,0.1,0.1,0.1))
    
    trout <- (1-logistic(post$ap))
    dens(trout , xlim=c(0, 1) , col="NA" , xlab="probability of tool use in a follow", yaxt="n" , ylab="" , cex.lab=1.5)
    shade( density(trout) ,
           PCI(trout , prob=0.9999) ,
           col=col.alpha("purple1" , alpha=0.75) )
    blog=3
    
    for(i in 1:17){
        squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
        points( median(squanch) , (i-.1)/blog , pch=18 , cex=0.75 )
    }
    # points(data_list_height$height_m ,rep(-0.01, length(data_list_height$height_m)) ,
    #        pch=19 , cex=0.5 , col=col.alpha("purple1" , alpha=0.05))
dev.off()




######mzag 1#################################################basic varef id terr rats
post <- extract.samples(mzag_1)
median(logistic(post$ap)) #overall terrest rate
HPDI(logistic(post$ap))

wamp <- precis(mzag_1 , depth=2, digits=2)
write.csv ( wamp[1:4] , file="tables/mzag1.csv")

precis(mzag_1 , depth=2)

for( i in 1 :17){ (print(median(logistic(post$ap + post$a_id[,i,1]) ) )) }

#general prob terr
png(file = "plots/prob_terr_mzag1.png",  width = 6, height = 6,units = "in" , res=300) 
trout <- logistic(post$ap)
par(mar=c(3,0.1,0.1,0.1))
dens(trout , xlim=c(0, 0.4) , col="NA" , xlab="probability of being terrestrial", yaxt="n" , ylab="" , cex.lab=1.5)
shade( density(trout) ,
       PCI(trout , prob=0.9999) ,
       col=col.alpha("purple1" , alpha=0.75) )
blog=1.5
for(i in 1:17){
    squanch <- logistic(post$ap + post$a_id[,i,1] )
    points( median(squanch) , (i-.1)/blog , pch=18 )
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
}
dev.off()

## height overall
median((1-logistic(post$ap))*exp(post$al))
HPDI ( (1-logistic(post$ap))*exp(post$al) )

png(file = "plots/height_overall_mzag1.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(3,0.1,0.1,0.1))

trout <- (1-logistic(post$ap))*exp(post$al)
dens(trout , xlim=c(0, 20) , col="NA" , xlab="height (m)", yaxt="n" , ylab="" , cex.lab=1.5)
shade( density(trout) ,
       PCI(trout , prob=0.9999) ,
       col=col.alpha("purple1" , alpha=0.75) )
blog=20

for(i in 1:17){
    squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,2])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=18 , cex=0.75 )
}
# points(data_list_height$height_m ,rep(-0.01, length(data_list_height$height_m)) ,
#        pch=19 , cex=0.5 , col=col.alpha("purple1" , alpha=0.05))

dev.off()

#####################mzag2########################
post <- extract.samples(mzag_2)
median(logistic(post$ap)) #overall terrest rate
HPDI(logistic(post$ap))

median(logistic(post$ap + post$bNp)) #overall terrest rate
HPDI(logistic(post$ap + post$bNp))

wamp <- precis(mzag_2 , depth=2, digits=2)
write.csv ( wamp[1:4] , file="tables/mzag2.csv")

#general prob terr tu ntu
png(file = "plots/prob_general_terr_tu_ntu.png",  width = 6, height = 6,units = "in" , res=300) 
trout <- logistic(post$ap)
par(mar=c(3,0.1,0.1,0.1))
dens(trout , xlim=c(0, 0.75) , col="NA" , xlab="probability of being terrestrial", yaxt="n" , ylab="" , cex.lab=1.5)
shade( density(trout) ,
       PCI(trout , prob=0.9999) ,
       col=col.alpha("purple1" , alpha=0.75) )

pike <- logistic(post$ap + post$bNp)
shade( density(pike) ,
       PCI(pike , prob=0.9999) ,
       col=col.alpha("orange2" , alpha=0.75) )

blog=1.5
for(i in 1:17){
    squanch <- logistic(post$ap + post$a_id[,i,1] )
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=23 ,col="black", bg="purple1" , cex=0.75)
}

for(i in 1:17){
    squanch <- logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=23, col="black", bg="orange2", cex=0.75)
}

legend("topright" , legend=c("no tool-use" , "tool-use") , fill=c("purple1" , "orange2") , bty='n' )


dev.off()

median ((1- logistic(post$ap) )*exp(post$al))
HPDI ((1- logistic(post$ap) )*exp(post$al))
pike <- (1- logistic(post$ap + post$bNp) )*exp(post$al + post$bNl)
median(pike)
HPDI(pike)
#general prob terr tu ntu
png(file = "plots/turate_general_terr_tu_ntu.png",  width = 6, height = 6,units = "in" , res=300)

trout <- (1- logistic(post$ap) )*exp(post$al)
dens(trout , xlim=c(0, 16) , ylim=c(0,3) ,col="NA" , xlab="average height (m) per follow", yaxt="n" , ylab="" , cex.lab=1.5)
shade( density(trout) ,
       PCI(trout , prob=0.9999) ,
       col=col.alpha("purple1" , alpha=0.75) )

pike <- (1- logistic(post$ap + post$bNp) )*exp(post$al + post$bNl)
shade( density(pike) ,
       PCI(pike , prob=0.9999) ,
       col=col.alpha("orange2" , alpha=0.75) )

blog=15
for(i in 1:17){
    squanch <- (1- logistic(post$ap + post$a_id[,i,1]) )*exp(post$al + post$a_id[,i,3])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=23 ,col="black", bg="purple1" , cex=0.75)
}

for(i in 1:17){
    squanch <- (1- logistic(post$ap + post$a_id[,i,1] + post$bNp + post$a_id[,i,2]) )*exp(post$al + post$a_id[,i,3] + post$bNl + post$a_id[,i,4])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=23, col="black", bg="orange2", cex=0.75)
}
legend("topright" , legend=c("no tool-use" , "tool-use") , fill=c("purple1" , "orange2") , bty='n' )

dev.off()

# covariance of tools and space
png(file = "plots/covary_prob_terr_space_gp_mzag_2.png", width = 6, height = 6,units = "in" , res=300) 
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

png(file = "plots/covary_prob_height|not_terr_space_gp_mzag_2.png") 
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


####height overall nutcracking follows

##height against distance

##make preds of post mean
# 
# ##plot preds stones as a function of height
# falafel <- range(apply (post$z_s ,2, mean))
# stone_seq <-  seq(from=min(falafel) , to=max(falafel) , length=30)
# # generate preds to plot
# p_preds <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bPTp*0) ) 
# p_med <- apply( p_preds , 2 , median )
# lambda_preds <- sapply( stone_seq, function(x) exp(  post$al + post$bSl*x + post$bPTl*0) ) 
# lambda_med <- apply( lambda_preds , 2 , median )
# joint_preds <- (1-p_preds)*(lambda_preds)
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/height_sd_maineff.png") 
# plot(data_list_height2$height_m~data_list_height2$log_avg_stone_std , col=col.alpha("black",0.5)) # flag above 30 m
# plot(data_list_height2$height_m~data_list_height2$log_avg_stone_std , 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="log(stone density per 110 m square grid)" ,
#      ylab="height (m)") 
# 
# #plot post median
# lines(joint_med ~ stone_seq , lw=2, col="slateblue4" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue4", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# data_list_height2$terr <- ifelse(data_list_height2$height_m==0,1,0) 
# 
# 
# png(file = "plots/terr_sd_maineff.png") 
# plot(data_list_height2$terr~data_list_height2$log_avg_stone_std , 
#      col=col.alpha("black",0.01) , ylim=c(0,1) , xlab="log(stone density per 110 m square grid)" ,
#      ylab="probability of being terrestrial" , pch=19) 
# #plot post median
# lines(p_med ~ stone_seq , lw=2, col="slateblue4" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( p_preds[j,] ~ stone_seq , lw=3, col=col.alpha("slateblue4", alpha=0.05) , lty=1)
# }
# dev.off()

############# m1f PALM DENSITY TU/NTU GP AND TERRESTRIALITY###################

post <- extract.samples(m1f)
wamp <- precis(m1f , depth=2, digits=2)
write.csv ( wamp[1:4] , file="tables/m1f.csv")

palm_seq <-  seq(from=min(data_list_height2$count_palm_std) , to=max(data_list_height2$count_palm_std) , length=30)
# generate preds to plot

p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*0 + post$bNxPTp*0*x ) )
p_med <- apply( p_preds , 2 , median )
lambda_preds <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*0 + post$bNxPTl*0*x ) )
lambda_med <- apply( lambda_preds , 2 , median )
joint_preds <- (1-p_preds)*(lambda_preds)
joint_med <- apply( joint_preds , 2 , median )

png(file = "plots/height_palm_maineff_ntu_m1f.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==0], 
     col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm tree density (standardized)" ,
     ylab="height (m)" , cex.lab=1.5) 

#plot post median
lines(joint_med ~ palm_seq , lw=2, col="purple1" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
}
text(-0.7,29.5,"b.")

dev.off()

p_preds2 <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*1 + post$bNxPTp*1*x ) )
p_med2 <- apply( p_preds2 , 2 , median )
lambda_preds2 <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*1 + post$bNxPTl*1*x ) )
lambda_med2 <- apply( lambda_preds , 2 , median )
joint_preds2 <- (1-p_preds2)*(lambda_preds2)
joint_med2 <- apply( joint_preds2 , 2 , median )

png(file = "plots/height_palm_maineff_tu_m1f.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
    plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==1], 
         col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm tree density (standardized)" ,
         ylab="height (m)" , cex.lab=1.5) 
    lines(joint_med2 ~ palm_seq , lw=2, col="orange2" , lty=1)
    for (j in sample( c(1:2000) , 100) ){
        lines( joint_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
    }
    text(-0.7,29.5,"d.")
dev.off()

####prob
p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*0 + post$bNxPTp*0*x ) )
p_med <- apply( p_preds , 2 , median )
joint_preds <- p_preds
joint_med <- apply( joint_preds , 2 , median )

png(file = "plots/prob_terr_palm_maineff_ntu_m1f.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==0], 
     col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="palm tree density (standardized)" ,
     ylab="probability of terrestriality" , cex.lab=1.5) 

#plot post median
lines(joint_med ~ palm_seq , lw=2, col="purple1" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
}
text(-0.7,1,"a.")

dev.off()

p_preds2 <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*1 + post$bNxPTp*1*x ) )
p_med2 <- apply( p_preds2 , 2 , median )
joint_preds2 <- p_preds2
joint_med2 <- apply( joint_preds2 , 2 , median )

png(file = "plots/prob_height_palm_maineff_tu_m1f.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==1], 
     col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="palm tree density (standardized)" ,
     ylab="probability of terrestriality" , cex.lab=1.5) 
lines(joint_med2 ~ palm_seq , lw=2, col="orange2" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
}
text(-0.7,1,"c.")
dev.off()
# 
# ########SIMPLER NON GP MODEL.    m1fsimp.   #############
# post <- extract.samples(m1fsimp)
# wamp <- precis(m1fsimp , depth=2, digits=2)
# write.csv ( wamp[1:4] , file="tables/m1fsimp.csv")
# 
# palm_seq <-  seq(from=min(data_list_height2$count_palm_std) , to=max(data_list_height2$count_palm_std) , length=30)
# # generate preds to plot
# 
# p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*0 + post$bNxPTp*0*x ) )
# p_med <- apply( p_preds , 2 , median )
# lambda_preds <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*0 + post$bNxPTl*0*x ) )
# lambda_med <- apply( lambda_preds , 2 , median )
# joint_preds <- (1-p_preds)*(lambda_preds)
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/height_palm_maineff_ntu_m1fsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==0], 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm tree density/ 22 m^2 grid (standardized)" ,
#      ylab="height (m)" , cex.lab=1.5) 
# 
# #plot post median
# lines(joint_med ~ palm_seq , lw=2, col="purple1" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# p_preds2 <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*1 + post$bNxPTp*1*x ) )
# p_med2 <- apply( p_preds2 , 2 , median )
# lambda_preds2 <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*1 + post$bNxPTl*1*x ) )
# lambda_med2 <- apply( lambda_preds , 2 , median )
# joint_preds2 <- (1-p_preds2)*(lambda_preds2)
# joint_med2 <- apply( joint_preds2 , 2 , median )
# 
# png(file = "plots/height_palm_maineff_tu_m1fsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==1], 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm tree density/ 22 m^2 grid (standardized)" ,
#      ylab="height (m)" , cex.lab=1.5) 
# lines(joint_med2 ~ palm_seq , lw=2, col="orange2" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# ####prob
# p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*0 + post$bNxPTp*0*x ) )
# p_med <- apply( p_preds , 2 , median )
# joint_preds <- p_preds
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/prob_terr_palm_maineff_ntu_m1fsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==0], 
#      col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="palm tree density/ 22 m^2 grid (standardized)" ,
#      ylab="probability of terrestriality" , cex.lab=1.5) 
# 
# #plot post median
# lines(joint_med ~ palm_seq , lw=2, col="purple1" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# p_preds2 <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*1 + post$bNxPTp*1*x ) )
# p_med2 <- apply( p_preds2 , 2 , median )
# joint_preds2 <- p_preds2
# joint_med2 <- apply( joint_preds2 , 2 , median )
# 
# png(file = "plots/prob_height_palm_maineff_tu_m1fsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==1], 
#      col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="palm tree density/ 22 m^2 grid (standardized)" ,
#      ylab="probability of terrestriality" , cex.lab=1.5) 
# lines(joint_med2 ~ palm_seq , lw=2, col="orange2" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
# }
# dev.off()

########
# palm_seq <-  seq(from=min(data_list_height2$count_palm_std) , to=max(data_list_height2$count_palm_std) , length=30)
# # generate preds to plot
# 
# p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*0 + post$bNxPTp*0*x ) )
# p_med <- apply( p_preds , 2 , median )
# lambda_preds <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*0 + post$bNxPTl*0*x ) )
# lambda_med <- apply( lambda_preds , 2 , median )
# joint_preds <- (1-p_preds)*(lambda_preds)
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/height_palm_maineff_ntu_m1fsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==0], 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm tree density/ 22 m^2 grid (standardized)" ,
#      ylab="height (m)" , cex.lab=1.5) 
# 
# #plot post median
# lines(joint_med ~ palm_seq , lw=2, col="purple1" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# p_preds2 <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*1 + post$bNxPTp*1*x ) )
# p_med2 <- apply( p_preds2 , 2 , median )
# lambda_preds2 <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*1 + post$bNxPTl*1*x ) )
# lambda_med2 <- apply( lambda_preds , 2 , median )
# joint_preds2 <- (1-p_preds2)*(lambda_preds2)
# joint_med2 <- apply( joint_preds2 , 2 , median )
# 
# 
# 
# 
# png(file = "plots/height_palm_maineff_tu_m1fsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$count_palm_std[data_list_height2$nutcracking==1], 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm tree density/ 22 m^2 grid (standardized)" ,
#      ylab="height (m)" , cex.lab=1.5) 
# lines(joint_med2 ~ palm_seq , lw=2, col="orange2" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
# }
# dev.off()




##OLDER PLOTS

# png(file = "plots/terr_palm_maineff.png") 
# plot(data_list_height_zero$terr~data_list_height_zero$count_palm_std, 
#      col=col.alpha("black",0.01) , ylim=c(0,1) , xlab="palm density per 22 m square grid  (standardized)" ,
#      ylab="probability of being terrestrial" , pch=19) 
# 
# #plot post median
# lines(p_med ~ palm_seq , lw=2, col="black" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( p_preds[j,] ~ palm_seq , lw=3, col=col.alpha("black", alpha=0.05) , lty=1)
# }
# 
# lines(p_med2 ~ palm_seq , lw=2, col="purple1" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( p_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
# }
# dev.off()


##plot palm density as a function of height in follows with and without tool use
# logit(p) <- ap + a_id[id,1] + a_gp[grid_id_follow] + (bNp + a_id[id,3])*nutcracking 
# + (bPTp + a_id[id,5])*count_palm_std + (bNxPTp + a_id[id,7])*count_palm_std*nutcracking,
# 
# log(mu) <- al + a_id[id,2] + a_gl[grid_id_follow] + (bNl + a_id[id,4])*nutcracking 
# + (bPTl + a_id[id,6])*count_palm_std +  (bNxPTl + a_id[id,8])*count_palm_std*nutcracking,



# ##plot preds stones as a function of height
# palm_seq <-  seq(from=min(data_list_height2$count_palm_std) , to=max(data_list_height2$count_palm_std) , length=30)
# # generate preds to plot
# p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bSp*0 + post$bPTp*x) ) 
# p_med <- apply( p_preds , 2 , median )
# lambda_preds <- sapply( palm_seq, function(x) exp(  post$al + post$bSl*0 + post$bPTl*x) ) 
# lambda_med <- apply( lambda_preds , 2 , median )
# joint_preds <- (1-p_preds)*(lambda_preds)
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/height_palm_maineff.png") 
# plot(data_list_height2$height_m~data_list_height2$count_palm_std, 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm density per 22 m square grid  (standardized)" ,
#      ylab="height (m)") 
# 
# #plot post median
# lines(joint_med ~ palm_seq , lw=2, col="forestgreen" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("forestgreen", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# png(file = "plots/terr_palm_maineff.png") 
# plot(data_list_height2$terr~data_list_height2$count_palm_std, 
#      col=col.alpha("black",0.01) , ylim=c(0,1) , xlab="palm density per 22 m square grid  (standardized)" ,
#      ylab="probability of being terrestrial" , pch=19) 
# 
# #plot post median
# lines(p_med ~ palm_seq , lw=2, col="forestgreen" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( p_preds[j,] ~ palm_seq , lw=3, col=col.alpha("forestgreen", alpha=0.05) , lty=1)
# }
# dev.off()

#################m1d height in first point scan from tool use #####################
plot(precis(m1d))
post <- extract.samples(m1d)
str(post)
trout <- logistic(post$ap)
median(trout)
HPDI(trout)
pike <- logistic(post$ap + post$bNp)
median(pike)
HPDI(pike)
wamp <- precis(m1d , depth=2, digits=2)
write.csv ( wamp[1:4] , file="tables/m1d.csv")

png(file = "plots/prob_terr_first_scan_m1d.png",  width = 6, height = 6,units = "in" , res=300) 
    trout <- logistic(post$ap)
    par(mar=c(3,0.1,0.1,0.1))
    dens(trout , xlim=c(0, 1) , col="NA" , xlab="probability of being terrestrial in first scan", yaxt="n" , ylab="" , cex.lab=1.5)
    shade( density(trout) ,
           PCI(trout , prob=0.9999) ,
           col=col.alpha("purple1" , alpha=0.75) )
    
    pike <- logistic(post$ap + post$bNp)
    shade( density(pike) ,
           PCI(pike , prob=0.9999) ,
           col=col.alpha("orange2" , alpha=0.75) )
    
    blog=1.5
    for(i in 1:17){
        squanch <- logistic(post$ap + post$a_id[,i,1] )
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
        points( median(squanch) , (i-.1)/blog , pch=23 ,col="black", bg="purple1" , cex=0.75)
    }
    
    for(i in 1:17){
        squanch <- logistic(post$ap + post$a_id[,i,1]  + post$bNp + post$a_id[,i,2])
        segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
        points( median(squanch) , (i-.1)/blog , pch=23, col="black", bg="orange2", cex=0.75)
    }
    
    legend("topright" , legend=c("no tool-use" , "tool-use") , fill=c("purple1" , "orange2") , bty='n' )

dev.off()

trout <- (1- logistic(post$ap) )*exp(post$al)
median(trout)
HPDI(trout)
pike <- (1- logistic(post$ap + post$bNp) )*exp(post$al + post$bNl)
median(pike)
HPDI(pike)

png(file = "plots/height_first_scan_m1d.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(3,0.1,0.1,0.1))

trout <- (1- logistic(post$ap) )*exp(post$al)
dens(trout , xlim=c(0, 30) , col="NA" , xlab="height (m) at beginning of focal follow", yaxt="n" , ylab="" , cex.lab=1.5)
shade( density(trout) ,
       PCI(trout , prob=0.9999) ,
       col=col.alpha("purple1" , alpha=0.75) )

pike <- (1- logistic(post$ap + post$bNp) )*exp(post$al + post$bNl)
shade( density(pike) ,
       PCI(pike , prob=0.9999) ,
       col=col.alpha("orange2" , alpha=0.75) )

blog=30
for(i in 1:17){
    squanch <- (1- logistic(post$ap + post$a_id[,i,1]) )*exp(post$al + post$a_id[,i,3])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=23 ,col="black", bg="purple1" , cex=0.75)
}

for(i in 1:17){
    squanch <- (1- logistic(post$ap + post$a_id[,i,1] + post$bNp + post$a_id[,i,2]) )*exp(post$al + post$a_id[,i,3] + post$bNl + post$a_id[,i,4])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=23, col="black", bg="orange2", cex=0.75)
}

dev.off()

## height overall
median((1-logistic(post$ap))*exp(post$al))
HPDI ( (1-logistic(post$ap))*exp(post$al) )

png(file = "plots/height_overall.png",  width = 6, height = 6,units = "in" , res=300) 

trout <- (1-logistic(post$ap))*exp(post$al)
par(mar=c(3,0.1,0.1,0.1))
dens(trout , xlim=c(0, 20) , col="NA" , xlab="height (m)", yaxt="n" , ylab="" , cex.lab=1.5)
shade( density(trout) ,
       PCI(trout , prob=0.9999) ,
       col=col.alpha("purple1" , alpha=0.75) )
blog=20

for(i in 1:17){
    squanch <- (1-logistic(post$ap + post$a_id[,i,1] ))*exp(post$al+ post$a_id[,i,2])
    segments(x0=HPDI(squanch)[1] , x1=HPDI(squanch)[2] , y0=(i-.1)/blog , y1=(i-.1)/blog  )
    points( median(squanch) , (i-.1)/blog , pch=18 , cex=0.75 )
}
# points(data_list_height$height_m ,rep(-0.01, length(data_list_height$height_m)) ,
#        pch=19 , cex=0.5 , col=col.alpha("purple1" , alpha=0.05))

dev.off()

##############m1f#######
post <- extract.samples(m1f)
data_list_height_zero$terr <- ifelse(data_list_height_zero$height_m==0,1,0) 

palm_seq <-  seq(from=min(data_list_height_zero$count_palm_std) , to=max(data_list_height_zero$count_palm_std) , length=30)
# generate preds to plot

p_preds <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*0 + post$bNxPTp*0*x ) )
p_med <- apply( p_preds , 2 , median )
lambda_preds <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*0 + post$bNxPTl*0*x ) )
lambda_med <- apply( lambda_preds , 2 , median )
joint_preds <- (1-p_preds)*(lambda_preds)
joint_med <- apply( joint_preds , 2 , median )

#png(file = "plots/height_palm_maineff.png") 
plot(data_list_height_zero$height_m~data_list_height_zero$count_palm_std, 
     col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="palm density per 22 m square grid  (standardized)" ,
     ylab="height (m)") 

#plot post median
lines(joint_med ~ palm_seq , lw=2, col="black" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds[j,] ~ palm_seq , lw=3, col=col.alpha("black", alpha=0.05) , lty=1)
}
#dev.off()

p_preds2 <- sapply( palm_seq, function(x) logistic(  post$ap + post$bPTp*x + post$bNp*1 + post$bNxPTp*1*x ) )
p_med2 <- apply( p_preds2 , 2 , median )
lambda_preds2 <- sapply( palm_seq, function(x) exp(  post$al  + post$bPTl*x + post$bNl*1 + post$bNxPTl*1*x ) )
lambda_med2 <- apply( lambda_preds , 2 , median )
joint_preds2 <- (1-p_preds2)*(lambda_preds2)
joint_med2 <- apply( joint_preds2 , 2 , median )

lines(joint_med2 ~ palm_seq , lw=2, col="slateblue" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("slateblue", alpha=0.05) , lty=1)
}




#png(file = "plots/terr_palm_maineff.png") 
plot(data_list_height_zero$terr~data_list_height_zero$count_palm_std, 
     col=col.alpha("black",0.01) , ylim=c(0,1) , xlab="palm density per 22 m square grid  (standardized)" ,
     ylab="probability of being terrestrial" , pch=19) 

#plot post median
lines(p_med ~ palm_seq , lw=2, col="black" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( p_preds[j,] ~ palm_seq , lw=3, col=col.alpha("black", alpha=0.05) , lty=1)
}

lines(p_med2 ~ palm_seq , lw=2, col="slateblue" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( p_preds2[j,] ~ palm_seq , lw=3, col=col.alpha("slateblue", alpha=0.05) , lty=1)
}
#dev.off()


library(dagitty)
myDAG1 <- dagitty( 'dag {
Terrest -> ToolUse
PalmTrees -> ToolUse
Stone -> ToolUse
Predation -> ToolUse
Predation -> Terrest
Predation[unobserved]
}' )



plot(myDAG1)

#can't answer this
adjustmentSets(
    myDAG1,
    exposure = "Terrest",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

# adjustmentSets(
#     myDAG1,
#     exposure = "PalmTrees",
#     outcome = "ToolUse",
#     type =  "canonical",
#     effect = "total"
# )  

myDAG2 <- dagitty( 'dag {
PalmTrees -> ToolUse
Stone -> ToolUse
Predation -> ToolUse
Predation -> Terrest
Predation[unobserved]
}' )

plot(myDAG2)
adjustmentSets(
    myDAG2,
    exposure = "Terrest",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

myDAG3 <- dagitty( 'dag {
Terrest -> ToolUse
PalmTrees -> ToolUse
Stone -> ToolUse
Predation -> ToolUse
Predation -> Terrest
Predation[unobserved]
PalmTrees -> Terrest
Stone -> Terrest
}' )

plot(myDAG3)
adjustmentSets(
    myDAG3,
    exposure = "Terrest",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

adjustmentSets(
    myDAG3,
    exposure = "PalmTrees",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

adjustmentSets(
    myDAG3,
    exposure = "Stone",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

myDAG4 <- dagitty( 'dag {
PalmTrees -> ToolUse
Stone -> ToolUse
Predation -> ToolUse
Predation -> Terrest
Predation[unobserved]
PalmTrees -> Terrest
Stone -> Terrest
}' )

plot(myDAG4)
adjustmentSets(
    myDAG4,
    exposure = "Terrest",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

adjustmentSets(
    myDAG4,
    exposure = "PalmTrees",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 

adjustmentSets(
    myDAG4,
    exposure = "Stone",
    outcome = "ToolUse",
    type =  "canonical",
    effect = "total"
) 


#################m1gsimp#####################
# 
# post <- extract.samples(m1gsimp)
# precis(m1gsimp, pars="z_s" , depth=3)
# wamp <- precis(m1gsimp , depth=2, digits=2)
# write.csv ( wamp[1:4] , file="tables/m1gsimp.csv")
# 
# # data_list_daily_2$log_avg_stone_std
# 
# stone_seq <-  seq(from=min(data_list_height2$log_avg_stone_std) , to=max(data_list_height2$log_avg_stone_std) , length=30)
# # generate preds to plot
# p_preds <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*0 + post$bNxSp*0*x ) )
# p_med <- apply( p_preds , 2 , median )
# lambda_preds <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*0 + post$bNxSl*0*x ) )
# lambda_med <- apply( lambda_preds , 2 , median )
# joint_preds <- (1-p_preds)*(lambda_preds)
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/height_stone_maineff_ntu_m1gsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==0], 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="log stone density (standardized))" ,
#      ylab="height (m)" , cex.lab=1.5) 
# 
# #plot post median
# lines(joint_med ~ stone_seq , lw=2, col="purple1" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ stone_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# p_preds2 <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*1 + post$bNxSp*1*x ) )
# p_med2 <- apply( p_preds2 , 2 , median )
# lambda_preds2 <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*1 + post$bNxSl*1*x ) )
# lambda_med2 <- apply( lambda_preds , 2 , median )
# joint_preds2 <- (1-p_preds2)*(lambda_preds2)
# joint_med2 <- apply( joint_preds2 , 2 , median )
# 
# png(file = "plots/height_stone_maineff_tu_m1gsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==1], 
#      col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="log stone density (standardized)" ,
#      ylab="height (m)" , cex.lab=1.5) 
# lines(joint_med2 ~ stone_seq , lw=2, col="orange2" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds2[j,] ~ stone_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# 
# post <- extract.samples(m1gsimp)
# precis(m1gsimp, pars="z_s" , depth=3)
# blarg <- log(data_list_height2$stone_raw_grid) 
# max(blarg)
# precis(m1f, pars="z_id" , depth=2)
# 
# # data_list_daily_2$log_avg_stone_std
# 
# stone_seq <-  seq(from=min(data_list_height2$log_avg_stone_std) , to=max(data_list_height2$log_avg_stone_std) , length=30)
# # generate preds to plot
# p_preds <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*0 + post$bNxSp*0*x ) )
# p_med <- apply( p_preds , 2 , median )
# lambda_preds <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*0 + post$bNxSl*0*x ) )
# lambda_med <- apply( lambda_preds , 2 , median )
# joint_preds <- p_preds
# joint_med <- apply( joint_preds , 2 , median )
# 
# png(file = "plots/probterr_stone_maineff_ntu_m1gsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==0], 
#      col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="log stone density (standardized))" ,
#      ylab="probability of terrestriality" , cex.lab=1.5) 
# 
# #plot post median
# lines(joint_med ~ stone_seq , lw=2, col="purple1" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds[j,] ~ stone_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
# }
# dev.off()
# 
# p_preds2 <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*1 + post$bNxSp*1*x ) )
# p_med2 <- apply( p_preds2 , 2 , median )
# lambda_preds2 <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*1 + post$bNxSl*1*x ) )
# lambda_med2 <- apply( lambda_preds , 2 , median )
# joint_preds2 <- p_preds2
# joint_med2 <- apply( joint_preds2 , 2 , median )
# 
# png(file = "plots/probterr_stone_maineff_tu_m1gsimp.png",  width = 6, height = 6,units = "in" , res=300) 
# par(mar=c(4.5,4.5,0.25,0.25))
# plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==1], 
#      col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="log stone density (standardized)" ,
#      ylab="probability of terrestriality" , cex.lab=1.5) 
# lines(joint_med2 ~ stone_seq , lw=2, col="orange2" , lty=1)
# for (j in sample( c(1:2000) , 100) ){
#     lines( joint_preds2[j,] ~ stone_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
# }
# dev.off()

#####full m1g
post <- extract.samples(m1g)
wamp <- precis(m1g , depth=2, digits=2)
write.csv ( wamp[1:4] , file="tables/m1g.csv")

# data_list_daily_2$log_avg_stone_std

stone_seq <-  seq(from=min(data_list_height2$log_avg_stone_std) , to=max(data_list_height2$log_avg_stone_std) , length=30)
# generate preds to plot
p_preds <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*0 + post$bNxSp*0*x ) )
p_med <- apply( p_preds , 2 , median )
lambda_preds <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*0 + post$bNxSl*0*x ) )
lambda_med <- apply( lambda_preds , 2 , median )
joint_preds <- (1-p_preds)*(lambda_preds)
joint_med <- apply( joint_preds , 2 , median )

png(file = "plots/height_stone_maineff_ntu_m1g.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==0], 
     col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="log stone density (standardized))" ,
     ylab="height (m)" , cex.lab=1.5) 
text(-2.5,29,"b.")
#plot post median
lines(joint_med ~ stone_seq , lw=2, col="purple1" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds[j,] ~ stone_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
}
dev.off()

p_preds2 <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*1 + post$bNxSp*1*x ) )
p_med2 <- apply( p_preds2 , 2 , median )
lambda_preds2 <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*1 + post$bNxSl*1*x ) )
lambda_med2 <- apply( lambda_preds , 2 , median )
joint_preds2 <- (1-p_preds2)*(lambda_preds2)
joint_med2 <- apply( joint_preds2 , 2 , median )

png(file = "plots/height_stone_maineff_tu_m1g.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==1], 
     col=col.alpha("black",0.05) , ylim=c(0,30) , xlab="log stone density (standardized)" ,
     ylab="height (m)" , cex.lab=1.5) 
lines(joint_med2 ~ stone_seq , lw=2, col="orange2" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds2[j,] ~ stone_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
}
text(-2.5,29,"d.")
dev.off()

stone_seq <-  seq(from=min(data_list_height2$log_avg_stone_std) , to=max(data_list_height2$log_avg_stone_std) , length=30)
# generate preds to plot
p_preds <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*0 + post$bNxSp*0*x ) )
p_med <- apply( p_preds , 2 , median )
lambda_preds <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*0 + post$bNxSl*0*x ) )
lambda_med <- apply( lambda_preds , 2 , median )
joint_preds <- p_preds
joint_med <- apply( joint_preds , 2 , median )

png(file = "plots/probterr_stone_maineff_ntu_m1g.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==0] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==0], 
     col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="log stone density (standardized))" ,
     ylab="probability of terrestriality" , cex.lab=1.5) 

#plot post median
lines(joint_med ~ stone_seq , lw=2, col="purple1" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds[j,] ~ stone_seq , lw=3, col=col.alpha("purple1", alpha=0.05) , lty=1)
}
text(-2.5,1,"a.")

dev.off()

p_preds2 <- sapply( stone_seq, function(x) logistic(  post$ap + post$bSp*x + post$bNp*1 + post$bNxSp*1*x ) )
p_med2 <- apply( p_preds2 , 2 , median )
lambda_preds2 <- sapply( stone_seq, function(x) exp(  post$al  + post$bSl*x + post$bNl*1 + post$bNxSl*1*x ) )
lambda_med2 <- apply( lambda_preds , 2 , median )
joint_preds2 <- p_preds2
joint_med2 <- apply( joint_preds2 , 2 , median )

png(file = "plots/probterr_stone_maineff_tu_m1g.png",  width = 6, height = 6,units = "in" , res=300) 
par(mar=c(4.5,4.5,0.25,0.25))
plot(data_list_height2$height_m[data_list_height2$nutcracking==1] ~ data_list_height2$log_avg_stone_std[data_list_height2$nutcracking==1], 
     col=col.alpha("black",0.00) , ylim=c(0,1) , xlab="log stone density (standardized)" ,
     ylab="probability of terrestriality" , cex.lab=1.5) 
lines(joint_med2 ~ stone_seq , lw=2, col="orange2" , lty=1)
for (j in sample( c(1:2000) , 100) ){
    lines( joint_preds2[j,] ~ stone_seq , lw=3, col=col.alpha("orange2", alpha=0.05) , lty=1)
}
text(-2.5,1,"c.")
dev.off()
