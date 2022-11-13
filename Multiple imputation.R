
# Function for LOCF
coalesce <- function(...) {
  apply(cbind(...), 1, function(x) {
    x[which(!is.na(x))[1]]
  })
}

# Function for subset
subset_Trt <- function(data, value) {
  data1 <- data %>% filter(Trt==value)
  return(data1)
}


size <- 161 

rho1 <- 0.7
s1 <- 1
s2 <- 1
base <- 0.7
t_pbase <- 1.029 
c_pbase <- 0.77 
expected <- (((t_pbase - base)/base) - ((c_pbase - base)/base))*100
v1 <- NULL
samplemean <- data.frame(v1)
N <- 10
  simmn<-function(ss,N) #m1,m2,sigma,sigma1,n1,n2,N,miss_pcnt,imput_n)
  # for (ss in size)
{
for (i in  1:N)
{
  set.seed(123789+i)
  # test <- data.frame(cbind(rtmvnorm(n=800, mean=c(compeffect, treateffect), 
  #                                   sigma=matrix(c(s1, s1*s2*rho, s1*s2*rho, s2),2 )
  #                                   # , lower=c(-0.5, -0.3), upper=c(1.5, 1.7), algorithm="rejection"
  # ), rho))
  bvn1 <- data.frame(mvrnorm(ss, mu = c(base, t_pbase), Sigma = matrix(c(s1, s1*s2*rho1, s1*s2*rho1, s2),2 ))) # from MASS package
  colnames(bvn1) <- c("X1","X2")
  bvn1$pch <- ((bvn1$X2/bvn1$X1)-1)*100
  Median_t<-median(bvn1$pch);
  bvn1$rand10 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.8), 1)
  bvn1$rand15 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.7), 1)
  bvn1$rand20 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.6), 1)
  bvn1$rand25 <- ifelse(bvn1$pch < Median_t, rbinom(ss,1,0.5), 1)
  bvn1$pchmis10 <- ifelse (bvn1$pch < Median_t & bvn1$rand10 == 0,NA, bvn1$pch)
  bvn1$pchmis15 <- ifelse (bvn1$pch < Median_t & bvn1$rand15 == 0,NA, bvn1$pch)
  bvn1$pchmis20 <- ifelse (bvn1$pch < Median_t & bvn1$rand20 == 0,NA, bvn1$pch)
  bvn1$pchmis25 <- ifelse (bvn1$pch < Median_t & bvn1$rand25 == 0,NA, bvn1$pch)
  bvn1$X2_10 <- ifelse (bvn1$pch < Median_t & bvn1$rand10 == 0,NA, bvn1$X2)
  bvn1$X2_15 <- ifelse (bvn1$pch < Median_t & bvn1$rand15 == 0,NA, bvn1$X2)
  bvn1$X2_20 <- ifelse (bvn1$pch < Median_t & bvn1$rand20 == 0,NA, bvn1$X2)
  bvn1$X2_25 <- ifelse (bvn1$pch < Median_t & bvn1$rand25 == 0,NA, bvn1$X2)
  
  set.seed(123789+i)
  bvn2 <- data.frame(mvrnorm(ss, mu = c(base, c_pbase), Sigma = matrix(c(s1, s1*s2*rho1, s1*s2*rho1, s2),2 ))) # from MASS package
  colnames(bvn2) <- c("X1","X2")
  bvn2$pch <- ((bvn2$X2/bvn2$X1)-1)*100
  Median_c<-median(bvn2$pch);
  bvn2$rand10 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.8), 1)
  bvn2$rand15 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.7), 1)
  bvn2$rand20 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.6), 1)
  bvn2$rand25 <- ifelse(bvn2$pch < Median_c, rbinom(ss,1,0.5), 1)
  bvn2$pchmis10 <- ifelse (bvn2$pch < Median_c & bvn2$rand10 == 0,NA, bvn2$pch)
  bvn2$pchmis15 <- ifelse (bvn2$pch < Median_c & bvn2$rand15 == 0,NA, bvn2$pch)
  bvn2$pchmis20 <- ifelse (bvn2$pch < Median_c & bvn2$rand20 == 0,NA, bvn2$pch)
  bvn2$pchmis25 <- ifelse (bvn2$pch < Median_c & bvn2$rand25 == 0,NA, bvn2$pch)
  bvn2$X2_10 <- ifelse (bvn2$pch < Median_c & bvn2$rand10 == 0,NA, bvn2$X2)
  bvn2$X2_15 <- ifelse (bvn2$pch < Median_c & bvn2$rand15 == 0,NA, bvn2$X2)
  bvn2$X2_20 <- ifelse (bvn2$pch < Median_c & bvn2$rand20 == 0,NA, bvn2$X2)
  bvn2$X2_25 <- ifelse (bvn2$pch < Median_c & bvn2$rand25 == 0,NA, bvn2$X2)
  # sum(is.na(bvn1$pchmis25))
  
  mean <- cbind(
    # summary from comparator drug
    meant=mean(bvn1$pch), sqrt(var(bvn1$pch)), ss,mean(bvn1$X1),mean(bvn1$X2), 
    deltat=((mean(bvn1$X2)/mean(bvn1$X1))-1)*100, 
    meant_10 = mean(bvn1$pchmis10, na.rm=TRUE),
    meant_15 = mean(bvn1$pchmis15, na.rm=TRUE),
    meant_20 = mean(bvn1$pchmis20, na.rm=TRUE),
    meant_25 = mean(bvn1$pchmis25, na.rm=TRUE),
    # summary from comparator drug
    meanc=mean(bvn2$pch), sqrt(var(bvn2$pch)), mean(bvn2$X1),mean(bvn2$X2), 
    deltac=((mean(bvn2$X2)/mean(bvn2$X1))-1)*100, 
    meanc_10 = mean(bvn2$pchmis10, na.rm=TRUE),
    meanc_15 = mean(bvn2$pchmis15, na.rm=TRUE),
    meanc_20 = mean(bvn2$pchmis20, na.rm=TRUE),
    meanc_25 = mean(bvn2$pchmis25, na.rm=TRUE),
    #treatment difference
    effect=mean(bvn1$pch) - mean(bvn2$pch), 
    effect_10=mean(bvn1$pchmis10, na.rm=TRUE) - mean(bvn2$pchmis10, na.rm=TRUE),
    effect_15=mean(bvn1$pchmis15, na.rm=TRUE) - mean(bvn2$pchmis15, na.rm=TRUE),
    effect_20=mean(bvn1$pchmis20, na.rm=TRUE) - mean(bvn2$pchmis20, na.rm=TRUE),
    effect_25=mean(bvn1$pchmis25, na.rm=TRUE) - mean(bvn2$pchmis25, na.rm=TRUE),
    # Delta method
    d_effect=(((mean(bvn1$X2)/mean(bvn1$X1))-1) - ((mean(bvn2$X2)/mean(bvn2$X1))-1))*100,
    d_effect_10=(((mean(bvn1$X2_10, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_10, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
    d_effect_15=(((mean(bvn1$X2_15, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_15, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
    d_effect_20=(((mean(bvn1$X2_20, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_20, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
    d_effect_25=(((mean(bvn1$X2_25, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_25, na.rm=TRUE)/mean(bvn2$X1))-1))*100
  )
  samplemean <- rbind(samplemean, mean)
  # samplemean <- samplemean %>% slice(501:5000)
  
  predata0 <- rbind(cbind(bvn1 %>% dplyr::select(X1, X2, pch),Treat = 1),cbind(bvn2 %>% dplyr::select(X1, X2, pch),Treat = 0))
  colnames(predata0) <- c("baseline", "postbase", "response", "Trt")
  
  for (j in  c(20,30,40,50))
  {
    if (j == 20){predata <- rbind(cbind(bvn1 %>% dplyr::select(X1, pchmis10, X2_10),Treat = 1),cbind(bvn2 %>% dplyr::select(X1, pchmis10, X2_10),Treat = 0)) 
    colnames(predata) <- c("baseline", "response", "postbase", "Trt")}
    if (j == 30){predata <- rbind(cbind(bvn1 %>% dplyr::select(X1, pchmis15, X2_15),Treat = 1),cbind(bvn2 %>% dplyr::select(X1, pchmis15, X2_15),Treat = 0)) 
    colnames(predata) <- c("baseline", "response", "postbase", "Trt")}
    if (j == 40){predata <- rbind(cbind(bvn1 %>% dplyr::select(X1, pchmis20, X2_20),Treat = 1),cbind(bvn2 %>% dplyr::select(X1, pchmis20, X2_20),Treat = 0)) 
    colnames(predata) <- c("baseline", "response", "postbase", "Trt")}
    if (j == 50){predata <- rbind(cbind(bvn1 %>% dplyr::select(X1, pchmis25, X2_25),Treat = 1),cbind(bvn2 %>% dplyr::select(X1, pchmis25, X2_25),Treat = 0)) 
    colnames(predata) <- c("baseline", "response", "postbase", "Trt")}
    
    #multiple imputation
    imp <- mice(data = predata, m = 10, method = "pmm", maxit = 50, print=FALSE);
    # First, turn the datasets into long format
    imp_long <- mice::complete(imp, action="long", include = TRUE)
    # Convert treatment variable to factor
    imp_long$Trt <- with(imp_long, 
                         as.factor(imp_long$Trt))
    # Convert back to mids type - mice can work with this type
    imp_long_mids<-as.mids(imp_long)
    # Regression 
    fitimp <- with(imp_long_mids,
                   lm(response ~ Trt))
    pooled <- summary(pool(fitimp))
    delta_estimatet <- imp_long %>% filter(.imp > 0 & Trt==1) %>% group_by(factor(.imp)) %>% 
      summarize(post_t = mean(postbase),
                base_t = mean(baseline))
    delta_estimatec <- imp_long %>% filter(.imp > 0 & Trt==0) %>% group_by(factor(.imp)) %>% 
      summarize(post_c = mean(postbase),
                base_c = mean(baseline))
    
    if (i < 2) {
      if (j == 20) {flagdat1_20 <- cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                         ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                           (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                         )
      colnames(flagdat1_20) <- c("Estimate", "Std_Error", "D_estimate")}
      else if (j == 30) {flagdat1_30 <- cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                              ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                 (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                              )
      colnames(flagdat1_30) <- c("Estimate", "Std_Error", "D_estimate")}
      else if (j == 40) {flagdat1_40 <- cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                              ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                 (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                              )
      colnames(flagdat1_40) <- c("Estimate", "Std_Error", "D_estimate")}
      else if (j == 50) {flagdat1_50 <- cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                              ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                 (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                              )
      colnames(flagdat1_50) <- c("Estimate", "Std_Error", "D_estimate")}
      
      # if (j == 20) {flagdat1_20=summary(pooled)}
      # else if (j == 30) {flagdat1_30=summary(pooled)}
      # else if (j == 40) {flagdat1_40=summary(pooled)}
      # else if (j == 50) {flagdat1_50=summary(pooled)}
    } else
    {
      if (j == 20) {flagdat1_20 <- rbind(flagdat1_20,cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                                           ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                              (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                                           ))}
      else if (j == 30) {flagdat1_30 <- rbind(flagdat1_30,cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                                                ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                                   (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                                                ))}
      else if (j == 40) {flagdat1_40 <- rbind(flagdat1_40,cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                                                ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                                   (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                                                ))}
      else if (j == 50) {flagdat1_50 <- rbind(flagdat1_50,cbind(summary(pool(fitimp))[2,2], summary(pool(fitimp))[2,3],
                                                                ((mean(delta_estimatet$post_t)/mean(delta_estimatet$base_t))-
                                                                   (mean(delta_estimatec$post_c)/mean(delta_estimatec$base_c)))*100
                                                                ))}
      
      # if (j == 20) {flagdat1_20=rbind(flagdat1_20,summary(pooled))}
      # else if (j == 30) {flagdat1_30=rbind(flagdat1_30,summary(pooled))}
      # else if (j == 40) {flagdat1_40=rbind(flagdat1_40,summary(pooled))}
      # else if (j == 50) {flagdat1_50=rbind(flagdat1_50,summary(pooled))}
    }
    
    #completers set
    predata1 <- predata[!(predata$response == "NA"),]
    fit2 <- aov(response~Trt, predata1)
    
    if (i < 2) {
      if (j == 20) {flagdat2_20 <- cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2])
                    colnames(flagdat2_20) <- c("Estimate", "Std_Error")}
      else if (j == 30) {flagdat2_30 <- cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2])
      colnames(flagdat2_30) <- c("Estimate", "Std_Error")}
      else if (j == 40) {flagdat2_40 <- cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2])
      colnames(flagdat2_40) <- c("Estimate", "Std_Error")}
      else if (j == 50) {flagdat2_50 <- cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2])
      colnames(flagdat2_50) <- c("Estimate", "Std_Error")}
      
      # if (j == 20) {flagdat2_20 =coef(summary.lm(fit2))}
      # else if (j == 30) {flagdat2_30=coef(summary.lm(fit2))}
      # else if (j == 40) {flagdat2_40=coef(summary.lm(fit2))}
      # else if (j == 50) {flagdat2_50=coef(summary.lm(fit2))}
    } else
    {
      if (j == 20) {flagdat2_20 <- rbind(flagdat2_20,cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2]))}
      else if (j == 30) {flagdat2_30 <- rbind(flagdat2_30,cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2]))}
      else if (j == 40) {flagdat2_40 <- rbind(flagdat2_40,cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2]))}
      else if (j == 50) {flagdat2_50 <- rbind(flagdat2_50,cbind(coef(summary.lm(fit2))[2,1], coef(summary.lm(fit2))[2,2]))}
      
      # d_effect_10=(((mean(bvn1$X2_10, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_10, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
      # d_effect_15=(((mean(bvn1$X2_15, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_15, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
      # d_effect_20=(((mean(bvn1$X2_20, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_20, na.rm=TRUE)/mean(bvn2$X1))-1))*100,
      # d_effect_25=(((mean(bvn1$X2_25, na.rm=TRUE)/mean(bvn1$X1))-1) - ((mean(bvn2$X2_25, na.rm=TRUE)/mean(bvn2$X1))-1))*100
      
      # if (j == 20) {flagdat2_20=rbind(flagdat2_20,coef(summary.lm(fit2)))}
      # else if (j == 30) {flagdat2_30=rbind(flagdat2_30,coef(summary.lm(fit2)))}
      # else if (j == 40) {flagdat2_40=rbind(flagdat2_40,coef(summary.lm(fit2)))}
      # else if (j == 50) {flagdat2_50=rbind(flagdat2_50,coef(summary.lm(fit2)))}
    }
    
    #LOCF
    predata2 <- predata
    predata2$response1 <- coalesce(predata$response, 0)
    predata2$postbase1 <- coalesce(predata$postbase, baseline)
    fit3 <- aov(response1~Trt, predata2)
    
    if (i < 2) {
      # if (j == 20) {flagdat3_20=coef(summary.lm(fit3))}
      # else if (j == 30) {flagdat3_30=coef(summary.lm(fit3))}
      # else if (j == 40) {flagdat3_40=coef(summary.lm(fit3))}
      # else if (j == 50) {flagdat3_50=coef(summary.lm(fit3))}
      if (j == 20) {flagdat3_20 <- cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2],
                                         ((mean(subset_Trt(predata2,1)$postbase1)/mean(subset_Trt(predata2,1)$baseline)) - 
                                            (mean(subset_Trt(predata2,0)$postbase1)/mean(subset_Trt(predata2,0)$baseline)))*100
                                         )
                    colnames(flagdat3_20) <- c("Estimate", "Std_Error","D_estimate")
                    }
      else if (j == 30) {flagdat3_30 <- cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2])
                        colnames(flagdat3_30) <- c("Estimate", "Std_Error","D_estimate")
                        }
      else if (j == 40) {flagdat3_40 <- cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2])
                        colnames(flagdat3_40) <- c("Estimate", "Std_Error","D_estimate")
                        }
      else if (j == 50) {flagdat3_50 <- cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2])
                        colnames(flagdat3_50) <- c("Estimate", "Std_Error","D_estimate")
                        }
    } else
    {
      # if (j == 20) {flagdat3_20=rbind(flagdat3_20,coef(summary.lm(fit3)))}
      # else if (j == 30) {flagdat3_30=rbind(flagdat3_30,coef(summary.lm(fit3)))}
      # else if (j == 40) {flagdat3_40=rbind(flagdat3_40,coef(summary.lm(fit3)))}
      # else if (j == 50) {flagdat3_50=rbind(flagdat3_50,coef(summary.lm(fit3)))}
      if (j == 20) {flagdat3_20 <- rbind(flagdat3_20,cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2]))}
      else if (j == 30) {flagdat3_30 <- rbind(flagdat3_30,cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2]))}
      else if (j == 40) {flagdat3_40 <- rbind(flagdat3_40,cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2]))}
      else if (j == 50) {flagdat3_50 <- rbind(flagdat3_50,cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2]))}
    }
  }
  
  #whole data
  fit1 <- aov(response~Trt, predata0)
  if (i < 2) {flagdat <- cbind(coef(summary.lm(fit1))[2,1], 
                               coef(summary.lm(fit1))[2,2],
                               ((mean(subset_Trt(predata0,1)$postbase)/mean(subset_Trt(predata0,1)$baseline)) - 
                                  (mean(subset_Trt(predata0,0)$postbase)/mean(subset_Trt(predata0,0)$baseline)))*100
                               )
              colnames(flagdat) <- c("Estimate", "Std_Error","D_estimate")} else
      # flagdat=coef(summary.lm(fit1))} else
  {
    flagdat=rbind(flagdat,cbind(coef(summary.lm(fit3))[2,1], coef(summary.lm(fit3))[2,2],
                                ((mean(subset_Trt(predata0,1)$postbase)/mean(subset_Trt(predata0,1)$baseline)) - 
                                   (mean(subset_Trt(predata0,0)$postbase)/mean(subset_Trt(predata0,0)$baseline)))*100
                                ))}
    # flagdat=rbind(flagdat,coef(summary.lm(fit1)))}
}

# return(list(flagdat))
  return(list(#Actual complete data estimates
              Orig_mean=mean(data.frame(flagdat)$Estimate),
              Orig_stderr=sqrt(mean(data.frame(flagdat)$Std_Error**2)),
              Orig_dmean=mean(data.frame(flagdat)$D_estimate),
              Orig_dstderr=sqrt(mean(data.frame(flagdat)$Std_Error**2)),
              #Completers set estimates ignoring the missing
              COMP_10_mean=mean(data.frame(flagdat2_20)$Estimate),
              COMP_10_stderr=sqrt(mean(data.frame(flagdat2_20)$Std_Error**2)),
              COMP_15_mean=mean(data.frame(flagdat2_30)$Estimate),
              COMP_15_stderr=sqrt(mean(data.frame(flagdat2_30)$Std_Error**2)),
              COMP_20_mean=mean(data.frame(flagdat2_40)$Estimate),
              COMP_20_stderr=sqrt(mean(data.frame(flagdat2_40)$Std_Error**2)),
              COMP_25_mean=mean(data.frame(flagdat2_50)$Estimate),
              COMP_25_stderr=sqrt(mean(data.frame(flagdat2_50)$Std_Error**2)),
              #Estimates imputing the missing with LOCF method
              LOCF_10_mean=mean(data.frame(flagdat3_20)$Estimate),
              LOCF_10_stderr=sqrt(mean(data.frame(flagdat3_20)$Std_Error**2)),
              LOCF_15_mean=mean(data.frame(flagdat3_30)$Estimate),
              LOCF_15_stderr=sqrt(mean(data.frame(flagdat3_30)$Std_Error**2)),
              LOCF_20_mean=mean(data.frame(flagdat3_40)$Estimate),
              LOCF_20_stderr=sqrt(mean(data.frame(flagdat3_40)$Std_Error**2)),
              LOCF_25_mean=mean(data.frame(flagdat3_50)$Estimate),
              LOCF_25_stderr=sqrt(mean(data.frame(flagdat3_50)$Std_Error**2)),
              #Estimates imputing the missing with MI method
              MI_10_mean=mean(data.frame(flagdat1_20)$Estimate),
              MI_10_stderr=sqrt(mean(data.frame(flagdat1_20)$Std_Error**2)),
              MI_10_dmean=mean(data.frame(flagdat1_20)$D_estimate),
              MI_15_mean=mean(data.frame(flagdat1_30)$Estimate),
              MI_15_stderr=sqrt(mean(data.frame(flagdat1_30)$Std_Error**2)),
              MI_15_dmean=mean(data.frame(flagdat1_30)$D_estimate),
              MI_20_mean=mean(data.frame(flagdat1_40)$Estimate),
              MI_20_stderr=sqrt(mean(data.frame(flagdat1_40)$Std_Error**2)),
              MI_20_dmean=mean(data.frame(flagdat1_40)$D_estimate),
              MI_25_mean=mean(data.frame(flagdat1_50)$Estimate),
              MI_25_stderr=sqrt(mean(data.frame(flagdat1_50)$Std_Error**2)),
              MI_25_dmean=mean(data.frame(flagdat1_50)$D_estimate)
              )
         )
}
 out.matrix <- data.frame(simmn(ss=size,N=10)) #mu7,mu8,sigma7,sigma8,791,791,1,20,10)
 out.matrix