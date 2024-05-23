########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################################################
# Function to do the population level estimations for Aim 1, i.e., model 1/model 2/ model 3


#################################################################################################################
# LM model / t-test between ASD & TD
# Model 1: blinkcoef_iv1 ~ \alpha_v1 + \gamma_v1 * ASD + \nu_v1
# Model 2: fixcoef_iv2 ~ \alpha_v2 + \gamma_v2 * ASD + \nu_v2
# Model 3: (fixcoef_iv2 - blinkcoef_iv1) ~ (\alpha_v2 - \alpha_v1) + (\gamma_v2 - \gamma_v1) + (\nu_v2 - \nu_v1)

#######################
####### Model 1 #######
#######################
# for loop to do t-test for each region
LM.pop.model <- function(combined_data_frame.asdtd, outcome, model3 = FALSE) {
    # Initialize two dataframe
    ASD.TD.lm.nonASD <- data.frame()
    ASD.TD.lm.ASD <- data.frame()
    if (model3) {
      for (i in 1:100) {
      subset.dat <- subset(combined_data_frame.asdtd, Region == i)
      test <- lm( (subset.dat[["fixcoef"]] - subset.dat[["blinkcoef"]]) ~ ASD, data = subset.dat)
      pvalue.nonASD <- summary(test)$coefficient[1,] ##!!!!!### 1 is beta0, i.e., non.ASD only. 2 is beta1, the difference between non.ASD and ASD 
      pvalue.ASD <- summary(test)$coefficient[2,] 
      ASD.TD.lm.nonASD <- rbind(ASD.TD.lm.nonASD, pvalue.nonASD)
      colnames(ASD.TD.lm.nonASD) <- c("coefs", "se", "tstat","pval")
      ASD.TD.lm.ASD <- rbind(ASD.TD.lm.ASD, pvalue.ASD)
      colnames(ASD.TD.lm.ASD) <- c("coefs", "se", "tstat","pval") 
      }
    } else {
      for (i in 1:100) {
        subset.dat <- subset(combined_data_frame.asdtd, Region == i)
        test <- lm( (subset.dat[[outcome]]) ~ ASD, data = subset.dat)
        pvalue.nonASD <- summary(test)$coefficient[1,] ##!!!!!### 1 is beta0, i.e., non.ASD only. 2 is beta1, the difference between non.ASD and ASD 
        pvalue.ASD <- summary(test)$coefficient[2,] 
        ASD.TD.lm.nonASD <- rbind(ASD.TD.lm.nonASD, pvalue.nonASD)
        colnames(ASD.TD.lm.nonASD) <- c("coefs", "se", "tstat","pval")
        ASD.TD.lm.ASD <- rbind(ASD.TD.lm.ASD, pvalue.ASD)
        colnames(ASD.TD.lm.ASD) <- c("coefs", "se", "tstat","pval") 
    }
    }
     return(list(ASD.TD.lm.nonASD, ASD.TD.lm.ASD))
}

########################################################################################################################################################################################################################################################################################################################
# Making brain mapping plots
# Function for coefficients
Brainmap.coefs <- function( mappingdat = ASD.TD.lm.nonASD[,1], xii = xii) {
  # mappingdat: V length data to be mapped
  brain.plot <- fMRIBrain_Mapping(dtseries_data = xii, mapping_data = mappingdat)
  return(brain.plot)
}

# Function for FDR pvalue
Brainmap.pval <- function(mappingdat = ASD.TD.lm.nonASD[,4], xii = xii, fdr.threshold = 0.05) {
  # do FDR adjustment
  fdr.pvalue <- p.adjust(mappingdat, method = "fdr")
  # set the pvals beyond the threshold to NAs
  fdr.pvalue[which(fdr.pvalue >= fdr.threshold)] <- NA
  # plot -log10(fdr.pval)
  inputformapping.fdr <- -log10(fdr.pvalue)
  brain.plot.fdr <- fMRIBrain_Mapping(dtseries_data = xii, mapping_data = inputformapping.fdr)
  return(brain.plot.fdr)
}



