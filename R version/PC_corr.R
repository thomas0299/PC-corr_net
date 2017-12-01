###### PC-corr algorithm associate to PCA analysis #####

### Released under MIT License
### Copyright (c) 2017 Sara Ciucci, Yan Ge, Claudio Durán and Carlo Vittorio Cannistraci

# Please cite:
# Enlightening discriminative network functional modules behind Principal Component Analysis separation in differential-omic science studies.
# Sara Ciucci, Yan Ge, Claudio Durán, Alessandra Palladini, Víctor Jiménez Jiménez, Luisa María Martínez Sánchez, 
# Yuting Wang, Susanne Sales, Andrej Shevchenko, Steven W. Poser, Maik Herbig, Oliver Otto, Andreas Androutsellis-Theotokis, 
# Jochen Guck, Mathias J. Gerl and Carlo Vittorio Cannistraci 
# Scientific Reports, 2017

# INPUT
#   x => (Numeric matrix MxN) Dataset with samples on the rows and features on the columns
#   sample_labels => (Character vector Mx1) Labels of the samples
#   feat_names => (Character vector Nx1) Names of the features
#   sample_names => (Character vector Mx1) Names of the samples
#   dis =>'yes' or 'no', depending if you want to display the sample names
#   in the scatterplot. Default: 'no'

# OUTPUT
#
#   -For a single cut-off: 
#
#   Edges => (Data frame Ox3) Edges in the PC-corr constructed network: 
#          -first/second columns: two nodes (features) that are connected by an edge 
#          -third column: PC-corr edge weight.
#   Nodes => (Data frame Px3) Nodes in the constructed network: 
#           -first column: nodes (features) in the network
#           -second column: corresponding node colors
#           -third column: correspondind loading value (V)
#
#   -For multiple cut-offs (say C cut-offs):
#
#   Edges => (Matrix of lists Cx2)
#            -first column: cut-offs
#            -second column: respective edge table in the PC-corr network for each cut-off 
#            (that is a data frame with the same structure as in the single cut-off case)
#   Nodes => (Matrix of lists Cx2)
#            -first column: cut-offs
#            -second column: respective node table in the PC-corr network  for each cut-off 
#            (that is a data frame with the same structure as in the single cut-off case)

PC_corr_19_5_2017<-function(x,sample_labels,feat_names, sample_names,dis) {
  
  # initialisation and default options --------------------------------------
  
  library(base)
  library(stats)
  library(xlsx)
  library(knitr)
  library(preprocessCore)
  library(ggplot2)
  library(gridExtra)
  library(pROC)
  library(ROCR)
  library(caTools)
  library(ggrepel)
  
  n<-nargs()
  if (n<4){
    stop('\nNot Enough Input Arguments')
  }
  
  if (n<5){
    dis<-'no';
  }
  
  
  if (sum(is.nan(x))==1) {
    stop(sprintf('\nThere is %d NaN value in your data matrix.\nPlease replace it.',sum(is.nan(x))))
  } else if (sum(is.nan(x))>1) {
    stop(sprintf('\nThere are %d NaN values in your data matrix. \nPlease replace them.',sum(is.nan(x))))
  }
  
  
  #  ------------------------------------------------------------------------
  
  labels <- sample_labels;
  nameLabels <- unique(labels) # name of the groups
  numbLabels <- length(nameLabels) # number of groups
  
  # Remove features with same identical values across all the samples
  x1 <- x 
  x1 <- x1-matrix(rep(colMeans(x1),each=nrow(x1)),nrow=dim(x1)[1],ncol=dim(x1)[2])
  ma <- colSums(x1==0)==nrow(x1)
  remov_feat <- sum(ma) #number of removed fatures
  x <- x[,!ma]
  feat_names<-feat_names[t(!ma)]
  
  
  # Normalizations of the dataset ------------------------------------------- 
  
  norm <- matrix(list(),12,1)
  norms <-  character()
  
  norm[[1]] <- x/matrix(rep(colSums(x),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2]) #dividing by the sum over the samples
  norms[1] <- 'DSOR'
  
  norm[[2]] <- x/matrix(rep(rowSums(x),each = ncol(x)),nrow = dim(x)[1],ncol = dim(x)[2],byrow=TRUE) #dividing by the sum over the features
  norms[2] <- 'DSOC'
  
  norm[[3]] <- log10(1+x)
  norms[3] <- 'LOG'
  
  norm[[4]] <- scale(x)
  norms[4] <- 'ZSCORE'
  
  norm[[5]] <- t(normalize.quantiles(t(x)))
  norms[5] <- 'QUANTILE T'
  
  norm[[6]] <- normalize.quantiles(x)
  norms[6] <- 'QUANTILE'
  
  norm[[7]] <- t(scale(t(x)))
  norms[7] <- 'ZSCORE T'
  
  norm[[8]] <- x + abs(min(x))
  norms[8] <- 'PLUS(ABS(MIN))'
  
  x_center <- x-matrix(rep(colMeans(x),each=nrow(x)),nrow=dim(x)[1],ncol=dim(x)[2])
  norm[[9]] <- x_center/matrix(rep(sqrt(apply(x_center, 2, sd)),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2])
  norms[9] <- 'PARETO SCALING'
  
  norm[[10]] <- sqrt(x)
  norms[10] <- 'SQRT'
  
  norm[[11]] <- x/matrix(rep(colMeans(x),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2])
  norms[11] <- 'MANORM'
  
  norm[[12]] <- x  #No normalization
  norms[12] <- '-' 
  
  #  ------------------------------------------------------------------------
  u_lab <- c( )
  flag <- 0
  while (flag == 0){
    cat("\nIs your data represented by\n[r]ranked labels(labels that are organized according to a progressive order. e.g. different stages of a disease, where Stage 1 < Stage 2 < Stage 3)\n[c]class labels (labels that are not necessary organized in a progressive order e.g. Condition A, Condition B, Condition C)? [r/c]:\n\n")
    u_lab <-readline(prompt="-> ")
    flag <- 1
    if ((u_lab != 'r') & (u_lab != 'c')){
      flag <- 0
      cat("Please introduce either 'r' for ranked labels or c' for class labels.\n")
    }
  }
  
  if (u_lab =='r') {
    flag <- 0 
    while (flag == 0){
      cat("\nAre the values of your ranked labels\n[d] discrete (Stage 1 < Stage 2 < Stage 3)\nor [con] continuous (different times of development of a cell line)? [d/con]:\n\n")
      u_lab <- readline(prompt="-> ")
      flag <- 1
      if ((u_lab != 'd') & (u_lab != 'con')){
        flag <- 0
        cat("Please introduce either 'd' for discrete labels or 'con' for continuous labels.\n")
      }
    }
  }
  
  labl_numb <- c()
  if (u_lab =='d') {
    #for correlation evaluators, labels are turned into numbers
    for (i in 1:numbLabels){
      labl_numb[labels %in% nameLabels[i]] <- i
    } 
    } else if (u_lab == 'con'){
      labl_numb <- as.numeric(labels)
    }

  
  #  ------------------------------------------------------------------------
  pc_nc<- matrix(list(),12,1)
  pc_c<- matrix(list(),12,1)
  ncPCA<- pc_nc
  cPCA<- pc_c
  explained_nc <- matrix(list(),12,1)
  explained_c <- matrix(list(),12,1)
  mw_ncPCA <- array(NA,dim=c(length(norm),choose(numbLabels,2),length(labels)))
  mw_cPCA <- mw_ncPCA
  AUC_nc <- array(NA,dim=c(length(norm),choose(numbLabels,2),length(labels)))
  AUC_c <- AUC_nc
  response <- c()
  AUPR_nc <- array(NA,dim=c(length(norm),choose(numbLabels,2),length(labels)))
  AUPR_c <- AUPR_nc
  rank_pears_corr_ncPCA <- array(NA,dim=c(length(norm),length(labels),1))
  rank_pears_corr_cPCA <- rank_pears_corr_ncPCA
  rank_spear_corr_ncPCA <- array(NA,dim=c(length(norm),length(labels),1))
  rank_spear_corr_cPCA <- rank_spear_corr_ncPCA

  flag <- 0 
  for (i in 1:length(norm)){
    
    ###non-centred PCA
    res_svd_nc <- svd(norm[[i]])
    snc <- res_svd_nc$d
    pc_nc[[i]] <- res_svd_nc$v 
    ncPCA[[i]]<- norm[[i]]%*%pc_nc[[i]] 
    
    ###centred PCA
    normCenter <- norm[[i]]-matrix(rep(colMeans(norm[[i]]),each=nrow(norm[[i]])),nrow=dim(norm[[i]])[1],ncol=dim(norm[[i]])[2])
    res_svd_c <- svd(normCenter)
    sc <- res_svd_c$d
    pc_c[[i]] <- res_svd_c$v 
    cPCA[[i]]<- normCenter%*%pc_c[[i]] 
    
    latent_nc <-(snc^2)/dim(x)[1] # component variance 
    explained_nc[[i]] <- 100*latent_nc/sum(latent_nc) # explained variance (%)
    
    latent_c <-(sc^2)/(dim(x)[1]-1) # component variance 
    explained_c[[i]] <- 100*latent_c/sum(latent_c) # explained variance (%)
    
    for (k in  1:dim(ncPCA[[i]])[2]){ #dimension
      
      if ((u_lab == 'c')|(u_lab == 'd')){
        n <- 1
        m <- 2
        for (j in 1:choose(numbLabels,2)){ # two-group comparison
          # Compute p-value of Mann-Whitney test
          mw_ncPCA[i,j,k] <- wilcox.test(ncPCA[[i]][is.element(labels,nameLabels[n]),k],ncPCA[[i]][is.element(labels,nameLabels[m]),k],paired=FALSE)$p.value
          mw_cPCA[i,j,k] <- wilcox.test(cPCA[[i]][is.element(labels,nameLabels[n]),k],cPCA[[i]][is.element(labels,nameLabels[m]),k],paired=FALSE)$p.value
          
          samp_lab <- c(sample_labels[is.element(labels,nameLabels[n])],sample_labels[is.element(labels,nameLabels[m])])
          scores_nc <- c(ncPCA[[i]][is.element(labels,nameLabels[n]),k],ncPCA[[i]][is.element(labels,nameLabels[m]),k])
          scores_c <- c(cPCA[[i]][is.element(labels,nameLabels[n]),k],cPCA[[i]][is.element(labels,nameLabels[m]),k])
          
          response <- c()
          # Compute AUC
          response[samp_lab == nameLabels[n]] <- 1
          response[samp_lab == nameLabels[m]] <- 0
          
          res_ROC_nc <- roc(response,scores_nc,direction = "<")
          res_ROC_c <- roc(response,scores_c,direction = "<")
          AUC_nc[i,j,k] <- res_ROC_nc$auc[1]
          AUC_c[i,j,k] <- res_ROC_c$auc[1]
          
          
          if (AUC_nc[i,j,k] < 0.5){
            flag <- 1
            AUC_nc[i,j,k] <- 1-AUC_nc[i,j,k]
            #Compute AUPR
            pred_nc<- prediction(scores_nc,samp_lab,label.ordering=c(nameLabels[n],nameLabels[m]))
            RP.perf_nc <- performance(pred_nc, "prec", "rec")
            AUPR_nc[i,j,k] <- trapz(unlist(RP.perf_nc@x.values)[2:length(unlist(RP.perf_nc@x.values))], unlist(RP.perf_nc@y.values)[2:length(unlist(RP.perf_nc@y.values))])
          } else {
            #Compute AUPR
            pred_nc<- prediction(scores_nc,samp_lab,label.ordering=c(nameLabels[m],nameLabels[n]))
            RP.perf_nc <- performance(pred_nc, "prec", "rec")
            AUPR_nc[i,j,k] <- trapz(unlist(RP.perf_nc@x.values)[2:length(unlist(RP.perf_nc@x.values))], unlist(RP.perf_nc@y.values)[2:length(unlist(RP.perf_nc@y.values))])
          }

          
          if (AUC_c[i,j,k] < 0.5){
            flag <- 1
            AUC_c[i,j,k] <- 1-AUC_c[i,j,k]
            #Compute AUPR
            pred_c<- prediction(scores_c,samp_lab,label.ordering=c(nameLabels[n],nameLabels[m]))
            RP.perf_c <- performance(pred_c, "prec", "rec")
            AUPR_c[i,j,k] <- trapz(unlist(RP.perf_c@x.values)[2:length(unlist(RP.perf_c@x.values))], unlist(RP.perf_c@y.values)[2:length(unlist(RP.perf_c@y.values))])
          } else {
            # Compute AUPR
            pred_c <- prediction(scores_c,samp_lab, label.ordering=c(nameLabels[m],nameLabels[n]))
            RP.perf_c <- performance(pred_c, "prec", "rec")
            AUPR_c[i,j,k] <- trapz(unlist(RP.perf_c@x.values)[2:length(unlist(RP.perf_c@x.values))], unlist(RP.perf_c@y.values)[2:length(unlist(RP.perf_c@y.values))])
          }
          
          

          m <- m + 1
          if (m > numbLabels){
            n <- n + 1
            m <- n + 1
          }
        }
        
        if (u_lab=='d'){
          rank_pears_corr_ncPCA[i,k,1] <- cor(ncPCA[[i]][,k],labl_numb, method="pearson")
          rank_pears_corr_cPCA[i,k,1] <- cor(cPCA[[i]][,k],labl_numb, method="pearson")
          rank_spear_corr_ncPCA[i,k,1] <- cor(ncPCA[[i]][,k],labl_numb, method="spearman")
          rank_spear_corr_cPCA[i,k,1] <- cor(cPCA[[i]][,k],labl_numb, method="spearman")
        }
      } else {
        

        rank_pears_corr_ncPCA[i,k,1] <- cor(ncPCA[[i]][,k],labl_numb, method="pearson")
        rank_pears_corr_cPCA[i,k,1] <- cor(cPCA[[i]][,k],labl_numb, method="pearson")
        rank_spear_corr_ncPCA[i,k,1] <- cor(ncPCA[[i]][,k],labl_numb, method="spearman")
        rank_spear_corr_cPCA[i,k,1] <- cor(cPCA[[i]][,k],labl_numb, method="spearman")
      }
    }
  }
  
  if (flag == 1){
    warning('Since some AUC-ROC were < 0.5, we inverted their labels and now AUC-ROC are >=0.5')
  }
  
  # Constructing the table of results ---------------------------------------
  norms <- rep(norms,times=length(labels)*2)
  centred <- rep('yes',times=length(labels)*12)
  non_centred <- rep('no',times=length(labels)*12)
  centr <- c(non_centred,centred)
  dim <- 1:length(labels)
  dim <- rep(dim, each=12)
  dim <- c(dim,dim)
  variance_c <- matrix(unlist(explained_c),nrow = 12, ncol = length(labels) , byrow=TRUE)
  variance_c <- as.vector(variance_c)
  variance_nc <- matrix(unlist(explained_nc),nrow = 12, ncol = length(labels) , byrow=TRUE)
  variance_nc <- as.vector(variance_nc)
  variance <- c(variance_nc,variance_c)
  
  if (numbLabels==2){
    if (u_lab=='c'){
      header <- c("P-value","AUC","AUPR","Norm","Centering","Dim","expl Var")
      pvals <- c(mw_ncPCA,mw_cPCA)
      AUCs <- c(AUC_nc,AUC_c)
      AUPRs <- c(AUPR_nc,AUPR_c)
      results <- data.frame(pvals,AUCs,AUPRs,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
      
    } else if (u_lab=='d'){
      header <- c("P-value","AUC","AUPR", "pears","spear","Norm","Centering","Dim","expl Var")
      pvals <- c(mw_ncPCA,mw_cPCA)
      AUCs <- c(AUC_nc,AUC_c)
      AUPRs <- c(AUPR_nc,AUPR_c)
      pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
      spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
      results <- data.frame(pvals,AUCs,AUPRs,pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
      
    } else if (u_lab=='con'){
      header <- c("pears","spear","Norm","Centering","Dim","expl Var")
      pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
      spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
      results <- data.frame(pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
      
    }
    
    flag <- 0
    while (flag==0){
      if ((u_lab=='c')|(u_lab =='d')){
        if (u_lab=='c'){
          cat("\nWould you like to rank the PCA results by P-value, AUC or AUPR [p/auc/aupr]:\n\n")
          u_rank <- readline(prompt="-> ")
        } else {
          cat("\nWould you like to rank the PCA results by P-value, AUC, AUPR, Pearson correlation or Spearman correlation [p/auc/aupr/pc/sc]:\n\n")
          u_rank <- readline(prompt="-> ")
        }
        
        flag <- 1
        
        if (u_rank=='p'){
          
          results <- results[order(results[,1],decreasing=FALSE),]
          #Only some results are shown on the screen
          if (nrow(results[results[,1]<0.05,])==0){
            print(kable(results[results[,1]<0.25,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[results[,1]<0.05,], format="markdown",row.names=FALSE))
          }
          
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
          
        } else if (u_rank =='auc'){
          
          results <- results[order(results[,2],decreasing=TRUE),]
          #Only some results are shown on the screen
          if (nrow(results[results[,2]>=0.7,])==0){
            print(kable(results[results[,2]>=0.6,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[results[,2]>=0.7,], format="markdown",row.names=FALSE))
          }
          
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank =='aupr'){
          
          results <- results[order(results[,3],decreasing=TRUE),]
          #Only some results are shown on the screen
          if (nrow(results[results[,3]>=0.7,])==0){
            print(kable(results[results[,3]>=0.6,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[results[,3]>=0.7,], format="markdown",row.names=FALSE))
          }
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank == 'pc'){
          if (u_lab=='c'){
            flag <- 0
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else {
            results <- results[order(abs(results[,4]),decreasing=TRUE),]
            #Only some results are shown on the screen
            if (nrow(results[abs(results[,4])>=0.6,])==0){
              print(kable(results[abs(results[,4])>=0.5,], format="markdown",row.names=FALSE))
            } else {
              print(kable(results[abs(results[,4])>=0.6,], format="markdown",row.names=FALSE))
            }
            #Creation of the excel file with all the results
            filename <-'results.xlsx'
            if (file.exists(filename)){
              a_xlsx <- file.remove(filename)
              write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
            } else {
              write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
            }
          }
        } else if (u_rank =='sc'){
          if (u_lab =='c'){
            flag <- 0
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else {
            results <- results[order(abs(results[,5]),decreasing=TRUE),]
            #Only some results are shown on the screen
            if (nrow(results[abs(results[,5])>=0.6,])==0){
              print(kable(results[abs(results[,5])>=0.5,], format="markdown",row.names=FALSE))
            } else {
              print(kable(results[abs(results[,5])>=0.6,], format="markdown",row.names=FALSE))
            }
            #Creation of the excel file with all the results
            filename <-'results.xlsx'
            if (file.exists(filename)){
              a_xlsx <- file.remove(filename)
              write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
            } else {
              write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
            }
          }
        } else {
          if (u_lab=='c'){
            flag <- 0 
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else {
            flag <- 0 
            cat("Please introduce either 'p', 'auc', 'aupr', 'pc' or 'sc'\n")
          }
        }
      } else {
        cat("\nWould you like to rank the PCA results by Pearson or Spearman correlation? [pc/sc]:\n\n")
        u_rank <- readline(prompt="-> ")
        flag <- 1
        if (u_rank =='pc'){
          results <- results[order(abs(results[,1]),decreasing=TRUE),]
          #Only some results are shown on the screen
          if (nrow(results[abs(results[,1])>=0.6,])==0){
            print(kable(results[abs(results[,1])>=0.5,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[abs(results[,1])>=0.6,], format="markdown",row.names=FALSE))
          }
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank=='sc'){
          results <- results[order(abs(results[,2]),decreasing=TRUE),]
          #Only some results are shown on the screen
          if (nrow(results[abs(results[,2])>=0.6,])==0){
            print(kable(results[abs(results[,2])>=0.5,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[abs(results[,2])>=0.6,], format="markdown",row.names=FALSE))
          }
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else {
          flag <- 0
          cat("Please introduce either 'pc' or 'sc'\n'")
        }
      }
    }
  } else {
    if ((u_lab=='c')|(u_lab=='d')){
      pval1 <- aperm(mw_ncPCA, c(1,3,2))
      pval1 <- matrix(pval1,12*dim(mw_ncPCA)[3],choose(numbLabels,2))
      storage.mode(pval1)<- "numeric"
      pval2 <- aperm(mw_cPCA,c(1,3,2))
      pval2 <- matrix(pval2,12*dim(mw_cPCA)[3],choose(numbLabels,2))
      storage.mode(pval2)<- "numeric"
      pvals <- rbind(pval1,pval2)
      avg <- apply(pvals, 1, mean)
      
      AUC1 <- aperm(AUC_nc,c(1,3,2))
      AUC1 <- matrix(AUC1,12*dim(mw_ncPCA)[3],choose(numbLabels,2))
      storage.mode(AUC1)<- "numeric"
      AUC2 <- aperm(AUC_c,c(1,3,2))
      AUC2 <- matrix(AUC2,12*dim(mw_cPCA)[3],choose(numbLabels,2))
      storage.mode(AUC2)<- "numeric"
      AUCs <- rbind(AUC1,AUC2)
      avg_AUC <- apply(AUCs, 1, mean)
      
      AUPR1 <- aperm(AUPR_nc,c(1,3,2))
      AUPR1 <- matrix(AUPR1,12*dim(mw_ncPCA)[3],choose(numbLabels,2))
      storage.mode(AUPR1)<- "numeric"
      AUPR2 <- aperm(AUPR_c,c(1,3,2))
      AUPR2 <- matrix(AUPR2,12*dim(mw_cPCA)[3],choose(numbLabels,2))
      storage.mode(AUPR2)<- "numeric"
      AUPRs <- rbind(AUPR1, AUPR2)
      avg_AUPR <- apply(AUPRs, 1, mean)
      
      o <- 1
      p <- 2
      
      group_head <-  character()
      for (i in 1:choose(numbLabels,2)){
        group_head[i] <- paste(nameLabels[o]," vs ",nameLabels[p], sep="")
        
        p <- p + 1
        if (p > numbLabels){
          o <- o + 1
          p <- o + 1
        }
      }
      
      if (u_lab =='c'){
        header_xls <- c('Avg P-val',group_head,'Avg AUC',group_head,'Avg AUPR',group_head,'Norm','Centering','Dimension','explained Variance')
        results_xls<- data.frame(avg,pvals,avg_AUC,AUCs,avg_AUPR, AUPRs, norms,centr,dim,variance, row.names = NULL)
        colnames(results_xls) <- header_xls
      }
      
      if (u_lab=='d'){
        header_xls <- c('Avg P-val',group_head,'Avg AUC',group_head,'Avg AUPR',group_head,'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance')
        pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
        spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
        results_xls<- data.frame(avg,pvals,avg_AUC,AUCs,avg_AUPR, AUPRs,pears_corr, spear_corr,norms,centr,dim,variance, row.names = NULL)
        colnames(results_xls) <- header_xls
      } 
      
    } else {
      pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
      spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
      header_xls <- c('pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance')
      results_xls <- data.frame(pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results_xls) <- header_xls
    }
    
    if (u_lab =='c'){
      header <- c('Avg P-val','Avg AUC','Avg AUPR','Norm','Centering','Dimension','expl Var')
      results <- data.frame(avg,avg_AUC,avg_AUPR,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
    } else if (u_lab =='d'){
      header <- c('Avg Pval','Avg AUC','Avg AUPR','pears','spear','Norm','Centering','Dim','expl Var')
      results <- data.frame(avg,avg_AUC,avg_AUPR,pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
    } else {
      header <- c('pears','spear','Norm','Centering','Dimension','expl Var')
      results <- data.frame(pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
    }
    
    flag <- 0
    while (flag == 0){
      if ((u_lab =='c')|(u_lab =='d')){
        if (u_lab =='c'){ 
          cat("\nWould you like to rank the PCA results by P-value, AUC or AUPR [p/auc/aupr]:\n\n")
          u_rank <- readline(prompt="-> ")
        } else {
          cat("\nWould you like to rank the PCA results by P-value, AUC, AUPR, Pearson correlation or Spearman correlation [p/auc/aupr/pc/sc]:\n\n")
          u_rank <- readline(prompt="-> ")
        }
        
        flag <- 1
        
        if (u_rank =='p'){
          
          idx <- order(results[,1],decreasing=FALSE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,1]<0.05,])==0){
            print(kable(results[results[,1]<0.25,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[results[,1]<0.05,], format="markdown",row.names=FALSE))
          }
          
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
          
        } else if (u_rank == 'auc'){
          
          idx <- order(results[,2],decreasing=TRUE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,2]>=0.7,])==0){
            print(kable(results[results[,2]>=0.6,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[results[,2]>=0.7,], format="markdown",row.names=FALSE))
          }
          
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank == 'aupr'){
          
          idx <- order(results[,3],decreasing=TRUE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,3]>=0.7,])==0){
            print(kable(results[results[,3]>=0.6,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[results[,3]>=0.7,], format="markdown",row.names=FALSE))
          }
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank =='pc'){
          
          if (u_lab =='c'){
            flag <- 0
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else {
            
            idx <- order(abs(results[,4]),decreasing=TRUE)
            results <- results[idx,]
            results_xls <- results_xls[idx,]
            #Only some results are shown on the screen
            if (nrow(results[abs(results[,4])>=0.6,])==0){
              print(kable(results[abs(results[,4])>=0.5,], format="markdown",row.names=FALSE))
            } else {
              print(kable(results[abs(results[,4])>=0.6,], format="markdown",row.names=FALSE))
            }
            #Creation of the excel file with all the results
            filename <-'results.xlsx'
            if (file.exists(filename)){
              a_xlsx <- file.remove(filename)
              write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
            } else {
              write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
            }
          }
        } else if (u_rank =='sc'){
          if (u_lab =='c'){
            flag <- 0
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else {
            
            idx <- order(abs(results[,5]),decreasing=TRUE)
            results <- results[idx,]
            results_xls <- results_xls[idx,]
            #Only some results are shown on the screen
            if (nrow(results[abs(results[,5])>=0.6,])==0){
              print(kable(results[abs(results[,5])>=0.5,], format="markdown",row.names=FALSE))
            } else {
              print(kable(results[abs(results[,5])>=0.6,], format="markdown",row.names=FALSE))
            }
            #Creation of the excel file with all the results
            filename <-'results.xlsx'
            if (file.exists(filename)){
              a_xlsx <- file.remove(filename)
              write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
            } else {
              write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
            }
          }
        } else {
          if (u_lab=='c'){
            flag <- 0
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else{
            flag <- 0
            cat("Please introduce either 'p', 'auc', 'aupr', 'pc' or 'sc'\n")
          }
        }
      } else {
        cat("\nWould you like to rank the PCA results by Pearson or Spearman correlation? [pc/sc]:\n\n'")
        u_rank <- readline(prompt="-> ")
        
        flag <- 1
        
        if (u_rank == 'pc'){
          idx <- order(abs(results[,1]),decreasing=TRUE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[abs(results[,1])>=0.6,])==0){
            print(kable(results[abs(results[,1])>=0.5,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[abs(results[,1])>=0.6,], format="markdown",row.names=FALSE))
          }
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank == 'sc'){
          idx <- order(abs(results[,2]),decreasing=TRUE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[abs(results[,2])>=0.6,])==0){
            print(kable(results[abs(results[,2])>=0.5,], format="markdown",row.names=FALSE))
          } else {
            print(kable(results[abs(results[,2])>=0.6,], format="markdown",row.names=FALSE))
          }
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else {
          flag <- 0
          cat("Please introduce either 'pc' or 'sc'\n'")
        }
      }
    }
  }
  
  
  # User interaction --------------------------------------------------------
  flag <- 0
  while (flag == 0){
    cat("\nSelect the normalization:\n\n")
    u_norm_n <- readline(prompt="-> ")
    flag <- 1
    normal <- c('DSOR', 'DSOC', 'LOG', 'ZSCORE', 'QUANTILE T', 'QUANTILE', 'ZSCORE T', 'PLUS(ABS(MIN))', 'PARETO SCALING', 'SQRT', 'MANORM',  '-')
    if (sum(u_norm_n==normal)==1){
      u_norm <- switch(u_norm_n, 'DSOR'= 1, 'DSOC'= 2, 'LOG'= 3, 'ZSCORE'= 4, 'QUANTILE T' = 5, 'QUANTILE' = 6,
                       'ZSCORE T'= 7, 'PLUS(ABS(MIN))' = 8, 'PARETO SCALING'= 9, 'SQRT' = 10, 'MANORM' = 11, '-' = 12) 
    } else {
      flag <- 0
      cat("Please introduce the exact name of the normalization.\n")
    }
  }
  
  
  flag <- 0
  while (flag == 0){
    cat('\nCentering version? [y/n]:\n\n')
    u_cent <- readline(prompt="-> ")
    if ((u_cent=='n')|(u_cent=='y')){
      flag <- 1
    } else{
      cat("Please introduce just 'y' or 'n'\n")
    }
  }
  
  
  flag <- 0 
  while (flag == 0){
    cat('\nSelect the dimension for generating the PC-corr network:\n\n')
    u_dim <- readline(prompt="-> ")
    storage.mode(u_dim) <- "numeric"
    if ((u_dim > 0) & (u_dim <= length(labels))){
      flag <- 1
    } else{
      cat("Please introduce an existing dimension.\n")
    }
  }
  
  
  flag <- 0 
  while (flag == 0){
    cat('\nSelect a cut-off or a set of cut-offs for generating the PC-corr network [number between 0 and 1]:\n\n')
    cat("Example: c(0.6, 0.65, 0.7)\n\n")
    cutoff <- readline(prompt="-> ")
    cutoff <- eval(parse(text=cutoff))
    if ((max(cutoff) >= 0) & (max(cutoff) <= 1)){
      flag <- 1
    } else{
      cat("Please introduce a correct cut-off or a set of cut-offs (in [0,1]).\n")
    }
  }
  
  if (u_cent=='y'){
    PCA <- cPCA
    if (u_lab == 'c'){
      if (u_rank == 'p'){
        mw_PCA <- mw_cPCA
        mw_PCA1 <- AUC_c
        mw_PCA2 <- AUPR_c
        evaluat <- 'P-value'
      } else if (u_rank == 'auc'){
        mw_PCA <- AUC_c
        mw_PCA1 <- mw_cPCA
        mw_PCA2 <- AUPR_c
        evaluat <- 'AUC'
      } else if (u_rank == 'aupr'){
        mw_PCA <- AUPR_c
        mw_PCA1 <- mw_cPCA
        mw_PCA2 <- AUC_c
        evaluat <- 'AUPR'
      }
    } else if (u_lab == 'd'){
      if (u_rank == 'p'){
        mw_PCA <- mw_cPCA
        mw_PCA1 <- AUC_c
        mw_PCA2 <- AUPR_c
        mw_PCA3 <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        mw_PCA4 <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        evaluat <- 'P-value'
      } else if (u_rank == 'auc'){
        mw_PCA <- AUC_c
        mw_PCA1 <- mw_cPCA
        mw_PCA2 <- AUPR_c
        mw_PCA3 <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        mw_PCA4 <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        evaluat <- 'AUC'
      } else if (u_rank == 'aupr'){
        mw_PCA <- AUPR_c
        mw_PCA1 <- mw_cPCA
        mw_PCA2 <- AUC_c
        mw_PCA3 <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        mw_PCA4 <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        evaluat <- 'AUPR'
      } else if (u_rank == 'pc'){
        mw_PCA <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        mw_PCA1 <- mw_cPCA
        mw_PCA2 <- AUC_c
        mw_PCA3 <- AUPR_c
        mw_PCA4 <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        evaluat <- 'Pearson correlation'
      } else if (u_rank == 'sc'){
        mw_PCA <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        mw_PCA1 <- mw_cPCA
        mw_PCA2 <- AUC_c
        mw_PCA3 <- AUPR_c
        mw_PCA4 <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        evaluat <- 'Spearman correlation'
      }
    } else {
      if (u_rank == 'pc'){
        mw_PCA <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        mw_PCA1 <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        evaluat <- 'Pearson correlation'
      } else {
        mw_PCA <- aperm(rank_spear_corr_cPCA, c(1,3,2))
        mw_PCA1 <- aperm(rank_pears_corr_cPCA, c(1,3,2))
        evaluat <- 'Spearman correlation'
      }
    }
    
    pc <- pc_c
    ttl <- 'centred PCA'
    explained <- explained_c
  } else if (u_cent == 'n'){
    PCA <- ncPCA
    if (u_lab == 'c'){
      if (u_rank == 'p'){
        mw_PCA <- mw_ncPCA
        mw_PCA1 <- AUC_nc
        mw_PCA2 <- AUPR_nc
        evaluat <- 'P-value'
      } else if (u_rank == 'auc'){
        mw_PCA <- AUC_nc
        mw_PCA1 <- mw_ncPCA
        mw_PCA2 <- AUPR_nc
        evaluat <- 'AUC'
      } else if (u_rank == 'aupr'){
        mw_PCA <- AUPR_nc
        mw_PCA1 <- mw_ncPCA
        mw_PCA2 <- AUC_nc
        evaluat <- 'AUPR'
      }
    } else if (u_lab == 'd'){
      if (u_rank == 'p'){
        mw_PCA <- mw_ncPCA
        mw_PCA1 <- AUC_nc
        mw_PCA2 <- AUPR_nc
        mw_PCA3 <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        mw_PCA4 <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        evaluat <- 'P-value'
      } else if (u_rank == 'auc'){
        mw_PCA <- AUC_nc
        mw_PCA1 <- mw_ncPCA
        mw_PCA2 <- AUPR_nc
        mw_PCA3 <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        mw_PCA4 <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        evaluat <- 'AUC'
      } else if (u_rank == 'aupr'){
        mw_PCA <- AUPR_nc
        mw_PCA1 <- mw_ncPCA
        mw_PCA2 <- AUC_nc
        mw_PCA3 <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        mw_PCA4 <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        evaluat <- 'AUPR'
      } else if (u_rank == 'pc'){
        mw_PCA <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        mw_PCA1 <- mw_ncPCA
        mw_PCA2 <- AUC_nc
        mw_PCA3 <- AUPR_nc
        mw_PCA4 <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        evaluat <- 'Pearson correlation'
      } else if (u_rank == 'sc'){
        mw_PCA <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        mw_PCA1 <- mw_ncPCA
        mw_PCA2 <- AUC_nc
        mw_PCA3 <- AUPR_nc
        mw_PCA4 <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        evaluat <- 'Spearman correlation'
      }
    } else {
      if (u_rank == 'pc'){
        mw_PCA <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        mw_PCA1 <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        evaluat <- 'Pearson correlation'
      } else {
        mw_PCA <- aperm(rank_spear_corr_ncPCA, c(1,3,2))
        mw_PCA1 <- aperm(rank_pears_corr_ncPCA, c(1,3,2))
        evaluat <- 'Spearman correlation'
      }
    }
    
    pc <- pc_nc
    ttl <- 'non-centred PCA'
    explained <- explained_nc
  }
  
  ind2 <- c()
  if ((u_lab == 'c')|(u_lab == 'd')){
    if (u_rank == 'p') {
      
      #P-value
      if (numbLabels == 2){
        val <- min(mw_PCA[u_norm,1,])
        ind1 <- which.min(mw_PCA[u_norm,1,])
        if (ind1 == u_dim){
          if (ind1 == length(labels)){
            val <- min(mw_PCA[u_norm,1,1:(ind1-1)])
            ind2 <- which.min(mw_PCA[u_norm,1,1:(ind1-1)])
          } else if (ind1 ==1){
            val <- min(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])
            ind2 <- which.min(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])
          } else {
            val <- min(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])
            ind2 <- which.min(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])
          }
        }
      } else {
        val <- min(colMeans(apply(mw_PCA[u_norm, ,], c(1,2), as.numeric)))
        ind1 <- which.min(colMeans(apply(mw_PCA[u_norm, ,], c(1,2), as.numeric)))
        if (ind1 == u_dim){
          if (ind1 == length(labels)){
            val <- min(colMeans(apply(mw_PCA[u_norm, ,1:(ind1-1)], c(1,2), as.numeric)))
            ind2 <- which.min(colMeans(apply(mw_PCA[u_norm, ,1:(ind1-1)], c(1,2), as.numeric)))
          } else if (ind1 == 1) {
            val <- min(colMeans(apply(mw_PCA[u_norm, ,(ind1+1):dim(mw_PCA)[3]], c(1,2), as.numeric)))
            ind2 <- which.min(colMeans(apply(mw_PCA[u_norm, ,(ind1+1):dim(mw_PCA)[3]], c(1,2), as.numeric)))
          } else {
            val <- min(colMeans(apply(mw_PCA[u_norm, ,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])], c(1,2), as.numeric)))
            ind2 <- which.min(colMeans(apply(mw_PCA[u_norm, ,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])], c(1,2), as.numeric)))
          }
        }
      }
    } else if ((u_rank == 'pc')|(u_rank == 'sc')){
      
      #Pearson correlation and Spearman correlation
      val <- max(abs(unlist(mw_PCA[u_norm,1,])))
      ind1 <- which.max(abs(unlist(mw_PCA[u_norm,1,])))
      if (ind1 == u_dim){
        if (ind1 == length(labels)){
          val <- max(abs(unlist(mw_PCA[u_norm,1,1:(ind1-1)])))
          ind2 <- which.max(abs(unlist(mw_PCA[u_norm,1,1:(ind1-1)])))
        } else if (ind1 == 1){
          val <- max(abs(unlist(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])))
          ind2 <- which.max(abs(unlist(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])))
        } else {
          val <- max(abs(unlist(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])))
          ind2 <- which.max(abs(unlist(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])))
        }
      }
    } else if ((u_rank == 'auc')|(u_rank == 'aupr')) {
      
      #AUC and AUPR
      if (numbLabels == 2){
        val <- max(mw_PCA[u_norm,1,])
        ind1 <- which.max(mw_PCA[u_norm,1,])
        if (ind1 == u_dim){
          if (ind1 ==length(labels)){
            val <- max(mw_PCA[u_norm,1,1:(ind1-1)])
            ind2 <- which.max(mw_PCA[u_norm,1,1:(ind1-1)])
          } else if (ind1 == 1){
            val <- max(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])
            ind2 <- which.max(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])
          } else {
            val <- max(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])
            ind2 <- which.max(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])
          }
        }
      } else {
        val <- max(colMeans(apply(mw_PCA[u_norm, ,], c(1,2), as.numeric)))
        ind1 <- which.max(colMeans(apply(mw_PCA[u_norm, ,], c(1,2), as.numeric)))
        if (ind1 == u_dim){
          if (ind1 == length(labels)){
            val <- max(colMeans(apply(mw_PCA[u_norm, ,1:(ind1-1)], c(1,2), as.numeric)))
            ind2 <- which.max(colMeans(apply(mw_PCA[u_norm, ,1:(ind1-1)], c(1,2), as.numeric)))
          } else if (ind1 == 1){
            val <- max(colMeans(apply(mw_PCA[u_norm, ,(ind1+1):dim(mw_PCA)[3]], c(1,2), as.numeric)))
            ind2 <- which.max(colMeans(apply(mw_PCA[u_norm, ,(ind1+1):dim(mw_PCA)[3]], c(1,2), as.numeric)))
          } else { 
            val <- max(colMeans(apply(mw_PCA[u_norm, ,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])], c(1,2), as.numeric)))
            ind2 <- which.max(colMeans(apply(mw_PCA[u_norm, ,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])], c(1,2), as.numeric)))
          }
        }
      }
    }
  } else {
    val <- max(abs(mw_PCA[u_norm,1,]))
    ind1 <- which.max(abs(mw_PCA[u_norm,1,]))
    if (ind1 == u_dim){
      if (ind1 == length(labels)){
        val <- max(abs(unlist(mw_PCA[u_norm,1,1:(ind1-1)])))
        ind2 <- which.max(abs(unlist(mw_PCA[u_norm,1,1:(ind1-1)])))
      } else if (ind1 == 1) {
        val <- max(abs(unlist(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])))
        ind2 <- which.max(abs(unlist(mw_PCA[u_norm,1,(ind1+1):dim(mw_PCA)[3]])))
      } else {
        val <- max(abs(unlist(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])))
        ind2 <- which.max(abs(unlist(mw_PCA[u_norm,1,c(1:(ind1-1),(ind1+1):dim(mw_PCA)[3])])))
      }
    }
  }
  
  if (length(ind2)==0) {
    ind <- ind1
  } else {
    if (ind2 >= ind1){
      ind <- ind2 + 1
    } else {
      ind <- ind2 
    }
  }
  
  
  
  # Scatter plot of the desired PCA -----------------------------------------
  #Dimension 1: Best discriminating dimension according to the evaluator
  #Dimension 2: Chosen dimension
  
  if (numbLabels == 2){
    if (u_rank == 'p'){
      if (val <= mw_PCA[u_norm,1,u_dim]){
        dim1 <- ind
        dim2 <- u_dim
      } else {
        dim1 <- u_dim
        dim2 <- ind
      }
    } else {
      if (abs(val) >= abs(mw_PCA[u_norm,1,u_dim])){
        dim1 <- ind
        dim2 <- u_dim
      } else {
        dim1 <- u_dim
        dim2 <- ind
      }
    }
  } else {
    if (u_rank == 'p'){
      if (mean(unlist(mw_PCA[u_norm,,ind])) <= mean(unlist(mw_PCA[u_norm,,u_dim]))){
        dim1 <- ind
        dim2 <- u_dim
      } else {
        dim1 <- u_dim
        dim2 <- ind
      }
    } else if ((u_rank == 'auc')|(u_rank =='aupr')){
      if (mean(unlist(mw_PCA[u_norm,,ind])) >= mean(unlist(mw_PCA[u_norm,,u_dim]))){
        dim1 <- ind
        dim2 <- u_dim
      } else {
        dim1 <- u_dim
        dim2 <- ind
      }
    } else {
      if (abs(val) >= abs(mw_PCA[u_norm,1,u_dim])){
        dim1 <- ind
        dim2 <- u_dim
      } else {
        dim1 <- u_dim
        dim2 <- ind
      }
    }
  }
  
  # Colours of the groups ---------------------------------------------------
  #Group on the left of the most discriminative dimension in the PCA plot (median value) : black
  #Group on the right of the most discriminative dimension in the PCA plot (median value) : red
  #If there are more than two groups, the colours of the other groups are random
  
  a <- vector('numeric')
  
  for (i in 1:numbLabels){
    a[i] <- median(PCA[[u_norm]][labels %in% nameLabels[i],dim1])
  }
  I <- order(a)
  col <- c()
  col[I[1]]<- 'black'
  col[I[numbLabels]] <- 'red'
  if (numbLabels > 2){
    for (i in 2:(numbLabels-1)){
      runif_numb <- runif(3,0.05,0.95)
      col[I[i]] <- rgb(runif_numb[1],runif_numb[2],runif_numb[3])
    }
  }
  
  # Density plots
  # x axis: density plot
  
  grp <- matrix(list(),numbLabels,1)
  f <- matrix(list(),numbLabels,1)
  xi <- matrix(list(),numbLabels,1)
  for (i in 1:numbLabels){
    grp[[i]] <- PCA[[u_norm]][labels %in% nameLabels[i],dim1]
    denst <- density(grp[[i]],kernel='gaussian',bw='nrd0', n=100)
    xi[[i]] <- denst$x
    f[[i]] <- denst$y
  }
  
  gx_Bottom <- ggplot()
  for (i in 1:numbLabels){
    dat_dens <- data.frame(cbind(xi[[i]],f[[i]]))
    colnames(dat_dens) <- c("xi","f")
    gx_Bottom <- gx_Bottom + geom_line(data=dat_dens,aes(x = xi, y = f), colour = col[i]) + theme(legend.position="none") 
  }
  gx_Bottom <- gx_Bottom + xlim(min(unlist(xi)), max(unlist(xi))) +xlab(" ") + theme(axis.title.y=element_blank()) + ylab(" ")
  
  # y axis: density plot
  grp1 <- matrix(list(),numbLabels,1)
  f1 <- matrix(list(),numbLabels,1)
  xi1 <- matrix(list(),numbLabels,1)
  for (i in 1:numbLabels){
    grp1[[i]] <- PCA[[u_norm]][labels %in% nameLabels[i],dim2]
    denst1 <- density(grp1[[i]],kernel='gaussian',bw='nrd0', n=100)
    xi1[[i]] <- denst1$x
    f1[[i]] <- denst1$y
  }
  
  gx_Left <- ggplot()
  for (i in 1:numbLabels){
    dat_dens1 <- data.frame(cbind(xi1[[i]],f1[[i]]))
    colnames(dat_dens1) <- c("xi1","f1")
    gx_Left <- gx_Left + geom_line(data=dat_dens1,aes(x = xi1, y = f1), colour = col[i]) + coord_flip() + theme(legend.position="none") 
  }
  gx_Left <- gx_Left + xlim(min(unlist(xi1)), max(unlist(xi1))) + ylab(" ") + xlab(" ") +labs(title="")
  
  # Scatterplot
  data_plt_all <- data.frame()
  col_plot <- c()
  for (i in 1:numbLabels){ 
    data_plt<- data.frame(PCA[[u_norm]][labels %in% nameLabels[i],c(dim1,dim2)],nameLabels[i],sample_names[labels %in% nameLabels[i]])
    colnames(data_plt) <- c("dim1","dim2","gr","label")
    data_plt_all <- rbind(data_plt_all,data_plt)
    col_plot <- c(col_plot,col[i])
  }
  
  if (dis == 'yes'){
    gx_Main <- ggplot(data_plt_all, aes(dim1,dim2)) + geom_point(aes(fill = gr),color="black",pch=21, size=5)+ geom_text_repel(aes(data_plt_all$dim1,data_plt_all$dim2, label = data_plt_all$label))+ scale_fill_manual(values = col_plot) + theme(legend.title=element_blank())+
      theme(legend.position=c(1,1),legend.justification=c(0.95,0.95),legend.background = element_rect(colour = 'black', fill = 'gray90', size = 0.5, linetype='solid'))
  } else if (dis == 'no'){
    gx_Main <- ggplot(data_plt_all, aes(dim1,dim2)) + geom_point(aes(fill = gr),color="black",pch=21, size=5)+ scale_fill_manual(values = col_plot) + theme(legend.title=element_blank())+
      theme(legend.position=c(1,1),legend.justification=c(0.95,0.95),legend.background = element_rect(colour = 'black', fill = 'gray90', size = 0.5, linetype='solid'))
    
  }
  
  # Labels in the axes ------------------------------------------------------
  # 2 Groups of labels ------------------------------------------------------
  if (numbLabels == 2){
    if (u_lab == 'c'){
      
      #class labels
      if (u_rank == 'p'){
        #p-value
        if (mw_PCA[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'auc'){
        #AUC
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'aupr'){
        #AUPR
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
        }
      }
    } else if (u_lab == 'd'){
      
      #discrete labels
      if (u_rank == 'p'){
        #p-value
        if (mw_PCA[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'auc'){
        #AUC
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'aupr'){
        #AUPR
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank =='pc'){
        #Pearson correlation
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'sc'){
        #Spearman correlation
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
        }
      }
    } else {
      if (u_rank == 'pc'){
        str_xlab <- c("PC",as.character(dim1)," ( pears = ", as.character(round(mw_PCA[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA1[u_norm,1,dim1],3))," )")
        str_ylab <- c("PC",as.character(dim2)," ( pears = ", as.character(round(mw_PCA[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA1[u_norm,1,dim2],3))," )")
      } else if (u_rank == 'sc'){
        str_xlab <- c("PC",as.character(dim1)," ( pears = ", as.character(round(mw_PCA1[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        str_ylab <- c("PC",as.character(dim2)," ( pears = ", as.character(round(mw_PCA1[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
      }
    }
  }
  
  
  
  # More than 2 Groups of labels ------------------------------------------------------
  if (numbLabels > 2){
    if (u_lab == 'c'){
      
      #class labels
      if (u_rank == 'p'){
        #p-value
        if (mean(mw_PCA[u_norm, ,dim1])>= 0.01){
          pdis1 <- as.character(round(mean(mw_PCA[u_norm, ,dim1]),3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mean(mw_PCA1[u_norm, ,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm, ,dim1]),3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mean(mw_PCA1[u_norm, ,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm, ,dim1]),3))," )")
        }
        if (mean(mw_PCA[u_norm, ,dim2])>= 0.01){
          pdis2 <- as.character(round(mean(mw_PCA[u_norm, ,dim2]),3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mean(mw_PCA1[u_norm, ,dim2]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm, ,dim2]),3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mean(mw_PCA1[u_norm, ,dim2]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm, ,dim2]),3))," )")
        }
      } else if (u_rank == 'auc'){
        #AUC
        if (mean(mw_PCA1[u_norm,,dim1])>= 0.01){
          pdis1 <- as.character(round(mean(mw_PCA1[u_norm,,dim1]),3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mean(mw_PCA[u_norm,,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm,,dim1]),3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mean(mw_PCA[u_norm,,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm,,dim1]),3))," )")
        }
        if (mean(mw_PCA1[u_norm,,dim2])>= 0.01){
          pdis2 <- as.character(round(mean(mw_PCA1[u_norm,,dim2]),3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mean(mw_PCA[u_norm,,dim2]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm,,dim2]),3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mean(mw_PCA[u_norm,,dim2]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm,,dim2]),3))," )")
        }
      } else if (u_rank == 'aupr'){
        #AUPR
        if (mean(mw_PCA1[u_norm,,dim1])>= 0.01){
          pdis1 <- as.character(round(mean(mw_PCA1[u_norm,,dim1]),3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mean(mw_PCA2[u_norm,,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA[u_norm,,dim1]),3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mean(mw_PCA2[u_norm,,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA[u_norm,,dim1]),3))," )")
        }
        if (mean(mw_PCA1[u_norm,,dim2])>= 0.01){
          pdis2 <- as.character(round(mean(mw_PCA1[u_norm,,dim2]),3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mean(mw_PCA2[u_norm,,dim2]),3)),", AUPR = ",as.character(round(mean(mw_PCA[u_norm,,dim2]),3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mean(mw_PCA2[u_norm,,dim2]),3)),", AUPR = ",as.character(round(mean(mw_PCA[u_norm,,dim2]),3))," )")
        }
      }
    } else if (u_lab == 'd'){
      
      #discrete labels
      if (u_rank == 'p'){
        #p-value
        if (mean(mw_PCA[u_norm,,dim1])>= 0.01){
          pdis1 <- as.character(round(mean(mw_PCA[u_norm,,dim1]),3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mean(mw_PCA1[u_norm,,dim1]),3)),", AUPR = ",as.character(round(mean(mw_PCA2[u_norm,,dim1]),3)),", pears = ", as.character(round(mean(mw_PCA3[u_norm,,dim1]),3)),", spear = ", as.character(round(mean(mw_PCA4[u_norm,,dim1]),3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA1[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'auc'){
        #AUC
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'aupr'){
        #AUPR
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA3[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank =='pc'){
        #Pearson correlation
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA4[u_norm,1,dim2],3))," )")
        }
      } else if (u_rank == 'sc'){
        #Spearman correlation
        if (mw_PCA1[u_norm,1,dim1]>= 0.01){
          pdis1 <- as.character(round(mw_PCA1[u_norm,1,dim1],3))
          str_xlab <- c("PC",as.character(dim1)," ( p = ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        } else {
          pdis1 <- "p < 0.01"
          str_xlab <- c("PC",as.character(dim1)," ( ",pdis1,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim1],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim1],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        }
        if (mw_PCA1[u_norm,1,dim2]>= 0.01){
          pdis2 <- as.character(round(mw_PCA1[u_norm,1,dim2],3))
          str_ylab <- c("PC",as.character(dim2)," ( p = ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
        } else {
          pdis2 <- "p < 0.01"
          str_ylab <- c("PC",as.character(dim2)," ( ",pdis2,", AUC = ",as.character(round(mw_PCA2[u_norm,1,dim2],3)),", AUPR = ",as.character(round(mw_PCA3[u_norm,1,dim2],3)),", pears = ", as.character(round(mw_PCA4[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
        }
      }
    } else {
      if (u_rank == 'pc'){
        str_xlab <- c("PC",as.character(dim1)," ( pears = ", as.character(round(mw_PCA[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA1[u_norm,1,dim1],3))," )")
        str_ylab <- c("PC",as.character(dim2)," ( pears = ", as.character(round(mw_PCA[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA1[u_norm,1,dim2],3))," )")
      } else if (u_rank == 'sc'){
        str_xlab <- c("PC",as.character(dim1)," ( pears = ", as.character(round(mw_PCA1[u_norm,1,dim1],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim1],3))," )")
        str_ylab <- c("PC",as.character(dim2)," ( pears = ", as.character(round(mw_PCA1[u_norm,1,dim2],3)),", spear = ", as.character(round(mw_PCA[u_norm,1,dim2],3))," )")
      }
    }
  }
  
  
  str_title <- c(ttl, " for norm: ",u_norm_n)
  gx_Main <- gx_Main + xlab(paste(str_xlab,collapse="")) + ylab(paste(str_ylab,collapse="")) + labs(title=paste(str_title, collapse = "")) +theme(plot.title = element_text(hjust = 0.5))+
    xlim(min(unlist(xi)), max(unlist(xi))) + ylim(min(unlist(xi1)), max(unlist(xi1))) 
  empty <- ggplot()+geom_point(aes(1,1), colour="white") +
    theme(                              
      plot.background = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )
  
  dev.new()
  #grid.arrange(gx_Left, gx_Main, empty, gx_Bottom, ncol=2, nrow=2, widths=c(2, 4), heights=c(4, 2))
  grid.arrange( gx_Left, gx_Main, empty, gx_Bottom,  layout_matrix= matrix(c(1,1,1,3,2,2,2,4,2,2,2,4,2,2,2,4),nrow=4))

  # bar plots ----------------------------------------------------------------
  if (numbLabels == 2){
    if ((u_rank !='pc') && (u_rank!='sc')){
      bar_data <- data.frame(dimens=factor(1:length(labels)),evaluator=mw_PCA[u_norm,1,])
      g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator)) + geom_bar(stat='identity',fill='steelblue')
    } else {
      a1 <- which(mw_PCA[u_norm,1,] >= 0)
      a2 <- which(mw_PCA[u_norm,1,] < 0)
      lev <- c()
      lev[a1] <- "pos"
      lev[a2] <- "neg"
      bar_data <- data.frame(dimens=factor(1:length(labels)),evaluator=abs(mw_PCA[u_norm,1,]),lev)
      g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=lev)) + geom_bar(stat='identity') 
      if (u_rank == 'pc'){
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Pear. corr < 0","Pear. corr \u2265 0"),values=c("#808080", "#000000"))
      } else {
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Spear. corr < 0","Spear. corr \u2265 0"),values=c("#808080", "#000000"))
      }
    }
  } else {
    if ((u_lab == 'c')||(u_lab == 'd')){
      if ((u_rank !='pc') && (u_rank!='sc')){
        bar_data <- data.frame(dimens=factor(1:length(labels)),evaluator=colMeans(mw_PCA[u_norm,,]))
        g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator)) + geom_bar(stat='identity',fill='steelblue')
      } else {
        a1 <- which(mw_PCA[u_norm,1,] >= 0)
        a2 <- which(mw_PCA[u_norm,1,] < 0)
        lev <- c()
        lev[a1] <- "pos"
        lev[a2] <- "neg"
        bar_data <- data.frame(dimens=factor(1:length(labels)),evaluator=abs(mw_PCA[u_norm,1,]),lev)
        g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=lev)) + geom_bar(stat='identity') 
        if (u_rank == 'pc'){
          g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Pear. corr < 0","Pear. corr \u2265 0"),values=c("#808080", "#000000"))
        } else {
          g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Spear. corr < 0","Spear. corr \u2265 0"),values=c("#808080", "#000000"))
        }
      }
    } else {
      a1 <- which(mw_PCA[u_norm,1,] >= 0)
      a2 <- which(mw_PCA[u_norm,1,] < 0)
      lev <- c()
      lev[a1] <- "pos"
      lev[a2] <- "neg"
      bar_data <- data.frame(dimens=factor(1:length(labels)),evaluator=abs(mw_PCA[u_norm,1,]),lev)
      g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=lev)) + geom_bar(stat='identity') 
      if (u_rank == 'pc'){
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Pear. corr < 0","Pear. corr \u2265 0"),values=c("#808080", "#000000"))
      } else {
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Spear. corr < 0","Spear. corr \u2265 0"),values=c("#808080", "#000000"))
      }
    }
  }
  eval_title <- c(evaluat,'s over principal components in ',ttl,' for norm: ',u_norm_n)
  g1 <- g1 + xlab("PC") + ylab(paste(str_ylab,collapse="")) + labs(title=paste(eval_title, collapse = "")) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, vjust= 0.5, hjust = 1),legend.position = c(0.92,0.9))
  if (( u_rank!='pc') && (u_rank!='sc')){
    g1 <- g1 + ylab(evaluat)
  } else {
    g1 <- g1 + ylab(paste(c('|',evaluat,'|'), collapse=""))
  }
  
  
  explain_var <- data.frame(dimens=factor(1:length(labels)),explain=explained[[u_norm]])
  g2 <- ggplot(data=explain_var, aes(x=dimens, y= explain))+ geom_bar(stat='identity', fill='steelblue') + labs(title= 'Explained variance for the respective principal components', x= 'PC', y =  'Explained Variance (%)') + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, vjust= 0.5, hjust = 1))
  
  dev.new()
  grid.arrange(g1,g2, ncol=1, nrow=2)   
  
  
  # PC-corr -----------------------------------------------------------------
  C_corr <- function(x, V, feat_names, cutoff) {
    #x: normalized data matrix (NxM), that was given as an input of PCA. 
    #The N rows contain the samples, while the M columns contain the features.
    #V: loadings of one PCA dimension (i.e. one column of the third output matrix of SVD function)
    #feat_names: the feature names  
    
    n<-nargs()
    if (n<3){
      stop("number of inputs must be at least three")
    }
    
    if(n<4){
      cutoff <- 0.7
    }
    
    V <- sign(V) * log10(1+abs(V)/mean(abs(V)))
    V <- sign(V) * ((abs(V)-min(abs(V)))/(max(abs(V))-min(abs(V))))
    
    # PC-corr on all the features
    
    x_temp <- x
    V_temp <- V
    
    
    index <- abs(V) > cutoff
    n <- length(V)-sum(index)

    
    
    if (n>0){
      ns <- toString(n) 
      if (n == 1){
        print(paste(ns,"feature was deleted because |V_feature|<cutoff",sep = " "));
      } else if (n > 1){
        print(paste(ns,"features were deleted because |V_feature|<cutoff",sep = " "));
      }
      x <- x[,index]
      V <- V[index]
      feat_names <- feat_names[index]
    }   
    

    
    # - Normal Pearson Corrleation calculation
    c <- cor(x)
    
    # - Apply the C-Corr formula
    pc_corr <- matrix(0,length(V),length(V))
    
    for (i in 1:(length(V)-1)) {
      for (j in (i+1):length(V)) {
        pc_corr[i,j] <- sign(c[i,j]) * min(abs(c(c[i,j],V[i],V[j])))
      }
    }
    
    if(max(abs(pc_corr[,])) <= cutoff) {
      c_temp <- cor(x_temp)
      pc_corr_temp <- matrix(0,length(V_temp),length(V_temp))
      for (i in 1:(length(V_temp)-1)){
        for (j in (i+1):length(V_temp)){
          pc_corr_temp[i,j] <- sign(c_temp[i,j]) * min(abs(c(c_temp[i,j],V_temp[i],V_temp[j])))
        }
      }
      warning('With this cut-off, there are no edges that have |PC_corr(i,j)| > cutoff')
      warning(paste('Try with another cut-off less than',toString(max(abs(pc_corr_temp[,]))), sep = " "))
    }
    
    idx <- (pc_corr < (-cutoff)) | (pc_corr > cutoff)
    pc_corr[idx==0] <- 0
    
    # Remove the singletons
    pc_corr <- pc_corr + t(pc_corr)
    ma <- (colSums(pc_corr[,]==0)) == dim(pc_corr)[1]
    pc_corr<- pc_corr[,!ma]
    pc_corr<- pc_corr[!ma,]
    sig_PC_corr <- pc_corr
    sig_PC_corr_Name <- feat_names[!ma]
    
    #Datamatrix with only the features in the network
    x1 <- x
    x1 <- x1[,!ma]
    
    # Create the edges and nodes variables
    # Edges
    node_i <- vector()
    node_j <- vector()
    pc_corr_ <- matrix()
    
    k <- 1;
    for (i in 1:(dim(sig_PC_corr)[1]-1)) {
      for (j in (i+1):dim(sig_PC_corr)[1]) {
        if(sig_PC_corr[i,j]!= 0) {
          node_i[k] <- sig_PC_corr_Name[i]
          node_j[k] <- sig_PC_corr_Name[j]
          pc_corr_[k] <- sig_PC_corr[i,j]
          
          k <- k+1
        }
      }
    }
    
    edges <- data.frame(node_i,node_j,pc_corr_)
    colnames(edges) <- c("Node i", "Node j", "PC-corr(i,j)")
    
    # Nodes
    node <- vector()
    loadings_V <- vector()
    
    Lcob <- is.element(feat_names,sig_PC_corr_Name)
    for (i in 1:length(sig_PC_corr_Name)) {
      node[i] <- sig_PC_corr_Name[i]
    }
    loadings_V <- V[Lcob]
    
    nodes <- data.frame(node,loadings_V)
    colnames(nodes) <- c("Node","Loading (V)")
    
    return(list(Edges=edges,Nodes=nodes,pc_corr=pc_corr, x1=x1))
    
  }
  
  
  match_V_samplCol <- function(col,x1,labels,nameLabels,Nodes){
    n <- which(col == 'red')
    m <- which(col == 'black')
    nodes_feat <- Nodes[,]
    NodeColor <- c()
    men1 <- c()
    men2 <- c()
    men3 <- c()
    
    for (i in 1:dim(x1)[2]){
      men1[i] <- nodes_feat[i,1]
      men2[i]<- mean(x1[labels %in% nameLabels[n],i]) #red group
      men3[i]<- mean(x1[labels %in% nameLabels[m],i]) #black group
    }
    
    men <- data.frame(men1,men2,men3)
    men_dif <- c()
    men_dif <- men[,2] - men[,3]
    m1 <- men_dif >= 0
    m2 <- men_dif <= 0
    
    b <- c()
    b <- nodes_feat[,2] <= 0 # 1 where it is a negative loading
    l1 <- sum(b==0)
    l2 <- sum(b==1)
    t1 <- sum(t(m1)*t(!b))/sum(b==0)
    t2 <- sum(t(m2)*t(b))/sum(b==1) 
    t3 <- sum(t(m1)*t(b))/sum(b==1)
    t4 <- sum(t(m2)*t(!b))/sum(b==0)
    
    if ((l1 == 0) && (t2 >= t3)){
      t <- t2
      NodeColor[b==1] <- 'Black'
      n1_f <- NaN
      n2_f <- 1-t
    } else if ((l1 == 0) && (t3 >= t2)){
      t <- t3
      NodeColor[b==1] <- 'Red'
      n1_f <- NaN
      n2_f <- 1-t
    } else if ((l2 == 0) && (t1 >= t4)){
      t <- t1
      NodeColor[b==0] <- 'Red'
      n1_f <- 1-t
      n2_f <- NaN
    } else if ((l2==0) && (t4>= t1)){
      t <- t4
      NodeColor[b==0]  <- 'Black'
      n1_f <- 1-t
      n2_f <- NaN
    } else if ((l2 != 0) && (l1 != 0) && (t1> t3)) {
      t <- t1
      NodeColor[b==0]  <- 'Red'
      NodeColor[b==1]  <- 'Black'
      n1_f <- 1-t
      n2_f <- 1-t2
    } else if ((l2 != 0) && (l1 != 0) && (t3> t1)) {
      t <- t3
      NodeColor[b==0]  <- 'Black'
      NodeColor[b==1]  <- 'Red'
      n1_f <- 1-t
      n2_f <- 1-t4
    }
    
    if (!is.na(t1) && (!is.na(t3))){
      if ((t1 ==0) && (t1 == t3)) {
        if ( sum(t(m2) * t(b)) >= sum(t(m2) * t(!b))){
          NodeColor[b==0] <- 'Red'
          NodeColor[b==1] <- 'Black'
          n1_f <- 1-t1
          n2_f <- 1-t2
        } else if ( sum(t(m2) * t(!b)) > sum(t(m2) * t(b))){
          NodeColor[b==0] <- 'Black'
          NodeColor[b==1] <- 'Red'
          n1_f <- 1-t3
          n2_f <- 1-t4
        }
      }
    }
    
    if (!is.na(t2) && (!is.na(t4))){
      if ((t2 ==0) && (t2 == t4)) {
        if ( sum(t(m1) * t(!b)) >= sum(t(m1) * t(b))){
          NodeColor[b==0] <- 'Red'
          NodeColor[b==1] <- 'Black'
          n1_f <- 1-t1
          n2_f <- 1-t2
        } else if ( sum(t(m1) * t(b)) > sum(t(m1) * t(!b))) {
          NodeColor[b==0] <- 'Black'
          NodeColor[b==1] <- 'Red'
          n1_f <- 1-t3
          n2_f <- 1-t4
        }
      }
    }
    
    return(list(NodeColor=NodeColor,n1_f=n1_f,n2_f=n2_f))
  }
  
  
  if (length(cutoff) == 1){
    pc_corr_res <- C_corr(norm[[u_norm]],pc[[u_norm]][,u_dim],feat_names,cutoff)
    pc_corr <- pc_corr_res$pc_corr
    x1 <- pc_corr_res$x1
    Edges <- pc_corr_res$Edges
    samplCol <- match_V_samplCol(col,pc_corr_res$x1,labels,nameLabels,pc_corr_res$Nodes)
    n1_f <- samplCol$n1_f
    n2_f <- samplCol$n2_f
    Nodes <- data.frame(pc_corr_res$Nodes[,1], samplCol$NodeColor,pc_corr_res$Nodes[,2])
    colnames(Nodes) <- c("Node","Colour","Loading (V)")
    
    edges <- Edges
    nodes <- Nodes
    # print(edges)
    # print(nodes)
  } else {
    Edges <- matrix(list(),length(cutoff),2)
    Nodes <- matrix(list(),length(cutoff),2)
    pc_corr <- matrix(list(),length(cutoff),2)
    x2 <- matrix(list(),length(cutoff),2)
    for (i in 1:length(cutoff)){
      Edges[[i,1]] <- cutoff[i]
      Nodes[[i,1]] <- cutoff[i]
      pc_corr[[i,1]] <- cutoff[i]
      x2[[i,1]] <- cutoff[i]
      
      pc_corr_res <- C_corr(norm[[u_norm]],pc[[u_norm]][,u_dim],feat_names,cutoff[i])
      Edges[[i,2]] <- pc_corr_res$Edges
      pc_corr[[i,2]] <- pc_corr_res$pc_corr
      x2[[i,2]] <- pc_corr_res$x1
       
      samplCol <- match_V_samplCol(col,pc_corr_res$x1,labels,nameLabels,pc_corr_res$Nodes)
      
      Nodes[[i,2]] <- data.frame(pc_corr_res$Nodes[,1], samplCol$NodeColor,pc_corr_res$Nodes[,2])
      colnames(Nodes[[i,2]]) <- c("Node","Colour","Loading (V)")
      
      
      cutoff1 <-  cutoff[i]
      edges <- Edges[[i,2]]
      nodes <- Nodes[[i,2]]
      # print(cutoff1)
      # print(edges)
      # print(nodes)
    }
  }

  filename <- 'PC-corr_net_edges-nodes_results.xlsx'
  if (file.exists(filename)){
    a_xlsx <- file.remove(filename)
  } 
  
  
  if (length(cutoff) == 1){
    write.xlsx(Edges, file=filename, sheetName ='Edges', row.names = FALSE,append=TRUE)
    write.xlsx(Nodes, file=filename, sheetName ='Nodes', row.names = FALSE,append=TRUE)
  } else {
    for (i in 1:length(cutoff)){
      write.xlsx(Edges[[i,2]], file=filename, sheetName = paste(c('Edges-cutoff ',toString(cutoff[i])),collapse=""), row.names = FALSE, append=TRUE)
      write.xlsx(Nodes[[i,2]], file=filename, sheetName = paste(c('Nodes-cutoff ',toString(cutoff[i])),collapse=""), row.names = FALSE, append=TRUE)
      
    }
  }  
  
  
  
  plot_graph <- function(pc_corr1,nodes1,cutoff){
    library(igraph)
    nodes11 <- unlist(lapply(nodes1[,1],as.character)) 
    pc_corr1 <- data.frame(pc_corr1)
    colnames(pc_corr1) <- nodes11
    rownames(pc_corr1) <- nodes11
    pc_corr1 <- as.matrix(pc_corr1)
    g <- graph_from_adjacency_matrix(pc_corr1, weighted=TRUE, mode="undirected")
    for (i in 1:length(V(g))){
      if (nodes1[i,2] == 'Red'){
        V(g)$color[i] <- rgb(1,0,0)
        V(g)$frame.color[i]<- rgb(1,0,0)
      } else {
        V(g)$color[i] <- rgb(0,0,0)
        V(g)$frame.color[i] <- rgb(0,0,0)
      }
    }
    
    val_edg_col <- matrix(nrow = length(E(g)),ncol = 3 )
    for (i in 1:length(E(g))){
      max_edg_w <- max(E(g)$weight)
      min_edg_w <- min(E(g)$weight)
      if (sign(E(g)[i]$weight)==1) {
        color_edg <- 1 - E(g)[i]$weight/max_edg_w
        E(g)[i]$color <- rgb(1,color_edg,color_edg)
        val_edg_col[i,] <- c(1,color_edg,color_edg)
      } else {
        color_edg <- E(g)[i]$weight/abs(min_edg_w) + 1
        E(g)[i]$color <- rgb(color_edg,color_edg,color_edg)
        val_edg_col[i,] <- c(color_edg,color_edg,color_edg)
      }
    }
    
    
    E(g)$width <- 2
    E(g)$lty <- 1
    #Replace the grey edges (edges under frustration) in the figure with dashed grey edges 
    edg_frust <- 0
    for (i in 1:(length(V(g))-1)){
      for (j in (i+1):length(V(g))){
        eh1 <- E(g)[i %--% j]
        if (length(eh1)==0){
          next
        } else{
          n <- (V(g)[i]$color == "#FF0000") # rgb of red is "#FF0000"
          m <- (V(g)[j]$color == "#FF0000") # rgb of red is "#FF0000"
          if (((n*m) == 1) & ( sign(eh1$weight) == -1)){
            E(g)[i %--% j]$color <- rgb(0.6,0.6,0.6)
            E(g)[i %--% j]$lty <- 2
            edg_frust <- edg_frust + 1
          } else if ((n+m == 0)& ( sign(eh1$weight) == -1)){
            E(g)[i %--% j]$color <- rgb(0.6,0.6,0.6)
            E(g)[i %--% j]$lty <- 2
            edg_frust <- edg_frust + 1
          } else if ((n+m == 1)& ( sign(eh1$weight) == 1)){
            E(g)[i %--% j]$color <- rgb(0.6,0.6,0.6)
            E(g)[i %--% j]$lty <- 2
            edg_frust <- edg_frust + 1
          }
        }
      }
    }
    
    
    val_red <- c()
    val_black <- c()
    for (i in 1:length(E(g))){
      if ((val_edg_col[i,1]==1)&(val_edg_col[i,2]!=0)){
        val_red <- c(val_red, val_edg_col[i,2])
      } else {
        val_red <- c(val_red)
      }
      if ((val_edg_col[i,1]!=1) & (val_edg_col[i,1]!=0.9)){
        val_black <- c(val_black,val_edg_col[i,1])
      } else{
        val_black <- c(val_black)
      }
    }
    
    vector.is.empty <- function(x) return(length(x) ==0 )
    if (vector.is.empty(val_red)){
      min_red <- c()
    }else {
      min_red <- min(val_red)
    }
    
    if (vector.is.empty(val_black)){
      min_black <- c()
    }else {
      min_black <- min(val_black)
    }
    
    
    fr_edg_frust <- edg_frust/length(E(g))#Percentage of edges under frustation
    X11(height=23,width=44, pointsize=12) 
    plot.igraph(g,vertex.shape="sphere", vertex.size=11, vertex.label.color=rgb(0.5,0.5,0.57), vertex.label.cex=.8, edge.curved=0)
    usr <- par( "usr" )
    fr_edg_frust_rd <- round( fr_edg_frust,3)*100
    nb_edges <- length(E(g))
    frust_txt1 <- do.call("expression", list(substitute("cutoff "== cutoff, list(cutoff=cutoff))))
    frust_txt2 <- do.call("expression", list(substitute("frustration "== frac(edg_frust,nb_edges) ,list(edg_frust=edg_frust,nb_edges=nb_edges))))
    frust_txt3 <- do.call("expression", list(substitute(" "== fr_edg_frust_rd~"%", list(fr_edg_frust_rd=fr_edg_frust_rd))))
    text( usr[ 1 ]-0.17, usr[ 4 ]+0.2,frust_txt1, adj = c( 0, 1 ) )
    text( usr[ 1 ]-0.17, usr[ 4 ]+0.1, frust_txt2, adj = c( 0, 1 ) )
    text( usr[ 1 ]+0.24, usr[ 4 ]+0.055, frust_txt3, adj = c( 0, 1 ) )
    
    plotHandle <- ((E(g)$color == rgb(0.6,0.6,0.6)) & (E(g)$lty == 2) )
    plotHandle1 <- ((sign(E(g)$weight)== 1) & (E(g)$color != rgb(0.6,0.6,0.6)) & (E(g)$lty != 2)) 
    plotHandle2 <- ((sign(E(g)$weight)== -1)& (E(g)$color != rgb(0.6,0.6,0.6))& (E(g)$lty != 2))
    plotHandle3 <- (V(g)$color == "#FF0000" )
    plotHandle4 <- (V(g)$color == "#000000" )
    
    plotHandle_1 <- sum(plotHandle==TRUE)
    plotHandle1_1 <- sum(plotHandle1==TRUE)
    plotHandle2_1 <- sum(plotHandle2==TRUE)
    plotHandle3_1 <- sum(plotHandle3==TRUE)
    plotHandle4_1 <- sum(plotHandle4==TRUE)
    
    plotHandle1m <- ((sign(E(g)$weight)== 1) & (E(g)$color != rgb(0.6,0.6,0.6)) & (E(g)$lty != 2)) & (E(g)$color == rgb(1,0,0))
    plotHandle1_1m <- sum(plotHandle1m==TRUE)  
    if (plotHandle1_1 == plotHandle1_1m){
      min_red <- 0
    } else {
      min_red <- min_red
    }
    
    plotHandle2m <- ((sign(E(g)$weight)== -1) & (E(g)$color != rgb(0.6,0.6,0.6)) & (E(g)$lty != 2)) & (E(g)$color == rgb(0,0,0))
    plotHandle2_1m <- sum(plotHandle2m==TRUE)  
    if (plotHandle2_1 == plotHandle2_1m){
      min_black <- 0
    } else {
      min_black <- min_black
    }
    
    
    if ((plotHandle1_1 ==0) & (plotHandle2_1 ==0) & (plotHandle3_1 !=0) & (plotHandle4_1 !=0) & (plotHandle_1 !=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('Edge under frustation', '\u2191 Red sample group', '\u2191 Black sample group'),lty=c(2,NA,NA),pch=c(NA,19,19),col=c(rgb(.6,.6,.6),"red","black"),cex=1.2)
    } else if ((plotHandle1_1 ==0) & (plotHandle2_1 ==0) & (plotHandle3_1 ==0) & (plotHandle4_1 !=0) & (plotHandle_1 !=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('Edge under frustation','\u2191 Black sample group'),lty=c(2,NA),pch=c(NA,19),col=c(rgb(.6,.6,.6),"black"),cex=1.2)
    } else if ((plotHandle1_1 ==0) & (plotHandle2_1 ==0) & (plotHandle3_1 !=0) & (plotHandle4_1 ==0) & (plotHandle_1 !=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('Edge under frustation','\u2191 Red sample group'),lty=c(2,NA),pch=c(NA,19),col=c(rgb(.6,.6,.6),"red"),cex=1.2)
    } else if ((plotHandle4_1==0) & (plotHandle2_1==0) & (plotHandle_1!=0) & (plotHandle1_1!=0) & (plotHandle3_1!=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','Edge under frustation', '\u2191 Red sample group'),lty=c(1,2,NA),pch=c(NA,NA,19),col=c(rgb(1,min_red,min_red),rgb(.6,.6,.6),"red"),cex=1.2)
    } else if (( plotHandle1_1==0) & (plotHandle3_1==0) & ( plotHandle_1!=0) & (plotHandle2_1!=0) & (plotHandle4_1!=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr < 0','Edge under frustation', '\u2191 Black sample group'),lty=c(1,2,NA),pch=c(NA,NA,19),col=c(rgb(min_black,min_black,min_black),rgb(.6,.6,.6),"black"),cex=1.2)
    } else if ((plotHandle2_1==0) & (plotHandle3_1!=0) & (plotHandle_1 !=0) & (plotHandle1_1!=0) & (plotHandle4_1!=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','Edge under frustation', '\u2191 Red sample group','\u2191 Black sample group'),lty=c(1,2,NA,NA),pch=c(NA,NA,19,19),col=c(rgb(1,min_red,min_red),rgb(.6,.6,.6),"red","black"),cex=1.2)
    } else if (( plotHandle1_1==0) & (plotHandle3_1!=0) & (plotHandle_1 !=0) & (plotHandle2_1!=0) & (plotHandle4_1!=0)){
      legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr < 0','Edge under frustation', '\u2191 Red sample group','\u2191 Black sample group'),lty=c(1,2,NA,NA),pch=c(NA,NA,19,19),col=c(rgb(min_black,min_black,min_black),rgb(.6,.6,.6),"red","black"),cex=1.2)
    } else if ((plotHandle3_1==0) & (plotHandle2_1==0) & (plotHandle1_1 !=0) & (plotHandle4_1 !=0) ){
      if (plotHandle_1 ==0){
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0', '\u2191 Black sample group'),lty=c(1,NA),pch=c(NA,19),col=c(rgb(1,min_red,min_red),"black"),cex=1) #only black nodes
      } else {
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','Edge under frustation', '\u2191 Black sample group'),lty=c(1,2,NA),pch=c(NA,NA,19),col=c(rgb(1,min_red,min_red),rgb(.6,.6,.6),"black"),cex=1.2) #only black nodes
      }
    } else if ((plotHandle4_1==0) & (plotHandle2_1==0) & (plotHandle1_1 !=0) & (plotHandle3_1 !=0)) {
      if (plotHandle_1 ==0){
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0', '\u2191 Red sample group'),lty=c(1,NA),pch=c(NA,19),col=c(rgb(1,min_red,min_red),"red"),cex=1.2) #only red nodes
      } else {
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','Edge under frustation', '\u2191 Red sample group'),lty=c(1,2,NA),pch=c(NA,NA,19),col=c(rgb(1,min_red,min_red),rgb(.6,.6,.6),"red"),cex=1.2) #only red nodes
      }
    } else if ((plotHandle2_1==0)& (plotHandle4_1!=0) & (plotHandle1_1 !=0) & (plotHandle3_1 !=0)) {
      if (plotHandle_1 ==0){
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0', '\u2191 Red sample group', '\u2191 Black sample group'),lty=c(1,NA,NA),pch=c(NA,19,19),col=c(rgb(1,min_red,min_red),"red","black"),cex=1.2) #only positive PC-corr
      } else {
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','Edge under frustation', '\u2191 Red sample group', '\u2191 Black sample group'),lty=c(1,2,NA,NA),pch=c(NA,NA,19,19),col=c(rgb(1,min_red,min_red),rgb(.6,.6,.6),"red","black"),cex=1.2) #only positive PC-corr
      }
    } else if ((plotHandle1_1==0)& (plotHandle2_1!=0) & (plotHandle4_1 !=0) & (plotHandle3_1 !=0)){
      if (plotHandle_1 ==0){
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr < 0', '\u2191 Red sample group', '\u2191 Black sample group'),lty=c(1,NA,NA),pch=c(NA,19,19),col=c(rgb(min_black,min_black,min_black),"red","black"),cex=1.2) #only negative PC-corr
      } else {
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr < 0','Edge under frustation', '\u2191 Red sample group', '\u2191 Black sample group'),lty=c(1,2,NA,NA),pch=c(NA,NA,19,19),col=c(rgb(min_black,min_black,min_black),rgb(.6,.6,.6),"red","black"),cex=1.2) #only negative PC-corr
      }
    } else if (( plotHandle1_1!=0) & (plotHandle2_1!=0) & (plotHandle3_1!=0) & (plotHandle4_1!=0)){
      if (plotHandle_1 ==0){
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','PC-corr < 0','\u2191 Red sample group', '\u2191 Black sample group'),lty=c(1,1,NA,NA),pch=c(NA,NA,19,19),col=c(rgb(1,min_red,min_red),rgb(min_black,min_black,min_black),"red","black"),cex=1.2) #no edges under frustation
      } else {
        legend(usr[2]-1,usr[ 4 ]+0.2,c('PC-corr > 0','PC-corr < 0','Edge under frustation','\u2191 Red sample group', '\u2191 Black sample group'),lty=c(1,1,2,NA,NA),pch=c(NA,NA,NA,19,19),col=c(rgb(1,min_red,min_red),rgb(min_black,min_black,min_black),rgb(.6,.6,.6),"red","black"),cex=1.2) 
      }
    }
    

    #tkplot(g)
  }
  
  
  #Graph plot
  for (k in 1:length(cutoff)){
    if (length(cutoff)==1){
      pc_corr1 <- pc_corr
      nodes1 <- Nodes
    } else {
      pc_corr1 <- pc_corr[[k,2]]
      nodes1 <- Nodes[[k,2]]
    }
    if ((dim(pc_corr1)[1] == 0) & (dim(pc_corr1)[2] == 0)){
      next
    } else{
      plot_graph(pc_corr1,nodes1,cutoff[k])
    }
  }

  list2env(list(Edges=Edges,Nodes=Nodes),.GlobalEnv)
  
}
  
  
  





































