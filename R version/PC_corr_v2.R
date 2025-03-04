###### PC-corr algorithm associate to PCA analysis #####

### Released under MIT License
### Copyright (c) 16 Dec 2017 Sara Ciucci, Yan Ge, Claudio Dur?n and Carlo Vittorio Cannistraci

# Please cite:
# Enlightening discriminative network functional modules behind Principal Component Analysis separation in differential-omic science studies.
# Sara Ciucci, Yan Ge, Claudio Dur?n, Alessandra Palladini, V?ctor Jim?nez Jim?nez, Luisa Mar?a Mart?nez S?nchez, 
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

PC_corr_v2<-function(x,sample_labels,feat_names, sample_names,dis) {
  

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
  library(tictoc)
  
  cat("\f")  
  
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
  
  cat("\nWelcome to PC-corr, a simple algorithm that associates to any PCA segregation a discriminative network of features.\n")
  cat("For more details, see https://www.nature.com/articles/srep43946\n\n")
  
  #  ------------------------------------------------------------------------
  
  labels <- sample_labels
  nameLabels <- unique(labels) # name of the groups
  numbLabels <- length(nameLabels) # number of groups
  
  # Remove features with same identical values across all the samples
  x1 <- x 
  x1 <- x1-matrix(rep(colMeans(x1),each=nrow(x1)),nrow=dim(x1)[1],ncol=dim(x1)[2])
  ma <- colSums(x1==0)==nrow(x1)
  remov_feat <- sum(ma) #number of removed fatures
  x <- x[,!ma]
  feat_names<-feat_names[t(!ma)]
  
  #  user interaction: type of labels ------------------------------------------------------------------------
  u_lab <- c( )
  flag <- 0
  while (flag == 0){
    cat("\nIs your data represented by\n[r] ranked labels(labels that are organized according to a progressive order. e.g. different stages of a disease, where Stage 1 < Stage 2 < Stage 3)\n[c] class labels (labels that are not necessary organized in a progressive order e.g. Condition A, Condition B, Condition C)? [r/c]:\n\n")
    u_lab <-readline(prompt="-> ")
    flag <- 1
    if ((u_lab != 'r') & (u_lab != 'c')){
      flag <- 0
      cat("\nPlease introduce either 'r' for ranked labels or 'c' for class labels\n")
    }
  }
  
  if (u_lab =='r') {
    flag <- 0 
    while (flag == 0){
      cat("\nAre the values of your ranked labels\n[d]   discrete (Stage 1 < Stage 2 < Stage 3)\n[con] continuous (different times of development of a cell line)? [d/con]:\n\n")
      u_lab <- readline(prompt="-> ")
      flag <- 1
      if ((u_lab != 'd') & (u_lab != 'con')){
        flag <- 0
        cat("\nPlease introduce either 'd' for discrete labels or 'con' for continuous labels\n")
      }
    }
  }
  
  labl_numb <- c()
  if (u_lab == 'd') {
    #for correlation evaluators, labels are turned into numbers
    for (i in 1:numbLabels){
      labl_numb[labels %in% nameLabels[i]] <- i
    } 
  } else if (u_lab == 'con'){
    labl_numb <- as.numeric(labels)
  }
  
  # Normalizations of the dataset ------------------------------------------- 
  
  flag <- 0
  while ( flag == 0){
    cat("\nThe analysis starts by normalizing or not the dataset.\n\nDo you want to apply:\n[1] no normalization \n[2] a preferred normalization \n[3] automatically all the set of available normalizations? [1/2/3]\n\n");
    u_norm_opt <- readline(prompt="-> ")
    
    flag <- 1
    
    if ((u_norm_opt !=1) & (u_norm_opt !=2) & (u_norm_opt !=3)){
      flag <- 0
      cat("\nPlease introduce either 1 for no normalization, 2 for a preferred normalization or 3 for all the set of available normalizations. \n")
    } else {
      flag <- 1
    }
  }
  
  
  if (u_norm_opt == 1) {
    norm <- matrix(list(),1,1)
    norms <-  character()
    norm[[1]] <- x  #No normalization
    norms[1] <- '-' 
  } else if (u_norm_opt == 3){
    norm <- matrix(list(),12,1)
    norms <-  character()
    norm[[1]] <- x/matrix(rep(colSums(x),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2]) #dividing by the column sum
    norms[1] <- 'DCS'
    
    norm[[2]] <- x/matrix(rep(rowSums(x),each = ncol(x)),nrow = dim(x)[1],ncol = dim(x)[2],byrow=TRUE) #dividing by the row sum
    norms[2] <- 'DRS'
    
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
  } else if (u_norm_opt == 2){
     
    flag <- 0
    while (flag == 0){
      norms_list_1 <- c('DCS','DRS','LOG','ZSCORE','QUANTILE T','QUANTILE','ZSCORE T','PLUS(ABS(MIN))','PARETO SCALING','SQRT','MANORM')
      cat("\nInput a preferred normalization from the following list:\n")
      cat("DCS, DRS, LOG, ZSCORE, QUANTILE T, QUANTILE, ZSCORE T, PLUS(ABS(MIN)), PARETO SCALING, SQRT, MANORM\n")
      cat("(For detailed information on the type of normalization, see the User guide)\n\n")
      cat("Example: LOG\n\n")
      u_norm_choice <- readline(prompt="-> ")
      flag <- 1
      u_norm_choice1 <- which(norms_list_1 %in% u_norm_choice)
      if (length(u_norm_choice1) == 0){
        flag <- 0
        cat("Please introduce the exact name of the normalization.\n")
      }
    }
    
    norm <- matrix(list(),1,1)
    norms <-  character()
    norms[1] <- u_norm_choice
    if (u_norm_choice == 'DCS'){
      norm[[1]] <- x/matrix(rep(colSums(x),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2]) #dividing by the column sum
    } else if (u_norm_choice == 'DRS'){
      norm[[1]] <- x/matrix(rep(rowSums(x),each = ncol(x)),nrow = dim(x)[1],ncol = dim(x)[2],byrow=TRUE) #dividing by the row sum
    } else if (u_norm_choice == 'LOG'){
      norm[[1]] <- log10(1+x)
    } else if (u_norm_choice == 'ZSCORE'){
      norm[[1]] <- scale(x)
    } else if (u_norm_choice == 'QUANTILE T'){
      norm[[1]] <- t(normalize.quantiles(t(x)))
    } else if (u_norm_choice == 'QUANTILE'){
      norm[[1]] <- normalize.quantiles(x)
    } else if (u_norm_choice == 'ZSCORE T'){
      norm[[1]] <- t(scale(t(x)))
    } else if (u_norm_choice == 'PLUS(ABS(MIN))'){
      norm[[1]] <- x + abs(min(x))
    } else if (u_norm_choice == 'PARETO SCALING'){
      x_center <- x-matrix(rep(colMeans(x),each=nrow(x)),nrow=dim(x)[1],ncol=dim(x)[2])
      norm[[1]] <- x_center/matrix(rep(sqrt(apply(x_center, 2, sd)),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2])
    } else if (u_norm_choice == 'SQRT'){
      norm[[1]] <- sqrt(x)
    } else if (u_norm_choice == 'MANORM'){
      norm[[1]] <- x/matrix(rep(colMeans(x),each = nrow(x)),nrow = dim(x)[1],ncol = dim(x)[2])
    } 
  }
  

  inf_norm <- c()
  nan_norm <- c()
  nonreal_norm <- c()
  for (i in 1:length(norm)){
    inf_norm[i] <- sum(is.infinite(norm[[i]]))
    nan_norm[i] <- sum(is.nan(norm[[i]]))
    nonreal_norm[i] <- is.complex(norm[[i]])
  }
  
  probl_norm <- inf_norm + nan_norm + nonreal_norm
  norm <- norm[!probl_norm]
  norms <- norms[!probl_norm]
  norms_list <- norms
 
  

  
  # number of elements in each group 
  number_el_group <-c()
  for (i in 1:length(nameLabels)){
    number_el_group[i] <- sum(labels==nameLabels[i])
  }
 
  

  flag <- 0 
  u_aupr <- c()
  
  if ((sum(duplicated(number_el_group))==0) & (numbLabels >2 )){
    while (flag == 0) {

      cat("\nFor the calculation of Area Under the ROC-Curve (AUC) and Area Under the Precision-Recall curve (AUPR), you need to provide a positive label for each pairwise group comparison.\n")
      cat("\nDo you want to calculate the AUC and AUPR values considering:\n[s] as positive label, the label of the smallest sample group in each pairwise group comparison\n[l] as positive label, the label of the largest sample group in each pairwise group comparison\n[r] a ranked list of possible positive labels \n\n")
      u_aupr <- readline(prompt="-> ")
      
      flag <- 1
      
      if ((u_aupr != 's') & (u_aupr != 'l') & (u_aupr != 'r')){
        flag <- 0
        
        cat("\nPlease introduce either 's' for the label of the smallest sample group as positive label, \n 'l' for the labels of the largest sample group as positive label, 'r' for a ranked list of possible positive labels.\n")
        
      } 
    }
  } else {
    flag <- 1
    u_aupr <- 'r'
  }
  
  u_aupr_r <- c()
  if (u_aupr == 'r'){
    flag <- 0
    while (flag == 0){
      if (numbLabels ==2){
        cat("\nInput the positive label for the calculation of AUC and AUPR values: \n");
      }else {
        cat("\nInput the ranked list of possible positive labels for the calculation of AUC and AUPR values: \n");
      }
      randperm_nameLabels <- sample(nameLabels, replace=FALSE)
      if (length(nameLabels)>2){
        cat("Example: c(",paste(encodeString(randperm_nameLabels[1:(length(nameLabels)-1)], quote = '"'),collapse = ","),")\n\n",sep="")
      } else if(length(nameLabels)==2) {
        cat("Example: ",encodeString(randperm_nameLabels[1], quote = '"'),"\n\n",sep="")
      }
      u_aupr_r <- readline(prompt="-> ")
      u_aupr_r <- eval(parse(text=u_aupr_r))
      flag <- 1
      if (sum(is.element(u_aupr_r,nameLabels))==0){
        flag <- 0
        if (numbLabels >2){
          cat("\nPlease introduce a correct ranked list of positive labels.\n")
        } else if (numbLabels == 2){
          cat("\nPlease introduce a correct positive label.\n")
        }
      }
    }
  }
  
  if ((u_aupr == 's')|(u_aupr == 'l')){
    u_aupr_r <- c()
  }
  
  # positive label depending on the chosen option-----------------------------------------------------------------
  positive_label_opt <- function(u_aupr,nameLabels1,nameLabels2,labels,u_aupr_r){
    
    if (u_aupr == 's'){
      len_n <- length(labels[labels==nameLabels1])
      len_m <- length(labels[labels==nameLabels2])
      if (len_n < len_m){
        possClass <- nameLabels1
      } else if (len_m < len_n){
        possClass <- nameLabels2
      }
    } else if (u_aupr == 'l'){
      len_n <- length(labels[labels==nameLabels1])
      len_m <- length(labels[labels==nameLabels2])
      if (len_n < len_m){
        possClass <- nameLabels2
      } else if (len_m < len_n){
        possClass <- nameLabels1
      }
    } else if (u_aupr == 'r'){
      idx_n <- which(u_aupr_r %in% nameLabels1) 
      idx_m <- which(u_aupr_r %in% nameLabels2) 
      if (length(idx_m)==0){
        possClass <- nameLabels1
      } else if (length(idx_n)==0){
        possClass <- nameLabels2
      } else {
        if (idx_n < idx_m){
          possClass <- nameLabels1
        } else if (idx_m < idx_n){
          possClass <- nameLabels2
        }
      }

    }
    return(possClass)
  }
  
  #  aupr evaluation ------------------------------------------------------------------------
  
  aupr_evaluation <- function(samp_lab,scores, possClass,negClass){
    pred<- prediction(scores,samp_lab,label.ordering=c(negClass,possClass))
    RP.perf <- performance(pred, "prec", "rec")
    recall <- unlist(RP.perf@x.values)[1:length(unlist(RP.perf@x.values))]
    precision <- unlist(RP.perf@y.values)[1:length(unlist(RP.perf@y.values))]
    if (precision[2] ==1){
      precision[1] <- 1
    } else {
      precision[1] <- 0
    }
    AUPR<- trapz(recall,precision)
    return(AUPR)
  }    
  
  #  ------------------------------------------------------------------------
  if (length(labels) <= length(feat_names)){
    dimension_PCA <- length(labels)
  } else {
    dimension_PCA <- length(feat_names)
  }
  pc_nc<- matrix(list(),length(norm),1)
  pc_c<- matrix(list(),length(norm),1)
  ncPCA<- pc_nc
  cPCA<- pc_c
  explained_nc <- matrix(list(),length(norm),1)
  explained_c <- matrix(list(),length(norm),1)
  mw_ncPCA <- array(NA,dim=c(length(norm),choose(numbLabels,2),dimension_PCA))
  mw_cPCA <- mw_ncPCA
  AUC_nc <- array(NA,dim=c(length(norm),choose(numbLabels,2),dimension_PCA))
  AUC_c <- AUC_nc
  response <- c()
  AUPR_nc <- array(NA,dim=c(length(norm),choose(numbLabels,2),dimension_PCA))
  AUPR_c <- AUPR_nc
  rank_pears_corr_ncPCA <- array(NA,dim=c(length(norm),dimension_PCA,1))
  rank_pears_corr_cPCA <- rank_pears_corr_ncPCA
  rank_spear_corr_ncPCA <- array(NA,dim=c(length(norm),dimension_PCA,1))
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
          
          samp_lab <- c(labels[is.element(labels,nameLabels[n])],labels[is.element(labels,nameLabels[m])])
          scores_nc <- c(ncPCA[[i]][is.element(labels,nameLabels[n]),k],ncPCA[[i]][is.element(labels,nameLabels[m]),k])
          scores_c <- c(cPCA[[i]][is.element(labels,nameLabels[n]),k],cPCA[[i]][is.element(labels,nameLabels[m]),k])
          
          possClass <- positive_label_opt(u_aupr, nameLabels[n],nameLabels[m],labels,u_aupr_r)
          negClass <- c(nameLabels[n],nameLabels[m])
          negClass <- negClass[negClass!=possClass]
            
          response <- c()
          # Compute AUC
          response[samp_lab == possClass] <- 1
          response[samp_lab == negClass] <- 0
          
          res_ROC_nc <- roc(response,scores_nc,direction = "<")
          res_ROC_c <- roc(response,scores_c,direction = "<")
          AUC_nc[i,j,k] <- res_ROC_nc$auc[1]
          AUC_c[i,j,k] <- res_ROC_c$auc[1]
          
          
          if (AUC_nc[i,j,k] < 0.5){
 
            AUC_nc[i,j,k] <- 1-AUC_nc[i,j,k]
            #Compute AUPR
            flip_scores_nc <- 2*mean(scores_nc)-scores_nc
            AUPR_nc[i,j,k] <- aupr_evaluation(samp_lab,flip_scores_nc, possClass,negClass)
          } else {
            #Compute AUPR
            AUPR_nc[i,j,k] <- aupr_evaluation(samp_lab,scores_nc, possClass,negClass)
          }

          
          if (AUC_c[i,j,k] < 0.5){
            
            AUC_c[i,j,k] <- 1-AUC_c[i,j,k]
            #Compute AUPR
            flip_scores_c <- 2*mean(scores_c)-scores_c
            AUPR_c[i,j,k] <- aupr_evaluation(samp_lab,flip_scores_c, possClass,negClass)
          } else {
            # Compute AUPR
            AUPR_c[i,j,k] <- aupr_evaluation(samp_lab,scores_c, possClass,negClass)
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
  
  # Certain number of permutations of the labels ---------------------------------------
  random_permutation_labels <- function(numb_rand,labels){
    # This function returns numb_rand random permutation of the labels
    #
    # INPUT
    # numb_rand => (integer number) number of permutations 
    # labels => (Character vector Mx1) labels of the samples 
    #
    # OUTPUT
    # labels_rp => (matrix of list numb_rand x 1) numb_rand random permutation of the sample labels:
    #              each element of the list contains a character vector of the random permutation of the sample labels 
    labels_rp <- matrix(list(),numb_rand,1)
    for (i in 1:numb_rand){
      labels_rp[[i]] <- sample(sample_labels,replace=FALSE)
    }
    
    return(labels_rp)
  }
  
  # Significance of segragation - Trustworthiness----------------------------
  significance_segragation <- function(numb_rand,PCA,dim,labels_rp,nameLabels,numbLabels,segr_meas,type,u_aupr,u_aupr_r){
    # This function evaluates the significane of the segragation
    # measure (P-value, AUC or AUPR). Is the dimension discriminative
    # because it is random or capture the main variability?
      
    # INPUT
    #   numb_rand => number of permutation
    #   PCA=> PCA results either centred or not centred, with a certain
    #   normalization or not
    #   dim => dimension of sample segragation
    #   labels_rp => numb_rand random permutation of the sample labels
    #   nameLabels => unique name of labels
    #   numbLabels => number of unique sample labels
    #   segr_meas => segragation measure on the
    #   input dimension dim:
    #   - for 'p','auc','aupr', just one value: either segragation
    #   measure according to P-value, AUC or AUPR
    #   - for 'all': segragation measures according to P-value (in pos. 1), 
    #   AUC (in pos.2)and AUPR (in pos. 3)
    #   type => type of segragation measure, depending on the input
    #   segr_meas: 'p' (for MW p-value), 'auc' (for AUC), 'aupr' (for
    #   AUPR) or 'all' (if you want to compute for all: MW p-value, AUC
    #   and AUPR)
    
    
    # OUTPUT
    #   p_Value => significance (p-value) of the input segragation
    #   measure:
    #   - a number if computed for just one input segragation measure
    #   (type= 'p', 'auc' or 'aupr')
    #   - a numeric vector if computed for just one input segragation measure
    #   (type= 'all'): in idx 1 for P-value, in 2 for AUC and in 3 for
    #   AUPR
    
    rand_segr_meas <- matrix(data=NA, nrow=choose(numbLabels,2),ncol=numb_rand)
    if (type =='all'){
      rand_segr_meas_mw <- matrix(data=NA, nrow=choose(numbLabels,2),ncol=numb_rand)
      rand_segr_meas_auc <- matrix(data=NA, nrow=choose(numbLabels,2),ncol=numb_rand)
      rand_segr_meas_aupr <- matrix(data=NA, nrow=choose(numbLabels,2),ncol=numb_rand)
    }
    for (pr in 1:numb_rand){
      
      labels <- labels_rp[[pr]]
      
      n <- 1
      m <- 2
      
      for (j in 1:choose(numbLabels,2)) { # two-group comparison
        
        if (type == 'p'){
          # Compute p-value of Mann-Whitney test
          rand_segr_meas[j,pr] <- wilcox.test(PCA[is.element(labels,nameLabels[n]),dim],PCA[is.element(labels,nameLabels[m]),dim],paired=FALSE)$p.value 
        } else if (type == 'auc'){
          samp_lab <- c(labels[is.element(labels,nameLabels[n])],labels[is.element(labels,nameLabels[m])])
          scores <- c(PCA[is.element(labels,nameLabels[n]),dim],PCA[is.element(labels,nameLabels[m]),dim])

          possClass <- positive_label_opt(u_aupr, nameLabels[n],nameLabels[m],labels,u_aupr_r)
          negClass <- c(nameLabels[n],nameLabels[m])
          negClass <- negClass[negClass!=possClass]
          
          response <- c()
          # Compute AUC
          response[samp_lab == possClass] <- 1
          response[samp_lab == negClass] <- 0
          
          res_ROC <- roc(response,scores,direction = "<")
          rand_segr_meas[j,pr]<- res_ROC$auc[1]
          if (rand_segr_meas[j,pr] < 0.5){
            rand_segr_meas[j,pr] <- 1-rand_segr_meas[j,pr]
          }
        } else if (type =="aupr"){
          samp_lab <- c(labels[is.element(labels,nameLabels[n])],labels[is.element(labels,nameLabels[m])])
          scores <- c(PCA[is.element(labels,nameLabels[n]),dim],PCA[is.element(labels,nameLabels[m]),dim])
          
          possClass <- positive_label_opt(u_aupr, nameLabels[n],nameLabels[m],labels,u_aupr_r)
          negClass <- c(nameLabels[n],nameLabels[m])
          negClass <- negClass[negClass!=possClass]
          
          response <- c()
          # Compute AUC
          response[samp_lab == possClass] <- 1
          response[samp_lab == negClass] <- 0
          
          res_ROC <- roc(response,scores,direction = "<")
          auc<- res_ROC$auc[1]
          if (auc < 0.5){
            #Compute AUPR
            flip_scores <- 2*mean(scores)-scores
            rand_segr_meas[j,pr] <- aupr_evaluation(samp_lab,flip_scores, possClass,negClass)
          } else {
            rand_segr_meas[j,pr] <- aupr_evaluation(samp_lab,scores, possClass,negClass)
          }
        } else if (type == 'all'){
          # Compute p-value of Mann-Whitney test
          rand_segr_meas_mw[j,pr] <- wilcox.test(PCA[is.element(labels,nameLabels[n]),dim],PCA[is.element(labels,nameLabels[m]),dim],paired=FALSE)$p.value 
          
          samp_lab <- c(labels[is.element(labels,nameLabels[n])],labels[is.element(labels,nameLabels[m])])
          scores <- c(PCA[is.element(labels,nameLabels[n]),dim],PCA[is.element(labels,nameLabels[m]),dim])
          
          possClass <- positive_label_opt(u_aupr, nameLabels[n],nameLabels[m],labels,u_aupr_r)
          negClass <- c(nameLabels[n],nameLabels[m])
          negClass <- negClass[negClass!=possClass]
          
          response <- c()
          # Compute AUC
          response[samp_lab == possClass] <- 1
          response[samp_lab == negClass] <- 0
          
          res_ROC <- roc(response,scores,direction = "<")
          rand_segr_meas_auc[j,pr]<- res_ROC$auc[1]
          if (rand_segr_meas_auc[j,pr] < 0.5){
            rand_segr_meas_auc[j,pr] <- 1- rand_segr_meas_auc[j,pr] 
            #Compute AUPR
            flip_scores <- 2*mean(scores)-scores
            rand_segr_meas_aupr[j,pr] <- aupr_evaluation(samp_lab,flip_scores, possClass,negClass)
          } else {
            rand_segr_meas_aupr[j,pr] <- aupr_evaluation(samp_lab,scores, possClass,negClass)
          }
        }
        
        m <- m + 1
        if (m > numbLabels){
          n <- n + 1
          m <- n + 1
        }
      }
    }
    
    if ((type == 'p') | (type == 'auc') | (type == 'aupr')){
          Rand_segr_meas <- colMeans(rand_segr_meas,1)
    } else if (type == 'all'){
      Rand_segr_meas_mw <- colMeans(rand_segr_meas_mw,1)
      Rand_segr_meas_auc <- colMeans(rand_segr_meas_auc,1)
      Rand_segr_meas_aupr <- colMeans(rand_segr_meas_aupr,1)
    }
    # Calculation of the significance (p-value) of the input segragation measure
    if (type == 'p'){
      p_Value <- (sum(Rand_segr_meas<segr_meas)+1)/(numb_rand+1)
      return(p_Value)
    } else if ((type =='auc')|(type=='aupr')){
      p_Value <- (sum(Rand_segr_meas>segr_meas)+1)/(numb_rand+1)
      return(p_Value)
    } else if (type == 'all') {
      p_Value <- c( )
      p_Value[1]<- (sum(Rand_segr_meas_mw<segr_meas[1])+1)/(numb_rand+1)
      p_Value[2] <- (sum(Rand_segr_meas_auc>segr_meas[2])+1)/(numb_rand+1)
      p_Value[3] <- (sum(Rand_segr_meas_aupr>segr_meas[3])+1)/(numb_rand+1)
      return(p_Value)
    }

  }
  
  
  # Constructing the table of results ---------------------------------------
  norms <- rep(norms,times=dim(ncPCA[[1]])[2]*2)
  centred <- rep('yes',times=dim(ncPCA[[1]])[2]*length(norm))
  non_centred <- rep('no',times=dim(ncPCA[[1]])[2]*length(norm))
  centr <- c(non_centred,centred)
  dim <- 1:dim(ncPCA[[1]])[2]
  dim <- rep(dim, each=length(norm))
  dim <- c(dim,dim)
  variance_c <- matrix(unlist(explained_c),nrow = length(norm), ncol = dim(ncPCA[[1]])[2] , byrow=TRUE)
  variance_c <- as.vector(variance_c)
  variance_nc <- matrix(unlist(explained_nc),nrow = length(norm), ncol = dim(ncPCA[[1]])[2] , byrow=TRUE)
  variance_nc <- as.vector(variance_nc)
  variance <- c(variance_nc,variance_c)
  
  
  if (numbLabels==2){
    if (u_lab=='c'){
      header <- c("P-value","AUC","AUPR","Norm","Centering","Dim","expl Var")
      header_xls <- c("P-value","AUC","AUPR","Norm","Centering","Dim","expl Var","Trustworthiness(p-value)","Trustworthiness(AUC)","Trustworthiness(AUPR)")
      pvals <- c(mw_ncPCA,mw_cPCA)
      AUCs <- c(AUC_nc,AUC_c)
      AUPRs <- c(AUPR_nc,AUPR_c)
      results <- data.frame(pvals,AUCs,AUPRs,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
      results_xls <- results
      
    } else if (u_lab=='d'){
      header <- c("P-value","AUC","AUPR", "pears","spear","Norm","Centering","Dim","expl Var")
      header_xls <- c("P-value","AUC","AUPR", "pears","spear","Norm","Centering","Dim","expl Var","Trustworthiness(p-value)","Trustworthiness(AUC)","Trustworthiness(AUPR)")
      pvals <- c(mw_ncPCA,mw_cPCA)
      AUCs <- c(AUC_nc,AUC_c)
      AUPRs <- c(AUPR_nc,AUPR_c)
      pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
      spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
      results <- data.frame(pvals,AUCs,AUPRs,pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
      results_xls <- results
      
    } else if (u_lab=='con'){
      header <- c("pears","spear","Norm","Centering","Dim","expl Var")
      pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
      spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
      results <- data.frame(pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      colnames(results) <- header
      results_xls <- results
    }
    
    cat("\nThe best discrimination in a PCA result (combination of normalization and centering) and along one dimension (principal component, PCn) can be assessed by different evaluators.\n");
    
    flag <- 0
    while (flag==0){
      if ((u_lab=='c')|(u_lab =='d')){
        if (u_lab=='c'){
          cat("\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR? [p/auc/aupr]:\n\n")
          u_rank <- readline(prompt="-> ")
        } else {
          cat("\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR\n[pc]   Pearson correlation\n[sc]   Spearman correlation? [p/auc/aupr/pc/sc]:\n\n")
          u_rank <- readline(prompt="-> ")
        }
        
        flag <- 1
        
        if (u_rank=='p'){
          idx <- order(results[,1],decreasing=FALSE)
          results <- results[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,1]<0.05,])==0){
            cat("\nThe table below shows the best results (p-value<0.25), automatically ranked with respect to P-value from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")  
            print(kable(format(results[results[,1]<0.25,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below shows the best results (p-value<0.05), automatically ranked with respect to P-value from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[results[,1]<0.05,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to P-value.\n\n")
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          results_xls <- results_xls[idx,]
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
          
        } else if (u_rank =='auc'){
          idx <- order(results[,2],decreasing=TRUE)
          results <- results[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,2]>=0.7,])==0){
            cat("\nThe table below shows the best results (AUC>=0.6), automatically ranked with respect to AUC from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[results[,2]>=0.6,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below shows the best results (AUC>=0.7), automatically ranked with respect to AUC from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[results[,2]>=0.7,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to AUC.\n\n")
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          results_xls <- results_xls[idx,]
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank =='aupr'){
          idx <- order(results[,3],decreasing=TRUE)
          results <- results[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,3]>=0.7,])==0){
            cat("\nThe table below shows the best results (AUPR>=0.6), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[results[,3]>=0.6,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below shows the best results (AUPR>=0.7), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[results[,3]>=0.7,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to AUPR.\n\n")
          #Creation of the excel file with all the results
          filename <-'results.xlsx'
          results_xls <- results_xls[idx,]
          if (file.exists(filename)){
            a_xlsx <- file.remove(filename)
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          } else {
            write.xlsx(results_xls, file=filename, sheetName ='PCA results', row.names = FALSE)
          }
        } else if (u_rank == 'pc'){
          if (u_lab=='c'){
            flag <- 0
            cat("Please introduce either 'p', 'auc' or 'aupr'\n")
          } else {
            idx <- order(abs(results[,4]),decreasing=TRUE)
            results <- results[idx,]
            #Only some results are shown on the screen
            if (nrow(results[abs(results[,4])>=0.6,])==0){
              cat("\nThe table below shows the best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
              print(kable(format(results[abs(results[,4])>=0.5,],digits=7), format="markdown",row.names=FALSE))
            } else {
              cat("\nThe table below shows the best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
              print(kable(format(results[abs(results[,4])>=0.6,],digits=7), format="markdown",row.names=FALSE))
            }
            cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Pearson correlation.\n\n")
            #Creation of the excel file with all the results
            filename <-'results.xlsx'
            results_xls <- results_xls[idx,]
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
            #Only some results are shown on the screen
            if (nrow(results[abs(results[,5])>=0.6,])==0){
              cat("\nThe table below shows the best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
              print(kable(format(results[abs(results[,5])>=0.5,],digits=7), format="markdown",row.names=FALSE))
            } else {
              cat("\nThe table below shows the best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
              print(kable(format(results[abs(results[,5])>=0.6,],digits=7), format="markdown",row.names=FALSE))
            }
            cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Spearman correlation.\n\n")
            #Creation of the excel file with all the results
            filename <-'results.xlsx'
            results_xls <- results_xls[idx,]
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
            cat("\nThe table below shows the best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[abs(results[,1])>=0.5,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below shows the best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[abs(results[,1])>=0.6,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Pearson correlation.\n\n")
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
            cat("\nThe table below shows the best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[abs(results[,2])>=0.5,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below shows the best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[abs(results[,2])>=0.6,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Spearman correlation.\n\n")
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
  } else {  # more than two groups
    if ((u_lab=='c')|(u_lab=='d')){
      pval1 <- aperm(mw_ncPCA, c(1,3,2))
      pval1 <- matrix(pval1,length(norm)*dim(mw_ncPCA)[3],choose(numbLabels,2))
      storage.mode(pval1)<- "numeric"
      pval2 <- aperm(mw_cPCA,c(1,3,2))
      pval2 <- matrix(pval2,length(norm)*dim(mw_cPCA)[3],choose(numbLabels,2))
      storage.mode(pval2)<- "numeric"
      pvals <- rbind(pval1,pval2)
      avg <- apply(pvals, 1, mean)
      
      AUC1 <- aperm(AUC_nc,c(1,3,2))
      AUC1 <- matrix(AUC1,length(norm)*dim(mw_ncPCA)[3],choose(numbLabels,2))
      storage.mode(AUC1)<- "numeric"
      AUC2 <- aperm(AUC_c,c(1,3,2))
      AUC2 <- matrix(AUC2,length(norm)*dim(mw_cPCA)[3],choose(numbLabels,2))
      storage.mode(AUC2)<- "numeric"
      AUCs <- rbind(AUC1,AUC2)
      avg_AUC <- apply(AUCs, 1, mean)
      
      AUPR1 <- aperm(AUPR_nc,c(1,3,2))
      AUPR1 <- matrix(AUPR1,length(norm)*dim(mw_ncPCA)[3],choose(numbLabels,2))
      storage.mode(AUPR1)<- "numeric"
      AUPR2 <- aperm(AUPR_c,c(1,3,2))
      AUPR2 <- matrix(AUPR2,length(norm)*dim(mw_cPCA)[3],choose(numbLabels,2))
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
        header <- c('Avg P-val','Avg AUC','Avg AUPR','Norm','Centering','Dimension','expl Var')
        results <- data.frame(avg,avg_AUC,avg_AUPR,norms,centr,dim,variance, row.names = NULL)
        colnames(results) <- header
        
      } else if (u_lab =='d'){
        header <- c('Avg Pval','Avg AUC','Avg AUPR','pears','spear','Norm','Centering','Dim','expl Var')
        results <- data.frame(avg,avg_AUC,avg_AUPR,pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
        colnames(results) <- header
        
      } 
      
      if (u_lab =='c'){
        header_xls1 <- c('Avg P-val',group_head,'Avg AUC',group_head,'Avg AUPR',group_head,'Norm','Centering','Dimension','explained Variance')
        header_xls <- c('Avg P-val',group_head,'Avg AUC',group_head,'Avg AUPR',group_head,'Norm','Centering','Dimension','explained Variance','Trustworthiness(p-value)','Trustworthiness(AUC)','Trustworthiness(AUPR)')
        results_xls<- data.frame(avg,pvals,avg_AUC,AUCs,avg_AUPR, AUPRs, norms,centr,dim,variance, row.names = NULL)
        colnames(results_xls) <- header_xls1
      }
      
      if (u_lab=='d'){
        header_xls1 <- c('Avg P-val',group_head,'Avg AUC',group_head,'Avg AUPR',group_head,'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance')
        header_xls <- c('Avg P-val',group_head,'Avg AUC',group_head,'Avg AUPR',group_head,'pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance','Trustworthiness(p-value)','Trustworthiness(AUC)','Trustworthiness(AUPR)')
        pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
        spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
        results_xls<- data.frame(avg,pvals,avg_AUC,AUCs,avg_AUPR, AUPRs,pears_corr, spear_corr,norms,centr,dim,variance, row.names = NULL)
        colnames(results_xls) <- header_xls1
      } 
      
    } else {
      header <- c('pears','spear','Norm','Centering','Dimension','expl Var')
      pears_corr <- c(rank_pears_corr_ncPCA,rank_pears_corr_cPCA)
      spear_corr <- c(rank_spear_corr_ncPCA,rank_spear_corr_cPCA)
      header_xls <- c('pearson-correlation','spearman-correlation','Norm','Centering','Dimension','explained Variance')
      results_xls <- data.frame(pears_corr,spear_corr,norms,centr,dim,variance, row.names = NULL)
      results <- results_xls
      colnames(results_xls) <- header_xls
      colnames(results) <- header
    }
    
    
    flag <- 0
    while (flag == 0){
      if ((u_lab =='c')|(u_lab =='d')){
        if (u_lab =='c'){ 
          cat("\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR? [p/auc/aupr]:\n\n")
          u_rank <- readline(prompt="-> ")
        } else {
          cat("\nWould you like to rank the PCA results by \n[p]    P-value\n[auc]  AUC\n[aupr] AUPR\n[pc]   Pearson correlation\n[sc]   Spearman correlation? [p/auc/aupr/pc/sc]:\n\n")
          u_rank <- readline(prompt="-> ")
        }
        
        flag <- 1
        
        if (u_rank =='p'){
          
          idx <- order(results[,1],decreasing=FALSE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[results[,1]<0.05,])==0){
            cat("\nThe table below show best results (p-value<0.25), automatically ranked with respect to P-value  from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.");
            cat("\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n")
            print(kable(format(results[results[,1]<0.25,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below show best results (p-value<0.05), automatically ranked with respect to P-value  from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.");
            cat("\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n")
            print(kable(format(results[results[,1]<0.05,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to P-value.\n\n")
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
            cat("\nThe table below show best results (AUC>=0.6), automatically ranked with respect to AUC from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.");
            cat("\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n")
            print(kable(format(results[results[,2]>=0.6,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below show best results (AUC>=0.7), automatically ranked with respect to AUC from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.");
            cat("\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n")
            print(kable(format(results[results[,2]>=0.7,],digits=7), format="markdown",row.names=FALSE))
          }
          cat('\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to AUC.\n\n')
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
            cat("\nThe table below show best results (AUPR>=0.6), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.");
            cat("\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n")
            print(kable(format(results[results[,3]>=0.6,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below show best results (AUPR>=0.7), automatically ranked with respect to AUPR from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.");
            cat("\nSince there are more than two groups, only the average value (and not the values in each pairwise comparison) for each evaluator is present.\n\n")
            print(kable(format(results[results[,3]>=0.7,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to AUPR.\n\n")
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
              cat("\nThe table below show best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n");
              print(kable(format(results[abs(results[,4])>=0.5,],digits=7), format="markdown",row.names=FALSE))
            } else {
              cat("\nThe table below show best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n");
              print(kable(format(results[abs(results[,4])>=0.6,],digits=7), format="markdown",row.names=FALSE))
            }
            cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Pearson correlation.\n\n")
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
              cat("\nThe table below show best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
              print(kable(format(results[abs(results[,5])>=0.5,],digits=7), format="markdown",row.names=FALSE))
            } else {
              cat("\nThe table below show best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
              print(kable(format(results[abs(results[,5])>=0.6,],digits=7), format="markdown",row.names=FALSE))
            }
            cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Spearman correlation.\n\n")
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
        cat("\nWould you like to rank the PCA results by \n[pc] Pearson correlation \n[sc] Spearman correlation? [pc/sc]\n\n")
        u_rank <- readline(prompt="-> ")
        
        flag <- 1
        
        if (u_rank == 'pc'){
          idx <- order(abs(results[,1]),decreasing=TRUE)
          results <- results[idx,]
          results_xls <- results_xls[idx,]
          #Only some results are shown on the screen
          if (nrow(results[abs(results[,1])>=0.6,])==0){
            cat("\nThe table below show best results (|pears|>=0.5), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n")
            print(kable(format(results[abs(results[,1])>=0.5,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below show best results (|pears|>=0.6), automatically ranked with respect to Pearson correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.\n\n");
            print(kable(format(results[abs(results[,1])>=0.6,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Pearson correlation.\n\n")
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
            cat("\nThe table below show best results (|spear|>=0.5), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.")
            print(kable(format(results[abs(results[,2])>=0.5,],digits=7), format="markdown",row.names=FALSE))
          } else {
            cat("\nThe table below show best results (|spear|>=0.6), automatically ranked with respect to the Spearman correlation from the most discriminative to the least discriminative, in order to assist the user in the creation of the network from the most discriminative PCn.")
            print(kable(format(results[abs(results[,2])>=0.6,],digits=7), format="markdown",row.names=FALSE))
          }
          cat("\nAll the results are returned in your current folder in an Excel table, named results.xlsx in the spreadsheet PCA results, automatically ranked with respect to Spearman correlation.\n\n")
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
  
  cat("\nAfter seeing the ranked results, you can choose any combination of normalization, centring,dimension (PCn) and cut-off for the network to obtain the PC-corr network  according to your interest (or need).\n")
  # User interaction --------------------------------------------------------
  if ((u_norm_opt == 1)|(u_norm_opt == 2)){
    u_norm <- 1
    u_norm_n <- norms_list[1]
  } else {
    
    flag <- 0
    while (flag == 0){
      cat("\n\nSelect the normalization (for detailed information see the User guide):\n\n")
      poss_norm_case <- sample(norms_list,replace=FALSE)
      while (poss_norm_case[1] == '-'){
        poss_norm_case <- sample(norms_list,replace=FALSE)
      }
      cat(paste("Examples: ",poss_norm_case[1],"or - \n",sep=" "))
      cat("where - stands for no normalization.\n\n")
      u_norm_n <- readline(prompt="-> ")
      flag <- 1
      u_norm <- which(norms_list %in% u_norm_n)
      if (length(u_norm) == 0){
        flag <- 0
        cat("Please introduce the exact name of the normalization.\n")
      }
    }
  }
  
  flag <- 0
  while (flag == 0){
    cat("\nCentering version?[y/n]\n[y] yes, centred PCA\n[n] no, non-centred PCA\n\n")
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
    if ((u_dim > 0) & (u_dim <= dim(ncPCA[[1]])[2])){
      flag <- 1
    } else{
      cat("Please introduce an existing dimension.\n")
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
          if (ind1 == dim(ncPCA[[1]])[2]){
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
          if (ind1 == dim(ncPCA[[1]])[2]){
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
        if (ind1 == dim(ncPCA[[1]])[2]){
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
          if (ind1 ==dim(ncPCA[[1]])[2]){
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
          if (ind1 == dim(ncPCA[[1]])[2]){
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
      if (ind1 == dim(ncPCA[[1]])[2]){
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
      runif_numb <- runif(3,0.15,0.85)
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
      bar_data <- data.frame(dimens=factor(1:dim(ncPCA[[1]])[2]),evaluator=mw_PCA[u_norm,1,])
      if (u_rank == 'p'){
        indexmin <- which(bar_data$evaluator %in% min(bar_data$evaluator))
        color_bars <- rep('steelblue', dim(bar_data)[1]) 
        color_bars[indexmin] <- "skyblue"
        g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=dimens)) + geom_bar(stat='identity') + scale_fill_manual(values=color_bars) +theme(legend.position='none')+ geom_hline(aes(yintercept=0.05),linetype="dashed", color = "red", size=1)+ scale_y_continuous(breaks=c(0.0,0.05,0.2,0.4,0.6,0.8,1.0),limits=c(0,1)) +   theme(axis.text.y = element_text(colour = c('black', "#FF0000",'black', 'black', 'black', 'black','black')))
        save(g1,indexmin,color_bars,bar_data,file='plot.RData')
        for (i in 1:length(indexmin)){
          g1 <- g1 + annotate(geom="text", x= indexmin[i], y=min(bar_data$evaluator)+0.02,label = round(min(bar_data$evaluator),digits=2), color = "#FF0000",size=3.5)
        }
      } else{
        indexmax <- which(bar_data$evaluator %in% max(bar_data$evaluator))
        color_bars <- rep('steelblue',dim(bar_data)[1]) 
        color_bars[indexmax] <- "skyblue"
        g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=dimens)) + geom_bar(stat='identity') + scale_fill_manual(values=color_bars) +theme(legend.position='none') + scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),limits=c(0,1)) + geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)
        for (i in 1:length(indexmax)){
          g1 <- g1 + annotate(geom="text", x= indexmax[i], y=max(bar_data$evaluator)+0.03,label = round(max(bar_data$evaluator),digits=2), color = "#FF0000",size=3.5)
        }

        
      }
    } else {
      a1 <- which(mw_PCA[u_norm,1,] >= 0)
      a2 <- which(mw_PCA[u_norm,1,] < 0)
      lev <- c()
      indexmax <- which(abs(mw_PCA[u_norm,1,]) == max(abs(mw_PCA[u_norm,1,])))
      lev[a1] <- "pos"
      lev[a2] <- "neg"
      bar_data <- data.frame(dimens=factor(1:dim(ncPCA[[1]])[2]),evaluator=abs(mw_PCA[u_norm,1,]),lev)
      g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=lev)) + geom_bar(stat='identity') 
      if (u_rank == 'pc'){
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Pear. corr < 0","Pear. corr \u2265 0"),values=c("#808080", "#000000"))+ geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)+scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),limits=c(0,1))  
      } else {
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Spear. corr < 0","Spear. corr \u2265 0"),values=c("#808080", "#000000"))+geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)+scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),limits=c(0,1))  
      }
      for (i in 1:length(indexmax)){
        g1 <- g1 + annotate(geom="text", x= indexmax[i], y=max(abs(mw_PCA[u_norm,1,]))+0.03,label = round(max(abs(mw_PCA[u_norm,1,])),digits=2), color = "#FF0000",size=3.5)
      }

    }
  } else {
    if ((u_lab == 'c')||(u_lab == 'd')){
      if ((u_rank !='pc') && (u_rank!='sc')){
        bar_data <- data.frame(dimens=factor(1:dim(ncPCA[[1]])[2]),evaluator=colMeans(mw_PCA[u_norm,,]))
        if (u_rank == 'p'){
          indexmin <- which(bar_data$evaluator %in% min(bar_data$evaluator))
          color_bars <- rep('steelblue',dim(bar_data)[1]) 
          color_bars[indexmin] <- "skyblue"
          g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=dimens)) + geom_bar(stat='identity') + scale_fill_manual(values=color_bars) + theme(legend.position='none')+ geom_hline(aes(yintercept=0.05),linetype="dashed", color = "red", size=1)+  scale_y_continuous(breaks=c(0.0,0.05,0.2,0.4,0.6,0.8,1.0),limits=c(0,1)) +   theme(axis.text.y = element_text(colour = c('black', "#FF0000",'black', 'black', 'black', 'black','black')))
          for (i in 1:length(indexmin)){
            g1 <- g1 + annotate(geom="text", x= indexmin[i], y=min(bar_data$evaluator)+0.02,label = round(min(bar_data$evaluator),digits=2), color = "#FF0000",size=3.5)
          }

        } else{
          indexmax <- which(bar_data$evaluator %in% max(bar_data$evaluator))
          color_bars <- rep('steelblue',dim(bar_data)[1]) 
          color_bars[indexmax] <- 'skyblue'
          g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=dimens)) + geom_bar(stat='identity') + scale_fill_manual(values=color_bars) +theme(legend.position='none') +scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),limits=c(0,1))  +geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)
          for (i in 1:length(indexmax)){
            g1 <- g1 + annotate(geom="text", x= indexmax[i], y=max(bar_data$evaluator)+0.03,label = round(max(bar_data$evaluator),digits=2), color = "#FF0000",size=3.5)
          }
          
        }
      } else {
        a1 <- which(mw_PCA[u_norm,1,] >= 0)
        a2 <- which(mw_PCA[u_norm,1,] < 0)
        indexmax <- which(abs(mw_PCA[u_norm,1,]) %in% max(abs(mw_PCA[u_norm,1,])))
        lev <- c()
        lev[a1] <- "pos"
        lev[a2] <- "neg"
        bar_data <- data.frame(dimens=factor(1:dim(ncPCA[[1]])[2]),evaluator=abs(mw_PCA[u_norm,1,]),lev)
        g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=lev)) + geom_bar(stat='identity') 
        if (u_rank == 'pc'){
          g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Pear. corr < 0","Pear. corr \u2265 0"),values=c("#808080", "#000000"))+ geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1) +scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),limits=c(0,1))  
        } else {
          g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Spear. corr < 0","Spear. corr \u2265 0"),values=c("#808080", "#000000"))+ geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)+scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8,1.0),limits=c(0,1))  
        }
        for (i in 1:length(indexmax)){
          g1 <- g1 + annotate(geom="text", x= indexmax[i], y=max(abs(mw_PCA[u_norm,1,]))+0.03,label = round(max(abs(mw_PCA[u_norm,1,])),digits=2), color = "#FF0000",size=3.5)
        }

      }
    } else {
      a1 <- which(mw_PCA[u_norm,1,] >= 0)
      a2 <- which(mw_PCA[u_norm,1,] < 0)
      lev <- c()
      indexmax <- which(abs(mw_PCA[u_norm,1,]) %in% max(abs(mw_PCA[u_norm,1,])))
      lev[a1] <- "pos"
      lev[a2] <- "neg"
      bar_data <- data.frame(dimens=factor(1:dim(ncPCA[[1]])[2]),evaluator=abs(mw_PCA[u_norm,1,]),lev)
      g1 <- ggplot(data=bar_data, aes(x=dimens,y=evaluator,fill=lev)) + geom_bar(stat='identity') 
      if (u_rank == 'pc'){
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Pear. corr < 0","Pear. corr \u2265 0"),values=c("#808080", "#000000")) + geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)
      } else {
        g1 <- g1 +  theme(legend.title = element_blank())+ scale_fill_manual(labels = c("Spear. corr < 0","Spear. corr \u2265 0"),values=c("#808080", "#000000"))+ geom_hline(aes(yintercept=0.5),linetype="dashed", color = "red", size=1)
      }
      
      for (i in 1:length(indexmax)){
        g1 <- g1 + annotate(geom="text", x= indexmax[i], y=max(abs(mw_PCA[u_norm,1,]))+0.03,label = round(max(abs(mw_PCA[u_norm,1,])),digits=2), color = "#FF0000",size=3.5)
      }

    }
  }
  eval_title <- c(evaluat,'s over principal components in ',ttl,' for norm: ',u_norm_n)
  if (( u_rank!='pc') && (u_rank!='sc')){
    g1 <- g1 + xlab("PC") + ylab(paste(str_ylab,collapse="")) + labs(title=paste(eval_title, collapse = "")) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, vjust= 0.5, hjust = 1))
  } else {
    g1 <- g1 + xlab("PC") + ylab(paste(str_ylab,collapse="")) + labs(title=paste(eval_title, collapse = "")) + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, vjust= 0.5, hjust = 1),legend.position = c(0.92,0.9))
  }
  if (( u_rank!='pc') && (u_rank!='sc')){
    g1 <- g1 + ylab(evaluat)
  } else {
    g1 <- g1 + ylab(paste(c('|',evaluat,'|'), collapse=""))
  }
  
  
  explain_var <- data.frame(dimens=factor(1:dim(ncPCA[[1]])[2]),explain=explained[[u_norm]])
  if (( u_rank!='pc') && (u_rank!='sc')){
    g2 <- ggplot(data=explain_var, aes(x=dimens, y= explain,fill=dimens))+ geom_bar(stat='identity') + scale_fill_manual(values=color_bars) +theme(legend.position='none')+labs(title= 'Explained variance for the respective principal components', x= 'PC', y =  'Explained Variance (%)') + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, vjust= 0.5, hjust = 1))
  } else {
    g2 <- ggplot(data=explain_var, aes(x=dimens, y= explain))+ geom_bar(stat='identity', fill='steelblue') + labs(title= 'Explained variance for the respective principal components', x= 'PC', y =  'Explained Variance (%)') + theme(plot.title = element_text(hjust = 0.5), axis.text.x=element_text(angle=90, vjust= 0.5, hjust = 1))
  }
  
  dev.new()
  g1 <- g1 + theme(plot.title = element_text(size = 14))  
  g2 <- g2 + theme(plot.title = element_text(size = 14))  
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
    feat_names_temp <- feat_names
    
    index <- abs(V) > cutoff
    n <- length(V)-sum(index)

    
    
    if (n>0){
      ns <- toString(n) 
      if (n == 1){
        cat(paste("\n",ns,"feature was deleted because |V_feature|<cutoff.\n",sep = ""))
      } else if (n > 1){
        cat(paste("\n",ns," features were deleted because |V_feature|<cutoff",sep = ""))
      }
      x <- x[,index]
      V <- V[index]
      feat_names <- feat_names[index]
    }   
    

    if (length(V) > 1){
      # - Normal Pearson Corrleation calculation
      c <- cor(x)
      
      # - Apply the C-Corr formula
      pc_corr <- matrix(0,length(V),length(V))
      
      for (i in 1:(length(V)-1)) {
        for (j in (i+1):length(V)) {
          pc_corr[i,j] <- sign(c[i,j]) * min(abs(c(c[i,j],V[i],V[j])))
        }
      } 
    } else {
      pc_corr <- matrix(0,length(V),length(V))
    }

    
    if((max(abs(pc_corr[,])) <= cutoff)|(length(V)<=1)) {
      c_temp <- cor(x_temp)
      pc_corr_temp <- matrix(0,length(V_temp),length(V_temp))
      for (i in 1:(length(V_temp)-1)){
        for (j in (i+1):length(V_temp)){
          pc_corr_temp[i,j] <- sign(c_temp[i,j]) * min(abs(c(c_temp[i,j],V_temp[i],V_temp[j])))
        }
      }
      cat('\nWith this cut-off, there are no edges that have |PC_corr(i,j)| > cutoff.\n\n')
    
      max_cutoff <- max(abs(pc_corr_temp[,]))
      max_cutoff_d2 <- max_cutoff - round(max_cutoff, digits = 2)
      if ((sign(max_cutoff_d2) == 1) | (sign(max_cutoff_d2) == 0)){
        max_possibl_cutoff <- round(max_cutoff, digits = 2)
      } else {
        max_possibl_cutoff <- round(max_cutoff, digits = 2) - 0.01
      }
                        
        
      cat(paste('\nTry with another cut-off less than ',toString(max_possibl_cutoff),'.\n\n', sep = ""))
      
      flag <- 0 
      while (flag == 0){
        cat(paste("\nThen select another cut-off for generating the PC-corr network [number between 0 and ",toString(max_possibl_cutoff),"]:\n\n", sep = ""))
        cutoff <- readline(prompt="-> ")
        cutoff <- eval(parse(text=cutoff))
        if ((cutoff >= 0) & (cutoff <= max_possibl_cutoff)){
          flag <- 1
        } else {
          cat(paste("Please introduce a correct cut-off (in [0,",toString(max_possibl_cutoff),"]", sep =""))
        }
      }
      
      if (flag == 1){
        index <- abs(V_temp) > cutoff
        n <- length(V_temp)-sum(index)
        if (n>0){
          if (n == 1){
            cat(paste("\n",ns,"feature was deleted because |V_feature|<cutoff.\n",sep = ""))
          } else if (n > 1){
            cat(paste("\n",ns," features were deleted because |V_feature|<cutoff",sep = ""))
          }
          x <- x_temp[,index]
          V <- V_temp[index]
          feat_names <- feat_names_temp[index]
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
      }
    }
    
    
    idx <- (pc_corr < (-cutoff)) | (pc_corr > cutoff)
    pc_corr[idx==0] <- 0
    
    # Remove the singletons
    pc_corr <- pc_corr + t(pc_corr) # symmetric
    ma <- (colSums(pc_corr[,]==0)) == dim(pc_corr)[1]
    pc_corr<- pc_corr[,!ma]
    pc_corr<- pc_corr[!ma,]
    sig_PC_corr <- pc_corr
    sig_PC_corr_Name <- feat_names[!ma]
    
    #Datamatrix with only the features in the network
    x1 <- x
    x1 <- x1[,!ma]
    
    cutoff_f <- cutoff
    
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
    
    return(list(Edges=edges,Nodes=nodes,pc_corr=pc_corr, x1=x1, cutoff_f=cutoff_f))
    
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
  
  # User interaction for cut-off -------------------------------------------------------------------------------------------
  flag <- 0 
  while (flag == 0){
    cat("\nSelect a cut-off or a set of cut-offs for generating the PC-corr network [number between 0 and 1]:\n\n")
    cat("Examples: 0.6 or c(0.6, 0.65, 0.7)\n\n")
    cutoff <- readline(prompt="-> ")
    cutoff <- eval(parse(text=cutoff))
    if ((max(cutoff) >= 0) & (max(cutoff) <= 1)){
      flag <- 1
    } else{
      cat("Please introduce a correct cut-off or a set of cut-offs (in [0,1]).\n")
    }
  }
  
  cut_off <-c()
  if (length(cutoff) == 1){
    pc_corr_res <- C_corr(norm[[u_norm]],pc[[u_norm]][,u_dim],feat_names,cutoff)
    cutoff_f <- pc_corr_res$cutoff_f
    pc_corr <- pc_corr_res$pc_corr
    x1 <- pc_corr_res$x1
    Edges <- pc_corr_res$Edges
    samplCol <- match_V_samplCol(col,pc_corr_res$x1,labels,nameLabels,pc_corr_res$Nodes)
    n1_f <- samplCol$n1_f
    n2_f <- samplCol$n2_f
    Nodes <- data.frame(pc_corr_res$Nodes[,1], samplCol$NodeColor,pc_corr_res$Nodes[,2])
    colnames(Nodes) <- c("Node","Colour","Loading (V)")
    
    cut_off <- cutoff_f
    edges <- Edges
    nodes <- Nodes
    # print(cut_off)
    # print(edges)
    # print(nodes)
    if ((dim(edges)[1]==1)) {
      cat(paste('\nAt cut-off',toString(round(cut_off,digits=2)),'the PC-corr network has', toString(dim(nodes)[1]),'nodes and',toString(dim(edges)[1]),'edge.\n\n',sep=" "))
    } else {
      cat(paste('\nAt cut-off',toString(round(cut_off,digits=2)),'the PC-corr network has', toString(dim(nodes)[1]),'nodes and',toString(dim(edges)[1]),'edges.\n\n',sep=" "))
    }
    
  } else {
    Edges <- matrix(list(),length(cutoff),2)
    Nodes <- matrix(list(),length(cutoff),2)
    pc_corr <- matrix(list(),length(cutoff),2)
    x2 <- matrix(list(),length(cutoff),2)
    cutoff_f <- c()
    for (i in 1:length(cutoff)){
      Edges[[i,1]] <- cutoff[i]
      Nodes[[i,1]] <- cutoff[i]
      pc_corr[[i,1]] <- cutoff[i]
      x2[[i,1]] <- cutoff[i]
      
      pc_corr_res <- C_corr(norm[[u_norm]],pc[[u_norm]][,u_dim],feat_names,cutoff[i])
      Edges[[i,2]] <- pc_corr_res$Edges
      pc_corr[[i,2]] <- pc_corr_res$pc_corr
      x2[[i,2]] <- pc_corr_res$x1
      cutoff_f[i] <- pc_corr_res$cutoff_f
        
      samplCol <- match_V_samplCol(col,pc_corr_res$x1,labels,nameLabels,pc_corr_res$Nodes)
      
      Nodes[[i,2]] <- data.frame(pc_corr_res$Nodes[,1], samplCol$NodeColor,pc_corr_res$Nodes[,2])
      colnames(Nodes[[i,2]]) <- c("Node","Colour","Loading (V)")
      
      
      cut_off[i] <-  cutoff_f[i]
      edges <- Edges[[i,2]]
      nodes <- Nodes[[i,2]]
      # print(cut_off[i])
      # print(edges)
      # print(nodes)
      if ((dim(edges)[1]==1)) {
        cat(paste('\nAt cut-off,',toString(round(cut_off[i],digits=2)),'the PC-corr network has', toString(dim(nodes)[1]),'nodes and',toString(dim(edges)[1]),'edge.\n\n',sep=" "))
      } else {
        cat(paste('\nAt cut-off,',toString(round(cut_off[i],digits=2)),'the PC-corr network has', toString(dim(nodes)[1]),'nodes and',toString(dim(edges)[1]),'edges.\n\n',sep=" "))
      }
  
    }
  }

  filename <- 'PC-corr_net_edges-nodes_results.xlsx'
  if (file.exists(filename)){
    a_xlsx <- file.remove(filename)
  } 
  
  
  if (length(cut_off) == 1){
    write.xlsx(Edges, file=filename, sheetName = paste(c('Edges-cutoff ',toString(cut_off)),collapse=""), row.names = FALSE,append =TRUE)
    write.xlsx(Nodes, file=filename, sheetName = paste(c('Nodes-cutoff ',toString(cut_off)),collapse=""), row.names = FALSE,append=TRUE)
  } else {
    for (i in 1:length(cut_off)){
      write.xlsx(Edges[[i,2]], file=filename, sheetName = paste(c('Edges-cutoff ',toString(cut_off[i])),collapse=""), row.names = FALSE, append=TRUE)
      write.xlsx(Nodes[[i,2]], file=filename, sheetName = paste(c('Nodes-cutoff ',toString(cut_off[i])),collapse=""), row.names = FALSE, append=TRUE)
      
    }
  }  
  
  cat("\nYou can find the tables of edge and node values of the PC-corr network in your current folder in an Excel file (in separate spreadsheet), named PC-corr_net_edges-nodes_results.xlsx, so that you can easily visualize the graph with another network visualization program.\n\n")
  
  
  # Plot PC-corr network -------------------------------------------------------------------------------------------------------------------------------------
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
    frust_txt1 <- bquote(cut-off == .(cutoff))
    frust_txt2 <- bquote(paste(frustration == frac(.(edg_frust),.(nb_edges)),phantom()==.(fr_edg_frust_rd),"%"))
    text( usr[ 1 ]-0.17, usr[ 4 ]+0.2,frust_txt1, adj = c( 0, 1 ), font = 14)
    text( usr[ 1 ]-0.17, usr[ 4 ]+0.1, frust_txt2, adj = c( 0, 1 ), font = 14)

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
  
  
  # User interaction for the visualization of the PC-corr network -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  flag <- 0
  
  while (flag == 0){
    if (length(cut_off)==1){
      cat("\nDo you want to visualize the PC-corr network? [y/n]:\n\n")
    } else {
      cat("\nDo you want to visualize the PC-corr networks? [y/n]:\n\n")
    }
    u_vis_net <- readline(prompt="-> ")
    if ((u_vis_net == 'n') | (u_vis_net =='y')){
      flag <- 1
    } else {
      cat("Please introduce just 'y' or 'n'\n")
    }
  }
  
  #Graph plot ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (u_vis_net == 'y'){
    # Visualize the PC-corr networks
    for (k in 1:length(cut_off)){
      if (length(cut_off)==1){
        pc_corr1 <- pc_corr
        nodes1 <- Nodes
      } else {
        pc_corr1 <- pc_corr[[k,2]]
        nodes1 <- Nodes[[k,2]]
      }
      if ((dim(pc_corr1)[1] == 0) & (dim(pc_corr1)[2] == 0)){
        next
      } else{
        
        na_rows <- which(rowSums(is.na(nodes1)) > 0)
        indices_to_remove <- nodes1$index[na_rows]
        nodes1 <- nodes1[-na_rows, ]
        pc_corr1 <- pc_corr1[-indices_to_remove, -indices_to_remove]
        
        plot_graph(pc_corr1,nodes1,cut_off[k])
      }
    }
  }
  
  
  
  # User interaction for the significance of measure of segregation -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  if ((u_lab == 'c') | (u_lab =='d')){
    flag <- 0
    cat("\nAdditionally, you can compute the trustworthiness of p-value, AUC and AUPR results, that is a measure of significance of segregation.\n")
    cat("The trustworthiness (of p-value/AUC/AUPR results) evaluates if the results are discriminative because it is random (trustworthiness p-value>0.05) or it captures main variability (trustworthiness p-value<=0.05).\n")
    while (flag == 0){
      cat("\nWould you like to compute the trustworthiness of measure of segregation (p-value,AUC and AUPR)? [1/2/3]\n[1] no \n[2] yes, but just for the selected case (normalization, dimension and centering)\n[3] yes, for all the cases\n\n")
      meas_seg <- readline(prompt="-> ")
      if ((meas_seg != 1) & (meas_seg != 2) & (meas_seg != 3)){
        flag <- 0
        cat("\nPlease introduce either 1 for no calculation of trustworthiness, 2 for its calculation for the the selected case (normalization, dimension and centering) or 3 for its calculation for all the cases. \n")
      } else {
        flag <- 1
      }
    }
  }

  
  if (meas_seg == 3){
    flag <- 0
    while (flag == 0){
      if ((dim(results)[1] *choose(numbLabels,2)) > 150){
        cat("\nAre you sure? Do you want to compute for all? It will take some time (hours). [y/n]")
        cat("\n(We suggest running it overnight)")
      } else {
        cat("\nAre you sure? Do you want to compute for all? It will take some time (minutes). [y/n]")
      }

      cat("\n[y] yes, I want to compute for all the cases\n[n] no, I don't want to compute for all the cases\n\n")
      meas_seg_all <- readline(prompt="-> ")
      if ((meas_seg_all == 'y') || (meas_seg_all == 'n')){
        flag <- 1
      } else{
        cat("\nPlease introduce just 'y' or 'n'\n")
      }
    }
  }
  
  if ((meas_seg == 2) || (meas_seg == 3)){
    labels_rp <- random_permutation_labels(1000,labels)
  }
  
  
  # Calculation of trustworthiness (pvalue,AUC,AUPR) for the selected options: norm, dimension and centred/non-centred version of PCA----------------------------------------------------------------------------------------------------
  if (meas_seg == 2){
    if (numbLabels ==2) {
      if (u_lab =='c'){
        if (u_cent == 'y'){
          u_cent_n <- 'yes'
        } else if (u_cent == 'n'){
          u_cent_n <- 'no'
        }
        chosen_opt <- (results_xls[,4] == u_norm_n) & (results_xls[,5] == u_cent_n) & (results_xls[,6] == u_dim)
        idx_opt <- which(chosen_opt)
        mw_opt <- results_xls[idx_opt,1]
        auc_opt <- results_xls[idx_opt,2]
        aupr_opt <- results_xls[idx_opt,3]
        numb_rand <- 1000
        # Significance of segragation (pvalue,AUC,AUPR)
        pVal_opt <- significance_segragation(numb_rand,PCA[[u_norm]],u_dim,labels_rp,nameLabels,numbLabels,c(mw_opt,auc_opt,aupr_opt),'all',u_aupr,u_aupr_r)
        pVal_mw_opt <- pVal_opt[1]
        pVal_auc_opt <- pVal_opt[2]
        pVal_aupr_opt <- pVal_opt[3]
        idx_excel <- which(header_xls %in% 'Trustworthiness(p-value)')
        p_Val_sign_p <- character(length=dim(results)[1])
        p_Val_sign_auc <-  p_Val_sign_p 
        p_Val_sign_aupr <-  p_Val_sign_p 
        results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names = NULL),row.names = NULL)
        results_xls[,idx_excel:(idx_excel +2)]<-lapply(results_xls[,idx_excel:(idx_excel +2)],as.character)
        results_xls[idx_opt,idx_excel:(idx_excel +2)] <- as.character(c(pVal_mw_opt,pVal_auc_opt,pVal_aupr_opt))
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)        
        colnames(results_xls) <- header_xls
        
      } else if (u_lab=='d'){
        if (u_cent == 'y'){
          u_cent_n <- 'yes'
        } else if (u_cent == 'n'){
          u_cent_n <- 'no'
        }
        chosen_opt <- (results_xls[,6] == u_norm_n) & (results_xls[,7] == u_cent_n) & (results_xls[,8] == u_dim)
        idx_opt <- which(chosen_opt)
        mw_opt <- results_xls[idx_opt,1]
        auc_opt <- results_xls[idx_opt,2]
        aupr_opt <- results_xls[idx_opt,3]
        numb_rand <- 1000
        # Significance of segragation (pvalue,AUC,AUPR)
        pVal_opt <- significance_segragation(numb_rand,PCA[[u_norm]],u_dim,labels_rp,nameLabels,numbLabels,c(mw_opt,auc_opt,aupr_opt),'all',u_aupr,u_aupr_r)
        pVal_mw_opt <- pVal_opt[1]
        pVal_auc_opt <- pVal_opt[2]
        pVal_aupr_opt <- pVal_opt[3]
        idx_excel <- which(header_xls %in% 'Trustworthiness(p-value)')
        p_Val_sign_p <- character(length=dim(results)[1])
        p_Val_sign_auc <-  p_Val_sign_p 
        p_Val_sign_aupr <-  p_Val_sign_p
        results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr, row.names = NULL),row.names = NULL)
        results_xls[,idx_excel:(idx_excel +2)]<-lapply(results_xls[,idx_excel:(idx_excel +2)],as.character)
        results_xls[idx_opt,idx_excel:(idx_excel +2)] <- as.character(c(pVal_mw_opt,pVal_auc_opt,pVal_aupr_opt))
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)        
        colnames(results_xls) <- header_xls
        
        
      } 
    } else { #more than two groups
      if (u_lab =='c'){
        if (u_cent == 'y'){
          u_cent_n <- 'yes'
        } else if (u_cent == 'n'){
          u_cent_n <- 'no'
        }
        chosen_opt <- (results[,4] == u_norm_n) & (results[,5] == u_cent_n) & (results[,6] == u_dim)
        idx_opt <- which(chosen_opt)
        mw_opt <- results[idx_opt,1]
        auc_opt <- results[idx_opt,2]
        aupr_opt <- results[idx_opt,3]
        numb_rand <- 1000
        # Significance of segragation (pvalue,AUC,AUPR)
        pVal_opt <- significance_segragation(numb_rand,PCA[[u_norm]],u_dim,labels_rp,nameLabels,numbLabels,c(mw_opt,auc_opt,aupr_opt),'all',u_aupr,u_aupr_r)
        pVal_mw_opt <- pVal_opt[1]
        pVal_auc_opt <- pVal_opt[2]
        pVal_aupr_opt <- pVal_opt[3]
        idx_excel <- which(header_xls %in% 'Trustworthiness(p-value)')
        p_Val_sign_p <- character(length=dim(results)[1])
        p_Val_sign_auc <-  p_Val_sign_p 
        p_Val_sign_aupr <-  p_Val_sign_p 
        results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names=NULL),row.names = NULL)
        results_xls[,idx_excel:(idx_excel +2)]<-lapply(results_xls[,idx_excel:(idx_excel +2)],as.character)
        results_xls[idx_opt,idx_excel:(idx_excel +2)] <- as.character(c(pVal_mw_opt,pVal_auc_opt,pVal_aupr_opt))
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)        
        colnames(results_xls) <- header_xls
        
      } else if (u_lab=='d'){
        if (u_cent == 'y'){
          u_cent_n <- 'yes'
        } else if (u_cent == 'n'){
          u_cent_n <- 'no'
        }
        chosen_opt <- (results[,6] == u_norm_n) & (results[,7] == u_cent_n) & (results[,8] == u_dim)
        idx_opt <- which(chosen_opt)
        mw_opt <- results[idx_opt,1]
        auc_opt <- results[idx_opt,2]
        aupr_opt <- results[idx_opt,3]
        numb_rand <- 1000
        # Significance of segragation (pvalue,AUC,AUPR)
        pVal_opt <- significance_segragation(numb_rand,PCA[[u_norm]],u_dim,labels_rp,nameLabels,numbLabels,c(mw_opt,auc_opt,aupr_opt),'all',u_aupr,u_aupr_r)
        pVal_mw_opt <- pVal_opt[1]
        pVal_auc_opt <- pVal_opt[2]
        pVal_aupr_opt <- pVal_opt[3]
        idx_excel_all <- which(header_xls %in% 'Trustworthiness(p-value)')
        p_Val_sign_p <- character(length=dim(results)[1])
        p_Val_sign_auc <-  p_Val_sign_p 
        p_Val_sign_aupr <-  p_Val_sign_p 
        results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names=NULL),row.names = NULL)
        results_xls[,idx_excel_all:(idx_excel_all +2)]<-lapply(results_xls[,idx_excel_all:(idx_excel_all +2)],as.character)
        results_xls[idx_opt,idx_excel_all:(idx_excel_all +2)] <- as.character(c(pVal_mw_opt,pVal_auc_opt,pVal_aupr_opt))
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)        
        colnames(results_xls) <- header_xls
        
      }
      
    } 
    
  }

  
  
# Calculation of trustworthiness (p-value,AUC and AUPR) for each case ------------------------------------------------------------------------------------------
  if (meas_seg == 3){
    tic("Total")
    if (meas_seg_all =='y'){
      
      if ( numbLabels == 2){ # two groups
        if (u_lab == 'c'){
          idx_all <- 1:dim(results_xls)[1]
          A <- round(idx_all/length(idx_all),digits=2)*100
          #A <- round(idx_all/length(idx_all),digits=1)*100
          unik <- !duplicated(A)
          ia2 <- seq_along(A)[unik]
          ia1 <- A[unik]
          #ia2 <- ia2[(ia1 %% 10) == 0]
          pVal_all <- matrix(data=NA,nrow=length(idx_all),ncol=3)
          pVal_mw_all <- c()
          pVal_auc_all <- c()
          pVal_aupr_all <- c()
          p_Val_sign_p <- character(length=dim(results)[1])
          p_Val_sign_auc <-  p_Val_sign_p 
          p_Val_sign_aupr <-  p_Val_sign_p 
          results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names = NULL),row.names = NULL)
          idx_excel_all <- which(header_xls %in% 'Trustworthiness(p-value)')
          results_xls[,idx_excel_all:(idx_excel_all +2)]<-lapply(results_xls[,idx_excel_all:(idx_excel_all +2)],as.numeric)
          cat("\n")
          for (i in 1:length(idx_all)){
            numb_rand <- 1000
            nm <- which(norms_list %in% results_xls[idx_all[i],4])
            if (results_xls[idx_all[i],5] == 'yes'){
              PCA_all <- cPCA[[nm]]
              pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results_xls[idx_all[i],6],labels_rp,nameLabels,numbLabels,c(results_xls[idx_all[i],1],results_xls[idx_all[i],2], results_xls[idx_all[i],3]),'all',u_aupr,u_aupr_r)
              pVal_mw_all[i] <-pVal_all[i,1] 
              pVal_auc_all[i] <- pVal_all[i,2] 
              pVal_aupr_all[i] <- pVal_all[i,3]  
              results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
            } else if (results_xls[idx_all[i],5] == 'no'){
              PCA_all <- ncPCA[[nm]]
              pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results_xls[idx_all[i],6],labels_rp,nameLabels,numbLabels,c(results_xls[idx_all[i],1],results_xls[idx_all[i],2], results_xls[idx_all[i],3]),'all',u_aupr,u_aupr_r)
              pVal_mw_all[i] <-pVal_all[i,1] 
              pVal_auc_all[i] <- pVal_all[i,2] 
              pVal_aupr_all[i] <- pVal_all[i,3]  
              results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
            }
            if (idx_all[i] %in% ia2){
              if (round(idx_all[i]/length(idx_all),digits=2)!=1){
                cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed")," \r")
              } else {
                cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed - Done"),"\r")
                
              }
              # cat(paste(" > Progress: ",toString(round(idx_all[i]/length(idx_all),digits=1)*100),"%\n",sep=""))
              
            }
          }
          colnames(results_xls) <- header_xls
          write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)        
          
        } else if (u_lab =='d'){
          tic("Total")
          idx_all <- 1:dim(results_xls)[1]
          A <- round(idx_all/length(idx_all),digits=2)*100
          #A <- round(idx_all/length(idx_all),digits=1)*100
          unik <- !duplicated(A)
          ia2 <- seq_along(A)[unik]
          ia1 <- A[unik]
          #ia2 <- ia2[(ia1 %% 10) == 0]
          pVal_all <- matrix(data=NA,nrow=length(idx_all),ncol=3)
          pVal_mw_all <- c()
          pVal_auc_all <- c()
          pVal_aupr_all <- c()
          p_Val_sign_p <- character(length=dim(results)[1])
          p_Val_sign_auc <-  p_Val_sign_p 
          p_Val_sign_aupr <-  p_Val_sign_p 
          results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names = NULL),row.names = NULL)
          idx_excel_all <- which(header_xls %in% 'Trustworthiness(p-value)')
          results_xls[,idx_excel_all:(idx_excel_all +2)]<-lapply(results_xls[,idx_excel_all:(idx_excel_all +2)],as.numeric)
          cat("\n")
          for (i in 1:length(idx_all)){
            numb_rand <- 1000
            nm <- which(norms_list %in% results_xls[idx_all[i],6])
            if (results_xls[idx_all[i],7] == 'yes'){
              PCA_all <- cPCA[[nm]]
              pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results_xls[idx_all[i],8],labels_rp,nameLabels,numbLabels,c(results_xls[idx_all[i],1],results_xls[idx_all[i],2], results_xls[idx_all[i],3]),'all',u_aupr,u_aupr_r)
              pVal_mw_all[i] <-pVal_all[i,1] 
              pVal_auc_all[i] <- pVal_all[i,2] 
              pVal_aupr_all[i] <- pVal_all[i,3]  
              results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
            } else if (results_xls[idx_all[i],7] == 'no'){
              PCA_all <- ncPCA[[nm]]
              pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results_xls[idx_all[i],8],labels_rp,nameLabels,numbLabels,c(results_xls[idx_all[i],1],results_xls[idx_all[i],2], results_xls[idx_all[i],3]),'all',u_aupr,u_aupr_r)
              pVal_mw_all[i] <-pVal_all[i,1] 
              pVal_auc_all[i] <- pVal_all[i,2] 
              pVal_aupr_all[i] <- pVal_all[i,3]  
              results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
            }
            if (idx_all[i] %in% ia2){
              if (idx_all != length(idx_all)){
                cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed")," \r")
              } else {
                cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed - Done")," \r")
                
              }
              # cat(paste(" > Progress: ",toString(round(idx_all[i]/length(idx_all),digits=1)*100),"%\n",sep=""))
              
            }
          }
        }
        colnames(results_xls) <- header_xls
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)   
      }
    } else { #more than two groups
      if (u_lab == 'c'){
        tic("Total")
        idx_all <- 1:dim(results)[1]
        A <- round(idx_all/length(idx_all),digits=2)*100
        #A <- round(idx_all/length(idx_all),digits=1)*100
        unik <- !duplicated(A)
        ia2 <- seq_along(A)[unik]
        ia1 <- A[unik]
        #ia2 <- ia2[(ia1 %% 10) == 0]
        pVal_all <- matrix(data=NA,nrow=length(idx_all),ncol=3)
        pVal_mw_all <- c()
        pVal_auc_all <- c()
        pVal_aupr_all <- c()
        p_Val_sign_p <- character(length=dim(results)[1])
        p_Val_sign_auc <-  p_Val_sign_p 
        p_Val_sign_aupr <-  p_Val_sign_p 
        results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names = NULL),row.names = NULL)
        idx_excel_all <- which(header_xls %in% 'Trustworthiness(p-value)')
        results_xls[,idx_excel_all:(idx_excel_all +2)]<-lapply(results_xls[,idx_excel_all:(idx_excel_all +2)],as.numeric)
        cat("\n")
        for (i in 1:length(idx_all)){
          numb_rand <- 1000
          nm <- which(norms_list %in% results[idx_all[i],4])
          if (results[idx_all[i],5] == 'yes'){
            PCA_all <- cPCA[[nm]]
            pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results[idx_all[i],6],labels_rp,nameLabels,numbLabels,c(results[idx_all[i],1],results[idx_all[i],2], results[idx_all[i],3]),'all',u_aupr,u_aupr_r)
            pVal_mw_all[i] <-pVal_all[i,1] 
            pVal_auc_all[i] <- pVal_all[i,2] 
            pVal_aupr_all[i] <- pVal_all[i,3]  
            results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
          } else if (results[idx_all[i],5] == 'no'){
            PCA_all <- ncPCA[[nm]]
            pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results[idx_all[i],6],labels_rp,nameLabels,numbLabels,c(results[idx_all[i],1],results[idx_all[i],2], results[idx_all[i],3]),'all',u_aupr,u_aupr_r)
            pVal_mw_all[i] <-pVal_all[i,1] 
            pVal_auc_all[i] <- pVal_all[i,2] 
            pVal_aupr_all[i] <- pVal_all[i,3]  
            results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
          }
          if (idx_all[i] %in% ia2){
            if (idx_all != length(idx_all)){
              cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed")," \r")
            } else {
              cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed - Done")," \r")
              
            }
            # cat(paste(" > Progress: ",toString(round(idx_all[i]/length(idx_all),digits=1)*100),"%\n",sep=""))
            
          }
        }
        colnames(results_xls) <- header_xls
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)   
      } else if (u_lab == 'd'){
        tic("Total")
        idx_all <- 1:dim(results)[1]
        #A <- round(idx_all/length(idx_all),digits=2)*100
        A <- round(idx_all/length(idx_all),digits=1)*100
        unik <- !duplicated(A)
        ia2 <- seq_along(A)[unik]
        ia1 <- A[unik]
        #ia2 <- ia2[(ia1 %% 10) == 0]
        pVal_all <- matrix(data=NA,nrow=length(idx_all),ncol=3)
        pVal_mw_all <- c()
        pVal_auc_all <- c()
        pVal_aupr_all <- c()
        p_Val_sign_p <- character(length=dim(results)[1])
        p_Val_sign_auc <-  p_Val_sign_p 
        p_Val_sign_aupr <-  p_Val_sign_p 
        results_xls <- cbind(results_xls, data.frame(p_Val_sign_p,p_Val_sign_auc,p_Val_sign_aupr,row.names = NULL),row.names = NULL)
        idx_excel_all <- which(header_xls %in% 'Trustworthiness(p-value)')
        results_xls[,idx_excel_all:(idx_excel_all +2)]<-lapply(results_xls[,idx_excel_all:(idx_excel_all +2)],as.numeric)
        cat("\n")
        for (i in 1:length(idx_all)){
          numb_rand <- 1000
          nm <- which(norms_list %in% results[idx_all[i],6])
          if (results[idx_all[i],7] == 'yes'){
            PCA_all <- cPCA[[nm]]
            pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results[idx_all[i],8],labels_rp,nameLabels,numbLabels,c(results[idx_all[i],1],results[idx_all[i],2], results[idx_all[i],3]),'all',u_aupr,u_aupr_r)
            pVal_mw_all[i] <-pVal_all[i,1] 
            pVal_auc_all[i] <- pVal_all[i,2] 
            pVal_aupr_all[i] <- pVal_all[i,3]  
            results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
          } else if (results[idx_all[i],7] == 'no'){
            PCA_all <- ncPCA[[nm]]
            pVal_all[i,] <- significance_segragation(numb_rand,PCA_all,results[idx_all[i],8],labels_rp,nameLabels,numbLabels,c(results[idx_all[i],1],results[idx_all[i],2], results[idx_all[i],3]),'all',u_aupr,u_aupr_r)
            pVal_mw_all[i] <-pVal_all[i,1] 
            pVal_auc_all[i] <- pVal_all[i,2] 
            pVal_aupr_all[i] <- pVal_all[i,3]  
            results_xls[idx_all[i],idx_excel_all:(idx_excel_all +2)] <- c(pVal_mw_all[i],pVal_auc_all[i],pVal_aupr_all[i])
          }
          if (idx_all[i] %in% ia2){
            if (idx_all != length(idx_all)){
              cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed")," \r")
            } else {
              cat(paste0("> calculating trustworthiness: ",round(idx_all[i]/length(idx_all),digits=2)*100,"% completed")," - Done\r")
              
            }
            # cat(paste(" > Progress: ",toString(round(idx_all[i]/length(idx_all),digits=1)*100),"%\n",sep=""))

          }
        }
        colnames(results_xls) <- header_xls
        write.xlsx(results_xls, file='results.xlsx', sheetName ='PCA results', row.names = FALSE)   
      }
    }
    cat("\n")
    toc()
  }
  
list2env(list(Edges=Edges,Nodes=Nodes),.GlobalEnv)

}

  
  





































