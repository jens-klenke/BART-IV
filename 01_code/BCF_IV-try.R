# bcf_iv <- function(
# Data 
data <- readRDS(data$path_in[1])

y <- data$y
w <- data$w
z <- data$z
x <- data$X

# Parameters
binary = FALSE
n_burn = 500
n_sim = 500
inference_ratio = 0.5
max_depth = 2
cp = 0.01
minsplit = 10
adj_method = "holm"
seed = 42

x_names <- paste0('x', 1:ncol(x))
    
  ######################################################
  ####         Step 0: Initialize the Data          ####
  ######################################################
  
  # Split data into Discovery and Inference
  set.seed(seed)
  index <- sample(nrow(x), nrow(x)*inference_ratio, replace=FALSE)
  
  # Initialize total dataset
  iv.data <- as.data.frame(cbind(y, w, z, x))
  names(iv.data) <- c("y", "w", "z", x_names)
  
  # Discovery and Inference Samples
  discovery <- iv.data[-index,]
  inference <- iv.data[index,]
  
  ######################################################
  ####  Step 1: Compute the Bayesian Causal Forest  ####
  ######################################################
  
  # Compute the Propensity Score though a Logistic Regression
  p.score <- glm(z ~ x[-index,],
                 family = binomial,
                 data = discovery)
  pihat <- predict(p.score, as.data.frame(x[-index,]))
  
  # Perform the Bayesian Causal Forest  to calculate the Proportion of Compliers (pic)
  bcf_pic_tree <- quiet(bartCause::bartc(w[-index], z[-index], x[-index,],
                                         n.samples = n_sim, n.burn = n_burn, 
                                         n.chains = 2L))
  
  bcf_tau_pic <- bartCause::extract(bcf_pic_tree, type = "ite")
  bcf_pic <- apply(bcf_tau_pic, 2, mean)
  # workaround! 
  bcf_pic[bcf_pic == 0] <- 1e-2 
  # mean(pic_bcf) / median(pic) == compliance
  
  ##### sparse BCF -> BCF https://github.com/albicaron/SparseBCF
  s_bcf_pic_tress <- quiet(SparseBCF::SparseBCF(w[-index], z[-index], x[-index,], pihat = pihat, 
                                      nsim = n_sim, nburn = n_burn)
                 )
  
  ######################################################
  ####     Continuous and Discrete Outcomes         ####
  ######################################################
  
#  if (binary == FALSE){
    
    # Perform the Bayesian Causal Forest for the ITT
  # is that \mu and ITT ? 
    'bcf - sparse bcf' 
    tictoc::tic()
    bcf_itt.tree <- quiet(bcf::bcf(y[-index], z[-index], x[-index,], x[-index,], pihat, nburn=n_burn, nsim=n_sim, save_tree_directory = NULL))
    bcf_tau_itt <- bcf_itt.tree$tau
    bcf_itt <- colMeans(bcf_tau_itt)
    tictoc::toc()
    
    # Get posterior of treatment effects
    bcf_tauhat <- bcf_itt/bcf_pic
    # driven by outliers
    # bcf_exp %>% dplyr::filter(x1 == 0 & x2 == 0) %>% dplyr::summarise(median(bcf_tauhat))
    bcf_exp <- as.data.frame(cbind(bcf_tauhat, x[-index,]))
    
    # repair names? !! by me 
    names(bcf_exp)[2:length(bcf_exp)] <- names(inference)[-(1:3)]
    
    
    'sparse BCF'
    ##### sparse BCF -> BCF https://github.com/albicaron/SparseBCF
    tictoc::tic()
    s_bcf_itt.tree <- quiet(
      SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                           nsim = n_sim, nburn = n_burn)
      )
    
    s_bcf_tau_itt <- s_bcf_itt.tree$tau
    s_bcf_itt <- colMeans(s_bcf_tau_itt)
    tictoc::toc()
    
    # Get posterior of treatment effects
    s_bcf_tauhat <- s_bcf_itt/bcf_pic
    # driven by outliers
    # s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
    s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))
    
    # repair names? !! by me 
    names(s_bcf_exp)[2:length(s_bcf_exp)] <- names(inference)[-(1:3)]
    
    ######################################################
    ####  Step 2: Build a CART on the Unit Level CITT ####
    ######################################################
    ' binary decision tree to discover, in an interpretable manner, the drivers of the heterogeneity??'
    bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                      data = bcf_exp,
                      maxdepth = max_depth,
                      cp=cp,
                      minsplit=minsplit)
    
    # plot tree
    rpart.plot::rpart.plot(bcf_fit.tree)
    
    # binary tree for sparse trees 
    s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                          data = s_bcf_exp,
                          maxdepth = max_depth,
                          cp=cp,
                          minsplit=minsplit)
    
    # plot tree
    rpart.plot::rpart.plot(s_bcf_fit.tree)
    
    
    ######################################################
    ####  Step 3: Extract the Causal Rules (Nodes)    ####
    ######################################################
    
    #  rules end terminal nodes?
    rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
    
    # Initialize Outputs
    bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol=7))
    names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", "Pi_obs", "ITT", "Pi_compliers")
    
    # Generate Leaves Indicator
    lvs <- leaves <- numeric(length(rules)) 
    lvs[unique(fit.tree$where)] <- 1
    leaves[rules[lvs==1]] <- 1
    
    ######################################################
    ####  Step 4: Run an IV Regression on each Node   ####
    ######################################################
    
    # Run an IV Regression on the Root
    iv.root <- ivreg(y ~ w | z,  
                     data = inference) # inference dataset
    summary <- summary(iv.root, diagnostics = TRUE)
    iv.effect.root <-  summary$coef[2,1]
    p.value.root <- summary$coef[2,4]
    p.value.weak.iv.root <- summary$diagnostics[1,4]
    proportion.root <- 1
    # share of compilers in root
    compliers.root <- length(which(inference$z==inference$w))/nrow(inference)
    itt.root <- iv.effect.root*compliers.root
    
    # Store Results for Root
    bcfivMat[1,] <- c( NA , round(iv.effect.root, 4), round(p.value.root, 4), round(p.value.weak.iv.root, 4), round(proportion.root, 4), round(itt.root, 4), round(compliers.root, 4))
    
    # Initialize New Data
    names(inference) <- paste(names(inference), sep="")
    
    # Run a loop to get the rules (sub-populations)
    for (i in rules[-1]){
      # Create a Vector to Store all the Dimensions of a Rule
      sub <- as.data.frame(matrix(NA, nrow = 1,
                                  ncol = nrow(as.data.frame(path.rpart(fit.tree, node=i, print.it = FALSE)))-1))
      quiet(capture.output(for (j in 1:ncol(sub)){
        # Store each Rule as a Sub-population
        sub[,j] <- as.character(print(as.data.frame(path.rpart(fit.tree,node=i,print.it=FALSE))[j+1,1]))
        sub_pop <- noquote(paste(sub , collapse = " & "))
      }))
      
      subset <- with(inference, inference[which(eval(parse(text=sub_pop))),])
      
      # Run the IV Regression
      if (length(unique(subset$w))!= 1 | length(unique(subset$z))!= 1){
        iv.reg <- ivreg(y ~ w | z,  
                        data = subset)
        summary <- summary(iv.reg, diagnostics = TRUE)
        iv.effect <-  summary$coef[2,1]
        p.value <- summary$coef[2,4]
        p.value.weak.iv <- summary$diagnostics[1,4]
        compliers <- length(which(subset$z==subset$w))/nrow(subset)
        itt <- iv.effect*compliers
        
        # Proportion of observations in the node
        proportion.node <- nrow(subset)/nrow(inference)
        
        ######################################################
        ####   Step 5: Output the Values of each CCACE   ####
        ######################################################
        
        bcfivMat[i,] <- c(sub_pop, round(iv.effect, 4), round(p.value, 4), round(p.value.weak.iv, 4), round(proportion.node, 4), round(itt, 4), round(compliers, 4))
      }
      
      # Delete data
      rm(subset)
    }
    
    # Adjust P.values 
    bcfiv_correction <- cbind(as.data.frame(bcfivMat), leaves)
    adj <- round(p.adjust( as.numeric(bcfiv_correction$pvalue[which(bcfiv_correction$leaves==1)]) ,  paste(adj_method)), 5)
    Adj_pvalue <- rep(NA, length(rules)) 
    Adj_pvalue[which(bcfiv_correction$leaves==1)] <- adj
    
    # Store Results
    bcfivResults <- cbind(as.data.frame(bcfivMat), Adj_pvalue)


#    } # end Continuous and Discrete Outcomes 
  
  # Return Results
  bcfivResults
