# used insite BCF_IV-own estimation
heterogeneous_treatment_estimation <- function(
    fit.tree, inference, adj_method,
    stan_model_first_stage = get("stan_model_first_stage", envir = .GlobalEnv),
    stan_model_second_stage = get("stan_model_second_stage", envir = .GlobalEnv), ...){
  # rules end terminal nodes?
  rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
  
  # Initialize Outputs
  bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol =7))
  names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", 
                       "Pi_obs", "ITT", "Pi_compliers")
  
  bayes_ivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol =7))
  names(bayes_ivMat) <- c("node", "CCACE", "CCACE l-95% CI", "CCACE u-95% CI", 
                          "Pi_obs", "ITT", "Pi_compliers")
  
  
  # Generate Leaves (end notes) Indicator
  lvs <- leaves <- numeric(length(rules)) 
  lvs[unique(fit.tree$where)] <- 1
  leaves[rules[lvs==1]] <- 1
  
  ####  Step 4: Run an IV Regression on each Node   ####
  
  # Run an IV Regression on the Root
  iv.root <- ivreg(y ~ w | z,  
                   data = inference) # inference dataset
  
  # bayes IV estimation
  bayes_iv.root <- brms_iv_function(inference, stan_model_first_stage, stan_model_second_stage)
                                 #get("stan_model_first_stage", envir = .GlobalEnv), 
                                 #get("stan_model_second_stage", envir = .GlobalEnv))
  
  # Store Results for Root
  # Store Results for Root
  bcfivMat[1, ] <- iv_summary_func(iv.root, inference, inference, sub_pop = 'root')
  bayes_ivMat[1, ] <- iv_summary_func(bayes_iv.root, inference, inference, 
                                      sub_pop = 'root', bayes = TRUE)

  # Initialize New Data
  names(inference) <- paste(names(inference), sep="")
  
  # Run a loop to get the rules (sub-populations)
  for (i in rules[-1]){
    # Create a Vector to Store all the Dimensions of a Rule
    sub <- as.data.frame(matrix(NA, nrow = 1,
                                ncol = nrow(as.data.frame(
                                  path.rpart(fit.tree, node = i, print.it = FALSE)
                                  )
                                  )-1)
                         )
    
    quiet(capture.output(for (j in 1:ncol(sub)){
      # Store each Rule as a Sub-population
      sub[,j] <- as.character(
        print(
          as.data.frame(
            path.rpart(fit.tree, node = i, print.it = FALSE))[j+1,1]
          )
        )
      # combine rule to one path 
      sub_pop <- noquote(paste(sub , collapse = " & "))
    }))
    
    # get subset 
    subset <- with(inference, inference[which(eval(parse(text = sub_pop))),])
    
    # Run the IV Regression
    if (length(unique(subset$w))!= 1 & length(unique(subset$z))!= 1 & nrow(subset) >2){
      
      # freq
      iv.reg <- ivreg(y ~ w | z, data = subset)
      
      # Bayes 
      bayes_iv <- brms_iv_function(subset, stan_model_first_stage, stan_model_second_stage)
                                   # get("stan_model_first_stage", envir = .GlobalEnv), 
                                   # get("stan_model_second_stage", envir = .GlobalEnv)) #
      
      ####   Step 5: Output the Values of each CCACE   ####
      
      bcfivMat[i,] <- iv_summary_func(iv.reg, subset, inference, sub_pop)
      
      bayes_ivMat[i, ] <- iv_summary_func(bayes_iv, subset, inference, sub_pop, bayes = T)
      
    }
    
    if (!(length(unique(subset$w))!= 1 & length(unique(subset$z))!= 1 & nrow(subset) >2)){
      bcfivMat[i,] <- c(sub_pop, 'est_problem', NA, NA, NA, NA, NA)
      bayes_ivMat[i, ] <- c(sub_pop, 'est_problem', NA, NA, NA, NA, NA)
    }
    
    # Delete data
    rm(subset)
  }
  
  # Adjust P.values 
  bcfiv_correction <- cbind(as.data.frame(bcfivMat), leaves)
  adj <-stats::p.adjust(as.numeric(bcfiv_correction$pvalue[which(bcfiv_correction$leaves==1)]),
                        paste(adj_method))
  Adj_pvalue <- rep(NA, length(rules))
  Adj_pvalue[which(bcfiv_correction$leaves==1)] <- adj
  
  # Store Results
  # changed by Jens, leaves denote End nodes
  bcfivResults <- cbind(as.data.frame(bcfivMat), leaves, Adj_pvalue)
  
  bayes_ivResults <- cbind(as.data.frame(bayes_ivMat), leaves)

  #### Return Results ####
  return(list('bcfivResults' = bcfivResults,
              'bayes_ivResults' = bayes_ivResults))
}
