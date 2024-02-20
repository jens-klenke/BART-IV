# used insite BCF_IV-own estimation
extract_causal_rules <- function(fit.tree, inference){
  # rules end terminal nodes?
  rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))
  
  # Initialize Outputs
  bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol =7))
  names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", 
                       "Pi_obs", "ITT", "Pi_compliers")
  
  # Generate Leaves (end notes) Indicator
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
  bcfivMat[1,] <- c(NA , round(iv.effect.root, 4), p.value.root,
                    p.value.weak.iv.root, round(proportion.root, 4), 
                    round(itt.root, 4), round(compliers.root, 4))
  
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
      
      bcfivMat[i,] <- c(sub_pop, round(iv.effect, 4), p.value, 
                        p.value.weak.iv, round(proportion.node, 4), 
                        round(itt, 4), round(compliers, 4))
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
  
  # Return Results
  return(bcfivResults)
}