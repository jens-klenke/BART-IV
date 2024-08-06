iv_summary_func <- function(obj, subset, inference, sub_pop, bayes = FALSE, ...){
  
  proportion <- nrow(subset)/nrow(inference)
  compliers <- length(which(subset$z==subset$w))/nrow(inference)
  
  # frequencies 
  if(!bayes){
    summary <- summary(obj, diagnostics = TRUE)
    iv.effect <-  summary$coef[2,1]
    p.value <- summary$coef[2,4]
    p.value.weak.iv <- summary$diagnostics[1,4]
    itt <- iv.effect*compliers
    
    # Store Results for Root
    summary_vec <- c(sub_pop , round(iv.effect, 4), p.value,
                     p.value.weak.iv, round(proportion, 4),
                     round(itt, 4), round(compliers, 4))
  }
  # frequencies 
  if(bayes){
    summary_vec <- c(sub_pop, round(obj[2, 1], 4), round(obj[2, 3], 4), 
                     round(obj[2, 4], 4), round(proportion, 4), round(obj[2, 1]*compliers, 4), 
                     round(compliers, 4))
  }
  
  # return summary vector
  return(summary_vec)
  
}

## summary IV function
#iv_summary_func(iv.root, inference, inference, sub_pop = 'root', bayes = TRUE)
