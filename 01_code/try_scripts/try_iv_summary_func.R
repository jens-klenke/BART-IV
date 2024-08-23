#iv_summary_func <- function(obj, subset, inference, sub_pop, pred_df, bayes = FALSE, ...){
obj <- iv.root # IVREG
obj <- bayes_iv.root # BAYES
subset <- inference
inference <- inference 
sub_pop <- 'root'
pred_df <- pred_root
  
  proportion <- nrow(subset)/nrow(inference)
  compliers <- length(which(subset$z==subset$w))/nrow(inference)
  
  # frequencies 
  if(!bayes){
    summary <- summary(obj, diagnostics = TRUE)
    iv.effect <-  summary$coef[2,1]
    p.value <- summary$coef[2,4]
    p.value.weak.iv <- summary$diagnostics[1,4]
    itt <- iv.effect*compliers
    
    pred_df %<>%
      dplyr::mutate(tau_pred = iv.effect)
      
    # Store Results
    summary_vec <- tibble::tibble(
      'node' = sub_pop, 
      'CCACE' = iv.effect, 
      'pvalue' = p.value,
      'Weak_IV_test' = p.value.weak.iv,
      'Pi_obs' = proportion, 
      'ITT' = itt,
      'Pi_compliers' = compliers, 
      'pred_df' = list(pred_df)
    )
    
  }
  # frequencies 
  if(bayes){
    
    pred_df %<>%
      dplyr::mutate(tau_pred = obj[2, 1])
    
    summary_vec <- tibble::tibble(
      "node" = sub_pop,
      "CCACE" = obj[2, 1],
      "CCACE_l-95%_CI" = obj[2, 3],
      "CCACE_u-95%_CI" = obj[2, 4],
      "Pi_obs" = proportion,
      "ITT" = obj[2, 1]*compliers,
      "Pi_compliers" = compliers,
      'pred_df' = list(pred_df)
    )
    
  }
  
  # return summary vector
  #return(summary_vec)
  
#}
