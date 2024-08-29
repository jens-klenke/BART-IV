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
    
    # prediction and coverage 
    pred_df %<>%
      dplyr::mutate(tau_pred = iv.effect, 
                    coverage = coverage_95_fun(summary, .$tau_true))
      
    # node coverage 
    cove <- pred_df %>%
      dplyr::summarise(mean(coverage))
    
    print(cove)
    # Store Results
    summary_vec <- tibble::tibble(
      'node' = sub_pop, 
      'CCACE' = iv.effect, 
      'pvalue' = p.value,
      'Weak_IV_test' = p.value.weak.iv,
      'Pi_obs' = proportion, 
      'ITT' = itt,
      'Pi_compliers' = compliers, 
      'pred_df' = list(pred_df), 
      'node_coverage' = unlist(pred_df$coverage)
    )
    
  }
  # frequencies 
  if(bayes){
    
    # metrics
    pred_df %<>%
      dplyr::mutate(tau_pred = obj[2, 1],
                    coverage = coverage_95_fun(obj, .$tau_true, bayes = TRUE))
    
    # node coverage 
    cove <- pred_df %>%
      dplyr::summarise(mean(coverage))
    
    summary_vec <- tibble::tibble(
      "node" = sub_pop,
      "CCACE" = obj[2, 1],
      "CCACE_l-95%_CI" = obj[2, 3],
      "CCACE_u-95%_CI" = obj[2, 4],
      "Pi_obs" = proportion,
      "ITT" = obj[2, 1]*compliers,
      "Pi_compliers" = compliers,
      'pred_df' = list(pred_df),
      'node_coverage' = mean(pred_df$coverage)
    )
    
  }
  
  # compute metrics
  summary_vec$node_pehe <- PEHE_fun(pred_df$tau_pred, pred_df$tau_true)
  summary_vec$node_bias <- bias_fun(pred_df$tau_pred, pred_df$tau_true)
  summary_vec$node_abs_bias <- abs_bias_fun(pred_df$tau_pred, pred_df$tau_true)
  
  # return summary vector
  #return(summary_vec)
  
#}



  
# Oracl Algorithm 