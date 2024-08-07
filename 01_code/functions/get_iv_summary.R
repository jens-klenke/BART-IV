get_iv_summary <- function(mod, inference, ...){
  
  # get summary
  summary <- summary(mod, diagnostics = TRUE)
  iv.effect <-  summary$coef[2,1]
  p.value <- summary$coef[2,4]
  p.value.weak.iv <- summary$diagnostics[1,4]
  proportion <- nrow(mod$model)/nrow(inference)
  # share of compilers
  compliers <- length(which(mod$model$z==mod$model$w))/nrow(mod$model)
  itt <- iv.effect*compliers
  
  # Store Results for Root
  return(
    c(deparse(substitute(mod)), round(iv.effect, 4), p.value, p.value.weak.iv,
      round(proportion, 4),round(itt, 4), round(compliers, 4))
  )
  
}