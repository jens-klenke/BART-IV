iv_summary_func <- function(obj, subset, inference, sub_pop, ...){
  summary <- summary(obj, diagnostics = TRUE)
  iv.effect <-  summary$coef[2,1]
  p.value <- summary$coef[2,4]
  p.value.weak.iv <- summary$diagnostics[1,4]
  proportion <- nrow(subset)/nrow(inference)
  # share of compilers in root
  compliers <- length(which(subset$z==subset$w))/nrow(inference)
  itt <- iv.effect*compliers
  
  # Store Results for Root
  c(sub_pop , round(iv.effect, 4), p.value,
    p.value.weak.iv, round(proportion, 4), 
                    round(itt, 4), round(compliers, 4))
}


## summary IV function
iv_summary_func(iv.root, inference, inference, sub_pop = 'root')
