## try script for function heterogeneous_treatment_estimation  in '01_code/functions/heterogeneous_treatment_estimation.R' ----

# fit.tree, inference, adj_method, inference_tau

# load stan models ----
stan_model_first_stage <- readRDS(here::here("05_stan_code/brms_first_stage.rds"))
stan_model_second_stage <- readRDS(here::here("05_stan_code/brms_second_stage.rds"))

# Parameters ----
fit.tree <- bcf_fit.tree
fit.tree <- s_bcf_fit.tree
s_bcf_fit.tree
inference
adj_method
tau_true

# function body ----
# function call
#heterogeneous_treatment_estimation <- function(fit.tree, inference, adj_method)

# rules end terminal nodes?
rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))

# Initialize Outputs # NEW
bcfivMat <- tibble::tibble(
  "node" = rep(NA_character_, length(rules)),
  "CCACE" = rep(NA_real_, length(rules)),
  "pvalue" = rep(NA_real_, length(rules)),
  "Weak_IV_test" = rep(NA_real_, length(rules)),
  "Pi_obs" = rep(NA_real_, length(rules)),
  "ITT" = rep(NA_real_, length(rules)),
  "Pi_compliers" = rep(NA_real_, length(rules)),
  "pred" = rep(NA, length(rules))
)

bayes_ivMat <- tibble::tibble(
  "node" = rep(NA_character_, length(rules)),
  "CCACE" = rep(NA_real_, length(rules)),
  "CCACE_l-95%_CI" = rep(NA_real_, length(rules)),
  "CCACE_u-95%_CI" = rep(NA_real_, length(rules)),
  "Pi_obs" = rep(NA_real_, length(rules)),
  "ITT" = rep(NA_real_, length(rules)),
  "Pi_compliers" = rep(NA_real_, length(rules)),
  "pred" = rep(NA, length(rules))
)

pred_df <- tibble::tibble(
  index = as.numeric(row.names(inference)),
  tau_true = tau_true[as.numeric(row.names(inference))], 
  tau_pred = NA_real_)

# Generate Leaves (end notes) Indicator
lvs <- leaves <- numeric(length(rules)) 
lvs[unique(fit.tree$where)] <- 1
leaves[rules[lvs==1]] <- 1

#### Step 4: Run an IV Regression on each Node ####

# Run an IV Regression on the Root
iv.root <- ivreg(y ~ w | z,  
                 data = inference) # inference dataset

# bayes IV estimation
bayes_iv.root <- brms_iv_function(inference, stan_model_first_stage, stan_model_second_stage, i = 'root')
# weak instrument? 
# https://github.com/zeileis/ivreg/blob/main/R/ivregMethods.R

pred_root <- pred_df

# Store Results for Root
bcfivMat[1, ] <- iv_summary_func(iv.root, inference, inference, sub_pop = 'root', pred_root) # NEW
bayes_ivMat[1, ] <- iv_summary_func(bayes_iv.root, inference, inference, sub_pop = 'root', pred_root, bayes = T) # NEW

# Initialize New Data
names(inference) <- paste(names(inference), sep="")

# i = 1
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
  
  pred_subset <- pred_df %>%
    dplyr::filter(index %in% as.numeric(row.names(subset))) # get the right taus
  
  # Run the IV Regression
  if (length(unique(subset$w))!= 1 & length(unique(subset$z))!= 1 & nrow(subset) >2){
    # OLS 
    iv.reg <- ivreg(y ~ w | z,  
                    data = subset)
    
    # Bayes 
    bayes_iv <- brms_iv_function(subset, stan_model_first_stage, stan_model_second_stage, i = i)
    
    #### Step 5: Output the Values of each CCACE   ####
    bcfivMat[i,] <- iv_summary_func(iv.reg, subset, inference, sub_pop, pred_subset)
    
    bayes_ivMat[i, ] <- iv_summary_func(bayes_iv, subset, inference, sub_pop, pred_subset, bayes = T)
    
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
bcfivResults <- bcfivMat %>%
  dplyr::mutate('leaves' = leaves,
                'adj_pvalue' = Adj_pvalue)

bayes_ivResults <- bayes_ivMat %>%
  dplyr::mutate('leaves' = leaves)


#### Return Results ####
#  return(list('bcfivResults' = bcfivResults,
#              'bayes_ivResults' = bayes_ivResults))