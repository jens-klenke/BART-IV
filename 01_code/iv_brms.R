
#### first test ####

# Load necessary library
# if not already installed
install.packages("brms")
install.packages("AER")  
library(AER)
library(brms)

set.seed(123)  # For reproducibility

# Sample size
n <- 1000

# Generate the instrument (Z)
Z <- rnorm(n)

# Generate confounder (U)
U <- rnorm(n)

# Generate the treatment (X), affected by Z and U
X <- 0.5 * Z + 0.3 * U + rnorm(n)

# Generate the outcome (Y), affected by X and U
Y <- 2 * X + 0.5 * U + rnorm(n)

# Combine into a data frame
data <- data.frame(Y, X, Z, U)

# Estimate the treatment effect using IV (Z as instrument for X)
iv_model <- ivreg(Y ~ X | Z, data = data)

# Summarize the IV model
summary(iv_model)

# Perform Bayesian IV regression using brms
# First stage: regress X on Z
first_stage <- brm(X ~ Z, data = data, chains = 2, iter = 2000, warmup = 1000)

# Extract fitted values (predicted X)
data$X_hat <- fitted(first_stage)[,1]

# Second stage: regress Y on X_hat
second_stage <- brm(Y ~ X_hat, data = data, chains = 2, iter = 2000, warmup = 1000)

# Summarize the Bayesian IV model
summary(second_stage)


#### second test  ####

# bcf_iv <- function(
    # Data

# load packages 
############## Packages ################
source(here::here('01_code/packages.R'))
# some packages have to be implemented by hand (from github)

# source all files in the functions folder
invisible(sapply(list.files(here::here('01_code/functions'), full.names = TRUE), 
                 source))

dataset <- generate_dataset()

y <- dataset$y
w <- dataset$w
w1 <- dataset$w1
w0 <- dataset$w0
z <- dataset$z
x <- dataset$X

# https://evalf20.classes.andrewheiss.com/example/cace/ 
# visualize compliance structure - oracle (if we would know status through w1 w0)
dataset %>%
  as_tibble() %>%
  mutate(status=case_when((w1==1 & w0==0) ~ "Complier",
                          # because we have one-sided compliance assumption, 
                          # we do not have always-takers!!
                          (w1==1 & w0==1) ~ "Always-taker",
                          (w1==0 & w0==0) ~ "Never-taker", 
                          (w1==0 & w0==1) ~ "Defier")
  ) %>%
  ggplot(aes(y = y, x = as.factor(z))) + 
  geom_point(aes(color = status, shape=as.factor(w)),
             position = position_jitter(height = NULL, width = 0.25)) + 
  facet_wrap(vars(status)) + 
  labs(color = "Type of person", shape = "Observed compliance by w",
       x = "Control/Treatment assignment by random instrument z", y = "Outcome") +
  scale_color_viridis_d(option = "plasma", end = 0.85) +
  theme_bw()

# visualize compliance structure - observed
dataset %>%
  as_tibble() %>%
  mutate(status=case_when((w1==1 & w0==0) ~ "Complier",
                          # because we have one-sided compliance assumption, 
                          # we do not have always-takers!!
                          (w1==1 & w0==1) ~ "Always-taker",
                          (w1==0 & w0==0) ~ "Never-taker", 
                          (w1==0 & w0==1) ~ "Defier")
  ) %>%
  ggplot(aes(y = y, x = as.factor(w))) + 
  geom_point(aes(color = status, shape=as.factor(w)),
             position = position_jitter(height = NULL, width = 0.25)) + 
  facet_wrap(vars(as.factor(z))) + 
  labs(color = "Type of person", shape = "Observed compliance by w",
       x = "Observed compliance by w", y = "Outcome") +
  scale_color_viridis_d(option = "plasma", end = 0.85) +
  theme_bw()

# Parameters
binary = FALSE
n_burn = 1000
n_sim = 1000
inference_ratio = 0.5
max_depth = 2
cp = 0.01
minsplit = 30
adj_method = "holm"
seed = 42
cost = 'both'

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

######################################################
####     Continuous and Discrete Outcomes         ####
######################################################

#  if (binary == FALSE){

# Perform the Bayesian Causal Forest for the ITT
# is that \mu and ITT ? 
'bcf - sparse bcf' 
tictoc::tic()
bcf_itt.tree <- quiet(
  SparseBCF::SparseBCF(y[-index], z[-index], x[-index,], pihat = pihat,
                       nsim = n_sim, nburn = n_burn, sparse = F))
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
# s_bcf_exp %>% dplyr::filter(V2 == 0 & V3 == 0) %>% dplyr::summarise(median(s_bcf_tauhat))
s_bcf_exp <- as.data.frame(cbind(s_bcf_tauhat, x[-index,]))

# repair names? !! by me 
names(s_bcf_exp)[2:length(s_bcf_exp)] <- names(inference)[-(1:3)]

######################################################
####  Step 2: Build a CART on the Unit Level CITT ####
######################################################
'binary decision tree to discover, in an interpretable manner, the drivers of the heterogeneity??'
if(!cost){
  bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                        data = bcf_exp,
                        maxdepth = max_depth,
                        cp = cp,
                        minsplit = minsplit
  )
  # plot tree
  rpart.plot::rpart.plot(bcf_fit.tree)
  
  # binary tree for sparse trees 
  s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                          data = s_bcf_exp,
                          maxdepth = max_depth,
                          cp = cp,
                          minsplit = minsplit
  )
  # plot tree
  rpart.plot::rpart.plot(s_bcf_fit.tree)
}

# with cost function
if(cost){
  bcf_fit.tree <- rpart(bcf_tauhat ~ .,
                        data = bcf_exp,
                        maxdepth = max_depth,
                        cp = cp,
                        minsplit = minsplit,
                        cost = scales::rescale(
                          colMeans(bcf_itt.tree$varprb_tau), 
                          to = c(10, 1))
  )
  # plot tree
  rpart.plot::rpart.plot(bcf_fit.tree)
  
  # binary tree for sparse trees
  s_bcf_fit.tree <- rpart(s_bcf_tauhat ~ .,
                          data = s_bcf_exp,
                          maxdepth = max_depth,
                          cp = cp,
                          minsplit = minsplit,
                          cost = scales::rescale(
                            colMeans(s_bcf_itt.tree$varprb_tau), 
                            to = c(10, 1))
  )
  # plot tree
  rpart.plot::rpart.plot(s_bcf_fit.tree)
}


######################################################
####  Step 3: Extract the Causal Rules (Nodes)    ####
######################################################

fit.tree <- s_bcf_fit.tree
inference <- inference
adj_method <- adj_method

# rules end terminal nodes?
rules <- as.numeric(row.names(fit.tree$frame[fit.tree$numresp]))

# Initialize Outputs
bcfivMat <- as.data.frame(matrix(NA, nrow = length(rules), ncol =7))
names(bcfivMat) <- c("node", "CCACE", "pvalue", "Weak_IV_test", 
                     "Pi_obs", "ITT", "Pi_compliers")

# Initialize Outputs for Bayesian IV
bcfivMat_bayesian <- as.data.frame(matrix(NA, nrow = length(rules), ncol = 11))
names(bcfivMat_bayesian) <- c("node", "CCACE", "Est.Error", "l-95_CI", 
                     "u-95_CI", "Rhat", "Bulk_ESS", "Tail_ESS", "Pi_obs", "ITT", "Pi_compliers")

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
# in a bayesian manner
# Perform Bayesian IV regression using brms
# First stage: regress X on Z
first_stage <- brm(w ~ z, data = inference, chains = 2, iter = 2000, warmup = 1000)
# Extract fitted values (predicted X)
inference$w_hat <- fitted(first_stage)[,1]
# Second stage: regress Y on X_hat
second_stage <- brm(y ~ w_hat, data = inference, chains = 2, iter = 2000, warmup = 1000)
# Summarize the Bayesian IV model

summary <- summary(iv.root, diagnostics = TRUE)
iv.effect.root <-  summary$coef[2,1]
iv.effect.root.bayes <- fixef(second_stage)[2,1]
p.value.root <- summary$coef[2,4]
p.value.weak.iv.root <- summary$diagnostics[1,4]
proportion.root <- 1
# share of compilers in root
compliers.root <- length(which(inference$z==inference$w))/nrow(inference)
itt.root <- iv.effect.root*compliers.root
itt.root.bayes <- iv.effect.root.bayes*compliers.root

# Store Results for Root
bcfivMat[1,] <- c(NA , round(iv.effect.root, 4), p.value.root,
                  p.value.weak.iv.root, round(proportion.root, 4), 
                  round(itt.root, 4), round(compliers.root, 4))
bcfivMat_bayesian[1,] <- c(NA,
                           summary(second_stage)$fixed[2,],
                           round(proportion.root, 4), 
                           round(itt.root.bayes, 4),
                           round(compliers.root, 4))
                           


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
  # changes as both have to have more than just one value 
  if (length(unique(subset$w))!= 1 & length(unique(subset$z))!= 1){
    iv.reg <- ivreg(y ~ w | z,     # problem when alle are z = 1
                    data = subset)
    # Perform Bayesian IV regression using brms
    # First stage: regress X on Z
    first_stage <- brm(w ~ z, data = subset, chains = 2, iter = 2000, warmup = 1000)
    # Extract fitted values (predicted X)
    subset$w_hat <- fitted(first_stage)[,1]
    # Second stage: regress Y on X_hat
    second_stage <- brm(y ~ w_hat, data = subset, chains = 2, iter = 2000, warmup = 1000)
  
    summary <- summary(iv.reg, diagnostics = TRUE)
    iv.effect <-  summary$coef[2,1]
    iv.effect.bayes <- fixef(second_stage)[2,1]
    p.value <- summary$coef[2,4]
    p.value.weak.iv <- summary$diagnostics[1,4]
    compliers <- length(which(subset$z==subset$w))/nrow(subset)
    itt <- iv.effect*compliers
    itt.bayes <- iv.effect.bayes*compliers
    
    
    # Proportion of observations in the node
    proportion.node <- nrow(subset)/nrow(inference)
    
    ######################################################
    ####   Step 5: Output the Values of each CCACE   ####
    ######################################################
    
    bcfivMat[i,] <- c(sub_pop, round(iv.effect, 4), p.value, 
                      p.value.weak.iv, round(proportion.node, 4), 
                      round(itt, 4), round(compliers, 4))
    bcfivMat_bayesian[i,] <- c(sub_pop,
                               summary(second_stage)$fixed[2,],
                               round(proportion.node, 4), 
                               round(itt.bayes, 4),
                               round(compliers, 4))
    
  }
  
  # changed to display estimation problems
  if (!(length(unique(subset$w))!= 1 & length(unique(subset$z))!= 1)){
    bcfivMat[i,] <- c(sub_pop, 'est_problem', NA, NA, NA, NA, NA)
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
bcfivResults_bayes <- as.data.frame(bcfivMat_bayesian)











