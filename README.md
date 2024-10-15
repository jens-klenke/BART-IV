# High dimensional BART-IV

Based on [Bayesian Causal Forest with Instrumental Variable [BayesIV]](https://github.com/fbargaglistoffi/BCF-IV) by F.J. Bargagli-Stoffi, K. De Witte and G.

# Chrstophs remarks (meeting Oct 15 2024)

- What is actually computed in $\pi_c(x) = \mathbb{E} \left[W_i | Z_i = 1, X_i = x \right] - \mathbb{E} \left[W_i | Z_i = 0, X_i = x \right]$?
    -   Should the prediction mainly be 0 and 1s, or close to the compiler rate for each individual?
- Is there a benefit for Bayes-IV in the inference step?
- Why are Bayes-IV and frequentistic-IV so similar? Prior? Programming mistake?
- Where is our paper in the literature? Why do we need our paper? Is there a practical case?
