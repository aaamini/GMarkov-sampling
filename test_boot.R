library(clime)
library(glasso)
library(tidyr)
library(tibble)
library(ggplot2)

source("viz.R")
source("utils.R")

# This is assumes Gam_seq exists in the environment -- e.g. run test_rand_sampling.R

# This is just the vector of lambdas for computing the regularization path (glasso calls lambda rho)
rholist = 10^seq(-2.5,0, len=50)

### Many regularization paths and the average
n = 50
Sigh_seq = lapply(Gam_seq, function(Gam) gen_data_from_Gam(n, Gam)$Sigh)  # generate new data from Gam_seq

reg_path = do.call(
  rbind, 
  lapply(seq_along(Sigh_seq[1:50]), function(rep) {
    Sigh = Sigh_seq[[rep]]
    Gam = Gam_seq[[rep]]
    Omegalist_glasso = glassopath(as.matrix(Sigh), rholist, trace = 0)$wi
    clime_out = clime(Sigh, lambda = rholist, sigma = T)
    
    data.frame(
      lambda = rholist,
      glasso = apply(Omegalist_glasso, 3, function(G) norm(G - Gam)), 
      clime = sapply(clime_out$Omegalist, function(G) norm(G - Gam)),
      n = n, 
      rep = rep
    )
  })
)

reg_path %>% 
  pivot_longer(c(glasso, clime), names_to = "method", values_to = "error") %>% 
  plot_paths_and_avg(error, alpha_range = c(0.2,1)) + 
  scale_x_continuous(trans = "log10") +
  ylim(c(0.25, 2.5)) +
  inline_legend(0.2,0.2) +
  xlab("Regularization parameter") +
  ylab("Operator norm error") 
  
ggsave(sprintf("reg_path_many_reps_n = %d.pdf", n), width = 5, height = 5)  


### Plot a single reg_path
Sigh = Sigh_seq[[1]]
Omegalist_glasso = glassopath(as.matrix(Sigh), rholist, trace = 0)$wi
clime_out = clime(Sigh, lambda = rholist, sigma = T)   

reg_path1 = data.frame(
  lambda = rholist,
  glasso = apply(Omegalist_glasso, 3, function(G) norm(G - Gam)), 
  clime = sapply(clime_out$Omegalist, function(G) norm(G - Gam))
)
reg_path1 %>%
  pivot_longer(-lambda, names_to = "method", values_to = "error") %>%
  ggplot(aes(lambda, error, color = method)) +
  geom_line(size = 1.2) +
  scale_x_continuous(trans = "log10") +
  ylab("Operator norm error") +
  theme_minimal()

ggsave("reg_path_single_rep.pdf", width = 5, height = 5)


# Optimal lambda
rholist[which.min(reg_path1$glasso)]
rholist[which.min(reg_path1$clime)]

### Simple test
lambda = 0.12
Gamh_glasso = glasso(as.matrix(Sigh), lambda)$wi
Gamh_clime = clime(Sigh, lambda = lambda, sigma = T)$Omega[[1]]

Gh_glasso = extract_graph(Gamh_glasso)
Gh_clime = extract_graph(Gamh_clime)

image(Matrix(Gamh_glasso))
image(Matrix(Gamh_clime))
 
image(Gh_glasso)
image(Gh_clime)
image(G)

### Systematic test
methods = list()
methods[["glasso"]] = function(Sigh) {
  extract_graph( glasso(as.matrix(Sigh), lambda)$wi  )
}

methods[["clime"]] = function(Sigh) {
  extract_graph( clime(Sigh, lambda = lambda, sigma = T)$Omega[[1]] )
}

n_methods = length(methods)


# n = 50
# # generate new data from Gam_seq
# Sigh_seq = lapply(Gam_seq, function(Gam) gen_data_from_Gam(n, Gam)$Sigh)  
# # apply the methods to get estimated graphs Gh[[method]][[i]]
# Gh_seq = lapply(methods, function(method) lapply(Sigh_seq, function(Sigh) method(Sigh)))
# 
# hist(comp_pairwise_dist(Gh_seq$glasso))
# # comp_dist_to_G0(Gh_seq$clime[-1], Gh_seq$clime[[1]])

runs = expand.grid(n = c(50, 100, 200))

res = do.call(
  rbind, 
  lapply(1:nrow(runs), function(j) {
    n = runs[j, "n"]
    
    # generate new data from Gam_seq
    Sigh_seq = lapply(Gam_seq, function(Gam) gen_data_from_Gam(n, Gam)$Sigh) 
    
    # apply the methods to get estimated graphs Gh[[method]][[i]]
    Gh_seq = lapply(methods, function(method) lapply(Sigh_seq, function(Sigh) method(Sigh)))
    
    do.call(
      rbind, 
      lapply(names(methods), function(method) 
        data.frame(pdist = comp_pairwise_dist(Gh_seq[[method]]), n = n, method = method)
    ))
  })
)

res %>% ggplot(aes(factor(n), pdist, fill = method)) + 
  geom_violin(color = NA) + 
  ylab("Pairwise Normalized Hamming Dist.") +
  xlab("Sample Size") + inline_legend(0.9,0.9)

ggsave(sprintf("var_plot_d = %d, lambda = %2.2f.pdf", d, lambda), width = 5, height = 5)
