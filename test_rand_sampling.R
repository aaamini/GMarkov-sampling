library(CVXR)
library(Matrix)

source("utils.R")

# Not very efficient to construct the problem every time using CVX, but should do for now
# Might be better to call the SCS solver directly with the data.
# get_problem_data(prob, 'SCS')$data 
find_epsmax = function(Del) {
  d = nrow(Del)
  eps = Variable(1)
  S = Variable(d, d, PSD =T) 

  con = list(S == Del*eps + diag(1, d))
  prob = Problem(Maximize(eps), con)
  res = solve(prob)
  res$value
}

normalize_infnorm = function(x) x / max(abs(x))

sample_Del = function(G) {
  sG = summary(G)
  g = length(sG$i)
  Del_values = normalize_infnorm(runif(g, min = -1, max = 1))
  Del = sparseMatrix(i = sG$i, j = sG$j, x = Del_values, symmetric = T)
  as(Del, "dgCMatrix")
}
make_prec_mat = function(eps, Del) {
  A = eps*Del
  diag(A) = 1
  A
}

sample_prec_mat = function(G) {
  Del = sample_Del(G)
  epsmax = find_epsmax(Del)
  eps = runif(1, min = 0, max = epsmax)
  Gam = make_prec_mat(eps, Del)
  list(Gam = Gam, Del = Del, eps = eps, epsmax = epsmax) 
}

set.seed(15)
d = 10
G = as(rsparsematrix(d, d, density = 0.2, symmetric = T), "nsparseMatrix")
diag(G) = F
image(G)

out = sample_prec_mat(G)
out$Gam
out$eps
out$epsmax
image(out$Del)
image(out$Gam)
eigen(out$Gam)$values
plot(eigen(out$Gam)$values)

# This is slow -- Generate Gam_seq (sequence of G-Markov precision matrices)
nreps = 100
Gam_seq = lapply(1:nreps, function(rep) sample_prec_mat(G)$Gam)
Gam_mean = Reduce(`+`, Gam_seq)/nreps
diag(Gam_mean) = 0
image(Gam_mean)

# summary(as(Gam_seq[[1]],"symmetricMatrix"))

Gam = Gam_seq[[1]]
Sigh = gen_data_from_Gam(50, Gam)$Sigh
image(Sigh)
image(Sig)
norm(Sigh - Sig)/norm(Sig)

# verify that the thresholding produces the same graph for all population level Gamma's
G_seq = lapply(Gam_seq, extract_graph, tol = 1e-4)
unique(G_seq)

