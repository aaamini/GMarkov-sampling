gen_data_from_Gam = function(n, Gam) {
  Sig = solve(Gam)
  R = chol(Sig)
  
  X = matrix(rnorm(n*d), nrow = n) %*% R
  Sigh = t(X) %*% X / n
  list(X = X, Sigh = Sigh)
}


n_hammaing = function(X,Y) {
  d = nrow(X)
  sum(abs(X-Y)) / (d*(d-1))
}

comp_dist_to_G0 = function(G_seq, G0) {
  sapply(G_seq, function(G) n_hammaing(G,G0))    
}

comp_pairwise_dist = function(G_seq) {
  N = length(G_seq)
  # sapply(G_seq, function(G0) comp_dist_to_G0(G_seq, G0))
  unlist(lapply(1:(N-1), function(i) comp_dist_to_G0(G_seq[(i+1):N], G_seq[[i]])))
}

extract_graph = function(Gamh, tol = 1e-3) {
  Gh = as(abs(Gamh) > tol, "sparseMatrix")  
  diag(Gh) = F
  # as(Gh, "symmetricMatrix")
  Gh
}
