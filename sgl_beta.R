sgl_beta <- function(X, y, group) {
  dat = list(x = X, y = y)
  mod <- cvSGL(dat, group, type="linear", min.frac = 0.005)
  min_lam = mod$lambdas[which.min(mod$lldiff)]
  fit_sgl = SGL(dat, grps, type="linear", lambdas = c(min_lam))
  beta_sgl = fit_sgl$beta
  return(beta_sgl)
}