library(pROC)
evaluate_beta <- function(beta, beta_true,y_hat,y,tol=1e-6){
  beta_close = sum(abs(beta-beta_true)<tol)/length(beta)
  myROC = roc(y, y_hat)
  return(list(prop = beta_close, sens = myROC$sensitivities, 
              spec = myROC$specificities,auc = myROC$auc))
}