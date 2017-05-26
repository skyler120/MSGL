lasso <- function(X, y, group=NA) {
  glmnet::glmnet(X,y, lambda=glmnet::cv.glmnet(X,y)$lambda.min)$beta
}