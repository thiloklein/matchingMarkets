# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) Karlson-Holm-Breen method for comparing probit coefficients
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Karlson-Holm-Breen method for comparing probit coefficients
#'
#' @description Significance test for confounding; that is, the difference between regression 
#' coefficients from same-sample nested logit and probit models. The test procedure follows
#' Karlson et al (2012), Section 3.4.
#'
#' @param X data frame comprising independent variables including confounding variable.
#' @param y vector of dependent variable.
#' @param z character string giving the name of the confounding variable in \code{X}.
#' 
#' @export
#' 
#' @import stats
#' 
#' @author Thilo Klein 
#' 
#' @keywords summary
#' 
#' @references Karlson, K.B., A. Holm and R. Breen (2012). Comparing regression coefficients between same-sample nested models using logit and probit: A new method. \emph{Sociological Methodology}, 42(1):286--313.
#' 
#' @examples
#' ## 1. load results from Klein (2015a)
#'  data(klein15a)
#'  M <- klein15a$model.list
#' 
#' ## 2. extract variables
#'  X <- do.call(rbind.data.frame, M$X)
#'  eta <- c(klein15a$coefs$eta, rep(0, length(M$X)-length(M$W)))
#'  X <- cbind(X,eta)
#'  y <- unlist(M$R)
#' 
#' ## 3. apply KHB method
#' khb(X=X, y=y, z="eta")
khb <- function(X,y,z){

  ## -------------------------------------------------------
  ## Karlson-Holm-Breen (2012). Sociological Methodology:  
  ## Comparing Regression Coefficients Between Same-sample 
  ## Nested Models Using Logit and Probit: A New Method. 
  
  ## Arguments:
  ## X : design matrix comprising independent variables including confounding variable (z)
  ## y : vector of dependent variable
  ## z : confounding variable
  ## -------------------------------------------------------
  
  ## --- Reduced and Full model, Equations (1) and (2) on page 289 ---
  glmR <- glm(y ~ -1 + ., family=binomial(link="probit"),data=X[,-which(names(X)==z)])
  glmF <- glm(y ~ -1 + ., family=binomial(link="probit"),data=X)

  ## --- Auxiliary regression, Equation (8) on page 292 ---
  lmA <- lm(X$eta ~ -1 + ., data=X)
  glmFs <- glm(y ~ -1 + . + lmA$resid, family=binomial(link="probit"),data=X[,-which(names(X)==z)])

  ## --- Recovery of parameters from Full model and Auxiliary regression ---
  b.yx.zt <- glmFs$coef[!names(glmFs$coef)%in%c("(Intercept)","lmA$resid")]
  b.yx.z <- glmF$coef[!names(glmF$coef)%in%c("(Intercept)",z)]
  b.yz.x <- glmF$coef[names(glmF$coef)==z]
  sigma.b.yz.x <- diag(vcov(glmF))[names(diag(vcov(glmF)))==z]
  
  if(sum(X[,1])!=dim(X)[1]){ ## X doesn't have an intercept
    t.zx <- lmA$coef
    sigma.t.zx <- diag(vcov(lmA))
  } else{ ## X has an intercept
    t.zx <- lmA$coef[-1]
    sigma.t.zx <- diag(vcov(lmA))[-1]
  }
  
  ## --- Significance Test (page 295-296) ---
  
  ## Equation (17)
  round(b.yx.zt - b.yx.z,5) == round(b.yz.x * t.zx,5)
  
  ## Test statistic, Equation (20)
  Z <- (b.yx.zt - b.yx.z)/sqrt((b.yz.x^2 * sigma.t.zx^2 + t.zx^2 * sigma.b.yz.x^2))
  p.val <- round((1-pnorm(Z)),4)
  
  cat("\nKarlson-Holm-Breen method\nNull hypothesis: Change in coefficient is not attributable to confounding by z.\n\n") 
  data.frame(p.value=ifelse(p.val<1e-04,"<1e-04",p.val))
}
