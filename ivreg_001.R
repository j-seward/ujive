################################################################################
# Program: This is essentially a copy of Kolesar's Matlab function without 
#            the sparsity considerations.
# Author: Jonathan Seward
# Affiliation: Baylor University
# Date Created: 10/18/2021
# Last Modified: 3/9/2022
################################################################################

# pracma gives the matlab like matrix manipulations so I don't have to write my 
#   own. it is necessary for the function to work
# I think tidyverse is not necessary for the function, but I can't remember. I 
#   should check that out in the future.

pacman::p_load(pracma, tidyverse)

ivreg <- function (y, T, Z, W, noConstant = TRUE) {

  # need to test before coercion
  y <- as.matrix(y)
  T <- as.matrix(T)
  Z <- as.matrix(Z)
  W <- as.matrix(W)
  
  # Check that the matrices (T, Z, W) and left hand side (y) have compatible dimensions
  n = dim(Z)[1]
  K = dim(Z)[2]
  L = dim(W)[2]
  
  if (dim(y)[2] != 1 | dim(T)[2] != 1) {
    stop('Y and T must be vectors.')
  } else if (dim(y)[1] != n | dim(T)[1] != n | dim(W)[1] != n) {
    stop('The number of rows in Y, W, and T must equal the number of rows in Z.')
  } else if (dim(y)[2] == n) {
    warning('Y is a row vector--needs to be a column vector.')
    y = t(y)
  } else if (dim(T)[2] == n) {
    warning('T should be a column vector.')
    T = t(T)
  }
  
  # If matrix W does not include a constant, add it
  if (noConstant == TRUE) {
    Wc <- cbind(W, 1) # matrix W with a constant added
  
    # ## I dropped the sparsity logic -- may need to put back in if shift to sparse matrices
    # if ( Rank(as.matrix(t(Wc) %*% Wc)) > Rank(as.matrix(t(W) %*% W)) ) {
      W = Wc
      print('Added a row of ones to the matrix of covariates.')
#    }
  }
  
  # 2. Point estimates
  # annihilator
  MW <- function (x) {
    x - W %*% ( mldivide(W, x) )
  }
  
  Yp <- MW(cbind(y, T))
  Zp <- MW(Z)
  
  YY  <- as.matrix( t(Yp) %*% Yp) # [Y T]'*M_{W}*[Y T]: need this to be full to compute eigenvalues
  YPY <- (t(Yp) %*% Zp) %*% ( mldivide(Zp, Yp) ) # [Y T]'*H_{W}*[Y T]:
  YMY <- YY - YPY # ditto
  
  # 2.1 k-class: OLS, TSLS, LIML, MBTLS
  k <- rbind( 0, 1, min(eig( mrdivide(YY, YMY) )), 
              mrdivide( (1- mrdivide(L,n) ),
              ( mrdivide(1-(K - 1) , n - mrdivide(L, n) ) ) ) ) 
  
  beta <- ( YY[1,2]-k %*% YMY[1,2] )/( YY[2,2]-k %*% YMY[2,2])
  
  
  # 2.2 JIVE, UJIVE
  ZW  <- cbind(W, Z)
  MZW <- function(x) {
    x - ZW %*% ( mldivide(ZW, x) )
  }
  D <- rowSums( mrdivide( ZW, (t(ZW) %*% ZW ) ) * ZW ) # D=diag(P_ZW) as a vector
  iIDZW <- matrix(1, n) / (matrix(1, n)-D) # (I-D)^{-1}
  D <- rowSums( mrdivide(W , (t(W) %*% W) ) * W ) # D=diag(P_W) as a vector
  iIDW <- matrix(1, n) / (matrix(1, n)-D) # (I-D)^{-1}
  
  hatTujive <- T - iIDZW * MZW(T)
  hatPjive <- MW(hatTujive)
  hatPujive <- hatTujive - (T - iIDW * MW(T) )
  
  betaLabels <- c('OLS', 'TSLS', 'LIML', 'MBTSLS', 'JIVE', 'UJIVE', 'RTSLS')
  beta <- rbind(beta, 
               (t(hatPjive) %*% y) / (t(hatPjive) %*% T), 
               (t(hatPujive) %*% y) / (t(hatPujive) %*% T), 
               mrdivide(YPY[1, 1], YPY[1, 2]) )
  betas <- data.frame(t(beta))
  names(betas) <- betaLabels
  
  # 5. Standard Errors
    se <- matrix(NaN, 4, 7)
    epsilon <- 
      function (x) {
        Yp[,1] - Yp[,2]*x;
      }
    
    # 5.1 Homoscedastic
      sig <-
        function (x) {
          mrdivide(t(epsilon(x)) %*% epsilon(x), n)
        }

      se[1, 1:6] <-
        sqrt( c(mrdivide( sig(beta[1]), ( t(Yp[,2]) %*% Yp[,2] )),
                mrdivide( sig(beta[2]), YPY[2,2] ),
                mrdivide( sig(beta[3]), YPY[2,2] ),
                mrdivide( sig(beta[4]), YPY[2,2] ),
                mrdivide( sig(beta[5]) %*% (t(hatPjive)  %*% hatPjive),  (t(hatPjive)  %*% T)^2 ),
                mrdivide( sig(beta[6]) %*% (t(hatPujive) %*% hatPujive), (t(hatPujive) %*% T)^2)
                )
              )

      # 5.2 Heteroscedastic
      # ols
      se[2,1] <- mrdivide(sqrt(sum( (epsilon(beta[1]) * Yp[,2])^2 )), YY[2,2])

      # tsls, liml, mbtsls
      hatP <- Zp %*% mldivide(Zp, Yp[,2])
      sekclass <- function (x) {
        mrdivide(sqrt(sum( (epsilon(x) * hatP)^2 )), YPY[2,2])
      }
      se[2, 2:4] <- c( sekclass(beta[2]), sekclass(beta[3]), sekclass(beta[4]) )


      # jive
      se[2, 5] <- mrdivide( sqrt(sum((epsilon(beta[5]) * hatPjive)^2)), t(hatPjive) %*% T )
      se[2, 6] <- mrdivide( sqrt(sum((epsilon(beta[6]) * hatPujive)^2)), t(hatPujive) %*% T )

      # 5.3 Many instruments

      # Notation
      Sp <- YMY/(n-K-L) # S_{perp}
      S  <- YPY/n
      mmin <- min(eig(mldivide(Sp,S)))

      # Hessian of random-effects
      lamre <- max(eig(mldivide(Sp,S))) - K/n
      a <- rbind(beta[3], 1)
      b <- rbind(1, -beta[3])

      Omre <- (n-K-L)*Sp/(n-L) + n*(S - ( (lamre * a %*% t(a)) / c(mrdivide(t(a),Sp)  %*% a)))/(n-L)
      Qs <- mrdivide(t(b) %*% S %*% b , (t(b) %*% Omre %*% b))
      c  <- (lamre*Qs) /  ((1-L/n) * (K/n+lamre))
      se[3, 3] <-
        sqrt(
          mrdivide(
            (-t(b) %*% Omre %*% b) / ((n*lamre) * (lamre+K/n)),
            (Qs %*% Omre[2,2]-S[2,2] +
               mrdivide(c, (1-c) %*% mrdivide(Qs, (mrdivide(t(a),Omre)%*%a)))
             )
            )
          )


      # mbtsls, using maximum URE likelihood plug-in estimator
      b <- rbind(1, -beta[4]) # b_mbtsls
      Lam11 <- max(c(0, t(b) %*% (S-K/n * Sp) %*% b))
      if (mmin > K/n) {
        Lam22 <- S[2,2] - K/n * Sp[2,2]
        Omure <- Sp
      } else {
        Lam22 <- mrdivide(lamre, (mrdivide(t(a), Omre) %*% a) )
        Omure <- Omre
      }

      Gamma <- function (x) {
        matrix(rbind(c(1, 0),c(-x, 1)), 2,2)
      }
      Sig <- t(Gamma(beta[4])) %*% Omure %*% Gamma(beta[4])
      h <- (1-L/n) * (K-1)/n / (1-L/n-(K-1)/n)
      Vvalid <- Sig[1,1]/Lam22 + h*(Sig[1,1]*Sig[2,2] + Sig[1,2]^2)/Lam22^2
      Vinvalid <- Vvalid + (Lam11*Omure[2,2] + Lam11*Lam22*n/K)/Lam22^2
      se[3:4, 4] <- sqrt(c(Vvalid, Vinvalid)/n)



    #3. Other outputs
    F <- YPY[2, 2] / (K * Sp[2,2]) # first-stage F
    Xi <- YPY/n - (K/n) * Sp # Xi

    overid <- c(NaN, NaN)
    pvalue <- c(NaN, NaN)
    if (dim(Z)[2] > 1) {
      overid[1] <- n*mmin/(1-K/n-L/n+mmin) # n* J_sargan
      pvalue[1] <- 1 - pchisq(overid[1], K-1) # p-value for Sargan

      overid[2] <- n*mmin # Cragg-Donald
      pvalue[2] <- 1 - pnorm(sqrt((n-K-L)/(n-L))*
                               qnorm(pchisq(overid[2], K-1)))
    }

    stats <- list('BetaLabs' = betaLabels,
                  'Beta' = beta,
                  'SE' = se,
                  'F' = F,
                  'Omega' = Sp,
                  'Xi' = Xi,
                  'Sargan' = c(overid[1], pvalue[1]),
                  'CD' = c(overid[2], pvalue[2]))

    return(stats)
}

# outcome is outcome variable of interest
# fs_out is first stage outcome variable
# demogs is demographics/characteristics
# clin_fe is instrument dummies
# time_fe is time fixed effects

# this lumps all of the Xs together
# fe <- as.matrix(bind_cols(demogs, time_fe))

# ivreg(y, T, Z, W): the function returns a list of stats
#stats  <- ivreg(outcome, fs_out, clin_fe, time_fe)
#stats2 <- ivreg(outcome, fs_out, clin_fe, fe)
