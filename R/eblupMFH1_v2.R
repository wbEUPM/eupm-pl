df2matR <- 
  function (var.df, r) 
{
  if (dim(var.df)[2] != sum(1:r)) 
    stop("Length of vardir is not appropiate with data")
  var.df <- as.data.frame(var.df)
  n <- nrow(var.df)
  R <- lapply(var.df, diag)
  R_1n <- matrix()
  for (i in 1:r) {
    R.row <- R[[i]]
    for (j in i:r) {
      if (i != j) 
        R.row <- cbind(R.row, R[[sum((r - i):r) + j - 
                                   r]])
    }
    if (i == 1) {
      R_1n <- R.row
    }
    else {
      tmp <- matrix(rep(0, n * n * (i - 1)), n, n * (i - 
                                                       1))
      R.row <- cbind(tmp, R.row)
      R_1n <- rbind(R_1n, R.row)
    }
  }
  for (i in 1:(r * n)) {
    for (j in i:(r * n)) {
      if (R_1n[j, i] != R_1n[i, j]) {
        R_1n[j, i] <- R_1n[i, j]
      }
    }
  }
  return(R_1n)
}


adiag_v2 <- 
  
  function (..., pad = as.integer(0), do.dimnames = TRUE) 
{
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1], list(pad = pad)))
    return(do.call("Recall", c(list(args[[1]]), list(jj), 
                               list(pad = pad))))
  }
  a <- args[[1]]
  b <- args[[2]]
  if (is.null(b)) {
    return(a)
  }
  if (is.null(dim(a)) & is.null(dim(b))) {
    dim(a) <- rep(1, 2)
    dim(b) <- rep(1, 2)
  }
  if (is.null(dim(a)) & length(a) == 1) {
    dim(a) <- rep(1, length(dim(b)))
  }
  if (is.null(dim(b)) & length(b) == 1) {
    dim(b) <- rep(1, length(dim(a)))
  }
  if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
    stop("a and b must have identical number of dimensions")
  }
  s <- array(pad, dim.a + dim.b)
  s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
  ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) + 
                  dim.a[[i]])
  out <- do.call("[<-", c(list(s), ind, list(b)))
  n.a <- dimnames(a)
  n.b <- dimnames(b)
  if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
    dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
    names(dimnames(out)) <- names(n.a)
  }
  return(out)
}

eblupMFH1_v2 <- function (formula,
          vardir,
          samevar = FALSE,
          MAXITER = 100,
          PRECISION = 1e-04,
          data)
{
  r = length(formula)
  if (!missing(data)) {
    formula.matrix = lapply(formula, function(x) {
      model.frame(x, na.action = na.omit, data)
    })
    y.matrix = unlist(lapply(formula, function(x) {
      model.frame(x, na.action = na.omit, data)
    }[[1]]))
    x.matrix = Reduce(adiag_v2, lapply(formula, function(x) {
      model.matrix(x, data)
    }))
    n = length(y.matrix) / r
    if (!all(vardir %in% names(data)))
      stop("Object vardir is not appropiate with data")
    if (length(vardir) != sum(1:r))
      stop("Length of vardir is not appropiate with data")
    if (any(is.na(data[, vardir])))
      stop("Object vardir contains NA values.")
    R = df2matR(data[, vardir], r)
  }
  else {
    formula.matrix = lapply(formula, function(x) {
      model.frame(x, na.action = na.omit)
    })
    y.matrix = unlist(lapply(formula, function(x) {
      model.frame(x, na.action = na.omit)
    }[[1]]))
    x.matrix = Reduce(adiag_v2, lapply(formula, function(x) {
      model.matrix(x)
    }))
    n = length(y.matrix) / r
    if (dim(vardir)[2] != sum(1:r)) {
      stop("Object vardir is not appropiate with data")
    }
    if (any(is.na(vardir)))
      stop("Object vardir contains NA values.")
    R = df2matR(vardir, r)
  }
  for (i in 1:r) {
    if (attr(attributes(formula.matrix[[i]])$terms, "response") ==
        1)
      textformula = paste(formula[[i]][2], formula[[i]][1], formula[[i]][3])
    else
      textformula = paste(formula[[i]][1], formula[[i]][2])
    if (length(na.action(formula.matrix[[i]])) > 0) {
      stop("Argument formula= ", textformula, " contains NA values.")
    }
  }
  y.var = sapply(formula, "[[", 2)
  I = diag(n)
  Id = diag(r)
  omega = matrix(0, r, r)
  d.omega = list()
  d.Omega = list()
  for (i in 1:r) {
    d.omega[[i]] = omega
    d.omega[[i]][i, i] = 1
    d.Omega[[i]] = kronecker(d.omega[[i]], I)
  }
  convergence = TRUE
  if (samevar) {
    sigmau = median(diag(R))
    kit = 0
    diff = rep(PRECISION + 1, r)
    while (any(diff > PRECISION) & (kit < MAXITER)) {
      kit = kit + 1
      sigmau1 = sigmau
      G = kronecker(sigmau1, Id)
      GI = kronecker(G, I)
      Omega = solve(GI + R)
      Xto = t(Omega %*% x.matrix)
      Q = solve(Xto %*% x.matrix)
      P = Omega - t(Xto) %*% Q %*% Xto
      Py = P %*% y.matrix
      s = (-0.5) * sum(diag(P)) + 0.5 * (t(Py) %*% Py)
      iF = 0.5 * sum(diag(P %*% P))
      sigmau = sigmau1 + solve(iF) %*% s
      diff = abs((sigmau - sigmau1) / sigmau1)
    }
    sigmau = as.vector(mapply(max, sigmau, rep(0, r)))
    names(sigmau) = y.var
    if (kit >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }
    GI = kronecker(diag(sigmau), I)
    Omega = solve(GI + R)
    Xto = t(Omega %*% x.matrix)
    Qh = solve(Xto %*% x.matrix)
    P = Omega - t(Xto) %*% Qh %*% Xto
    Py = P %*% y.matrix
    beta = Qh %*% Xto %*% y.matrix
    res = y.matrix - x.matrix %*% beta
    eblup = x.matrix %*% beta + GI %*% Omega %*% res
    eblup.df = data.frame(matrix(eblup, n, r))
    names(eblup.df) = y.var
    se.b = sqrt(diag(Qh))
    t.val = beta / se.b
    pv = 2 * pnorm(abs(t.val), lower.tail = FALSE)
    coef = cbind(beta, se.b, t.val, pv)
    colnames(coef) = c("beta", "std.error", "t.statistics", "p.value")
    d = kronecker(Id, I) - GI %*% Omega
    gg1 = diag(GI %*% Omega %*% R)
    gg2 = diag(d %*% x.matrix %*% Qh %*% t(x.matrix) %*%
                 t(d))
    dg = Omega - GI %*% Omega %*% Omega
    g3 = (dg %*% (GI + R) %*% t(dg)) / iF
    gg3 = diag(g3)
    mse = gg1 + gg2 + 2 * gg3
    mse.df = data.frame(matrix(0, n, r))
    names(mse.df) = y.var
    for (i in 1:r)
      mse.df[, i] = mse[((i - 1) * n + 1):(i *
                                             n)]
  }
  else {
    sigmau = apply(matrix(diag(R), n, r), 2, median)
    kit = 0
    diff = rep(PRECISION + 1, r)
    while (any(diff > rep(PRECISION, r)) & (kit < MAXITER)) {
      kit = kit + 1
      sigmau1 = sigmau
      if (r == 1) {
        G = sigmau1
      }
      else {
        G = diag(as.vector(sigmau1))
      }
      GI = kronecker(G, I)
      Omega = solve(GI + R)
      Xto = t(Omega %*% x.matrix)
      Q = solve(Xto %*% x.matrix)
      P = Omega - t(Xto) %*% Q %*% Xto
      Py = P %*% y.matrix
      s = sapply(d.Omega, function(x)
        (-0.5) * sum(diag(P %*% x)) + 0.5 * (t(Py) %*% x %*% Py))
      iF = matrix(unlist(lapply(d.Omega, function(x)
        lapply(d.Omega, function(y)
          0.5 * sum(diag(
            P %*% x %*% P %*%
              y
          ))))), r)
      sigmau = sigmau1 + solve(iF) %*% s
      diff = abs((sigmau - sigmau1) / sigmau1)
      cat('Iteration: ', kit, ". Diff: ", as.vector(diff), "\n")
    }
    sigmau = as.vector(mapply(max, sigmau, rep(0, r)))
    names(sigmau) = y.var
    if (kit >= MAXITER && diff >= PRECISION) {
      convergence = FALSE
    }
    if (r == 1) {
      G = sigmau
    }
    else {
      G = diag(as.vector(sigmau))
    }
    GI = kronecker(G, I)
    Omega = solve(GI + R)
    Xto = t(Omega %*% x.matrix)
    Qh = solve(Xto %*% x.matrix)
    P = Omega - t(Xto) %*% Qh %*% Xto
    Py = P %*% y.matrix
    beta = Qh %*% Xto %*% y.matrix
    res = y.matrix - x.matrix %*% beta
    eblup = x.matrix %*% beta + GI %*% Omega %*% res
    eblup.df = data.frame(matrix(eblup, n, r))
    names(eblup.df) = y.var
    se.b = sqrt(diag(Qh))
    t.val = beta / se.b
    pv = 2 * pnorm(abs(t.val), lower.tail = FALSE)
    coef = cbind(beta, se.b, t.val, pv)
    colnames(coef) = c("beta", "std.error", "t.statistics", "p.value")
    FI = solve(iF)
    d = kronecker(Id, I) - GI %*% Omega
    gg1 = diag(GI %*% Omega %*% R)
    gg2 = diag(d %*% x.matrix %*% Qh %*% t(x.matrix) %*%
                 t(d))
    dg = lapply(d.Omega, function(x)
      x %*% Omega - GI %*%
        Omega %*% x %*% Omega)
    g3 = list()
    for (i in 1:r) {
      for (j in 1:r) {
        g3[[(i - 1) * r + j]] = FI[i, j] * (dg[[i]] %*%
                                              (GI + R) %*% t(dg[[j]]))
      }
    }
    gg3 = diag(Reduce("+", g3))
    mse = gg1 + gg2 + 2 * gg3
    mse.df = data.frame(matrix(0, n, r))
    names(mse.df) = y.var
    for (i in 1:r)
      mse.df[, i] = mse[((i - 1) * n + 1):(i *
                                             n)]
  }
  u.cap = GI %*% Omega %*% res
  u.cap.df = as.data.frame(matrix(u.cap, n, r))
  names(u.cap.df) = y.var
  result = list(
    eblup = NA,
    MSE = NA,
    randomEffect = NA,
    Rmatrix = NA,
    fit = list(
      method = NA,
      convergence = NA,
      iterations = NA,
      estcoef = NA,
      refvar = NA,
      informationFisher = NA
    )
  )
  result$eblup = signif(eblup.df, digits = 5)
  result$MSE = signif(mse.df, digits = 5)
  result$randomEffect = signif(u.cap.df, digits = 5)
  result$Rmatrix = signif(R, digits = 5)
  result$fit$method = "REML"
  result$fit$convergence = convergence
  result$fit$iterations = kit
  result$fit$estcoef = signif(coef, digits = 5)
  result$fit$refvar = signif(data.frame(t(sigmau)), digits = 5)
  result$fit$informationFisher = signif(iF, digits = 5)
  return(result)
}


eblupMFH2_v2 <- 
  
function (formula, vardir, MAXITER = 100, PRECISION = 1e-04, 
          data) 
{
  r = length(formula)
  if (!missing(data)) {
    formula.matrix = lapply(formula, function(x) {
      model.frame(x, na.action = na.omit, data)
    })
    y.matrix = unlist(lapply(formula, function(x) {
      model.frame(x, na.action = na.omit, data)
    }[[1]]))
    x.matrix = Reduce(adiag_v2, lapply(formula, function(x) {
      model.matrix(x, data)
    }))
    n = length(y.matrix)/r
    if (!all(vardir %in% names(data))) 
      stop("Object vardir is not appropiate with data")
    if (length(vardir) != sum(1:r)) 
      stop("Length of vardir is not appropiate with data")
    if (any(is.na(data[, vardir]))) 
      stop("Object vardir contains NA values.")
    vardir = data[, vardir]
    R = df2matR(vardir, r)
  }
  else {
    formula.matrix = lapply(formula, function(x) {
      model.frame(x, na.action = na.omit)
    })
    y.matrix = unlist(lapply(formula, function(x) {
      model.frame(x, na.action = na.omit)
    }[[1]]))
    x.matrix = Reduce(adiag_v2, lapply(formula, function(x) {
      model.matrix(x)
    }))
    n = length(y.matrix)/r
    if (dim(vardir)[2] != sum(1:r)) {
      stop("Object vardir is not appropiate with data")
    }
    if (any(is.na(vardir))) 
      stop("Object vardir contains NA values.")
    R = df2matR(vardir, r)
  }
  for (i in 1:r) {
    if (attr(attributes(formula.matrix[[i]])$terms, "response") == 
        1) 
      textformula = paste(formula[[i]][2], formula[[i]][1], 
                          formula[[i]][3])
    else textformula = paste(formula[[i]][1], formula[[i]][2])
    if (length(na.action(formula.matrix[[i]])) > 0) {
      stop("Argument formula= ", textformula, " contains NA values.")
    }
  }
  y.var = sapply(formula, "[[", 2)
  I = diag(n)
  Id = diag(r)
  d.Omega = list()
  sigmau = c(mean(sapply(vardir[, 1:r], median)), 1 - 1e-04)
  kit = 0
  diff = rep(PRECISION + 1, 2)
  convergence = TRUE
  while (any(diff > rep(PRECISION, 2)) & (kit < MAXITER)) {
    kit = kit + 1
    sigmau1 = sigmau
    Omega_AR = matrix(0, r, r)
    if (r == 1) {
      G = sigmau1[1]/(1 - sigmau1[2]^2)
    }
    else {
      for (i in 1:r) {
        for (j in 1:r) {
          Omega_AR[i, j] = sigmau1[2]^(abs(i - j))/(1 - 
                                                      sigmau1[2]^2)
        }
      }
      G = sigmau1[1] * Omega_AR
    }
    GI = kronecker(G, I)
    Omega = solve(GI + R)
    Xto = t(Omega %*% x.matrix)
    Q = solve(Xto %*% x.matrix)
    P = Omega - t(Xto) %*% Q %*% Xto
    d.Omega[[1]] = kronecker(Omega_AR, I)
    d.Omega[[2]] = matrix(NA, r, r)
    for (i in 1:r) {
      for (j in 1:r) {
        k = abs(i - j)
        d.Omega[[2]][i, j] = sigmau1[1] * (k * sigmau1[2]^(k - 
                                                             1) + (2 - k) * (sigmau1[2]^(k + 1)))/(1 - 
                                                                                                     sigmau1[2]^2)^2
      }
    }
    d.Omega[[2]] = kronecker(d.Omega[[2]], I)
    Py = P %*% y.matrix
    s = sapply(d.Omega, function(x) (-0.5) * sum(diag(P %*% 
                                                        x)) + 0.5 * (t(Py) %*% x %*% Py))
    iF = matrix(unlist(lapply(d.Omega, function(x) lapply(d.Omega, 
                                                          function(y) 0.5 * sum(diag(P %*% x %*% P %*% y))))), 
                2)
    sigmau = sigmau1 + solve(iF) %*% s
    if (abs(sigmau[2]) > 1) {
      sigmau[2] = sigmau1[2]
    }
    diff = abs((sigmau - sigmau1)/sigmau1)
    cat('Iteration: ', kit, ". Diff: ", as.vector(diff), "\n")
  }
  sigmau[1] = max(sigmau[1], 0)
  if (kit >= MAXITER && diff >= PRECISION) {
    convergence = FALSE
  }
  if (r == 1) {
    G = sigmau[1]/(1 - sigmau[2]^2)
  }
  else {
    for (i in 1:r) {
      for (j in 1:r) {
        Omega_AR[i, j] = sigmau[2]^(abs(i - j))/(1 - 
                                                   sigmau[2]^2)
      }
    }
    G = sigmau[1] * Omega_AR
  }
  GI = kronecker(G, I)
  Omega = solve(GI + R)
  Xto = t(Omega %*% x.matrix)
  Qh = solve(Xto %*% x.matrix)
  P = Omega - t(Xto) %*% Qh %*% Xto
  Py = P %*% y.matrix
  d.Omega[[1]] = kronecker(Omega_AR, I)
  d.Omega[[2]] = matrix(NA, r, r)
  for (i in 1:r) {
    for (j in 1:r) {
      k = abs(i - j)
      d.Omega[[2]][i, j] = sigmau[1] * (k * sigmau[2]^(k - 
                                                         1) + (2 - k) * (sigmau[2]^(k + 1)))/((1 - sigmau[2]^2)^2)
    }
  }
  d.Omega[[2]] = kronecker(d.Omega[[2]], I)
  beta = Qh %*% Xto %*% y.matrix
  res = y.matrix - x.matrix %*% beta
  eblup = x.matrix %*% beta + GI %*% Omega %*% res
  eblup.df = data.frame(matrix(eblup, n, r))
  names(eblup.df) = y.var
  se.b = sqrt(diag(Qh))
  t.val = beta/se.b
  pv = 2 * pnorm(abs(t.val), lower.tail = FALSE)
  coef = cbind(beta, se.b, t.val, pv)
  colnames(coef) = c("beta", "std.error", "t.statistics", 
                     "p.value")
  FI = solve(iF)
  d = kronecker(Id, I) - GI %*% Omega
  gg1 = diag(GI %*% Omega %*% R)
  gg2 = diag(d %*% x.matrix %*% Qh %*% t(x.matrix) %*% t(d))
  dg = lapply(d.Omega, function(x) x %*% Omega - GI %*% Omega %*% 
                x %*% Omega)
  g3 = list()
  for (i in 1:2) {
    for (j in 1:2) {
      g3[[(i - 1) * 2 + j]] = FI[i, j] * (dg[[i]] %*% 
                                            (GI + R) %*% t(dg[[j]]))
    }
  }
  gg3 = diag(Reduce("+", g3))
  mse = gg1 + gg2 + 2 * gg3
  mse.df = data.frame(matrix(0, n, r))
  names(mse.df) = y.var
  for (i in 1:r) mse.df[, i] = mse[((i - 1) * n + 1):(i * 
                                                        n)]
  u.cap = GI %*% Omega %*% res
  u.cap.df = as.data.frame(matrix(u.cap, n, r))
  names(u.cap.df) = y.var
  T.test = (sigmau[2, ])/(sqrt(FI[2, 2]))
  p.val = 2 * pnorm(abs(T.test), lower.tail = FALSE)
  rho.df = data.frame(signif(data.frame(sigmau[2], T.test, 
                                        p.val), digits = 5))
  names(rho.df) = c("rho", "T.test", "p-value")
  result = list(eblup = NA, MSE = NA, randomEffect = NA, Rmatrix = NA, 
                fit = list(method = NA, convergence = NA, iterations = NA, 
                           estcoef = NA, refvar = NA, rho = NA, informationFisher = NA))
  result$eblup = signif(eblup.df, digits = 5)
  result$MSE = signif(mse.df, digits = 5)
  result$randomEffect = signif(u.cap.df, digits = 5)
  result$Rmatrix = signif(R, digits = 5)
  result$fit$method = "REML"
  result$fit$convergence = convergence
  result$fit$iterations = kit
  result$fit$estcoef = coef
  result$fit$refvar = signif(sigmau[1, ], digits = 5)
  result$fit$rho = rho.df
  result$fit$informationFisher = signif(iF, digits = 5)
  return(result)
}
