no.scalefn <- function(x, ...){
  return(x)
}
no.rescalefn <- function(x, zs, ...){
  if ( is.null(zs) )
    return(matrix(NA, nrow = 0, ncol = 0))
  return(zs)
}



intpm<-function(ri){
(t(ri)%*%ri)}

norm2.Frnormfn <- function(m) {


  return(sum(apply(m,2,intpm))) # 
}





Larchetypes<-function (data, k, weights = NULL, maxIterations = 100, minImprovement = sqrt(.Machine$double.eps), 
    maxKappa = 1000, verbose = FALSE, saveHistory = FALSE, family = archetypesFamily("original"), ...) 
{
    mycall <- match.call()
    famargs <- list(...)
    memento <- NULL
    snapshot <- function(i) {
        a <- list(archetypes = as.archetypes(t(family$rescalefn(x, 
            family$undummyfn(x, zs))), k, alphas = t(alphas), 
            betas = t(betas), rss = rss, kappas = kappas, zas = t(family$rescalefn(x, 
                family$undummyfn(x, zas))), residuals = resid, 
            reweights = reweights, weights = weights, family = list(class = family$class)))
        memento$save(i, a)
    }
    printIter <- function(i) {
        cat(i, ": rss = ", formatC(rss, 8, format = "f"), ", improvement = ", 
            formatC(imp, 8, format = "f"), "\n", sep = "")
    }
    x1 <- t(data)
    x1 <- family$scalefn(x1, ...)
    x1 <- family$dummyfn(x1, ...)
    x0 <- family$globweightfn(x1, weights, ...)
    x <- x0
    n <- ncol(x)
    m <- nrow(x)
    init <- family$initfn(x, k, ...)
    betas <- init$betas
    alphas <- init$alphas
    zas <- NULL
    zs <- x %*% betas
    #resid <- zs %*% alphas - x
  resid <- zs[1: (nrow(zs)-1),] %*% alphas - x[1: (nrow(x)-1),]
    rss <- family$normfn(resid,...)/n
    reweights <- rep(1, n)
    kappas <- c(alphas = kappa(alphas), betas = kappa(betas), 
        zas = -Inf, zs = kappa(zs))
    isIll <- c(kappas) > maxKappa
    errormsg <- NULL
    if (saveHistory) {
        memento <- new.memento()
        snapshot(0)
    }
    i <- 1
    imp <- +Inf
    tryCatch(while ((i <= maxIterations) & (imp >= minImprovement)) {
        reweights <- family$reweightsfn(resid, reweights, ...)
        x <- family$weightfn(x0, reweights, ...)
        alphas <- family$alphasfn(alphas, zs, x, ...)
        zas <- family$zalphasfn(alphas, x, ...)
  resid1n <- zas[1: (nrow(zas)-1),] %*% alphas - x[1: (nrow(x)-1),]

        rss1 <- family$normfn(resid1n,...)/n
        kappas[c("alphas", "zas")] <- c(kappa(alphas), kappa(zas))
        betas <- family$betasfn(betas, x, zas, ...)
        zs <- x %*% betas
        kappas[c("betas", "zs")] <- c(kappa(betas), kappa(zs))
        alphas0 <- family$alphasfn(alphas, zs, x0, ...)
        #resid <- zs %*% alphas0 - x0
  resid <- zs[1: (nrow(zs)-1),] %*% alphas0 - x0[1: (nrow(x0)-1),]

        rss2 <- family$normfn(resid, ...)/n
        imp <- rss - rss2
        rss <- rss2
        kappas <- c(alphas = kappa(alphas), betas = kappa(betas), 
            zas = kappa(zas), zs = kappa(zs))
        isIll <- isIll & (kappas > maxKappa)
        if (verbose) 
            printIter(i)
        if (saveHistory) 
            snapshot(i)
        i <- i + 1
    }, error = function(e) errormsg <<- e)
    if (!is.null(errormsg)) {
        warning("k=", k, ": ", errormsg)
        return(as.archetypes(NULL, k, NULL, NA, iters = i, call = mycall, 
            history = history, kappas = kappas))
    }
    if (any(isIll)) 
        warning("k=", k, ": ", paste(names(isIll)[isIll], collapse = ", "), 
            " > maxKappa", sep = "")
    alphas <- family$alphasfn(alphas, zs, x1)
    betas <- family$betasfn(betas, x1, zs)
    zs <- family$undummyfn(x1, zs)
    zs <- family$rescalefn(x1, zs)
    resid <- zs %*% alphas - t(data)
    return(as.archetypes(t(zs), k, t(alphas), rss, iters = (i - 
        1), call = mycall, history = memento, kappas = kappas, 
        betas = t(betas), family = family, familyArgs = famargs, 
        residuals = t(resid), weights = weights, reweights = reweights, 
        scaling = attr(x1, ".Meta")))
}



LstepArchetypesMod <- function (data, k, nrep = 3, verbose = TRUE,saveHistory=FALSE){
  mycall <- match.call()
  as <- list()
  for (i in 1:length(k)) {
    as[[i]] <- list()
    class(as[[i]]) <- "repArchetypes"
    for (j in seq_len(nrep)) {
      if (verbose) 
        cat("\n*** k=", k[i], ", rep=", j, ":\n", sep = "")
      as[[i]][[j]] <- Larchetypes(data, k = k[i],saveHistory=FALSE, family = archetypesFamily("original",scalefn = no.scalefn, rescalefn = no.rescalefn,normfn=norm2.Frnormfn))
    }
  }
  return(structure(as, class = "stepArchetypes", call = mycall))
}
