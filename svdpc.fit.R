### Modifications Francesca Rizzardi, 2016-07-14
### svdpc.fit.R: SVD PC fit algorithm
### $Id: ORIG-svdpc.fit.R 108 2007-03-19-2016-07-03 17:46:06Z bhm $

MYsvdpc.fit <- function(X, Y, ncomp, lbd=-Inf, stripped = FALSE, ...)
{
    Y <- as.matrix(Y)
    if (!stripped) {
        ## Save dimnames:
        dnX <- dimnames(X)
        dnY <- dimnames(Y)
    }
    ## Remove dimnames during calculation  (doesn't seem to matter; in fact,
    ## as far as it has any effect, it hurts a tiny bit in most situations).
    ## dimnames(X) <- dimnames(Y) <- NULL

    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]

    B <- array(0, dim = c(npred, nresp, ncomp))
    bsave <- array(0, dim = c(npred, nresp, ncomp))
    if (!stripped) fitted <- array(0, dim = c(nobj, nresp, ncomp))

    ## Center variables:
    Xmeans <- colMeans(X)
    X <- X - rep(Xmeans, each = nobj)
    Ymeans <- colMeans(Y)
    Y <- Y - rep(Ymeans, each = nobj)

    huhn <- La.svd(X)
    D <- huhn$d[1:ncomp]
    TT <- huhn$u[,1:ncomp, drop=FALSE] %*% diag(D, nrow = ncomp)
    Vt = huhn$vt[1:ncomp,, drop=FALSE] # ncomp x npred
    P <- t(Vt) # npred x ncomp
    tQ <- crossprod(TT, Y) / D^2

    for (a in 1:ncomp) {
        B[,,a] <- P[,1:a, drop=FALSE] %*% tQ[1:a,]
	if (!stripped) {
	    b <- (B[,,a, drop=FALSE] >= lbd)*1 # matrix, trick: converts T/F -> 1/0
	    # NOTE that drop=FALSE for B means that b is npred x resp x 1
	    bsave[,,a] = b[,,1]
	    B[,,a] <- B[,,a] * b
	    for (i in 1:nresp) {
		bv <- diag(b[,i,1]) # npred x npred
		Vb = Vt[1:a,] %*% bv %*% P[,1:a]
		fitted[,i,a] <- (TT[,1:a, drop=FALSE] %*% Vb %*% tQ[1:a,])[,i]
		}
	}
    }

    if (stripped) {
        ## Return as quickly as possible
        list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
    } else {
        residuals <- c(Y) - fitted
        fitted <- fitted + rep(Ymeans, each = nobj) # Add mean

        ## Add dimnames and classes:
        objnames <- dnX[[1]]
        if (is.null(objnames)) objnames <- dnY[[1]]
        prednames <- dnX[[2]]
        respnames <- dnY[[2]]
        compnames <- paste("Comp", 1:ncomp)
        nCompnames <- paste(1:ncomp, "comps")
        dimnames(TT) <- list(objnames, compnames)
        dimnames(P) <- list(prednames, compnames)
        dimnames(tQ) <- list(compnames, respnames)
        dimnames(B) <- list(prednames, respnames, nCompnames)
        dimnames(fitted) <- dimnames(residuals) <-
            list(objnames, respnames, nCompnames)
        names(D) <- compnames
        class(TT) <- "scores"
        R <- P                          # To avoid class "loadings" on projection
        class(P) <- class(tQ) <- "loadings"

        list(coefficients = B, b = bsave,
             scores = TT, loadings = P,
             Yloadings = t(tQ),
             projection = R,
             Xmeans = Xmeans, Ymeans = Ymeans,
             fitted.values = fitted, residuals = residuals,
             Xvar = D^2, Xtotvar = sum(X * X))
    }
}
