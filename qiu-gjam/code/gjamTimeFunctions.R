colF <- colorRampPalette( c('#8c510a','#d8b365','#c7eae5','#5ab4ac','#01665e','#2166ac') )

colT <- colorRampPalette( c('#8c510a','#d8b365','#c7eae5','#5ab4ac','#01665e','#2166ac') )

mergeList <- function( list1, list2 ){
  
  # update elements of list1 if contained list2
  
  stopifnot(is.list(list1), is.list(list2))
  
  n2 <- names(list2)
  
  for(k in n2){
    if(!k %in% names(list1)){
      list1 <- append( list1, list2[k] )
    }else{
      list1[k] <- list2[k]
    }
  }
  list1
}

.pasteCols <- function(mm){
  tmp <- apply(mm,1,paste0,collapse='-')
  names(tmp) <- NULL
  tmp
}

buildGuildTraits <- function(ydata, traitTable, traitNames, 
                             traitCode = 'code6', minCooccur=5,
                             diagTrait = -1){
  # competitive guilds
  # traitCode  - column in traitTable that matches colnames of ydata
  # traitNames - columns in traitTable holding traits used to define guilds
  #            - columns are characters, factors, or integers
  
  
  nt   <- length(traitNames)
  levs <- vector('list', nt )   # levels for each trait
  names(levs) <- traitNames
  ynames <- colnames(ydata)
  
  mm <- match(ynames, traitTable[,traitCode])
  wna <- which(is.na(mm))
  if(length(wna) > 0)stop( paste0('missing from traitTable: ', ynames[wna],collapse=', ') )
  
  traitTab <- traitTable[mm ,]
  
  for(k in 1:nt){
    
    tk <- traitTab[,traitNames[k]]
    tk[is.na(tk)] <- 'UNKN'
    
    levs[[k]] <- sort(unique(tk))
    traitTab[,traitNames[k]] <- tk
  }
  groups <- expand.grid( levs )
  
  gnames <- apply( groups, 1, paste0, collapse = '_' )
  mm     <- match(ynames, traitTab[,traitCode])
  tdata  <- traitTab[mm, traitNames]
  
  tnames <- apply( tdata, 1, paste0, collapse = '_' )
  tcodes <- match(tnames, gnames)
  
  nall   <- sort(unique(tnames))
  tall   <- sort(unique(tcodes))
  
  ng <- length(nall)
  guildList <- vector( 'list', ng  )
  names(guildList) <- nall
  
  # below co-occurence threshold
  ty <- as.matrix( ydata ) 
  ty[ ty != 0 ] <- 1
  tt <- crossprod(ty)
  tt[ tt < minCooccur ] <- 0
  tt[ tt > 0 ] <- 1
  
  for(i in 1:length(tall)){          # not in same guild
    ii <- tall[i]
    wi <- which(tcodes == ii)
    wj <- which(tcodes != ii)
    #  tii <- tt[wi,wi, drop=F]
    #  tij <- tt[wi,wj, drop=F]
    
    tt[wi,wj] <- 0
    tt[wj,wi] <- 0
    guildList[[i]] <- wi
  }
  #  diag(tt) <- diagTrait
  zeros <- which(tt == 0, arr.ind=T)
  ones  <- which(tt == 1, arr.ind=T)
  fromTo <- rbind( ones, ones[,2:1])
  
  g  <- sort(unique(ones[,1]))
  gl <- numeric(0)
  for(k in 1:length(g)){
    ik <- which(ones[,1] == g[k])
    jk <- ones[ik, 2]
    jk <- sort( unique( c(ik, jk) ) )
    gl <- append(gl, list(jk))
  }
  
  attr(guildList, 'traitNames') <- traitNames
  
  list(alpha = tt, zeroAlpha = zeros, guildList = guildList,
       fromTo = fromTo)
}

gjamPredTime <- function(output, nsim = 10, quant = .05, minw = .1){
  
  # evaluates carrying capacity for TIME
  # must supply either 'groups' or 'ngrid'
  # groups - xdata$groups column to predict for groups
  # ngrid  - no. values per covariate if !BYGROUP
  # nsim   - simulation steps per covariate combination
  # NOTE: x is centered and standardized
  
  #  breaks <- optim <- rhoStandXmu <- sensAlpha <- sensBeta <- sensRho <- 
  #    wA <- wL <- xStand <- xUnstand <- NULL
  
  x          <- output$inputs$xStand
  xl         <- output$inputs$xRho
  #xnames     <- output$inputs$xnames
  xlnames    <- output$inputs$xlnames
  factorBeta <- output$inputs$factorBeta
  factorRho  <- output$inputs$factorRho
  ydata      <- output$inputs$y
  effMat     <- output$inputs$effMat
  uindex     <- output$parameters$uindex
  gindex     <- output$parameters$gindex
  other      <- output$inputs$other
  notOther   <- output$inputs$notOther
  bg         <- output$parameters$betaMu
  xnames     <- rownames(bg)
  Rmat       <- output$parameters$RmatStandXmu
  Amat       <- output$parameters$Amu
  sigMu      <- output$parameters$sigMu
  wB         <- output$parameters$wB
  wL         <- output$parameters$wL
  wA         <- output$parameters$wA
  bgibbs     <- output$chains$bgibbs
  lgibbs     <- output$chains$lgibbs
  alphaGibbs <- output$chains$alphaGibbs
  groups     <- output$inputs$xdata$groups
  ng         <- output$modelList$ng
  burnin     <- output$modelList$burnin
  timeList   <- output$inputs$timeList
  groupCol   <- output$inputs$timeList$groups
  timeCol    <- output$inputs$timeList$times
  group      <- output$inputs$xdata[,groupCol]
  # time       <- output$inputs$xdata[,timeCol]
  assign(timeCol, time)
  time       <- output$inputs$xdata$times
  timeZero   <- output$inputs$timeList$timeZero
  # year       <- output$inputs$xdata[,groupCol]
  
  cx <- colnames(xl)[!colnames(xl) %in% colnames(x)]
  if(length(cx) > 0)x <- cbind(x, xl[,cx])
  
  S <- ncol(ydata)
  snames <- colnames(ydata)
  termB  <- termR <- termA <- FALSE
  
  if(!is.null(bgibbs))termB <- TRUE
  if(!is.null(lgibbs))termR <- TRUE
  if(!is.null(alphaGibbs))termA <- TRUE
  
  w <- ydata/effMat
  
  # wstart <- apply( ydata/effMat, 2, mean, na.rm=T)
  
  fn <- function(wt, termA, termB, termR, xmu, xlmu,
                 gindex, uindex, notOther,
                 amat, bmat, rmat, epsilon){
    wz <- wt
    wz[ wz < 0 ] <- 0
    ff <- wz*0
    
    if(termB)ff[,notOther] <- t(xmu[,rownames(bmat)]%*%bmat[,notOther])
    if(termR){
      Vmat <- wz[ drop=F, ,gindex[,'colW']]*xlmu[ drop=F, ,gindex[,'colX']]
      colnames(Vmat) <- rownames(gindex)
      ff[,notOther] <- ff[,notOther] + Vmat%*%rmat[colnames(Vmat),notOther]
    }
    if(termA){
      Umat <- wz[ drop=F, ,uindex[,1]]*wz[ drop=F, ,uindex[,2]]
      ff[,notOther] <- ff[,notOther] + Umat%*%amat[,notOther]
    }
    ff[,notOther] <-  ff[,notOther] + wz[,notOther] + epsilon
    ff
  }
  
  
  ii <- burnin:ng
  
  times  <- sort(unique(time))
  groups <- sort(unique(group))
  nt     <- length(times)
  ng     <- length(groups)
  
  pbar <- txtProgressBar(min = 1,max = nt, style=1)
  
  ylo <- yhi <- ymu <- w*0
  
  yj <- array(NA, dim = c( ng, ncol(w), nsim) )
  dimnames(yj)[[1]] <- groups
  dimnames(yj)[[2]] <- colnames(w)
  
  wj <- w[timeZero+1,] + .01
  for(k in 1:nsim)yj[,,k] <- wj
  
  for(j in 1:nt){
    
    gj <- which(time == j)  # x[gj,] is one step ahead of w[timeZero,]
    gi <- group[gj]
    
    zj <- yj <- yj[gi,,] 
    zj[ zj < 0 ] <- 0
    zj <- zj + minw
    
    if(j> 1){
      
      for(k in 1:nsim){
        
        ig <- sample(ii, 1)   
        
        muw <- w*0
        
        if(termB){
          bmat <- bg*0
          bmat[ 1:length(bg) ] <- bgibbs[ig,]
        }
        
        if(termA){
          amat <- Amat*0
          amat[ wA ] <- alphaGibbs[ig,]
        }
        
        if(termR){
          rmat <- Rmat*0
          rmat[ wL ] <- lgibbs[ig,]
        }
        
        epsilon <- .rMVN(1, 0, sigma = sigMu )[notOther]
        
        yj[,,k] <- fn(wt = zj[ , ,k], termA, termB, termR, 
                      xmu = x[drop=F, gj, ], xlmu = xl[drop=F, gj, ],
                      gindex, uindex, notOther,
                      amat, bmat, rmat, epsilon)
      }
    }
    
    yj[ yj < 0 ] <- 0
    ymu[gj,] <- apply( yj, c(1,2), mean )
    ylo[gj,] <- apply( yj, c(1,2), quantile, quant )
    yhi[gj,] <- apply( yj, c(1,2), quantile, 1 - quant )
    
    setTxtProgressBar(pbar,j)
    
  }
  
  ymu <- ymu*effMat[,notOther]
  ylo <- ylo*effMat[,notOther]
  yhi <- yhi*effMat[,notOther]
  
  list(ymu = ymu, ylo = ylo, yhi = yhi)
}
##########################


##########################

.wrapperEquilAbund <- function(output, covars = NULL, nsim = 10, ngrid = NULL, 
                               BYFACTOR = FALSE, BYGROUP = TRUE, verbose = FALSE){
  
  # evaluates equil abundances for TIME
  # covars - names of covariates in x to use as gradients
  # must supply 'ngrid' or BYGROUP
  # groups - xdata$groups column to predict for groups
  # ngrid  - no. values per covariate if !BYGROUP
  # nsim   - simulation steps per covariate combination
  # NOTE: x is centered and standardized
  
  #  breaks <- optim <- rhoStandXmu <- sensAlpha <- sensBeta <- sensRho <- 
  #    wA <- wL <- xStand <- xUnstand <- NULL
  
  x          <- output$inputs$xStand
  xl         <- output$inputs$xRho
  xnames     <- output$inputs$xnames
  xlnames    <- output$inputs$xlnames
  factorBeta <- output$inputs$factorBeta
  factorRho  <- output$inputs$factorRho
  ydata      <- output$inputs$y
  effMat     <- output$inputs$effMat
  uindex     <- output$parameters$uindex
  gindex     <- output$parameters$gindex
  other      <- output$inputs$other
  notOther   <- output$inputs$notOther
  bg         <- output$parameters$betaMu
  Rmat       <- output$parameters$RmatStandXmu
  Amat       <- output$parameters$Amu
  sigMu      <- output$parameters$sigMu
  wB         <- output$parameters$wB
  wL         <- output$parameters$wL
  wA         <- output$parameters$wA
  bgibbs     <- output$chains$bgibbs
  lgibbs     <- output$chains$lgibbs
  alphaGibbs <- output$chains$alphaGibbs
  groups     <- output$inputs$xdata$groups
  ng         <- output$modelList$ng
  burnin     <- output$modelList$burnin
  
  if( !is.null(ngrid) )BYGROUP <- F
  
  cx <- colnames(xl)[!colnames(xl) %in% colnames(x)]
  if(length(cx) > 0)x <- cbind(x, xl[,cx])
  
  S <- ncol(ydata)
  snames <- colnames(ydata)
  termB  <- termR <- termA <- FALSE
  
  if(!is.null(bgibbs))termB <- TRUE
  if(!is.null(lgibbs))termR <- TRUE
  if(!is.null(alphaGibbs))termA <- TRUE
  
  
  xnames  <- colnames(x)
  vxnames <- xnames[ !xnames %in% factorBeta$isFactor ][-1] # exclude intercept
  vxnames <- vxnames[ !vxnames %in% factorRho$isFactor ]
  
  gnames <- vxnames
  
  if( !is.null(covars) )gnames <- vxnames[ vxnames %in% covars ]
  
  ix <- grep(':', xnames)
  intb <- xnames[ix]
  if (length(ix) > 0)xnames <- xnames[ -ix ]
  
  ix <- grep(':', xlnames)
  intr <- xlnames[ix]
  if (length(ix) > 0)xlnames <- xlnames[ -ix ]
  
  ix <- grep(':', vxnames)
  intv <- vxnames[ix]
  if (length(ix) > 0)vxnames <- vxnames[ -ix ]
  
  if( !BYGROUP ){   # a prediction grid for covariates
    
    sgrid <- 0
    if(ngrid > 1)sgrid <- seq(-2, 2, length = ngrid)  # std devs for standardized x
    
    grid <- vector( 'list', length(vxnames) )
    for(k in 1:length(vxnames)){
      if(vxnames[k] %in% gnames){
        grid[[k]] <- sgrid
      }else{
        grid[[k]] <- 0
      }
    }
    names(grid) <- vxnames
    xgrid <- expand.grid(grid)
    
  }else{                # predict for groups
    ii <- rep( groups, ncol(x) )
    jj <- rep( colnames(x), each = nrow(x) )
    xgrid <- tapply( as.vector(x), list(groups = ii, x = jj), mean, na.rm=T)
    xgrid <- xgrid[,vxnames, drop=F]
  }
  
  fnames <- unique( c(factorBeta$facNames,  factorRho$facNames) )
  
  nfact <- length(fnames)
  factorColumns <- character(0)
  
  if(nfact > 0){
    
    flist  <- vector( 'list', length(fnames))
    
    for(k in 1:length(fnames)){ # find factor in either factorBeta or factorRho
      
      
      fk <- fl <- numeric(0)
      knames <- lnames <- character(0)
      
      if ( !is.null(factorBeta$factorList) ) {
        wf <- grep(fnames[k], names(factorBeta$factorList))
        if(length(wf) > 0){
          fk <- factorBeta$contrast[fnames[k]][[wf]]
          knames <- factorBeta$factorList[fnames[k]][[1]]
          fcc   <- factorBeta$factorList[[wf]]
        }
      }
      if (!is.null(factorRho$factorList) ) {
        wf <- grep(fnames[k], names(factorRho$factorList))
        if(length(wf) > 0){
          fk <- factorRho$contrast[fnames[k]][[1]]
          knames <- factorRho$factorList[fnames[k]][[1]]
          fcc   <- factorRho$factorList[[wf]]
        }
      }
      if(!BYFACTOR)fk <- fk[drop = F, 1, ]
      
      xindex <- rep( 1:nrow(xgrid), each = nrow(fk) )
      kindex <- rep( 1:nrow(fk), nrow(xgrid) )
      
      colnames(fk) <- fcc
      
      xgrid <- cbind(xgrid[drop=F, xindex,], fk[drop=F, kindex,])
      
      factorColumns <- c( factorColumns, fcc )
    }
  }
  
  
  intercept <- 1
  xgrid <- as.matrix( cbind(intercept, xgrid) )
  xgrid <- xgrid[drop=F, ,xnames]
  
  attr(xgrid,'factors') <- factorColumns
  
  wstart <- apply( ydata/effMat, 2, mean, na.rm=T)
  
  fn <- function(w, termA, termB, termR, xmu, xlmu,
                 gindex, uindex, notOther,
                 amat, bmat, rmat, epsilon){
    wz <- w
    wz[ wz < 0 ] <- 0
    ff <- wz*0
    
    if(termB)ff[notOther] <- t(xmu[,rownames(bmat)]%*%bmat[,notOther])
    if(termR){
      Vmat <- wz[ gindex[,'colW']]*xlmu[ gindex[,'colX']]
      names(Vmat) <- rownames(gindex)
      ff[notOther] <- ff[notOther] + t(Vmat%*%rmat[names(Vmat),notOther])
    }
    if(termA){
      Umat <- matrix( wz[uindex[,1]]*wz[uindex[,2]], 1)
      ff[notOther] <- ff[notOther] + t(Umat%*%amat[,notOther])
    }
    ff[notOther] <-  ff[notOther] - epsilon
    sum( ff^2 )
  }
  
  lo <- rep(0, S)
  hi <- 2*max(ydata/effMat)
  
  if( verbose ) print(dim(xgrid))
  
  if(nrow(xgrid) > 10)pbar <- txtProgressBar(min = 1,max = nrow(xgrid), style=1)
  
  ccMu <- matrix(0, nrow(xgrid), S)
  colnames(ccMu) <- snames
  rownames(ccMu) <- rownames(xgrid)
  ccSd <- ccMu
  
  xnames <- rownames(bg)
  
  ii <- burnin:ng
  
  for(j in 1:nrow(xgrid)){
    
    xmu   <- xgrid[drop=F, j, xnames]
    xlmu  <- xgrid[drop=F, j, xlnames]
    wstar <- matrix(NA, nsim, S)
    
    for(k in 1:nsim){
      
      ig <- sample(ii, 1)
      
      if(termB){
        bmat <- bg*0
        bmat[ 1:length(bg) ] <- bgibbs[ig,]
        
        cnn <- colnames(xmu)
        if(length(intb) > 0){
          for(m in 1:length(intb)){
            cii <- unlist( strsplit(intb[m], ':') )
            xmu <- cbind(xmu, xmu[1,cii[1]]*xmu[1,cii[2]])
          }
          colnames(xmu) <- c(cnn,intb)
        }
      }
      
      if(termA){
        amat <- Amat*0
        amat[ wA ] <- alphaGibbs[ig,]
      }
      
      if(termR){
        rmat <- Rmat*0
        rmat[ wL ] <- lgibbs[ig,]
        cnn <- colnames(xlmu)
        if(length(intr) > 0){
          for(m in 1:length(intr)){
            cii <- unlist( strsplit(intr[m], ':') )
            xlmu <- cbind(xlmu, xlmu[1,cii[1]]*xlmu[1,cii[2]])
          }
          colnames(xlmu) <- c(cnn,intr)
        }
      }
      
      epsilon <- .rMVN(1, 0, sigma = sigMu )[notOther]
      
      wstar[k,] <- optim(wstart, fn, method = "L-BFGS-B", 
                         termA = termA, termB = termB, termR = termR, 
                         xmu = xmu, xlmu = xlmu,
                         gindex = gindex, uindex = uindex, 
                         notOther = notOther,
                         amat = amat, bmat = bmat, rmat = rmat, 
                         epsilon = epsilon, lower = lo, upper  = hi )$par
    }
    ccMu[j,] <- colMeans(wstar)
    ccSd[j,] <- apply(wstar, 2, sd)
    
    if(nrow(xgrid) > 10)setTxtProgressBar(pbar,j)
    
  }
  ccMu[,other] <- ccSd[,other] <- NA
  
  list(x = xgrid, ccMu = ccMu, ccSd = ccSd)
}

plotEquilAbund <- function(output, nsim = 20, ngrid = NULL, BYFACTOR = FALSE, 
                           verbose = T){
  
  wstar <- .wrapperEquilAbund(output, nsim, ngrid, BYFACTOR, 
                              verbose = verbose)
  ccMu <- wstar$ccMu[,notOther]
  ccSd <- wstar$ccSd[,notOther]
  ccx  <- wstar$x
  ccx  <- ccx[, !colnames(ccx) %in% attributes(ccx)$factors, drop=F]
  ccx  <- ccx[,-1, drop=F]
  
  np <- ncol(ccMu)
  
  npage <- 1
  o   <- 1:np
  if(np > 16){
    npage <- ceiling(np/16)
    np    <- 16
  }
  
  mfrow <- .getPlotLayout(np)
  
  nbin <- 12
  ylimit <- c(0, max(ccMu) )
  
  for(m in 1:ncol(ccx)){   # loop over predictors
    
    xm  <- colnames(ccx)[m]
    xx  <- ccx[,m]
    atx <- quantile(xx,seq(0, 1, length=nbin))
    
    xlimit <- c( sum( c(.7, .3)*atx[1:2] ), sum( c(.3, .7)*atx[(nbin-1):nbin] ) )
    
    k   <- 0
    add <- F
    o   <- 1:np
    o   <- o[o <= 16]
    
    for(p in 1:npage){
      
      file <- paste('equilAbund_', xm, '_', p,'.pdf',sep='')
      
      if(SAVEPLOTS)pdf( file=.outFile(outFolder,file) )
      
      npp <- ncol(ccMu) - k
      if(npp > np)npp <- np
      mfrow <- .getPlotLayout(np)
      par(mfrow=mfrow, bty='n', omi=c(.3,.3,0,0), mar=c(3,3,2,1), 
          tcl= tcl, mgp=mgp)
      
      for(j in o){
        
        yy  <- ccMu[,j]
        
        minbin <- 5
        xmids  <- 0
        while( length(xmids) == 1){
          minbin <- minbin - 1
          tt  <- .getBin(xx, yy, minbin = minbin)
          xbin <- tt$xbin
          xmids <- tt$xmids
        }
        
        c95 <- tapply( yy, list(bin = xbin), quantile, pnorm( c(-1.96, -1, 1, 1.96) ))
        ci <- matrix( unlist( c95 ), ncol = 4, byrow = T )
        rownames(ci) <- names(c95)
        
        xmids <- xmids[ as.numeric(names(c95)) ]
        
        .shadeInterval( xmids, loHi = ci[,c(1, 4)],col=specColor[j],PLOT=T,add=F,
                        xlab=' ',ylab=' ', xlim = xlimit,  
                        LOG=F, trans = .3)
        .shadeInterval( xmids, loHi = ci[,c(2, 3)],col=specColor[j], PLOT=T, add=T, 
                        trans = .3)
        mu <- tapply( yy, xbin, mean )
        lines(xmids,  mu, lty=2, lwd=2, col = specColor[j])
        
        
        k <- k + 1
        if(k > 26)k <- 1
        
        lab <- colnames(ccMu)[j]
        
        .plotLabel( lab,above=T )
      }
      mtext(xm, 1, outer=T)
      mtext('Equilibrium abundance', 2, outer=T)
      
      if(!SAVEPLOTS){
        readline('equilibrium abundance -- return to continue ')
      } else {
        dev.off()
      }
      o <- o + 16
      o <- o[o <= SO]
    }
  }
}

.rMVN <- function (nn, mu, sigma = NULL, sinv = NULL){
  
  # nn - no. samples from one mu vector or nrow(mu) for matrix
  
  if(!is.null(sigma)){
    m <- ncol(sigma)
  }else if(!is.null(sinv)){
    m <- ncol(sinv)
  }else{
    stop( '.rMNV requires either sigma or sinv' )
  }
  
  if(length(mu) > 1){
    if( !is.matrix(mu) ) mu <- matrix( mu, nn, length(mu) )  # mu is a vector of length m
    if( ncol(mu) == 1 & nn == 1 )  mu <- t(mu)
    if( length(mu) == m & nn > 1) mu <- matrix( mu, nn, length(mu), byrow=T )
  }
  
  if(is.null(sinv)){          # from sigma
    
    vv <- try(svd(sigma),T)
    
    if( inherits(vv,'try-error') ){
      ev <- eigen(sigma, symmetric = TRUE)
      rr <- t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    } else {
      rr <- vv$v %*% (t(vv$u) * sqrt(vv$d)) 
    }
    
  }else{     # from sinv
    
    L  <- chol(sinv)  
    rr <- backsolve(t(L), diag(m), upper.tri = F) 
  }
  ps <- matrix(rnorm(nn * m), nn) %*% rr
  ps + mu 
}

gjamSimTime <- function(S, Q = 0, nsite, ntime = 50, termB, termR, termA, obsEffort = 100,
                        predPrey = NULL, zeroAlpha = NULL, PLOT = FALSE){
  
  # observations are counts ('DA' in gjam)
  # S           - no. species in Y
  # Q           - no. predictors in X
  # nsite       - no. groups/plots/sites
  # ntime       - mean length of time series
  # termB, termR, termA are logical for inclusion of immigration/emigration X%*%B, 
  #               DI growth VL, DD growth UA
  # obsEffort   - effort used in gjam
  # predPrey    - interactions assumed to be competition, unless predPrey = c(pred, prey), where
  #               prey and pred are the rows/columns in alpha 
  # zeroAlpha   - second column does not affect first column
  # PLOT        - logical to plot series
  
  cols <- colF( S )
  
  wdata <- NULL
  if(Q < 1)Q <- 1
  
  if(Q <= 1 | !termB){
    form <- '~ 1'
  }else{
    form <- paste0( 'x', c(2:Q), collapse = '+')
    form <- as.formula( paste(' ~', form ) )
  }
  
  if(!is.null(predPrey) & length(predPrey) == 2)   predPrey <- matrix(predPrey, nrow=1 )
  if(!is.null(zeroAlpha) & length(zeroAlpha) == 2)zeroAlpha <- matrix(zeroAlpha, nrow=1 )
  
  beta <- rhoTrue <- alphaTrue <- wstar <- NULL
  
  nt     <- 2 + rpois(nsite, ntime - 2)
  ntot   <- sum(nt)
  groups <- rep(1:nsite, nt)
  
  times  <- 1:nt[1]
  if(nsite > 1){
    for(k in 2:nsite)times <- c(times, c(1:nt[k]))
  }
  
  w <- matrix(1000,ntot,S)
  colnames(w) <- paste('s', 1:S, sep='')
  
  if(termB){ # environmental immigration/emigration
    
    bsd <- 2
    if(termR)bsd <- 2
    if(termA)bsd <- 10
    
    x     <- matrix( 0, ntot, Q)
    beta  <- matrix( rnorm(Q*S, 0, bsd), Q, S )
    beta[1,] <- rnorm(S, 0, 1)
    colnames(x) <- rownames(beta) <- paste('x',1:Q,sep='')
    rownames(beta)[1] <- 'intercept'
    colnames(beta) <- paste('s', 1:S, sep='')
    
  }else{
    
    x <- matrix(1, ntot, 1)
    colnames(x) <- 'intercept'
    
  }
  
  if(termR){ # growth rate rho
    
    gam <- runif(S, -.03, .03)       
    if(termB){
      gam <- runif(S, -.03, .05)
    }
    if(termA){
      gam <- runif(S, .01, .1)
    }
    rhoTrue <- gam
  }
  
  if(termA){ # alpha matrix
    
    intrWt <- 1.2
    if(S > 10)intrWt <- 10
    
    allPos <- F
    np     <- 0
    
    while(!allPos){
      daa <- runif(S,-1/80000,-1/150000)  # competition
      ir  <- runif(S^2, min(daa), 0)
      aa  <- matrix(ir, S, S)
      diag(aa) <- intrWt*daa
      
      if(!is.null(predPrey))aa[ predPrey ] <- -aa[ predPrey ] # prey has positive effect on pred
      if(!is.null(zeroAlpha))aa[ zeroAlpha ] <- 0             # no effect
      alphaTrue <- aa
      
      wstar <- -solve(crossprod(aa))%*%t(aa)%*%gam # carrying capacity
      if( all(wstar > 0) )allPos <- T
      np <- np + 1
    }
    print( paste(np, 'iterations for termA') )
  }
  
  # residual covariance
  sigma <- diag(1, S)
  XB    <- 0
  
  for(k in 1:nsite){
    
    ww     <- matrix(10, nt[k], S)
    ww[1,] <- runif(S, 200, 400)
    if(!termR & !termA)ww[1,] <- beta[1,]
    if(termR & termA)  ww[1,] <- runif(S, wstar*.1, wstar*2)
    
    xx <- matrix( rnorm(nt[k]*Q), nt[k], Q)
    xx[,1] <- 1
    if(termB)XB <- xx%*%beta
    
    for(t in 2:nt[k]){
      
      ep    <- .rMVN(1, 0, sigma)
      ww[t,] <- ww[t-1,] + t( ep )
      if(termB) ww[t,] <- ww[t,] + XB[t,] 
      if(termR) ww[t,] <- ww[t,] + ww[t-1,]*gam
      if(termA) ww[t,] <- ww[t,] + diag(ww[t-1,])%*%aa%*%t(ww[drop=F,t-1,]) 
      
      if(sum(ww[t,]) > 1000000)stop('try again')
    }
    wk <- which(groups == k)
    w[wk,] <- ww
    if(termB)x[wk,] <- xx
  }
  
  wkeep <- which( is.finite(rowSums(w)) & apply(w, 1, min) > 0 )
  w <- w[wkeep,]
  if(termB){
    x <- x[wkeep,]
    x <- x[,!colnames(x) == 'x']
  }
  times <- times[wkeep]
  groups <- groups[wkeep]
  
  y <- round( w*obsEffort )
  colnames(y) <- colnames(w)
  
  if(PLOT){
    ylim <- c(0, 1.5*max(y))
    xlim <- c(0, max(nt))
    par(bty='n', cex=1.5)
    xlab <- expression( paste("Time ", italic(t), sep=''))
    ylab <- expression( paste("Count ", italic(y[s][t]), sep=''))
    
    wk <- which(groups == 1)
    plot(times[wk], y[wk,1],type='l', xlim = xlim, ylim=ylim, xlab=xlab, ylab=ylab)
    
    for(k in 1:nsite){
      wk <- which(groups == k)
      for(s in 1:S){
        lines(times[wk], y[wk,s],col = cols[s],lwd=2)
        if(termA)abline(h=wstar[s]*obsEffort, lty=2, col= cols[s], lwd=2)
      }
    }
  }
  rownames(w) <- paste(groups, times, sep='-')
  
  if(termA){
    print( 'eigenvalues and carrying capacities' )
    wdata <- data.frame( eigenvalues = eigen(aa)$values, carryingCapacity = wstar )
    print( wdata )
  }
  
  rho <- matrix(0, S, S)
  diag(rho) <- rhoTrue
  
  trueValues <- list(beta = beta, rho = rho, alpha = alphaTrue, 
                     sigma = sigma, w = w)
  
  wt <- which( sapply(trueValues, length) > 0 )
  trueValues <- trueValues[ wt ]
  
  xdata <- data.frame( groups = groups, times = times, x, stringsAsFactors = F)
  
  list(xdata = xdata, ydata = y, edata = y*0 + obsEffort, formula = form, 
       groups = groups, times = times, trueValues = trueValues,
       wdata = wdata)
}

.cleanNames <- function(xx){
  
  xx <- .replaceString(xx,'-','')
  xx <- .replaceString(xx,'_','')
  xx <- .replaceString(xx,' ','')
  xx <- .replaceString(xx,"'",'')
  xx
}


.replaceString <- function(xx,now='_',new=' '){  #replace now string in vector with new
  
  ww <- grep(now,xx,fixed=T)
  if(length(ww) == 0)return(xx)
  
  for(k in ww){
    s  <- unlist( strsplit(xx[k],now,fixed=T) )
    ss <- s[1]
    if(length(s) == 1)ss <- paste( ss,new,sep='')
    if(length(s) > 1)for(kk in 2:length(s)) ss <- paste( ss,s[kk],sep=new)
    xx[k] <- ss
  }
  xx
}


getDesign <- function(form, xdata){
  
  xm <- xs <- NULL
  
  tmp <- model.frame(form, data=xdata, na.action=NULL) # standardized design
  x   <- model.matrix(form, data=tmp)
  
  # do not center intercept or factors
  wf <- names( which( attr(attributes(tmp)$terms, 'dataClasses') == 'factor' ) )
  wf <- which( startsWith(colnames(x), wf ) )
  wf <- c(1, wf)
  wp <- c(1:ncol(x))[-wf]
  
  if(length(wp) > 0){
    xm  <- colMeans(x, na.rm=T)
    xs  <- apply(x, 2, sd, na.rm=T)
    
    w0 <- which(xs[-1] == 0)
    if(length(w0) > 0)stop( paste('no variation in predictor', names(xs)[w0-1]) )
    x[,wp]   <- t( (t(x[,wp]) - xm[wp])/xs[wp]  )
  }
  
  list(x = x, xmean = xm, xsd = xs)
  
}

gjamTimePrior <- function( xdata, ydata, edata, priorList, minSign = 5, 
                           betaMax = 10, rhoMax = 1 ){
  
  # minSign - minimum number of co-occurrences to estimate alpha
  
  bp <- lp <- ap <- NULL
  formulaBeta <- formulaRho <- alphaSign <- NULL
  termB <- termR <- termA <- FALSE
  
  betaPrior <- rhoPrior <- alphaPrior <- NULL
  
  for(k in 1:length(priorList))assign( names(priorList)[k], priorList[[k]] )
  
  S <- ncol(ydata)
  w <- as.matrix(ydata)/as.matrix(edata)
  n <- nrow(w)
  ynames <- colnames(ydata)
  
  colnames(w) <- .cleanNames( colnames(w ) )
  
  maxw <- apply(w, 2, max)
  quaw <- apply(w,2, quantile, c(.1, .9), na.rm=T)
  ranw <- apply(quaw, 2, diff)
  ranw[ ranw == 0 ] <- maxw[ ranw == 0 ]/5
  
  wmin <- quaw[1,]
  wmin[ wmin == 0] <- ranw[wmin == 0]/20
  
  if(!is.null(betaPrior))termB <- TRUE
  if(!is.null(rhoPrior)) termR <- TRUE
  if(!is.null(alphaSign))termA <- TRUE
  
  if( 'formulaBeta' %in% names(priorList) )termB <- TRUE
  
  # variance in dy
  gr  <- unique(xdata$groups)
  nr  <- length(gr)
  
  timeZero <- grep('-0', rownames(ydata)) # rownames from gjamFillMissing
  timeLast <- c( timeZero - 1, nrow(ydata))[-1]
  
  if(termB){  # only XB
    
    x <- getDesign(formulaBeta, xdata)$x
    
    kname <- as.vector( outer( colnames(x), colnames(w), FUN = paste, sep='-') )
    
    xnames <- attributes(x)$dimnames[[2]]
    blo <- matrix(Inf, ncol(x), ncol(ydata))
    rownames(blo) <- colnames(x)
    colnames(blo) <- colnames(ydata)
    bhi <- -blo
    
    sumw <- sumw2 <- rep(0, ncol(w))
    sumx <- sumx2 <- rep(0, ncol(x))
    
    xxj <- ddw <- numeric(0)
    
    for(j in gr){
      wj <- which(xdata[,'groups'] == j & is.finite(rowSums(x)) &
                    is.finite(rowSums(w)))
      wj <- wj[ !wj %in% timeZero ]
      if(length(wj) < 2)next
      
      dw <- w[ drop = F, wj[ -1 ],] - w[ drop = F, wj[ -length(wj) ], ]
      rw <- matrix( ranw, nrow(dw), ncol(dw), byrow=T)
      dw[dw > rw]  <- rw[dw > rw]
      dw[dw < -rw] <- rw[dw < -rw]
      
      xj <- x[drop=F, wj,][-1,, drop=F]
      
      xxj <- rbind(xxj, xj)
      ddw <- rbind(ddw, dw)
    }
    
    IXX <- solve( crossprod(xxj) )
    WX  <- crossprod(xxj,ddw)
    bb <- IXX%*%WX
    
    sig <- crossprod( ddw - xxj%*%bb )
    bvr  <- kronecker(sig,IXX)
    
    bsd  <- matrix( sqrt(diag(bvr)), nrow(bb), ncol(bb) )
    
    bl <- bb - 10*bsd
    bh <- bb + 10*bsd
    
    rownames(bl)[1] <- rownames(bh)[1] <- 'intercept'
    
    blo[ blo < -betaMax ] <- -betaMax
    bhi[ bhi > betaMax ]  <- betaMax
    
    if( is.list(betaPrior) ){
      gg  <- gjamPriorTemplate(formulaBeta, xdata, ydata = ydata, 
                               lo = betaPrior$lo, hi = betaPrior$hi)
      blo <- gg[[1]]
      bhi <- gg[[2]]
      
      blo[ !rownames(blo)%in% names(betaPrior$lo), colnames(blo) ] <- 
        bl[ !rownames(blo)%in% names(betaPrior$lo),  ]
      
      bhi[ !rownames(bhi)%in% names(betaPrior$hi), colnames(bhi) ] <- 
        bh[ !rownames(bhi)%in% names(betaPrior$hi),  ]
      
    }else{
      blo <- bl
      attr(blo,'formula') <- formulaBeta
      bhi <- bh
    }
    
    bp  <- list(lo = blo, hi = bhi )
  }
  
  if(termR){  # rho
    
    x   <- getDesign(formulaRho, xdata)$x
    colnames(x)[1] <- 'intercept'
    
    mlo <- matrix(-10, ncol(x), S)
    colnames(mlo) <- ynames
    rownames(mlo) <- colnames(x)
    mhi <- -mlo
    
    dw <- w[-timeZero,] - w[-timeLast,]
    xj <- x[-timeLast,]
    XX <- crossprod(xj)
    XI <- try( solve(XX), T)
    if( !inherits(XI,'try-error') ){
      bb <- XI%*%crossprod(xj, dw)
      rownames(bb) <- colnames(x)
      mlo <- -2*abs(bb)
      mhi <- 2*abs(bb)
    }
    
    for(k in 1:length(rhoPrior$lo)){
      rlo <- rhoPrior$lo[[k]]
      if(length(rlo) == 1)rlo <- rep(rlo, S)
      mlo[ names(rhoPrior$lo)[k], ] <- rlo
    }
    for(k in 1:length(rhoPrior$hi)){
      rhi <- rhoPrior$hi[[k]]
      if(length(rhi) == 1)rhi <- rep(rhi, S)
      mhi[ names(rhoPrior$hi)[k], ] <- rhi
    }
    
    sumb <- (mlo + mhi)/2
    sumn <- mlo*0 + 1
    
    for(j in gr){
      
      wj <- which(xdata[,'groups'] == j & is.finite(rowSums(x)) &
                    is.finite(rowSums(w)))
      if(length(wj) < 2)next
      wj <- wj[ !wj %in% timeZero ]
      w0 <- wj - 1
      
      yj <- w[ wj, ] + matrix(wmin,length(wj),S, byrow=T)
      
      dw <- w[drop=F, wj,] - w[drop=F, w0, ]
      rw <- matrix( ranw, nrow(dw), ncol(dw), byrow=T)
      dw[dw > rw]  <- rw[dw > rw]
      dw[dw < -rw] <- rw[dw < -rw]
      dw  <- dw/yj                   # per capita
      SS  <- F
      
      if( length(wj) > ncol(x) ){     
        
        xj <- x[drop=F,w0,]
        
        wm <- which( is.na(xj), arr.ind=T )
        if(length(wm) > 0){
          ic <- unique(wm[,2])
          for( m in 1:nrow(wm) ){
            xm  <- xj[, wm[m,2] ]
            tt  <- 1:length(xm)
            
            xsm <- lm(xm ~ tt, na.action = na.omit)
            xmiss <- predict( xsm, newdata = data.frame(tt = tt) )
            
            xm[ is.na(xm) ] <- xmiss[ is.na(xm) ]
            xj[, wm[m,2] ] <- xm
          }
          x[w0,] <- xj
        }
        
        XX <- crossprod(xj)
        XI <- try( solve(XX), T)
        if( inherits(XI,'try-error') )next
        bb <- XI%*%crossprod(xj, dw)
        
        ww <- which(bb[1,] > mhi[1,])   # intercept outside prior
        if(length(ww) > 0){
          bb[1,ww] <- mhi[1,ww]
          SS <- T
        }
        ww <- which(bb[1,] < mlo[1,])
        if(length(ww) > 0){
          bb[1,ww] <- mlo[1,ww]
          SS <- T
        }
        if(SS & ncol(x) > 1){
          dx <- dw - matrix(bb[1,],nrow(dw),S)
          bx <- solve( crossprod(xj[,-1]) )%*%crossprod(xj[,-1], dx)
          bb[2:nrow(bb),] <- bx
        }
        
        sumb[ bb != 0 ] <- sumb[ bb != 0 ] + bb[ bb != 0 ]
        sumn[ bb != 0 ] <- sumn[ bb != 0 ] + 1
      }
    }
    
    if( is.na(range(mlo)[1]) | is.na(range(mhi)[1]) ){
      
      dw <- w[-timeZero,] - w[-timeLast,]
      xj <- x[-timeLast,]
      XX <- crossprod(xj)
      XI <- try( solve(XX), T)
      if( !inherits(XI,'try-error') ){
        bb <- XI%*%crossprod(xj, dw)
        rownames(bb) <- colnames(x)
        mlo <- -2*abs(bb)
        mhi <- 2*abs(bb)
      }
    }
    
    rm <- sumb/sumn
    rm[ !is.finite(rm) ] <- 0
    #  rl <- -1.5*abs(rm)
    rl <- mlo
    rh <- -rl
    rl[ rl[1,] > mlo[1,] ] <- mlo[rl[1,] > mlo[1,] ]
    rh[ rh[1,] < mhi[1,] ] <- mhi[rh[1,] < mhi[1,] ]
    
    lp <- list(lo = rl, hi = rh)
    
  }
  
  if(termA){  # alpha
    
    if( !is.list(rhoPrior) ){
      if(is.null(lp))stop(' must have rhoPrior if there is alphaSign' )
      ap <- NULL
    }else{
      
      #must co-occur
      yind <- ydata
      rownames(alphaSign) <- colnames(alphaSign) <- 
        colnames(yind) <- .cleanNames( colnames(alphaSign ) )
      
      yind[yind > 1] <- 1
      yind <- crossprod(yind)
      
      yind[yind < minSign] <- 0  # minimum co-occurence
      yind[yind > 1] <- 1
      
      alphaSign <- alphaSign*yind[rownames(alphaSign),colnames(alphaSign)]
      
      rho <- (lp$lo['intercept',] + lp$hi['intercept',])/2 + .01
      
      timeZero <- which(xdata$times == 0)
      timeLast <- (timeZero - 1)[-1]
      wt <- unique( c(timeZero, timeLast) )
      
      quanw <- apply(w, 2, quantile, .9, na.rm=T)
      maxw  <- apply(w, 2, max, na.rm=T)/2
      maxw[quanw > 0] <- quanw[ quanw > 0]
      minw  <- maxw/20
      
      wmu    <- colMeans( w, na.rm=T )
      wdelta <- apply( w, 2, diff )  # growth rates
      wdelta[!is.finite(wdelta) ] <- 0
      wdelta <- wdelta[-wt,]         # exclude first and last
      
      wplus <- w + matrix(minw, nrow(w), S, byrow=T)
      
      ww <- n/crossprod(wplus)
      wrange <- apply(wdelta, 2, quantile, c(.05,.95), na.rm=T)
      
      w0 <- which( wrange[1,] == wrange[2,] )
      if(length(w0) > 0){
        wrange[1,w0] <- wrange[1,w0] - 1
        wrange[2,w0] <- wrange[2,w0] + 1
      }
      
      wrange[1, wrange[1,] >= 0 ] <- mean( wrange[1, wrange[1,] < 0 ] )
      wrange[2, wrange[2,] <= 0 ] <- mean( wrange[2, wrange[2,] > 0 ] )
      
      wlo <- matrix(wrange[1,], S, S, byrow=T)*ww   # E[dw]/E[w_s * w_s']
      whi <- matrix(wrange[2,], S, S, byrow=T)*ww   # E[dw]/E[w_s * w_s']
      
      rlh <- matrix( rho/maxw, S, S, byrow= T) # rho/w_s
      
      whi[ whi > rlh ]  <- rlh[ whi > rlh ] 
      wlo[ wlo < -rlh ] <- -rlh[ wlo < -rlh ] 
      
      alo <- ahi <- wlo*0
      ww  <- which(alphaSign < 0)
      
      scale <- 100
      if(max(whi) > .1)scale <- 1
      if(length(ww) > 0){
        alo[ww] <- scale*wlo[ww]
        ahi[ww] <- 0
      }
      ww <- which(alphaSign > 0)
      if(length(ww) > 0){
        alo[ww] <- 0
        ahi[ww] <- scale*whi[ww]
      }
      ww <- which(alphaSign == 0)
      if(length(ww) > 0)alo[ww] <- ahi[ww] <- NA
      
      ap <- list(lo = alo, hi = ahi)
    }
  }
  
  list(betaPrior = bp, rhoPrior = lp, alphaPrior = ap, formulaBeta = formulaBeta,
       formulaRho = formulaRho)
}

foodWebDiagram <- function(S, guildList = NULL, predPrey = NULL, zeroAlpha = NULL,
                           intraComp = 1:S, label = NULL, PLOT = TRUE, layout = 'rr'){
  
  # S - no. species
  # default interaction is negative, arrows only for negative interactions
  # guildList - overrides default, only members of the same guild compete
  # predPrey  - matrix with 2 columns, second column is prey of first column
  # zeroAlpha - matrix with 2 columns, second column does not affect first column
  # intraComp - needed for intraspecific comp if guildList is specified
  # layout can be 'tree', 'rr', ...
  
  require( DiagrammeR )
  
  pp <- numeric(0)
  
  fromTo <- as.matrix( expand.grid( c(1:S), c(1:S) ) )
  ft <- .pasteCols( fromTo )
  
  qq <- numeric(0)
  if( !is.null(guildList) ){
    
    fromTo <- numeric(0)
    for(k in 1:length(guildList)){
      if( is.character(guildList[[k]]) )guildList[[k]] <- match(guildList[[k]],snames)
      ft <- as.matrix(expand.grid(guildList[[k]], guildList[[k]]))
      fromTo <- rbind(fromTo, ft)
    }
    if(!is.null(intraComp)){
      ft <- cbind(1:S, 1:S)
      fromTo <- rbind(fromTo, ft)
    }
    ft <- .pasteCols( fromTo )
    ww <- which(!duplicated(ft))
    ft <- ft[ww]
    fromTo <- fromTo[ww,]
    #   zeroAlpha <- NULL
  }
  
  if(!is.null(predPrey)){
    pp <- .pasteCols( predPrey[drop=F,,c(2,1)] )
    qq <- .pasteCols( predPrey)
    #   fromTo <- fromTo[ !ft %in% pp, ]
    fromTo <- rbind(fromTo, predPrey, predPrey[drop=F,,c(2,1)])
    ft <- .pasteCols( fromTo )
    ww <- which(!duplicated(ft))
    fromTo <- fromTo[ww,]
    ft <- ft[ww]
  }
  
  if(!is.null(zeroAlpha)){
    za <- .pasteCols( zeroAlpha[drop=F,,c(1,2)])
    
    fromTo <- fromTo[ !ft %in% za, ]
    ft <- ft[ !ft %in% za ]
  }
  
  if( length(qq) > 0 ){
    ft <- .pasteCols( fromTo )
    fromTo <- rbind( fromTo[ft %in% qq,], fromTo[!ft %in% qq, ] )
  }
  ft <- .pasteCols( fromTo )
  pp <- pp[ pp %in% ft ]
  
  if(!PLOT){
    return( fromTo )
  }else{
    ecol  <- rep( 'tan', length(ft) )
    #  ncol  <- rep( 'tan', S )
    ncol  <- colF(S)
    shape <- rep( 'rectangle', S)
    
    if( length(qq) > 0 ){
      
      wt <- which( ft %in% qq )
      ecol[ wt ] <- 'brown'
      
    }
    if( length(pp) > 0 ){
      
      wt <- which( ft %in% pp )
      ecol[ wt ] <- 'blue'
    }
    
    if( is.null(layout) )layout <- "tree"
    if(is.null(label))label <- paste('s', 1:S, sep='')
    nodes <- create_node_df(n = S, label = label, style = "filled", color = ncol, shape = shape)
    edges <- create_edge_df(from = fromTo[,1], to = fromTo[,2], color = ecol )
    graph <- create_graph(nodes_df = nodes, edges_df = edges)
    render_graph( graph,  layout = layout )
  }
}

lowerFirstLetter <- function(xx){
  s <- unlist(strsplit(xx, " "))
  s <- paste(tolower(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

upperFirstLetter <- function(xx, FIRSTONLY = F){
  
  # FIRSTONLY - only first letter of first word when a string has multiple words
  
  if(FIRSTONLY){
    return( paste(toupper(substring(xx, 1, 1)), substring(xx, 2),
                  sep = "") )
  }
  s <- unlist(strsplit(xx, " "))
  s <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
             sep = "", collapse = " ")
  unlist(strsplit(s, " "))
}

columnSplit <- function(vec, sep='_', ASFACTOR = F, ASNUMERIC = FALSE,
                        LASTONLY = FALSE){
  
  vec <- as.character(vec)
  nc  <- length( strsplit(vec[1], sep, fixed  = TRUE)[[1]] )
  
  mat <- matrix( unlist( strsplit(vec, sep, fixed  = TRUE) ), ncol=nc, byrow  = TRUE )
  if(LASTONLY & ncol(mat) > 2){
    rnn <- mat[,1]
    for(k in 2:(ncol(mat)-1)){
      rnn <- columnPaste(rnn,mat[,k])
    }
    mat <- cbind(rnn,mat[,ncol(mat)])
  }
  if(ASNUMERIC){
    mat <- matrix( as.numeric(mat), ncol=nc )
  }
  if(ASFACTOR){
    mat <- data.frame(mat)
  }
  if(LASTONLY)mat <- mat[,2]
  mat
}

columnPaste <- function(c1, c2, sep='-', NOSPACE = FALSE){
  
  c1    <- as.character(c1)
  c2    <- as.character(c2)
  if(NOSPACE){
    c1   <- .replaceString(c1, ' ', '')
    c2   <- .replaceString(c2, ' ', '')
  }
  c12   <- apply( cbind(c1, c2) , 1, paste0, collapse=sep)
  
  c12
}


.plotLabel <- function(label,location='topleft',cex=1.3,font=1,
                       above=F,below=F,bg=NULL){
  
  if(above){
    adj <- 0
    if(location == 'topright')adj=1
    title(label,adj=adj, font.main = font, font.lab =font)
    return()
  }
  if(below){
    adj <- 0
    if(location == 'bottomright')adj=1
    mtext(label,side=1,adj=adj, outer=F,font.main = font, font.lab =font,cex=cex)
    return()
  }
  
  if(is.null(bg)){
    tmp <- legend(location,legend=' ',bty='n')
  } else {
    tmp <- legend(location,legend=label,bg=bg,border=bg,text.col=bg,bty='o')
  }
  
  xt <- tmp$rect$left # + tmp$rect$w
  yt <- tmp$text$y
  
  pos <- 4
  tmp <- grep('right',location)
  if(length(tmp) > 0)pos <- 2
  
  XX <- par()$xlog
  YY <- par()$ylog
  
  if(XX)xt <- 10^xt
  if(YY)yt <- 10^yt
  
  text(xt,yt,label,cex=cex,font=font,pos=pos)
}

.getColor <- function(col,trans){
  
  # trans - transparency fraction [0, 1]
  
  tmp <- col2rgb(col)
  rgb(tmp[1,], tmp[2,], tmp[3,], maxColorValue = 255, 
      alpha = 255*trans, names = paste(col,trans,sep='_'))
}