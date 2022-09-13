#*
#*


DCDKineticModel <- function(Mat1, Mat2, pTime1=NULL, pTime2=NULL, GT=NULL
                            , maineffect.models = FALSE, perBin=20){

  if(!require(AMCBGeneUtils))
    devtools::install_github("AndresMCB/AMCBGeneUtils")

# if not provided, we use Cosmic Cancer Gene Census (CGC) as conservative GT
  if(is.null(GT)){
    GT <- AMCBGeneUtils::CGC.driverNames
  }

  CD.baseline <- intersect(GT,putativeCD)



  # Ordering datasets following ptime ascending order

  ptime1 <- ptime1[order(ptime1),, drop=F]
  D.normal <- Mat1[row.names(ptime1),, drop = F]

  ptime2 <- ptime2[order(ptime2),, drop=F]
  D.tumour <- Mat2[row.names(ptime2),, drop = F]

  # binning Gene Expression
  bins <- floor(min(nrow(D.normal),nrow(D.tumour))/perBin)

  D.normal <- apply(D.normal, MARGIN = 2
                    , function(Feature,bins){
                      aux <- vector(length = bins)
                      delta <- round(length(Feature)/bins)
                      for (j in 1:(bins-1)) {
                        aux[j] <- mean(Feature[(delta*(j-1)+1):(delta*j)])
                      }
                      aux[bins] <- mean(Feature[(delta*j+1):length(Feature)])
                      return(aux)
                    }
                    ,bins)

  D.tumour <- apply(D.tumour, MARGIN = 2
                    , function(Feature,bins){
                      aux <- vector(length = bins)
                      delta <- round(length(Feature)/bins)
                      for (j in 1:(bins-1)) {
                        aux[j] <- mean(Feature[(delta*(j-1)+1):(delta*j)])
                      }
                      aux[bins] <- mean(Feature[(delta*j+1):length(Feature)])
                      return(aux)
                    }
                    ,bins)


  D <- rbind(c(D.normal),c(D.tumour))
  time <- 0:(bins-1)
  env <- 1:2
  ck.fit <- CausalKinetiX(D ,time, env, target = 1,
                          pars=list(expsize=1
                                    ,maineffect.models = maineffect.models
                                    ,average.reps=TRUE))

  ck.fit$variable.scores <- data.frame(variable=colnames(Mat1)[-1]
                                       ,ck.score = ck.fit$variable.scores[-1]
                                       ,ck.ranking = ck.fit$ranking[-1])

  Ranked <- ck.fit$variable.scores%>%
    arrange(ck.ranking)

  top <- matrix(ncol = 5, nrow = 2)
  colnames(top) <-  c("50","100","150","200","250")
  row.names(top) <- c("CDinGT","p.val")
  ntop <- c(50,100,150,200,250)

  for (i in 1:5) {
    if(!are_miR){
      top[1,i] <- length(intersect(GT,Ranked$variable[1:ntop[i]]))
    }else{
      top[1,i] <- length(intersect(GT
                                   ,get_miR (Ranked$variable[1:ntop[i]])))
    }
    top[2,i] <- 1 - phyper(q = top[1,i]-1
                           , m = length(CD.baseline)
                           , n = length(putativeCD)-length(CD.baseline)
                           , k = ntop[i], lower.tail = T, log.p = FALSE)
  }


  ck.fit$CD.baseline <- CD.baseline
  ck.fit$top <- top
  return(list(ck.fit))
}
