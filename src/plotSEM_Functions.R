#Plot SEM for log2 normalized data

plotSEM <- function(basedir,TFname,error=FALSE, cexError=FALSE, colError=FALSE, reverse=FALSE, cex=1) {
  score <- read.delim (paste(basedir, "/", TFname, ".sem", sep=""), header = T)
  
  baseline <- read.delim(paste(basedir, "/BASELINE/baseline.maximums", sep=""), header=FALSE)
  baseline.sig.mean <- baseline[1,2]
  baseline.sig.count <- baseline[1,3]
  baseline.sig.std <- baseline[1,4]
  baseline.sig.stderr <- baseline[1,5]
  baseline.rnd.mean <- baseline[2,2]
  baseline.rnd.count <- baseline[2,3]
  baseline.rnd.std <- baseline[2,4]
  baseline.rnd.stderr <- baseline[2,5]
  
  if(error) {
    stderr <- read.delim (paste(basedir, "/", TFname, ".sterr", sep=""), header = T)
  }
  
  if(reverse) {
    temp <- score
    score$A <- rev(temp$T)
    score$T <- rev(temp$A)
    score$C <- rev(temp$G)
    score$G <- rev(temp$C)
    if(error) {
      temp <- stderr
      stderr$A <- rev(temp$T)
      stderr$T <- rev(temp$A)
      stderr$C <- rev(temp$G)
      stderr$G <- rev(temp$C)
    }
  }

  baseline.rnd.mean <- log2(baseline.rnd.mean/baseline.sig.mean)
  baseline.rnd.stderr <- baseline.rnd.stderr/baseline.sig.mean
  baseline.sig.stderr <- baseline.sig.stderr/baseline.sig.mean
  
  maxY<-max(score$A, score$C, score$G, score$T, 0)
  minY<-min(score$A, score$C, score$G, score$T, -1)
  
  plot (score[,1], score$A, ylim = c(minY, maxY), xlim = c(1, nrow(score)), xaxp = c(1, nrow(score), nrow(score)-1), type = "n", font = 2, main = paste("SNP Effect Matrix of ", TFname, sep=""), ylab = "SEM Score", xlab = paste("Location on ", TFname, " motif",sep=""))
  abline (h = 0, lwd = 4, col = "gray")
  abline (h = baseline.rnd.mean, lwd = 4, col = "gray", lty=2)
  
  #size scaling
  cexA = cex
  cexT = cex
  cexG = cex
  cexC = cex
  if(cexError) {
    cexA = semCexScale(stderr$A)
    cexT = semCexScale(stderr$T)
    cexC = semCexScale(stderr$C)
    cexG = semCexScale(stderr$G)
  }
  
  #colors
  colA = "green"
  colT = "red"
  colC = "blue"
  colG = "darkgoldenrod3"
  if(colError) {
    colA = semColScale(stderr$A, color=colA)
    colT = semColScale(stderr$T, color=colT)
    colC = semColScale(stderr$C, color=colC)
    colG = semColScale(stderr$G, color=colG)
  }
  
  text (score[,1], score$A, col = colA, labels = "A", font = 2, cex=cexA)
  text (score[,1], score$T, col = colT, labels = "T", font = 2, cex=cexT)
  text (score[,1], score$C, col = colC, labels = "C", font = 2, cex=cexC)
  text (score[,1], score$G, col = colG, labels = "G", font = 2, cex=cexG)
  
  if(error) {
    abline (h = baseline.rnd.mean-baseline.rnd.stderr, lwd = 1, col = "lightgray", lty = 3)
    abline (h = baseline.rnd.mean+baseline.rnd.stderr, lwd = 1, col = "lightgray", lty = 3)
    abline (h = 0-baseline.sig.stderr, lwd = 1, col = "lightgray", lty = 3)
    abline (h = 0+baseline.sig.stderr, lwd = 1, col = "lightgray", lty = 3)
    arrows(score[,1], score$A-stderr$A, score[,1], score$A+stderr$A, length = 0.05, angle = 90, code=3, col = "slategray")
    arrows(score[,1], score$T-stderr$T, score[,1], score$T+stderr$T, length = 0.05, angle = 90, code=3, col = "slategray")
    arrows(score[,1], score$C-stderr$C, score[,1], score$C+stderr$C, length = 0.05, angle = 90, code=3, col = "slategray")
    arrows(score[,1], score$G-stderr$G, score[,1], score$G+stderr$G, length = 0.05, angle = 90, code=3, col = "slategray")
  }
  
}

semCexScale <- function(counts, min=.3, max=4) {
  #linear scale
  scale <- ((counts/max(counts)) * (max-min)) + min
  return(scale)
}

semColScale <- function(counts, color, min=.3, max=1) {
  x<- col2rgb(color)/255
  #linear scale
  alphascale <- ((counts/max(counts)) * (max-min)) + min
  for(i in 1:length(alphascale)) {
    counts[i] <- rgb(x[1], x[2], x[3], alpha=alphascale[i])
  }
  return(counts)
}