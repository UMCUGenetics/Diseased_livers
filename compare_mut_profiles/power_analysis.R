##
ttestEffSizes <- function(n1=2:8, n2=2:6){
   out <- expand.grid(n1=n1, n2=n2)
   out$eff_size_metric <- "cohens_d"
   out$eff_size <- sapply(1:nrow(out), function(i){
      n1 <- out$n1[i]
      n2 <- out$n2[i]
      pwr::pwr.t2n.test(n1=n1, n2=n2, sig.level=0.05, power=0.8, alternative='greater')$d
   })
   return(out)
}

ttestEffSizes()

##
chiSquareEffSizes <- function(N=2:10, df=2:10){
   out <- expand.grid(N=N, df=df)
   out$eff_size_metric <- "cohens_w"
   out$eff_size <- sapply(1:nrow(out), function(i){
      #i=1
      pwr::pwr.chisq.test(N=out$N[i], df=out$df[i], sig.level=0.05, power=0.8)$w
   })
   return(out)
}

chiSquareEffSizes(N=c(13,11), df=9)
chiSquareEffSizes(N=c(13,11), df=6)

