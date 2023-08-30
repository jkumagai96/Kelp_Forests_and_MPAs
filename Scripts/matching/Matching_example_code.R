###################################################
### chunk number 1: setup
###################################################
library("Matching")
data("lalonde")
attach(lalonde)


###################################################
### chunk number 2: 
###################################################
Y  <- lalonde$re78   
Tr <- lalonde$treat  


###################################################
### chunk number 3: pscore1
###################################################
glm1  <- glm(Tr~age + educ + black +
             hisp + married + nodegr + re74  + re75,
             family=binomial, data=lalonde)


###################################################
### chunk number 4: match1
###################################################
rr1 <- Match(Y=Y, Tr=Tr, X=glm1$fitted)


###################################################
### chunk number 5: 
###################################################
    m1 = Match(Y=Y, Tr=Tr, X=glm1$fitted) 


###################################################
### chunk number 6: 
###################################################
    m1 = Match(Y=Y, Tr=Tr, X=glm1$fitted, estimand="ATT", M=1, ties=TRUE, 
               replace=TRUE) 


###################################################
### chunk number 7: matchbalance1 eval=FALSE
###################################################
## MatchBalance(Tr~age + I(age^2) + educ + I(educ^2) + black +
##              hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
##              u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
##              match.out=rr1, nboots=1000, data=lalonde)


###################################################
### chunk number 8: 
###################################################
MatchBalance(Tr~nodegr, match.out=rr1, nboots=1000, data=lalonde)


###################################################
### chunk number 9: 
###################################################
MatchBalance(Tr~re74, match.out=rr1, nboots=1000, data=lalonde)


###################################################
### chunk number 10: 
###################################################
  qqplot(lalonde$re74[rr1$index.control], lalonde$re74[rr1$index.treated])
  abline(coef=c(0,1), col=2)


###################################################
### chunk number 11: dw.pscore
###################################################
dw.pscore <- glm(Tr~age + I(age^2) + educ + I(educ^2) + black +
                 hisp + married + nodegr + re74  + I(re74^2) + re75 + 
                 I(re75^2) + u74 + u75, family=binomial, data=lalonde)
rr.dw <- Match(Y=Y, Tr=Tr, X=dw.pscore$fitted)


###################################################
### chunk number 12: dw.matchbalance eval=FALSE
###################################################
## MatchBalance(Tr~age + I(age^2) + educ + I(educ^2) + black +
##              hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
##              u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
##              data=lalonde, match.out=rr.dw, nboots=1000)  


###################################################
### chunk number 13: 
###################################################
MatchBalance(Tr~nodegr+re74+I(re74^2), match.out=rr.dw, nboots=1000, data=lalonde)


###################################################
### chunk number 14: genoud1
###################################################
X <- cbind(age, educ, black, hisp, married, nodegr, re74, re75, u74, u75)

BalanceMatrix <- cbind(age, I(age^2), educ, I(educ^2), black,
                       hisp, married, nodegr, re74 , I(re74^2), re75, 
                       I(re75^2), u74, u75, I(re74*re75), I(age*nodegr), 
                       I(educ*re74), I(educ*re75))

gen1 <- GenMatch(Tr=Tr, X=X, BalanceMatrix=BalanceMatrix, pop.size=1000)


###################################################
### chunk number 15: genoud2
###################################################
  mgen1 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen1)


###################################################
### chunk number 16: genoud2 eval=FALSE
###################################################
##   MatchBalance(Tr~age + I(age^2) + educ + I(educ^2) + black +
##              hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
##              u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
##              data=lalonde, match.out=mgen1, nboots=1000)  


###################################################
### chunk number 17: genmatchbalance
###################################################
MatchBalance(Tr~nodegr+re74+I(re74^2), match.out=mgen1, nboots=1000, data=lalonde)


###################################################
### chunk number 18: output
###################################################
  summary(mgen1)


###################################################
### chunk number 19:  eval=FALSE
###################################################
##   c("localhost","localhost","musil","musil","deckard")


###################################################
### chunk number 20: parallel1 eval=FALSE
###################################################
## library("snow")
## library("Matching")
## data("lalonde")
## attach(lalonde)
## 
## cl <- makeCluster(c("musil","quetelet","quetelet"), type="SOCK")
## 
## X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74)
## 
## genout <- GenMatch(Tr=treat, X=X, cluster=cl)
## 
## stopCluster(cl)


###################################################
### chunk number 21: parallel2 eval=FALSE
###################################################
## library("Matching")
## data("lalonde")
## attach(lalonde)
## 
##  X = cbind(age, educ, black, hisp, married, nodegr, u74, u75, re75, re74)
##  Xbig <- rbind(X, X, X, X)
## 
##  Ybig <- c(treat, treat, treat, treat)
## 
##  GenMatch(Tr=Ybig, X=Xbig, BalanceMatrix=Xbig, estimand="ATE", M=1,
##           pop.size=1000, max.generations=10, wait.generations=1,
##           int.seed=3818, unif.seed=3527,
##           cluster=c("localhost","localhost", 
##                     "localhost","localhost"))


###################################################
### chunk number 22: appendix.pscore1
###################################################
glm1  <- glm(Tr~age + educ + black +
             hisp + married + nodegr + re74  + re75,
             family=binomial, data=lalonde)
rr1 <- Match(Y=Y, Tr=Tr, X=glm1$fitted)
MatchBalance(Tr~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
             match.out=rr1, nboots=1000, data=lalonde)


###################################################
### chunk number 23: <appendix.dw
###################################################
dw.pscore <- glm(Tr~age + I(age^2) + educ + I(educ^2) + black +
                 hisp + married + nodegr + re74  + I(re74^2) + re75 + 
                 I(re75^2) + u74 + u75, family=binomial, data=lalonde)
rr.dw <- Match(Y=Y, Tr=Tr, X=dw.pscore$fitted)
MatchBalance(Tr~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
             data=lalonde, match.out=rr.dw, nboots=1000)  


###################################################
### chunk number 24: <appendix.genout.prelim eval=FALSE
###################################################
## X <- cbind(age, educ, black, hisp, married, nodegr, re74, re75, u74, u75)
## 
## BalanceMatrix <- cbind(age, I(age^2), educ, I(educ^2), black,
##                        hisp, married, nodegr, re74 , I(re74^2), re75, 
##                        I(re75^2), u74, u75, I(re74*re75), I(age*nodegr), 
##                        I(educ*re74), I(educ*re75))
## 
## gen1 <- GenMatch(Tr=Tr, X=X, BalanceMatrix=BalanceMatrix, pop.size=1000)


###################################################
### chunk number 25: <appendix.genmatch
###################################################
  mgen1 <- Match(Y=Y, Tr=Tr, X=X, Weight.matrix=gen1)

  MatchBalance(Tr~age + I(age^2) + educ + I(educ^2) + black +
             hisp + married + nodegr + re74  + I(re74^2) + re75 + I(re75^2) +
             u74 + u75 + I(re74*re75) + I(age*nodegr) + I(educ*re74) + I(educ*re75),
             data=lalonde, match.out=mgen1, nboots=1000)  


