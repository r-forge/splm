pkgname <- "splm"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('splm')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("bsjktest")
### * bsjktest

flush(stderr()); flush(stdout())

### Name: bsjktest
### Title: Baltagi, Song, Jung and Koh LM test for spatial panels
### Aliases: bsjktest.formula bsjktest
### Keywords: htest

### ** Examples

data(Produc, package="Ecdat")
Produc <- Produc[Produc$year<1975, ]
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
test1<-bsjktest(fm,data=Produc, w=usaww,
  test="C.1")
test1



cleanEx()
nameEx("bsktest")
### * bsktest

flush(stderr()); flush(stdout())

### Name: bsktest
### Title: Baltagi, Song and Koh LM test for spatial panels
### Aliases: bsktest bsktest.formula bsktest.lm bsktest.splm
### Keywords: htest

### ** Examples

data(Produc, package="Ecdat")
Produc <- Produc[Produc$year<1975, ]
data(usaww)
fm <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
test1<-bsktest(fm,data=Produc, w=mat2listw(usaww),
  test="SLM1")
class(test1)
test1
ml2 <- spfeml(fm, data = Produc, , mat2listw(usaww), model = "error", effects = "pooled")
class(ml2)
test5bis<-bsktest(ml2, w=mat2listw(usaww),index=Produc[,c(1,2)] ,test="CLMmu")
summary(test5bis)



cleanEx()
nameEx("effects.splm")
### * effects.splm

flush(stderr()); flush(stdout())

### Name: effects.splm
### Title: method for extracting fixed effects
### Aliases: effects.splm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE) 
summary(lag)
eff <- effects(lag) 
print(eff)
err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)
summary(err)
eff <- effects(err) 
print(eff)



cleanEx()
nameEx("listw2dgCMatrix")
### * listw2dgCMatrix

flush(stderr()); flush(stdout())

### Name: listw2dgCMatrix
### Title: Interface between Matrix class objects and weights list
### Aliases: listw2dgCMatrix
### Keywords: spatial

### ** Examples

library(spdep)
data(columbus)
listw<-nb2listw(col.gal.nb)
spW<-listw2dgCMatrix(listw)



cleanEx()
nameEx("print.effects.splm")
### * print.effects.splm

flush(stderr()); flush(stdout())

### Name: print.effects.splm
### Title: method for printing fixed effects from objects of class
###   effects.splm
### Aliases: print.effects.splm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE) 
summary(lag)
eff <- effects(lag) 
print(eff)
err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)
summary(err)
eff <- effects(err) 
print(eff)



cleanEx()
nameEx("print.splm")
### * print.splm

flush(stderr()); flush(stdout())

### Name: print.splm
### Title: print method for class splm
### Aliases: print.splm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
GM<-spregm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,w=usaww,method="fulweigh")
summary(GM)



cleanEx()
nameEx("spfegm")
### * spfegm

flush(stderr()); flush(stdout())

### Name: spfegm
### Title: GM estimator for spatial panel data models
### Aliases: spfegm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
GM<-spfegm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp-1, data=Produc,w=usaww,method="fulweigh", effects='spfe', lag=FALSE)
summary(GM)



cleanEx()
nameEx("spfeml")
### * spfeml

flush(stderr()); flush(stdout())

### Name: spfeml
### Title: Spatial Panel Fixed Effects Models Estimation
### Aliases: spfeml
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE)
summary(lag)
eff <- effects(lag)
err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)
summary(err)
eff <- effects(err)
print(eff)
write.effects.splm(eff)



cleanEx()
nameEx("splm-package")
### * splm-package

flush(stderr()); flush(stdout())

### Name: splm-package
### Title: Spatial panel models: estimation and testing
### Aliases: splm-package splm
### Keywords: package

### ** Examples

data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
GM<-spregm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,w=usaww,method="fulweigh")
summary(GM)
respaterr <- spreml(fm, data = Produc, w = usaww, errors="semre")
summary(respaterr)




cleanEx()
nameEx("spregm")
### * spregm

flush(stderr()); flush(stdout())

### Name: spregm
### Title: GM estimator for spatial panel data models
### Aliases: spregm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
GM<-spregm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,w=usaww,method="fulweigh")
summary(GM)



cleanEx()
nameEx("spreml")
### * spreml

flush(stderr()); flush(stdout())

### Name: spreml
### Title: Spatial Panel Random Effects Model Estimation
### Aliases: spreml
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1974, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
respaterr <- spreml(fm, data = Produc, w = usaww, errors="semre")
summary(respaterr)



cleanEx()
nameEx("spsegm")
### * spsegm

flush(stderr()); flush(stdout())

### Name: spsegm
### Title: Spatial Simultaneous Equations
### Aliases: spsegm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1973, ]
eq1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq2 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq3 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
formula<-list(tp1 = eq1, tp2 = eq2, tp3=eq3)
w<-mat2listw(usaww)
se<-spsegm(formula, data=Produc, w=w, panel= TRUE,lags=list(c(TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE),c(TRUE,TRUE,TRUE)), errors=list(FALSE,TRUE,FALSE),endogenous=list(c(FALSE,TRUE,FALSE),c(TRUE,FALSE,FALSE),c(TRUE,FALSE,FALSE)))
summary(se)



cleanEx()
nameEx("spseml")
### * spseml

flush(stderr()); flush(stdout())

### Name: spseml
### Title: Spatial SUR - Lag and Error
### Aliases: spseml similar.listw similar.listw_spam similar.listw_Matrix
###   can.be.simmed llsurerror llsurlag
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1973, ]
eq1 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq2 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
eq3 <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp 
formula<-list(tp1 = eq1, tp2 = eq2, tp3=eq3)
listw<-mat2listw(usaww)
sur.error<-spseml(formula, data = Produc, w = listw, model = "error", method = "eigen", quiet = TRUE)
summary(sur.error)
sur.lag<-spseml(formula, data = Produc, w = listw, model = "lag", method = "eigen", quiet = FALSE)
summary(sur.lag)



cleanEx()
nameEx("summary.splm")
### * summary.splm

flush(stderr()); flush(stdout())

### Name: summary.splm
### Title: summary method for class splm
### Aliases: summary.splm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat") 
data(usaww)
Produc <- Produc[Produc$year<1975, ] 
GM<-spregm(log(gsp)~log(pcap)+log(pc)+log(emp)+unemp, data=Produc,w=usaww,method="fulweigh")
summary(GM)



cleanEx()
nameEx("write.effects.splm")
### * write.effects.splm

flush(stderr()); flush(stdout())

### Name: write.effects.splm
### Title: method for writing a table with fixed effects
### Aliases: write.effects.splm
### Keywords: spatial

### ** Examples

data(Produc, package = "Ecdat")
data(usaww)
Produc <- Produc[Produc$year<1975, ]
fm <- log(gsp) ~ log(pcap) + log(pc) + log(emp) + unemp
lag <- spfeml(fm, data = Produc, listw = mat2listw(usaww), effects = "sptpfe", method = "eigen", quiet = TRUE)
summary(lag)
eff <- effects(lag) 
err <- spfeml(fm, data = Produc, listw = mat2listw(usaww), model = "error", effects = "tpfe", method = "eigen", quiet = TRUE)
summary(err)
eff <- effects(err)
summary(eff)
write.effects.splm(eff)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
