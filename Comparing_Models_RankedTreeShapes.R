install.packages("ape")
install.packages("apTreeshape")
install.packages("distory")
install.packages("MASS")
library("ape")
library("apTreeshape")
library("distory")
library("MASS")
set.seed(1)#To generate same results each time
numtrees=25#Change for number of trees of each model
numtips=5#Change depending on number of tips you want in each tree
alltreeslist7=rtreeshape(n=4*numtrees, tip.number=numtips, model="yule")
yuletreeslist7=rtreeshape(n=numtrees, tip.number=numtips, model="yule")
for (i in 1:numtrees) {
  yuletreeslist7[[i]] = ladderize(as.phylo(rtreeshape(n=1, tip.number=numtips, model="yule")[[1]]))
  yuletreeslist7[[i]]$tip.label = yuletreeslist7[[1]]$tip.label #Standardizing tip labels
  alltreeslist7[[i]]=yuletreeslist7[[i]]#Adding to full tree list
}
PDAtreeslist7=rtreeshape(n=numtrees, tip.number=numtips, model="pda")
for (i in 1:numtrees) {
  PDAtreeslist7[[i]] = ladderize(as.phylo(rtreeshape(n=1, tip.number=numtips, model="pda")[[1]]))
  PDAtreeslist7[[i]]$tip.label = yuletreeslist7[[1]]$tip.label 
  alltreeslist7[[numtrees+i]]=PDAtreeslist7[[i]]
}
Aldoustreeslist7=rtreeshape(n=numtrees, tip.number=numtips, model="aldous")
for (i in 1:numtrees) {
  Aldoustreeslist7[[i]] = ladderize(as.phylo(rtreeshape(n=1, tip.number=numtips, model="aldous")[[1]]))
  Aldoustreeslist7[[i]]$tip.label = yuletreeslist7[[1]]$tip.label 
  alltreeslist7[[2*numtrees+i]]=Aldoustreeslist7[[i]]
}
Coaltreeslist7 = rtreeshape(n=numtrees, tip.number=numtips, model="aldous")
for (i in 1:numtrees) {
  treee=ladderize(rcoal(n=numtips,br=rep(1,numtips-1)))
  treee$tip.label = yuletreeslist7[[1]]$tip.label
  alltreeslist7[[3*numtrees+i]]=treee
}
distanceval = dist.multiPhylo(alltreeslist7, 
                              method ="geodesic", force.multi2di=FALSE, outgroup=NULL,convert.multifurcating=FALSE, 
                              use.random.resolution=FALSE, scale=NULL, verbose=FALSE)
distanceval=as.matrix(distanceval)
namevec=rep(0,100)
for(i in 1:4){
  for(j in 1:numtrees){
    namevec[j+(i-1)*numtrees]=i
  }
}
fit <- cmdscale(distanceval, eig=TRUE, k=2) # k is the number of dim (needs MASS PACKAGE)
fit # view results
row.names(distanceval) = namevec
# plot solution 
x3 <- fit$points[,1]
y3 <- fit$points[,2]
plot(x3, y3, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="MDS Comparing Tree Models, n=5", type="n")
text(x3, y3, labels = row.names(distanceval), cex=1.0)
text(x3[1:25], y3[1:25], labels = row.names(distanceval)[1:25], col="blue", cex=1.0)
text(x3[26:50], y3[26:50], labels = row.names(distanceval)[26:50], col="red", cex=1.0)
text(x3[51:75], y3[51:75], labels = row.names(distanceval)[51:75], col="green", cex=1.0)
text(x3[76:100], y3[76:100], labels = row.names(distanceval)[76:100], col="orange", cex=1.0)
#Adding Labels for Centroids
text(mean(x3[1:25]), mean(y3[1:25]), col="blue", cex=1.0, "Yule")
text(mean(x3[26:50]), mean(y3[26:50]), col="red", cex=1.0, "PDA")
text(mean(x3[51:75]), mean(y3[51:75]), col="green", cex=1.0, "Aldous")
text(mean(x3[76:100]), mean(y3[76:100]), col="orange", cex=1.0, "Coalescent")