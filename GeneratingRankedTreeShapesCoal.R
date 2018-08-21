install.packages("ape")
install.packages("apTreeshape")
install.packages("distory")
install.packages("MASS")
library("ape")
library("apTreeshape")
library("distory")
library("MASS")
create_F<-function(tree){
  
  #n is the number of individual samples 
  edges<-tree$edge
  x<-tree$Nnode+1
  n_samples<-x
  values1<-max(ape::node.depth.edgelength(tree))-ape::node.depth.edgelength(tree)
  values1<-values1[(x+1):length(values1)]
  values<-sort(values1,decreasing=T)
  correct.label<-seq(x+1,2*x-1)
  ##do the crossreference
  reference<-matrix(c(values1,seq(x+1,2*x-1),rep(0,length(values1))),ncol=3)
  
  for (j in 1:length(values1)){
    reference[reference[,1]==values[j],3]<-correct.label[j]
  }
  newedges<-matrix(0,nrow=nrow(edges),ncol=4)
  newedges[,1:2]<-edges
  newedges[,3:4]<-edges
  
  for (j in (x+1):(2*x-1)){
    newedges[newedges[,1]==j,3]<-reference[reference[,2]==j,3]
    newedges[newedges[,2]==j,4]<-reference[reference[,2]==j,3]
  }
  edges<-newedges[,3:4]
  #tree.edge is edges correcting the labels
  F<-matrix(0,nrow=n_samples-1,ncol=n_samples-1)
  diag(F)<-n_samples:2
  F[cbind(2:(n_samples-1),1:(n_samples-2))]<-diag(F)[-(n_samples-1)]-2
  x<-nrow(edges)
  y<-max(edges)
  
  
  for (cini in 1:(n_samples-2)){
    ini<-cini+1
    for (j in y:(n_samples+2) ){
      F[ini,cini]<-F[(ini-1),cini]-sum(edges[edges[,1]==j,2]<=n_samples)
      ini<-ini+1
    }
    #remove those two rows
    find<-seq(1,nrow(edges))[edges[,1]==y]
    edges<-edges[-find,]
    edges[edges[,2]==y,2]<-1
    y<-y-1
  }
  return(F)
}
#simulate a tree from the standard coalescent model
library("ape","phylodyn")
tree<-rcoal(n=10,br=rep(1,9))
plot(ladderize(tree),show.tip.label=FALSE)
axisPhylo()
F<-create_F(tree)
#this function takes the F matrix and generate the Adjacency information
parentlist<-function(F){
  n<-F[1,1]
  Fdiff<-F[1:(n-2),]-F[2:(n-1),]
  parentlist<-matrix(0,nrow=2*(n-1),ncol=2)
  parentlist[,1]<-sort(rep(seq(1,n-1),2))
  for (i in 2:(n-2)){
    if (Fdiff[i,1]==1){parentlist[parentlist[,1]==i,2]<-c(0,min(seq(0,n-2)[Fdiff[i,]==2]))}
    if (Fdiff[i,1]==0){parentlist[parentlist[,1]==i,2]<-c(min(seq(0,n-2)[Fdiff[i,]==1]),min(seq(0,n-2)[Fdiff[i,]==2]))}
  }
  parentlist[parentlist[,1]==(n-1),2]<-c(min(seq(0,n-2)[F[n-1,]==1]),min(seq(0,n-2)[F[n-1,]==2]))
  return(parentlist)
}

adje<-parentlist(F)

#this function takes adjacency data and generates the F matrix

Ffromparentlist<-function(parent.list){
  n_samples<-max(parent.list[,1])+1
  F <- matrix(0, nrow = n_samples - 1, ncol = n_samples - 1)
  diag(F) <- n_samples:2
  for (j in 2:(n_samples-1)){
    numbersing<-sum(parent.list[parent.list[,1]==j-1,2]==0)
    if (numbersing==2) {F[j,1:(j-1)]<-F[j-1,1:(j-1)]-2}
    if (numbersing==1) {
      vin<-max(parent.list[parent.list[,1]==j-1,2])
      F[j,1:vin]<-F[j-1,1:vin]-1
      F[j,(vin+1):(j-1)]<-F[j-1,(vin+1):(j-1)]-2
    }
    if (numbersing==0) {
      vinm<-max(parent.list[parent.list[,1]==j-1,2])
      vin<-min(parent.list[parent.list[,1]==j-1,2])
      F[j,1:vin]<-F[j-1,1:vin]
      F[j,(vin+1):(vinm)]<-F[j-1,(vin+1):(vinm)]-1
      F[j,(vinm+1):(j-1)]<-F[j-1,(vinm+1):(j-1)]-2
    }
  }
  return(F)  
}

Fback<-Ffromparentlist(adje)

sum(abs(Fback-F)) ##same ranked tree shape?

#This function takes a matrix F and a vector of times and generates a phylo tree

tree_from_F <- function(matF, coal_times){
  #generate an ape tree 
  #F is the actual form used in the code that differs from the paper's notation
  
  n= dim(matF)[1]
  edge=matrix(rep(0,6*n-6),ncol=3)
  edge[,2]= (1:(2*n-2))
  vintages = c()
  times=c(rep(0,n),coal_times)
  for (j in n:2){
    new_node = 2*n-j+1
    next_leaf = intersect(which(edge[,1]==0),1:n)[1]
    F_difference = rev(matF[,j]-matF[,j-1])
    
    if (F_difference[1]==2){
      edge[next_leaf:(next_leaf+1),1]=new_node
      vintages=c(vintages, new_node)
    }
    else if (F_difference[1]==1){
      selected_vintage = which(F_difference == 2)[1]+n-1
      edge[selected_vintage,1]=new_node
      edge[next_leaf,1]=new_node
      vintages = c(vintages[vintages != selected_vintage],new_node)
    }
    else {
      selected_vintage1 =which(F_difference == 1)[1]+n-1
      selected_vintage2 =which(F_difference == 2)[1]+n-1
      edge[selected_vintage1,1]=new_node
      edge[selected_vintage2,1]=new_node
      vintages = vintages[vintages!=selected_vintage1]
      vintages = vintages[vintages!=selected_vintage2]
      vintages<-c(vintages,new_node)
    }
  }
  #edge[5:8,2]=c(6,7,8,5)
  edge[1:n,]=edge[order(edge[1:n,2]),]
  #edge=edge[order(edge[,1]),]
  
  for (j in 1:(2*n-2)) {
    #I don't understand this
    edge[j,3]=times[edge[j,1]]-times[edge[j,2]]
  }
  edge[,1]=3*n-edge[,1]
  edge[-(1:n),2]=3*n-edge[-(1:n),2]
  
  final_tree=rcoal(n,br=coal_times)
  final_tree$edge=edge[,-3]
  final_tree$edge.length=edge[,3]
  final_tree$Nnode=n-1
  class(final_tree) <- "phylo"
  final_tree <- reorder(final_tree,"postorder")
  final_tree$edge[final_tree$edge[, 2] <= n, 2] <- 1:n
  return(final_tree)
}

#auxiliary function, you can ignore it
consolidateF<-function(Fju){
  ##generates an F matrix consistent with paper notation
  newF<-matrix(0,nrow=nrow(Fju),ncol=nrow(Fju))
  for (j in 1:nrow(newF)){
    newF[nrow(Fju)-j+1,]<-rev(Fju[,j])
  }
  newF2<-matrix(0,nrow=Fju[1,1],ncol=Fju[1,1])
  newF2[2:Fju[1,1],2:Fju[1,1]]<-newF
  return(newF2)
}

originalTree<-tree_from_F(consolidateF(Fback),seq(1,9))
par(mfrow=c(1,2))
plot(ladderize(originalTree),show.tip.label=FALSE)
plot(ladderize(tree),show.tip.label=FALSE)

#simulate a list of trees from the standard coalescent model
numtrees=500
numtips=5
library("apTreeshape","ape","phylodyn")
simulationtreeslist=rtreeshape(n=numtrees, tip.number=numtips, model="yule")
for (i in 1:numtrees) {
  simtree=rcoal(n=numtips,br=rep(1,numtips-1))
  simulationtreeslist[[i]] = ladderize(simtree)
} 

rankedlist <- list()#going to store all the trees of the same ranked tree shape in a list (list of lists)
list1 <- list()
list1[[length(list1)+1]] = (simulationtreeslist[[1]])
rankedlist[[1]]= list1
unique <- TRUE

for (i in 2:numtrees) {
  for(j in 1:length(rankedlist)){
    if(sum(abs(parentlist(create_F(simulationtreeslist[[i]]))-parentlist(create_F(rankedlist[[j]][[1]]))))==0){#if trees are the same
      #append(rankedlist[j],simulationtreeslist[[i]])
      unique=FALSE
      rankedlist[[j]][[length(rankedlist[[j]])+1]]=simulationtreeslist[[i]] #add current tree to this list
      break
    }
  }#if going into this code this tree should have a unique shape
  if(unique==TRUE){
    newlist = list()#making new list
    newlist[[length(newlist)+1]] = simulationtreeslist[[i]]
    rankedlist[[length(rankedlist)+1]] = newlist #adding new list with unique tree to rankedlist
    #rankedlist[length(rankedlist)+1]=newlist  
  }
  unique=TRUE
}
heightvector = rep(0,5)
for(i in 1:5){
  heightvector[i] = length(rankedlist[[i]])
}
barplot(heightvector, main="Ranked Tree Shape n=5, 500 Trees", 
        xlab="Tree Type", ylab="Frequency")

#BHV MIN distance
dimmat = length(rankedlist)
mindist=10000
distmatrix3 = matrix(0, nrow=dimmat, ncol=dimmat)
for(i in 1:dimmat){
  for(j in i:dimmat){
    for(k in 1: length(rankedlist[[i]])){
      for(m in 1: length(rankedlist[[j]])){
        distanceval = dist.multiPhylo(list(rankedlist[[i]][[k]],rankedlist[[j]][[m]]), 
                                      method ="geodesic", force.multi2di=FALSE, outgroup=NULL,convert.multifurcating=FALSE, 
                                      use.random.resolution=FALSE, scale=NULL, verbose=FALSE)
        if(distanceval<mindist){
          mindist=distanceval
          distmatrix3[i,j]= distanceval
          distmatrix3[j,i]= distanceval
        }
      }
    }
    mindist=10000
  }
}
for(i in 1:dimmat){
  distmatrix3[i,i]=0
}

# Classical MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

fit <- isoMDS(distmatrix3, k=2) # k is the number of dim
fit # view results
row.names(distmatrix3) = (1:dimmat)
# plot solution 
x3 <- fit$points[,1]
y3 <- fit$points[,2]
plot(x3, y3, xlab="Coordinate 1", ylab="Coordinate 2", 
     main="Non-metric_MDS - Ranked Tree Shapes n=5 ", type="n")
text(x3, y3, labels = row.names(distmatrix3), cex=.7)