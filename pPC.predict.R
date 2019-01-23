pPC.predict<-function(ORIGINAL_TREE, ORIGINAL_DATA_MATRIX, NEW_DATA_MATRIX){ #adapted from phytools phyl.pca code; TAKES DATA USED TO CONSTRUCT pPC AXES AND A NEW SET OF DATA AND PROJECTS pPC SCORES FOR THOSE NEW DATA
	if(dim(ORIGINAL_DATA_MATRIX)[2]!=dim(NEW_DATA_MATRIX)[2]){stop("Error: new data has a different number of variables from original data")}
	
	C.a<-vcv.phylo(ORIGINAL_TREE)[rownames(ORIGINAL_DATA_MATRIX),rownames(ORIGINAL_DATA_MATRIX)]

	temp<-phyl.vcv(ORIGINAL_DATA_MATRIX,C.a,lambda=1)
	V <- temp$R
	a <- t(temp$alpha)
	#C <- temp$C

	n<-dim(ORIGINAL_DATA_MATRIX)[1]
	m<-dim(ORIGINAL_DATA_MATRIX)[2]

	es = eigen(V)
	Evec <- es$vectors[, 1:min(n - 1, m)]
    #A <- matrix(rep(a, n), n, m, byrow = T)
	#original scores should be given by: (ORIGINAL_DATA_MATRIX - A) %*% Evec

	A.new<-matrix(rep(a, dim(NEW_DATA_MATRIX)[1]), dim(NEW_DATA_MATRIX)[1], m, byrow = T)
	simulated.pPC<-(NEW_DATA_MATRIX-A.new)%*%Evec

	return(simulated.pPC)
}
