##########################
#setClass("AWFE4",
#	representation(SetShapes="list",
#			IntersectionShapes="list",
#			IntersectionMidpoints="list",
#			universe.x="numeric","universe.y"="numeric",
#			"Venn",areas="data.frame"))


###############################
# the Edwards construction for 4 sets; used by the Chow Ruskey algorithm
makeAWFE.all <- function(n,V=Venn()) { 
	if(n>6) { stop("Not implemented for n>6")}

	data(AWFEto6) # nb this data set is created by the AWFE vignette
	An <- AWFEto$A[[n]]; Dn <- AWFEto$D[[n]]; S <- AWFEto$S
	new("AWFE",A=An,D=Dn,S=S[1:n],V)
}

makeAWFE <- function(n=4) {
	AWFE <- makeAWFE.all(n)
	AWFE@A
}

makeAWFE4 <- function() {
	nodes = paste("p",1:14,sep="")
	AWFE44.V <- new("graphNEL",nodes=nodes,
		edgemode="directed"
		)
	edgeDataDefaults(AWFE44.V,"Set") <- NA
	edgeDataDefaults(AWFE44.V,"Fixed") <- FALSE

	addSet <- function(SetNumber,s1,gr) {
		for (ix in seq_along(s1)) {
			s1a <- paste("p",c(s1,s1[1]),sep="")
			gr <- addEdge(s1a[ix],s1a[ix+1],gr)
			edgeData(gr,s1a[ix],s1a[ix+1],"Set") <- SetNumber
		}
		gr
	}

	AWFE44.V <- addSet(1,c(1,2,3,4,5,6),AWFE44.V )
	AWFE44.V <- addSet(2,c(1,7,8,4,9,10),AWFE44.V )
	AWFE44.V <- addSet(3,c(2,11,8,12,6,13,9,14),AWFE44.V )
	AWFE44.V <- addSet(4,c(11,7,12,5,13,10,14,3),AWFE44.V )
	
	AWFE44.V
}

# a function for debuggging
.plotAWFE4 <- function() {
	AWFE44.V <- makeAWFE4()

	en <- edgeNames(AWFE44.V)
	ed <- (edgeData(AWFE44.V,,attr="Set"))
	ed <- lapply(ed,function(Set){c("red","green","black","blue")[Set]})
	names(ed) <- en
	eattrs <- list("color"=ed)
	library(Rgraphviz)

	plot(AWFE44.V,edgeAttrs=eattrs) #colour does
}


cutAWFE4 <- function(AWFE4.V) {
	# produce a cut in the graph
	CK.cutG <- AWFE4.V
	CK.cutG <- removeEdge("p9","p14",CK.cutG)
	CK.cutG <- removeEdge("p10","p14",CK.cutG)
	CK.cutG <- removeEdge("p10","p1",CK.cutG)
	CK.cutG <- removeEdge("p6","p1",CK.cutG)
	CK.cutG
}



####################
# decorate the dual of an n-cube graph with the topology of the Edwards graph
makeAWFEstar <- function(n=4) { 
	AWFE <- makeAWFE.all(n)
	AWFE@D
}


makeAWFE4star <- function(n=4) {
	stopifnot(n==4)
	AWFE44.Vs <- makeQn(4)
	nodeDataDefaults(AWFE44.Vs,"Region") <- new("graphNEL")
	AWFE44.V <- makeAWFE4() 
	nodeData(graph=AWFE44.Vs,"1111","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(4,9,14,3)))
	nodeData(graph=AWFE44.Vs,"0111","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(3,11,4,8)))
	nodeData(graph=AWFE44.Vs,"1011","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(4,5,9,13)))
	nodeData(graph=AWFE44.Vs,"1101","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(14,9,10)))
	nodeData(graph=AWFE44.Vs,"1110","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(2,3,14)))
	nodeData(graph=AWFE44.Vs,"0011","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(4,5,8,12)))
	nodeData(graph=AWFE44.Vs,"0101","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(7,8,11)))
	nodeData(graph=AWFE44.Vs,"0110","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(2,3,11)))
	nodeData(graph=AWFE44.Vs,"1010","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(5,6,13)))
	nodeData(graph=AWFE44.Vs,"1001","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(9,10,13)))
	nodeData(graph=AWFE44.Vs,"1100","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(1,2,10,14)))
	nodeData(graph=AWFE44.Vs,"1000","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(1,6,10,13)))
	nodeData(graph=AWFE44.Vs,"0100","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(1,2,7,11)))
	nodeData(graph=AWFE44.Vs,"0010","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(5,6,12)))
	nodeData(graph=AWFE44.Vs,"0001","Region") <- subGraph(graph=AWFE44.V,paste(sep="","p",c(7,8,12)))
	AWFE44.Vs
}

cutAWFE.old <- function(AWFE4.V) { cutAWFE4 (AWFE4.V)}
cutAWFE <- function(AWFE) { AWFE <- scythe.AWFE(AWFE); AWFE@A }
scythe.AWFE <- function(AWFE) { 
	# actually only need graph for algorithm but calc dual for debug display
	graph <- AWFE@A; dual=AWFE@D
	dn <- nodes(dual)
	nSets <- nchar(dn[1]) 
	source.sink.path <- sapply(nSets:0,function(ns){
		paste(paste(rep("0",nSets-ns),collapse=""),paste(rep("1",ns),collapse=""), sep="")
		})
	for (six in 1:nSets) {
		from.face.name <- source.sink.path[six]
		to.face.name <- source.sink.path[six+1]
		from.face <- nodeData(dual,from.face.name,attr="face")[[1]]
		to.face <- nodeData(dual,to.face.name,attr="face")[[1]]
		ed <- edgeData(from.face,attr="Set")
		edix <- ed[(sapply(ed,function(x)x==six))]
		scythe.edge <- names(edix)[1]
		fromto.node.name <- strsplit(scythe.edge,split="|",fixed=TRUE)[[1]]
		from.node.name <- fromto.node.name[1]
		to.node.name <- fromto.node.name[2]
		graph <- removeEdge(from.node.name,to.node.name,graph)
		from.face <- removeEdge(from.node.name,to.node.name,from.face)
		to.face <- removeEdge(from.node.name,to.node.name,to.face)
		nodeData(dual,c(from.face.name,to.face.name),"face") <- list(from.face,to.face)
	}
	AWFE@A <- graph
	AWFE@D <- dual
	AWFE
}

###########################################
#
setClass("AWFE",	representation(A="graphNEL",D="graphNEL",	S="list","Venn"))

compute.AWFE <- function(V,doWeights=FALSE) {
	if(doWeights) { warning("Weights ignored by AWFE algorithm") }
	n <- NumberOfSets(V)
	makeAWFE.all(n,V)
}



setMethod("VisibleRange","AWFE",function(object) {
	UniverseRange(object)
})

setMethod("UniverseRange","AWFE",function(object) {
	xrange <- 5*c(-1,1)
	yrange <- xrange
	list(x=xrange,y=yrange)
})

		
setMethod("SetLabels","AWFE",function(object){
		data.frame()
})


twopi <- 2* pi


projection.thetaphi <- function(thetaphi,projection="PS") {
	theta <- thetaphi[,1]; phi <- (thetaphi[,2])
	if (projection=="PS") {
		rho <- cos(phi)/(1-sin(phi))* 2
		x <- rho * cos(theta)
		y <- rho * sin(theta)
		xy <- cbind(x,y)
	} else if (projection=="EC") {
		y <- sin(phi)
		x <- theta
		xy <- cbind(x,y)
	}
	xy
}

setMethod("PlotSetBoundaries","AWFE",function(object,gp){
	S=object@S
	projection="PS"
	s <- seq(0,1,length=100)

	for (Six in  1:length(S)) {
			thexy <- S[[Six]](s,projection)
			grid.lines(thexy[,1],thexy[,2],gp=gp[Six],default.units="native")
		}

})

setMethod("PlotUniverse","AWFE",function(object,gp){
	if(missing(gp)) { gp <- NULL }
	})


thetah.to.xy <- function(thetah,projection) {
	xy <- if (projection=="EC") { 
		thetah
	} 	else { # PS
		phi <- asin(thetah[,2])
		thetaphi <- cbind(thetah[,1],phi)
		projection.thetaphi(thetaphi,projection)
	}
	xy
}

setMethod("IntersectionMidpoints","AWFE",function(object){
	dual <- object@D
	xy <- t(sapply(nodes(dual),function(dn){
		thetah <- nodeData(dual,dn,"thetah")[[1]]
		if (!is.null(thetah)) {
			xy <- thetah.to.xy(thetah,projection="PS")
			xy
			}
	}))

	# pull out the outer triangle vertices
	VLabels <- data.frame(Label=nodes(dual),x=NA,y=NA,hjust=I("center"),vjust=I("center"))
	VLabels[,c("x","y")] <- xy
	names(VLabels)[2:3] <- paste("Midpoint.",sep="",names(VLabels)[2:3])
	VLabels
})

#setMethod("Areas","AWFE4",function(object){
#	Areas <- sapply(object@IntersectionShapes,.polygon.area)
#	dark.matter <- paste(rep("0",NumberOfSets(object)),collapse="")
#	Areas[dark.matter] <- Areas[dark.matter]-sum(Areas[names(Areas)!=dark.matter])
#	Areas
#})


setMethod("PlotUniverse","AWFE",function(object,gp){
	ur <- UniverseRange(object)
	grid.rect(x=mean(ur$x),y=mean(ur$y),width=diff(ur$x),height=diff(ur$y),default.units="native")
})

setMethod("PlotNodes","AWFE",function(object,gp){
	graph <- object@A
	nn <- nodes(graph)
	thetah <- do.call("rbind",nodeData(graph,nn,"thetah"))
	xy <- thetah.to.xy(thetah,projection="PS")
	grid.text(x=xy[,1],y=xy[,2],label=nn,default.units="native",just=c("right","top")) 
})

setMethod("PlotFaces","AWFE",function(object,gp,arrow){
	if (missing(gp)) { gp  <- gpar(col=trellis.par.get("superpose.symbol")$col) }

	if(missing(arrow)) { 	arrow <- grid::arrow(length=unit(2,"mm"),type="closed")} 
	graph <- object@A
	dual <- object@D
	for (dn in nodes(dual)) {

	face <- nodeData(dual,dn,"face")[[1]]
	for (face.from in nodes(face)) {
		for (face.to in edges(face,face.from)[[1]]) {
			earc <- edgeData(graph,face.from,face.to,"arc")[[1]]
			Six<- edgeData(graph,face.from,face.to,"Set")[[1]]

			xy <- arc.to.xy(earc,nintervals=100,projection="PS")
	 		grid.lines(x=xy[,1],y=xy[,2],default.units="native",gp=gpar(gp[Six]),arrow=arrow) 
			}
		}
	}
})


