warning("Entering ThreeSquares")


	
######################

# TODO need to implement DrawnVenn for these classes and use that
# to provide Areas and PlotFaces methods
setClass("TwoSquares",representation(x="numeric",y="numeric",width="numeric",height="numeric",
	universe.x="numeric",universe.y="numeric","DrawnVenn","Venn" ))

compute.S2 <- function(V,doWeights=TRUE,doEuler=FALSE) {
	stopifnot(NumberOfSets(V)==2)
	Weight <- Weights(V)
	if (!doWeights){
		Weight <- 1+0*Weight
		if (doEuler) {
			Weight[Weights(V)==0] <- 0
		} 
	}
	wab <- Weight["00"]
	wAb <- Weight["10"]
	waB <- Weight["01"]
	wAB <- Weight["11"]

	wA <- wAb+ wAB
	wB <- waB + wAB

	s1 <- sqrt(wA)
	s2 <- sqrt(wB)

	# squares are centred at (-+d/2,0) and have sides s1, s2
	if (wAb==0) {
		#d <- (s2-s1)/4 to have it completely inside
		d <- (s2-s1)/2
		if (!doEuler) {
			d <- (s2-s1)/2+ 0.1*s1
		}
	} else if (waB==0) {
		#d <- (s1-s2)/4
		d <- (s1-s2)/2
		if (!doEuler) {
			d <- (s1-s2)/2+ 0.1*s2
		}
	} else {
		d <- (s1+ s2)/2 -   wAB / min(s1,s2)
		if(wAB==0 & !doEuler) {
			d <- d - .1* min(s1,s2)
		}
	}

	l1 <- -d/2-s1/2; l2 <- d/2-s2/2
	r1 <- -d/2+s1/2; r2 <- d/2+s2/2

	# first we construct as though no intersection
	diagram <- new("graphNEL",nodes=c("p1","p2"),edgemode="directed")
	region <- new("graphNEL",nodes=IndicatorString(V),edgemode="directed")
	emptyxy <- matrix(NA,ncol=2,nrow=1)

	nodeDataDefaults(diagram,"xy") <- emptyxy 
	nodeDataDefaults(region,"xy") <- emptyxy 
	edgeDataDefaults(diagram,"arc") <- new("edgeLines")
	if (FALSE) { # code that is not finished...
	r01 <- diagram
	r11 <- diagram
	r10 <- diagram
	arcs.1 <- new("edgeLines",Set=1,xy=cbind(c(l1,r1,r1,l1),c(-s1/2,-s1/2,s1/2,s1/2)))
	arcs.2 <- new("edgeLines",Set=2,xy=cbind(c(l2,r2,r2,l2),c(-s2/2,-s2/2,s2/2,s2/2)))
	diagram <- addEdge("p1","p1",diagram)
	edgeData(diagram,"p1","p1","arc") <- list(arcs.1)
	diagram <- addEdge("p2","p2",diagram)
	edgeData(diagram,"p2","p2","arc") <- list(arcs.2)
	nodeData(diagram,"p1","xy") <- list(arcs.1@xy[1,,drop=FALSE])
	nodeData(diagram,"p2","xy") <- list(arcs.2@xy[1,,drop=FALSE])

	if ( l2<=l1 & r2>= r1) { # 1 in 2
		diagram <- addNode("q1",diagram)
		nodeData(diagram,"q1","xy") <- list(cbind(l2,s1/2))
	} else if ( l1<=l2 & r1>= r2) { # 2 in 1
		diagram <- addNode(c("q1","q2"),diagram)
		nodeData(diagram,c("q1","q2"),"xy") <- list(cbind(r2,s2/2),cbind(r2,-s2/2))
	} else if (l1 < l2 & l2 <= r2) { # 1 intersect 2 
		diagram <- addNode(c("q1","q2"),diagram)
		if (s1 > s2) { xix <- r1; yix <- s2/2 } else { xix <- l2; yix <- s1/2}
		nodeData(diagram,c("q1","q2"),"xy") <- list(cbind(xix,yix),cbind(xix,-yix))
	}	
	} # unfinished code


	
	drawnVenn <- new("DrawnVenn",regions=region,diagram=diagram)

	S2 <- new("TwoSquares",x=c(-d/2,d/2),y=c(0,0),width=c(s1,s2),height=c(s1,s2),
		drawnVenn,V)
	S2 <- .square.universe (S2,doWeights=doWeights)
	S2
}
setMethod("PlotSetBoundaries","TwoSquares",function(object,gp){
	if(missing(gp)) { gp <- NULL }
	grid.rect(x=object@x,y=object@y,width=object@width,height=object@height,default.unit="native",gp=gp)
	}
 )
setMethod("PlotUniverse","TwoSquares",.plotuniverse.square)
setMethod("PlotNodes","TwoSquares",function(object,gp){
	warning("PlotNodes not implemented for TwoSquares")
})

setMethod("VisibleRange","TwoSquares",function(object) {
	centres.x <- object@x
	centres.y <- object@y
	width <- object@width
	height <- object@height
	xrange <- c(min(centres.x-width/2),max(centres.x+width/2))
	yrange <- c(min(centres.y-height/2),max(centres.y+height/2))
		list(x=xrange ,y=yrange)

})
setMethod("UniverseRange","TwoSquares",.universe.range)

setMethod("IntersectionMidpoints","TwoSquares",function(object) {
	V <- as(object,"Venn")
	VI <- data.frame(Indicator(V))
	VS <- IndicatorString(V)
	VI$Midpoint.x <- NA
	VI$Midpoint.y <- NA
	VI$hjust <- "center"
	VI$vjust <- "bottom"
	s1 <- object@width[1]; s2 <- object@width[2]
#	 choose the side with the smaller circle
	goleft <- s1 < s2
	d <- object@x[2]-object@x[1]
	smidge.x <- 0.01 * diff(UniverseRange(object)$x)
	smidge.y <- 0.01 * diff(UniverseRange(object)$y)
	ybottom <- min(UniverseRange(object)$y) + smidge.y
	if (goleft) {
		dark.matter.x <- min(UniverseRange(object)$x)+smidge.x
		dark.matter.hjust <- "left"
	} else {
		dark.matter.x <- max(UniverseRange(object)$x)-smidge.x
		dark.matter.hjust <- "right"
	}

	VI[VS=="00" ,c("Midpoint.x","Midpoint.y")] <- c(dark.matter.x,ybottom)
	x10 <- if(d+s1/2<s2/2) {  # 1 inside 2
		NA
	} else {
		object@x[1]-((s1+s2)/2-d)/2
	}
	VI[VS=="10" ,c("Midpoint.x","Midpoint.y")] <- c(x10,0)

	x01 <- if(d+s2/2<s1/2) {  # 2 inside 1
		NA
	} else {
		object@x[2]+((s1+s2)/2-d)/2
	}
	VI[VS=="01" ,c("Midpoint.x","Midpoint.y")] <- c(x01,0)
	
	x11 <- if (d>=(s1+s2)/2) { NA } else { object@x[1]+s1/2-(s1/2+s2/2-d)/2 } 
	VI[VS=="11" ,c("Midpoint.x","Midpoint.y")] <- c(x11,0)
	VI[VS=="00","hjust"] <- dark.matter.hjust

	VI
})

	

setMethod("SetLabels","TwoSquares",function(object){
	# unlikely to be right for general circles
	# also assumes particular order
	y=object@y;x=object@x;s=object@width
	VLabels <- data.frame(Label=rep("unset",2),x=NA,y=NA,hjust=I("center"),vjust=I("top"))
	VLabels[1,2:3] <- c( x[1]-s[1]/2,y[1]+s[1]/2 );  VLabels$hjust[1] <- "left"
	VLabels[2,2:3] <- c( x[2]+s[2]/2, y[2]+s[2]/2 );VLabels$hjust[2] <- "right"
	VLabels
})


setMethod("PlotFaces","TwoSquares",function(object,gp,arrow){
	warning("PlotFaces not implemented for class TwoSquares")
})

setMethod("Areas","TwoSquares",function(object){
	warning("Areas not implemented for class TwoSquares")
})





###########################
setClass("ThreeSquares",
	representation(SetShapes="list",
			IntersectionShapes="list",
			universe.x="numeric","universe.y"="numeric",
			"DrawnVenn", # implementation only to define position of nodes
			"Venn",areas="data.frame"))

################################

compute.S3 <- function(V,doWeights=TRUE) {
	stopifnot(NumberOfSets(V)==3)
	wght <- Weights(V)
	VS <- IndicatorString(V)
		
	if (!doWeights) { 
		x <- 1; y <- 3
		wght <- y + 0 * wght
		wght["111"] <- wght["101"] <- wght["011"] <- wght["001"] <- x
	}

	# currently easier thant TwoSqures because we don't worry about zero weights
	diagram <- new("graphNEL",nodes=c("p1","p2"),edgemode="directed")
	region <- new("graphNEL",nodes=IndicatorString(V),edgemode="directed")
	emptyxy <- matrix(NA,ncol=2,nrow=1)

	nodeDataDefaults(diagram,"xy") <- emptyxy 
	nodeDataDefaults(region,"xy") <- emptyxy 
	edgeDataDefaults(diagram,"arc") <- new("edgeLines")
	diagram <- addNode(c("b1","b2","c1","c2","c3","c4"),diagram) # named cf AWFE nodes

	IntersectionShapes <- vector("list",length(VS))
	names(IntersectionShapes) <- VS
	
	wabc <- wght[VS=="111"] 
	d1 <- 1
	wabc.width <- sqrt(wabc)/d1
	wabc.height <- sqrt(wabc)*d1
	wabc.x <- c(-wabc.width/2,wabc.width/2,wabc.width/2,-wabc.width/2)
	wabc.y <- c(-wabc.height/2,-wabc.height/2,wabc.height/2,+wabc.height/2)

	IntersectionShapes[["111"]] <- matrix(c(wabc.x,wabc.y),byrow=FALSE,ncol=2)
	nodeData(diagram,c("b1","c4","c1"),"xy") <- list(IntersectionShapes[["111"]][3,],
									IntersectionShapes[["111"]][2,],
									IntersectionShapes[["111"]][4,])	
	wac <- wght[VS=="101"] 
	wac.height <- wac/wabc.width
	wac.x <- wabc.x
	wac.y <-  wabc.height/2+c(0,0,wac.height,wac.height)
	IntersectionShapes[["101"]] <- matrix(c(wac.x,wac.y),byrow=FALSE,ncol=2)

	wbc <- wght[VS=="011"] 
	wbc.height<- wabc.height
	wbc.width <- wbc/wbc.height
	wbc.y <- wabc.y
	wbc.x <-  wabc.width/2+c(0,wbc.width,wbc.width,0)
	IntersectionShapes[["011"]] <- matrix(c(wbc.x,wbc.y),byrow=FALSE,ncol=2)

	wab <- wght[VS=="110"]
	# (d +wabc.height)^2=wab+wabc
	#d <- 0.5 * wab/wabc.width
	# 
	d<- sqrt(wab+wabc)-wabc.height
	stopifnot(d>=0)
	wab.width <- (wab- d* wabc.width)/(wabc.height+d) + wabc.width
	wab.height <- wabc.height+d
	wab.x <-  wabc.width/2 + c(-wab.width,0,0,-wabc.width,-wabc.width,-wab.width)
	wab.y <- wabc.height/2+c(-wab.height,-wab.height,-wabc.height,-wabc.height,0,0) 
	IntersectionShapes[["110"]] <- matrix(c(wab.x,wab.y),byrow=FALSE,ncol=2)

	nodeData(diagram,c("b2","c2","c3"),"xy") <- list(IntersectionShapes[["110"]][6,],
									IntersectionShapes[["101"]][3,],
									IntersectionShapes[["011"]][2,])	

	wc <- wght[VS=="001"]
	wc.height <- (wac.height+wabc.height)
	if (wc < wbc.width* wc.height - wbc ) { # not square
		dbc <- wc/(wbc.width+wc.height)
		wc.x <- wabc.width/2+  c(wbc.width,wbc.width+dbc, wbc.width+dbc,           dbc,      dbc,        0,         0,wbc.width)
		wc.y <- -wabc.height/2+c(        0,            0,wbc.height+dbc,wbc.height+dbc,wc.height,wc.height,wbc.height,wbc.height)
	} else {
		wc.width <- (wbc+wc)/wc.height # > wbc.width
		wc.x <- wabc.width/2+  c(wbc.width,wc.width,wc.width,        0,         0,wbc.width)
		wc.y <- -wabc.height/2+c(        0,                 0,          wc.height,wc.height,wbc.height,wbc.height)
	}
	IntersectionShapes[["001"]] <- matrix(c(wc.x,wc.y),byrow=FALSE,ncol=2)

	wb <- wght[VS=="010"]
	if (!doWeights) { # special is case to send 010 the other way
		wb.height <- wab.height
		wb.width <- ( wbc+wb)/wb.height
		wb.x <- wabc.width/2+c(0,wb.width,wb.width,wbc.width,wbc.width,0)
		wb.y <- wabc.height/2-wb.height+c(0,0,wb.height,wb.height,wb.height-wbc.height,wb.height-wbc.height)
	} else {
	wb.width <- wab.width + wbc.width
	if ( wb < wbc.width*(wab.height-wabc.height) ) {
		xa <- 2; xb <- wab.width + wbc.width +wab.height; xc <- - wb
		db <- (-xb+sqrt(xb^2-4*xa*xc))/(2*xa)
		wb.width <- wb.width+db
		wb.height <- wb.height+db
		wb.x <- -wabc.width/2-(wab.width-wabc.width)-db+
			c(0,wab.width+2*db,wab.width+2*db,wb.width,wb.width,wab.width+db,wab.width+db,db,db,0)
		wb.y <- wabc.height/2-wb.height+
			c(0,0,wab.height-wabc.height,wab.height-wabc.height,
			wb.height-wabc.height,wb.height-wabc.height,db,db,wb.height,wb.height)
	} else {
		wb.height <- (wabc+wbc+wab+wb)/wb.width
		# that's the height giving no overlap with a
		sf <- (wb.height/wab.height) # > 1, the extent below
		db <- sqrt(sf)
		wb.height <- wb.height/db
		wb.width <- wb.width*db
		wb.x <- +wabc.width/2+wbc.width-wb.width+
			c(0,wb.width,wb.width,wb.width-wbc.width,wb.width-wbc.width,wb.width-wbc.width-wab.width,wb.width-wbc.width-wab.width,0)
		wb.y <- wabc.height/2-wb.height+
			c(0,0,wb.height-wabc.height,wb.height-wabc.height,wb.height-wab.height,wb.height-wab.height,wb.height,wb.height)
	}
	} # else (!doWeights)
	IntersectionShapes[["010"]] <- matrix(c(wb.x,wb.y),byrow=FALSE,ncol=2)

	wa <- wght[VS=="100"]

	wa.width <- wab.width 
	if ( wa < (wab.width-wabc.width)*wab.height ) {
		da <- wa/(wa.width+ wac.height)
		wa.x <- -wabc.width/2-(wab.width-wabc.width)+
			c(0,wab.width-wabc.width,wab.width-wabc.width,wab.width,wab.width,wab.width-wabc.width-da,wab.width-wabc.width-da,0)
		wa.y <- wabc.height/2+
			c(0,0,wac.height,wac.height,wac.height+da,wac.height+da,da,da)
	} else {
		wa.height <- (wac+wb)/wa.width
		wa.x <- -wabc.width/2-(wab.width-wabc.width)+
			c(0,wab.width-wabc.width,wab.width-wabc.width,wab.width,wab.width,0)
		wa.y <- wabc.height/2+
			c(0,0,wac.height,wac.height,wa.height,wa.height)
	}
	IntersectionShapes[["100"]] <- matrix(c(wa.x,wa.y),byrow=FALSE,ncol=2)

		SetShapes <- list()
	SetShapes[[1]] <- do.call(rbind,list(
		IntersectionShapes[["110"]][1:2,],IntersectionShapes[["100"]][5:nrow(IntersectionShapes[["100"]]),]
		))
	SetShapes[[2]] <- if (doWeights) {
		do.call(rbind,list(
		IntersectionShapes[["010"]][1:(nrow(IntersectionShapes[["010"]])-6),],IntersectionShapes[["011"]][3,],
			IntersectionShapes[["010"]][nrow(IntersectionShapes[["010"]]),]
		))
		} else {
		do.call(rbind,list(
		IntersectionShapes[["110"]][1,],
		IntersectionShapes[["010"]][2:3,],
		IntersectionShapes[["110"]][6,]
		))
		}
	SetShapes[[3]] <- do.call(rbind,list(
		IntersectionShapes[["111"]][1,],IntersectionShapes[["001"]][2:(nrow(IntersectionShapes[["001"]])-3),],
			IntersectionShapes[["101"]][4,]
		))
	drawnVenn <- new("DrawnVenn",regions=region,diagram=diagram)
	res <- 	new("ThreeSquares",SetShapes=SetShapes,
		IntersectionShapes=IntersectionShapes,
		drawnVenn ,
		V)
	res <- .square.universe(res,doWeights=doWeights)
	res@IntersectionShapes[["000"]] <- matrix( c(res@universe.x[c(1,2,2,1)],res@universe.y[c(1,1,2,2)]),byrow=FALSE,ncol=2)
	res
}

setMethod("VisibleRange","ThreeSquares",function(object) {
	object@IntersectionShapes[["000"]] <- NULL
	xrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,1])}))	
	yrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,2])}))
	list(x=xrange,y=yrange)
})

setMethod("UniverseRange","ThreeSquares",function(object){
	xrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,1])}))	
	yrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,2])}))
	list(x=xrange,y=yrange)
})


		
setMethod("IntersectionMidpoints","ThreeSquares",function(object){
	# three regions are squares and we want the mean
	VI <- Indicator(as(object,"Venn"))
	mids <- 	data.frame(t(sapply(object@IntersectionShapes,function(x)apply(x,2,mean))))
	colnames(mids) <- c("Midpoint.x","Midpoint.y")
	mids[sapply(object@IntersectionShapes,function(x)nrow(x)!=4),] <- NA
	mids$IndicatorString <- names(object@IntersectionShapes)
	mids <- cbind(VI,mids)

	IS <- object@IntersectionShapes
	smidge <- 0.05 * sapply(UniverseRange(object),diff)
	mids[mids$IndicatorString=="000",c("Midpoint.x","Midpoint.y")] <-c(max(object@universe.x),max(object@universe.y))-smidge
	mids[mids$IndicatorString=="110",c("Midpoint.x","Midpoint.y")] <-apply(IS[["110"]][c(1,4),],2,mean)
	mids[mids$IndicatorString=="100",c("Midpoint.x","Midpoint.y")] <- apply(rbind(IS[["100"]][6,],IS[["101"]][4,]),2,mean)
	lh <- if(nrow(IS[["001"]])==6) { 4 } else {3}
	mids[mids$IndicatorString=="001",c("Midpoint.x","Midpoint.y")] <- apply(rbind(IS[["011"]][3,],IS[["001"]][lh,]),2,mean)
	mids[mids$IndicatorString=="010",c("Midpoint.x","Midpoint.y")] <- apply(rbind(IS[["010"]][2,],IS[["111"]][2,]),2,mean)
	mids <- subset(mids,,-IndicatorString)
	mids 
})
setMethod("PlotSetBoundaries","ThreeSquares",function(object,gp){
	grid.polygon(x=object@SetShapes[[1]][,1],y=object@SetShapes[[1]][,2],default.units="native",gp=gp[1])
	grid.polygon(x=object@SetShapes[[2]][,1],y=object@SetShapes[[2]][,2],default.units="native",gp=gp[2])
	grid.polygon(x=object@SetShapes[[3]][,1],y=object@SetShapes[[3]][,2],default.units="native",gp=gp[3])
})

setMethod("PlotUniverse","ThreeSquares",function(object,gp){
	if(missing(gp)) { gp <- NULL }
	outest <- object@IntersectionShapes[["000"]]
	grid.polygon(x=outest[,1],y=outest[,2],default.units="native",gp=gp)
})


setMethod("SetLabels","ThreeSquares",function(object){
	# pull out the outer triangle vertices
	VLabels <- data.frame(Label=rep("unset",3),x=NA,y=NA,hjust=I("center"),vjust=I("center"))
	VLabels[1,c("x","y")] <- object@SetShapes[[1]][2,]
	VLabels[2,c("x","y")] <- object@IntersectionShapes[["011"]][3,]
	VLabels[3,c("x","y")] <- object@SetShapes[[3]][2,]
	VLabels$vjust <- c("bottom","top","bottom")
	VLabels$hjust <- c("right","right","right")

	VLabels
})


setMethod("Areas","ThreeSquares",function(object){
	Areas <- sapply(object@IntersectionShapes,.polygon.area)
	dark.matter <- paste(rep("0",NumberOfSets(object)),collapse="")
	Areas[dark.matter] <- Areas[dark.matter]-sum(Areas[names(Areas)!=dark.matter])
	Areas
})

#######################
# Four squares
setClass("FourSquares",
	 representation(
	"DrawnVenn","Venn"))


compute.S4 <- function(V,doWeights=FALSE,s=.25,likeSquares=TRUE) {
	if (doWeights) { warning("Cant do a weighted S4") }
	if (NumberOfSets(V) != 4) { stop("fournotfour")}
	diagram <- makeAWFE4()
	nodeDataDefaults(diagram,"xy") <- c(NA,NA)
	top.2 <- 2+s
	px <- c(0,0,0,0, 0, 0,2,1,-1,-2,1, 1,-1,-1)
	py <- c(top.2,2,1,0,-1,-2,0,0, 0, 0,1,-1,-1, 1)

	nn  <- paste("p",1:14,sep="")
	pxyl <-  lapply(seq_along(px),function(i) c(px[i],py[i]))

	nodeData(diagram,nn,"xy") <- pxyl      

	edgeDataDefaults(diagram,"arc") <- new("edgeLines")
	# (nearly) all the edges are straight lines
	ee <- edges(diagram)
	for (from in names(ee)) {
		for (to in ee[[from]])  {
			eSet <- edgeData(diagram,from,to,"Set")[[1]]
			exy <- do.call(rbind,nodeData(diagram,c(from,to),"xy"))
			edgeData(diagram,from,to,"arc") <- new("edgeLines",Set=eSet,xy=exy)
		}
	}
	edgeData(diagram,"p6","p1","arc") <- new("edgeLines",Set=1,xy=matrix(
		c(0,-2,   0,-2-s,   -2 - 2*s,-2-s,   -2-2*s, top.2+s,0,top.2+s,0,top.2),ncol=2,byrow=TRUE))
	edgeData(diagram,"p1","p7","arc") <- new("edgeLines",Set=2,xy=matrix(
		c(0,top.2, 2+s,top.2, 2+s,0, 2,0),ncol=2,byrow=TRUE))
	edgeData(diagram,"p10","p1","arc") <- new("edgeLines",Set=2,xy=matrix(
		c(-2,0, -2-s,0, -2-s,top.2, 0,top.2),ncol=2,byrow=TRUE))
	diag.to.cityblock <- function(xy,first="x") {
		block.point <- if(first=="x") { c(xy[2,1],xy[1,2])} else { c(xy[1,1],xy[2,2]) }
		xy <- rbind(xy[1,], block.point,		xy[2,])
	}
	edge.diag.to.cityblock <- function(diagram,from,to,first) {
		object <- edgeData(diagram,from,to,"arc")[[1]]
		object@xy <- diag.to.cityblock(object@xy,first)
		edgeData(diagram,from,to,"arc") <- list(object)
		diagram
	}

	if (likeSquares) {
		diagram <- edge.diag.to.cityblock(diagram,"p2","p11","x")
		diagram <- edge.diag.to.cityblock(diagram,"p11","p7","x")
		diagram <- edge.diag.to.cityblock(diagram,"p7","p12","y")
		diagram <- edge.diag.to.cityblock(diagram,"p12","p6","y")
		diagram <- edge.diag.to.cityblock(diagram,"p6","p13","x")
		diagram <- edge.diag.to.cityblock(diagram,"p13","p10","x")
		diagram <- edge.diag.to.cityblock(diagram,"p10","p14","y")
		diagram <- edge.diag.to.cityblock(diagram,"p14","p2","y")

	}
	
	# we might want to make 1 and 2 not overlap, 
	region <- makeAWFE4star ()
	nodeDataDefaults(region,attr="xy") <- matrix(NA,nrow=0,ncol=2)

	for (face.name in nodes(region)) {
	#	cat(face.name,"\n")
		face <- nodeData(region,face.name,attr="Region")[[1]]
		edgeDataDefaults(face,"arc" ) <- new("edgeLines")
		for (from in nodes(face)) {
			for (to in edges(face,from)[[1]]) {
				edgeData(face,from,to,"arc") <- edgeData(diagram,from,to,"arc")
			}
		}
		face <- rectifyEdges(face)

		nodeData(region,face.name,attr="Region") <- list(face)
		if (face.name != "0000") {
			nodeData(region,face.name,attr="xy") <- list(face.centroid(face)) 
		}
	
	}
	nodeData(region,"1000",attr="xy") <- list(c(-1.5,-1.5))
	nodeData(region,"1100",attr="xy") <- list(c(-1.5, 1.5))
	nodeData(region,"0100",attr="xy") <- list(c( 1.5, 1.5))
	nodeData(region,"0000",attr="xy") <- list(c( 1.5,-1.5))

	drawnVenn <- new("DrawnVenn",regions=region,diagram=diagram)
	new("FourSquares",drawnVenn,V)
}	

face.centroid <- function(face) {
		# given straight line edges, the centroid is the mean of the 
		# vertices
		# pull out the start and end position of every arc
		ed <- lapply(edgeData(face,attr="arc"),function(edgeline){edgeline@xy})
		ed <- do.call(rbind,ed)
		# that has starts and ends all mixed up
		ed <- unique(ed)
		
		apply(ed,2,mean)
}
rectifyEdges <- function(face) {
	# better be a cycle when seen undirected

	if (numNodes(face)==0) { return(face)}
	rec <- face
	top.node <- nodes(rec)[1]
	this.node <- top.node
	
	attrnames <- names(edgeDataDefaults(rec))
	repeat {
		next.nodes <- edges(rec)[[this.node]]
		if (length(next.nodes)==0) {
			next.in.nodes <- inEdges(this.node,rec)[[1]]
			if (length(next.in.nodes)==0) {stop("Nowhere to go") }
			next.node <- next.in.nodes[1]
			rec <- addEdge(this.node,next.node,rec)
			for (attr in setdiff(attrnames,"arc")) {
				edgeData(rec,this.node,next.node,attr ) <- edgeData(rec,next.node,this.node,attr)  # 
			}
			if ("arc" %in% attrnames) {
				edgeData(rec,this.node,next.node,"arc") <- lapply(edgeData(rec,next.node,this.node,"arc"),reverseArc)
			}
			rec <- removeEdge(next.node,this.node,rec)	
		} else {
			next.node <- next.nodes[1]
		}
		this.node <- next.node
		if (this.node==top.node) { break }
	}

	rec
}
	
setMethod("SetLabels","FourSquares",function(object){
	S1arc <-  edgeData(object@diagram,"p6","p1",attr="arc")[[1]]
	S1top <- max(S1arc@xy[,2])
	S2arc <-  edgeData(object@diagram,"p1","p7",attr="arc")[[1]]
	S2top <- max(S2arc@xy[,2]);S2right <- max(S2arc@xy[,1])
	data.frame(Label="unset",x=c(0,S2right,2,1),y=c(S1top,S2top,-1,-2),hjust=c("right","right","right","right"),vjust=c("top","top","bottom","bottom"),
		stringsAsFactors=FALSE)
})
setMethod("PlotUniverse","FourSquares",function(object,gp){
	ur <- UniverseRange(object)
	grid.rect(x=mean(ur$x),y=mean(ur$y),width=diff(ur$x),height=diff(ur$y),default.units="native")
})

setMethod("UniverseRange","FourSquares",function(object) {
	vr <- VisibleRange(object)
	vr
})
setMethod("VisibleRange","FourSquares",function(object) {
	S1arc <-  edgeData(object@diagram,"p6","p1",attr="arc")[[1]]
	yrange <- range(S1arc@xy[,2])
	S2arc <-  edgeData(object@diagram,"p1","p7",attr="arc")[[1]]
	xrange <- range(c(S1arc@xy[,1],S2arc@xy[,2]))
	list(x=xrange,y=yrange)
})
