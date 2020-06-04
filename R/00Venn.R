warning("Entering 00Venn")
# need to integrate DrawnVenn and AWFE classes
#**************************************************************************
#
# # $Id: 00Venn.R,v 1.18 2007/04/23 22:01:46 js229 Exp $

# $Log: 00Venn.R,v $
# Revision 1.18  2007/04/23 22:01:46  js229
# swap funcdef
#
# Revision 1.17  2007/03/20 18:50:17  js229
# no message
#
# Revision 1.15  2007/03/18 12:10:07  js229
# prepare for area check
#
# Revision 1.13  2007/03/17 19:06:16  js229
# tidy up vignette and split triangles out
#
# Revision 1.12  2007/03/16 00:08:08  js229
# mainly improvements to the vignette
#
# Revision 1.10  2007/03/08 23:53:12  js229
# post us and first AWFE implementation
#
# Revision 1.9  2007/02/23 08:15:38  js229
# awfe nearly there
#
# Revision 1.8  2007/02/17 20:34:44  js229
# four squares done; three ellipses started
#
# Revision 1.6  2007/02/13 23:09:37  js229
# unbroken but faces not working
#
# Revision 1.4  2007/02/09 19:43:41  js229
# fiddle with circles
#
# Revision 1.3  2007/02/08 21:44:59  js229
# doEuler; big tidy; rearrange vignette
#
# Revision 1.1  2007/02/08 14:59:57  js229
# working package
#
# Revision 1.1.1.1  2007/02/08 13:51:02  js229
# no message
#

setClass("Venn",
	representation( IndicatorWeight="matrix",IntersectionSets="character")
)

Venn <- function(Indicator,Weight,SetNames,NumberOfSets) {
	if(missing(NumberOfSets)) {
		NumberOfSets <- 0
		if (!missing(Indicator)) {
			NumberOfSets <- ncol(Indicator)
		}
		if (!missing(SetNames)) {
			NumberOfSets <- length(SetNames)
		}
	}
		
	if (missing(Indicator)) {
		if(NumberOfSets==0) {
			Indicator <- matrix(nrow=0,ncol=0)
		} else {
			Indicator <- (data.matrix(do.call(expand.grid,lapply(seq(1,length=NumberOfSets),function(x){c(0,1)}))))==1
		}
	}
	if(!missing(SetNames)) {
		if (!is.null(SetNames)) {
			colnames(Indicator) <- SetNames
		}
	}
	if(missing(Weight)) {
		Weight <- rep(1,nrow(Indicator))
	} else if( length(Weight) >0 & length(Weight) != nrow(Indicator)) {
		stop("Weight length does not match number of intersections")
	}
	
	IndicatorWeight <- cbind(Indicator,.Weight=Weight)

	new("Venn",IndicatorWeight =IndicatorWeight )
}

VennFromSets <- function(setList,collapse="") {
	SetNames <- names(setList) 
	if (is.null(SetNames)) { SetNames <- seq_along(setList) }
	VN <- Venn(SetNames=SetNames)
	VIsig <- IndicatorString(VN)
	IntersectionSets <- sapply(1:length(VIsig),function(x)NULL)
	names(IntersectionSets) <- VIsig
	universe <- NULL
	for ( iset in setList) { universe <- union (universe,iset) }
	for (element in universe) {
		sig <- as.numeric(!is.na(sapply(setList,match,x=element)))
		sig <- paste(sig,collapse="")
		IntersectionSets[[sig]] <- union(IntersectionSets[[sig]],element)
		}
	Weights <- sapply(IntersectionSets,length)
	
	Vres <- Venn(SetNames=SetNames,Weight=Weights)
	IntersectionSets <- sapply(IntersectionSets,paste,collapse=collapse)
	Vres@IntersectionSets <- IntersectionSets
	Vres
}



setGeneric("NumberOfSets",function(object){standardGeneric("NumberOfSets")})
setMethod("NumberOfSets","Venn",function(object){ncol(object@IndicatorWeight)-1})
Indicator <- function(object){
	object <- as(object,"Venn")
	object@IndicatorWeight[,-ncol(object@IndicatorWeight),drop=FALSE]}
setGeneric("SetNames",function(object){standardGeneric("SetNames")})
IndicatorString <- function(object){
	ind <- Indicator(object)
	inn <- apply(ind,1,paste,collapse="")
	inn
}
Weights <- function(object) {
	V <- as(object,"Venn")
	wght <- V@IndicatorWeight[,".Weight"]
	names(wght) <- IndicatorString(V)
	wght
}

"Weights<-" <- function(object,value) {
	V <- as(object,"Venn")
	VS <- IndicatorString(V)
	if (is.null(names(value))) {
		names(value) <- VS
	}
	value <- value[match(VS,names(value))]
	object@IndicatorWeight[,".Weight"] <- as.numeric(value)
	object
}

dark.matter.signature <- function(object) {
	V <- as(object,"Venn")
	VS <- IndicatorString(V)
	VS[regexpr("1",VS)<0]
}	


setMethod("SetNames","Venn", function(object) {cn <- colnames(object@IndicatorWeight);cn[cn!=".Weight"] })
setMethod("[","Venn", function(x,i,j,...,drop) {
	if (!missing(i)) {
		stop("Can't subset on rows")
	}
	if (!missing(j)) {
		Indicator <- x@IndicatorWeight
		newIndicator <- aggregate(Indicator[,".Weight",drop=FALSE],
			by=data.frame(Indicator[,j,drop=FALSE]),FUN=sum)
		for (col in setdiff(colnames(newIndicator ),".Weight")) {
			newIndicator [,col] <- as.numeric(as.character(newIndicator [,col])) # sigh
		}
		x@IndicatorWeight <- data.matrix(newIndicator)	

		IntersectionSets <- x@IntersectionSets
		if (length(IntersectionSets )>0) {
			newIntersectionSets <- aggregate(list(IntersectionSets=IntersectionSets),data.frame(Indicator[,j,drop=FALSE]),FUN=paste,collapse="")
			x@IntersectionSets <- as.character(newIntersectionSets$IntersectionSets)
		}
	}
	x
	})

.WeightVisible <- function(V) {# excluding  elements not in any sets
	wght <- V@IndicatorWeight
	wghtsum <- apply(wght[,-ncol(wght)],1,sum)
	sum(wght[wghtsum>0,ncol(wght)]) 
}
.WeightUniverse <- function(V) {
	wght <- Weights(V)
	sum(wght) 
}


################################################
setGeneric("PlotSetBoundaries",function(object,gp){standardGeneric("PlotSetBoundaries")})
setGeneric("PlotFaces",function(object,gp,arrow){standardGeneric("PlotFaces")})
setGeneric("PlotNodes",function(object,gp){standardGeneric("PlotNodes")})
setGeneric("PlotUniverse",function(object,gp){standardGeneric("PlotUniverse")})
setGeneric("IntersectionMidpoints",function(object){standardGeneric("IntersectionMidpoints")}) 
setGeneric("SetLabels",function(object){standardGeneric("SetLabels")}) 
setGeneric("Areas",function(object){standardGeneric("Areas")}) 
setGeneric("ComputeAreas",function(object,nintervals){standardGeneric("ComputeAreas")}) 
setGeneric("VisibleRange",function(object){standardGeneric("VisibleRange")}) 
setGeneric("UniverseRange",function(object){standardGeneric("UniverseRange")}) 


setClass("DrawnVenn",representation(diagram="graph",regions="graph"),
	prototype=list(diagram=new("graphNEL"),regions=new("graphNEL")))

setClass("edgeSector",representation(Set="numeric",centre="numeric",hand="numeric"))
setClass("edgeLines",representation(Set="numeric",xy="matrix"))
setClass("edgeFunction",representation(Set="numeric",s="numeric",f="function"))
setClassUnion("edgeDrawn",c("edgeSector","edgeLines","edgeFunction"))

setGeneric("Arc.to.xy",function(diagram,from,to,ftarc,nintervals){standardGeneric("Arc.to.xy")})
setGeneric("reverseArc",function(object){standardGeneric("reverseArc")})
setGeneric("plotArc",function(object,gp){standardGeneric("plotArc")})

setMethod("reverseArc","edgeLines",function(object){
	object@xy <- object@xy[ rev(seq_len(nrow(object@xy))),]
	object
})
setMethod("plotArc","edgeLines",function(object,gp){
	Set <- object@Set
	gp <- gp[Set]
	grid.lines(x=object@xy[,1],y=object@xy[,2],gp=gp,default.units="native")
})

setMethod("PlotSetBoundaries","DrawnVenn",function(object,gp){
	if (missing(gp)) {
		gp <- gpar(col=trellis.par.get("superpose.symbol")$col)
	}
	dv <- as(object,"DrawnVenn")
	diagram <- dv@diagram
	edgeArcs <- edgeData(diagram,attr="arc")
	lapply(edgeArcs,plotArc,gp=gp)
})

setMethod("PlotNodes","DrawnVenn",function(object,gp){
	dv <- as(object,"DrawnVenn")
	diagram <- dv@diagram
	pxy <- nodeData(diagram,attr="xy")
	xy <- do.call(rbind,pxy)
	grid.text(x=xy[,1],y=xy[,2],label=names(pxy),default.units="native",
		just=c("left","bottom"))
})


setMethod("IntersectionMidpoints","DrawnVenn",function(object){
	dv <- as(object,"DrawnVenn")
	regions <- dv@regions
	cxy <- nodeData(regions,attr="xy")
	df <- data.frame(do.call(rbind,cxy))
	colnames(df) <- c("Midpoint.x","Midpoint.y")
	df$VS <- nodes(regions)
#	df <- df[match(df$VS,IndicatorString(as(object,"Venn"))),]
	df$hjust <- "center"
	df$vjust <- "center"

	df
})

sector.to.xy <-  function(diagram,from,to,ftarc,nintervals=100) {
	# draw the segment of a circle from point from to point to
	p.xy <- nodeData(diagram,c(from,to),"xy")
	centre <- ftarc@centre
	p.xy.0 <- sapply(p.xy,function(xy){xy-centre})
	
	# they better both be on the circle
	rads <- as.numeric(apply(p.xy.0,2,function(x)sqrt(sum(x^2))))
	stopifnot ( all.equal(rads,rep(rads[1],length(rads))))
	r <- rads[1]
	hand <- ftarc@hand

	if (from==to | isTRUE(all.equal(p.xy[[1]],p.xy[[2]]))) { 
		if (hand>0) { thetafrom <- 0; thetato <- 2* pi } else { thetafrom <- 2*pi; thetato <- 0}
	} else {
		thetafrom <- -atan2(p.xy.0[2,1],p.xy.0[1,1])
		thetato <- -atan2(p.xy.0[2,2],p.xy.0[1,2])
		# if hand> 1 we always go clockwise, increasing thetafrom ( reverse of direction atan2 goes)
		# so thetato must be bigger than thetafrom
		if (hand>0) {	if (thetato<thetafrom)  { thetato <- thetato + 2* pi } }
		if (hand<0) {	if (thetato>thetafrom)  { thetato <- thetato - 2* pi } }
	}
	theta <- seq(from=thetafrom,to=thetato,length=nintervals)
	x <- r * cos(theta)+centre[1]; y <- r* sin(theta)+centre[2]
	cbind(x,y)
}

setMethod("Arc.to.xy",c("graph","character","character","edgeSector","numeric"),
	sector.to.xy)
setMethod("Arc.to.xy",c("graph","character","character","edgeLines","numeric"),
	function(diagram,from,to,ftarc,nintervals) {})


.universe.range <- function (object) 
{
    list(x = object@universe.x, y = object@universe.y)
}
.plotuniverse.square <- function(object,gp) {
	if(missing(gp)) { gp <- NULL }
		grid.rect(x=mean(object@universe.x),y=mean(object@universe.y),
		width=diff(object@universe.x),height=diff(object@universe.y),default.unit="native",gp=gp)
	}

.square.universe <- function(object,doWeights,smudge=0.05) {
	if (doWeights) {
		# minimal square box 
		minimal.square.universe.area <- diff(VisibleRange(object)$x)*diff(VisibleRange(object)$y)
		V <- as(object,"Venn")
		visible.area <- .WeightVisible(V)
		dark.matter.area <- .WeightUniverse(V) - .WeightVisible(V)
		dark.matter.scale.squared <- (dark.matter.area + visible.area)/minimal.square.universe.area
		if (dark.matter.scale.squared  < 1 + smudge) {
	#		warning("Square box is too large for the weight")
			dark.matter.scale.squared = 1 + smudge
		}
	} else {
		dark.matter.scale.squared = 1.2
	}

	object@universe.x <- VisibleRange(object)$x * sqrt(dark.matter.scale.squared)
	object@universe.y <- VisibleRange(object)$y * sqrt(dark.matter.scale.squared)
	
	object
}



####################

compute.Venn <- function(V,doWeights=TRUE,doEuler=FALSE,type="circles",nintervals) {
	nSets <- NumberOfSets(V)
	if (nSets < 2) {
		stop("Not enough sets")
	} 
	if (type=="AWFE"  ) {
		C3 <- compute.AWFE(V)
	} else if (type=="ChowRuskey") { 
		C3 <- compute.CR(V,doWeights)
	} else if (nSets==2) {
		if (type=="circles") {
			C3 <- compute.C2(V,doWeights,doEuler)
		} else if (type=="squares") {
			C3 <- compute.S2(V,doWeights,doEuler)
		}
	} else if (nSets==3){
		if (type=="circles") {
			C3 <- compute.C3(V,doWeights,nintervals)
	} else if (type=="squares") {
			C3 <- compute.S3(V,doWeights)
		} else if (type=="triangles") {
			C3 <- compute.T3(V,doWeights)
		}
	} else if (nSets==4) {
		if (type=="ellipses") {
			C3 <- compute.E4(V,doWeights)
		} else {
			C3 <- compute.S4(V,doWeights)
		}
	} else {
		stop(sprintf("Unimplemented type %s for %d sets\n",type,nSets))
	}
	C3
}

plot.Venn <- function(V,doWeights=TRUE,doEuler=FALSE,type="circles",add=FALSE,show=list(intersection.weight=TRUE),
		gp.V.set.colour){
	C3 <- compute.Venn(V,doWeights=doWeights,doEuler=doEuler,type=type)
	if (!add) {
		grid.newpage()
	}
	PlotVennGeometry(C3,gp.V.set.colour=gp.V.set.colour,show=show)
}

PlotVennGeometry <- function(C3,gp.V.set.colour,show=list(intersection.weight=TRUE)) {
	show.default <- list(universe=TRUE,Sets=TRUE,SetLabels=TRUE,
		dark.matter=TRUE,
		Faces=FALSE,
		intersection.weight=TRUE,intersection.elements=FALSE,indicator.string=FALSE)
	unshown <- names(show)[! names(show) %in% names(show.default)]
	if (length(unshown)>0) {
		warning(sprintf("Unknown show parameters %s",paste(unshown,collapse=",")))
	}
	show <- show[names(show) %in% names(show.default)]

	show.default[names(show)] <- show
	if (missing(gp.V.set.colour)) {
		tpg <- trellis.par.get("superpose.polygon")
		gp.V.set.colour <- gpar(col=tpg$col)
	}
	
	CreateViewport(C3)
	
	if(show.default$universe) {
		PlotUniverse(C3)
	}
	if (show.default$Faces) {
		if (is(C3,"DrawnVenn")|is(C3,"AWFE")) {
			PlotFaces(C3)
		} else {
			warning("Can't plot faces for this object")
		}
	}
	if(show.default$Sets) {
		PlotSetBoundaries(C3,gp=gp.V.set.colour)
	}
	if(show.default$SetLabels) {
		PlotSetLabels (C3,gp=gp.V.set.colour) 
	}

	if(show.default$intersection.weight) {
		PlotIntersectionText(C3,element.plot="weight",show.dark.matter=show.default$dark.matter)	
	}
	if(show.default$intersection.elements) {
		PlotIntersectionText(C3,element.plot="elements",show.dark.matter=show.default$dark.matter)	
	}
	if(show.default$indicator.string) {
		PlotIntersectionText(C3,element.plot="inn",show.dark.matter=show.default$dark.matter)	
	}


	UpViewports()	
}

CreateViewport <- function(object) {
	xData <- UniverseRange(object)$x
	yData <- UniverseRange(object)$y
	makevp.eqsc(xData,yData)
}
UpViewports <- function() {
	upViewport()
	upViewport()
	upViewport()
}
makevp.eqsc <- function(xrange,yrange) {
	# cf Fig 7.4 of Murrell R Graphics
	pushViewport(plotViewport(name="Vennmar",c(1,1,1,1)))
	pushViewport(viewport(name="Vennlay",layout=grid.layout(1,1,widths=diff(xrange),heights=diff(yrange),respect=TRUE)))
	pushViewport(viewport(name="Vennvp",layout.pos.row=1,layout.pos.col=1,xscale=xrange,yscale=yrange))
	}

PlotIntersectionText <- function(object,element.plot="weight",show.dark.matter=TRUE) {
	V <- as(object,"Venn")
	nSets <- NumberOfSets(V)
	VI <- IntersectionMidpoints(object);
	if( element.plot=="weight" ) {
		VINames <- Weights(V)
	} else if (element.plot=="inn" ) {
		VINames <-  IndicatorString(V)
	} else if (element.plot=="elements") {
		elements <- V@IntersectionSets
		if (is.null(elements)) {
			warning("No intersection sets elements known")
			VINames <- ""
		} else {
			VINames <- unlist(elements) # only work for scalar elements
		}
	}
	VI$Names <- VINames
	if (!show.dark.matter) {
		is.dark <- IndicatorString(V)==dark.matter.signature(V)
		VI <- VI[!is.dark,]
	}
	if (!"hjust" %in% colnames(VI)) { VI$hjust <- "centre" }
	if (!"vjust" %in% colnames(VI)) { VI$vjust <- "centre" }
	hj <-sapply( VI$hjust,function(EXPR){switch(EXPR,left=0,right=1,center=,centre=0.5)})
	vj <-sapply( VI$vjust,function(EXPR){switch(EXPR,top=1,bottom=0,center=,centre=0.5)})
	for (ij in 1:nrow(VI)) {
		grid.text(x=VI$Midpoint.x[ij],y=VI$Midpoint.y[ij],hjust=hj[ij],
		vjust=vj[ij],label=VI$Names[ij],default.unit="native")
	}

	#grid.text(x=VI$Midpoint.x,y=VI$Midpoint.y,label=VINames,default.unit="native") 
}


PlotSetLabels <- function(object,gp) {
	VLabels <- SetLabels(object)
	if(nrow(VLabels)==0){ warning("Can't show Set labels"); return()}
#	print(VLabels)
	VLabels$Label <- SetNames(as(object,"Venn"))
	# just may not be vectorised...
	hj <-sapply( VLabels$hjust,function(EXPR){switch(EXPR,left=0,right=1,center=,centre=0.5)})
	vj <-sapply( VLabels$vjust,function(EXPR){switch(EXPR,top=1,bottom=0,center=,centre=0.5)})

	for (ij in 1:nrow(VLabels)) {
		grid.text(x=VLabels$x[ij],y=VLabels$y[ij],hjust=hj[ij],
		vjust=vj[ij],gp=gp[ij],label=as.character(VLabels$Label[ij]),default.unit="native")
	}
}



sector.to.xy <-  function(diagram,from,to,ftarc,nintervals=100) {
	# draw the segment of a circle from point from to point to
	p.xy <- nodeData(diagram,c(from,to),"xy")
	centre <- ftarc$centre
	p.xy.0 <- sapply(p.xy,function(xy){xy-centre})
	
	# they better both be on the circle
	rads <- as.numeric(apply(p.xy.0,2,function(x)sqrt(sum(x^2))))
	stopifnot ( all.equal(rads,rep(rads[1],length(rads))))
	r <- rads[1]
	hand <- ftarc$hand

	if (from==to | isTRUE(all.equal(p.xy[[1]],p.xy[[2]]))) { 
		if (hand>0) { thetafrom <- 0; thetato <- 2* pi } else { thetafrom <- 2*pi; thetato <- 0}
	} else {
		thetafrom <- -atan2(p.xy.0[2,1],p.xy.0[1,1])
		thetato <- -atan2(p.xy.0[2,2],p.xy.0[1,2])
		# if hand> 1 we always go clockwise, increasing thetafrom ( reverse of direction atan2 goes)
		# so thetato must be bigger than thetafrom
		if (hand>0) {	if (thetato<thetafrom)  { thetato <- thetato + 2* pi } }
		if (hand<0) {	if (thetato>thetafrom)  { thetato <- thetato - 2* pi } }
	}
	theta <- seq(from=thetafrom,to=thetato,length=nintervals)
	x <- r * cos(theta)+centre[1]; y <- r* sin(theta)+centre[2]
	cbind(x,y)
}

numEdges <- function(G) {
	if(numNodes(G)==0) return(0)
	sum(sapply(edges(G),length))
}
face.to.xy <- function(diagram,face,nintervals=100) {
	if (numEdges(face)==0) return(list())
	edgeCount <- sapply(edges(face),length)
	free.nodes <- names(edgeCount)[edgeCount==0]
	if(length(free.nodes)>0) {face <- removeNode(free.nodes,face) }
	first.node <- nodes(face)[1]
	to.first <- inEdges(first.node,face)
	cut.face <- removeEdge(to.first[[1]],first.node,face)
	face.nodes <- dfs(cut.face)$discovered
	lapply(face.nodes,function(from){
		to <- edges(face)[[from]]
		stopifnot(length(to)==1)
		to <- to[1]
		arc <- edgeData(face,from,to,"arc")[[1]]
		arc.Set <- arc@Set
		arc.xy <- Arc.to.xy(diagram,from,to,arc,nintervals)
		list(Set=arc.Set,xy=arc.xy)
	})
}


.PlotFaces.drawnVenn  <- function(object,gp,arrow){
	if (missing(gp)){ gp <- gpar()}
	if ( !all(c("diagram","regions") %in% slotNames(object))
		| numNodes(object@regions)==0) {
		warning("Cant plot faces for this object")
		return()
	}
	diagram <- object@diagram
	regions <- object@regions
	faces=nodeData(regions,attr="Region")
	for (face.name in names(faces)) {
		PlotFace(face.name,faces,diagram,gps=gp,gpl=gp)
	}
}
setMethod("PlotFaces","DrawnVenn",.PlotFaces.drawnVenn)



PlotFace <- function(face.name,faces,diagram,nintervals=100,gps=gpar(),gpl=gpar()) {
#cat(face.name,"\n")
	if (is.null(gpl$col)){ 
		gpl$col <- trellis.par.get("superpose.symbol")$col
	}
	if (is.null(gpl$lty)) {
		gpl$lty <- 0
	} # region boundaries usually put in as sets
	fix <- match(face.name,names(faces))
	if (is.null(gps$fill)){
		gps$fill <- trellis.par.get("superpose.polygon")$col[fix]
	}	
	if (is.null(gps$lty)){
		gps$lty <- 0
	}	
	face <- faces[[face.name]]
	face.Sxy <- face.to.xy(diagram,face,nintervals=nintervals)
	
	all.xy <- do.call(rbind,lapply(face.Sxy,function(x)x$xy))
	if (!is.null(all.xy)) {
		grid.polygon(all.xy[,1],all.xy[,2],default.units="native",gp=gps)
	}
	lapply(face.Sxy,function(x){
		if (length(x)>0) {
		Set=x$Set
		xy <- x$xy
		gpss <- gps
		gpss$col <- gps$col[Set]
		grid.lines(xy[,1],xy[,2],
			default.units="native",gp=gpl)
		}
	})
}


.Areas.DrawnVenn  <- function(object){
	if ( !all(c("regions") %in% slotNames(object))
		| numNodes(object@regions)==0) {
		warning("Cant plot faces for this object")
		return()
	}
	diagram <- object@diagram
	regions <- object@regions
	faces=nodeData(regions,attr="Region")
	face.names <- as.vector(names(faces))
	res <- sapply(face.names, function(face.name) {
		.Area.face(face.name,faces,diagram)
	})
	res
}
setMethod("Areas","DrawnVenn",.Areas.DrawnVenn )


.polygon.area <- function(Shape) {
	if (is.null(nrow(Shape))) return(NA)
	if (nrow(Shape)<3) {
		return(0)
	}
	xy1 <- Shape; xy2 <- Shape[ c(2:nrow(Shape),1),]
	x1 <- xy1[,1];y1 <-xy1[,2];x2<-xy2[,1];y2<- xy2[,2]
	det <- x1*y2 - x2*y1
	abs(sum(det)/2)
}
# this function also used by ChowRuskey.R at least

.Area.face<- function(face.name,faces,diagram,nintervals=100) {
	face <- faces[[face.name]]
	face.Sxy <- face.to.xy(diagram,face,nintervals=nintervals)
	all.xy <- do.call(rbind,lapply(face.Sxy,function(x)x$xy))
	.polygon.area(all.xy)
}


