warning("Entering Circles.R")

################################################
# some geometries 
# r is radius of circle
# d distance from origin

CircleRepresentation <- representation(
	r="numeric",x="numeric",y="numeric",type="character",	areas="numeric",
	centroids.x="numeric",centroids.y="numeric",
	universe.x="numeric",universe.y="numeric",
	nintervals="numeric","DrawnVenn","Venn")
setClass("ThreeCircles",CircleRepresentation )
setClass("TwoCircles",CircleRepresentation )

.merge.graphs <- function(G1,G2) {
	n1 <- nodes(G1)
	n2 <- nodes(G2)
	n12 <- setdiff(n1,n2)
	if (length(n12)>0) {	G2 <- addNodes(n12,G2)  }
	n21 <- setdiff(n2,n1)
	if (length(n21)>0) {	G1 <- addNodes(n21,G1)  }
	G<- union(G1,G2)
}


TwoCircles <- function(r,d,V) {
	if (length(r) !=2 ) {
		if (length(r)==1 ) {
			r <- rep(r,2)
		} else {
			stop("Need two circle radii")
		}
	}
	# Two circles, radius r1 and r2, distance d apart
	d.origin  <-  d/2
	VS <- IndicatorString(V)
	diagram <- new("graphNEL",nodes=c("p1","p2"),edgemode="directed")
	nodeDataDefaults(diagram,"xy") <- c(NA,NA)
	edgeDataDefaults(diagram,"arc") <- new("edgeSector")

	arcs.1 <- new("edgeSector",Set=1,centre=c(-d.origin,0),hand=1)
	arcs.2 <- new("edgeSector",Set=2,centre=c(+d.origin,0),hand=1)
	reverse.arc <- function(arc) { arc@hand <- - arc@hand; arc }


	x <- (d.origin) * c(-1,1)
	y <- c(0,0)
	c10 <- NA;c11 <- NA; c01 <- NA
	areas <- rep(0,length(VS))

	haveIntersection <- (d>abs(r[1]-r[2]) & d < r[1]+r[2])  # relies on numerics...
	r01 <- diagram
	r11 <- diagram
	r10 <- diagram
	aed <- function(face,n1,n2,arcnum) {
	 p1 <- paste("p",n1,sep=""); p2 <- paste("p",n2,sep="");
	 face <- addEdge(p1,p2,face)
	 arc <- if(abs(arcnum)==1){arcs.1}else{arcs.2}
	 if (arcnum<0) { arc <- reverse.arc(arc)}
	 edgeData(face,p1,p2,"arc") <- list(arc)
	 face
	}
	if (haveIntersection) {
		# distance from each circle centre to the chord of the intersection
		d1 <- (d^2 - r[2]^2+ r[1]^2) /( 2* d)
	 	d2 <- d - d1
		# half-height from the chord to the intersection
		yh <- (1/(2*d))* sqrt(4*d^2*r[1]^2-(d^2-r[2]^2+r[1]^2)^2)
		p1xy <- c(-d.origin+d1,yh)
		p2xy <- c(-d.origin+d1,-yh)

		r10 <- aed(r10,2,1,1); r10 <- aed(r10,1,2,-2)
		r01 <- aed(r01,1,2,2); r01 <- aed(r01,2,1,-1)
		r11 <- aed(r11,1,2,1); r11 <- aed(r11,2,1,2)

		c10 <-  ( (x[1]-r[1]) + (x[2]-r[2]))/2
		c11 <-  ( (x[2]-r[2]) + (x[1]+r[1]))/2
		c01 <-  ( (x[1]+r[1]) + (x[2]+r[2]))/2

		areas[VS=="10"] <- pi * r[1]^2
		areas[VS=="01"] <- pi * r[2]^2

		# area of a circular segment
		Ap <- function(r,d) {r^2 * acos(d/r) - d * sqrt(r^2-d^2)}
		# triangle heights
		lensarea <- Ap(r[1],d1)+Ap(r[2],d2)
		areas[VS=="10"] <- areas[VS=="10"] - lensarea
		areas[VS=="01"] <- areas[VS=="01"] - lensarea
		areas[VS=="11"] <- lensarea

	} else {	
		if (d==r[1]+r[2]) {
			p1xy <- c(-d.origin+r[1],0)
			p2xy <- p1xy
			r01 <- aed(r01,1,1,1)
			r10 <- aed(r10,2,2,2)
			c10 <-  x[1]
			c01 <-  x[2]
			areas[VS=="10"] <- pi * r[1]^2
			areas[VS=="01"] <- pi * r[2]^2

		} else if ( r[2]>=r[1]) {
			p1xy <- c(-d.origin-r[1],0)
			p2xy <- p1xy
			r01 <- aed(r01,1,2,2);r01 <- aed(r01,2,1,-1)
			r11 <- aed(r11,1,1,1)		
			c01 <- ( (x[1]+r[1]) + (x[2]+r[2]))/2
			c11 <- x[1]
			areas[VS=="01"] <- pi * r[2]^2 -pi * r[1]^2
			areas[VS=="11"] <- pi * r[1]^2
		} else if ( r[1]>=r[2]) {
			p1xy <- c(d.origin+r[2],0)
			p2xy <- p1xy
			r10 <- aed(r10,1,2,1); r10 <- aed(r10,2,1,-2)
			r11 <- aed(r11,2,2,2)	
			c10 <- ( (x[1]-r[1]) + (x[2]-r[2]))/2
			c11 <- x[2]
			areas[VS=="10"] <- pi * r[1]^2 -pi * r[2]^2
			areas[VS=="11"] <- pi * r[2]^2

		} else { stop("Whoops")}
	}
	nodeData(diagram,c("p1","p2"),"xy") <- list(p1xy,p2xy)
	

	region <- makeQn(2)
	nodeDataDefaults(region,"Region") <- new("graphNEL")

	nodeData(region,"10",attr="Region") <- r10
	nodeData(region,"01",attr="Region") <- r01
	nodeData(region,"11",attr="Region") <- r11

		

	centroids.y <- rep(0,length(VS))
	# only approximate..
	centroids.x <- rep(NA,length(VS))
	centroids.x[VS=="10"] <- c10
	centroids.x[VS=="01"] <- c01
	centroids.x[VS=="11"] <- c11


	new("TwoCircles",r=r,x=x,y=y,areas=areas,centroids.x=centroids.x,centroids.y=centroids.y,
		diagram=diagram,regions=region,V)
}


# 2-d diagrams
compute.C2 <- function(V,doWeights=TRUE,doEuler=FALSE) {
	Vcalc <- V
	if(!doWeights) {
		if (!doEuler) {
		# want each area, even if empty
			wght <- Weights(Vcalc )
			wght <- 1+ 0*wght
			Weights(Vcalc ) <- wght
		} else {
			wght <- Weights(Vcalc )
			wght [wght !=0] <- c(3,2,2,1)[wght !=0]
			Weights(Vcalc ) <- wght
		}
	}
	dList <- .Venn.2.weighted.distance (Vcalc,doEuler) # returns radii of two circles and their distance
	r1 <- dList$r1;r2 <- dList$r2; d <- dList$d; 
	C2 <- TwoCircles(r=c(r1,r2),d=d,V) # d in TwoCircles is distance of centre from origin
	C2 <- .square.universe(C2,doWeights)
	
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


#####################################
# given four weights we seek to find two circles all of whose
# interesction areas are proportional to those weights
.Venn.2.weighted.distance <- function(V,doEuler) {
	if (NumberOfSets(V) != 2) {
		stop(sprintf("Wrong number of sets (%d) in .Venn.2.weighted.distance",NumberOfSets(V)))
	}
	# cf Chow and Ruskey, 200?
	Weight <- Weights(V)
	wab <- Weight["00"] #not used 
	wAb <- as.numeric(Weight["10"])
	waB <- as.numeric(Weight["01"])
	wAB <- as.numeric(Weight["11"])

	inneroff <- 1 # set = 2 to put inner circles completely inside
	r1 <- sqrt( (wAb+wAB)/pi)
	r2 <- sqrt( (waB+wAB)/pi) # area proportional to weights

	if (wAb==0) { #inside the other one 
		if (doEuler) {
			d <- (r2-r1)/inneroff
		} else { # bodge it outside to ensure the Venntersection..doesnt really make much sense
			d <- (r2-r1)+0.1*r1
		}
	} else if (waB==0) { # also
		if (doEuler) {
			d <- (r1-r2)/inneroff
		} else { # bodge it outside to ensure the Venntersection..doesnt really make much sense
			d <- (r1-r2)+0.1*r1
		}
	} else if (wAB==0) {
		d <- inneroff *  (r1+r2)  # no intersection
		if (!doEuler) {
			d <- 0.95 * d
		}
		
	} else  { # all nonzero
	alpha  <- function(d,r1,r2) {
		alphaD <- (d^2+r1^2-r2^2)
		alphaD <- ifelse(alphaD ==0,alphaD ,alphaD / (2 * r1 * d) )
		alphaD <- ifelse(alphaD>1,1,alphaD)
		alphaD <- ifelse(alphaD< -1,-1,alphaD)
		alpha <-  2 *acos( alphaD )
		alpha
	}
	sigma <- function(d,r1,r2) { 
		alpha <- alpha(d,r1,r2)
		betaD <- (d^2+r2^2-r1^2); 
		betaD <- ifelse(betaD==0,betaD,betaD/ (2 * r2 * d) )
		betaD <- ifelse(betaD< (-1),-1,betaD )
		betaD <- ifelse(betaD> 1,1,betaD )
		beta <-   2 *acos( betaD )
		sigma <- 0.5 * r1^2*(alpha-sin(alpha))+ 0.5 * r2^2 * (beta-sin(beta))
		sigma
	}
	sigmaminuswAB <- function(d,r1,r2,wAB) {
		sigma(d,r1,r2) - wAB
	}
	d <-	uniroot(sigmaminuswAB ,lower=r1-r2,upper=r1+r2,r1=r1,r2=r2,wAB=wAB)$root
	} # wAB != 0 
	list(r1=r1,r2=r2,d=d )

}



ThreeCircles <- function(r,x,y,d,angles,V,nintervals=100) {
	exact <- FALSE
	if (missing(x) | missing(y)) {
		if (missing(d)) {
			stop("Need x and y or d")
		}
		exact <- (length(d)==1 & length(r)==1 & missing(angles)) # actually also need angles uniform...
		if (missing(angles)) {
			angles <- pi/2-c( 0, 2*pi/3, 4 * pi/3)
		}
		x <- d*cos(angles)
		y <- d*sin(angles)
	}
	if (length(r) !=3 ) {
			if (length(r)==1 ) {
				r <- rep(r,3)
			} else {
				error("Need three circle radii")
			}
		}
	DV<- new("DrawnVenn")
	C3 <- new("ThreeCircles",r=r,	x=x,y=y,DV,V)

	if (exact) {
		r <- r[1]
		intersections<- sqrt((r^2-3*d^2/4))-d/2
		outersector <- ( d+r + intersections)/2
	
		outersections <- sqrt((r^2-3*d^2/4)) +d/2
		innersector <- (r-d+outersections)/2
		VI <- data.frame(Indicator(V),check.names=FALSE)
		VIsum <- apply(VI,1,sum)
		VI$Midpoint.x <- NA
		VI$Midpoint.y <- NA
		VI[VIsum==1,c("Midpoint.x","Midpoint.y")] <- c(x=outersector *cos(angles),y=outersector *sin(angles))
		VI[VIsum==2,c("Midpoint.x","Midpoint.y")] <- c(x=innersector *cos(-angles+2*pi/3),y=innersector *sin(-angles+2*pi/3))
		VI[VIsum==3,c("Midpoint.x","Midpoint.y")] <- c(0,0)
		C3@centroids.x <- VI$Midpoint.x
		C3@centroids.y <- VI$Midpoint.y
		C3@areas <- rep(as.numeric(NA),8);names(C3@areas) <- IndicatorString(V)
	} else {
	
		C3 <- ComputeAreas(C3,nintervals=nintervals)
	}

	C3
}


.pairwise.overlaps <- function(V) {
	VI <- Indicator(V)
	VIsum <- apply(VI,1,sum)
	Vpairs <- unique(VI[VIsum==2,])
	Vwhich <- which(Vpairs==1,arr.ind=TRUE)
	Vindex <- Vpairs
	Vindex [Vwhich] <- Vwhich[,2]
	isdisjoint <- function(vrow) {
		vsub <- V[,vrow]
		vsum <- apply(Indicator(vsub),1,sum)
		overweight <- Weights(vsub)[vsum==2]
		overweight==0
	}
	Vpairs <- data.frame(Vpairs,check.names=FALSE)
	Vpairs$Disjoint <- apply(Vindex,1,isdisjoint)
	Vpairs
}


compute.C3 <- function(V,doWeights=TRUE,nintervals=100) {
	doEuler <- TRUE
if (doWeights) {	
		overlaps <- .pairwise.overlaps(V)
		dList12 <- .Venn.2.weighted.distance (V[,c(1,2)],doEuler ) # returns radii of two circles and their distance
		dList23 <- .Venn.2.weighted.distance (V[,c(2,3)],doEuler ) # 
		dList31 <- .Venn.2.weighted.distance (V[,c(3,1)],doEuler ) #
		
		disjointcount <- sum(overlaps$Disjoint)
	#	cp <- dList12$d; b <- dList23$d; a <- dList31$d
		dp <- c( cp= dList12$d, b = dList23$d, a = dList31$d)
		smidge <- 1 - 1e-4
		if (disjointcount==3) {
			# all disjoint, arrange with centres in equilateral triangle
			dmax <- max(dp)
			dp <- dp*0+dmax
		} else if (disjointcount==2) {
			#one is disjoint from both the others
			# set it at the same distance from both
			# and far enough that it will be a triangle 
			conjoint <- overlaps[!overlaps$Disjoint,-ncol(overlaps)]
			conjointix <-  ( which(conjoint==0)%%3 +1) 
			conjointdistance <-  dp[conjointix ]
			disjointdistances <- dp[-conjointix]
			disjointdistances <- max(disjointdistances,conjointdistance/2)
			dp[-conjointix] <- disjointdistances 
		} else if (disjointcount==1) {
			# only one missing intersection; we set its distance to
			# force a (near) straightline
			# nb this will still fail if the intersections are large enough
			disjoint <- overlaps[overlaps$Disjoint,-ncol(overlaps)]
			disjointix <- ( which(disjoint ==0)%%3 +1) 
			conjointdistances <-  dp[-disjointix ]
		 	dp[disjointix] <- sum(conjointdistances)*smidge
		} 
		cp= dp["cp"]; b = dp["b"]; a=dp["a"]

		# can we satisfy the triangle inequality? if not, bodge it
		if (a+b < cp) {
			cp <-  (a+b)*smidge
		}
		if (b+cp < a) {
			a <-  (b+cp)*smidge
		}
		if (cp+a < b ) {
			b <- (cp+a)*smidge
		}

		# use the SSS rule to compute angles in the triangle
		CP <- acos( (a^2+b^2-cp^2)/(2*a*b))
		B <- acos( (a^2+cp^2-b^2)/(2*a*cp))
		A <- acos( (cp^2+b^2-a^2)/(2*b*cp)) 
		# arbitrarily bisect one and calculate xy offsets from first centre
		theta <- B/2
		o21 <- cp * c( sin(theta), - cos(theta))
		o31 <- a * c( -sin(B-theta),	-cos(B-theta))
		# pick a centre for the first one
		c1 <- c(0, dList12$r1)
		c2 <- c1 + o21
		c3 <- c1 + o31
		C3 <- ThreeCircles(r=c(dList12$r1,dList12$r2,dList31$r1),
			x=c(c1[1],c2[1],c3[1]),y=c(c1[2],c2[2],c3[2]),V=V,nintervals=nintervals)
	} else {
		d <- 0.4
		C3 <- ThreeCircles(r=0.6,d=0.4,V=V)
	}
	C3 <- .square.universe(C3,doWeights=doWeights)
	C3
}


.CircleVisibleRange <- function(object){
	centres.x <- object@x
	centres.y <- object@y
	radius <- object@r
	xrange <- c(min(centres.x-radius),max(centres.x+radius))
	yrange <- c(min(centres.y-radius),max(centres.y+radius))
	list(x=xrange,y=yrange)
}

.universe.range <- function(object){list(x=object@universe.x,y=object@universe.y)}

setMethod("VisibleRange","ThreeCircles",.CircleVisibleRange)
setMethod("VisibleRange","TwoCircles",.CircleVisibleRange)
setMethod("UniverseRange","ThreeCircles",.universe.range)
setMethod("UniverseRange","TwoCircles",.universe.range)


setMethod("PlotNodes","ThreeCircles",function(object,gp) {
		warning("PlotNodes not implemented for ThreeCircles")
})

.plotboundaries.circles <- function(object,gp) {
	if(missing(gp)) { gp <- NULL }
		grid.circle(x=object@x,y=object@y,r=object@r,default.unit="native",gp=gp)
	}

setMethod("PlotSetBoundaries","ThreeCircles",.plotboundaries.circles )
setMethod("PlotSetBoundaries","TwoCircles",.plotboundaries.circles )

setMethod("PlotUniverse","ThreeCircles",.plotuniverse.square)
setMethod("PlotUniverse","TwoCircles",.plotuniverse.square)

setMethod("IntersectionMidpoints","TwoCircles",function(object){
	VI <- data.frame(	Midpoint.x =object@centroids.x,Midpoint.y= object@centroids.y)
	VI$hjust <- "center"
	VI$vjust <- "center"
	VS <- IndicatorString(as(object,"Venn"))

	# choose the side with the smaller circle
	goleft <- object@r[1]<=object@r[2]
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
	VI[VS=="00",c("Midpoint.x","Midpoint.y")] <- c(dark.matter.x,ybottom)
	VI[VS=="00","hjust"] <- dark.matter.hjust
	VI[VS=="00","vjust"] <- "bottom"
	VI
})

setMethod("IntersectionMidpoints","ThreeCircles",function(object){
	V <- as(object,"Venn")
	VI <- data.frame(	Midpoint.x =object@centroids.x,Midpoint.y= object@centroids.y)
	VI$hjust <- "center"
	VI$vjust <- "center"

	goleft <- object@r[2]<=object@r[3]
	smidge.x <- 0.01 * diff(UniverseRange(object)$x)
	smidge.y <- 0.01 * diff(UniverseRange(object)$y)
	ytop <- max(UniverseRange(object)$y) - smidge.y

	if (goleft) {
		dark.matter.x <- min(UniverseRange(object)$x)+smidge.x
		dark.matter.hjust <- "left"
	} else {
		dark.matter.x <- max(UniverseRange(object)$x)-smidge.x
		dark.matter.hjust <- "right"
	}

	VS <- IndicatorString(V)
	
	ybottom <- min(VisibleRange(object)$y) - 0.05 * abs(diff(VisibleRange(object)$y))
	VI[VS=="000",c("Midpoint.x","Midpoint.y")] <- c(dark.matter.x,ytop)
	VI[VS=="000","hjust"] <- dark.matter.hjust
	VI[VS=="000","vjust"] <- "top"
	VI
})
setMethod("SetLabels","ThreeCircles",function(object){
	y=object@y;x=object@x
	r <- object@r
	VLabels <- data.frame(Label=rep("unset",3),x=NA,y=NA,hjust=I("center"),vjust=I("center"))
	yscale <- diff(VisibleRange(object)$y); smidge <- 0.01*yscale
	VLabels[1,2:3] <- c( x[1], y[1] +r[1]-smidge );  VLabels$vjust[1] <- "top"
	VLabels[2,2:3] <- c( x[2], y[2]+r[2]-smidge );VLabels$vjust[2] <- "top"
	VLabels[3,2:3] <- c( x[3], y[3]+r[3]-smidge );VLabels$vjust[3] <- "top"
	
	VLabels
})
setMethod("SetLabels","TwoCircles",function(object){
	y=object@y;x=object@x
	r <- object@r
	yscale <- diff(VisibleRange(object)$y); smidge <- 0.01*yscale
	VLabels <- data.frame(Label=rep("unset",2),x=NA,y=NA,hjust=I("center"),vjust=I("center"))
	VLabels[1,2:3] <- c( x[1], y[1]+r[1]-smidge );  VLabels$vjust[1] <- "top"
	VLabels[2,2:3] <- c( x[2], y[2]+r[2]-smidge);VLabels$vjust[2] <- "top"
	VLabels
})

if (FALSE) {
# now computed via drawnVenn
setMethod("Areas","TwoCircles",function(object){
	res <- object@areas
	names(res) <- IndicatorString(object)
	res
})
}

setMethod("Areas","ThreeCircles",function(object){
	res <- object@areas
	names(res) <- IndicatorString(object)
	res
})


.ComputeAreaCircles <- function(object,nintervals) {
	V <- as(object,"Venn")
	setcount <- NumberOfSets(V) 
	if (!setcount %in% c(2,3)) {
		stop("Wrong number of circles")
	}
	if(missing(nintervals)) {
		nintervals <- 1000
		nintervals <- 100
	}
	centres.x <- object@x
	centres.y <- object@y
	radius <- object@r
	xrange <- VisibleRange(object)$x
	yrange <- VisibleRange(object)$y

	agrid <- expand.grid(x=seq(xrange[1],xrange[2],length=nintervals),
					y=seq(yrange[1],yrange[2],length=nintervals))

	agrid$Count <- 0
	agrid$Var1 <- as.numeric((agrid$x-centres.x[1])^2+(agrid$y-centres.y[1])^2 < radius[1]^2)
	agrid$Var2 <- as.numeric((agrid$x-centres.x[2])^2+(agrid$y-centres.y[2])^2 < radius[2]^2)
	if (setcount==3) {
		agrid$Var3 <- as.numeric((agrid$x-centres.x[3])^2+(agrid$y-centres.y[3])^2 < radius[3]^2)
	}
	VI <- data.frame(Indicator(V),check.names=FALSE)
	VI$roworder <- 1:nrow(VI)
	varSetNames <- SetNames(V)
	colnames(agrid)[-c(1:3)] <- varSetNames
	V.count <- aggregate(agrid[,"Count",drop=FALSE],by=agrid[,varSetNames],FUN=length)
	VI <- merge(VI,V.count,all.x=TRUE)
	VI$Count[is.na(VI$Count)] <- 0
	V.count <- aggregate(agrid[,"x",drop=FALSE],by=agrid[,varSetNames],FUN=mean)
	VI <- merge(VI,V.count,all.x=TRUE)
	V.count <- aggregate(agrid[,"y",drop=FALSE],by=agrid[,varSetNames],FUN=mean)
	VI <- merge(VI,V.count,all.x=TRUE)
	VI$Area <- VI$Count/sum(VI$Count) * (xrange[2]-xrange[1]) * (yrange[2]-yrange[1])
	VI <- VI[order(VI$roworder),]
	VI <- VI[,c(varSetNames ,"Area","x","y")]
	colnames(VI) <- c(varSetNames ,"Area","Centroid.x","Centroid.y")
	# should be ok but check VI is still in the right order
		
	object@areas <- VI$Area
	object@centroids.x <- VI$Centroid.x
	object@centroids.y <- VI$Centroid.y
	object@nintervals <- nintervals
	object
}

setMethod("ComputeAreas","ThreeCircles",.ComputeAreaCircles)
setMethod("ComputeAreas","TwoCircles",.ComputeAreaCircles)

