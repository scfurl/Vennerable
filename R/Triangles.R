warning("Entering Triangles.R")
###########################
setClass("ThreeTriangles",
	representation(SetShapes="list",
			IntersectionShapes="list",
			"DrawnVenn",
			"Venn",areas="data.frame"))

.inscribetriangle.feasible <- function(wghts) {
	w0 <- 1- sum(wghts)
	stopifnot(all(wghts <= 1) & all(wghts>=0) & w0>=0)
	wa <- wghts[1];wb <- wghts[2]; wc <- wghts[3]
	Delta <- w0^2 - 4 * wa * wb * wc
	return (Delta>=0)
}
.inscribetriangle.compute <- function (wghts) {
	wa <- wghts[1];wb <- wghts[2]; wc <- wghts[3]
	stopifnot(.inscribetriangle.feasible(wghts)) 
	pa <- (1-wc)
	pb <- (wb+wc-wa-1)
	pc <- wa * (1-wb)
	sc <- if (wa>0) { 
		(-pb-sqrt( pb^2 - 4 * pa * pc))/(2*pa) 
		} else if (wb+wc<1) { 
		(1-wb-wc)/(1-wc) 
		} else { 
		0 
		} 
	sb <- if (sc>0 ) { 1 - wa/sc } else { wc/(1-wb) }
	sa <- wb/(1-sc)
	c(sc,sa,sb) # nb order around triangle
}
.inscribetriangle.inscribe <- function(xy,wghts) {
	scalef <- NA
	
	isfeasible <- .inscribetriangle.feasible(wghts)
	if (!isfeasible) {
		scalef <- 4 * wghts[1]*wghts[2]*wghts[3]/(1-sum(wghts))^2
		scalef <- scalef^(1/3)
		wghts <- wghts / (scalef*1.001)
		isfeasible <- .inscribetriangle.feasible(wghts)
		stopifnot(!isfeasible)
	}
	if (!isfeasible) return(list(feasible=FALSE))
	scab <- .inscribetriangle.compute (wghts)
	inner.xy <- (1-scab)*xy + scab * (xy[c(2,3,1),])
	return(list(feasible=TRUE,inner.xy=inner.xy,scalef=scalef))
}

compute.T3 <- function(V,doWeights=TRUE) {
	stopifnot(NumberOfSets(V)==3)
	wght <- Weights(V)
	VS <- IndicatorString(V)
	Vorig <- V
	if (!doWeights) { 
		Wsum <- apply(Indicator(V),1,sum)
		wght[Wsum==0] <- 2
		wght[Wsum==1] <- 4
		wght[Wsum==2] <- 1
		wght[Wsum==3] <- 1
		V@IndicatorWeight[,".Weight"] <- wght
	}


	# with wa, wb, wc on the outside and w0 the universe.

	WeightUniverse <- .WeightUniverse(V)
	WeightVisible <- .WeightVisible(V)
	WeightInvisible <- WeightUniverse-WeightVisible 


	w0ratio <- WeightInvisible/WeightVisible 
	# we ignore w0 from now on and all the other weights sum to one
	wght <- wght/WeightVisible 
	
	# the inner triangle contains wab,wbc,wca and 
	# the innest contains wabc


	# the outer triangle
	wa <- wght[VS=="100"]
	wb <- wght[VS=="010"]
	wc <- wght[VS=="001"]
	outer.weights <- c(wa,wb,wc)
	outer.innerw <- 1 - sum(outer.weights)
	outer.inner.ratios <- outer.weights/outer.innerw #  ratio of each wa, wb,wc to pooled inner weights

	outer.feasible <- .inscribetriangle.feasible(outer.weights)

	# the inner triangle 
	wab <- wght[VS=="110"]
	wbc <- wght[VS=="011"]
	wca <- wght[VS=="101"]
	wabc <-  wght[VS=="111"]

	inner.weights <- c(wab,wbc,wca)
	inner.innerw <- wabc
	# we resclae the inner weights...
	sf <- (sum(inner.weights)+inner.innerw)
	Weight.Inner <- sf * WeightVisible
	if (sf>0) {	
		inner.weights <- inner.weights/sf 
		inner.feasible <- .inscribetriangle.feasible(inner.weights)
	} else {
		inner.feasible <- FALSE
	}
	
	if (inner.feasible & outer.feasible) {
		# whole triangle should have area in Weights
		side <- sqrt(4 * WeightVisible /(3*sqrt(3)))
		angles <- pi/2-c(0,2*pi/3,4*pi/3)
		outer.xy <- t(sapply(angles,function(a)c(x=side * cos(a),y= side * sin(a))))
		
		inner <- .inscribetriangle.inscribe(outer.xy,wghts=outer.weights)
		inner.xy <- inner$inner.xy
		innest <- .inscribetriangle.inscribe(inner.xy,wghts=inner.weights)
		innest.xy=innest$inner.xy
		# finally we construct the outside triangle
		# outer.xy is equilateral with centre at zero, so just scale 
		# if inner triangle has area A and rim has area A' then scaling 
		# is (A'+A)=s^2 A so s^2=1+A'/A. A'/A is the w0ratio calculated above
		outest.xy <- outer.xy * sqrt( 1+ w0ratio)
	} else { 
		if (inner.feasible) { # but not outer
			# so we make the inner equilateral of area Weight.Inner
			side <- sqrt(4 * Weight.Inner /(3*sqrt(3)))
			angles <- pi/6-c(0,2*pi/3,4*pi/3)
			inner.xy <- t(sapply(angles,function(a)c(x=side * cos(a),y= side * sin(a))))
			innest <- .inscribetriangle.inscribe(inner.xy,wghts=inner.weights)
			innest.xy=innest$inner.xy
			# 
			outer.heights <- 2 * outer.inner.ratios * Weight.Inner /side # TODO wrong?
			outer.distance <- side/(2*sqrt(3))+ outer.heights
			outer.angles <- pi/2-c(0,2*pi/3,4*pi/3)
			outer.x <- outer.distance * cos(outer.angles)
			outer.y <- outer.distance *sin(outer.angles)
			outer.xy <- matrix(c(outer.x,outer.y),ncol=2,byrow=FALSE)
			# TODO THIS IS QUITE WRONG
			outest.xy <- outer.xy * sqrt( 1+ w0ratio)
		} else { # inner and out infeasible
		stop("Not implemented yet")
		}
	}
	

	outline.a.xy  <- matrix(c(outer.xy[1,],inner.xy[1,],innest.xy[1,],innest.xy[2,],inner.xy[3,]),ncol=2,byrow=TRUE)
	outline.b.xy  <- matrix(c(outer.xy[2,],inner.xy[2,],innest.xy[2,],innest.xy[3,],inner.xy[1,]),ncol=2,byrow=TRUE)
	outline.c.xy  <- matrix(c(outer.xy[3,],inner.xy[3,],innest.xy[3,],innest.xy[1,],inner.xy[2,]),ncol=2,byrow=TRUE)

	df <- Indicator(V)
	dfIx <- apply(df,1,paste,collapse="")
	IntersectionShapes <- vector("list",length(dfIx))
	names(IntersectionShapes) <- dfIx
	IntersectionShapes[["000"]]<-  outest.xy
	IntersectionShapes[["100"]] <-  matrix(c(outer.xy[1,],inner.xy[1,],inner.xy[3,]),ncol=2,byrow=TRUE)
	IntersectionShapes[["010"]] <-  matrix(c(outer.xy[2,],inner.xy[2,],inner.xy[1,]),ncol=2,byrow=TRUE)
	IntersectionShapes[["001"]] <-  matrix(c(outer.xy[3,],inner.xy[3,],inner.xy[2,]),ncol=2,byrow=TRUE)
	IntersectionShapes[["110"]] <-  matrix(c(inner.xy[1,],innest.xy[1,],innest.xy[3,]),ncol=2,byrow=TRUE)
	IntersectionShapes[["011"]] <-  matrix(c(inner.xy[2,],innest.xy[2,],innest.xy[1,]),ncol=2,byrow=TRUE)
	IntersectionShapes[["101"]] <-  matrix(c(inner.xy[3,],innest.xy[3,],innest.xy[2,]),ncol=2,byrow=TRUE)
	IntersectionShapes[["111"]] <-  matrix(c(innest.xy[1,],innest.xy[2,],innest.xy[3,]),ncol=2,byrow=TRUE)

	SetShapes=list(outline.a.xy,outline.b.xy,outline.c.xy)
############# beginning of replacement by drawnvenn object
	diagram <- new("graphNEL") # this remains empty for now
	regions <- new("graphNEL")
	nodeDataDefaults(regions,"xy") <- matrix(as.numeric(NA),nrow=0,ncol=2)
	VS <- IndicatorString(Venn(NumberOfSets=3))
	regions <- addNode(VS,regions)
	mids <- data.frame(t(sapply(IntersectionShapes,function(x)apply(x,2,mean))))
	outerty <- min(sapply(SetShapes,function(x)min(x[,2])))
	outesty <- min(IntersectionShapes[["000"]][,2])
	mids[VS=="000",]  <- c(0,mean(outerty,outesty))

	for (ix in 1:nrow(mids)) {
		vs <- rownames(mids)[ix]
		nodeData(regions,vs,"xy") <- list(c(mids[ix,1],mids[ix,2]))
		if (wght[VS==vs]==0) {
			nodeData(regions,vs,"xy") <- list(c(NA,NA))
		}	
	}
	drawnVenn <- new("DrawnVenn",regions=regions,diagram=diagram)
	res <- 	new("ThreeTriangles",SetShapes=SetShapes,
		IntersectionShapes=IntersectionShapes,drawnVenn,
		Vorig)
	res
}

setMethod("VisibleRange","ThreeTriangles",function(object) {
	object@IntersectionShapes[["000"]] <- NULL
	xrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,1])}))	
	yrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,2])}))
	list(x=xrange,y=yrange)
})

setMethod("UniverseRange","ThreeTriangles",function(object) {
	xrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,1])}))	
	yrange <- range(sapply(object@IntersectionShapes,function(x){range(x[,2])}))
	list(x=xrange,y=yrange)
})

# should inherit from the method defined for "DrawnVenn"
if(FALSE) {
setMethod("IntersectionMidpoints","ThreeTriangles",function(object){
	V <- as(object,"Venn")
	VS <- IndicatorString(V)
	mids <- 	data.frame(t(sapply(object@IntersectionShapes,function(x)apply(x,2,mean))))
	colnames(mids) <- c("Midpoint.x","Midpoint.y")
	# only the 000 is not a whole triangle
	# find the base of the others 
	outerty <- min(sapply(object@SetShapes,function(x)min(x[,2])))
	outesty <- min(object@IntersectionShapes[["000"]][,2])
	mids[VS=="000",c("Midpoint.x","Midpoint.y")]  <- c(0,mean(outerty,outesty))
	mids
})
}

setMethod("PlotSetBoundaries","ThreeTriangles",function(object,gp){
	grid.polygon(x=object@SetShapes[[1]][,1],y=object@SetShapes[[1]][,2],default.units="native",gp=gp[1])
	grid.polygon(x=object@SetShapes[[2]][,1],y=object@SetShapes[[2]][,2],default.units="native",gp=gp[2])
	grid.polygon(x=object@SetShapes[[3]][,1],y=object@SetShapes[[3]][,2],default.units="native",gp=gp[3])
})

setMethod("PlotUniverse","ThreeTriangles",function(object,gp){
	if(missing(gp)) { gp <- NULL }
	outest <- object@IntersectionShapes[["000"]]
	grid.polygon(x=outest[,1],y=outest[,2],default.units="native",gp=gp)
})


setMethod("SetLabels","ThreeTriangles",function(object){
	# pull out the outer triangle vertices
	outer.xy <- do.call(rbind,list(object@IntersectionShapes[["100"]][1,],
			object@IntersectionShapes[["010"]][1,],
			object@IntersectionShapes[["001"]][1,]))
	# scale in our out
	sf <- 0.9
	outer.xy <- outer.xy * sf
	VLabels <- data.frame(Label=rep("unset",3),x=NA,y=NA,hjust=I("center"),vjust=I("center"))
	VLabels[,2:3] <-outer.xy; 
	VLabels$vjust <- c("top","bottom","bottom")
	VLabels$hjust <- c("center","right","left")

	VLabels
})



.triangle.area <- function(xy) {
	if(nrow(xy)!=3) stop("Not a triangle")
	if(ncol(xy)!=2) stop("Not xy coords")
	edgexy <- matrix(c(xy[1,]-xy[2,],xy[2,]-xy[3,],xy[3,]-xy[1,]),ncol=2,byrow=TRUE)
	edges <- sqrt(apply(edgexy^2,1,sum))
	s <- sum(edges)/2
	heron <- sqrt(s*(s-edges[1])*(s-edges[2])*(s-edges[3]))
	heron
}

# should be replaced by DrawnVenn method 
setMethod("Areas","ThreeTriangles",function(object){
	V <- as(object,"Venn")
	VS <- IndicatorString(V)
	mids <- 	data.frame(Area=sapply(object@IntersectionShapes,.triangle.area))
	# only the 000 is not a whole triangle
	# find the diff from the others 
	mids$Area[VS=="000"]  <- mids$Area[VS=="000"] - sum(mids$Area[VS!="000"])
	res <- as.numeric(mids$Area)
	names(res) <- rownames(mids)
	res
})
setMethod("PlotNodes","ThreeTriangles",function(object,gp){
	warning("PlotNodes not implemented for ThreeTriangles")
})


	

