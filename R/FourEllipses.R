warning("Entering FourEllipses")

#######################
# Four squares
setClass("FourEllipses",
	 representation(
	"DrawnVenn","Venn"))



compute.E4 <- function(V,doWeights=FALSE,s=.25) {
	if (doWeights) { warning("Cant do a weighted E4") }
	if (NumberOfSets(V) != 4) { stop("fournotfour")}
diagram <- new("graphNEL")
regions <- new("graphNEL")
VS <- IndicatorString(Venn(NumberOfSets=4))
regions <- addNode(VS,regions)
nodeDataDefaults(regions,"xy") <- matrix(as.numeric(NA),nrow=0,ncol=2)
nodeData(regions,"0000","xy") <- list(c(0,9))
nodeData(regions,"1000","xy") <- list(c(-3.5,6))
nodeData(regions,"0100","xy") <- list(c(3.5,6))
nodeData(regions,"0010","xy") <- list(c(-6,0))
nodeData(regions,"0001","xy") <- list(c(6,0))
nodeData(regions,"1100","xy") <- list(c(0,3))
nodeData(regions,"1010","xy") <- list(c(-4.5,3))
nodeData(regions,"1001","xy") <- list(c(3,-2))
nodeData(regions,"0110","xy") <- list(c(-3,-2))
nodeData(regions,"0101","xy") <- list(c(4,3))
nodeData(regions,"0011","xy") <- list(c(0,-4.5))
nodeData(regions,"1110","xy") <- list(c(-2,0))
nodeData(regions,"1101","xy") <- list(c(2,0))
nodeData(regions,"1011","xy") <- list(c(1,-2.75))
nodeData(regions,"0111","xy") <- list(c(-1,-2.75))
nodeData(regions,"1111","xy") <- list(c(0,-1.5))

	drawnVenn <- new("DrawnVenn",regions=regions,diagram=diagram)
	new("FourEllipses",drawnVenn,V)
}	

setMethod("SetLabels","FourEllipses",function(object){
	data.frame(Label=as.character(1:4),
		x=c(-5,5,-7.5,7.5),
		y=c(9,9,-3,-3),
		hjust=c("center"),vjust=c("center"),
		stringsAsFactors=FALSE)

})

setMethod("PlotNodes","FourEllipses",function(object,gp){
	warning("PlotNodes not implemented for FourEllipses")
})

setMethod("PlotUniverse","FourEllipses",function(object,gp){
	ur <- UniverseRange(object)
	grid.rect(x=mean(ur$x),y=mean(ur$y),width=diff(ur$x),height=diff(ur$y),default.units="native")
})
# should inherit from the method defined for "DrawnVenn"
#setMethod("IntersectionMidpoints","FourEllipses",function(object){
#	dv <- as(object,"DrawnVenn")
#	regions <- dv@regions
#	cxy <- nodeData(regions,attr="xy")
#	df <- data.frame(do.call(rbind,cxy))
#	colnames(df) <- c("Midpoint.x","Midpoint.y")
#	df$VS <- rownames(df)
#	df <- df[match(df$VS,IndicatorString(as(object,"Venn"))),]
#	df
#})
setMethod("UniverseRange","FourEllipses",function(object) {
	vr <- VisibleRange(object)
	vr
})
setMethod("VisibleRange","FourEllipses",function(object) {
		list(x=c(-10,10),y=c(-7,15))
})
setMethod("PlotSetBoundaries","FourEllipses",function(object,gp) {
	phi <- 0.8; dex <- 1.7;dey <- 2.5; a<- 7.6; e<- 0.9
	x0 <- c( -0.9, -5.0)
	
	nintervals <- 200
	ellipse <- function(f1,phi,e,a,nintervals=200,gp=gpar()) {
		twoc <- a* e* 2
		f2 <- f1+  twoc*c(-cos(phi),sin(phi))
		theta <- seq(0,2*pi,length=nintervals)
		r <- (a * (1-e^2))/(1+e*cos(theta+phi))
		x <- f1[1]+r*cos(theta)
		y <- f1[2]+r*sin(theta)
		grid.lines(x,y,default.units="native",gp=gp)
	}
	ellipse (x0+c(0,0),-phi ,e,-a,gp=gp[4])
	ellipse (x0+c(dex,0),phi ,e,a,gp=gp[3])
	ellipse (x0+c(-dey,dey),-phi ,e,-a,gp=gp[2])
	ellipse (x0+c(dex+dey,dey),phi ,e,a,gp=gp[1])

})