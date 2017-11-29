#makes a polygonal gate on a ggplot
#note that this needs some calibration that I haven't quite figured out, so the points you click on are not quite the correct locations on the graph
fcsGate = function(p,n=4, draw=T){
  require(grid)
  pData = ggplot_build(p)
  #xrange=pData$layout$panel_ranges[[1]]$x.range
  #yrange=pData$layout$panel_ranges[[1]]$y.range
  xrange = pData$plot$coordinates$limits$x;
  yrange = pData$plot$coordinates$limits$y;
  xName=as.character(pData$plot$mapping$x)
  yName=as.character(pData$plot$mapping$y)
  points = data.frame(n=1:n)
  points[[xName]]=NA;
  points[[yName]]=NA;
  seekViewport(grid.ls()[[1]])
  for(i in 1:n){
    y <-  grid.locator("npc")
    y <- as.numeric(substring(y, 1, nchar(y)-3))
    message(sprintf("x=%g, y=%g",y[1],y[2]))
    #points[[xName]][i] <- min(xrange) + y[1]*diff(xrange)
    #points[[yName]][i] <- min(yrange) + y[2]*diff(yrange)
    points[[xName]][i] <- min(xrange) + y[1]/0.75*diff(xrange)
    points[[yName]][i] <- min(yrange) + y[2]/0.75*diff(yrange)
    if(i>1 && draw){
      p = p+geom_segment(x=points[[xName]][i-1], y=points[[yName]][i-1], xend=points[[xName]][i], yend=points[[yName]][i], colour="red"); print(p)
    }
  }
  return(points);
}

#returns a logical vector containing whether each point in the dataframe is within the given polygon (made with fcsGate)
isWithin = function(data, polygon){
  require(sp)
  point.in.polygon(data[[names(polygon)[2]]], data[[names(polygon)[3]]], polygon[[2]], polygon[[3]])
}
