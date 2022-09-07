#################################################################################
# Author: Guillaume Boglioni Beaulieu
# Description: Function to generate the data and 3D Figures of PIBD examples
# Last update: 06/09/2022
#################################################################################

# Libraries
library('rgl')
library('latticeExtra')
library('gridExtra')
library('grid')
library('ggplot2')

# Load external functions
source("3D_visualisation_functions.R")
source("generate_PI.R")

# Arguments of function 'p.i.figures'
#   seed:             Seed for simulations
#   type:             Identifier of the type of dependence
#   n:                Sample size   
#   my.col:           Colour scale of the points on the scatterplot
#   alpha:            Parameter used in certain examples 
#                     ('triangles', 'bernstein', 'bernstein4D' and 'copula')
#   pt.size:          Size of the points on the plots
#   plot.type:        Type visualisation: '2D.black','2D.colour', '2D.matrix','3D'
#   var.ind:          For plot.type = '2D, black', '2D.colour' or '3D.matrix', 
#                     order of appearance of the variables on the scatterplots
#   emp.colors:       TRUE to have empirical estimation of point concentration
#                     on 3D scatterplots
#   scale.col:        TRUE for a RELATIVE color scale 
#                     (i.e. full color palette of 'my.col' is used)
#   rad:              Radius of balls in empirical estimation of points concentration
#   grid.method:      TRUE if we want the 3D color coding using a fixed grid
#   g:                Integer for the size of the grid
#   adjust.pts.size:  TRUE to have the size of points proportional
#                     to the concentration around them
#   type.3d:          Type of points: 'p' is for points, 's' is for 3D spheres
#   legend:           For 2D coloured scatterplots, should a legend be displayed
#   margins:          If "TRUE", then simulated data will have marginals given by
#                     arguments mar1, mar2, mar3 (specifying quantile functions)

p.i.figures <- function(seed=1, type = 'bernstein', n = 1000, 
                        my.col = c("blue", "white", "red"), alpha = 1, 
                        pt.size = 5, plot.type = '2D', var.ind = c(1,2,3),
                        emp.colors = F, scale.col = F, rad=1/10, grid.method=F, 
                        g=8, adjust.pts.size=F, type.3d = 'p', legend=T, 
                        margins=F, mar1 = qnorm, mar2 = qexp, mar3 = qlnorm){
  
  # Set seed for random data generation
  set.seed(seed)
  
  # Generate data (according to 'type')
  data <- generate.PI(n, type, alpha)
  # At first, we set the "colors" of points to be a single color
  color <- my.col[1]
  
  # Color for the '4D example' (Example 3.7)
  if (type == 'bernstein4D'){
    color <- myColorRamp(my.col, data$U4)
  }
  
  # Update margins
  if(margins){
    data <- as.data.frame(
      cbind(X1=mar1(data[,1]), X2=mar2(data[,2]), X3=mar3(data[,3]))
    )
  }
  
  # Empirical estimation of 'h' for all points
  # ('h' calculated with function 'pts.concentration')
  if (emp.colors){
    # Option to have 3D plots where points vary in size (based on 'h' values)
    if(adjust.pts.size){
      # Values of the concentration index
      h.values <- pts.concentration(data, rad, scaled=scale.col, 
                                    grid.method=grid.method, g=g, return.h = T)
      # Values on the [0,1] scale (for color coding)
      f.values <- pts.concentration(data, rad, scaled=scale.col, 
                                    grid.method=grid.method, g=g)
      # For the grid method, the 'data' needs to be the data points of the grid
      if(grid.method){data <- as.data.frame(create.grid(g=g))}
      n.balls <- length(data[,1])
      data <- cbind(data, h.values, f.values)
      # Create "factors" for the size of the points 
      # (with max value the max of "h.values" in the sample)
      size <- as.numeric(cut(c(h.values,0, max(h.values)), 100))
      dataList <- split(data, size[-c(n.balls+1,n.balls+2)])
      # Find colors for the different 'size' factors
      # (we take the average colour of points having the same 'size')
      colors <- rep(NA, length(dataList))
      for(i in seq_along(dataList)){
        colors[i] <- myColorRamp(colors=my.col, mean(dataList[[i]]$f.values)) 
      }
      
    }else{
      color  <- myColorRamp(colors=my.col, 
                            values=pts.concentration(data, rad, scaled=scale.col, 
                                                     grid.method=grid.method, g=g)) 
    }  
  }
  
  # Plot 2D scatterplot
  if(plot.type=='2D.black'){
      par(mfrow=c(1,2))
      par(mar = c(5, 4.5, 1.5, 1.5))
      plot(data[,var.ind[1]], data[,var.ind[3]], cex=pt.size, pch = 20, 
           xlab = colnames(data)[var.ind[1]], 
           ylab =colnames(data)[var.ind[3]], cex.axis = 1.25, cex.lab = 1.75)
      plot(data[,var.ind[2]], data[,var.ind[3]], cex=pt.size, pch = 20, 
           xlab = colnames(data)[var.ind[2]], 
           ylab =colnames(data)[var.ind[3]], cex.axis = 1.25, cex.lab = 1.75)
      }
  if(plot.type=='2D.colour'){
      return(scatter2D.plus(data[,var.ind[1]], data[,var.ind[2]], 
                            data[,var.ind[3]], my.colors = my.col,
                            lab.size=1.75, pt.size=pt.size, 
                            X.lab = colnames(data)[var.ind[1]], 
                            Y.lab =colnames(data)[var.ind[2]]))
  }
  
  # Plot matrix of 2D scatterplot
  if(plot.type=='2D.matrix'){
    
    # Transform data to the ECDF
    ecdf.data <- apply(data, 2, rank)/n
    
    p12 <- scatter2D.plus(ecdf.data[,var.ind[1]], ecdf.data[,var.ind[2]], 
           ecdf.data[,var.ind[3]], pt.size = pt.size, my.color=my.col, 
           legend=legend, X.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[1]])*")"), 
           Y.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[2]])*")"))
    p13 <- scatter2D.plus(ecdf.data[,var.ind[1]], ecdf.data[,var.ind[3]], 
           ecdf.data[,var.ind[2]], pt.size = pt.size, my.color=my.col, 
           legend=legend, X.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[1]])*")"), 
           Y.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[3]])*")"))
    p23 <- scatter2D.plus(ecdf.data[,var.ind[2]], ecdf.data[,var.ind[3]], 
           ecdf.data[,var.ind[1]], pt.size = pt.size, my.color=my.col, 
           legend=legend, X.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[2]])*")"), 
           Y.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[3]])*")"))
    
    p11 <- ggplot(data, aes(x=data[,var.ind[1]])) + geom_histogram(aes(y=..density..), 
           colour="black", fill="lightblue", binwidth=2*IQR(data[,var.ind[1]])/n^(1/3))+ 
           labs(x = colnames(data)[var.ind[1]], y = "") + theme_classic()
    p22 <- ggplot(data, aes(x=data[,var.ind[2]])) + geom_histogram(aes(y=..density..), 
           colour="black", fill="lightblue", binwidth=2*IQR(data[,var.ind[2]])/n^(1/3))+ 
           labs(x = colnames(data)[var.ind[2]], y = "") + theme_classic()
    p33 <- ggplot(data, aes(x=data[,var.ind[3]])) + geom_histogram(aes(y=..density..), 
           colour="black", fill="lightblue", binwidth=2*IQR(data[,var.ind[3]])/n^(1/3))+ 
           labs(x = colnames(data)[var.ind[3]], y = "") + theme_classic()
    
    p21 <- scatter2D.plus(ecdf.data[,var.ind[1]], ecdf.data[,var.ind[2]], 
           ecdf.data[,var.ind[3]], pt.size = pt.size/2, 
           X.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[1]])*")"), 
           Y.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[2]])*")"), 
           my.colors = "black", legend = F)
    p31 <- scatter2D.plus(ecdf.data[,var.ind[1]], ecdf.data[,var.ind[3]], 
           ecdf.data[,var.ind[2]], pt.size = pt.size/2, 
           X.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[1]])*")"), 
           Y.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[3]])*")"), 
           my.colors = "black", legend = F)
    p32 <- scatter2D.plus(ecdf.data[,var.ind[2]], ecdf.data[,var.ind[3]], 
           ecdf.data[,var.ind[1]], pt.size = pt.size/2, 
           X.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[2]])*")"), 
           Y.lab = bquote(hat(F)*"("*.(colnames(data)[var.ind[3]])*")"), 
           my.colors = "black", legend = F)
    
    print(grid.arrange(p11, p12,p13, p21, p22, p23, p31, p32, p33, ncol=3, nrow = 3))
    
  }
  
  # Plot 3D scatterplot
  if(plot.type=='3D'){
    
    # Option where we alter the size of each point
    if(adjust.pts.size){
      # Create plot (points are added in a second step)
      with(data, plot3d(NULL, NULL, NULL, 
                        colnames(data)[1], ylab = colnames(data)[2], 
                        zlab = colnames(data)[3], col='white', size=0))
      par3d(windowRect = 50 + c(0, 0, 500, 500)) # adjust size of window
      
      # Use separate calls of points3d() to plot points of each size
      for(i in seq_along(dataList)){
        if(type.3d == 'p'){with(dataList[[i]], points3d(U1, U2, U3,  col=colors[i], 
                           size=pt.size*n*mean(dataList[[i]]$h.values)))}
        if(type.3d == 's'){with(dataList[[i]], spheres3d(U1, U2, U3, col=colors[i], 
                           radius=pt.size*n*mean(dataList[[i]]$h.values)/500))}
      }
      
    } else{
    
      if (grid.method){
        data <- create.grid(g=g)
      }
      open3d()
      plot3d(data[,1], data[,2], data[,3], col=color, size=pt.size, type=type.3d,
             xlab=colnames(data)[1], ylab=colnames(data)[2], zlab=colnames(data)[3])
      par3d(windowRect = 50 + c(0, 0, 500, 500)) # adjust size of window
    } # End 'else' for adjust.pts.size 
  } # End 'else' for 3D scatterplot
} # End function


