# Load needed package
library('pdist')

# Functions to find fraction of volume of a ball sitting OUTSIDE the [0,1]^3 "cube"
# (when the ball extends outside the cube in either 1, 2 or 3 directions)
# The arguments d1,d2,d3 are as defined in the article by Freireich et al. (2010)

f.face <- function(d,r){
  dr <- d/r
  return((3*dr^2 - dr^3)/4)
}

f.edge <- function(d1,d2,r){
  a  <- 1-d1/r
  b  <- 1-d2/r
  x2 <- 1 - a^2 - b^2
  if(x2 <= 0){
    return(0)}
  else{
    x <- sqrt(x2)
    return((2*a*b*x-(3*a-a^3)*atan(x/b)-(3*b-b^3)*atan(x/a)+2*atan(x*a/b)
            +2*atan(x*b/a))/(4*pi))
  }
}

f.corner <- function(d1,d2,d3,r){
  a  <- 1-d1/r
  b  <- 1-d2/r
  c  <- 1-d3/r

  if((a^2+b^2+c^2) >= 1){
    return(0)}
  else{
    A  <- sqrt(1-a^2-c^2)
    B  <- sqrt(1-b^2-c^2)
    return(f.edge(d1,d2,r)/2 - 
          (6*a*b*c-2*a*A*c-2*b*B*c-(3*b-b^3)*atan(c/B)-(3*a-a^3)*atan(c/A)+
          (3*c-c^3)*(atan(A/a)-atan(b/B))+2*(atan(c*a/A)+atan(c*b/B)))/(8*pi))
  }
}

# Functions for the fraction of volume of the ball comprised inside the "cube"
# (different functions apply, depending on in how many directions (1,2 or 3) 
# the ball extends outside the cube)
v.1inter <- function(d,r){
  1-(f.face(d=d,r=r))
}

v.2inter <- function(d1,d2,r){
  1-(f.face(d=d1,r=r) + f.face(d=d2,r=r) - f.edge(d1=d1,d2=d2,r=r))
}

v.3inter <- function(d1,d2,d3,r){
  1-(f.face(d1,r)+f.face(d2,r)+f.face(d3,r)
     -f.edge(d1,d2,r)-f.edge(d1,d3,r)-f.edge(d2,d3,r)
     +f.corner(d1,d2,d3,r))
}

# Function 'v.star' returns the proportion (volume wise) of a Ball[p,r] 
# which is inside the unit cube [0,1]^3
# Arguments of function 'v.star':
#   p:  Vector of coordinates (x,y,z) of the point at the center of the ball
#   r:  Radius of the Ball

v.star <- function(p,r){
  p1 <- p[1]
  p2 <- p[2]
  p3 <- p[3]
  p_pm_r <- c(p1-r,p1+r,p2-r,p2+r,p3-r,p3+r)
  
  # Count in how many directions does the Ball[p,r] falls outside [0,1]^3
  count  <- sum(p_pm_r < 0) + sum(p_pm_r > 1)
  # Distances from p to edges of [0,1]^3 
  # (for the directions for which the Ball[p,r] extends outside [0,1]^3)
  bounds <- c(p[p < r], 1 - p[p + r > 1])
  # (d is (d1, d2, d3) with d1, d2, d3 corresponding to Freireich et al. (2010))
  d <- r - bounds
  
  if (count == 0){
    delta <- 1
  }
  # Case the Ball[p,r] extends outside [0,1]^3 in one direction
  if (count == 1){
    delta <- v.1inter(d=d[1], r)
  }
  # Case the Ball[p,r] extends outside [0,1]^3 in two directions
  if (count == 2){
    delta <- v.2inter(d1=d[1], d2=d[2], r)
  }
  # Case the Ball[p,r] extends outside [0,1]^3 in three directions
  if (count == 3){
    delta <- v.3inter(d1=d[1], d2=d[2], d3=d[3], r)
  }
  if (count > 3){
    print("Error: 'count' cannot be >3")
  }
  
  #Return result
  return(delta)
}

# Function to create a 3D grid in [0,1]^3
create.grid <- function(g=10){
U1 <- NULL
U2 <- NULL
for (i in 1:g){
  U1 <- c(U1, rep(i,g^2))
  U2 <- c(U2, rep(i,g))
}
U3 <- rep(1:g,g^2)
return(cbind(U1,U2,U3)/(g+1))
}

# Function that finds the 'concentration index' h() for all points of a dataset 
# (data assumed in [0,1]^3)
# Arguments:
#     data:         (nx3) matrix containing the data
#     rad:          radius of the small balls around each point
#     scaled:       TRUE to have 0-to-1 scaling on a RELATIVE basis 
#                   (i.e. full colour palette is used)
#     return.h:     TRUE to return the value of 'h', otherwise returns 'f(h)'
#     grid.method:  TRUE to use "fixed grid" method
#     g:            Size of grid, if "fixed grid" method is used

pts.concentration <- function(data, rad, scaled=F, return.h = F, 
                              grid.method = F, g=10){
  # Sample size
  n         <- length(data[,1])
  # Initial (not smoothed) concentration
  distances <- as.matrix(dist(data, upper = TRUE))
  # Count number of points within radius "rad" of any point
  hp      <- (rowSums(distances <= rad)-1)/(n-1) #"-1" excludes the point itself
  hp      <- hp / apply(data, 1, v.star, r = rad)
  
  # Grid (g*g*g)
  grid <- create.grid(g=g)
  
  # Values that depend on whether we use the 'fixed grid' method or not
  if (grid.method){
    # Number of balls is the number of points of the grid, g^3
    n.balls <- g^3
    # Distances between each point of the grid and each point of the sample
    grid.dist <- as.matrix(pdist(X = grid, Y = data))
  } else{
    n.balls <- n
    grid.dist <- distances
  }
  
  # Find the final estimates of "h"
  # (each point's value is a weighted average of its nearest neighbours)
  final.hp <- rep(NA, n.balls)
  for(i in 1:n.balls){
    closest.pt  <- grid.dist[i,] <= rad
    dist        <- grid.dist[i,][closest.pt]
    # If there are ZERO points in the ball, hp should just be 0
    if (length(dist) == 0){final.hp[i] <- 0}
    # If there is only ONE point in the ball, hp should be that of this point 
    if (length(dist) == 1){final.hp[i] <- hp[closest.pt]}
    # If there are at least TWO points
    if (length(dist)  > 1){
      # If one point is EXACTLY the center of the ball, 
      # we give it a "distance" equal to the distance with THE closest neighbour
      if(min(dist)==0){dist[which.min(dist)] <- sort(dist)[-1][1]}
      # Weights are simply the reverse of the distances
      weights     <- 1/dist
      # Compute final value: a weighted average of all points inside the ball
      final.hp[i] <- sum(hp[closest.pt] * weights / sum(weights))
    }
  }  
  
  # Parameters for scaling to a color scale (0-to-1)
  v = 4*pi*rad^3/3
  # Case where we want the entire color palette to be used (i.e. 'relative' scale)
  if (scaled){
    a <- min(final.hp)
    c <- max(final.hp)
    if (a > v){
      print("Warning: 'a' is bigger than 'v'. Replacing 'v' by (a+c)/2")
      v <- (a+c)/2
    }
    if (c < v){
      print("Error: 'c' must be larger than 'v'.")
    }
    # Case where the color palette is absolute (not relative to the given sample)  
  } else{
    a <- 0
    c <- 2 * rad/sqrt(3)
  }
  beta      <- -log(2)/log((v-a)/(c-a))
  
  # Put the points on a 0-to-1 scale
  pt.scaled <- ((final.hp-a)/(c-a))^beta
  # Cap values at '1'
  ind <- pt.scaled > 1
  pt.scaled[ind] <- 1
  # Return results
  ifelse(return.h, return(final.hp), return(pt.scaled))
}

# Function that assigns a range of colors to a range of values
# Values must be between 0 and 1
myColorRamp <- function(colors, values) {
  x <- colorRamp(colors)(values)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}

# Function 'scatter2D.plus' creates scatterplots of two variables,
# where a third variable is represented by a colour
scatter2D.plus <- function(x,y,z, legend = TRUE, 
                           my.colors = c("magenta", "white", "cornflowerblue"), 
                           X.lab = "X", Y.lab = "Y", nb.levels = 50, 
                           pt.size=3, lab.size=1){
  Rx <- diff(range(x))
  Ry <- diff(range(y))
  levelplot(z ~ x + y, colorkey = legend, panel = panel.levelplot.points,
            col.regions = colorRampPalette(my.colors)(nb.levels),  
            xlab=list(label = X.lab, cex=lab.size), 
            ylab=list(label = Y.lab, cex=lab.size), cex = pt.size, 
            col =  "black", pch = 21, xlim = c(min(x)-0.05*Rx, max(x)+0.05*Rx), 
            ylim = c(min(y)-0.05*Ry, max(y)+0.05*Ry))  
}

