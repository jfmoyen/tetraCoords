
#### Tetrahedral coordinates
tetraCoords <- function(A,B,C,D){
  #' Calculate the (x,y,z) plotting coordinates from the values of 4 components
  #' @param{A,B,C,D} Coordinates of a point in (a,b,c,d) form
  #' @import tibble

  #### Definition of tetrahedron apices
  a_A <- c(0,0,0)
  a_B <- c(1/2,sqrt(3)/2,0)
  a_C <- c(1,0,0)
  a_D <- c(1/2,sqrt(3)/6,sqrt(6)/3)

  #### Sum of coordinates
  SS <- A+B+C+D

  #### Normalize
  A <- A/SS
  B <- B/SS
  C <- C/SS
  D <- D/SS

  #### New coords

  x <- B/2 + C + D/2
  y <- sqrt(3)/2*B + sqrt(3)/6*D
  z <- sqrt(6)/3 * D

  return(tibble(x,y,z))

}


#### Assemble a matrix
prepareData <- function(data,coord_names,object_names,object_type="Nature"){
  #' Convenience function. Assemble and converts a matrix from a vector
  #' @param data a vector that holds the data
  #' @param coord_names Names of the coordinates (will be colnames)
  #' @param object_names Names of the objects (will be rownames)
  #' @param object_type String describing what the "objects" are (minerals, end-members...)
  #' Will be used in the tibble to identify them
  #' @import tibble

  ncols <- length(coord_names)

  ee <- matrix(data,byrow=T,ncol=ncols)
  colnames(ee) <- coord_names

  out <- tibble({{object_type}} := object_names) %>%
    bind_cols(ee)

  return(out)
}


#### Empty tetra canvas for RGL
drawggTetraPlot <-function(apicesNames  = c("A","B","C","D"), axes = F, zscale = 3.5, show = T){
  #' Draw an empty canvas with the edges of tetraedron
  #' @param apicesNames Names of the apices to be drawn on plot
  #' @param axes (F) Boolean, plot the cartesian axes?
  #' @param zscale A magical number to be passed to rgldev, adjusting the vertical scale.
  #' Apparently 3.5 looks nice and gives apparent square aspect ratio.
  #' @param show (T) Boolean. If true, immediately show the plot otherwise just return it.
  #' @import tibble dplyr rgl devout devoutrgl ggrgl ggplot2

  ## Define apices cordinates
  a_A <- c(0,0,0)
  a_B <- c(1/2,sqrt(3)/2,0)
  a_C <- c(1,0,0)
  a_D <- c(1/2,sqrt(3)/6,sqrt(6)/3)

  tetra <- rbind(a_A,a_B,a_C,a_D)
  colnames(tetra)<-c("x","y","z")

  # To tibble
  tetraVerts <- as_tibble(tetra)%>%
    mutate(idx = seq(n()))

  # Build the tetrahedron "mesh" by defining edges with the right connectivity
  tetraEdges <- tibble(start=c(1,1,1,2,2,3),end=c(2,3,4,3,4,4)) %>%
    left_join(tetraVerts,by=c(start="idx")) %>%
    left_join(tetraVerts,by=c(end="idx")) %>%
    rename(xstart=x.x,ystart=y.x,zstart=z.x,xend=x.y,yend=y.y,zend=z.y)

  # Plot an empty plot with just the tetra
  p <- ggplot(tetraEdges) +
    geom_segment_3d(aes(x=xstart, y=ystart, z=zstart, xend=xend, yend=yend, zend=zend),
                    size = 0.2) +
    geom_text_z(data=tetraVerts,aes(x=x,y=y,z=z,label=apicesNames ),size=18)+
    theme(legend.position = 'none') +
    coord_equal()

  # Adjustment
  if(axes){
    p <- p + theme_ggrgl()
  }else{
    p <- p +theme_void()
  }

  # Show the plat if required
  if(show){
    devoutrgl::rgldev(zscale=zscale)
    p
  }else{
    p
  }

}


##### Tibble to matrix #####
dataTib2Mat <- function(data,coord_names,object_type){
  #' Convert a tibble-formatted data into a matrix fit for coordinate mapping
  #' @param data a tibble that holds the data
  #' @param coord_names Names of the coordinates (will be colnames)
  #' @param object_type Name of the column from which rownames will be taken

  rn <- data[[object_type]]
  if(  any(duplicated(rn)) ){stop("Objects must be unique")}

  ee <- select(data,all_of(coord_names)) %>% as.matrix
  rownames(ee) <- rn

  return(ee)
}


##### Coordinate mapping #####
coordMap <- function(data,transMat){
  #' Perform coordinate mapping on a data matrix
  #' @param data the tibble containing data
  #' @param transMat Transformation matrix (a matrix named in rows and cols)


  # Transformation matrix must be square to be invertible
  # (well, amongst other things...)
  if(ncol(transMat) != nrow(transMat) ){stop("Transformation matrix must be square")}
  imat <- solve(transMat)

  # Tibbles don't particularly like matrix mulmtiplication, so...
  qq<-t(apply(
    select(data,colnames(transMat)),
    1,
    function(z){
      ee <- z%*%imat
      ee
    }
  ))
  colnames(qq)<-rownames(transMat)

  return(data %>% bind_cols(qq))
  }
