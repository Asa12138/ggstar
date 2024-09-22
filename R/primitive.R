#' @importFrom grid unit convertX convertY is.unit unit.c gpar
starGrob <- function(x=0.5, y=0.5,
                     starshape=1, 
                     angle=0, 
                     phase = 0,
                     gp = gpar(fill="black",
                               fontsize=2,
                               alpha=1,
                               col=NA,
                               lwd=0.5),
                     position.units = "npc", 
                     size.units="mm", ...){
    if (! all(starshape %in% seq_len(length(starshape_ntab)))){
        stop("the starshape should be one of 1 to 33 !")
    }
    N <- length(x)
    stopifnot(length(y)==N)
    angle <- deg2rad(x=angle)
    if (!is.unit(x)){x <- unit(x, position.units)}
    if (!is.unit(y)){y <- unit(y, position.units)}
    xv <- convertX(x, position.units, TRUE)
    yv <- convertY(y, position.units, TRUE)
    n <- match_n(starshape)
    if (is.null(gp)){
        size <- 2
    }else{
        size <- gp$fontsize
    }
    lnxy <- mapply(build_polygenxy_id.lengths, 
                       starshape=starshape, 
                       phase=phase, SIMPLIFY=FALSE)
    vertices <- unlist(lapply(lnxy, function(x)nrow(x)))
    # ar is the aspect ratio. It can control the 
    # height and width ratio of shapes.
    ar <- match_ar(starshape)
    lxy <- mapply(stretch_rotate_move, p=lnxy,
                  size=size, ar=ar, angle=angle, 
                  x=xv, y=yv, 
                  MoreArgs=list(position.units=position.units,
                                size.units=size.units),
                  SIMPLIFY = FALSE)
    allx <- do.call("unit.c", lapply(lxy, "[[", 1))
    ally <- do.call("unit.c", lapply(lxy, "[[", 2))
    
    grobs <- polygonGrob(allx, ally, id.lengths = vertices, gp = gp, ...)
    return(grobs)
}

deg2rad <- function(x){x * pi / 180}
rad2deg <- function(x){x * 180 / pi}

my_add_num=3
# index of starshape = numbers of edge (n)
starshape_ntab <- c(5, 6, 7, 8,
                    5, 6, 7, 8,
                    8, 8, 3, 4,
                    4, 4, 50, 0,
                    0, 0, 0, 0,
                    0, 3, 3, 6,
                    50, 3, 0, 4,
                    0, 0, 50, 50)
starshape_ntab<-c(starshape_ntab,rep(4,my_add_num))

names(starshape_ntab) <- seq_len(length(starshape_ntab))

match_n <- function(starshape){
    n <- starshape_ntab[match(starshape,names(starshape_ntab))]
    return(unname(n))
}

# index of starshape = aspect ratio (ar) 
starshape_artab <- c(rep(1, 9), 0.5, 1, 0.5, rep(1,12),0.5, 0.18, 1, 1, 1, 1,1,1,1)
starshape_artab<-c(starshape_artab,rep(1,my_add_num))
names(starshape_artab) <- seq_len(length(starshape_artab))

match_ar <- function(starshape){
    ar <- starshape_artab[match(starshape,names(starshape_artab))]
    return(unname(ar))
}

#' @importFrom gridExtra polygon_regular
build_polygenxy_id.lengths <- function(starshape, phase){
    # the edge numbers
    n <- match_n(starshape)
    if (starshape %in% c(1, 2, 3, 4, 9, 10, 14, 22, 24)){
        phase2 <- phase + pi/n
        tmpplxy <- mapply(polygon_regular, 
               phase=c(phase, phase2), 
               n=rep(n, 2), SIMPLIFY=FALSE)
        if (starshape==1){
            tmpplxy[[2]] <- 0.38 * tmpplxy[[2]]
        }else if (starshape==2){
            tmpplxy[[2]] <- 0.556 * tmpplxy[[2]]
        }else if (starshape==3){
            tmpplxy[[2]] <- 0.32 * tmpplxy[[2]]
        }else if (starshape==4){
            tmpplxy[[2]] <- 0.756 * tmpplxy[[2]]
        }else if (starshape==14){
            tmpplxy[[2]] <- 0.35 * tmpplxy[[2]]
        }else if (starshape==22){
            tmpplxy[[2]] <- 0.2 * tmpplxy[[2]]
        }else if (starshape==24){
            tmpplxy[[2]] <- 0.26 * tmpplxy[[2]]
        }else{
            tmpplxy[[2]] <- 0.5 * tmpplxy[[2]]
        }
        tmpplxy <- lapply(tmpplxy,function(x)data.frame(x))
        plxy <- as.matrix(mapply(function(x,y){rbind(x,y)},tmpplxy[[1]],tmpplxy[[2]]))
        colnames(plxy) <- c("x", "y")
    }else if (starshape==16){
        t <- seq(0, 2*pi, by=0.08)
        plxy <- 0.06 *as.matrix(data.frame(x=16 * sin(t)^3,
                     y=13*cos(t) - 5*cos(2*t) - 2*cos(3*t) - cos(4*t) + 2))
    }else if (starshape==13){
        phase <- phase + pi/n
        plxy <- polygon_regular(n=n, phase=phase)
    }else if (starshape==17){
        plxy <- 1.2*matrix(c(-1, -1, 1, 0.5, -0.5, -0.5),nrow=3)
    }else if (starshape==18){
        plxy <- 1.2*matrix(c(-1, 1, 1, 0.5, 0.5, -0.5),nrow=3) 
    }else if (starshape==19){
        plxy <- matrix(c(-1,-1,1,1,-1,-1), nrow=3)
    }else if (starshape==20){
        plxy <- matrix(c(-1,1,1,1,1,-1), nrow=3)
    }else if (starshape==21){
        plxy <- 0.8 * matrix(c(-1, 1, 1, -1,
	                           0.5, 0.5, -0.5, -0.5), nrow=4)
    }else if (starshape==23){
        phase <- phase + pi/n       
        plxy <- 0.8*polygon_regular(n=n, phase=phase)
    }else if (starshape==26){
        phase <- phase + pi/n
        plxy <- 0.7*polygon_regular(n=n, phase=phase)
    }else if (starshape==27){
        plxy <- 0.7*data.frame(x=c(0, -0.25, -0.65, -0.5, -1.1, -0.5, -0.65, 
                                   -0.25, 0, 0.25, 0.65, 0.5, 1.1, 0.5, 0.65, 0.25),
                               y=c(1.4, 0.5, 0.65, 0.25, 0, -0.25, -0.65, -0.5,
                                   -1.4, -0.5, -0.65, -0.25, 0, 0.25, 0.65, 0.5)) 
        plxy <- as.matrix(plxy)
    }else if (starshape==29){
        plxy <- 0.7 * data.frame(x=c(-0.2, -0.2, -1, -1, -0.2, -0.2, 0.2, 0.2, 1, 1, 0.2, 0.2),
                                 y=c(1, 0.2, 0.2, -0.2, -0.2, -1, -1, -0.2, -0.2, 0.2, 0.2, 1))
        plxy <- as.matrix(plxy)
    }else if (starshape==30){
        plxy <- 0.58 * data.frame(x=c(-1,-1.6, 1.6, 1),
                                 y=c(1, -1, -1, 1))
        plxy <- as.matrix(plxy)
    }else if (starshape==31){
        plxy <- polygon_regular(n=n+1, phase=phase)
        endxy <- matrix(apply(plxy[c(nrow(plxy)/2, nrow(plxy)/2+1),], 2, mean), nrow = 1)
        plxy <- plxy[1:(nrow(plxy)/2), ]
        plxy <- rbind(plxy, endxy)
    }else if (starshape == 32){
        phase <- pi/2
        plxy <- 0.85*polygon_regular(n = n+1, phase = phase)
        endxy1 <- matrix(apply(plxy[c(nrow(plxy)/2, nrow(plxy)/2+1),], 2, mean), nrow = 1)
        endxy2 <- matrix(c(max(plxy[,1]), max(plxy[,2]), min(plxy[,1]), max(plxy[,2])), ncol=2, byrow=T)
        plxy <- plxy[seq(nrow(plxy)/2),]
        plxy <- rbind(plxy, endxy1, endxy2)
    }else if (starshape == 33){
        x1<-seq(-1,1,by=0.05)
        y1<-sqrt(1-x1^2)
        y2<-seq(-1,1,by=0.05)
        x2<-sqrt(1-y2^2)
        
        data.frame(x=c(x1,x2+1,-x1,-x2-1),
                   y=c(y1+1,-y2,-y1-1,y2)) %>% as.matrix() -> plxy
    }else if (starshape == 34){
        x1<-seq(-1,1,by=0.5)
        y1<-sqrt(1-x1^2)
        y2<-seq(-1,1,by=0.5)
        x2<-sqrt(1-y2^2)
        
        data.frame(x=c(x1,x2+1,-x1,-x2-1),
                   y=c(y1,-y2,-y1-1,y2)) %>% as.matrix() -> plxy
    }else if (starshape == 35){
        coordinates <- data.frame(
                                  x = c(0.536,0.706,0.855,0.958,1,0.975,
                                        0.885,0.746,0.579,0.409,0.26,0.153,
                                        0.1,0.103,0.155,0.242,0.343,0.437,
                                        0.506,0.537,0.527,0.476,0.396,
                                        0.301,0.206,0.123,0.062,0.024,0.005,
                                        -0.003,-0.013,-0.036,-0.079,-0.146,
                                        -0.233,-0.33,-0.423,-0.496,-0.537,
                                        -0.538,-0.496,-0.419,-0.321,-0.222,
                                        -0.142,-0.101,-0.112,-0.179,-0.298,
                                        -0.455,-0.627,-0.789,-0.917,-0.991,-1,
                                        -0.942,-0.825,-0.669,-0.497,-0.335,
                                        -0.205,-0.123,-0.099,-0.128,-0.199,
                                        -0.296,-0.396,-0.48,-0.531,-0.542,
                                        -0.51,-0.443,-0.354,-0.257,-0.166,
                                        -0.094,-0.044,-0.017,-0.005,0.003,
                                        0.018,0.05,0.106,0.183,0.277,0.373,
                                        0.459,0.517,0.539,0.518,0.457,0.368,
                                        0.266,0.174,0.112,0.095,0.134,
                                        0.229,0.369,0.536),
                                  y = c(-0.164,-0.162,-0.132,-0.078,
                                        -0.009,0.062,0.121,0.158,0.167,0.15,
                                        0.117,0.083,0.064,0.076,0.13,0.23,
                                        0.369,0.534,0.702,0.851,0.956,1,0.976,
                                        0.885,0.742,0.57,0.396,0.25,0.155,
                                        0.128,0.173,0.282,0.438,0.614,
                                        0.782,0.913,0.988,0.995,0.935,0.817,
                                        0.661,0.492,0.331,0.201,0.112,0.069,
                                        0.066,0.091,0.126,0.156,0.167,
                                        0.151,0.108,0.045,-0.027,-0.094,-0.142,
                                        -0.165,-0.161,-0.135,-0.099,-0.07,
                                        -0.065,-0.097,-0.174,-0.295,-0.449,
                                        -0.619,-0.78,-0.909,-0.986,-0.996,
                                        -0.938,-0.818,-0.658,-0.481,-0.317,
                                        -0.194,-0.132,-0.141,-0.22,-0.356,
                                        -0.525,-0.7,-0.853,-0.958,-1,-0.972,
                                        -0.881,-0.742,-0.576,-0.408,-0.261,
                                        -0.15,-0.085,-0.063,-0.075,-0.108,
                                        -0.143,-0.164)
                       )
        as.matrix(coordinates) -> plxy
    }else{
        plxy <- 0.8*polygon_regular(n=n, phase=phase)
    }
    return (plxy)
}

#' @importFrom grid grid.draw
grid.star <- function(x=0.5, y=0.5,
                      starshape=1,
                      angle=0, 
                      phase = 0,
                      gp = NULL,
                      position.units = "npc",
                      size.units="mm",
                      draw = TRUE, vp = NULL, ...){
    sg <- starGrob(x = x, y = y, 
                   starshape = starshape,
                   angle = angle, 
                   gp = gp, 
                   position.units = position.units,
                   size.units = size.units,
                   vp = vp,...)
    if (draw){
        grid.draw(sg)
    }
    invisible(sg)
}

# reference the gridExtra
stretch_rotate_move <- function(p, size, 
                                ar, angle, x, 
                                y, position.units, size.units){
    central <- size * p %*%
    diag(c(sqrt(ar), 1/sqrt(ar))) %*%
         rbind(c(cos(angle), -sin(angle)),
         c(sin(angle),  cos(angle)))
    list(x = unit(central[,1], size.units) + unit(x, position.units),
	 y = unit(central[,2], size.units) + unit(y, position.units))
}

