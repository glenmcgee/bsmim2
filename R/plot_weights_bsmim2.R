#' Plot weights (rho, theta, theta*) for bsmim
#' 
#'
#' @param object An object of class bsmim
#'
#' @return a list of boxplots.
#' @improt ggplot2 tidyverse
#' @export
#'

plot_weights_bsmim2 <- function(object){
  
  if(class(object)!="bsmim"){
    stop("input must be of class bsmim")
  }
  
  ## theta
  plots_theta <- list()
  for(mm in 1:length(object$theta)){                 ## loop over m (indices)

    dfplot <- data.frame(object$theta[[mm]])
    colnames(dfplot) <- 1:ncol(dfplot)
    dfplot <- gather(dfplot,"l","theta")
    if(ncol(object$theta[[mm]])==1){dfplot$l <- ""}
    
    p <- ggplot(dfplot, aes(x=l, y=theta))+
      geom_boxplot()+##color=alpha("black",0.8)
      geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
      scale_y_continuous(expression(paste("Component Weights (",theta,")")),limits=c(-1,1)) + 
      scale_x_discrete("Exposure Component") + 
      ggtitle("") + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
    plots_theta[[mm]] <- p
  }
  
  ## thetastar
      ## first get the max and min
      dflist <- list()
      for(mm in 1:length(object$theta)){                 ## loop over m (indices)
        dflist[[mm]] <- data.frame(object$theta[[mm]])*sqrt(as.matrix(object$rho)[,mm])
      }
      maxy <- max(sapply(dflist,max))
      miny <- min(sapply(dflist,min))
  
  plots_thetastar <- list()
  for(mm in 1:length(object$theta)){                 ## loop over m (indices)
    
    dfplot <- data.frame(object$theta[[mm]])*sqrt(as.matrix(object$rho)[,mm])
    colnames(dfplot) <- 1:ncol(dfplot)
    dfplot <- gather(dfplot,"l","thetastar")
    if(ncol(object$theta[[mm]])==1){dfplot$l <- ""}

    p <- ggplot(dfplot, aes(x=l, y=thetastar))+
      geom_boxplot()+##color=alpha("black",0.8)
      geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
      scale_y_continuous(expression(paste("Exposure Importance (",theta,"*)")),limits=c(miny,maxy)) + 
      scale_x_discrete("Exposure Component") + 
      ggtitle("") + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
    plots_thetastar[[mm]] <- p
  }
  
  ## thetaPOS (i.e. weights as proportions, when theyre all positive)
  plots_thetaPOS <- list()
  if(!is.null(object$thetaPOS)){
    for(mm in 1:length(object$thetaPOS)){                 ## loop over m (indices)
      
      dfplot <- data.frame(object$thetaPOS[[mm]])
      colnames(dfplot) <- 1:ncol(dfplot)
      dfplot <- gather(dfplot,"l","thetaPOS")
      if(ncol(object$thetaPOS[[mm]])==1){dfplot$l <- ""}
      
      p <- ggplot(dfplot, aes(x=l, y=thetaPOS))+
        geom_boxplot()+##color=alpha("black",0.8)
        geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
        scale_y_continuous(expression(paste("Component Weights Proportions")),limits=c(0,1)) + 
        scale_x_discrete("Exposure Component") + 
        ggtitle("") + 
        theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
      
      plots_thetaPOS[[mm]] <- p
    }
  }
  
  ## w (raw exposure weights)
  plots_w <- list()
  for(mm in 1:length(object$w)){                 ## loop over m (indices)
    
    dfplot <- data.frame(object$w[[mm]])
    if(!is.null(colnames(object$x[[mm]]))){
      colnames(dfplot) <- colnames(object$x[[mm]])
    }
    dfplot <- gather(dfplot,"l","w")
    if(ncol(object$w[[mm]])==1){
      colnames(dfplot) <- "w"
      dfplot$l <- ""
    }
    
    p <- ggplot(dfplot, aes(x=l, y=w))+
      geom_boxplot()+##color=alpha("black",0.8)
      geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
      scale_y_continuous(expression(paste("Component Weights (",w[m],")")),limits=c(-1,1)) + 
      scale_x_discrete("Exposure Component",limits=colnames(object$x[[mm]])) + 
      ggtitle("") + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
    plots_w[[mm]] <- p
  }
  
  ## wstar
  ## first get the max and min
  dflist <- list()
  for(mm in 1:length(object$w)){                 ## loop over m (indices)
    dflist[[mm]] <- data.frame(object$w[[mm]])*sqrt(as.matrix(object$rho)[,mm])
  }
  maxy <- max(sapply(dflist,max))
  miny <- min(sapply(dflist,min))
  
  plots_wstar <- list()
  for(mm in 1:length(object$w)){                 ## loop over m (indices)
    
    dfplot <- data.frame(object$w[[mm]])*sqrt(as.matrix(object$rho)[,mm])
    if(!is.null(colnames(object$x[[mm]]))){
      colnames(dfplot) <- colnames(object$x[[mm]])
    }
    dfplot <- gather(dfplot,"l","wstar")
    if(ncol(object$w[[mm]])==1){
      colnames(dfplot) <- "wstar"
      dfplot$l <- ""
    }
    
    p <- ggplot(dfplot, aes(x=l, y=wstar))+
      geom_boxplot()+##color=alpha("black",0.8)
      geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
      scale_y_continuous(expression(paste("Exposure Importance (",w[m],"*)")),limits=c(miny,maxy)) + 
      scale_x_discrete("Exposure Component",limits=colnames(object$x[[mm]])) + 
      ggtitle("") + 
      theme_bw() +
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
    
    plots_wstar[[mm]] <- p
  }
  
  ## wPOS (raw exposure weights) as proportions (when positive)
  plots_wPOS <- list()
    if(!is.null(object$wPOS)){
    for(mm in 1:length(object$wPOS)){                 ## loop over m (indices)
      
      dfplot <- data.frame(object$wPOS[[mm]])
      if(!is.null(colnames(object$x[[mm]]))){
        colnames(dfplot) <- colnames(object$x[[mm]])
      }
      dfplot <- gather(dfplot,"l","w")
      if(ncol(object$w[[mm]])==1){
        colnames(dfplot) <- "w"
        dfplot$l <- ""
      }
      
      p <- ggplot(dfplot, aes(x=l, y=w))+
        geom_boxplot()+##color=alpha("black",0.8)
        geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
        scale_y_continuous(expression(paste("Component Proportions")),limits=c(0,1)) + 
        scale_x_discrete("Exposure Component",limits=colnames(object$x[[mm]])) + 
        ggtitle("") + 
        theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
      
      plots_wPOS[[mm]] <- p
    }
  }
  
  ## rho
  dfplot <- data.frame(as.matrix(object$rho))
  colnames(dfplot) <- 1:ncol(dfplot)
  dfplot <- gather(dfplot,"l","rho")

  p <- ggplot(dfplot, aes(x=l, y=rho))+
    geom_boxplot()+##color=alpha("black",0.8)
    geom_hline(yintercept=0,linetype=2,alpha=0.8)+ 
    scale_y_continuous(expression(paste("Index Importance (",rho,")")))+#,limits=c(-1,1)) + 
    scale_x_discrete("Exposure Index") + 
    ggtitle("") + 
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 
  
  plots_rho <- p
  
  
  

  
  return(list(theta=plots_theta,
              thetastar=plots_thetastar,
              thetaPOS=plots_thetaPOS,
              rho=plots_rho,
              w=plots_w,
              wstar=plots_wstar,
              wPOS=plots_wPOS))
}





