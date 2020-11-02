#' Model age-proxy-uncertainty (robust charcoal accumulation rate / CHAR) based on Dietze et al., 2019.
#' 
#' This function models the uncertainty imposed on a proxy data set with 
#' continuously sampled proxy values and their respective proxy error in stratigraphic order. 
#' cm depth and age (+ respective error) should mark the upper and lower 
#' boundary of a sample (not the center point). 
#' 
#' The approch consists of several succesive steps. First, a density 
#' estimate of the proxy flux (i.e., reciprocal of sedimentation rate) is 
#' calculated based on Monte Carlo simulations. 
#' Then an age density is calculated by combining the mean ages and sd of ages of 
#' upper and lower boundary of a samples.
#' Bins of a certain time resolution are created and in a combined MC run, 
#' the density distributions of fluxes are randomnly sampled from all samles that have 
#' an age probability in this age bin.
#' Finally, the resulting distributions per bins are characterized by median and interquartile ranges.
#' 
#' @param data \code{Data frame}, input data set. Must contain the following 
#' data in the specified order: \code{depth (cm)}, \code{age (a)}, 
#' \code{age uncertainty (a)}, \code{proxy value} and 
#' \code{proxy uncertainty}.
#' 
#' @param n \code{Numeric} value, number of Monte Carlo runs. Default is 
#' \code{100}.
#' 
#' @param resolution_density \code{Numeric} value, temporal resolution of the 
#' density estimates during calculation. Default is \code{0.1} years.
#' 
#' @param resolution_out \code{Numeric} value, temporal resolution of the 
#' output data set. Default is \code{3} years.
#' 
#' @param n_density \code{Numeric} value, number of values for density 
#' estimates. Default is \code{5000}. 
#' 
#' @param scale \code{Logical} value, option to scale (z-transform the data),
#' default is \code{FALSE}.
#' 
#' @param plot \code{Logical} value, option to create a plot.
#' 
#' @return something to still describe.
#' 
#' 
#' 
CHARrobust <- function(
  data, 
  n = 100,
  resolution_density = 0.1,
  resolution_out = 3,
  n_density = 5000,
  scale = FALSE,
  BP = FALSE, # for dates in age relative to 1950 AD
  plot = TRUE,
  ...
) {
  
  ## check input data
  if(sum(names(data) == c("depth", 
                          "age", 
                          "age_error", # should be 1 sd
                          "proxy", 
                          "proxy_error")) < 5) {
    
    stop("Input data structure is not correct! See function documentation.")
  }
  if(BP == TRUE) 
    {data$age <- 1950-data$age}
  
  ## PART 1 - Generate unit deposition time -----------------------------------
  
  t_unit_raw <- lapply(X = 1:n, FUN = function(i, data) {
    
    ## draw random age estimates assuming normal distribution
    age_i <- rnorm(n = nrow(data),
                   mean = data$age,
                   sd = data$age_error )
    
    ## calculate age differences (remove negative sign, if working in AD/BC)
    age_diff <- -c(NA,diff(age_i))
    
    ## set inverted ages to NA
    age_diff[age_diff < 0] <- NA
    
    ## calculate unit deposition time (= 1/sedimentation rate)
    t_unit <-  age_diff / diff(c(NA, data$depth)) 
    
    ## return result
    return(list(t_unit = t_unit,
                age = age_i))
    
  },
  data = data)
  
  ## convert list to matrices
  t_unit <- do.call(cbind, lapply(X = t_unit_raw, FUN = function(X) {
    X$t_unit
  }))
  
  t_age <- do.call(cbind, lapply(X = t_unit_raw, FUN = function(X) {
    X$age
  }))
  
  ## calculate mean unit deposition times
  t_unit_mean <- apply(X = t_unit,  
                       MARGIN = 1,
                       FUN = mean,
                       na.rm = TRUE)
  
  ## calculate sd unit deposition times
  t_unit_sd <- apply(X = t_unit,
                     MARGIN = 1,
                     FUN = sd,
                     na.rm = TRUE)
  
  ## PART 2 - Generate sample wise flux density distributions -----------------------------------
  flux_raw <- 
    lapply(X = 1:n, FUN = function(X, data, t_unit_mean, t_unit_sd) {
      
      
      ## draw random unit deposition times assuming normal distribution
      sed_i <- rnorm(n = nrow(data) - 1,
                     mean = t_unit_mean[2:nrow(data)],
                     sd = t_unit_sd[2:nrow(data)])
      
      ## identify NAs
      i_na <- is.na(data$proxy)
      
      ## replace NAs by zero
      data_0 <- data
      data_0[i_na, 4:5] <- 0
      
      ## draw random proxy value estimates assuming normal distribution
      proxy_i <- rnorm(n = nrow(data_0) - 1, 
                       mean = data_0$proxy[-1], 
                       sd = data_0$proxy_error[-1])
      
      ## replace zero data by NA
      proxy_i[i_na[-1]] <- NA
      
      ## calculate age differences 
      proxy_flux <- proxy_i / sed_i
      
      ## set inverted ages to NA
      proxy_flux[proxy_flux < 0] <- NA
      
      ## return result
      return(proxy_flux)
      
    },
    data = data, 
    t_unit_mean = t_unit_mean, 
    t_unit_sd = t_unit_sd)
  
  ## convert list to matrix
  flux <- do.call(cbind, flux_raw)
  
  ## calculate empiric density function
  flux_density <- apply(X = flux,
                     MARGIN = 1,
                     FUN = function(x) {
                       
                       d <- try(density(x = x, 
                                        na.rm = TRUE,
                                        from = 0, 
                                        to = quantile(x = x,
                                                      probs = 0.99,
                                                      na.rm = TRUE),
                                        n = n_density), 
                                silent = TRUE)
                       
                       if(class(d) == "try-error") {
                         
                         d <- NA
                       }
                       
                       return(d)
                     })
  
  
  ## PART 3 - Generate sample wise age density distribution  -----------------------------------
  ## sort age and sample order
  data_order <- data[order(data$age),]
  
  da_min <- data_order$age_error[data_order$age == min(data_order$age)]
  da_max <- data_order$age_error[data_order$age == max(data_order$age)]
  
  ## convert age and age uncertainty to pairwise list
  age_info <- as.list(as.data.frame(rbind(data_order$age, 
                                          data_order$age_error))) 
  
  ## define age density age vector
  age_index <- seq(from = min(data_order$age) - 5 * da_min, 
                   to = max(data_order$age) + 5 * da_max, 
                   by = resolution_density)
  
  ## calculate densities for each age value
  density_raw <- lapply(X = age_info, FUN = function(age_info, age_index) {
    
    dnorm(x = age_index, mean = age_info[1], sd = age_info[2])
  },
  age_index = age_index)
  
  ## generate output data set
  ii <- seq(from = 1, to = length(age_index))

  age_density <- vector(mode = "list",
                        length = length(age_info) - 1)

  ## generate all density estimates
  for(i in 1:(length(age_info) - 1)) {
    
    i_l <- ii[density_raw[[i]] == max(density_raw[[i]])][1]
    i_u <- ii[density_raw[[i + 1]] == max(density_raw[[i + 1]])][1]
    
    ## define lower tail age vector
    d_l <- density_raw[[i]][1:i_l]
    
    ## define upper tail age vector
    d_u <- density_raw[[i + 1]][i_u:length(ii)]
    
    ##normalise density vectors
    d_l_n <- d_l / max(d_l)
    d_u_n <- d_u / max(d_u)
    
    ## define central part density vector
    d_m <- rep(1, i_u - i_l -1)
    
    ## merge age and density vector
    d_merged <- c(d_l_n, d_m, d_u_n)
    
    ## normalise merged density vector
    d_merged_n <- d_merged / sum(d_merged)
    
    ## generate and assign output data set
    age_density[[i]] <- d_merged_n
  }
  
  ## revert order to original order
  age_density <- age_density[length(age_density):1]
  
  ## PART 4 - Combine age and flux density estimates -----------------------------------
  
   ## calculate fluxes
  flux_combined <- vector(mode = "list", 
                          length = length(age_density))
  
  for(i in 1:length(age_density)) {
    
    age <- sample(x = age_index, 
                  size = n, 
                  replace = TRUE, 
                  prob = age_density[[i]])
    
    if(is.na(data$proxy[i + 1])) {
      
      proxy <- rep(NA, times = n)
    } else {
      
      if(scale == TRUE) {
        
        flux_density[[i]]$x <- scale(flux_density[[i]]$x)
      }
      
      proxy <- sample(x = flux_density[[i]]$x, 
                      size = n, 
                      replace = TRUE, 
                      prob = flux_density[[i]]$y)
    }
    
     flux_combined[[i]]<- data.frame(age = age,
                                     proxy = proxy)
  }
  
  x_raw <- unlist(lapply(X = flux_combined,
                         FUN = function(x) x$age))
  
  y_raw <- unlist(lapply(X = flux_combined, 
                         FUN = function(x) x$proxy))

  bins <- seq(from = max(data$age), 
              to = min(data$age), #
              by = -resolution_out)
  
  data_out <- matrix(nrow = length(bins) - 1, 
                   ncol = 5)
  
  colnames(data_out) <- c("q_10","q_25", "q_50", "q_75","q_90")
  
  for(i in 1:(length(bins) - 1)) {
    
    x_i <- y_raw[x_raw <= bins[i] & x_raw > bins[i + 1]]
    
    data_out[i,] <- quantile(x = x_i, 
                             probs = c(0.1,0.25, 0.5, 0.75, 0.9), 
                             na.rm = TRUE)
  }
  
  if(plot == TRUE) {

    extraArgs <- list(...)
    
    if ("main" %in% names(extraArgs)) {
      main <- extraArgs$main
    }
    else {
      main <- ""
    }
    
    if ("xlab" %in% names(extraArgs)) {
      xlab <- extraArgs$xlab
    }
    else {
      xlab <- "Time"
    }
    
    if ("ylab" %in% names(extraArgs)) {
      ylab <- extraArgs$ylab
    }
    else {
      ylab <- "Proxy value"
    }
    
    if ("xlim" %in% names(extraArgs)) {
      xlim <- extraArgs$xlim
    }
    else {
      xlim <- range(bins[-1])
    }

    if ("ylim" %in% names(extraArgs)) {
      ylim <- extraArgs$ylim
    }
    else {
      ylim <- range(data_out[,2:4], na.rm = T)
    }
    
    plot(NA, 
         xlim = xlim,
         ylim = ylim,
         xlab = xlab,
         ylab = ylab,
         main = main)
    
    polygon(x = c(bins[-1], rev(bins[-1])), 
            y = c(data_out[,2], rev(data_out[,4])),
            col = "grey", 
            border = NA)
    
    lines(x = bins[-1],
          y = data_out[,3])
    
  }
  
  data_out <- data.frame(t_lower = bins[-length(bins)],
                         t_med = apply(cbind(bins[-length(bins)],bins[-1]),1,median),
                         t_upper = bins[-1],
                         q_10 = data_out[,1],
                         q_25 = data_out[,2],
                         q_50 = data_out[,3],
                         q_75 = data_out[,4],
                         q_90 = data_out[,5])
  
  return(data_out)
}


