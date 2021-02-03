### CharcoalFireReconstruction v1.0

### Script for the evaluation of sedimentary charcoal records

### Description
# This script can be used to assess fire signals within sedimentary charcoal records in both a "classic" and a "robust" approach.
# The classic approach (classic CHAR) calculates charcoal accumulation rates (CHAR) from charcoal counts per depth and unit of time. It separates 
# the influx record (calculated with the "paleofire" R package, Blarquez et al., 2014) into background and peak components (similar to Higuera et al., 2009, 2010), 
# obtains signal to noise index values for the peak component (after Kelly et al., 2011) and identifies peaks above a previously established global threshold as 
# fire episodes. The robust approach (robust CHAR) follows a similar overall procedure, but also includes a Monte Carlo-based method to incorporate age, proxy 
# and user input related uncertainties into the resulting timeseries (see Dietze et al., 2019). 
# Please refer to the readme file and references therein for more information.

### Included steps:
# 1. Run setup
# 2. Classic CHAR: Separation of background and peak components
# 3. Classic CHAR: Calculation of the signal to noise index (SNI)
# 4. Robust CHAR: Estimating CHAR with added age and proxy uncertainties
# 5. Robust CHAR: Separation of background and peak components
# 6. Generation of .txt files specifying run settings and results
# 7. Generation of a compiled plot of classic and robust CHAR

### Important:
# See the readme file for full references of methods and code used in this script. Please cite references if this script is used for a publication. 
# Please use this script at your own risk.



# # #
# 1.#   Run setup
# # #

    # Set working directory (.../CharcoalAnalysisR)
setwd("M:/PhD/Khamra/R/Combined script/Script for Github/Revised version/CharcoalAnalysisR-master")

    # Load packages
library(locfit)
library(mixtools)
library(paleofire)
library(stats)
library(zoo)

    # Load functions
source("./Function_RobustCHAR_Dietze2019.R")
source("./Function_SNI_Kelly2011.R")

    # Name of this record / run (a corresponding folder will be created in the wd to save upcoming files, if not already created) 
record_name = "RandomRecord" # Should be the name of the .csv input file and the folder it is located in, e.g. ./Records/record_name/record_name.csv
run_name = "_RunGLOBAL1" # Additional naming resource if multiple runs of the same record with different settings should be done
dir.create(path = paste0("./Records/", record_name, sep = ""))
   
    # Input file named after record_name, with these columns: depth_top, depth_bot, age_top, age_bot, age_uncert, vol, proxy, proxy_uncert
data <- read.csv(paste0("./Records/", record_name, "/", record_name, ".csv"), header = TRUE, sep = ",")

    # Define the proxy of this run (i.e. charcoal counts)
proxy = data[7]

    # Define window width to use in this run (as percentage of total record length; for LOWESS and SNI)
win_width = 0.25

    # Use a global (TRUE) or local threshold (FALSE)? -> Local threshold was used for tests only and should be regarded as experimental
use_global_tresh = TRUE

    # Parameters for robust CHAR (see Dietze et al., 2019)
n = 2000 # Number of MC runs, the more the longer it takes and the smoother the curve
range_s <- c(5, 25) # Percent of points to include in smoothing
n_s <- 1000 # Number of smoothings



# # #
# 2.#   Classic CHAR: Separation of background and peak components (see Higuera et al., 2009, 2010; Dietze et al., 2019)
# # #

df <- data.frame(depthtop = data[1],
                 depthbot = data[2],
                 agetop = data[3],
                 agebot = data[4],
                 volume = data[6], # In cm^3 
                 CHAR = proxy) # Absolute charcoal counts
df <- as.matrix(na.omit(df))

CHAR <- paleofire::pretreatment(df[,1:5],df[,6]) # Calculating CHAR using median time resolution (default) with the paleofire package (Blarquez et al., 2014)
x <- CHAR$ybpI # Interpolated ages
y <- CHAR$accI # Interpolated CHAR

    # Classic CHAR background component (relying on window width input above)
win <- (max(x)-min(x))*win_width
locbootA <- locfit(y ~ lp(x, deg = 1, h = win/2), 
                   maxk = 500, 
                   family = "qgauss")
predbootA <- predict(locbootA, newdata = x, se.fit = TRUE)

plot(x, y, type = "l", lty = 2)
lines(x, predbootA$fit)

    # Classic CHAR peak component
peak <- as.data.frame(cbind(peak = y-predbootA$fit,
                            x = x))

peak1 <- peak$peak[peak$peak > 0] # Consider only positive peaks

plot(peak$peak, type = "h")
hist(peak1, breaks = 20)


if (use_global_tresh == TRUE) { # Global threshold
  # Gaussian mixture model to differentiate between char fire and char noise, the two constituents of the peak component
  gm <- normalmixEM(peak1, k=2)
  plot(gm,which=2)
  
  # Usually the treshold is defined as the 99th percentile of the distribution of the noise component
  global_treshold <- gm$mu[1]+qnorm(0.99)*gm$sigma[1];global_treshold
  
  # Identifying positive peaks above the treshold and mark as fire episodes
  fire_pos <- which(peak$peak > global_treshold)
  fire_bp <- x[fire_pos]
  fire_ad <- 1950-fire_bp  #fire dates in AD
  
  treshold = global_treshold #for further use in this script
  peak = peak$peak
  
} else { # Local threshold, experimental
  # Function for calculation of threshold in moving window
  locthresh <- function(x) {
    gm <- normalmixEM(x, k = 2)
    treshold <- gm$mu[1]+qnorm(0.99)*gm$sigma[1]
    return(treshold)
  }
  
  # Window width of local treshold function
  thresh_win = length(peak1)*win_width
  
  tresh_raw <- rollapply(peak1, thresh_win, locthresh, align = "center", partial = TRUE)
  
  # Plot smoothed local threshold
  sep <- as.data.frame(cbind(peak1,
                             x = peak$x[peak$peak > 0],
                             tresh_raw))
  
  tresh_loess <- loess(sep$tresh_raw ~ sep$x, span = 0.7) # Smoothing of raw local treshold values, span can be adjusted
  local_treshold <- predict(tresh_loess)
  
  plot(sep$peak1 ~ sep$x, type = "h")
  lines(local_treshold ~ sep$x, col = "red", lwd = 2)
  
  # Identifying positive peaks above the treshold and mark as fire episodes
  sep <- cbind(sep,
               local_treshold)
  peak <- merge(sep, peak, all = TRUE)
  peak$local_treshold <- na.spline(peak$local_treshold, na.rm = FALSE) # Interpolation of threshold values to all sample depths, only used for visualization
  
  fire_pos <- which(peak$peak >= peak$local_treshold)
  fire_bp <- x[fire_pos]
  fire_ad <- 1950-fire_bp
  
  treshold = peak$local_treshold # For further use in this script
  peak = peak$peak
}



# # #
# 3.#   Classic CHAR: signal-to-noise index SNI (see Kelly et al., 2011)
# # #

    # Load data
SNIdat <- as.matrix(na.omit(data.frame(age = x,
                                       CHAR = peak,
                                       tresh = treshold)))

    # Calculate SNIs  
SNI <- CharSNI(SNIdat, win)



# # #
# 4.#   Robust CHAR (see Dietze et al., 2019)
# # #

    # Create input dataframe for robust CHAR
input <- data.frame(depth = data$depth_top, 
                   age = data$age_top,
                   age_error = data$age_uncert,
                   proxy = proxy/data$vol, # Using concentration (# cm^-3)
                   proxy_error = data$proxy_uncert) # Error/uncertainty could be measurement error from duplicates

    # Output resolution
res_out = median(diff(input$age)) # Median resolution (fit to individual record)

    # Running the model (for more extensive records this might take a while)
char.uncert <- CHARrobust(input, n = n, resolution_out = res_out, 
                          xlab="Age (CE)",
                          ylab = "CHAR full age uncertainty (# cm^-2 yr^-1)",
                          xlim = range(1950-data[3]),
                          BP = TRUE,
                          scale = FALSE)



# # #
# 5.#   Robust CHAR: Separation of background and peak components (see Dietze et al., 2019)
# # #

    # Separating background and peak components 
proxMC <- char.uncert
cutoff <- c(5, nrow(proxMC)-5)
ageMC <- proxMC$t_med

    # Robust CHAR background component
char_back <- do.call(cbind, lapply(X = 1:n_s, FUN = function(X, proxMC, n_s, range_s) {
  n_data <- nrow(proxMC)
  s_i <- runif(n = 1,
               min = range_s[1],
               max = range_s[2])
  span <- s_i / n_data
  
  lowess(x = scale(proxMC$q_50), f = span)$y
}, proxMC = proxMC, n_s = n_s, range_s = range_s))

char_back_stats <- rbind(q1 = apply(X = char_back,
                                    MARGIN = 1,
                                    FUN = mean),
                         q2 = apply(X = char_back,
                                    MARGIN = 1,
                                    FUN = sd))

y_back = c(char_back_stats[1,]+2*char_back_stats[2,],
           rev(char_back_stats[1,]-2*char_back_stats[2,]))

    # Robust CHAR peak component
char_peak <- do.call(cbind, lapply(X = 1:n_s, FUN = function(X, proxMC, n_s, range_s) {
  n_data <- nrow(proxMC)
  s_i <- runif(n = 1,
               min = range_s[1],
               max = range_s[2])
  span <- s_i / n_data
  scale(proxMC$q_50) - lowess(x = scale(proxMC$q_50), f = span)$y
}, proxMC = proxMC,n_s = n_s,range_s = range_s))

char_peak_stats <- rbind(q1 = apply(X = char_peak,
                                    MARGIN = 1,
                                    FUN = mean),
                         q2 = apply(X = char_peak,
                                    MARGIN = 1,
                                    FUN = sd))

y_peak = c(char_peak_stats[1,]+2*char_peak_stats[2,],
           rev(char_peak_stats[1,]-2*char_peak_stats[2,]))



# # #
# 6.#   Generation of .txt files specifying run settings and results
# # #

    # Results file
run_res <- rbind(c("--- Temporal resolution of the record ---", "\n",
                   "Median temporal resolution:", median(diff(data$age_top)), "\n",
                   "Standard deviation (1 sig):", sd(diff(data$age_top), na.rm = TRUE), "\n",
                   "Max temporal resolution:", max(diff(data$age_top)), "\n",
                   "Min temporal resolution:", min(diff(data$age_top)), "\n",
                   "\n",
                   "--- Signal-to-noise index (SNI) ---", "\n",
                   "SNI mean:", mean(SNI$SNI), "\n",
                   "SNI med:", median(SNI$SNI), "\n",
                   "SNI max:", max(SNI$SNI), "\n",
                   "SNI min:", min(SNI$SNI), "\n",
                   "Percentage of SNI > 3:", (sum(SNI$SNI > 3) / length(SNI$SNI)) * 100, "\n",
                   "Percentage of SNI < 3:", (sum(SNI$SNI < 3) / length(SNI$SNI)) * 100, "\n",
                   "\n",
                   "--- CHAR and fire episodes / fire return intervals (FRI) ---", "\n",
                   "Mean CHAR of the record:", mean(CHAR$accI),"\n",
                   "Standard deviation (1 sig):", sd(CHAR$accI), "\n",
                   "Number of fire episodes:", length(fire_ad), "\n",
                   "FRI mean:", mean(diff(fire_ad)*-1), "\n",
                   "Standard deviation (1 sig):", sd(diff(fire_ad)*-1), "\n",
                   "FRI med:", median(diff(fire_ad)*-1), "\n",
                   "FRI max:", max(diff(fire_ad)*-1), "\n",
                   "FRI min:", min(diff(fire_ad)*-1), "\n",
                   "Identified fire years:", fire_ad, "\n",
                   "Fire years with SNI > 3:", length(which(SNI$SNI[which(peak>treshold)]>3)), "\n",
                   "Precentage of fire years with SNI > 3:", (length(which(SNI$SNI[which(peak>treshold)]>3))/length(fire_ad))*100, "\n",
                   "Fire frequency per 100 years:", (length(fire_ad)/max(x)*100)))

cat(run_res, file = paste0("./Records/", record_name, "/", record_name, run_name, "_Results.txt"), sep = "\n", append = FALSE)

    # Settings file
run_sett <- rbind(c("--- Classic CHAR inputs ---",  "\n",
                    "Number of samples:", length(proxy$proxy), "\n",
                    "Number of non-zero values within these samples:", length(which(proxy!=0)), "\n",
                    "Median temporal resolution the record was interpolated to:", CHAR$yrInterp, "\n",
                    "Window width chosen:", win, "\n",
                    "Window width accounts for", win_width,"% of total record length", "\n",
                    "Applied global threshold:", use_global_tresh, "\n",
                    "(Mean) Treshold:", mean(treshold), "\n",
                    "\n",
                    "--- Robust CHAR inputs ---",  "\n",
                    "Output resolution:", res_out, "\n",
                    "Number of MC runs:", n, "\n",
                    "range_s:", range_s, "\n",
                    "n_s:", n_s, "\n",
                    "Cutoff:", cutoff))

cat(run_sett, file = paste0("./Records/", record_name, "/", record_name, run_name, "_Settings.txt"), sep = "\n", append = FALSE)



# # #
# 7.#   Generation of a complied plot of classic and robust CHAR
# # #

# Create output plot to folder
pdf(file = paste0("./Records/", record_name, "/", record_name, run_name, "_plotted.pdf", sep = ""), width = 12, height = 14)

    # Plot classic CHAR peak component and SNI
par(mfrow = c(6,1), mar = c(0,5,2,5))
plot(1950-x, peak, xaxt = "n", xlab = "", ylab = "particles cm^-2 yr^-1", las = 1, type="h", lend = 1, lwd = 3, 
     col = ifelse(peak < treshold,'darkgrey','grey30'), yaxt = "n")
abline(h = 0, col = "gray35")
title("Classic CHAR peak component", adj = 0.01, line = -1)
title(paste0("Run: ", record_name, run_name), adj = 0.01, line = 1)
axis(2, las = 1)
if (use_global_tresh == TRUE) {
  abline(h = treshold, lwd = 2, col = "black", lty = "longdash", lend = 1) # Draw global threshold
} else {
  lines(1950-x, treshold, lwd = 2, col = "black", lty = "longdash", lend = 1) # Draw local threshold
}
par(mar = c(0,5,0,5))
plot(1950-x, SNI$SNI, type="l", lwd = 2, ylab = "", las = 1, xlab = "", xaxt = "n", yaxt = "n")
title("SNI", adj = 0.01, line = -1)
abline(h = 3, lwd = 2, col = "tomato3", lty = "solid") # Draw SNI cutoff value
axis(4, las = 1)
mtext("Index", side = 4, line = 3, cex = 0.7)

    # Classic CHAR (sum) with background component and identified fire episodes
plot(1950-x, y, xlim = range(1950-data$age_top), xlab = "", ylim = range(CHAR$accI), ylab = "particles cm^-2 yr^-1", las = 1, type = "l",lty = 1, lwd = 2, 
     col = "black", xaxt = "n", yaxt = "n")
axis(2, las = 1)
lines(1950-x, predbootA$fit, type="l", lwd = 3, col = "deepskyblue3")
ypos1 <- c(rep(1.0*max(CHAR$accI),length(fire_ad))) # Positions of points and segments for each fire episode on the y-axis
ypos2 <- c(rep(-1.0*min(CHAR$accI),length(fire_ad)))
SNIplot <- data.frame(fire_ad = fire_ad,
                      SNI = SNI$SNI[which(peak>treshold)])
segments(x0 = SNIplot$fire_ad, y0 = ypos1, x1 = SNIplot$fire_ad, y1 = ypos2, lwd = 1, 
         col = ifelse(SNIplot$SNI > 3, "tomato3", "darkgrey"),lty = "solid") # Fire episodes with SNI >3 are drawn in red, otherwise in grey
points(x = SNIplot$fire_ad, y = ypos1, pch = 25 , bg = "0", 
       col=ifelse(SNIplot$SNI > 3, "tomato3", "darkgrey"), lwd = 3)
title("Classic CHAR", adj = 0.01, line = -1)

    # Robust CHAR background component
plot(NA, xlim = range(1950-data$age_top), ylim = range(y_back),
     xlab="", ylab = "", las = 1, xaxt = "n", yaxt = "n")
axis(4, las = 1)
polygon(x = c(ageMC, rev(ageMC)), y = y_back,
        col = "grey",  border = NA)
lines(ageMC[],char_back_stats[1,], lwd = 2)
abline(h = 0, col = "gray35")
title("Robust CHAR background component", adj = 0.01, line = -1)
mtext("particles cm^-2 yr^-1", side = 4, line = 3, cex = 0.7)

    # Robust CHAR peak component
plot(NA,  xlim = range(1950-data$age_top), ylim = range(y_peak),
     xlab="", ylab = "zscores", las = 1, xaxt = "n")
polygon(x = c(ageMC[], rev(ageMC[])), y = y_peak,
        col = "grey",  border = NA)
lines(ageMC[],char_peak_stats[1,], lwd = 2)
abline(h = 0, col = "gray35")
title("Robust CHAR peak component", adj = 0.01, line = -1)

    # Robust CHAR (sum)
par(mar = c(5,5,0,5))
plot(NA, ylim = c(min(char.uncert$q_25), max(range(char.uncert$q_75))),  xlim = range(1950-data$age_top), xlab="Age (CE)", 
     ylab = "", las = 1, yaxt = "n")
axis(4, las = 1)
polygon(x = c(char.uncert$t_med[], rev(char.uncert$t_med[])), y = c(char.uncert$q_25, rev(char.uncert$q_75)),
        col = "grey",  border = NA)
lines(char.uncert$t_med[], char.uncert$q_50, lwd = 2)
title("Robust CHAR", adj = 0.01, line = -1)
mtext("particles cm^-2 yr^-1", side = 4, line = 3, cex = 0.7)

dev.off()

