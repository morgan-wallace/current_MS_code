##### re-doing the last code... #####

##### installing the pre-rec libraries #####

library(dplyr)
library(tidyr)
 
###

## install package ##  

#install.packages("staRdom")

##

#install.packages("devtools") # Run this only, if devtools is not installed already.


###

##### Run just the libraries #####

library("dplyr")
library("tidyr")
library("staRdom")

##### end libaries

##### Cores #####

cores <- detectCores(logical = FALSE)


##### Raw Data #####


### eem data ###

system.file("extdata/EEMs", package = "staRdom")


### absorb data ###

system.file("extdata/absorbance", package = "staRdom")


### other raw data ,  table with data on diluted samples included ###

system.file("extdata/metatable_dreem.csv", package = "staRdom")


### calculated PARAFAC model ###

data(pf_models)

##### end of raw import data #####


##### import eems #####
folder <- system.file("extdata/EEMs/", package = "staRdom") # folder containing example EEMs 

eem_list <- eem_read(folder, recursive = TRUE, import_function = eem_csv) 

# in case you use your own data, just replace folder by a path. e.g. "C:/folder/another folder" and change import_function according to instrument.

#eem_list <- eem_read(folder, import_function = "cary")

# Horiba Fluoromax-4 ("fluoromax4")
#  generic csv files (eem_csv when excitation wavelengths are in columns).

##### end eem import #####


##### view eems #####

eem_overview_plot(eem_list, spp=9, contour = TRUE)


### end 


##### absorbence #####
absorbance_path = system.file("extdata/absorbance", package = "staRdom") 
# load example data, set a path without using system.file to use your own data e.g. "C:/folder/another folder"


## load the data ##

absorbance <- absorbance_read(absorbance_path, cores = cores) # load csv or txt tables in folder

### end absorb ### 



##### dilutions #####

metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom") # path to example data, can be replaced by a path to your own data
meta <- read.table(metatable, header = TRUE, sep = ",", dec = ".", row.names = 1) # load data


###


##### data and error checking #####

problem <- eem_checkdata(eem_list,absorbance,meta,metacolumns = c("dilution"),error=FALSE)


##### names fix #####

eem_list <- eem_name_replace(eem_list,c("\\(FD3\\)"),c(""))

###


##### Absorbance baseline correction #####
#The instrumental baseline drift in absorbance data can be corrected by subtracting the mean of the absorbance at high wavelengths (Li and Hur 2017). 
# The default is to use the spectrum between 680 and 700 nm but any other range can be set manually.
absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

###

##### Spectral correction #####
#done to remove instrument-specific influences on the EEMs

excorfile <- system.file("extdata/CorrectionFiles/xc06se06n.csv",package="staRdom")
Excor <- data.table::fread(excorfile)

emcorfile <- system.file("extdata/CorrectionFiles/mcorrs_4nm.csv",package="staRdom")
Emcor <- data.table::fread(emcorfile)

# adjust range of EEMs to cover correction vectors
eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1]))

eem_list <- eem_spectral_cor(eem_list,Excor,Emcor)

###

##### Blank subtraction #####
#Blanks are samples of ultrapure water that must contain either “nano”, “miliq”, “milliq”, “mq” or “blank” in their file names

# extending and interpolation data
eem_list <- eem_extend2largest(eem_list, interpolation = 1, extend = FALSE, cores = cores)

# blank subtraction
eem_list <- eem_remove_blank(eem_list)

###

### display ###
eem_overview_plot(eem_list, spp=9, contour = TRUE)



##### Inner-filter effect correction #####

eem_list <- eem_ife_correction(eem_list,absorbance, cuvl = 5)

###

##### Raman normalisation #####

eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")

# I think this is called normalization otherwise? ##

### display ###

eem_overview_plot(eem_list, spp=9, contour = TRUE)


##### Remove blanks from sample set #####

eem_list <- eem_extract(eem_list, c("nano", "miliq", "milliq", "mq", "blank"),ignore_case = TRUE)


###

absorbance <- dplyr::select(absorbance, -matches("nano|miliq|milliq|mq|blank", ignore.case = TRUE))

###

##### Remove and interpolate scattering #####

remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
remove_scatter_width <- c(15,15,15,15)

eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)

###

eem_overview_plot(eem_list, spp=9, contour = TRUE)

###


##### interp the blanks removed #####
# The types of interpolation are : (0) setting all NAs to 0, 
# (1) spline interpolation with mba.points (Lee, Wolberg, and Shin 1997), 
# (2) excitation and emission wavelength-wise interpolation with pchip (Moler 2004) and subsequent mean, 
# (3) excitation wavelength-wise interpolation with pchip and 
# (4) linear excitation and emission wavelength-wise interpolation with na.approx and again subsequent mean calculation.

# reccomended to start with type 1. #

eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)

#

eem_overview_plot(eem_list, spp=9, contour = TRUE)

#

##### Correct for dilution #####
#Each EEM is multiplied with a dilution factor (e.g. 10 if 1 part sample was diluted with 9 parts ultrapure water)
# mine are 50 sample to 50 ultra pure water (intended....)

dil_data <- meta["dilution"]

eem_list <- eem_dilution(eem_list,dil_data)

##

eem_overview_plot(eem_list, spp=9) # plot spared, due to no dilution it looks like the previous plot.

###


##### Smooth data #####
# Depending on your instrument, smoothing the data could be beneficial for peak picking. For PARAFAC analysis smoothing is not advised. 
# The parameter n specifies the moving average window size in nm.

eem4peaks <- eem_smooth(eem_list, n = 4, cores = cores)

###


##### Overview of samples #####

summary(eem_list)

###

##### Peak picking and indices #####
#Peaks and indices known from the literature (Huguet et al. 2009; Zepp, Sheldon, and Moran 2004; McKnight et al. 2001; Ohno 2002) can be calculated.
#
# If wavelength ranges needed for certain indices or peaks (e.g. the humification index (HIX) uses an excitation wavelength of 254 nm,
# but EEMs usually contain measurements at 250 and 255 nm) an interpolation is done automatically between existing measurements (no extrapolation!).


bix <- eem_biological_index(eem4peaks)
coble_peaks <- eem_coble_peaks(eem4peaks)
fi <- eem_fluorescence_index(eem4peaks)
hix <- eem_humification_index(eem4peaks, scale = TRUE)

indices_peaks <- bix %>%
  full_join(coble_peaks, by = "sample") %>%
  full_join(fi, by = "sample") %>%
  full_join(hix, by = "sample")

indices_peaks

###

##### Absorbance indices #####

# abs_parms can be used to calculate a254, a300, E2:E3, E4:E6, S275-295, S350-400, S300-700, SR and the wavelength distribution of absorption spectral slopes 
# (Helms et al. 2008; Twardowski et al. 2004; Massicotte 2016; Loiselle et al. 2009). 
# Li and Hur (2017) give a broad overview of possible values and their applications used in recent literature. 
# Missing wavelengths, needed for certain indices or ratios are interpolated automatically.


slope_parms <- abs_parms(absorbance, cuvl = 1, cores = cores)

slope_parms
## (above to print thesample to the Absorbance indices )##
###



######################################################################################################


##### Should this be it's own code page? #####



##### Creating a PARAFAC model #####

# Load drEEM example dataset #

dreem_raw <- tempfile()
download.file("http://models.life.ku.dk/sites/default/files/drEEM_dataset.zip",dreem_raw)
dreem_data <- unz(dreem_raw, filename="Backup/PortSurveyData_corrected.mat", open = "rb") %>%
  R.matlab::readMat()
unlink(dreem_raw)

eem_list <- lapply(dreem_data$filelist.eem, function(file){
  #file <- dreem_data$filelist.eem[1]
  n <- which(dreem_data$filelist.eem == file)
  file <- file %>%
    gsub("^\\s+|\\s+$", "", .) %>% # trim white spaces in filenames
    sub(pattern = "(.*)\\..*$", replacement = "\\1", .) # remove file extension from sample name
  eem <- list(file = paste0("drEEM/dataset/",file),sample = file,x = dreem_data$XcRU[n,,] %>% as.matrix(),ex = dreem_data$Ex %>% as.vector(), em = dreem_data$Em.in %>% as.vector(), location = "drEEM/dataset/")
  class(eem) <- "eem"
  attr(eem, "is_blank_corrected") <- TRUE
  attr(eem, "is_scatter_corrected") <- FALSE
  attr(eem, "is_ife_corrected") <- TRUE
  attr(eem, "is_raman_normalized") <- TRUE
  attr(eem, "manufacturer") <- "unknown"
  eem
}) %>%
  `class<-`("eemlist")

# add sample name suffix, R has sometimes troubles, when sample names start with a number.
eem_names(eem_list) <- paste0("d",eem_names(eem_list))

###

# In the drEEM tutorial, all samples containing “bl” or “0A” are removed from the set.#

ol <- function(x){x==("bl") | x == "0A"}
extract <- dreem_data$sites %>% unlist() %>% ol() %>% which()
eem_list <- eem_list %>% eem_extract(extract)

###


# Scattering has still to be removed and is done as described below. #

eem_list <- eem_rem_scat(eem_list, remove_scatter = c(TRUE, TRUE, TRUE, TRUE), remove_scatter_width = c(15,15,18,19), interpolation = FALSE, cores = cores)


### below is for special purposes #
# If you have worked with the basic analysis template, you can use the resulting data right away. 
# In the case you did the corrections with different sample sets separately and want to combine them, you can use eem_import_dir to combine EEM samples from several RData or RDa files.
# Put all files in one folder and run the following:

#  eem_list <- eem_import_dir(dir)

###



##### Sample set wavelength ranges #####

# eem_extend2largest adds NAs if values present in another sample are missing (e.g. different wavelength slits or spectra ranges were used). These values can be interpolated.

# eem_red2smallest removes wavelengths, that are missing in at least one sample form the whole set.

###

##### Find and remove noise in EEM data #####

##  eem_extract removes whole samples either by name or by number.

##  eem_range removes data outside the given wavelength ranges in all samples.

##  eem_exclude removes data from the sample set, provided by a list.

##   eem_rem_scat and eem_remove_scattering are used to set data in Raman and Rayleigh scattering of 1st and 2nd order to NA. While the later on removes one scattering at a time, the first one wraps it up to remove several scatterings in one step.

##   eem_setNA replaces data by NA in rectangular shape and in specific samples.

##   eem_matmult multiplies each EEM matrix by a certain matrix. This matrix can be used to set parts of the data to 0 or NA (e.g. the area where emission wavelength is shorter than excitation wavelength).

##   eem_interp is used to interpolate data previously set to NA.

### demo ###

# Sample “d667sf” find. match the language --> " ^d667sf$ ” is used. 
# ^ stands for the beginning of the string and $ for the end and ensure an exact match.


eem_list %>% 
  eem_extract(sample = "^d667sf$", keep = TRUE) %>%
  ggeem(contour = TRUE)


###


# The noisy range below 250 nm excitation and above 580 nm emission can be removed from the samples with the following command.
# As mentioned above, this was already removed in the example data.


eem_list <- eem_list %>% eem_range(ex = c(250,Inf), em = c(0,580))


###

### Demo - visual irregularities ###

# Visually found irregularities in patterns are manually replaced by NA and interpolated. 
# From the sample “d667sf” a band covering the excitation wavelengths from 245 to 350 is removed and a rectangle covering emission wavelengths 560 to 576 and excitation wavelengths 280 to 295 is removed in all samples. 

#

eem_list <- eem_list %>%
  eem_setNA(sample = 176, ex = 345:350, interpolate = FALSE) %>%
  eem_setNA(em = 560:576, ex = 280:295, interpolate = FALSE)

###


eem_list <- eem_interp(eem_list, type = 1, extend = FALSE, cores = cores)



#####





