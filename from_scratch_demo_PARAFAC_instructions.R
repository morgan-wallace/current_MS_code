##### Creating a PARAFAC model #####

##### installing the pre-rec libraries #####

library(dplyr)
library(tidyr)

###


##### Run just the libraries #####

library("dplyr")
library("tidyr")
library("staRdom")

##### end libraries #####

##### start instructions? #####

##### Loading data #####

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





