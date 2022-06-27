#C_bacchus - endemic species south africa 
setwd("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica")

require(raster)
require(rgdal)

#Variables processing 
setwd('Layers')

## names of variables
path <- "C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\Layers\\current\\bio_2-5m_bil\\"
clim_list <- list.files(path= path,  pattern=".bil", full.names=TRUE)

# stacking the bioclim variables to process them at one go 
variables <- raster::stack(clim_list) 
require(rgdal)
shape <- readOGR(dsn = "C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica",layer = "projection") #SHAPEFILE YOU CREATED FOR PROJECTION AREA
var_mask <- mask(crop(variables, shape), shape)

rnames <- paste0("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\current\\", names(variables), ".asc") # users select the format

sav <- lapply(1:nlayers(var_mask), function(x) {
  writeRaster(var_mask[[x]], filename = rnames[x], format = "ascii", overwrite=TRUE) # change format accordingly
})
 
# replace names

filez <- list.files('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\current')
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\current')
file.rename(from=filez, to=sub(pattern="bio", replacement="bio_", filez))

     
#Selecting species
#Number of records
Occ <- read.csv( 'C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\C_bacchus\\C_bacchus.csv')

#Thinning Dataset------


#Download worldclim variables version 2.1--------------
#
# argument options
# period: "historical", "future"
# res: 10m, 5m, 2.5m, 30s (30s only for historical)
# time: 2021-2040, 2041-2060, 2061-2080, 2081-2100 (only for future) 
# SSP: 126, 245, 370, 585 (only for future) 
# GCM: BCC-CSM2-MR, CNRM-CM6-1, CNRM-ESM2-1, CanESM5, IPSL-CM6A-LR, MIROC-ES2L, 
#      MIROC6, MRI-ESM2-0 (only for future) 
# output_dir: the directory the user wants to put the results in 


# current
require(raster)
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\current')

#source("Functions (1).R")

#Download climate data 
#clim <- get_NWC_bio(period = "historical", res = "10m", time = "2041-2060", SSP = "NULL",GCM = NULL) 

#create the full path for the raster you will split the raster stack you just created
#rnames <-  paste0("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Carmen\\Present20\\", names(clim), ".asc")

#sav <- lapply(1:nlayers(clim), function(x) {
#  writeRaster(clim[[x]], rnames[[x]], format = "ascii", overwrite=T) # change format accordingly
#})

#changing names ----
#setwd("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Carmen\\Present20")
#Files <- list.files(all.files = TRUE, pattern = ".asc")
#newName <- sub("wc2.1_10m_bio_", "bio_", Files)

#file.rename(Files, newName)


#changing names ----
#setwd("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Carmen\\Present20")
#Files <- list.files(all.files = TRUE, pattern = ".asc")
#newName <- sub("wc2.1_10m_bio_", "bio_", Files)

#file.rename(Files, newName)

#get_maxent(version = "latest", quiet = FALSE)

#Crop calibration area 
#M- calibration area---
install.packages('ellipsenm')
require(ellipsenm)

occ <- read.csv('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\C_bacchus.csv')

#Buffer distance is in Kilometers 
#For mandarina we will assum 2 buffer areas 
sp <- buffer_area(na.omit(occ), longitude = "Longitude", latitude = "Latitude", 
                    buffer_distance = 50)

path <- "C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\current"
Calib <- stack(list.files(path = path, pattern = ".asc", full.names = TRUE))
var_mask <- mask(crop(Calib, sp), sp)

## names for layers
rnames <- paste0("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\M_variables\\Set1\\", names(Calib), ".asc") 

## saving layers in new folder
sav <- lapply(1:nlayers(var_mask), function(x) {
  writeRaster(var_mask[[x]], filename = rnames[x], format = "ascii", overwrite=T) # change format accordingly
})
filez <- list.files('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\BC8570\\')
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\G_variables\\Set1\\BC8570')
file.rename(from=filez, to=sub(pattern="bc45bi70", replacement="bio_", filez))

#Start models------
require(kuenm)
require(raster)
require(devtools)

setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species')
spvec <- dir('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species')
exclude <- c( "G_variables" )
spvector <- spvec[!spvec %in% exclude]
i=1


all <- read.csv(paste(spvector[i], "/", paste(spvector[i], ".csv", sep = ""), sep = ""))
all <- unique(all)


all$check <- paste(all[,2], all[,3], sep = "_")
train <- all[sample(nrow(all), round((length(all[,1])/4 *3))), ]
test <- all[!all[,4] %in% train[,4], ]

all$check <- NULL
train$check <- NULL
test$check <- NULL

write.csv(all, paste(spvector[i], "/", paste(spvector[i], "_joint.csv", sep = ""), sep = ""),
          row.names = FALSE)
write.csv(train, paste(spvector[i], "/", paste(spvector[i], "_train.csv", sep = ""), sep = ""),
          row.names = FALSE)
write.csv(test, paste(spvector[i], "/", paste(spvector[i], "_test.csv", sep = ""), sep = ""),
          row.names = FALSE)

sp1 <- spvector

# occurrences
occ_joint <- paste(sp1, "/", paste(sp1, "_joint.csv", sep = ""), sep = "")
occ_tra <- paste(sp1, "/", paste(sp1, "_train.csv", sep = ""), sep = "")
occ_test <- paste(sp1, "/", paste(sp1, "_test.csv", sep = ""), sep = "")

# variables
M_var_dir <- paste(sp1, "/", "M_variables", sep = "")# m for each species
G_var_dir <- "G_variables" #g for all species

# other arguments
batch_cal <- paste(sp1, "/", "Candidate_models", sep = "")
out_dir <- paste(sp1, "/", "Candidate_Models", sep = "")
reg_mult <- c( 1, 1.5, 2, 3, 4, 5)
f_clas <- "no.t.h"
background <- 10000
maxent_path <- "C:\\Maxent"
wait <- TRUE
run <- TRUE
out_eval <- paste(sp1, "/", "Calibration_results", sep = "")
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- FALSE
selection <- "OR_AICc"
paral_proc <- FALSE
batch_fin <- paste(sp1, "/", "Final_models", sep = "")
mod_dir <- paste(sp1, "/", "Final_Models", sep = "")
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "logistic"
project <- TRUE
ext_type <- "all"
write_mess <- TRUE
write_clamp <- FALSE
args <- NULL
wait <-FALSE
run1 <- TRUE
i=1
# runing calibration and final models in loop

# create candidate models
kuenm_cal(occ.joint = occ_joint[i], occ.tra = occ_tra[i], M.var.dir = M_var_dir[i], 
          batch = batch_cal[i], out.dir = out_dir[i], reg.mult = reg_mult, f.clas = f_clas, 
          maxent.path = maxent_path, wait = wait, run = run)
# evaluate all candidate models and select the best based on pROC, Omission rate and complexity AICc
kuenm_ceval(path = out_dir[i], occ.joint = occ_joint[i], occ.tra = occ_tra[i], occ.test = occ_test[i],
            batch = batch_cal[i], out.eval = out_eval[i], threshold = threshold,
            rand.percent = rand_percent, iterations = iterations, kept = kept,
            selection = selection, parallel.proc = paral_proc)

#Calibration results thinning ---------
j=1
best <- list.files(path = paste(spvector[j], "Calibration_results", sep = "/"),
                   pattern = "best", full.names = TRUE)
bestr <- read.csv(best)

file.rename(best, paste(spvector[j], "Calibration_results", "tseb_candidate_models_OR_AICc.csv", sep = "/"))

bestr <- bestr[bestr[, 6] == 0, ]

if (dim(bestr)[2] > 1) {
  models <- bestr[, 1]
  
  sn <- ".*set2"
  stn <- gregexpr(sn, models)
  stan <- regmatches(models, stn)
  statn <- unlist(stan)
  
  if (length(statn) > 0) {
    bestb <- bestr[bestr[, 1] == statn, ][1, ]
  }else {
    bestb <- bestr[1, ]
  }
}else {
  bestb <- bestr
}

write.csv(bestb, file = best, row.names = FALSE)

i=1
# create final models using the parameterizations selected before

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = TRUE, write.clamp = FALSE,
          maxent.path = maxent_path, args = NULL, wait = TRUE, run =  TRUE) 

#model threshold------------
sp1 <- spvector 
sp2 <- gsub(" ", "_", sp1)

sp_name <- "bacchus"
fmod_dir <- paste(sp1, "Final_Models", sep = '/')
format <- "asc"
project <- TRUE
stats <- c("med", "range")
rep <- TRUE

#scen <- dir('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Maria_eduarda\\KUENM4\\G_variables\\Set1')
scenarios <- c("BC4550",  "BC4570",  "BC8550",  "BC8570" , "current", 
               "CC4550",  "CC4570",  "CC8550",  "CC8570" , 
               "MC4550",  "MC4570",  "MC8550",  "MC8570" ,
               "MR4550",  "MR4570",  "MR8550",  "MR8570" ,
               "LGM-MIROC",  "LGM-CCSM4")    
scenarios2 <- c("current", "LGM-MIROC",  "LGM-CCSM4")    
ext_type <- c("EC") # you can select only one type of extrapolation if needed
out_dir <- paste(sp1, "Final_Model_Stats2", sep = '/')

# other arguments were defined before
occ <- paste(sp1, paste(sp1, "joint.csv", sep = '_'), sep = "/")
thres <- 10
curr <- "current"
emi_scenarios<- c( '45', '85')
c_mods <- c('BC', 'MR', 'MC', 'CC','MIROC', 'CCSM4')
c_mods2 <- c('MIROC', 'CCSM4')
time_periods <- c('50', '70','LGM')
time_periods2 <- c('LGM')
out_dir1 <-  paste(sp1, "Projection_Changes1", sep = '/')

split <- 200 #partitioning the stack of models of different scenarios
out_dir2 <- paste(sp1, "Variation_from_sourcesPres", sep = '/')
out_dir3 <- paste(sp1, "Variation_from_sourcesLGM", sep = '/')
i=1

#sp_name <- always make sure to check  the name in the rasters of the models. Sometimes
#the name might be different and there will appear the error 
#Error in .local(x, ...) : no filenames supplied

#Model statistics -----
kuenm_modstats(sp.name = sp_name[i], fmod.dir = fmod_dir[i], format = format, project = project, 
               statistics = stats, replicated = rep, proj.scenarios = scenarios2, 
               ext.type = ext_type, out.dir = out_dir[i]) 

#Threshold models --------------
kuenm_projchanges(occ = occ[i], fmod.stats = out_dir[i], threshold = 5, 
                  current = curr, time.periods = time_periods, 
                  clim.models = c_mods, ext.type = ext_type, out.dir = out_dir1[i]) 

### Variance in model predictions---------
kuenm_modvar(sp.name = sp_name[i], fmod.dir = fmod_dir[i], is.swd = F,
             replicated = rep, format = format, project = project,
             current = curr, emi.scenarios = emi_scenarios,time.periods=	time_periods,
             clim.models = c_mods, ext.type = ext_type, split.length = split,
             out.dir = out_dir2[i])

kuenm_modvar(sp.name = sp_name[i], fmod.dir = fmod_dir[i], is.swd = F,
             replicated = rep, format = format, project = project,
             current = curr,  clim.models = c_mods2, ext.type = ext_type,
             split.length = split, out.dir = out_dir3[i])
#Mop----

percent <- 3 
paral <- FALSE 
comp_each <- 3000 
j=1

sets <- list.files(path = paste(spvector[j], "Calibration_results", sep = "/"),
                   pattern = "best", full.names = TRUE)
sets_var <- as.character(read.csv(sets)[1, 1])

sets_var <- gsub("^M.*_", "", sets_var)
M_var_dir <- paste(spvector[j], "M_variables", sep = "/")
out_mop <- paste(spvector[j], "MOP_results", sep = "/")

kuenm_mmop(G.var.dir = G_var_dir, M.var.dir = M_var_dir,
           sets.var = sets_var, out.mop = out_mop, is.swd = F, 
           percent = percent, parallel = paral, 
           comp.each = comp_each)

#Binary and Binary MOP - Bin
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\MOP_results\\Set1')
mops <- "MOP_\\d\\S*.tif$"
bin <-
files <- list.files(pattern = mops, full.names = T, recursive = T) #list MOPS

for (i in 1:length(files)) {
        mop <- raster(files[i])   
        mop <- mop > 0 # strict extrapolation
        #dir.create(MAIN[i]) #creating folder to save binary version of MOPS - MOP_results_bin
        #dir.create(SET[i])
        writeRaster(mop, filename = paste (files[i],"_bin"), format = "GTiff", overwrite=TRUE)
        cat(i, "of", length(files), "\n")
}
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\Projection_Changes')        

bins <- "binary.*if$"
bin_files <- list.files(pattern = bins, full.names = T, recursive = T)
bin_files

bins_exc <- "binary_c.*on.tif$"
bin_files_exc <- list.files(pattern = bins_exc, full.names = T, recursive = T)
bin_files <- bin_files[!bin_files %in% bin_files_exc]


#Manually transform in Binary minus MOP
setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\Projection_Changes\\Changes_EC\\Period_50\\Scenario_45\\Binary')
b <-raster('.\\binary_CC.tif')
m <-raster('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\MOP_results\\Set1\\MOP_3%_BC4550.tif _bin.tif')
bin1 <- b * m
writeRaster(bin1, '.\\binary_ME_BC.tif', format = "GTiff", overwrite=TRUE)


setwd("C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species2")

require(raster)
require(rgdal)

#cropping Study area
require(rgdal)
shape <- readOGR(dsn = "C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica",layer = "SALed_shape") #SHAPEFILE YOU CREATED FOR PROJECTION AREA
plot(shape)


setwd('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/Final_Model_StatsSALed/Statistics_EC')
dir <- list.files(full.names = T,  all.files = T, recursive = T)

setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\MOP-e SALed' )
dir2<- list.files(full.names = T,  all.files = T, recursive = T)
dir2

setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\MOP_resultsSALed\\Set1' )
dir3<- list.files(full.names = T,  all.files = T, recursive = T)
dir3

setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\Projection_ChangesSALed\\Changes_EC' )
dir2<- list.files(full.names = T,  all.files = T, recursive = T)
dir2

setwd('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Projects\\Beetle_SAfrica\\species\\C_bacchus\\Variation_from_sourcesSALed2\\Variation_EC' )
dir4<- list.files(full.names = T,  all.files = T, recursive = T)
dir4

i=1
for (i in 1:length(dir)) {
  x <- raster(dir[i])
  var_mask<- mask(crop(x, shape), shape)
  writeRaster(var_mask, filename = paste(dir[i]), overwrite=TRUE)
  }
}

#Multiply binary * MOP

##Binary
setwd('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/Projection_ChangesSA')

bins <- "binary.*if$"
bin_files <- list.files(pattern = bins, full.names = T, recursive = T)
bin_files
bins_exc <- "binary_c.*on.tif$"
bin_files_exc <- list.files(pattern = bins_exc, full.names = T, recursive = T)
bin_files <- bin_files[!bin_files %in% bin_files_exc]
bins_exc2 <- "binary_current.tif"
bin_files_exc2 <- list.files(pattern = bins_exc2, full.names = T, recursive = T)
bin_files2 <- bin_files[!bins_exc2 %in% bin_files_exc]
ex <- c(3, 8, 13, 18, 22 ) #Number of bin_files 2 bcs of binary_ and binary_current
bin_files3 <- bin_files2[-ex]
bin_files3 

setwd('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/MOP_resultsSALed/Set1')

mopr <- "_bin"
mop <- list.files(pattern = mopr, full.names = T, recursive = T)
ex <- c(9 ) #Number of bin_files 2 bcs of binary_ and binary_current
mop_ <- mop[-ex]
mop_
MOP-e
path <- 'C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/MOP-e' 

i=5
  setwd('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/MOP_resultsSALed/Set1')
  mopr <- "_bin"
  mop <- list.files(pattern = mopr, full.names = T, recursive = T)
  mop <- raster(mop[9])
  mop
  plot(mop)
  setwd('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/Projection_ChangesSALed')
  bin_files3
  bin <- raster(bin_files3[18])
  bin 
  plot(bin)
  bin1 <- bin * mop
  writeRaster(bin1, filename = paste(path, 'LGM_MIROC', sep = '/'), format = "GTiff", overwrite=TRUE)
  
  
 #Current -----------
  mop <- raster('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/MOP_resultsSALed/Set1/MOP_3%_current.tif _bin.tif ')
  setwd('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/Projection_ChangesSALed')
  bin_files3
  bin <- raster('C:/Users/mari1/Dropbox (Senckenberg)/Projects/Beetle_SAfrica/species/C_bacchus/Projection_ChangesSALed/Changes_EC/Period_50/Scenario_45/binary_comparison.tif')
  bin1 <- bin * mop
  writeRaster(bin1, filename = paste(path, 'current', sep = '/'), format = "GTiff", overwrite=TRUE)
 #----------
  
  
#in case we need to calculate the percentage
  # before
  vals <- na.omit(values(bin))
  uval <- sort(unique(vals))
  
  ## processing
  all <- length(vals)
  percents <- sapply(1:length(uval), function(x) {
    per <- round(sum(vals == uval[x]) * 100 / all, 2)
    return(per)
  })
  names(percents) <- uval
  
  # after
  vals <- na.omit(values(bin1))
  uval <- sort(unique(vals))
  
  ## processing
  all <- length(vals)
  percents1 <- sapply(1:length(uval), function(x) {
    per <- round(sum(vals == uval[x]) * 100 / all, 2)
    return(per)
  })
  names(percents1) <- uval
  
  results[[i]] <- list(before = percents, after = percents1)
  

  


#Calculation of Area: 
r<- raster('C:\\Users\\mari1\\Dropbox (Senckenberg)\\Maria_eduarda\\KUENM4\\waltheri\\Projection_Changes\\Changes_EC\\Period_50\\Scenario_45\\Binary\\binary_MRI.tif')

#get sizes of all cells in raster [km2]
cell_size<-area(r, na.rm=TRUE, weights=FALSE)

#delete NAs from vector of all raster cells
##NAs lie outside of the rastered region, can thus be omitted
cell_size<-cell_size[!is.na(cell_size)]

#compute area [km2] of all cells in geo_raster
raster_area<-length(cell_size)*median(cell_size)

#print area of Georgia according to raster object
print(paste("Area of Georgia (raster):",round(raster_area, digits=1),"km2"))


require(raster)
#Don't forget to indicate the species, GCMS, , RCP and time periods in the names you are using to save the maps you created
bin <- raster ('UPLOAD BINARY SUITABILITY  MAP') #These maps are in the folder Projection Changes
plot(bin) #try to always plot the maps you create, to avoid problems and mixed-ups
mop <- raster('UPLOAD MOP MAP')
Binmop <- mop>0 # creating Binary mop - transforming everything bigger than 0 to 1 - strict extrapolation
plot(Binmop)
writeRaster(Binmop, 'SAVE IN THE EXTRA FOLDER WE DISCUSSED')

NonExtrSuit <- bin * Binmop # This is the non extrapolative suitability final map
plot(NonExtrSuit)
writeRaster(NonExtrSuit, 'SAVE IN THE EXTRA FOLDER WITH NON EXTRAOLATIVE SUITABILITY MAPS')