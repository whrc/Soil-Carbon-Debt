## Derivation of potential soil carbon (https://github.com/whrc/Soil-Carbon-Debt)
## Code by: Tom.Hengl@isric.org
## Contributions by: J. (Jon) Sanderman (WHRC) and G. (Greg) Fiske (WHRC)
## Cite as: Sanderman, J., Hengl, T., Fiske, G., 2017? "The soil carbon debt of 12,000 years of human land use", sumbitted to PNAS.

list.of.packages <- c("raster", "rgdal", "nnet", "plyr", "ROCR", "randomForest", "R.utils", "plyr", "parallel", "psych", "mda", "dismo", "snowfall", "hexbin", "lattice", "ranger", "xgboost", "doParallel", "caret", "plotKML", "GSIF")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Packages in use:
setwd("/data/WHRC_soilcarbon/model10km")
load(".RData")
library(plyr)
library(aqp)
library(stringr)
library(sp)
library(rgdal)
library(devtools)
#devtools::install_github('dmlc/xgboost')
library(xgboost)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger")
library(ranger)
library(nnet)
library(caret)
library(hexbin)
library(snowfall)
library(utils)
library(plotKML)
library(GSIF)
library(raster)
library(R.utils)
library(doParallel)
library(foreign)
library(tools)
#library(doSNOW)
#library(doMC)
#library(randomForestSRC)
library(parallel)
#library(mxnet)
#options(rf.cores=detectCores(), mc.cores=detectCores())

load("/data/models/equi7t3.rda")
plotKML.env(convert="convert", show.env=FALSE)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate =  "/usr/bin/gdal_translate"
  gdalwarp =  "/usr/bin/gdalwarp"
  gdalbuildvrt = "/usr/bin/gdalbuildvrt"
  saga_cmd = "/usr/local/bin/saga_cmd"
}
system("gdal-config --version")
source("/data/models/saveRDS_functions.R")
source("/data/models/extract.equi7t3.R")
source("WHRC_functions.R")

## all processing done on ca 10 km grid
ncols = 4320
nrows = 2160
xllcorner = -180
yllcorner = -90
xurcorner = 180
yurcorner = 90
cellsize = 0.0833333
NODATA_value = -9999

library(maps)
library(maptools)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")

## 0. Input data -----------

## Climatic variables (https://crudata.uea.ac.uk/cru/data/hrg/)
## Harris, I., Jones, P.D., Osborn, T.J. and Lister, D.H. (2014), Updated high-resolution grids of monthly climatic observations – the CRU TS3.10 Dataset. Int. J. Climatol., 34: 623–642. doi: 10.1002/joc.3711
## Downloaded from: http://www.ipcc-data.org/observ/clim/cru_ts2_1.html
cl.zip.lst <- list.files(path="./climate", pattern = glob2rx("*.zip$"), full.names = TRUE)
## extract:
sapply(cl.zip.lst, function(x){system(paste("7za e ", x," -r -y"))})
cl.lst <- list.files(pattern=glob2rx("cru_*_*_*-*_*.tif")) ## "cru_*_*_1901-1930_*.tif"
## 324 layers at 50 km resolution
GDALinfo("cru_pre_clim_1961-1990_01.tif")
#r = readGDAL("cru_pre_clim_1961-1990_01.tif")
#plot(raster(r))
GDALinfo("cru_frs_clim_1961-1990_01.tif")
GDALinfo("cru_tmn_clim_1961-1990_01.tif")
GDALinfo("cru_tmx_clim_1961-1990_01.tif")
GDALinfo("cru_tmp_clim_1961-1990_01.tif")
GDALinfo("cru_vap_clim_1961-1990_01.tif")
GDALinfo("cru_wet_clim_1961-1990_01.tif")
GDALinfo("cru_cld_clim_1961-1990_01.tif")
GDALinfo("cru_dtr_clim_1961-1990_01.tif")
## Missing value flags are incorrect?
gdalwarp_clim <- function(x, srcnodata = "254", dstnodata = "255"){
  out = paste0("./stacked/", gsub(".tif", "_10km.tif", basename(x)))
  if(!file.exists(out)){ system(paste0(gdalwarp, ' ', x, ' ', out, ' -srcnodata \"', srcnodata,'\" -dstnodata \"', dstnodata,'\" -co \"COMPRESS=DEFLATE\" -r \"cubicspline\" -tr ', cellsize, ' ', cellsize))} ##  -t_srs \"+proj=longlat +datum=WGS84\" ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner
}
#unlink(paste0("./stacked/", gsub(".tif", "_10km.tif", basename(cl.lst[1]))))
#gdalwarp_clim(cl.lst[1])
## resample to 10 km resolution:
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalwarp", "gdalwarp_clim", "cl.lst", "cellsize", "xllcorner", "yllcorner", "xurcorner", "yurcorner")
out <- sfClusterApplyLB(cl.lst, gdalwarp_clim)
sfStop()
unlink(cl.lst)

## DEM parameters / Geology and landform
## from 250m to 10 km:
des <- read.csv("/data/models/SoilGrids250m_COVS250m.csv")
tcovs <- c(as.character(des$WORLDGRIDS_CODE[c(grep("L??USG5", des$WORLDGRIDS_CODE), grep("???MRG5", des$WORLDGRIDS_CODE), grep("L??USG5", des$WORLDGRIDS_CODE), grep("F??USG5", des$WORLDGRIDS_CODE), grep("C??GLC5", des$WORLDGRIDS_CODE))]), "MNGUSG")
#source("mosaick_function.R")
## 84 layers
#sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>25,25,length(tcovs)))
#sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdal_translate", "ext", "tcovs", "mosaick.equi7t3", "make_mosaick", "tile.names")
#out <- sfLapply(tcovs, function(x){try( make_mosaick(i="dominant", varn=x, ext=ext, in.path="/data/covs1t", tr=0.00833333, ot="Int16", dstnodata=-32768, tile.names=tile.names) )})
#sfStop()
#geog.lst <- c(list.files(path="./GEOG", pattern=glob2rx("F*_agg_ll.tif$"), full.names =TRUE), list.files(path="./GEOG", pattern=glob2rx("L*_agg_ll.tif$"), full.names =TRUE), list.files(path="./GEOG", pattern=glob2rx("*MRG5_agg_ll.tif$"), full.names =TRUE), "./GEOG/MNGUSG_agg_ll.tif")
geog.lst <- c(list.files(path="/data/stacked250m", pattern="USG5", full.names =TRUE), list.files(path="/data/stacked250m", pattern="MRG5", full.names =TRUE), paste0("/data/stacked250m/S0",3:9,"ESA4.tif"), paste0("/data/stacked250m/S10ESA4.tif"), "./GEOG/MNGUSG_agg_ll.tif", paste0("/data/EarthEnv/MODCF_monthlymean_0",c(1:9),".tif"), paste0("/data/EarthEnv/MODCF_monthlymean_",c(10:12),".tif"), "/data/DAAC/average_soil_and_sedimentary-deposit_thickness.tif")
## 50 layers
## resample to 10 km resolution:
sfInit(parallel=TRUE, type="SOCK", cpus=50)
sfExport("gdalwarp", "geog.lst", "cellsize", "xllcorner", "yllcorner", "xurcorner", "yurcorner")
out <- sfClusterApplyLB(geog.lst, function(x){ if(!file.exists(paste0('./stacked/', gsub(".tif", "_10km.tif", basename(x))))){ system(paste0(gdalwarp, ' ', x, ' ./stacked/', gsub(".tif", "_10km.tif", basename(x)), ' -co \"COMPRESS=DEFLATE\" -r \"average\" -t_srs \"+proj=longlat +datum=WGS84\"  -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner)) }})
sfStop()
#unlink(geog.lst)

## Intact areas / protected planet:
## 2014, IFL Mapping Team: Greenpeace, University of Maryland, Transparent World, World Resource Institute, WWF Russia. Results/reports can be viewed at www.intactforests.org
## Protected planet WPDA Data (http://www.protectedplanet.net/)
pp.zip.lst <- list.files(path="./intact", pattern = glob2rx("*.zip$"), full.names = TRUE)
sapply(pp.zip.lst, function(x){system(paste("7za x ", x," -r -y"))})
ogrInfo("WDPA_Apr2016-shapefile-polygons.shp", "WDPA_Apr2016-shapefile-polygons")
ogrInfo("ifl_2013.shp", "ifl_2013")
## rasterize:
rasterize_pol <- function(INPUT, FIELD, cellsize, xllcorner, yllcorner, xurcorner, yurcorner){
  out = paste0(strsplit(basename(INPUT), "\\.")[[1]][1], ".sdat")
  ## "./stacked10km/", 
  if(!file.exists(out)){
    system(paste0('/usr/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"', INPUT, '\" -FIELD \"', FIELD, '\" -GRID \"', gsub(".sdat", ".sgrd", out), '\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
  }
}
shp.lst = c("WDPA_Apr2016-shapefile-polygons.shp", "ifl_2013.shp", "ifl_2000.shp")
field.lst = c("STATUS_YR","IFL_ID","IFL_ID")
x = sapply(1:length(shp.lst), function(x){rasterize_pol(INPUT=shp.lst[x], FIELD=field.lst[x], cellsize, xllcorner, yllcorner, xurcorner, yurcorner)})
#plot(stack(gsub(".shp", ".sdat", basename(shp.lst[-1]))))

## Terrestrial ecoregions (http://maps.tnc.org/gis_data.html)
download.file("http://maps.tnc.org/files/shp/terr-ecoregions-TNC.zip", "terr-ecoregions-TNC.zip")
system(paste("7za x terr-ecoregions-TNC.zip -r -y"))
#ogrInfo("terr-ecoregions-TNC.shp", "terr-ecoregions-TNC")
ecoregions.db <- read.dbf("tnc_terr_ecoregions.dbf")
str(ecoregions.db)
str(levels(ecoregions.db$WWF_MHTNAM))
ecoregions.db$WWF_i = as.integer(ecoregions.db$WWF_MHTNAM)
write.dbf(ecoregions.db, "tnc_terr_ecoregions.dbf") 
rasterize_pol(INPUT="tnc_terr_ecoregions.shp", FIELD="WWF_i", cellsize, xllcorner, yllcorner, xurcorner, yurcorner)
ecoregions_leg = data.frame(Value=1:length(levels(ecoregions.db$WWF_MHTNAM)), Classes=levels(ecoregions.db$WWF_MHTNAM))
write.csv(ecoregions_leg, "ecoregions_leg.csv")
unlink("./stacked/tnc_terr_ecoregions_10km.tif")
system(paste0(gdalwarp, ' tnc_terr_ecoregions.sdat ./stacked/tnc_terr_ecoregions_10km.tif -co \"COMPRESS=DEFLATE\" -r \"near\" -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner))

## GAUL country borders:
unlink("./stacked/GAUL_COUNTRIES_10km.tif")
system(paste0(gdalwarp, ' /data/aggregated/GAUL_COUNTRIES_1km.tif ./stacked/GAUL_COUNTRIES_10km.tif -co \"COMPRESS=DEFLATE\" -t_srs \"+proj=longlat +datum=WGS84\" -r \"near\" -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner))

## Hyde data set:
system('wget -r -l1 --no-parent ftp://ftp.pbl.nl/hyde/hyde3.2/2016_beta_release/zip/')
hyde.lst <- list.files(path="./ftp.pbl.nl/hyde/hyde3.2/2016_beta_release/zip", pattern = glob2rx("*_lu.zip$"), full.names = TRUE)
## extract files and stack
sapply(hyde.lst, function(x){system(paste("7za e ", x," -r -y"))})
cropland.lst <- list.files(pattern = glob2rx("cropland*.asc$"), full.names = TRUE) 
## 74 slices
cropland.tbl <- data.frame(filename=cropland.lst)
## function to strip year:
strip_year = function(x, name, ext=".asc"){
  xn = sapply(paste(x), function(x){strsplit(strsplit(x, name)[[1]][2], ext)[[1]][1]})
  xi = rep(NA, length(xn))
  xi[grep("AD", xn)] <- as.numeric( sapply(paste(xn[grep("AD", xn)]), function(x){strsplit(x, "AD")[[1]][1]}) )
  xi[grep("BC", xn)] <- -as.numeric( sapply(paste(xn[grep("BC", xn)]), function(x){strsplit(x, "BC")[[1]][1]}) )
  return(xi)
}

cropland.tbl$Year <- strip_year(cropland.tbl$filename, name="cropland")
cropland.tbl <- cropland.tbl[order(cropland.tbl$Year),]
pasture.lst <- list.files(pattern = glob2rx("pasture*.asc$"), full.names = TRUE)
pasture.tbl <- data.frame(filename=pasture.lst)
pasture.tbl$Year <- strip_year(pasture.tbl$filename, name="pasture")
pasture.tbl <- pasture.tbl[order(pasture.tbl$Year),]

#library(animation)
plot_world10km <- function(i, tbl){
  out.file = paste0(tbl[i,"filename"], ".png")
  if(!file.exists(out.file)){
    png(file = out.file, width = 4320/2, height = 2160/2, type="cairo")
    par(mar=c(0,0,0,0), oma=c(0,0,0,0))
    image(raster(paste(tbl[i,"filename"])), asp=1, col=c(SAGA_pal[["SG_COLORS_YELLOW_RED"]], rep("#BF0000",20)), zlim=c(0,100))
    text(-65, paste(tbl[i,"Year"]), cex=6)
    lines(country, col="black")
    dev.off()
  }
}
## Create animation:
sapply(1:nrow(cropland.tbl), plot_world10km, cropland.tbl)
system(paste0('convert -delay 100 ', paste(gsub(".asc", ".asc.png", cropland.tbl$filename), collapse=" "), ' cropland_historic_Hyde.gif'))
sapply(1:nrow(pasture.tbl), plot_world10km, pasture.tbl)
system(paste0('convert -delay 100 ', paste(gsub(".asc", ".asc.png", pasture.tbl$filename), collapse=" "), ' pasture_historic_Hyde.gif'))
## Resample to 10 km:
asc.lst = list.files(pattern = glob2rx("*.asc$"), full.names = TRUE)
## many layers!
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalwarp", "asc.lst", "cellsize", "xllcorner", "yllcorner", "xurcorner", "yurcorner")
out <- sfClusterApplyLB(asc.lst, function(x){ system(paste0(gdalwarp, ' ', x, ' ./stacked/', gsub(".asc", "_10km.tif", basename(x)), ' -co \"COMPRESS=DEFLATE\" -r \"near\" -t_srs \"+proj=longlat +datum=WGS84\"  -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner)) })
sfStop()
unlink(asc.lst)

## MODIS land cover for years 2000-2015:
mcd.lst <- c(paste0("/data/MCD12Q1/LandCover_", 2001:2013, "001_L1_500m.tif"), paste0("/data/MCD12Q1/LandCover_", 2001:2013, "001_L5_500m.tif"))
#sapply(mcd.lst, function(x){system(paste0(gdalwarp, ' ', x, ' ', basename(x), ' -co \"COMPRESS=DEFLATE\" -r \"near\" -s_srs \"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +units=m +no_defs\" -t_srs \"+proj=longlat +datum=WGS84\"  -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner))})
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalwarp", "mcd.lst", "cellsize", "xllcorner", "yllcorner", "xurcorner", "yurcorner")
out <- sfClusterApplyLB(mcd.lst, function(x){ out.file = paste0("./stacked/", gsub("500m", "10km", basename(x))); if(!file.exists(out.file)){ system(paste0(gdalwarp, ' ', x, ' ', out.file, ' -co \"COMPRESS=DEFLATE\" -r \"near\" -s_srs \"+proj=sinu +R=6371007.181 +nadgrids=@null +wktext\" -t_srs \"+proj=longlat +datum=WGS84\"  -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner))} })
sfStop()

## Historic forest cover
## http://www.unep-wcmc.org/resources-and-data/generalised-original-and-current-forest
shpF.lst = c("./other/ofc_gen.shp", "./other/cfc_gen.shp")
shpF.db.lst <- lapply(gsub(".shp", ".dbf", shpF.lst), read.dbf)
for(i in 1:length(shpF.db.lst)){
  shpF.db.lst[[i]]$TYPE_INT = as.integer(shpF.db.lst[[i]]$TYPE)
  write.dbf(shpF.db.lst[[i]], gsub(".shp", ".dbf", shpF.lst[[i]])) 
}
fieldF.lst = c("TYPE_INT","TYPE_INT")
x = sapply(1:length(shpF.lst), function(x){rasterize_pol(INPUT=shpF.lst[x], FIELD=fieldF.lst[x], cellsize, xllcorner, yllcorner, xurcorner, yurcorner)})
forestcover_leg = data.frame(Value=1:length(levels(shpF.db.lst[[1]]$TYPE)), Classes=levels(shpF.db.lst[[1]]$TYPE))
write.csv(forestcover_leg, "forestcover_leg.csv")

## convert to geotifs:
sapply(list.files(pattern=glob2rx("*.sdat$")), function(x){ system(paste0(gdalwarp, ' ', x, ' ./stacked/', gsub(".sdat", "_10km.tif", basename(x)), ' -co \"COMPRESS=DEFLATE\" -r \"near\" -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner)) })

## Potential wetlands GIEMS (http://www.estellus.fr/index.php?static13/giems-d15):
unlink("./stacked/giems_d15_v10_10km.tif")
system(paste0(gdalwarp, ' /data/EarthEnv/giems_d15_v10.tif ./stacked/giems_d15_v10_10km.tif -co \"COMPRESS=DEFLATE\" -t_srs \"+proj=longlat +datum=WGS84\" -r \"average\" -tr ', cellsize, ' ', cellsize, ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner))
plot(raster("./stacked/giems_d15_v10_10km.tif"), col=SAGA_pal[[1]])

## function to fill in missing values:
source("/data/models/tiler.R")
missing.tifs = list.files(path="./stacked", pattern=glob2rx("MODCF_monthlymean_*_10km.tif"), full.names=TRUE)
#, list.files(path="./stacked", pattern="ESA4", full.names=TRUE))
for(i in 1:length(missing.tifs)){
  close.gaps(inputTile=missing.tifs[i], maskTile="landmask_10km.sgrd", outTile=missing.tifs[i], ot="Int16", nodata="32767", nodata_out="-32768", a_srs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", zmin=0, method="stepwise", cpus=56, fix.zmin=FALSE)
}

## Load soil profiles:
#SOCS_profs = readRDS("SOCS_global_profs_December_2016.rds")
OCD_profs = readRDS("OCD_global_profs_March_2017.rds")

## Add points from Remnant / Native DB
#rnDB <- list(site=read.csv("Remnant_native_SOC_DBv1_site.csv"), horizon=read.csv("Remnant_native_SOC_DBv1_horizon.csv"))
#rnDB$site$SOCS = rnDB$site$Reported.100.cm.SOC/10
#str(rnDB)

## IMPORT RASTERS INTO ONE BIG STACK -----------

hyde.cl = c("tot_rice", "tot_rainfed", "tot_irri", "rf_rice", "rf_norice", "rangeland", "pasture", "ir_rice", "grazing", "cropland")
cru.cl = c("cru_pre_clim", "cru_tmp_clim", "cru_tmx_clim","cru_tmn_clim")
## _1901-1930
selP = c("C05GLC5","DEMMRG5","SLPMRG5","TWIMRG5","MNGUSG","MODCF_monthlymean","USG5","ESA4","average_soil","LandCover_2013001","ifl_", sapply(hyde.cl, function(x){paste0(x,900)}), sapply(hyde.cl, function(x){paste0(x,1800)}), sapply(hyde.cl, function(x){paste0(x,1910)}), sapply(hyde.cl, function(x){paste0(x,1960)}), sapply(hyde.cl, function(x){paste0(x,1990)}), sapply(hyde.cl, function(x){paste0(x,1910)}), sapply(hyde.cl, function(x){paste0(x,2016)}), sapply(hyde.cl, function(x){paste0(x,2010)}), sapply(cru.cl, function(x){paste0(x,"_1901-1930")}), sapply(cru.cl, function(x){paste0(x,"_1961-1990")}), "WDPA", "ofc_gen_10km", "cfc_gen_10km", "giems_d15_v10")
## 103 types
tot.tif <- list.files(path="./stacked", pattern = glob2rx("*.tif$"), full.names = TRUE)
## 1236 images
tifP = tot.tif[unlist(sapply(selP, function(x){grep(x, tot.tif)}))]
str(tifP)
## 242 layers
g10kmP = readGDAL(tifP[1])
x = parallel::mclapply(1:length(tifP), function(j){ readGDAL(tifP[j], silent=TRUE)$band1 }, mc.cores=56)
g10kmP@data = data.frame(x)
## replace '-' symbols otherwise reports problems with column names:
names(g10kmP) = gsub("\\-", "\\.", basename(file_path_sans_ext(tifP)))
rm(x)
#names(g10kmP) = basename(file_path_sans_ext(tifP[1]))
#for(j in 2:length(tifP)){ 
#  g10kmP@data[,basename(file_path_sans_ext(tifP[j]))] = readGDAL(tifP[j], silent=TRUE)$band1 
#}
#plot(raster(g10kmP["cru_tmp_clim_1901-1930_03_10km"]), col=SAGA_pal[[1]])
#plot(raster(g10kmP["cfc_gen_10km"]), col=SAGA_pal[[1]])
g10kmP$cfc_gen_10km <- ifelse(is.na(g10kmP$cfc_gen_10km), 6, g10kmP$cfc_gen_10km)
g10kmP$ofc_gen_10km <- ifelse(is.na(g10kmP$ofc_gen_10km), 6, g10kmP$ofc_gen_10km)
g10kmP$cfc_gen_10km <- as.factor(g10kmP$cfc_gen_10km)
g10kmP$ofc_gen_10km <- as.factor(g10kmP$ofc_gen_10km)
summary(g10kmP$ofc_gen_10km)
#plot(raster(g10kmP["cfc_gen_10km"]), col=SAGA_pal[[1]])
## 10.6GB object
g10kmP$ifl_2000_10km <- readGDAL("./stacked/ifl_2000_10km.tif")$band1
#g10kmP$cropland2009AD_10km <- readGDAL("./stacked/cropland2009AD_10km.tif")$band1
summary(as.factor(g10kmP$LandCover_2013001_L1_10km))
g10kmP$MASK <- ifelse(is.na(g10kmP$LandCover_2013001_L1_10km) | g10kmP$LandCover_2013001_L1_10km == 15 | g10kmP$LandCover_2013001_L1_10km == 0, NA, 1)
summary(g10kmP$MASK)
#writeGDAL(g10kmP["MASK"], "landmask_10km.tif", type="Byte", mvFlag=255, options = "COMPRESS=DEFLATE")
writeGDAL(g10kmP["MASK"], "landmask_10km.sdat", type="Byte", mvFlag=255, drivername = "SAGA")
g10kmP$s.intact = !is.na(g10kmP$ifl_2000_10km) | (!is.na(g10kmP$'WDPA_Apr2016.shapefile.polygons_10km') & g10kmP$cropland2010AD_10km<10)
summary(g10kmP$s.intact)
## 602,000 pixels in mask
g10mP.mask <- g10kmP["s.intact"]
g10mP.mask <- as(g10mP.mask, "SpatialPixelsDataFrame")
g10mP.mask <- g10mP.mask[g10mP.mask$s.intact==TRUE,]
#plot(raster(g10mP.mask))
#lines(country)
g10mP.mask@data[,1] <- as.numeric(g10mP.mask@data[,1])
writeGDAL(g10mP.mask[1], "intact_mask_10km.tif", type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
gc(); gc()

## Filter out missing values (latitudes >65 degrees N)
rnd = spsample(g10kmP["MASK"], n=5e3, type = "random")
ovRND = over(y=g10kmP, x=spTransform(rnd, CRS(proj4string(g10kmP))))
ovRND = ovRND[!is.na(ovRND$MASK),]
sel.na <- colSums(sapply(ovRND, function(x){!is.na(x)}))
sel.na[which(sel.na<1000)]

## 1. Historic soil organic carbon stock  -----------
## modelled using point data
## m.OCS = f ( climate, relief, surface geology, land cover )

## overlay existing profiles and 10km covs
ovM2 <- over(y=g10kmP, x=spTransform(OCD_profs, CRS(proj4string(g10kmP))))
save.image()

#plot(raster(g10kmP["cru_tmn_clim_1901.1930_06_10km"]), col=SAGA_pal[[1]])
## select columns for model building (ONLY COVS representing CURRENT CONDITIONS):
pred.selC = names(g10kmP)[c(grep("MRG5", names(g10kmP)), grep("ESA4", names(g10kmP)), grep("MODCF_monthlymean", names(g10kmP)), grep("USG", names(g10kmP)), grep("2010AD", names(g10kmP)), grep("clim_1961", names(g10kmP)), grep("cfc_gen", names(g10kmP)), grep("giems_d15", names(g10kmP)))] 
#pred.selO = names(g10kmP)[c(grep("MRG5", names(g10kmP)), grep("USG", names(g10kmP)), grep("900AD", names(g10kmP)), grep("clim_1901", names(g10kmP)), grep("ofc_gen", names(g10kmP)))]
## without climate as dynamic variable
str(pred.selC)
## 106 in total
#ovMC2 <- cbind(as.data.frame(SOCS_profs[c("SOURCEID","dSOCS_30cm","dSOCS_100cm","dSOCS_200cm")]), ovM2[,pred.selC])
ovMC2 <- cbind(as.data.frame(OCD_profs[c("SOURCEID","DEPTH.f","OCDENS")]), ovM2[,pred.selC])
names(ovMC2) = gsub("2010AD", "", gsub("1961.1990_", "", names(ovMC2)))

## Replace factors with indicators:
in1km = data.frame(model.matrix(~cfc_gen_10km-1, ovMC2))
ovMC2 = cbind(ovMC2, in1km)
#str(ovMC2)
predT.sel = names(ovMC2)[grep(pattern="_10km", names(ovMC2))]

## For simulated points make sure that OCS is 0 for HYDE = 0:
sel.sim = grep("SIM_", paste(ovMC2$SOURCEID))
#View(ovMC2[sel.sim,(c("SOURCEID","OCDENS","grazing_10km", "cropland_10km"))])
ovMC2[sel.sim,"OCDENS"] = ifelse(ovMC2[sel.sim,"grazing_10km"]==0&ovMC2[sel.sim,"cropland_10km"]==0,0,NA)
saveRDS.gz(ovMC2, file="regMatrix_OCD_global_profs.rds")

unlink("mrf.OCD_2010.rds")
unlink("mgb.OCD_2010.rds")
library(xgboost)
library(ranger)
library(caret)
## Fit models / sub-sample to speed up model fitting ----
Nsub <- 1e4 
formulaString.OCD <- as.formula(paste0('OCDENS ~ DEPTH.f + ', paste(predT.sel[-which(predT.sel=="cfc_gen_10km")], collapse="+")))
## Initiate cluster
require(parallel)
cl <- makeCluster(56)
registerDoParallel(cl)
cat("Results of model fitting 'randomForest / XGBoost':\n\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"))
cat("\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
cat(paste("Variable:", all.vars(formulaString.OCD)[1]), file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
cat("\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
LOC_ID <- ovMC2$SOURCEID
## Caret training settings (reduce number of combinations to speed up):
ctrl <- trainControl(method="repeatedcv", number=3, repeats=1)
gb.tuneGrid <- expand.grid(eta = c(0.3,0.4,0.5), nrounds = c(50,100,150), max_depth = 2:3, gamma = 0, colsample_bytree = 0.8, min_child_weight = 1)
rf.tuneGrid <- expand.grid(mtry = seq(10,60,by=5))
out.rf <- paste0("mrf.OCD_2010.rds")
if(!file.exists(out.rf)){
  dfs <- ovMC2[,all.vars(formulaString.OCD)]
  sel <- complete.cases(dfs)
  dfs <- dfs[sel,]
  if(nrow(dfs)<Nsub){Nsub=nrow(dfs)}
  ## optimize mtry parameter:
  if(!file.exists(gsub("mrf","t.mrf",out.rf))){
    t.mrfX <- caret::train(formulaString.OCD, data=dfs[sample.int(nrow(dfs), Nsub),], method="ranger", trControl=ctrl, tuneGrid=rf.tuneGrid)
    saveRDS.gz(t.mrfX, file=gsub("mrf","t.mrf",out.rf))
  } else {
    t.mrfX <- readRDS.gz(gsub("mrf","t.mrf",out.rf))
  }
  ## fit RF model using 'ranger' (fully parallelized)
  mrfX <- ranger(formulaString.OCD, data=dfs, importance="impurity", write.forest=TRUE, mtry=t.mrfX$bestTune$mtry, num.trees=85)  
  saveRDS.gz(mrfX, file=out.rf)
  ## Top 15 covariates:
  sink(file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE, type="output")
  print(mrfX)
  cat("\n Variable importance:\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
  xl <- as.list(ranger::importance(mrfX))
  print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:40]])))
  ## save fitting success vectors:
  fit.df <- data.frame(LOC_ID=LOC_ID[sel], observed=dfs[,1], predicted=predictions(mrfX))
  unlink(paste0("RF_fit_OCD.csv.gz"))
  write.csv(fit.df, paste0("RF_fit_OCD.csv"))
  gzip(paste0("RF_fit_OCD.csv"))
  mg.out = paste0("mgb.OCD_2010.rds")
  if(!file.exists(mg.out)){
    ## fit XGBoost model (uses all points):
    mgbX <- caret::train(formulaString.OCD, data=dfs, method="xgbTree", trControl=ctrl, tuneGrid=gb.tuneGrid) 
    saveRDS.gz(mgbX, file=mg.out)
  }
  importance_matrix <- xgb.importance(mgbX$coefnames, model=mgbX$finalModel)
  cat("\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
  print(mgbX)
  cat("\n XGBoost variable importance:\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
  print(importance_matrix[1:40,])
  cat("--------------------------------------\n", file=paste0("Potential_OCD_resultsFit_LandUse_only.txt"), append=TRUE)
  sink()
}
stopCluster(cl); 
closeAllConnections()
#rm(mrfX); rm(mgbX)
mrfX
## 55%
mgbX
## 45%

## Predict OCD values using current and past climate/land cover/land use
cfc.levs = levels(ovMC2$cfc_gen_10km)

predict_e <- function(mrfX, mgbX=NULL, newdata, cfc.levs, depth=100){
  sel.comp = complete.cases(newdata@data) & !is.na(newdata$MASK)
  nc = which(names(newdata)=="ofc_gen_10km")
  if(length(nc)>0){ names(newdata)[nc] = "cfc_gen_10km" }
  for(k in cfc.levs){
    newdata@data[,paste0("cfc_gen_10km",k)] = ifelse(newdata@data[,"cfc_gen_10km"]==k, 1, 0)
  }
  newdata$DEPTH.f = depth
  v1 <- predict(mrfX, newdata@data[sel.comp,])$predictions
  m <- newdata["cfc_gen_10km"]
  m <- as(m, "SpatialPixelsDataFrame")
  if(!is.null(mgbX)){
    gm1.w = 1/mrfX$prediction.error
    gm2.w = 1/(min(mgbX$results$RMSE, na.rm=TRUE)^2)
    v2 <- predict(mgbX, newdata@data[sel.comp,])
    m@data[sel.comp,"predicted"] <- data.frame(Reduce("+", list(v1*gm1.w, v2*gm2.w)) / (gm1.w+gm2.w))*10 ## 10 x kg/m-cubic
  } else {
    m@data[sel.comp,"predicted"] <- v1
  }
  return(m["predicted"])
}

DepthI = c(0,30,100,200)
unlink(list.files(pattern=glob2rx(paste0("OCD_", DepthI, "cm_year_*_10km.tif$"))))
#mrfX = readRDS.gz(file="mrf.OCD_2010.rds")
#mgbX = readRDS.gz(file="mgb.OCD_2010.rds")
periods = c("NoLU","900AD","1800AD","1910AD","1960AD","1990AD","2010AD","2016AD")
ofc.lst = c("ofc_gen","ofc_gen","ofc_gen","cfc_gen","cfc_gen","cfc_gen","cfc_gen","cfc_gen")

## run in loop - TAKES CA 30 MINS
for(j in 1:length(periods)){
  pred.selC = c(names(g10kmP)[c(grep("MRG5", names(g10kmP)), grep("ESA4", names(g10kmP)), grep("MODCF_monthlymean", names(g10kmP)), grep("USG", names(g10kmP)), grep(periods[j], names(g10kmP)), grep("clim_1961", names(g10kmP)), grep(ofc.lst[j], names(g10kmP)), grep("giems_d15", names(g10kmP)))], "MASK")
  g10kmP_current <- g10kmP[pred.selC]
  names(g10kmP_current) = gsub(periods[j], "", gsub("1961.1990_", "", names(g10kmP_current)))
  ## Special case = reduce all Land use to 0
  if(periods[j]=="NoLU"){
    g10kmP_current@data[,c("tot_rice_10km","tot_rainfed_10km","tot_irri_10km", "rf_rice_10km","rf_norice_10km","rangeland_10km","pasture_10km","ir_rice_10km","grazing_10km","cropland_10km")] = 0
  } else {
    ## Filter out all LU values smaller than 3% (HYDE is often not that precise):
    for(i in c("tot_rice_10km","tot_rainfed_10km","tot_irri_10km", "rf_rice_10km","rf_norice_10km","rangeland_10km","pasture_10km","ir_rice_10km","grazing_10km","cropland_10km")){
      g10kmP_current@data[,i] <- ifelse(g10kmP_current@data[,i]<3, 0, g10kmP_current@data[,i])
    }
  }
  for(k in DepthI){
    out.tif = paste0("OCD_", k, "cm_year_", periods[j], "_10km.tif")
    if(!file.exists(out.tif)){
      ## predict OCD
      pr.OCD_c <- predict_e(mrfX, mgbX, newdata=g10kmP_current, cfc.levs=cfc.levs, depth=k)
      writeGDAL(pr.OCD_c[1], out.tif, type="Int16", mvFlag=-9999, options="COMPRESS=DEFLATE")
    }
  }
}

## Derive cumulative SOCS for 0-2 m:
sum_SOCS = function(tif.lst, depthT = c(30,70,100), year, depth.sel=200){
  out.tif = paste0("SOCS_0_", depth.sel, "cm_year_", year, "_10km.tif")
  s = stack(tif.lst)
  s = as(s, "SpatialGridDataFrame")
  for(i in 1:ncol(s)){ s@data[,i] = ifelse(s@data[,i]<0, 0, s@data[,i]) }
  x = list(NULL)
  for(i in 1:(length(tif.lst)-1)){
    x[[i]] = rowMeans(s@data[,i:(i+1)], na.rm=TRUE)*depthT[i]/100 
  }
  for(k in 1:length(depth.sel)){
    if(depth.sel[k]==200){
      s$SOCS = rowSums(as.data.frame(x[1:3]), na.rm=TRUE)
    }
    if(depth.sel[k]==100){
      s$SOCS = rowSums(as.data.frame(x[1:2]), na.rm=TRUE)
    }
    if(depth.sel[k]==30){
      s$SOCS = rowSums(as.data.frame(x[1]), na.rm=TRUE)
    }
    ## tones / ha
    writeGDAL(s["SOCS"], out.tif[k], type="Int16", mvFlag=0, options="COMPRESS=DEFLATE")
  }
}

tif.lst <- lapply(periods, function(x){paste0("OCD_",c(0,30,100,200),"cm_year_",x,"_10km.tif")})
sfInit(parallel=TRUE, cpus=length(periods))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("sum_SOCS", "tif.lst", "periods")
missing.lst <- sfLapply(1:length(periods), function(i){sum_SOCS(tif.lst[[i]], year=periods[i])})
sfStop()

## for NoLU and 2010:
sum_SOCS(tif.lst[[1]], year=periods[1], depth.sel=c(30,100,200))
sum_SOCS(tif.lst[[7]], year=periods[7], depth.sel=c(30,100,200))

## plot difference:
SOC_10km = stack(c("SOCS_0_200cm_year_NoLU_10km.tif","SOCS_0_200cm_year_2016AD_10km.tif"))
SOC_10km = as(SOC_10km, "SpatialGridDataFrame")
#plot(stack(SOC_10km))
SOC_10km$dif = (SOC_10km@data[,1]-SOC_10km@data[,2])
plot(SOC_10km["dif"], zlim=c(-50,200), col=SAGA_pal[[1]])
writeGDAL(SOC_10km["dif"], "SOCS_0_200cm_year_difference_10km.tif", options = c("COMPRESS=DEFLATE"))
unlink("SOCS_0_200cm_year_difference_10km_xy.tif")
system('gdalwarp SOCS_0_200cm_year_difference_10km.tif SOCS_0_200cm_year_difference_10km_xy.tif -t_srs \"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs\" -te -16810131 -6625155 16810131 8343004 -co \"COMPRESS=DEFLATE\"')
system('gdalwarp ./stacked/cropland2016AD_10km.tif cropland2016AD_10km_xy.tif -t_srs \"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs\" -te -16810131 -6625155 16810131 8343004 -co \"COMPRESS=DEFLATE\"')
system('gdalwarp ./stacked/pasture2016AD_10km.tif pasture2016AD_10km_xy.tif -t_srs \"+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs\" -te -16810131 -6625155 16810131 8343004 -co \"COMPRESS=DEFLATE\"')
download.file(url="http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip", "ne_110m_admin_0_countries.zip", "auto")
unzip("ne_110m_admin_0_countries.zip")
world <- readOGR(".", "ne_110m_admin_0_countries")
world <- spTransform(world, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
world@bbox
image(raster("SOCS_0_200cm_year_difference_10km_xy.tif"), zlim=c(-20,120), col=SAGA_pal[[1]])
lines(as(world, "SpatialLines"))

## Plot variable importance:
xl = as.list(ranger::importance(mrfX))
xl = t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[2:26]]))

pdf(file = "Fig_RF_importance_plot_200cm.pdf", width = 7, height = 7.5)
par(mar=c(2.5,9,2.5,0.5), oma=c(1,1,1,1))
plot(x=rev(xl)/max(xl)*100, y=1:25, pch = 19, col="blue", xlab="Importance (%)", xlim=c(0,105), ylim=c(0,26), yaxp=c(0,25,25), xaxs="i", yaxs="i", cex=1.4, yaxt="n", ylab="", main="SOC density model importance plot", cex.main=1)
abline(h=1:25, lty=2, col="grey")
#axis(2, at=1:25, labels=rev(attr(xl, "dimnames")[[1]]), las=2)
axis(2, at=1:25, labels=rev(c("Max. temp. September", "Max. temp. August", "MCF October", "MCF December", "Elevation", "MCF September", "MCF January", "MCF November", "Max. temp. October", "MCF June", "Wetness Index", "MCF May", "MCF July", "MCF February", "MCF March", "MCF August", "MCF April", "Snow prob. March", "Precipitation October", expression(bold("Grazing (HYDE)")), "Max. temp. April", "Hillands class", expression(bold("Cropland (HYDE)")), "Max. temp. February", expression(bold("Total rainfed (HYDE)")))), las=2)
dev.off()

## plot 1 to 1 relationships -----
library(scales)
library(hexbin)
pfun <- function(x,y, ...){
  panel.hexbinplot(x,y, ...)
  panel.loess(x, y, ..., col = "black",lty=1,lw=2,span=1/18)
}
pal = R_pal[["bpy_colors"]][1:18]
#par(mfcol=c(1,2), oma=c(1,1,1,1))
p3 = hexbinplot(ovMC2$OCDENS~ovMC2$cru_tmx_clim_10_10km, colramp=colorRampPalette(pal), ylab="Organic carbon density (kg/cubic-m)", xlab="Grazing (HYDE)", type="g", ylim=c(0,120), lwd=1, lcex=8, inner=.2, cex.labels=.8, xbins=30, asp=1, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)) ## , panel=pfun
p4 = hexbinplot(ovMC2$OCDENS~ovMC2$DEMMRG5_agg_ll_10km, colramp=colorRampPalette(pal), ylab="Organic carbon density (kg/cubic-m)", xlab="Grazing (HYDE)", type="g", ylim=c(0,120), lwd=1, lcex=8, inner=.2, cex.labels=.8, xbins=30, asp=1, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1))
p1 = hexbinplot(ovMC2$OCDENS~ovMC2$grazing_10km, colramp=colorRampPalette(pal), ylab="Organic carbon density (kg/cubic-m)", xlab="Grazing (HYDE)", type="g", ylim=c(0,120), xlim=c(0,100), lwd=1, lcex=8, inner=.2, cex.labels=.8, xbins=30, asp=1, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1)) 
p2 = hexbinplot(ovMC2$OCDENS~ovMC2$cropland_10km, colramp=colorRampPalette(pal), ylab="", xlab="Cropland (HYDE)", type="g", ylim=c(0,120), xlim=c(0,100), lwd=1, lcex=8, inner=.2, cex.labels=.8, xbins=30, asp=1, colorcut=c(0,0.01,0.03,0.07,0.15,0.25,0.5,0.75,1))
library(gridExtra)    
pdf(file = "Fig_hexbinplots_correlations_HYDE.pdf", width = 8, height = 4, pointsize=14)
par(oma=c(1,1,1,1))
#do.call(grid.arrange, c(list(p3,p4,p1,p2), ncol=2))
do.call(grid.arrange, c(list(p1,p2), ncol=2))
dev.off()

library(leaflet)
library(htmlwidgets)
unlink("SOCS_0_200cm_year_difference_10km_plt.tif")
system(paste0(gdalwarp, ' SOCS_0_200cm_year_difference_10km.tif SOCS_0_200cm_year_difference_10km_plt.tif -co \"COMPRESS=DEFLATE\" -overwrite -s_srs EPSG:4326 -t_srs EPSG:3857 -tr 10000 10000 -multi -of GTiff -te -20037508 -7706358 20022492 15506358')) ## \"+init=epsg:3857\"
unlink("SOCS_0_200cm_year_NoLU_10km_plt.tif")
system(paste0(gdalwarp, ' SOCS_0_200cm_year_NoLU_10km.tif SOCS_0_200cm_year_NoLU_10km_plt.tif -co \"COMPRESS=DEFLATE\" -overwrite -s_srs EPSG:4326 -t_srs EPSG:3857 -tr 10000 10000 -multi -of GTiff -te -20037508 -7706358 20022492 15506358'))
system(paste0(gdalwarp, ' SOCS_0_200cm_year_2016AD_10km.tif SOCS_0_200cm_year_2016AD_10km_plt.tif -co \"COMPRESS=DEFLATE\" -overwrite -s_srs EPSG:4326 -t_srs EPSG:3857 -tr 10000 10000 -multi -of GTiff -te -20037508 -7706358 20022492 15506358'))

#plot(raster("OCS_2m_10km.tif"))
r = readGDAL("SOCS_0_200cm_year_difference_10km_plt.tif")
r$band1 <- ifelse(r$band1< -27, -27, ifelse(r$band1>85, 85, r$band1))
#r$band1 <- ifelse(r$band1<0, 0, ifelse(r$band1>450, 450, r$band1))
summary(r$band1)
r = raster(r)
pal <- colorNumeric(SAGA_pal[[1]], values(r), na.color = "transparent")
m1 <- leaflet() %>% addTiles() %>% addRasterImage(r, colors=pal, opacity=0.6, project=FALSE, maxBytes = 6 * 1024 * 1024) %>% addLegend(pal=pal, values=values(r), title="SOCS 0--200 cm in t/ha (difference)")
saveWidget(m1, file="SOCS_0_200cm_year_difference_10km.html")

r2 = readGDAL("SOCS_0_200cm_year_NoLU_10km_plt.tif")
summary(r2$band1)
r2$band1 <- ifelse(r2$band1<0, 0, ifelse(r2$band1>680, 680, r2$band1))
r2 = raster(r2)
pal2 <- colorNumeric(SAGA_pal[[1]][5:20], values(r2), na.color = "transparent")
m2 <- leaflet() %>% addTiles() %>% addRasterImage(r2, colors=pal2, opacity=0.6, project=FALSE, maxBytes = 6 * 1024 * 1024) %>% addLegend(pal=pal2, values=values(r2), title="SOCS 0---200 cm in t/ha (no LU)")
saveWidget(m2, file="SOCS_0_200cm_year_NoLU_10km.html")
r3 = readGDAL("SOCS_0_200cm_year_2016AD_10km_plt.tif")
r3$band1 <- ifelse(r3$band1<0, 0, ifelse(r3$band1>680, 680, r3$band1))
r3 = raster(r3)
pal3 <- colorNumeric(SAGA_pal[[1]][5:20], values(r3), na.color = "transparent")
m3 <- leaflet() %>% addTiles() %>% addRasterImage(r3, colors=pal3, opacity=0.6, project=FALSE, maxBytes = 6 * 1024 * 1024) %>% addLegend(pal=pal3, values=values(r3), title="SOCS 0---200 cm in t/ha (2016AD)")
saveWidget(m3, file="SOCS_0_200cm_year_2016AD_10km.html")

save.image()

## comparison histograms:
ca.sg = stack(list("OCS_100cm_historic_10km_xy.tif","OCS_100cm_current_10km_xy.tif"))
ca.sg = as(ca.sg, "SpatialGridDataFrame")
str(ca.sg@data)
#lm(ca.sg$OCS_100cm_historic_10km~ca.sg$OCS_100cm_current_10km, ca.sg@data[!is.na(ca.sg$OCS_100cm_historic_10km),])
with(ca.sg@data[sample.int(length(ca.sg$OCS_100cm_historic_10km_xy),20000),], psych::scatter.hist(log1p(OCS_100cm_historic_10km_xy),log1p(OCS_100cm_current_10km_xy), xlab="historic", ylab="current", pch=19, title="Organic carbon stock (0--100 cm)", col=alpha("lightblue", 0.4), cex=1.5))
## difference in current and historic OCS total:
ca.sg = as(ca.sg, "SpatialPixelsDataFrame")
ratioOCS = ca.sg@data$OCS_100cm_historic_10km_xy/ca.sg@data$OCS_100cm_current_10km_xy*100
hist(ratioOCS[ratioOCS<300 & ratioOCS>0 & !is.na(ratioOCS)], breaks=seq(0,300,by=10), col="grey", xlim=c(0,300), xlab="Difference in percent", main="Historic / current SOCS ratio")
quantile(ratioOCS[ratioOCS<300 & ratioOCS>0 & !is.na(ratioOCS)], c(.1,.9))
mean(ratioOCS[ratioOCS<450 & ratioOCS>0 & !is.na(ratioOCS)])

## Prediction error for RF -----
## based on an empirical solution explained in: https://github.com/imbs-hl/ranger/issues/136

predict_cv_resid = function(formulaString, data, nfold, coords=c("LONWGS84", "LATWGS84")){
  data <- data[complete.cases(data[,all.vars(formulaString)]),]
  ## only one point per profile to avoid auto-correlation problems
  sel <- dismo::kfold(data, k=nfold, by=data$SOURCEID)
  out = list(NULL)
  for(j in 1:nfold){
    s.train <- data[!sel==j,]
    s.test <- data[sel==j,]
    m <- ranger(formulaString, data=s.train, write.forest=TRUE)
    pred <- predict(m, s.test, na.action = na.pass)$predictions
    obs.pred <- as.data.frame(list(s.test[,all.vars(formulaString)[1]], pred))
    names(obs.pred) = c("Observed", "Predicted")
    obs.pred[,"ID"] <- row.names(s.test)
    obs.pred$fold = j
    obs.pred[,coords] = s.test[,coords]
    out[[j]] = obs.pred
  }
  out <- plyr::rbind.fill(out)
  out <- out[order(as.numeric(out$ID)),]
  return(out)
}

## TAKES 30 mins
resid.OCDENS = predict_cv_resid(formulaString.OCD, data=ovMC2, nfold=4)
#xyplot(Predicted~Observed, resid.OCDENS, asp=1, par.settings=list(plot.symbol = list(col=alpha("black", 0.6), fill=alpha("red", 0.6), pch=21, cex=0.9)), xlab="measured", ylab="predicted (ranger)")
ovMC2$resid = NA
ovMC2[complete.cases(ovMC2[,all.vars(formulaString.OCD)]),"resid"] = abs(resid.OCDENS$Observed - resid.OCDENS$Predicted)

## model absolute residuals as function of covariates:
var.fm2 = as.formula(paste("resid ~ ", paste0(all.vars(formulaString.OCD)[-1], collapse = "+")))
var.OCDENS.rf <- ranger(var.fm2, data = ovMC2[complete.cases(ovMC2[,all.vars(var.fm2)]),], write.forest = TRUE, num.trees=85, mtry = mrfX$mtry)
var.OCDENS.rf
## R-square = 35% 

## predict uncertainty:
for(d in c(0,30,100,200)){
  pred.selC = c(names(g10kmP)[c(grep("MRG5", names(g10kmP)), grep("ESA4", names(g10kmP)), grep("MODCF_monthlymean", names(g10kmP)), grep("USG", names(g10kmP)), grep("2010AD", names(g10kmP)), grep("clim_1961", names(g10kmP)), grep(ofc.lst[j], names(g10kmP)), grep("giems_d15", names(g10kmP)))], "MASK")
  g10kmP_current <- g10kmP[pred.selC]
  names(g10kmP_current) = gsub("2010AD", "", gsub("1961.1990_", "", names(g10kmP_current)))
  pr.OCD_var <- predict_e(var.OCDENS.rf, newdata=g10kmP_current, cfc.levs=cfc.levs, depth=d)
  writeGDAL(pr.OCD_var[1], paste0("abs_error_OCDENS_",d,"cm_2010AD_10km_ll.tif"), type="Int16", mvFlag=-9999, options="COMPRESS=DEFLATE")
}

## All covariates:
#system("7za a SOCS_global_Covariates_10km.7z ./stacked/*.tif")
