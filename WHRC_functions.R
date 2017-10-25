## Derivation of potential soil carbon (https://github.com/whrc/Soil-Carbon-Debt)
## Code by: Tom.Hengl@isric.org

gdalwarp_clim <- function(x, srcnodata = "254", dstnodata = "255"){
  out = paste0("./stacked/", gsub(".tif", "_10km.tif", basename(x)))
  if(!file.exists(out)){ system(paste0(gdalwarp, ' ', x, ' ', out, ' -srcnodata \"', srcnodata,'\" -dstnodata \"', dstnodata,'\" -co \"COMPRESS=DEFLATE\" -r \"cubicspline\" -tr ', cellsize, ' ', cellsize))} ##  -t_srs \"+proj=longlat +datum=WGS84\" ' -te ', xllcorner,' ', yllcorner, ' ', xurcorner, ' ', yurcorner
}

## function to strip year:
strip_year = function(x, name, ext=".asc"){
  xn = sapply(paste(x), function(x){strsplit(strsplit(x, name)[[1]][2], ext)[[1]][1]})
  xi = rep(NA, length(xn))
  xi[grep("AD", xn)] <- as.numeric( sapply(paste(xn[grep("AD", xn)]), function(x){strsplit(x, "AD")[[1]][1]}) )
  xi[grep("BC", xn)] <- -as.numeric( sapply(paste(xn[grep("BC", xn)]), function(x){strsplit(x, "BC")[[1]][1]}) )
  return(xi)
}

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

## rasterize shapefiles:
rasterize_pol <- function(INPUT, FIELD, cellsize, xllcorner, yllcorner, xurcorner, yurcorner){
  out = paste0(strsplit(basename(INPUT), "\\.")[[1]][1], ".sdat")
  ## "./stacked10km/", 
  if(!file.exists(out)){
    system(paste0('/usr/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"', INPUT, '\" -FIELD \"', FIELD, '\" -GRID \"', gsub(".sdat", ".sgrd", out), '\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
  }
}

## Fill-in missing pixels
landmask_fix = function(OCS.lst, s.year, landmask.file="./OCD/landmask_10km.tif", out.dir="./SOCS/"){
  s = raster::stack(c(OCS.lst[grep(s.year, OCS.lst)], OCS.lst[grep("NoLU", OCS.lst)], landmask.file))
  s = as(s, "SpatialGridDataFrame")
  ## fill in all missing values based on the dominant value:
  for(i in c("30cm","100cm","200cm")){
    s@data[,"fix"] = ifelse(is.na(s$landmask_10km), NA, ifelse(s@data[,paste0("SOCS_0_",i ,"_year_",s.year,"_10km")]==0, s@data[,paste0("SOCS_0_",i ,"_year_NoLU_10km")], s@data[,paste0("SOCS_0_",i ,"_year_",s.year,"_10km")]))
    writeGDAL(s["fix"], paste0(out.dir, "SOCS_0_",i ,"_year_",s.year,"_10km.tif"), type="Int16", mvFlag=-32767, options="COMPRESS=DEFLATE")
  }
}

## Prediction error for RF ----
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


## Potential Vegetation
predict_PVG_tiles <- function(i, gm, tile.tbl, covs, levs, method="ranger", out.path="/data/WHRC_soilcarbon/model10km/tiled", lc1="/data/WHRC_soilcarbon/model10km/stacked1km/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1_1km.tif", des="/data/WHRC_soilcarbon/model10km/stacked1km/desertPR_sin_1km.tif", pvg.m){
  out.tifs <- gsub("-", ".", paste0(out.path, "/T", tile.tbl[i,"ID"], "/", normalizeFilename(levs, sub.sign = "."), "_T", tile.tbl[i,"ID"], ".tif"))
  if(any(!file.exists(out.tifs))){
    m = readGDAL(fname=lc1, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    names(m) = "landcover"
    m$desert = readGDAL(des, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
    m$desert = ifelse(is.na(m$desert), 0, m$desert)
    ## NO need to predict for permanent snow / shifting sands:
    m = as(m, "SpatialPixelsDataFrame")
    sel.p = !(m$landcover==220 | m$landcover==210 | m$desert==1)
    if(sum(sel.p)>2){
      m = m[sel.p,]
      for(j in 1:length(covs)){
        cname = gsub("_1km", "", gsub("-", "_", file_path_sans_ext(basename(covs[j]))))
        m@data[,cname] <- signif(readGDAL(covs[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index], 4)
      }
      ## Fill-in missing values (if necessary):
      sel.mis = sapply(m@data, function(x){sum(is.na(x))>0})
      if(sum(sel.mis)>0){
        for(j in which(sel.mis)){
          if(length(grep(pattern="USG5", names(m)[j]))>0){ 
            repn = 0 
          } else {
            repn = quantile(m@data[,j], probs=.5, na.rm=TRUE)
            if(is.na(repn)){
              repn = quantile(pvg.m[,names(m)[j]], probs=.5, na.rm=TRUE)
            }
          }
          m@data[,j] = ifelse(is.na(m@data[,j]), repn, m@data[,j])
        }
      }
      ## predict:
      if(method=="ranger"){
        m@data = data.frame(predict(gm, m@data, num.threads=1)$predictions)
        #plot(stack(m[6:8]), zlim=c(0,.35), col=SAGA_pal[[1]])
        for(j in 1:ncol(m)){
          out.tif <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/", names(m)[j], "_T", tile.tbl[i,"ID"], ".tif")
          m@data[,j] <- round(m@data[,j]*100)
          writeGDAL(m[j], out.tif, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
        }
        ## match most probable class
        m$cl <- apply(m@data,1,which.max)
        writeGDAL(m["cl"], paste0(out.path, "/T", tile.tbl[i,"ID"], "/PVGcl_T", tile.tbl[i,"ID"], ".tif"), type="Byte", mvFlag=255, options="COMPRESS=DEFLATE", catNames=list(gm$forest$levels))
      }
      ## Save output:
      saveRDS(m, file=paste0(out.path, "/T", tile.tbl[i,"ID"], "/m_T", tile.tbl[i,"ID"], ".rds"))
      gc()
    }
  }
}

## Predict SOCS using 2 models:
predict_SOCS_tiles <- function(i, tile.tbl, varn, n.tif, out.path="/data/WHRC_soilcarbon/model10km/tiled", gmF, gmX, from.tile=TRUE, lc.class=NULL, lc.name="ESACCI_LC_L4_LCCS_Map_300m_P5Y_2010_v1.6.1", gmF.w=.5, gmX.w=.5, lc1="/data/WHRC_soilcarbon/model10km/stacked1km/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1_1km.tif", des="/data/WHRC_soilcarbon/model10km/stacked1km/desertPR_sin_1km.tif", reg.m){
  out.tif <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/", varn, "_T", tile.tbl[i,"ID"], ".tif")
  if(!file.exists(out.tif)){
    if(from.tile==TRUE){
      m = readGDAL(paste0(out.path, "/T", tile.tbl[i,"ID"], "/", strsplit(basename(n.tif[1]),"_")[[1]][3], "_T", tile.tbl[i,"ID"], ".tif"), silent = TRUE)
      m = as(m, "SpatialPixelsDataFrame")
      names(m) = file_path_sans_ext(basename(n.tif[1]))
    } else {
      m = readGDAL(fname=lc1, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
      names(m) = "landcover"
      m$desert = readGDAL(fname=des, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
      m$desert = ifelse(is.na(m$desert), 0, m$desert)
      ## NO need to predict for permanent snow / shifting sands:
      m = as(m, "SpatialPixelsDataFrame")
      sel.p = !(m$landcover==220 | m$landcover==210 | m$desert==1)
      m = m[sel.p,]
    }
    for(j in 1:length(n.tif)){
      if(from.tile==TRUE){
        m@data[,file_path_sans_ext(basename(n.tif[j]))] <- readGDAL(paste0(out.path, "/T", tile.tbl[i,"ID"], "/", strsplit(basename(n.tif[j]),"_")[[1]][3], "_T", tile.tbl[i,"ID"], ".tif"), silent = TRUE)$band1[m@grid.index]
      } else {
        m@data[,file_path_sans_ext(basename(n.tif[j]))] <- readGDAL(fname=n.tif[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index]
      }
    }
    names(m) = gsub("-", "_", names(m))
    if(from.tile==FALSE){
      names(m) = gsub("_1km", "", names(m))
    }
    if(any(names(m) %in% lc.name)){
      m@data[,lc.name] = as.factor(m@data[,lc.name])
      for(k in lc.class){
        m@data[,paste0(lc.name,k)] <- ifelse(m@data[,lc.name]==k, 1, 0) 
      }
    }
    ## Fill-in missing values (if necessary):
    sel.mis = sapply(m@data, function(x){sum(is.na(x))>0})
    if(sum(sel.mis)>0){
      for(j in which(sel.mis)){
        if(length(grep(pattern="USG5", names(m)[j]))>0){ 
          repn = 0 
        } else {
          repn = quantile(m@data[,j], probs=.5, na.rm=TRUE)
          if(is.na(repn)){
            repn = quantile(reg.m[,grep(names(m)[j], names(reg.m))], probs=.5, na.rm=TRUE)
          }
        }
        m@data[,j] = ifelse(is.na(m@data[,j]), repn, m@data[,j])
      }
    }
    if(nrow(m)>2){
      if(is.null(gmX)){ 
        m$SOCS_m = predict(gmF, m@data)$predictions*10
        writeGDAL(m["SOCS_m"], out.tif, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      } else {
        ## Weighted mean / ensemble prediction (in tonnes / ha):
        if(missing(gmF.w)){ gmF.w = 1/gmF$prediction.error }
        if(missing(gmF.w)){ gmX.w = 1/(min(gmX$results$RMSE, na.rm=TRUE)^2) }
        m@data = data.frame(Reduce("+", list(predict(gmF, m@data)$predictions*gmF.w, predict(gmX, m@data)*gmX.w)) / (gmF.w+gmX.w)) * 10
        m@data[,1] <- ifelse(m@data[,1] < 0, 0, m@data[,1])
        writeGDAL(m[1], out.tif, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
      gc()
    }
  }
}


## Potential SOCS assuming potential vegetation
predict_potSOCS_tiles = function(i, tile.tbl, n.tif, out.path="/data/WHRC_soilcarbon/model10km/tiled", varn="OCS_0_100cm_PNV", pnv.df){
  out.tif <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/", varn, "_T", tile.tbl[i,"ID"], ".tif")
  if(!file.exists(out.tif)){
    m = readGDAL(paste0(out.path, "/T", tile.tbl[i,"ID"], "/", strsplit(basename(n.tif[1]),"_")[[1]][3], "_T", tile.tbl[i,"ID"], ".tif"), silent = TRUE)
    m = as(m, "SpatialPixelsDataFrame")
    names(m) = file_path_sans_ext(basename(n.tif[1]))
    for(j in 2:length(n.tif)){
      m@data[,file_path_sans_ext(basename(n.tif[j]))] <- readGDAL(paste0(out.path, "/T", tile.tbl[i,"ID"], "/", strsplit(basename(n.tif[j]),"_")[[1]][3], "_T", tile.tbl[i,"ID"], ".tif"), silent = TRUE)$band1[m@grid.index]
    }
    ## multiply by class centers (weighted average):
    m$SOCS_m = rowSums(t(t(m@data[,pnv.df$Biome_class]) * pnv.df$SOCS_m), na.rm=TRUE)/10
    writeGDAL(m["SOCS_m"], out.tif, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    gc()
  }
}

## aggregate per class / tile
summary_OCS_tiles <- function(i, tileS.tbl, out.path="./Stiled", lcs1="ESACCI_LC_L4_LCCS_2010_1km_sin.tif", sg1="OCS_0_100cm_current_1km_sin.tif", sg2="OCS_0_100cm_PNV_1km_sin.tif", cnt="GAUL_COUNTRIES_1km_sin.tif", lc.leg, countries, AREA = (1000*1000)/1e4){
  out.csv = paste0(out.path, "/T", tileS.tbl[i,"ID"], "/",c("SOCS_current_agg_LC2010","SOCS_PNV_agg_LC2010"),"_T", tileS.tbl[i,"ID"], ".csv") ## "SOCS_historic_agg_LC2010"
  out.tif = paste0(out.path, "/T", tileS.tbl[i,"ID"], "/OCS_0_100cm_",c("current","PNV"),"_T", tileS.tbl[i,"ID"], ".tif")
  if(any(!file.exists(out.csv))){
    m = readGDAL(fname=lcs1, offset=unlist(tileS.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    lst = list(sg1,sg2,cnt)
    for(j in 1:length(lst)){
      m@data[,j+1] = readGDAL(fname=lst[[j]], offset=unlist(tileS.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tileS.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
    }
    names(m) = c("LC2010","mSOCS_current","mSOCS_PNV","GAUL_COUNTRY")
    m = as(m, "SpatialPixelsDataFrame")
    m$LC2010 = factor(m$LC2010, levels=as.character(lc.leg$Value), labels=levels(lc.leg$NAME))
    m$GAUL_COUNTRY = factor(m$GAUL_COUNTRY, levels=as.character(1:nrow(countries)), labels=levels(countries$NAMES))
    m = m[!is.na(m$GAUL_COUNTRY),]
    saveRDS(m, paste0(out.path, "/T", tileS.tbl[i,"ID"], "/stacked_T", tileS.tbl[i,"ID"], ".rds"))
    writeGDAL(m["mSOCS_current"], out.tif[1], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    writeGDAL(m["mSOCS_PNV"], out.tif[2], type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
    ## Aggregate per country
    SOC_agg.LC2010 <- plyr::ddply(m@data, .(LC2010), summarize, Total_OCS_current_t=round(sum(mSOCS_current*AREA, na.rm = TRUE)/1e6,2), Sum_OCS=sum(mSOCS_current, na.rm=TRUE), N_OCS=sum(!is.na(mSOCS_current)))
    write.csv(SOC_agg.LC2010, out.csv[1])
    SOC_agg.LC2010 <- plyr::ddply(m@data, .(LC2010), summarize, Total_OCS_PNV_t=round(sum(mSOCS_PNV*AREA, na.rm = TRUE)/1e6,2), Sum_OCS=sum(mSOCS_PNV, na.rm=TRUE), N_OCS=sum(!is.na(mSOCS_PNV)))
    write.csv(SOC_agg.LC2010, out.csv[2])
  }
}

MRT_resample = function(INPUT_FILENAME, OUTPUT_FILENAME, SPECTRAL_SUBSET="1", OUTPUT_PROJECTION_TYPE = "SIN", OUTPUT_PROJECTION_PARAMETERS = "6371007.181 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0", RESAMPLING_TYPE = "BILINEAR", SPATIAL_SUBSET_UL_CORNER, SPATIAL_SUBSET_LR_CORNER, PIXEL_SIZE, MRT = '/data/MRT/MRT/bin/'){
  ## generate the prm file:
  prm = set.file.extension(INPUT_FILENAME, ".prm")
  filename = file(prm, open="wt")
  write(paste0('INPUT_FILENAME = \"', INPUT_FILENAME, '\"'), filename)
  write(paste0('SPECTRAL_SUBSET = (', SPECTRAL_SUBSET,')'), filename)
  write(paste0('SPATIAL_SUBSET_TYPE = INPUT_LAT_LONG'), filename)
  write(paste0('SPATIAL_SUBSET_UL_CORNER = (', paste0(SPATIAL_SUBSET_UL_CORNER, collapse=" "),')'), filename)
  write(paste0('SPATIAL_SUBSET_LR_CORNER = (', paste0(SPATIAL_SUBSET_LR_CORNER, collapse=" "),')'), filename)
  write(paste0('OUTPUT_FILENAME = \"', OUTPUT_FILENAME, '\"'), filename)
  write(paste0('RESAMPLING_TYPE = ', RESAMPLING_TYPE), filename)
  write(paste0('OUTPUT_PROJECTION_TYPE = ', OUTPUT_PROJECTION_TYPE), filename)
  ## -t projection_type [AEA ER GEO HAM IGH ISIN LA LCC MERCAT MOL PS SIN TM UTM]
  write(paste0('OUTPUT_PROJECTION_PARAMETERS = (', OUTPUT_PROJECTION_PARAMETERS, ')'), filename)
  write(paste0('DATUM = NoDatum'), filename)
  write(paste0('PIXEL_SIZE =  ', PIXEL_SIZE), filename)
  close(filename)
  system(paste0(MRT, 'resample -p ', prm))
}

latlon2sin = function(input.file, output.file, mod.grid, tmp.dir="./tmp/", proj, pixsize, cleanup.files=TRUE, te){
  ## reproject grid in tiles:
  out.files = paste0(tmp.dir, "T", mod.grid$ID, "_", basename(input.file))
  te.lst = apply(mod.grid@data[,1:4], 1, function(x){paste(x, collapse=" ")})
  if(missing(proj)){ proj = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs" }
  sfInit(parallel=TRUE, cpus=48)
  sfExport("mod.grid", "te.lst", "proj", "out.files")
  #sfLibrary(rgdal)
  x <- sfClusterApplyLB(1:length(out.files), function(i){ invisible( system(paste0('gdalwarp ', input.file, ' ', out.files[i], ' -r \"near\" -t_srs \"', proj, '\" -tr ', pixsize, ' ', pixsize, ' -te ', te.lst[i]), show.output.on.console = FALSE, intern = TRUE) ) }) ## -co \"COMPRESS=DEFLATE\"
  sfStop()
  ## mosaic:
  tmp.lst = list.files(path=tmp.dir, pattern=basename(input.file), full.names=TRUE)
  out.tmp <- tempfile(fileext = ".txt")
  vrt.tmp <- tempfile(fileext = ".vrt")
  cat(tmp.lst, sep="\n", file=out.tmp)
  system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
  if(missing(te)){
    system(paste0(gdalwarp, ' ', vrt.tmp, ' ', output.file, ' -ot \"Int16\" -dstnodata \"-32767\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\"'))
  } else {
    system(paste0(gdalwarp, ' ', vrt.tmp, ' ', output.file, ' -ot \"Int16\" -dstnodata \"-32767\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\" -te ', te))
  }
  if(cleanup.files==TRUE){ unlink(out.files) }
}

## Derive cumulative SOCS for 0-2 m ----
sum_SOCS = function(tifs, depthT = c(30,70,100), year, depth.sel=200){
  out.tif = paste0("./SOCS/SOCS_0_", c(depthT[1], sum(depthT[1:2]), sum(depthT[1:3])), "cm_year_", year, "_10km.tif")
  s = stack(tifs)
  s = as(s, "SpatialGridDataFrame")
  for(i in 1:ncol(s)){ s@data[,i] = ifelse(s@data[,i]<0, 0, s@data[,i]) }
  x = list(NULL)
  for(i in 1:(length(tifs)-1)){
    x[[i]] = rowMeans(s@data[,i:(i+1)], na.rm=TRUE)*depthT[i]/100 
  }
  for(k in 1:length(depthT)){
    if(depthT[k]==100){
      s$SOCS = rowSums(as.data.frame(x[1:3]), na.rm=TRUE)
    }
    if(depthT[k]==70){
      s$SOCS = rowSums(as.data.frame(x[1:2]), na.rm=TRUE)
    }
    if(depthT[k]==30){
      s$SOCS = rowSums(as.data.frame(x[1]), na.rm=TRUE)
    }
    ## tones / ha
    writeGDAL(s["SOCS"], out.tif[k], type="Int16", mvFlag=-32767, options="COMPRESS=DEFLATE")
  }
}


