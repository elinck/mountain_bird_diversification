# R script to generate mountain elevational transects
#
# Ignacio Quintero
#
# t(-_-t)
#
# 10 08 2020
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(data.table)

###
# NOW IN CLUSTER
###

# load DEM
r0 = raster('~/mn30_grd/mn30_grd/w001001.adf')

r0[r0==0] = NA

rasterOptions(tmpdir = "/mnt/data/personal/ignacioq/raster_tmp")

# read mountains
mshp = readShapePoly('~/MountSys/Global_GMBA.shp')

# crop raster for each mountain
rmt = list() # raster for only mountain
rex = list() # raster for the extent associated
for (i in 1:length(mshp)) {
  ext = extent(mshp[i,])
  ext@xmin = ext@xmin - 1
  ext@xmax = ext@xmax + 1
  ext@ymin = ext@ymin - 1
  ext@ymax = ext@ymax + 1
  rex[[i]] = rmt[[i]] = crop(r0, ext)
  exn = extract(rmt[[i]], mshp[i,], cellnumbers = TRUE)
  exn = exn[!sapply(exn, is.null)]
  exn = sapply(exn, function(x) x[,1])
  exn = na.omit(unlist(exn))
  rmt[[i]][setdiff(1:ncell(rmt[[i]]), exn)] = NA
  cat(i, '\n')
}


# save cropped rasters
save(rmt, rex, file = '~/mountain_rasters.rda')


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(data.table)
load('~/mountain_rasters.rda')

# make linear transects for all rectangles
source('~/source.r')


# make transects according to 50% of elevational range
# and lowlands defined as 10% of elevational range
mxv = sapply(rmt, maxValue)
mnv = sapply(rmt, minValue)

rgi = (mxv - mnv)*0.2

lss = list()
for (i in 1:length(rmt)) {
  lss[[i]] = 
    make_transects(rmt[[i]], rex[[i]],
       min_dist = 110000, buffer_width = 0.5, min_transect_length = rgi[i])
  cat(i, '\n')
}


save(lss, file = '~/data/mountain_bird_diversification/data/mount_transects.rda')

load('~/data/mountain_bird_diversification/data/mount_transects.rda')

# make array from list of lines
lssv = list()

for (j in 1:length(lss)) {

  lsi  = lss[[j]]
  if (length(lsi) == 0) {
    lssv[[j]] = lsi
    next
  }
  lsiv = array(NA_real_, dim = c(length(lsi), 4))

  lsiv[,1] = sapply(lsi, '[', 1, 1)
  lsiv[,2] = sapply(lsi, '[', 1, 2)
  lsiv[,3] = sapply(lsi, '[', 2, 1)
  lsiv[,4] = sapply(lsi, '[', 2, 2)

  lssv[[j]] = lsiv
}

# remove intersecting transects
for (j in 1:length(lss)) {
  lsi = lssv[[j]]
  if (length(lsi) == 0) {
    lssv[[j]] = lsi
    next
  }
  for (i1 in nrow(lsi):1) {
    l1 = lsi[i1,]
    if (anyNA(l1)) next
    if (any(intersectV(l1, lsi[setdiff(1:nrow(lsi),i1),]))) {
      lsi[i1,] = 0.0
    }
  }
  lsi = lsi[!(lsi[,1] == lsi[,3] & lsi[,2] == lsi[,4]),]
  # remove NAs
  lssv[[j]] = lsi
  cat(j, '\n')
}




library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(data.table)

r0 = raster('~/data/mountain_bird_diversification/data/mn30_grd/mn30_grd/w001001.adf')
mshp = readOGR('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')

NAvalue(r0) = 0

r0 = crop(r0, extent(-180.0001, 179.999, -60.00014, 83.99986))

pdf(file = '~/data/mountain_bird_diversification/plots/mount_transects.pdf',
  height = 10, width = 22)
  par(bty = 'n')
  plot(r0, col = brewer.pal(9,"Spectral"), axes = FALSE)
  plot(mshp, add = TRUE, col = rgb(128,0,128,150, max = 255))
  for (j in 1:length(lss)) {
    lsi = lssv[[j]]
    if (length(lsi) == 0) next
    segments(lsi[,1], lsi[,2], lsi[,3], lsi[,4], lwd = 2, col = "black")
  }
dev.off()

# BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# extract temperature
source('~/repos/mountain_bird_diversification/code/source.r')
load('~/data/mountain_bird_diversification/data/mount_transects.rda')
load('~/data/mountain_bird_diversification/data/mountain_rasters.rda')

# make stack of temperatures
fmin = list.files('~/data/mountain_bird_diversification/data/tmin/',
 full.names = TRUE)
fmax = list.files('~/data/mountain_bird_diversification/data/tmax/', 
  full.names = TRUE)
smin = stack(fmin)
smax = stack(fmax)

tr = stack(smin, smax)


# read mountains
mshp = readOGR('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')

# crop rasters
trr = list() # raster 
for (i in 1:length(mshp)) {
  ext = extent(mshp[i,])
  ext@xmin = ext@xmin - 1
  ext@xmax = ext@xmax + 1
  ext@ymin = ext@ymin - 1
  ext@ymax = ext@ymax + 1
  trr[[i]] = crop(tr, ext)
  cat(i, '\n')
}
trr[[40]] = NULL

# add elevation
rex = mapply(stack, rex, trr)

## extract Temperatures for each line

# make list from array of lines
lss = list()

for (j in 1:length(lssv)) {
  lsi  = lssv[[j]]

  if (length(lsi) == 0) {
    lss[[j]] = lsi
    next
  }

  li = list()
  for (i in 1:nrow(lsi)) {
    li[[i]] = matrix(lsi[i,], ncol = 2, nrow = 2, byrow = TRUE)
  }

  lss[[j]] = li
}

# make spatial lines and extract values
ex = list()
for (j in 1:length(lss)) {
  lsi = lss[[j]]
  if (length(lsi) == 0) {
    ex[[j]] = lsi
    next
  }
  lsi = lapply(lsi, Line)
  for (i in 1:length(lsi)) {
    lsi[[i]] = Lines(lsi[[i]], ID = i)
  }
  spl = SpatialLines(lsi)
  proj4string(spl) = proj4string(rex[[j]])
  ex[[j]] = extract(rex[[j]],spl)
  cat(j, '\n')
}

# Estimate temperature and elevation overlap/distance
to = list()
ed = list()
for (j in 1:length(ex)) {
  exi = ex[[j]]
  if (length(exi) == 0) {
    to[[j]] = list()
    ed[[j]] = list()
    next
  }
  toi = list()
  edi = list()
  for (i in 1:length(exi)) {
    dt = exi[[i]]
    toi[[i]] = toverlapM(dt[,2:25])
    edi[[i]] = dist(dt[,1])
  }
  to[[j]] = toi
  ed[[j]] = edi
  cat(j, '\n')
}

save(to, ed,
  file = '~/data/mountain_bird_diversification/results/mount_transects_toverlap.rda')

# make plots
load('~/data/mountain_bird_diversification/results/mount_transects_toverlap.rda')

mshp = mshp[-40,]

for (j in 1:length(to)) {

  toi = to[[j]]
  edi = ed[[j]]

  if (length(toi) == 0) next

  nam = mshp@data[j,'Name']
  if (nam =="US Great Basin/Sierra Nevada") {
    nam =  "US Great Basin-Sierra Nevada"
  }

  jpeg(file = paste0('~/data/mountain_bird_diversification/plots/temp_ele_', nam,'.jpeg'),
  width = 700, height = 500)

    plot(1, type = 'n', xlim = c(0,4000), ylim = c(0,12),
      xlab = 'Distance in elevation (m)', ylab = 'Temperature overlap', bty = 'n')

    for (i in 1:length(toi)) {
      d1 = to[[j]][[i]]
      d2 = ed[[j]][[i]]

      x = as.vector(d2)
      y = d1[lower.tri(d1)]

      ord = order(x)
      x = x[ord]
      y = y[ord]

      1/length(toi)

      lines(x, y, col = rgb(0,0,0,5/length(toi)), lwd = 2)
    }

  dev.off()
}




