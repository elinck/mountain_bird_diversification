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

# load DEM
r0 = raster('~/data/mountain_bird_diversification/data/mn30_grd/mn30_grd/w001001.adf')

# read mountains
mshp = readOGR('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')

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
  exn = sapply(exn, '[', , 1)
  exn = na.omit(unlist(exn))
  rmt[[i]][setdiff(1:ncell(rmt[[i]]), exn)] = NA
  cat(i, '\n')
}

# delete madeira
rmt[[40]] = NULL
rex[[40]] = NULL

# save cropped rasters
save(rmt, rex, file = '~/data/mountain_bird_diversification/data/mountain_rasters.rda')


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(raster)
library(sp)
library(rgeos)
library(rgdal)
library(data.table)
load('~/data/mountain_bird_diversification/data/mountain_rasters.rda')

# make linear transects for all rectangles
source('~/repos/mountain_bird_diversification/code/source.r')


# make transects according to 50% of elevational range
# and lowlands defined as 10% of elevational range
mxv = sapply(rmt, maxValue)
mnv = sapply(rmt, minValue)

rgi = (mxv - mnv)*0.2

lss = list()
for (i in 23:length(rmt)) {
  lss[[i]] = 
    make_transects(rmt[[i]], rex[[i]],
       min_dist = 110000, buffer_width = 0.5, min_transect_length = rgi[i])
  cat(i, '\n')
}


save(lss, file = '~/data/mountain_bird_diversification/data/mount_transects_23_46.rda')




save(r, lss, file = '~/data/mountain_bird_diversification/data/mount_transects.rda')

load('~/data/mountain_bird_diversification/data/mount_transects.rda')

# make array from list of lines
lssv = list()

for (j in 1:length(lss)) {
  lsi  = lss[[j]]
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



# plot results
for (j in 1:length(lss)) {
  j = 1
  plot(rmt[[j]])
  lsi = lssv[[j]]
  segments(lsi[,1], lsi[,2], lsi[,3], lsi[,4])
}


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# extract temperature
source('~/repos/mountain_bird_diversification/code/source.r')
load('~/data/mountain_bird_diversification/data/rect_transects.rda')


# make stack of temperatures
fmin = list.files('~/data/mountain_bird_diversification/data/tmin/',
 full.names = TRUE)
fmax = list.files('~/data/mountain_bird_diversification/data/tmax/', 
  full.names = TRUE)
smin = stack(fmin)
smax = stack(fmax)

tr = stack(smin, smax)

# load extents
e = list()
e[[1]] = extent(-83.94194,-82.54267,35.02879,35.71528)
e[[2]] = extent(-81.19186,-78.70011,37.08391,39.41059)
e[[3]] = extent(-74.8787,-69.28792,42.51517,45.42436)
e[[4]] = extent(-67.74567,-62.75549,-22.72413,-15.65367)
e[[5]] = extent(-76.00255,-68.58037,-18.61333,-14.04033)
e[[6]] = extent(-81.99395,-75.93518,-4.268851,2.284316)
e[[7]] = extent(-77.60973,-73.21778,1.448905 ,8.605978)

# crop rasters
tr = lapply(e, function(x) crop(tr,x))

# add elevation
r = mapply(stack, r, tr)

## extract Temperatures for each line

# make list from array of lines
lss = list()

for (j in 1:length(lssv)) {
  lsi  = lssv[[j]]

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
  lsi = lapply(lsi, Line)
  for (i in 1:length(lsi)) {
    lsi[[i]] = Lines(lsi[[i]], ID = i)
  }
  spl = SpatialLines(lsi)
  proj4string(spl) = proj4string(r[[j]])
  ex[[j]] = extract(r[[j]],spl)
  cat(j, '\n')
}

# Estimate temperature and elevation overlap/distance
to = list()
ed = list()
for (j in 1:length(ex)) {
  exi = ex[[j]]
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
  file = '~/data/mountain_bird_diversification/results/rect_transects_toverlap.rda')

# make plots
load('~/data/mountain_bird_diversification/results/rect_transects_toverlap.rda')


jpeg(file = '~/data/mountain_bird_diversification/plots/rect_transects_toverlap.jpeg',
  width = 1000, height = 500)
  par(mfrow = c(2,4), bty = 'n', las = 3)

  for (j in 1:length(to)) {

    toi = to[[j]]
    edi = ed[[j]]

    plot(1, type = 'n', xlim = c(0,4000), ylim = c(0,12),
      xlab = 'Distance in elevation (m)', ylab = 'Temperature overlap')

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
  }
dev.off()


