# R script to generate elevational transects
#
# Ignacio Quintero
#
# t(-_-t)
#
# 12 07 2020
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

library(raster)
library(sp)
library(rgeos)
library(data.table)

# load DEM
r0 = raster('~/data/mountain_bird_diversification/data/mn30_grd/mn30_grd/w001001.adf')

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
r = lapply(e, function(x) crop(r0,x))

# remove non-mountains
mshp = readShapePoly('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')

for (i in 1:length(r)) {
  exn = extract(r[[i]], mshp, cellnumbers = TRUE)
  exn = exn[!sapply(exn, is.null)]
  exn = sapply(exn, '[', , 1)
  exn = na.omit(unlist(exn))
  r[[i]][setdiff(1:ncell(r[[i]]), exn)] = NA
  cat(i, '\n')
}


# make linear transects for all rectangles
source('~/repos/mountain_bird_diversification/code/source.r')

lss = lapply(r, make_transects, 
  lowlands = 800, buffer_width = 0.05, min_transect_length = 500)

save(r, lss, file = '~/data/mountain_bird_diversification/data/rect_transects.rda')

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
par(mfrow = c(2,4))
for (j in 1:length(r)) {
  plot(r[[j]])
  lsi = lssv[[j]]
  segments(lsi[,1], lsi[,2], lsi[,3], lsi[,4])
}


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






