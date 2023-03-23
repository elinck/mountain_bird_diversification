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
library(sf)
library(data.table)



# load DEM
r0 = raster('~/data/mountain_bird_diversification/data/wc2.1_30s_elev.tif')

rasterOptions(tmpdir = "~/Desktop/tmp")

# read mountains
mshp = read_sf('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')

# crop raster for each mountain
rmt = list() # raster for only mountain
rex = list() # raster for the extent associated
for (i in 1:nrow(mshp)) {
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

for (i in seq_along(rmt)) {
  writeRaster(rmt[[i]], 
    paste0('/Users/quintero/data/mountain_bird_diversification/data/mountain_rasters/rmt_',i),
    overwrite = TRUE)
  writeRaster(rex[[i]], 
    paste0('/Users/quintero/data/mountain_bird_diversification/data/mountain_rasters/rex_',i),
  overwrite = TRUE)
}

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
library(raster)
library(sf)
library(data.table)


rasterOptions(tmpdir = "~/Desktop/tmp")

# make linear transects for all rectangles
source('~/repos/mountain_bird_diversification/code/source.r')


# make transects according to 40% of elevational range
# and lowlands defined as 10% of elevational range

lss = list()
for (i in seq_len(47)) {

  rmti = raster( 
    paste0('/Users/quintero/data/mountain_bird_diversification/data/mountain_rasters/rmt_',i))
  rexi = raster( 
    paste0('/Users/quintero/data/mountain_bird_diversification/data/mountain_rasters/rex_',i))

  mxv = maxValue(rmti)
  mnv = minValue(rmti)

  lss[[i]] = 
    make_transects(rmti, rexi,
       min_dist = 110000, buffer_width = 50000, min_transect_length = 1000)

  save(lss, file = '~/data/mountain_bird_diversification/data/mount_transects.rda')

  cat(i, '\n')
}

save(lss, file = '~/data/mountain_bird_diversification/data/mount_transects.rda')

# load('~/data/mountain_bird_diversification/data/mount_transects.rda')


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
library(sf)
library(data.table)
library(RColorBrewer)

r0 = raster('~/data/mountain_bird_diversification/data/wc2.1_30s_elev.tif')
mshp = read_sf('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')

r0 = crop(r0, extent(-180.0001, 179.999, -60.00014, 83.99986))

pdf(file = '~/data/mountain_bird_diversification/plots/mount_transects.pdf',
  height = 10, width = 22)
  par(bty = 'n')
  plot(r0, col = brewer.pal(9,"Spectral"), axes = FALSE)
  plot(mshp, add = TRUE, col = rgb(128,0,128,150, max = 255))
  for (j in 1:length(lss)) {
    lsi = lssv[[j]]
    lsi = matrix(lsi, ncol= 4)
    if (length(lsi) == 0) next
    segments(lsi[,1], lsi[,2], lsi[,3], lsi[,4], lwd = 2, col = "black")
    print(j)
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

r0 = raster('~/data/mountain_bird_diversification/data/wc2.1_30s_elev.tif')

# read mountains
mshp = read_sf('~/data/mountain_bird_diversification/data/MountSys/Global_GMBA.shp')


## extract Temperatures for each line

lssv[[40]] = as.matrix(lssv[[40]], ncol = 4)

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
el = list()
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
  ex[[j]] = extract(tr, spl)
  el[[j]] = extract(r0, spl)
  cat(j, '\n')
}


# remove NAs
for (j in 1:length(ex)) {
  exi = ex[[j]]
  eli = el[[j]]

  lna = sapply(eli, anyNA)
  if (length(lna) == 0) next

  wna = which(sapply(eli, anyNA))

  for (i in wna) {
    elj = eli[[i]]

    inc = TRUE
    if (is.na(elj[1])) {
      inc = TRUE 
    } else {
      if (is.na(elj[length(elj)])) {
        inc = FALSE
      } else {
        if (elj[1] > elj[length(elj)]) { 
            inc = FALSE
        }
      }
    }

    if (inc) {
      idx = (max(which(is.na(elj)))+1):length(elj)
    } else {
      idx = 1:(min(which(is.na(elj)))-1)
    }
    eli[[i]] = elj[idx]
    exi[[i]] = exi[[i]][idx]
  }
}




# Estimate temperature and elevation overlap/distance
to = list()
ed = list()
for (j in 1:length(ex)) {
  exi = ex[[j]]
  eli = el[[j]]
  if (length(exi) == 0) {
    to[[j]] = list()
    ed[[j]] = list()
    next
  }
  toi = list()
  edi = list()
  for (i in 1:length(exi)) {
    exii = exi[[i]]
    elii = eli[[i]]
    toi[[i]] = toverlapM(exii)
    edi[[i]] = dist(elii)
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




