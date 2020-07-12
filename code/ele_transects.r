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

# make linear transects for all rectangles
source('~/repos/mountain_bird_diversification/code/source.r')

lss = lapply(r, make_transects, lowlands = 1000, buffer_width = 0.05)

save(lss, file = '~/data/mountain_bird_diversification/data/rect_transects.rda')



plot(ri)

lapply(lss, lines)