# R source to generate elevational transects
#
# Ignacio Quintero
#
# t(-_-t)
#
# 12 07 2020
#
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

make_transects = function(ri, lowlands = 1000, buffer_width = 0.05) {

  lns = list()
  i   = 0

  # make line
  l = array(NA_real_, dim = c(2,2))

  while(TRUE) {

    # make counter
    i = i + 1

    # set highest point
    mxv = maxValue(ri)

    print(mxv)

    if (mxv <= lowlands) break

    l[1,] = xyFromCell(ri, cell = which(ri[] == mxv))[1,]
    dfp = distanceFromPoints(ri, l[1,])
    dfp[ri >= lowlands] = NA

    # set closest lowland point
    l[2,] = xyFromCell(ri, cell = which.min(dfp))[1,]

    # save line
    lns[[i]] = l

    # remove buffer around point
    rb = buffer(SpatialLines(list(Lines(Line(l), ID = 1))), 
      width = buffer_width)
    exn = extract(ri, rb, cellnumbers = TRUE)[[1]][,1]
    ri[exn] = NA

    cat('successfully made', i, 'transect(s) \n')
  }

  return(lns)
}




