# R source to generate elevational transects
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

make_transects = function(ri, lowlands = 1000, buffer_width = 0.05, min_transect_length = 500) {

  lns = list()
  i   = 0

  # make line
  l = array(NA_real_, dim = c(2,2))

  while(TRUE) {

    # make counter
    i = i + 1

    # set highest point
    mxv = maxValue(ri)

    if (mxv <= (lowlands + min_transect_length)) break

    l[1,] = xyFromCell(ri, cell = which(ri[] == mxv))[1,]
    dfp = distanceFromPoints(ri, l[1,])
    dfp[ri >= lowlands] = NA
    dfp[is.na(ri)] = NA

    # set closest lowland point
    wmncell = which.min(dfp)
    if (is.na(wmncell)) break
    l[2,] = xyFromCell(ri, cell = wmncell)[1,]

    mnv = ri[wmncell]

    # save line
    if ((mxv - mnv) > min_transect_length) {
      lns[[i]] = l
      cat('added transect (', mxv,'=>', mnv,') \n')
    }

    # remove buffer around point
    rb = buffer(SpatialLines(list(Lines(Line(l), ID = 1))), 
      width = buffer_width)
    exn = extract(ri, rb, cellnumbers = TRUE)[[1]][,1]
    ri[exn] = NA

    cat('successfully eval', i, 'transect(s) \n')
  }

  return(lns)
}




#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Line intersection functions
# if point r is on segment p->q
on_segment = function(px, py, qx, qy, rx, ry) {
  if (qx <= max(px, rx) && qx >= min(px, rx) && 
      qy <= max(py, ry) && qy >= min(py, ry)) {
    return(TRUE)  
  } else {
    return(FALSE)  
  }
}

on_segmentV = function(px, py, qx, qy, rx, ry) {
  qx <= ifelse(px > rx, px, rx) &
  qx >= ifelse(px < rx, px, rx) &
  qy <= ifelse(py > ry, py, ry) &
  qy >= ifelse(py < ry, py, ry)
}


# orientation of p->q
orientation = function(px, py, qx, qy, rx, ry) {
  val = (qy - py) * (rx - qx) - (qx - px) * (ry - qy)
    if (val == 0) return(0L) else {
      if (val > 0) return (1L) else return(2L)
    }
}

orientationV = function(px, py, qx, qy, rx, ry) {
  val = (qy - py) * (rx - qx) - (qx - px) * (ry - qy)
  ifelse(val == 0, 0, ifelse(val > 0, 1, 2))
}



intersectV = function(l1, l2) {
  p1x = l1[1] 
  p1y = l1[2]
  q1x = l1[3]
  q1y = l1[4]

  p2x = l2[,1]
  p2y = l2[,2]
  q2x = l2[,3]
  q2y = l2[,4]

  o1 = orientationV(p1x, p1y, q1x, q1y, p2x, p2y)
  o2 = orientationV(p1x, p1y, q1x, q1y, q2x, q2y)
  o3 = orientationV(p2x, p2y, q2x, q2y, p1x, p1y)
  o4 = orientationV(p2x, p2y, q2x, q2y, q1x, q1y)

  b = logical(nrow(l2))
  b = b + (o1 != o2 & o3 != o4)

  #Special Cases 
  # p1, q1 and p2 are colinear and p2 lies on segment p1q1 
  b = b + (o1 == 0 & on_segmentV(p1x, p1y, p2x, p2y, q1x, q1y))

  # p1, q1 and q2 are colinear and q2 lies on segment p1q1 
  b = b + (o2 == 0 & on_segmentV(p1x, p1y, q2x, q2y, q1x, q1y))

  # p2, q2 and p1 are colinear and p1 lies on segment p2q2 
  b = b + (o3 == 0 & on_segmentV(p2x, p2y, p1x, p1y, q2x, q2y))

  # p2, q2 and q1 are colinear and q1 lies on segment p2q2 
  b = b + (o4 == 0 & on_segmentV(p2x, p2y, q1x, q1y, q2x, q2y))

  return(b > 0)
}

















# returns TRUE if line l1 intersects with line l2, FALSE otherwise
intersect = function(l1, l2) {
  p1x = l1[1,1] 
  p1y = l1[1,2]
  q1x = l1[2,1]
  q1y = l1[2,2]
  p2x = l2[1,1]
  p2y = l2[1,2]
  q2x = l2[2,1]
  q2y = l2[2,2]

  o1 = orientation(p1x, p1y, q1x, q1y, p2x, p2y)
  o2 = orientation(p1x, p1y, q1x, q1y, q2x, q2y)
  o3 = orientation(p2x, p2y, q2x, q2y, p1x, p1y)
  o4 = orientation(p2x, p2y, q2x, q2y, q1x, q1y)

  if (o1 != o2 && o3 != o4) return(TRUE) 

  #Special Cases 
  # p1, q1 and p2 are colinear and p2 lies on segment p1q1 
  if (o1 == 0 && on_segment(p1x, p1y, p2x, p2y, q1x, q1y)) return(TRUE) 

  # p1, q1 and q2 are colinear and q2 lies on segment p1q1 
  if (o2 == 0 && on_segment(p1x, p1y, q2x, q2y, q1x, q1y)) return(TRUE) 

  # p2, q2 and p1 are colinear and p1 lies on segment p2q2 
  if (o3 == 0 && on_segment(p2x, p2y, p1x, p1y, q2x, q2y)) return(TRUE) 

  # p2, q2 and q1 are colinear and q1 lies on segment p2q2 
  if (o4 == 0 && on_segment(p2x, p2y, q1x, q1y, q2x, q2y)) return(TRUE) 

  return (FALSE)
}


#vectorized function to calculate overlap
# t1mn = vector of minimum monthly temperature for site 1
# t1mx = vector of maximum monthly temperature for site 1
# t2mn = vector of minimum monthly temperature for site 2
# t2mx = vector of maximum monthly temperature for site 2
toverlap = function(t1mn, t1mx, t2mn, t2mx) {

  r1 = t1mx - t1mn
  r2 = t2mx - t2mn

  max1 = pmax(t1mx, t1mn)
  min1 = pmin(t1mx, t1mn)
  max2 = pmax(t2mx, t2mn)
  min2 = pmin(t2mx, t2mn)

  ov1 = (pmin(max1,max2) - pmax(min1,min2))/r1
  ov2 = (pmin(max1,max2) - pmax(min1,min2))/r2
  wz = max1 < min2 | max2 < min1
  ov1[wz] = ov2[wz] = 0.0

  return(sum((ov1+ov2)/2))
}


# estimate toverlap matrix for a data frame with 24 columns
# arranged as minimum monthly temperature and then maximum monthly temperature
toverlapM = function(dt){

  ovM = array(NA_real_, dim = c(nrow(dt), nrow(dt)))

  for (j in 1:nrow(dt)) {
    for (i in 1:j) {
      if (i == j) {
        ovM[j,i] = 12.0
      } else {
        ovM[j,i] = toverlap(dt[i,1:12], dt[i,13:24], dt[j,1:12], dt[j,13:24])
      }
    }
  }

  return(ovM)
}



