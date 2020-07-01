#June 2009
#Juan Parra - juanluisparra@gmail.com
#script to generate temperature overlap and elevation difference for all pw combnations in a defined mountain

#load env data for all points

#calculate t range for each month

#calculate t overlap for any pw combination

#the following has to be done by mountain system

##########################################################
#########        Toverlap Function          ##############
##########################################################
#This function works to calculate Temperature overlap sensu 

#function to calculate overlap
toverlap <- function(data1,data2) {
range1   <- matrix(NA,nrow=12,ncol=1)
range2   <- matrix(NA,nrow=12,ncol=1)
overlap1 <- matrix(NA,nrow=12,ncol=1)
overlap2 <- matrix(NA,nrow=12,ncol=1)
overlap3 <- matrix(NA,nrow=12,ncol=1)
firstTmax = 15  #number of column before tmax jan. assumes that tmin and tmax are arranged sequentially for the 12 months.
firstTmin =3 #number of column before tmin jan. assumes that tmin and tmin are arranged sequentially for the 12 months.  
for (i in 1:12) {
	range1[i,1] <- data1[1,firstTmax + i] - data1[1,firstTmin + i]
	range2[i,1] <- data2[1,firstTmax + i] - data2[1,firstTmin + i]
	max1 <- max(c(data1[1,firstTmax + i], data1[1,firstTmin + i]))
	min1 <- min(c(data1[1,firstTmax + i], data1[1,firstTmin + i]))
	max2 <- max(c(data2[1,firstTmax + i], data2[1,firstTmin + i]))
	min2 <- min(c(data2[1,firstTmax + i], data2[1,firstTmin + i]))
		if ( max1 < min2 | max2 < min1) {
		overlap1[i,1] <- 0
		overlap2[i,1] <- 0
		}
		else {
		overlap1[i,1] <- (min(c(max1, max2)) - max(c(min1,min2)))/range1[i,1]
		overlap2[i,1] <- (min(c(max1, max2)) - max(c(min1,min2)))/range2[i,1]
		}		
	}
overlap4 <- sum((overlap1 + overlap2)/2)
return(overlap4)
}

##################################################################################
#########        Apply Toverlap function to mountain system         ##############
##################################################################################
#test should be a list, each element containing the following table in the following format:
#x (column with longitude), y(column with latitude), elev (column with elevation), tmin01, ... ,tmin12, tmax01, ... ,tmax12 (values of tmin and tmax for each month)
#each element of the list can represent a different mountain system

#out would be a list where each element contains the results for each mountain system

out <- vector("list", 2) 

for (o in 1:1) {
S <- dim(test[[o]]) 
out[[o]] <- matrix(NA, nrow=choose(S[1],2),ncol=9)
colnames(out[[o]]) <- c("mountainID","x1","y1","x2","y2", "alt1","alt2", "elevdiff","Toverlap")
u = 1                                                     
for (i in 1:(S[1]-1)) {                                         
	for (a in 1:i) {   
	out[[o]][u,1] <- o                    #id for mountain mountain
	
	out[[o]][u,2] <- test[[o]][i+1,1]     #x1 Careful that all of these ones match the actual column that they refer to
	out[[o]][u,3] <- test[[o]][i+1,2]     #y1 
	out[[o]][u,4] <- test[[o]][a,1]       #x2
	out[[o]][u,5] <- test[[o]][a,2]       #y2
	
	out[[o]][u,6] <- test[[o]][i+1,3]     #alt1
	out[[o]][u,7] <- test[[o]][a,3]       #alt2
	
	out[[o]][u,8] <- abs(test[[o]][i+1,3] - test[[o]][a,3]) #elevdiff
	
	out[[o]][u,9] <- toverlap(test[[o]][i+1,],test[[o]][a,])#toverlap 

	u = u+1                                                 
}}
}
