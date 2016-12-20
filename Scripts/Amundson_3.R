## Supplementary file 1 of Amundson et al. AUK

# C. Amundson camundson@usgs.gov 10_17_13

# Generating autocorrelated point count data with covariates for use with combined distance-sampling and time-
# removal BUGS model

# Will need the following packages: boot, gstat, lattice, latticeExtra, plyr, R2jags, unmarked.

# Data generated are based on the following
#  An imaginary bird 
# Grassland obligate, likes being near wetlands, dislikes agriculture, cares not about forest
# Lives in a world consisting of 4 habitat types: grass, wetland, ag, and trees
# Observers went out and conducted 300-m radius point-transect counts
# Observers measured proportion of each habitat within 300-m of the point

# 1.  Each covariate and transect effect are spatially autocorrelated - meaning they're not distributed
# randomly in the environment.  For the transect effect that means that birds are not evenly distributed
# independent of covariates - intrinsic autocorrelation
# 2.  In addition, habitat covariates are moderately correlated with each other (e.g., greater wetland area with greater grassland area)
# Furthermore, perceptibility is related to tree cover, and availability is related to date within season
# This appendix details how the data should be simulated for entry into JAGS (R2jags)

########################## START ###################################################################################

#################################### SOME INPUT FOR DATA SIMULATION ################################################
nsites<-100 # can change amount of replication
J<-3 # time periods
maxd<- 300 # maximum radial distance.  We can scale (e.g., maxd/10) to facilitate convergence-results are the same
# as long as you also scale breaks, delta, and interpret sigma on the same scale.
nbreaks<-5 # number of distance bins
breaks <- seq(0,maxd,by=maxd/nbreaks) # distance bins NOTE THAT THESE ARE OUTER BOUNDS OF EACH DISTANCE BIN AND NOT MIDPOINTS NEEDED FOR BUGS MODEL



####################################################################################################################
######## GENERATE HABITAT COVARIATES ###############################################################################
# generate a spatially autocorrelated landscape
# this takes a little processing time - I set the seed here so repeat simulations use the same landscape.
ncells = 10000 # number of cells in 100 x 100 grid (AKA: landscape)
set.seed(9116) # this will set the random seed - undo if replicating simulated datasets. Set here so output provided makes sense.
# fake, uncorrelated observations
X1 <- rnorm(ncells)
X2 <- rnorm(ncells)
X3 <- rnorm(ncells)
X4 <- rnorm(ncells)

###############################################
# Fake sigma... correlated decreases distance.
# Spatial autocorrelation WITHIN each covariate - clustered distributions of each habitat type in the landscape
sigma <- diag(ncells)
# Assign correlations for 'level' of spatial clustering within each habitat (e.g., wetlands are often patchily distributed in the landscape)
wet_corr<-0.5 # wetland distribution more clustered
ag_corr<-0.1 # agriculture is pretty uniform
grass_corr<-0.3
trees_corr<-0.4

sigma_wet <- wet_corr ^ abs(row(sigma)-col(sigma))
sigma_ag <- ag_corr ^ abs(row(sigma)-col(sigma))
sigma_grass <- grass_corr ^ abs(row(sigma)-col(sigma))
sigma_trees <- trees_corr ^ abs(row(sigma)-col(sigma))

###############################################
# Y is autocorrelated...centered on 0
Y_wet <- t(X1 %*% chol(sigma_wet))
Y_ag <- t(X2 %*% chol(sigma_ag))
Y_grass <- t(X3 %*% chol(sigma_grass))
Y_trees <- t(X4 %*% chol(sigma_trees))


plot(Y_wet)
plot(Y_ag)
plot(Y_grass)
plot(Y_trees)


# Skew in common and rare habitats, which is realistic.
hist(Y_ag)
hist(Y_wet)
hist(Y_grass)
hist(Y_trees)

# Setting up correlation among covariates: order = ag wet grass tree 
"[,1] [,2] [,3] [,4]
[1,]  1.0 -0.4  0.6 -0.1
[2,] -0.4  1.0 -0.2  0.1
[3,]  0.6 -0.2  1.0  0.2
[4,] -0.1  0.1  0.2  1.0"

# Correlation matrix - pick your poison - how are habitats correlated with each other? (e.g., more wetlands in grassland than farmland)
R<-matrix(cbind(1,-0.4,0.6,-0.1, -0.4,1,-0.2,0.1,0.6,-0.2,1,0.2,-0.1,0.1,0.2,1), nrow=4)
R # Make sure symmetric and positive definite or chol function won't work
U<-t(chol(R))
nvars<-dim(U)[1]

tran_vars<-t(cbind(Y_wet,Y_ag,Y_grass,Y_trees))
random_corr<-matrix(tran_vars, nrow=nvars)
X<-U%*% random_corr
newX<-t(X)
raw<-as.data.frame(newX)
orig.raw<-as.data.frame(t(random_corr))
names(raw)<-c("wet","ag","grass","tree")
cor(raw)
plot(head(raw,100))
plot(head(orig.raw,100))

# Unstandardizing the covariates to obtain real values with desired proportions. Lots of farmland, few wetlands.
wet<-raw$wet*0.2+0.15 # mean 15% wetland with SD=0.2
wet<-ifelse(wet<0,0,wet) # limit to > 0
summary(wet) # creates reasonable distribution
plot(wet)
ag<-raw$ag*0.3+0.7 # proportion of each site that's ag - here 70% with SD=0.3
ag<-ifelse(ag<0,0,ifelse(ag>1,1,ag)) # limit to between 0 - 1
summary(ag) # decent distribution
grass<-raw$grass*0.2+0.1 # landscape is 10% grass with SD=0.2
grass<-ifelse(grass<0,0,grass) # limit to > 0
summary(grass)
tree<-raw$tree*0.2+0.05 # here should be about 5% trees
tree<-ifelse(tree<0,0,tree)
summary(tree) # again, reasonable
plot(tree)

# Now scale so sum at each site = 1 (i.e., it makes the whole universe at each site - no more, no less)
comp<-rowSums(cbind(wet,ag,grass,tree))
hist(comp) # note a few sites maybe zero - these will produce NAs that need to be assigned values
newwet<-wet/comp
newag<-ag/comp
newgrass<-grass/comp
newtree<-tree/comp
# Now assign NAs a reasonable habitat composition.  For simplicity - we assume all NA sites are 100% ag
newwet[is.na(newwet)]=0
newag[is.na(newag)]=1
newgrass[is.na(newgrass)]=0
newtree[is.na(newtree)]=0
newcomp<-rowSums(cbind(newwet,newag,newgrass,newtree))
summary(newcomp)
cor(cbind(newwet,newag,newgrass,newtree)) # actual correlations among habitats in the landscape
# new correlation - final
"          newwet      newag    newgrass     newtree
newwet    1.0000000 -0.7709536  0.49846384 -0.20455434
newag    -0.7709536  1.0000000 -0.77564139 -0.34457040
newgrass  0.4984638 -0.7756414  1.00000000  0.02803595
newtree  -0.2045543 -0.3445704  0.02803595  1.00000000"

# plot landscape
library(gstat)
# create structure
xy <- expand.grid(1:100, 1:100)
grid<-xy[,1:2]

# Spatial.plot function to map landscape - created by M. Kery, Swiss Ornithological Institute, 4/2012
#####################################################################################################################
"spatial.plot" <-
  function(x,y){
    nc<-as.numeric(cut(y,20))
    plot(x,pch=" ")
    points(x,pch=20,col=topo.colors(20)[nc],cex=1)
    image.scale(y,col=topo.colors(20))
    
  }

"image.scale" <-
  function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks",
                                                                     "ranges"))
  {
    # sort out the location
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    if (missing(x))
      x <- mx + 1.05*dx/2  # default x to right of image
    else if (is.list(x)) {
      if (length(x$x) == 2)
        size <- c(diff(x$x), -diff(x$y)/n)
      y <- x$y[1]
      x <- x$x[1]
    } else x <- x[1]
    if (is.null(size))
      if (is.null(y)) {
        size <- 0.5*0.618*dy/n  # default size, golden ratio
        y <- my + 0.618*dy/2	# default y to give centred scale
      } else size <- (y-my)*2/n
    if (length(size)==1)
      size <- rep(size, 2)	# default square boxes
    if (is.null(y))
      y <- my + n*size[2]/2
    # draw the image scale
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2],
         col = rev(col), xpd = TRUE)
    # sort out the labels
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format="f", digits=digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
      ypts <- y - c(0, i) * size[2]
    else {
      bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
      ypts <- y - (i - 0.5) * size[2]
    }
  }
### END FUNCTION #########################################################################################################

# Now map proportion of each habitat in the landscape for illustration
spatial.plot(grid,newwet)
spatial.plot(grid,newag)
spatial.plot(grid,newgrass)
spatial.plot(grid,newtree)



# Now sample nsites (cells) from landscape and scale for analysis

set.seed(9116)
# Scale in second step here to facilitate plots of results
wet1<-sample(newwet,nsites,replace=FALSE)
wet<-as.numeric(scale(wet1))
ag1<-sample(newag,nsites,replace=FALSE)
ag<-as.numeric(scale(ag1))
grass1<-sample(newgrass,nsites,replace=FALSE)
grass<-as.numeric(scale(grass1))
tree1<-sample(newtree,nsites,replace=FALSE)
tree<-as.numeric(scale(tree1))

# convert proportions of square cells to circles because birds surveyed in a circle around the observer.
wet<-wet*(1/.78539)
ag<-ag*(1/.78539)
grass<-grass*(1/.78539)
tree<-tree*(1/.78539)


######################################################################################################################
# FOR JAGS ANALYSIS - SIMULATE POPULATION SAMPLE

# Associate covariates with lambda per survey point/site
# Coefficients
b1_trees<--0.05 # doesn't care about trees
b2_ag<--0.5 # hates ag
b3_grass<-1 # loves grass
b4_wet<-0.5 # likes wetland area
b5_wet2<--0.5 # but not too many - needs grass too

alpha<-5 # baseline abundance
#create transect specific relationships between counts - 100 transects, 10 pts each
# keep transect constant - like betas
tranie<-double() # create vector of transects - 10 survey points each
for(i in 1:10){
  trani=rep(runif(1,0.1,0.15),10) # this is trial and error -
  #put some amount of error that would create somewhat realistic results, but not outliers
  tranie = rbind(tranie, data.frame(trani))
}
beta0<-tranie$trani*alpha # transect-specific intercepts

# Abundance per transect or point
lambda1 <- exp(rowSums(cbind(beta0,b1_trees*tree,b2_ag*ag,b3_grass*grass,b4_wet*wet,b5_wet2*wet*wet))) 

summary(lambda1) # make sure simulated counts are reasonable before proceeding - no outliers or unrealistically high numbers
lambda<- lambda1*(1/.78539) # adjust for sim on square - now applies to circular count area
plot(lambda) # outliers?

# Population size per site
pop<-rpois(nsites,lambda)
summary(pop) # make sure simulated counts are reasonable before proceeding
sum(pop) # total population size of surveyed area

# Plot abundance by transect
transie<-as.factor(beta0)
plot(transie,pop,xaxt="n",ylab="Abundance",xlab="Transect") # how does lambda vary now by transect? Outliers?


#######################################

# Create site-level availability that varies by covariate  - season date
date<-rnorm(nsites,0,1)
b.a<--0.3 # availability declines throughout the season assuming surveys started during peak availability
library(boot)
p.a<-logit(0.9) #base availability
a.mean<-exp(p.a+b.a*date)/(1+exp(p.a+b.a*date)) # logit - site-level availability
mean(a.mean) # overall realized availability
mean.a<-1-(1-a.mean)^(1/J) # solve for meanpa from eq. p0=1-((1-meanpa)^J) Zippin 1958  
# meanpa is the average per time interval probability of availability

#######################################

# Create site-level perceptibility that varies by covariate - tree cover
sigma0 <- 137# For half-normal shape parameter PLAY WITH THIS UNTIL YOU REACH APPROXIMATE pdet YOU WANT

sigma<-exp(log(sigma0)+-0.3*tree) # here a negative association between tree cover and detection probability

summary(sigma) # This number used to calculate pdet - trial and error - play with this and covariate effects until you get desired pdet
# Given sigma0 of 137, maxd = 300 and covariate effects, true.pdet ~ 0.4

######### To help decipher sigma - this will backtransform sigma to pdet
library(unmarked)
ea <- 2*pi * integrate(grhn, 0,maxd, sigma=mean(sigma))$value # effective area
sqrt(ea / pi) # effective radius 
(true.pdet<-ea / (pi*maxd^2) ) # here is the perceptibility out to maxd

########################################################################################################################

## Now create observation level dataset that can be sampled
# NEED TO RESTRICT COUNT TO INDIVIDUALS IN POP THEN LATER ADD ZEROS
z<-data.frame(count=pop,siteid=1:nsites,sigma=sigma,mean.a=mean.a)
zz<-subset(z,count>0) # only individuals in population-sites with zero counts are removed temporarily
str(zz)
obsonly<-data.frame()
for (i in 1:length(zz$count)){ # basically make 1s for every individual at a site, then link to site id
  tmp=zz[i,]
  ind=rep(1,tmp$count)
  siteid=rep(tmp$siteid,tmp$count)
  sigma1=rep(tmp$sigma,tmp$count)
  mean.a1=rep(tmp$mean.a,tmp$count)
  obsonly=rbind(obsonly,data.frame(ind,siteid,sigma=sigma1,mean.a=mean.a1))
}

Z<-obsonly$ind  # create vector of only individuals in population (no zero counts)

# Create distances for each observation - translate them to pdet (perceptibility)
ncut<-sum(Z) # total number of observations
delta<-maxd/(ncut+1)
dclass<-seq(delta,maxd,,ncut)
probs<-  2*dclass/(maxd*maxd)
probs<-probs/sum(probs)
# Distances
r<- sample(dclass,length(Z),replace=TRUE,prob=probs)
angle<-runif(length(Z),0,360)
u<-r*cos(angle)
v<-r*sin(angle)
d<-sqrt(u^2+v^2) # random distances for every individual in the population
plot(u,v) # distances in a circle
hist(d) # histogram of distances - should have more obs in further distances because area is larger
# Creates distance classes; starts at 1 so no values are zero
d3<-cut(d,breaks,labels=1:nbreaks)  # break distance values into bins

pdet<- ifelse(d< maxd,1,0)*exp(-d*d/(2*(obsonly$sigma^2))) # this is perceptibility per site
mean(pdet) # average simulated pd
#- will be a little different than sampled pd

#################### START SAMPLING INDIVIDUALS FROM POPULATION #######################################  

# First, draw individuals for each time period with average availability 
# Then remove individuals captured in previous time interval

t1<-rbinom(Z,1,obsonly$mean.a) # capture some in time period 1
Z[t1==1]<-0 # remove the captured individuals from consideration for future time periods
t2<-rbinom(Z,1,obsonly$mean.a*Z)
Z[t2==1]<-0
t3<-rbinom(Z,1,obsonly$mean.a*Z)
tt<-data.frame(t1,t2,t3)
tint<-ifelse(tt$t1==1,1,ifelse(tt$t2==1,2,ifelse(tt$t3==1,3,0))) # create vector of time intervals for obs
Z1<-rowSums(tt)
mavail<-colSums(tt)/sum(pop) # remember Z is not = to count anymore, but rather unavailable individuals
sum(mavail) # simulated availability

newp<-pdet*Z1 # now figure out how many we observed (i.e., detected, pd) based on who was available

y<-rbinom(obsonly$ind,1,newp)  # generate count of birds based on combined probability of detection
# but remember here y is still observation-level

data<-NULL # create null dataset to insert values

## Create data frame with all our info
data<-data.frame(y,u,v,d=as.numeric(d3),tint,sigma=obsonly$sigma,siteid=obsonly$siteid)

# get rid of distance NAs
data$d[is.na(data$d)]=0

## Some examination of the dataset we created
summary(data)
data[1:10,]
plot(data[,c("u","v")])

# Now add in the zero counts from sites where no birds observed - see creation of 'z' above for reference
fulldata<-merge(z,data,all=TRUE) # observations plus sites with zero counts that are NA
x<-subset(data,y==1) # break it down to observations again, but with additional info from other columns

############################ FORMAT DATA FOR JAGS/BUGS ENTRY ##########################################################

# Subset to individuals with greater than zero count again for JAGS data entry
fulldata[is.na(fulldata)]=0 # change NAs to zeros; count doesn't make sense in this data frame - need to summarize by site into new.count
# get observed counts per site
library(plyr)
new.count<-ddply(fulldata,~siteid,summarize,y=sum(y),count=median(count),mtint=mean(tint),siteid=max(siteid))
str(new.count) # make sure it has R rows and sum of count = pop
sum(new.count$count) # must equal A - pop size
sum(new.count$y) # total number of individuals observed

y<-new.count$y # this is now the number of birds seen per site/point
y # make sure it looks appropriate

# Create transects
ntrans=length(unique(tranie$trani))
transect<-double() # create vector of transects - 10 survey points each
for(i in 1:ntrans){
  trans=rep(i,10)
  transect = rbind(transect, data.frame(trans))
}

tinterval<-x$tint # time interval for each observation ONLY (not sites with zero counts obviously)
nobs<-length(tinterval) # total count of individuals observed
nsites<-length(y) # number of sites
surveyid<-x$siteid # survey point/site ID only for individuals observed
delta<- maxd/nbreaks # bin size or width
mdpts<-seq(delta/2,maxd,delta) # generates midpoint distance of bins up to max distance
dclass<-x$d # distance class for each observation

########## DATA CREATION COMPLETE ####################################################################################

# To run combined distance-sampling and time-removal model from Amundson et al. AUK in JAGS
# See Supplemental file 2 for BUGS/JAGS code to create sim_model.txt file

library(R2jags)

# Bundle data

win.data<-list(y=y,surveyid=as.numeric(surveyid),dclass=as.numeric(dclass),nsurveys=nsites,nobs=nobs,delta=delta,nbreaks=nbreaks,
               mdpts=mdpts,maxd=maxd,J=max(unique(tinterval)),tinterval=as.numeric(tinterval),tree=tree,grass=grass,
               ag=ag,wet=wet,ntrans=ntrans,tran=transect$trans,date=date)

# create initial values for N and navail that are very close to actual values or model will not run!
Nst<-win.data$y+1

# Inits function
inits <- function(){list(N=Nst,navail=Nst,sigma.0=runif(1,120,150), beta0=runif(ntrans,-1,1),beta.a1=runif(1,-1,1),
                         beta1=runif(1,-1,1),beta2=runif(1,-1,1),mu.tran=runif(1,-1,1),sd.tran=runif(1,0,2),
                         beta3=runif(1,-1,1),beta4=runif(1,-1,1),beta5=runif(1,-1,1),beta.p1=runif(1,-1,1),
                         beta.a0=runif(1,-1,1))} #

# parameters to estimate
# careful printing pavail, pdet and N - will have nsites (e.g., in this example 100) values
params<-c("meansig","meanpdet","meanpavail","beta.a0","beta.a1","sigma.0","beta.p1","beta1","beta2","beta3","beta4","beta5",
          "meanN","mu.tran","sd.tran","totN","bayesp.pa","bayesp.pd","beta0","beta.tran","N")

# MCMC settings
# pavail can be subject to poor mixing in field data - keep thin high, burn-in long, and conduct sufficient number of iterations
nc<-3
ni<-15000
nb<-7500
nt<-50

## ONLY WORKS IN JAGS
# A bug fix for JAGS - model may produce error without this fix
set.factory("bugs::Conjugate", FALSE, type="sampler")

sim.fit<-jags(data=win.data, inits=inits, parameters.to.save=params, model.file="sim_model.txt",
              n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

# Examine results
print(sim.fit,dig=3)

# Create plot of loess trend of abundance vs. wetland with CIs and predicted data points

sim_att<-attach.jags(sim.fit) # attach jags model output
sim_N<-colMeans(sim_att$N) # mean abundance per site
sim_N.CI<-t(apply(sim_att$N,2,quantile,probs=c(0.025,0.975))) # create 95% CIs for each site/point


library(lattice)
library(latticeExtra)

# Functions for creating loess trends and plotting data.
panel.smoother<-function(x,y){
  panel.loess(x,y,col="black",lwd=2,span=1)
}
panel.smoother2<-function(x,y){
  panel.loess(x,y,col="black",lwd=1,lty=2,span=1)
}
panel.data<-function(x,y){
  panel.xyplot(x,y,col="darkgray")
}

xyplot(sim_N~grass1,panel=panel.smoother,xlab=list("Proportion grass",cex=2), 
       ylab=list("Predicted abundance per site",cex=2),scales=list(cex=2,
                                                                   col="black"))+ as.layer(xyplot(sim_N.CI[,1]~grass1,panel=panel.smoother2)) +
  as.layer(xyplot(sim_N.CI[,2]~grass1,panel=panel.smoother2))

# End code.
