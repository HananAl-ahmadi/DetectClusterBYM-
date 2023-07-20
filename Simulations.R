### Code from scratch # Jonatan
library(spatstat)
library(geoR)
library(splancs)
library(INLA)
library(maptools)
library(spdep)
library(doParallel)
library(Rgraphviz)

#### loading RF observations ###################################################
load('sim1000.RData')

####  Generating counties and their respectives graphs #########################
set.seed(13)
centroids <- rpoispp(50)
counties <- dirichlet(centroids)
map <- as(counties, "SpatialPolygons")
nb <- poly2nb(map)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

####  Generating the shape  #############################################
window <- ellipse(a = .3,b = .1, centre = c(.6,.6), phi = pi / 6)

###################### BYM and NAIVE functions #################################
#Computing confusion matrix
ConfMatrix <- function(TrueCluster, ObsCluster){
  TP <- TrueCluster * ObsCluster #only preserves one values
  TN <- eval.im(ifelse(TrueCluster == 0 & ObsCluster == 0, 1, 0))
  FP <- eval.im(ifelse(TrueCluster == 0 & ObsCluster == 1, 1, 0))
  FN <- eval.im(ifelse(TrueCluster == 1 & ObsCluster == 0, 1, 0))
  
  #Sensitivity (true positive rate)
  TPR <- sum(TP) / (sum(TP) + sum(FN))
  #Specificity
  TNR <- sum(TN) / (sum(TN) + sum(FP))
  #Positive predictive value
  PPV <- sum(TP) / (sum(TP) + sum(FP))
  #fall-out or false positive rate
  FPR <- sum(FP) / (sum(FP) + sum(TN))
  #Positive likelihood ratio
  PLR <- TPR / FPR #They should be computed before
  return(c(TPR = TPR, TNR = TNR, PPV = PPV,  FPR = FPR, PLR = PLR))
}

BYM <- function(PointPattern, population){
  cases <- PointPattern
  areasq <- quadratcount.ppp(cases, tess = counties)
  # Add the data to the map 
  d <- as.data.frame(areasq)
  d$counties <- d$tile
  d$Y <- d$Freq
  d$E <- npoints(PointPattern) * tile.areas(counties)
  d$SIR <- d$Y/d$E
  d <- subset(d, select = c("counties", "Y", "E", "SIR")) 
  rownames(d) <- d$counties
  map <- SpatialPolygonsDataFrame(map, d, match.ID = TRUE)
  
  # BYM method performed by INLA
  map$re_u <- map$re_v <-  1:nrow(map@data)
  
  formula <- Y ~ f(re_u, model = "besag", graph = g, scale.model = TRUE) + f(re_v, model = "iid")
  
  res <- inla(formula, family = "poisson", data = map@data, E = E, 
              control.compute = list(return.marginals.predictor = TRUE),
              num.threads = "1:1")
  
  ########################## The Exceedance probabilities BYM #########################
  exc <- sapply(res$marginals.fitted.values, FUN = function(marg){1 - inla.pmarginal(q =2 , marginal = marg)})
  map$exc <- exc
  ########################## convert the exceedance probabilities to zero and one #########################
  # calculate the 75 quantile of the exceedance probabilities
  q.75 <- quantile(map$exc, probs = 0.75)
  
  # Replace this by the INLA output transformed to 0s and 1's 
  Pinla <- ifelse(map$exc < q.75, 0, 1)
  Y <- counties
  Yim <- as.im(Y)
  obsCluster <- eval.im(t(Pinla)[Yim])
  #plot(obsCluster)
  bSquare <- as.im(unit.square())
  bSquare$v <- bSquare$v - 1 #Canvas of zeros 
  bClust <- bSquare[window, drop = F]
  bClust <- bClust + 1
  trueCluster <- bSquare + bClust
  trueCluster$v[is.na(trueCluster$v)] <- 0
  return(ConfMatrix(TrueCluster = trueCluster, ObsCluster = obsCluster))
}

Naive <- function(PointPattern){
  cases <- PointPattern
  #Counts by area
  areasq <- quadratcount.ppp(cases, tess = counties)
  #Transforming cluster shape in a binary image 
  bSquare <- as.im(unit.square())
  bSquare$v <- bSquare$v - 1 #Canvas of zeros 
  bClust <- bSquare[window, drop = F]
  bClust <- bClust + 1
  trueCluster <- bSquare + bClust
  trueCluster$v[is.na(trueCluster$v)] <- 0 
  
  # Going from counts to grid
  areasim <- intensity.quadratcount(areasq, image = T)

  #Artificial clustering for example
  cutcounts <- cut(areasim, 4, labels = F) #Basic clustering technique
  cutcounts$v <- ifelse(cutcounts$v == 4, 1, 0)
  obsCluster <- cutcounts
  return(ConfMatrix(TrueCluster = trueCluster, ObsCluster = obsCluster))
}

###################### OneBigSimulation ########################################
# We have all the ingredients except for the point pattern
onesimu <- function(RF, increasedRisk = 1.5, population = exp(9.2)){
  B <- as.im(data.frame(z = RF$data, x = RF[[1]][,1], y =  RF[[1]][,2])) #The whole window
  C <- B[window, drop = F] #Only the ellipsis
  C <- C + log(increasedRisk) # Transformation of only ellipse
  C$v[is.na(C$v)] <- 0 #filling the empty space with zeros
  D <- eval.im(exp(B + C))
  
  #Generating cases according to risk given by exp(RF)
  PP <- rpoispp(lambda = (population * D)) # 5 is the baseline
  
  bymCF <- BYM(PointPattern = PP, population = population)
  naiveCF <- Naive(PointPattern = PP)
  return(c(bymCF, naiveCF))
}

####################### RepeatingSimulations  ########################################
Results <- mclapply(sim1000, onesimu, mc.cores = 45, mc.preschedule = F)
avg <- colMeans(do.call(rbind, Results))
