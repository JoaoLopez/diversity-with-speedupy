Tracy Widom on CENTC distance matrix
========================================================

## Require stuff

```r
require("mclust")
```

```
## Loading required package: mclust Package 'mclust' version 4.2
```

```r
require("ape")
```

```
## Loading required package: ape
```

```r
require("RMTstat")
```

```
## Loading required package: RMTstat
```


## New Tracy Widom Plan - PCoA and then TW

## Functions needed for Tracy Widom

### Mmake function (only needed for non-distance data)

```r
Mmake <- function(data) {
    pcdat <- data
    mv <- apply(pcdat, 2, mean, na.rm = TRUE)
    # make vector of sqrt(p*(1-p))
    vv <- mv/2
    vv <- sqrt(vv * (1 - vv))  # need to use sd of data, not p(1-p)
    # redo data and vectors to take out monomophic sites
    # dat<-pcdat[,vv>0&is.na(vv)==FALSE] #don't need either but need something
    # named dat below mv<-mv[vv>0&is.na(vv)==FALSE] #we don't need?
    # vv<-vv[vv>0&is.na(vv)==FALSE] #we don't need?
    M <- t((t(pcdat) - mv)/vv)
    M[is.na(M)] <- 0
    return(M)
}
```


### TWcalc function (SM says ignore, as this is specific for SNP data)

```r
TWcalc <- function(dat, k) {
    x <- dat
    y <- x %*% t(x)
    eig <- eigen(y)
    eig <- eig$values
    # eig<-eig[1:(length(eig))] #Why do this? n<-ncol(x) # UNUSED?
    m0 = nrow(x)
    
    # I think this is TW residuals
    TWres_ <- c()
    for (i in 1:k) {
        m = (m0) - i + 1
        # m=m
        eiv <- eig[i:(m - 1)]
        # n_=n
        n_ = ((m + 1) * sum(eiv)^2)/((m - 1) * sum(eiv^2) - sum(eiv)^2)
        S = ((sqrt(n_ - 1) + sqrt(m)))/n_ * (1/sqrt(n_ - 1) + 1/sqrt(m))^(1/3)
        u_ = (sqrt(n_ - 1) + sqrt(m))^2/n_
        L1 <- (m - 1) * eiv[1]/sum(eiv)
        TW_ <- (L1 - u_)/S
        TWres_ <- c(TWres_, TW_)  # why cat these? 
    }
    
    # TWres_
    pres <- c()
    for (i in 1:length(TWres_)) {
        TW <- TWres_[i]
        dif <- abs(test[, 1] - TW)
        p <- test[dif == min(dif), 2]
        pres <- c(pres, p)  # why cat these? 
    }
    # pres
    res <- list(TWres_, pres)
    return(res)
}
```


### Do some test data

Read in some test data: JC distance matrix original 218 CENTC from genbank and Tracy-Widom table
(need to fix this so don't have to specify path)

```r
# read in table
dist_mat <- read.csv("~/src/tw_centc/test_matrix.csv", row.names = 1, header = T)
```

```
## Warning: cannot open file
## '/Users/paulbilinski/src/tw_centc/test_matrix.csv': No such file or
## directory
```

```
## Error: cannot open the connection
```

```r
test <- read.table("~/src/tw_centc/twtable.txt", header = FALSE)
```

```
## Warning: cannot open file '/Users/paulbilinski/src/tw_centc/twtable.txt':
## No such file or directory
```

```
## Error: cannot open the connection
```


Standardize matrix and do PC <- don't do this for distance matrix, only for individual observations (i.e. SNP data, phenotypes, environments, etc.)

```r
made_x <- Mmake(dist_mat)
```

```
## Error: object 'dist_mat' not found
```


PCOA from ape package. Check to make sure this is giving us eigenvalues 

```r
# pc<-prcomp(made_x) # for nondistance data
pc <- pcoa(dist_mat)  # get details on how this works
```

```
## Error: object 'dist_mat' not found
```


Skree plot (check to make sure in order)

```r
plot(pc$values[, 1])
```

```
## Error: object 'pc' not found
```


Do Tracy Widom. k is number of PCs to evaluate. Here, trying 150.
This TW is specific for SNP data from Patterson et al. 2006 (according to SM).  Need to redo using TracyWidom distribution from RMTstat package. Check details. 

```r
TW <- TWcalc(made_x, 150)
```

```
## Error: object 'made_x' not found
```

```r
km <- which(TW[[2]] > 0.05)[1]
```

```
## Error: object 'TW' not found
```


Cluster individuals based on first km significant principle coordinates

```r
pcx <- Mclust(pc$x[, 1:km], G = (km + 1))  # originally used VVI model
```

```
## Error: object 'pc' not found
```

```r
mainC <- pcx$classification
```

```
## Error: object 'pcx' not found
```


Plot stuff (if we can)

```r
col.regions = rainbow(km + 1)
```

```
## Error: object 'km' not found
```

```r
require(rgl)
```

```
## Loading required package: rgl
```

```
## Warning: there is no package called 'rgl'
```

```r
plot3d(pc$x[, 1], pc$x[, 2], pc$x[, 3], xlab = "PC1", ylab = "PC2", zlab = "PC3", 
    col = col.regions[mainC], type = "s", radius = 0.1)
```

```
## Error: could not find function "plot3d"
```

```r
plot(pc$x[, 1] ~ pc$x[, 3])
```

```
## Error: object 'pc' not found
```

```r

# the below fails w/o column names correctly dim<-1:3
# text3d(pc$x[,dim],text=rownames(M),col=cols,cex=0.6)
```


