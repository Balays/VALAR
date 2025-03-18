

# 3. Cluster Stability via pvclust to assess the statistical support for clusters:
library(pvclust)

#fit <- pvclust(countData, method.hclust="complete", method.dist="correlation", nboot=1000)
#plot(fit)
#pvrect(fit)

#fit2 <- pvclust(countData, method.hclust="complete", method.dist="euclidean", nboot=1000)
#plot(fit2)
#pvrect(fit2)

fit <- pvclust(countData, method.hclust=method.hclust, method.dist=method.dist, nboot=1000)
plot(fit)
pvrect(fit)
