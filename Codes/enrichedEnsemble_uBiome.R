#################################
# ENRICHED ENSEMBLE: MICROBIOME #
# Dea Putri, Davit Sargsyan     #
#################################

###--------------- Load libraries ---------------###
# Install qvalue package from bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")

require(data.table)
require(ggplot2)
require(parallel)
require(glmnet)





###--------------- Load data ---------------###
load("Data/otu_431_clr.Rda")
load("Data/group_crcad.Rda")

genus <- data.table(id=rownames(otu.clr),
                    otu.clr)
tmp <- data.table(id=names(group),
                  disease=group)

genus <- merge(tmp,
               genus,
               by="id",
               all=TRUE)

genus$disease <- factor(genus$disease)


# Data cleaning: remove non-zero variance
# Definition of non-zero variance: < 20% unique values
#                                : greater than 20 the ratio of most common and 2nd most common values

x <- genus$SVs9518

clean.func <- function(x) {
  out <- NULL
  n <- length(x)
  
  # Check the number of unique values
  n.u <- length(unique(x))
  chk.u <- n.u/n
  
  # Check the ratio of most & 2nd most common values
  if(n.u==1){
    out <- "removed"
  } else {
    chk.r <- sort(table(x), decreasing = TRUE)[1]/sort(table(x), decreasing = TRUE)[2]
    
    if(chk.u < .2 & chk.r > 20) {
      out <- "removed"
    } 
  }
  
  return(out)
}

rmv.genus <- lapply(genus[, 3:ncol(genus)],
                    FUN = function(a) {
                      # t.test(a ~ rna$disease,
                      #        alternative = "two.sided",
                      #        var.equal = FALSE)$p.value
                      chk <- clean.func(a)
                      return(chk)  #obtained the p-values
                    })
rmv.genus <- do.call("c", rmv.genus)

genus <- as.data.frame(genus)
genus <- genus[, !(names(genus) %in% names(rmv.genus))]





###--------------- Enriched ensemble method ---------------###
### Step 1: univariable p-values
univar <- lapply(genus[, 3:ncol(genus)],
                 FUN = function(a) {
                   # t.test(a ~ rna$disease,
                   #        alternative = "two.sided",
                   #        var.equal = FALSE)$p.value
                   s1 <- summary(glm(genus$disease ~ a,
                                     family="binomial"))
                   return(s1$coefficients[2, 4])  #obtained the p-values
                 })
univar <- do.call("c", univar)


# Problem with the margins (too large)
par(mar=c(1,1,1,1))
hist(univar, 20)
summary(univar)



### Step 2: weights
q <- qvalue::qvalue(univar)$qvalues
wgt <- -log(q)
wgt <- wgt/sum(wgt)

# Problem with the margins (too large)
par(mar=c(1,1,1,1))
hist(wgt, 20)



### Subset for features