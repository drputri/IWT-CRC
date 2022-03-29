#################################
# ENRICHED ENSEMBLE: MICROARRAY #
# Dea Putri, Davit Sargsyan     #
#################################

###--------------- Load libraries ---------------###
# Install qvalue package from bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("qvalue")
BiocManager::install("mixOmics")

require(data.table)
require(ggplot2)
require(parallel)
require(glmnet)
require(mixOmics)






###--------------- Load data ---------------###
load("Data/rna_crcad.Rda")
load("Data/group_crcad.Rda")

rna <- data.table(id=rownames(rna.crcad),
                  rna.crcad)

tmp <- data.table(id = names(group),
                  disease = group)

rna <- merge(tmp,
             rna,
             by="id",
             all=TRUE)
summary(rna$disease)

rna <- as.data.frame(rna)






###--------------- Enriched ensemble method ---------------###
### Step 1: univariable p-values
univar <- lapply(rna[, 3:ncol(rna)],
                 FUN = function(a) {
                   # t.test(a ~ rna$disease,
                   #        alternative = "two.sided",
                   #        var.equal = FALSE)$p.value
                   s1 <- summary(glm(rna$disease ~ a,
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




### Step 3: LASSO with cross-validation
#### Create a cluster of cores
ncores <- detectCores()
ncores

cl <- makeCluster(getOption("cl.cores",
                            4))
clusterExport(cl=cl, 
              varlist = c("rna",
                          "wgt"))

cl



#### Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:1000,
                   fun = function(i) {
                     require(data.table)
                     require(glmnet)
                     require(mixOmics)
                     
                     # sample 2/3 of the subjects at random
                     id_train <- sample(x=rna$id,
                                       size=floor(2*nrow(rna)/3),
                                       replace=FALSE)
                     
                     # sample sqrt(features), weighted
                     feat_keep <- sample(x=3:ncol(rna),
                                         size=floor(sqrt(ncol(rna) - 2)),
                                         replace=FALSE,
                                         prob=wgt)
                     
                     
                     # training data
                     dt1 <- rna[rna$id %in% id_train,
                                c(2,
                                  feat_keep)]
                     
                     # test data
                     dt2 <- rna[!rna$id %in% id_train,
                                c(2,
                                  feat_keep)]
                     
                     # glmnet
                     try({
                       # get lambda min
                       m0 <- cv.glmnet(x=as.matrix(dt1[, -1]),
                                       y=dt1$disease,
                                       family="binomial",
                                       alpha=1,
                                       nfolds=5)
                       
                       # use lambda min in glmnet
                       m1 <- glmnet(x=as.matrix(dt1[, -1]),
                                    y=dt1$disease,
                                    family="binomial",
                                    alpha=1,
                                    lambda=m0$lambda.min)
                       
                       feat_selected <- rownames(m1$beta)[as.numeric(m1$beta) != 0]
                       
                       # model evaluation
                       x.test <- model.matrix(disease ~., dt2)[,-1]
                       pred.out <- as.data.frame(predict(m1,
                                                         newx=x.test, 
                                                         s="lambda.min", 
                                                         type="class",
                                                         alpha=1))
                       colnames(pred.out) <- "pred"
                       pred.out$pred <- factor(pred.out$pred,
                                               levels=c("Adenoma", "CRC"))
                       conf <- caret::confusionMatrix(pred.out$pred, as.factor(dt2$disease))
                     })
                     
                     if (length(feat_selected) > 0) {
                       
                       
                       return(list(status="OK",
                                   feat_sampled=rownames(m1$beta),
                                   feat_sampled_all=data.table(feat=colnames(rna)[3:ncol(rna)],
                                                               feat_sampled=colnames(rna)[3:ncol(rna)] %in%
                                                                 rownames(m1$beta)),
                                   feat_selected=feat_selected,
                                   feat_selected_all=data.table(feat=colnames(rna)[3:ncol(rna)],
                                                                feat_selected=colnames(rna)[3:ncol(rna)] %in% feat_selected),
                                   dev_ratio=m1$dev.ratio,
                                   BER=get.BER(t(conf$table)),
                                   mce=1-conf$overall[["Accuracy"]],
                                   kappa=conf$overall[["Kappa"]],
                                   specificity=conf$byClass[["Specificity"]],
                                   sensitivity=conf$byClass[["Sensitivity"]],
                                   ppv=conf$byClass[["Pos Pred Value"]],
                                   npv=conf$byClass[["Neg Pred Value"]]))
                     } else {
                       return(status = NA)
                     }
                   })
})

stopCluster(cl)
gc()


# Remove empty iterations
tmp <- lapply(out,
              function(a) {
                is.na(a$status)
              })

tmp <- do.call("c", tmp)

out[tmp == TRUE] <- NULL



## Number of times each feature was sampled
# - Unweighted
tmp <- lapply(out,
              function(a){
                as.data.frame(a$feat_sampled_all) 
              })


dt_feat_sampled <- Reduce(function(x, y)  {merge(x, y, all = TRUE, by = "feat")}, tmp)

# - Weighted
tmp <- lapply(out,
              function(a){
                df1 <- as.data.frame(a$feat_sampled_all) 
                df1$feat_sampled_wgt <- df1$feat_sampled*a$dev_ratio
                df1$feat_sampled <- NULL
                return(df1)
              })


dt_feat_sampled_wgt <- Reduce(function(x, y)  {merge(x, y, all = TRUE, by = "feat")}, tmp)

# Merge
t1 <- data.table(id=dt_feat_sampled[, 1],
                 feat_sampled=rowSums(dt_feat_sampled[, -1]))

t2 <- data.table(id = dt_feat_sampled_wgt[, 1],
                 feat_sampled_wgt=rowSums(dt_feat_sampled_wgt[, -1]))

t12 <- merge(t1, t2, by = "id")


# Plot: feature by importance
res_uArray <- read_csv("Output/res_uArray.csv") 

dt <- res_uArray %>%
  filter(complete.cases(.),
         imp!=0) %>%
  mutate(name=droplevels(factor(name)))


ggplot(dt, aes(x=reorder(id, -imp), y=imp), binwidth=0) +
  geom_bar(stat="identity", fill="black", width=1, position=position_dodge()) +
  labs(x="Microarray Features",
       y="Importance") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Plot: feature by weighted importance
dt <- res_uArray %>%
  filter(complete.cases(.),
         imp_wgt!=0) %>%
  mutate(name=droplevels(factor(name)))


ggplot(dt, aes(x=reorder(id, -imp_wgt), y=imp_wgt), binwidth=0) +
  geom_bar(stat="identity", fill="black", width=1, position=position_dodge()) +
  labs(x="Microarray Features",
       y="Importance") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


