library(phytools)
library(geiger)
library(phylolm)
library(ggplot2)
library(forcats)
library(geomorph)
library(dplyr)
library(mvMORPH)

setwd("~/Documents/Projects/Project_Squirrels/squirrel_limbs")

load("~/Documents/Projects/Project_Squirrels/squirrel_limbs/Rdata/femur_data.RData")

#### femur data ####
femur_data <- read.csv("data/femur_combined_data.csv")
rownames(femur_data) <- femur_data[,1]
tree <- read.nexus("data/squirrel_mcctree.nex")

#remove species in phylogeny that aren't present in dataset
tree_pr <- treedata(phy = tree, data = femur_data, sort=TRUE)$phy

femur_data <- femur_data[tree_pr$tip.label,]

size <- (as.numeric(as.character(femur_data$geomean))) #save size as vector
names(size) <- rownames(femur_data) #give species names for size vector
lnsize <- log(size) #natural log transform size vector (for stats purposes)

ecotype <- as.factor(femur_data$ecotype) #save ecology as factor
names(ecotype) <- rownames(femur_data) #give species names for ecology factor

col_ecotype <- c("#073b4c", "#118ab2", "#ffd166", "#06d6a0") #set colors for ecological groups: blue = flying; orange = ground; green = tree

#internal structures data 
femur_Cg <- (as.numeric(as.character(femur_data$Cg))) #save size as vector
names(femur_Cg) <- rownames(femur_data) #give species names for size vector

femur_DE <- (as.numeric(as.character(femur_data$DE))) #save size as vector
names(femur_DE) <- rownames(femur_data) #give species names for size vector

femur_CSS <- (as.numeric(as.character(femur_data$CSS))) #save size as vector
names(femur_CSS) <- rownames(femur_data) #give species names for size vector

#external structures data
femur_rawlandmarks <- femur_data[,12:467]
femur_gpa <- gpagen(arrayspecs(femur_rawlandmarks, 152, 3))
femur_shape <- femur_gpa$coords
lnfemur_size <- log(femur_gpa$Csize)
femur_shape_2D <- two.d.array(femur_shape)
femur_shape_data <- data.frame(ecotype, lnfemur_size, femur_shape_2D)

##### data - individual ecotypes #####
#chipmunk data
femur_data_chipmunk <- filter(femur_data, ecotype == "chipmunk")
femur_Cg_chipmunk <- femur_data_chipmunk$Cg
names(femur_Cg_chipmunk) <- rownames(femur_data_chipmunk)
femur_DE_chipmunk <- femur_data_chipmunk$DE
names(femur_DE_chipmunk) <- rownames(femur_data_chipmunk)
femur_CSS_chipmunk <- femur_data_chipmunk$CSS
names(femur_CSS_chipmunk) <- rownames(femur_data_chipmunk)
femur_shape_data_chipmunk <- filter(femur_shape_data, ecotype == "chipmunk")
femur_lnsize_chipmunk <- femur_shape_data_chipmunk$lnfemur_size
names(femur_lnsize_chipmunk) <- rownames(femur_shape_data_chipmunk) 
femur_shape_chipmunk <- arrayspecs(femur_shape_data_chipmunk[,3:458], 152,3 )
femur_shape_chipmunk_2D <- two.d.array(femur_shape_chipmunk)
lnfemur_size_chipmunk <- femur_shape_data_chipmunk$lnfemur_size
tree_pr_chipmunk <- treedata(phy = tree, data = femur_shape_data_chipmunk, sort=TRUE)$phy

#gliding data
femur_data_gliding <- filter(femur_data, ecotype == "gliding")
femur_Cg_gliding <- femur_data_gliding$Cg
names(femur_Cg_gliding) <- rownames(femur_data_gliding)
femur_DE_gliding <- femur_data_gliding$DE
names(femur_DE_gliding) <- rownames(femur_data_gliding)
femur_CSS_gliding <- femur_data_gliding$CSS
names(femur_CSS_gliding) <- rownames(femur_data_gliding)
femur_shape_data_gliding <- filter(femur_shape_data, ecotype == "gliding")
femur_lnsize_gliding <- femur_shape_data_gliding$lnfemur_size
names(femur_lnsize_gliding) <- rownames(femur_shape_data_gliding) 
femur_shape_gliding <- arrayspecs(femur_shape_data_gliding[,3:458], 152,3 )
femur_shape_gliding_2D <- two.d.array(femur_shape_gliding)
lnfemur_size_gliding <- femur_shape_data_gliding$lnfemur_size
tree_pr_gliding <- treedata(phy = tree, data = femur_shape_data_gliding, sort=TRUE)$phy

#ground data
femur_data_ground <- filter(femur_data, ecotype == "ground")
femur_Cg_ground <- femur_data_ground$Cg
names(femur_Cg_ground) <- rownames(femur_data_ground)
femur_DE_ground <- femur_data_ground$DE
names(femur_DE_ground) <- rownames(femur_data_ground)
femur_CSS_ground <- femur_data_ground$CSS
names(femur_CSS_ground) <- rownames(femur_data_ground)
femur_shape_data_ground <- filter(femur_shape_data, ecotype == "ground")
femur_lnsize_ground <- femur_shape_data_ground$lnfemur_size
names(femur_lnsize_ground) <- rownames(femur_shape_data_ground) 
femur_shape_ground <- arrayspecs(femur_shape_data_ground[,3:458], 152,3 )
femur_shape_ground_2D <- two.d.array(femur_shape_ground)
lnfemur_size_ground <- femur_shape_data_ground$lnfemur_size
tree_pr_ground <- treedata(phy = tree, data = femur_shape_data_ground, sort=TRUE)$phy

#tree data
femur_data_tree <- filter(femur_data, ecotype == "tree")
femur_Cg_tree <- femur_data_tree$Cg
names(femur_Cg_tree) <- rownames(femur_data_tree)
femur_DE_tree <- femur_data_tree$DE
names(femur_DE_tree) <- rownames(femur_data_tree)
femur_CSS_tree <- femur_data_tree$CSS
names(femur_CSS_tree) <- rownames(femur_data_tree)
femur_shape_data_tree <- filter(femur_shape_data, ecotype == "tree")
femur_lnsize_tree <- femur_shape_data_tree$lnfemur_size
names(femur_lnsize_tree) <- rownames(femur_shape_data_tree) 
femur_shape_tree <- arrayspecs(femur_shape_data_tree[,3:458], 152,3 )
femur_shape_tree_2D <- two.d.array(femur_shape_tree)
lnfemur_size_tree <- femur_shape_data_tree$lnfemur_size
tree_pr_tree <- treedata(phy = tree, data = femur_shape_data_tree, sort=TRUE)$phy



#### femur_shape ####
##### phylogenetic PCA ##### 
femur_shape_pca <- gm.prcomp(femur_shape, phy = tree_pr, GLS = FALSE)
summary(femur_shape_pca)
plot(femur_shape_pca, phylo = TRUE, main = "phylo PCA")

col_ecotype_phylomorph <-c(col_ecotype[as.factor(ecotype)],rep("white",tree_pr$Nnode))
names(col_ecotype_phylomorph)<-1:(length(tree_pr$tip)+tree_pr$Nnode)
phylomorphospace(tree_pr, femur_shape_pca$x[,1:2], control = list(col.node = col_ecotype_phylomorph), label = "off")

# femur phy signal
femur_physig <- physignal(femur_shape, phy = tree_pr, iter = 999)

##### null model: shape ~ 1 #####
femur_shape_null_mv <- mvgls(femur_shape_2D ~ 1, tree = tree_pr, method = "LOOCV", model = "lambda")
save(femur_shape_null_mv, file= "femur_shape_null_mv.Rdata")
femur_shape_null_mv_EIC <- EIC(femur_shape_null_mv)
save(femur_shape_null_mv_EIC, file = "femur_shape_null_mv_EIC.Rdata")


##### size model: shape ~ size #####
femur_shape_allometry_mv <- mvgls(femur_shape_2D ~ lnfemur_size, tree = tree_pr, method = "LOOCV", model = "lambda")
save(femur_shape_allometry_mv, file= "femur_shape_allometry_mv.Rdata")
femur_shape_allometry_mv_EIC <- EIC(femur_shape_allometry_mv, nboot = 1000)
save(femur_shape_allometry_mv_EIC, file= "femur_shape_allometry_mv_EIC.Rdata")
femur_shape_allometry_mv_test <- manova.gls(femur_shape_allometry_mv, type = "II", test = "Pillai")
save(femur_shape_allometry_mv_test, file= "femur_shape_allometry_mv_test.Rdata")


##### ecotype*size model: shape ~ ecotype*size   #####
femur_ecotype_mvpancova <- mvgls(femur_shape_2D ~ lnfemur_size*ecotype, tree = tree_pr, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpancova, file= "femur_ecotype_mvpancova.Rdata")
femur_ecotype_mvpancova_EIC <- EIC(femur_ecotype_mvpancova)
save(femur_ecotype_mvpancova_EIC, file= "femur_ecotype_mvpancova_EIC.Rdata")
femur_ecotype_mvpancova_test <- manova.gls(femur_ecotype_mvpancova, type = "II", test = "Pillai")
save(femur_ecotype_mvpancova_test, file= "femur_ecotype_mvpancova_test.Rdata")
femur_ecotype_mvpancova_pw <- pairwise.glh(femur_ecotype_mvpancova, term="ecotype", test="Pillai", adjust="holm", nperm=1000, verbose=TRUE)
save(femur_ecotype_mvpancova_pw, file= "femur_ecotype_mvpancova_pw.Rdata")


##### ecotype+size model: shape ~ ecotype+size   #####
femur_ecotype_mvpancova_noint <- mvgls(femur_shape_2D ~ lnfemur_size+ecotype, tree = tree_pr, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpancova_noint, file= "femur_ecotype_mvpancova_noint.Rdata")
femur_ecotype_mvpancova_noint_EIC <- EIC(femur_ecotype_mvpancova_noint, nboot = 1000)
save(femur_ecotype_mvpancova_noint_EIC, file= "femur_ecotype_mvpancova_noint_EIC.Rdata")
femur_ecotype_mvpancova_noint_test <- manova.gls(femur_ecotype_mvpancova_noint, type = "II", test = "Pillai")
save(femur_ecotype_mvpancova_noint_test, file= "femur_ecotype_mvpancova_noint_test.Rdata")


##### ecotype model: shape ~ ecotype #####
femur_ecotype_mvpanova <- mvgls(femur_shape_2D ~ ecotype, tree = tree_pr, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpanova, file= "femur_ecotype_mvpanova.Rdata")
femur_ecotype_mvpanova_EIC <- EIC(femur_ecotype_mvpanova, nboot = 1000)
save(femur_ecotype_mvpanova_EIC, file= "femur_ecotype_mvpanova_EIC.Rdata")
femur_ecotype_mvpanova_test <- manova.gls(femur_ecotype_mvpanova, test = "Pillai")
save(femur_ecotype_mvpanova_test, file= "femur_ecotype_mvpanova_test.Rdata")

#panova
femur_ecotype_panova <- procD.pgls(femur_shape ~ ecotype, phy = tree_pr, iter = 999)
summary(femur_ecotype_panova)
femur_ecotype_pw <- pairwise(femur_ecotype_panova, groups = ecotype)
summary(femur_ecotype_pw)

#anova
femur_ecotype_anova <- procD.lm(femur_shape ~ ecotype, iter = 999)
summary(femur_ecotype_anova)
femur_ecotype_pancova_pw <- pairwise(femur_ecotype_anova, groups = ecotype)
summary(femur_ecotype_pancova_pw)

##### femur_shape AIC #####
geiger::aicw(c(femur_shape_allometry_mv_EIC$EIC,
               femur_ecotype_mvpanova_EIC$EIC,
               femur_ecotype_mvpancova_EIC$EIC,
               femur_ecotype_mvpancova_noint_EIC$EIC,
               femur_shape_null_mv_EIC$EIC))

##### within ecotypes allometry #####
###### chipmunk ######
femur_ecotype_mvpancova_chipmunk <- mvgls(femur_shape_chipmunk_2D ~ lnfemur_size_chipmunk, tree = tree_pr_chipmunk, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpancova_chipmunk, file= "femur_ecotype_mvpancova_chipmunk.Rdata")
femur_ecotype_mvpancova_test_chipmunk <- manova.gls(femur_ecotype_mvpancova_chipmunk, type = "II", test = "Pillai")
save(femur_ecotype_mvpancova_test_chipmunk, file= "femur_ecotype_mvpancova_test_chipmunk.Rdata")

femur_shape_allometry_chipmunk <- procD.pgls(femur_shape_chipmunk ~ lnfemur_size_chipmunk, phy = tree_pr_chipmunk,  iter = 999)
summary(femur_shape_allometry_chipmunk)
femur_allo_plot <- plotAllometry(femur_shape_allometry_chipmunk, size = lnfemur_size_chipmunk, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

femur_shape_allometry_chipmunk_np <- procD.lm(femur_shape_chipmunk ~ lnfemur_size_chipmunk, iter = 999)
summary(femur_shape_allometry_chipmunk_np)
femur_allo_plot_np <- plotAllometry(femur_shape_allometry_chipmunk_np, size = lnfemur_size_chipmunk, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

###### gliding ######
femur_shape_data_gliding <- filter(femur_shape_data, ecotype == "gliding")
femur_shape_gliding <- arrayspecs(femur_shape_data_gliding[,3:458], 152,3 )
femur_shape_gliding_2D <- two.d.array(femur_shape_gliding)
lnfemur_size_gliding <- femur_shape_data_gliding$lnfemur_size
tree_pr_gliding <- treedata(phy = tree, data = femur_shape_data_gliding, sort=TRUE)$phy

femur_ecotype_mvpancova_gliding <- mvgls(femur_shape_gliding_2D ~ lnfemur_size_gliding, tree = tree_pr_gliding, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpancova_gliding, file= "femur_ecotype_mvpancova_gliding.Rdata")
femur_ecotype_mvpancova_test_gliding <- manova.gls(femur_ecotype_mvpancova_gliding, type = "II", test = "Pillai")
save(femur_ecotype_mvpancova_test_gliding, file= "femur_ecotype_mvpancova_test_gliding.Rdata")

femur_shape_allometry_gliding <- procD.pgls(femur_shape_gliding ~ lnfemur_size_gliding, phy = tree_pr_gliding,  iter = 999)
summary(femur_shape_allometry_gliding)
femur_allo_plot <- plotAllometry(femur_shape_allometry_gliding, size = lnfemur_size_gliding, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

femur_shape_allometry_gliding_np <- procD.lm(femur_shape_gliding ~ lnfemur_size_gliding, iter = 999)
summary(femur_shape_allometry_gliding_np)
femur_allo_plot_np <- plotAllometry(femur_shape_allometry_gliding_np, size = lnfemur_size_gliding, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

###### ground ######
femur_shape_data_ground <- filter(femur_shape_data, ecotype == "ground")
femur_shape_ground <- arrayspecs(femur_shape_data_ground[,3:458], 152,3 )
femur_shape_ground_2D <- two.d.array(femur_shape_ground)
lnfemur_size_ground <- femur_shape_data_ground$lnfemur_size
tree_pr_ground <- treedata(phy = tree, data = femur_shape_data_ground, sort=TRUE)$phy

femur_ecotype_mvpancova_ground <- mvgls(femur_shape_ground_2D ~ lnfemur_size_ground, tree = tree_pr_ground, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpancova_ground, file= "femur_ecotype_mvpancova_ground.Rdata")
femur_ecotype_mvpancova_test_ground <- manova.gls(femur_ecotype_mvpancova_ground, type = "II", test = "Pillai")
save(femur_ecotype_mvpancova_test_ground, file= "femur_ecotype_mvpancova_test_ground.Rdata")

femur_shape_allometry_ground <- procD.pgls(femur_shape_ground ~ lnfemur_size_ground, phy = tree_pr_ground,  iter = 999)
summary(femur_shape_allometry_ground)
femur_allo_plot <- plotAllometry(femur_shape_allometry_ground, size = lnfemur_size_ground, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

femur_shape_allometry_ground_np <- procD.lm(femur_shape_ground ~ lnfemur_size_ground, iter = 999)
summary(femur_shape_allometry_ground_np)
femur_allo_plot_np <- plotAllometry(femur_shape_allometry_ground_np, size = lnfemur_size_ground, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

###### tree ######
femur_shape_data_tree <- filter(femur_shape_data, ecotype == "tree")
femur_shape_tree <- arrayspecs(femur_shape_data_tree[,3:458], 152,3 )
femur_shape_tree_2D <- two.d.array(femur_shape_tree)
lnfemur_size_tree <- femur_shape_data_tree$lnfemur_size
tree_pr_tree <- treedata(phy = tree, data = femur_shape_data_tree, sort=TRUE)$phy

femur_ecotype_mvpancova_tree <- mvgls(femur_shape_tree_2D ~ lnfemur_size_tree, tree = tree_pr_tree, method = "LOOCV", model = "lambda")
save(femur_ecotype_mvpancova_tree, file= "femur_ecotype_mvpancova_tree.Rdata")
femur_ecotype_mvpancova_test_tree <- manova.gls(femur_ecotype_mvpancova_tree, type = "II", test = "Pillai")
save(femur_ecotype_mvpancova_test_tree, file= "femur_ecotype_mvpancova_test_tree.Rdata")

femur_shape_allometry_tree <- procD.pgls(femur_shape_tree ~ lnfemur_size_tree, phy = tree_pr_tree,  iter = 999)
summary(femur_shape_allometry_tree)
femur_allo_plot <- plotAllometry(femur_shape_allometry_tree, size = lnfemur_size_tree, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

femur_shape_allometry_tree_np <- procD.lm(femur_shape_tree ~ lnfemur_size_tree, iter = 999)
summary(femur_shape_allometry_tree_np)
femur_allo_plot_np <- plotAllometry(femur_shape_allometry_tree_np, size = lnfemur_size_tree, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)




#### femur_Cg #### 
##### null model: femur_Cg ~ 1 #####
femur_Cg_pgls_null <- phylolm(femur_Cg~1, phy = tree_pr, model = "lambda", boot = 1000)


##### size model: femur_Cg ~ size #####
femur_Cg_pgls_size <- phylolm(femur_Cg~lnfemur_size, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_Cg_pgls_size)


##### size+ecotype model: femur_Cg ~ size + ecotype #####
femur_Cg_pANCOVAnoint_size_ecotype <- phylolm(femur_Cg~lnfemur_size+ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_Cg_pANCOVAnoint_size_ecotype)


##### size*ecotype model: femur_Cg ~ size * ecotype #####
femur_Cg_pANCOVA_size_ecotype <- phylolm(femur_Cg~lnfemur_size*ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_Cg_pANCOVA_size_ecotype)

######  Make table of coefficients (intercept and slope) ######  
femur_Cg_pgls_size_all_int <- femur_Cg_pgls_size$coefficients["(Intercept)"]
femur_Cg_pgls_size_all_int_boot <- femur_Cg_pgls_size$bootstrap[,"(Intercept)"]

femur_Cg_pANCOVA_size_ecotype_chip_int <- femur_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]
femur_Cg_pANCOVA_size_ecotype_glide_int <- femur_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_Cg_pANCOVA_size_ecotype$coefficients["ecotypegliding"]
femur_Cg_pANCOVA_size_ecotype_ground_int <- femur_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_Cg_pANCOVA_size_ecotype$coefficients["ecotypeground"]
femur_Cg_pANCOVA_size_ecotype_tree_int <- femur_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_Cg_pANCOVA_size_ecotype$coefficients["ecotypetree"]
femur_Cg_pANCOVA_size_ecotype_chip_int_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"]
femur_Cg_pANCOVA_size_ecotype_glide_int_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_Cg_pANCOVA_size_ecotype$bootstrap[,"ecotypegliding"]
femur_Cg_pANCOVA_size_ecotype_ground_int_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_Cg_pANCOVA_size_ecotype$bootstrap[,"ecotypeground"]
femur_Cg_pANCOVA_size_ecotype_tree_int_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_Cg_pANCOVA_size_ecotype$bootstrap[,"ecotypetree"]

femur_Cg_pgls_size_all_slope <- femur_Cg_pgls_size$coefficients["lnfemur_size"]
femur_Cg_pgls_size_all_slope_boot <- femur_Cg_pgls_size$bootstrap[,"lnfemur_size"]

femur_Cg_pANCOVA_size_ecotype_chip_slope <- femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size"]
femur_Cg_pANCOVA_size_ecotype_glide_slope <- femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypegliding"]
femur_Cg_pANCOVA_size_ecotype_ground_slope <- femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypeground"]
femur_Cg_pANCOVA_size_ecotype_tree_slope <- femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_Cg_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypetree"]
femur_Cg_pANCOVA_size_ecotype_chip_slope_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"]
femur_Cg_pANCOVA_size_ecotype_glide_slope_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypegliding"]
femur_Cg_pANCOVA_size_ecotype_ground_slope_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypeground"]
femur_Cg_pANCOVA_size_ecotype_tree_slope_boot <- femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_Cg_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypetree"]

femur_Cg_pANCOVA_size_ecotype_int_cof <- rbind(
  femur_Cg_pgls_size_all_int,
  femur_Cg_pANCOVA_size_ecotype_chip_int,
  femur_Cg_pANCOVA_size_ecotype_glide_int,
  femur_Cg_pANCOVA_size_ecotype_ground_int,
  femur_Cg_pANCOVA_size_ecotype_tree_int)
femur_Cg_pANCOVA_size_ecotype_int_95 <- rbind(
  quantile(femur_Cg_pgls_size_all_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_tree_int_boot, prob = c(0.025, 0.975)))
femur_Cg_pANCOVA_size_ecotype_slope_cof <- rbind(
  femur_Cg_pgls_size_all_slope,
  femur_Cg_pANCOVA_size_ecotype_chip_slope,
  femur_Cg_pANCOVA_size_ecotype_glide_slope,
  femur_Cg_pANCOVA_size_ecotype_ground_slope,
  femur_Cg_pANCOVA_size_ecotype_tree_slope)
femur_Cg_pANCOVA_size_ecotype_slope_95 <- rbind(
  quantile(femur_Cg_pgls_size_all_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot, prob = c(0.025, 0.975)))

femur_Cg_pANCOVA_size_ecotype_cof <- data.frame(femur_Cg_pANCOVA_size_ecotype_int_cof, femur_Cg_pANCOVA_size_ecotype_int_95, femur_Cg_pANCOVA_size_ecotype_slope_cof, femur_Cg_pANCOVA_size_ecotype_slope_95 )
colnames(femur_Cg_pANCOVA_size_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95", "Slope", "Slope_L95", "Slope_U95")
rownames(femur_Cg_pANCOVA_size_ecotype_cof) <- c("all", "chip", "gliding", "ground", "tree")
round(femur_Cg_pANCOVA_size_ecotype_cof, digits = 3)

######   plot femur_Cg ~ lnsize*ecology ######  
plot(femur_Cg ~ lnfemur_size, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)
abline(a = femur_Cg_pgls_size_all_int, b = femur_Cg_pgls_size_all_slope, col = "black")
abline(a = femur_Cg_pANCOVA_size_ecotype_chip_int, b = femur_Cg_pANCOVA_size_ecotype_chip_slope, col = "#073b4c", lty = 2)
abline(a = femur_Cg_pANCOVA_size_ecotype_glide_int, b = femur_Cg_pANCOVA_size_ecotype_glide_slope, col = "#118ab2", lty = 2)
abline(a = femur_Cg_pANCOVA_size_ecotype_ground_int, b = femur_Cg_pANCOVA_size_ecotype_ground_slope, col = "#ffd166", lty = 2)
abline(a = femur_Cg_pANCOVA_size_ecotype_tree_int, b = femur_Cg_pANCOVA_size_ecotype_tree_slope, col = "#06d6a0")

######   slope differences between ecotypes ######  
# chip - glide
femur_Cg_slope_size_chip_glide_SE <- sqrt((sd(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)))

femur_Cg_slope_size_chip_glide_meandiff <- mean(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)

femur_Cg_slope_size_chip_glide_CI <- mean(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot) + c(-1, 1) * 1.96 * femur_Cg_slope_size_chip_glide_SE 

femur_Cg_slope_size_chip_glide_table <- cbind(femur_Cg_slope_size_chip_glide_meandiff, femur_Cg_slope_size_chip_glide_CI[1], femur_Cg_slope_size_chip_glide_CI[2])
rownames(femur_Cg_slope_size_chip_glide_table) <- ("femur_Cg_slope_size_chip_glide")
colnames(femur_Cg_slope_size_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_Cg_slope_size_chip_glide_table

# chip - ground
femur_Cg_slope_size_chip_ground_SE <- sqrt((sd(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)))

femur_Cg_slope_size_chip_ground_meandiff <- mean(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)

femur_Cg_slope_size_chip_ground_CI <- mean(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * femur_Cg_slope_size_chip_ground_SE  

femur_Cg_slope_size_chip_ground_table <- cbind(femur_Cg_slope_size_chip_ground_meandiff, femur_Cg_slope_size_chip_ground_CI[1], femur_Cg_slope_size_chip_ground_CI[2])
rownames(femur_Cg_slope_size_chip_ground_table) <- ("femur_Cg_slope_size_chip_ground")
colnames(femur_Cg_slope_size_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_Cg_slope_size_chip_ground_table

# chip - tree
femur_Cg_slope_size_chip_tree_SE <- sqrt((sd(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)))

femur_Cg_slope_size_chip_tree_meandiff <- mean(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)

femur_Cg_slope_size_chip_tree_CI <- mean(femur_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_Cg_slope_size_chip_tree_SE  

femur_Cg_slope_size_chip_tree_table <- cbind(femur_Cg_slope_size_chip_tree_meandiff, femur_Cg_slope_size_chip_tree_CI[1], femur_Cg_slope_size_chip_tree_CI[2])
rownames(femur_Cg_slope_size_chip_tree_table) <- ("femur_Cg_slope_size_chip_tree")
colnames(femur_Cg_slope_size_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_Cg_slope_size_chip_tree_table

# glide - ground
femur_Cg_slope_size_glide_ground_SE <- sqrt((sd(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)))

femur_Cg_slope_size_glide_ground_meandiff <- mean(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)

femur_Cg_slope_size_glide_ground_CI <- mean(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * femur_Cg_slope_size_glide_ground_SE  

femur_Cg_slope_size_glide_ground_table <- cbind(femur_Cg_slope_size_glide_ground_meandiff, femur_Cg_slope_size_glide_ground_CI[1], femur_Cg_slope_size_glide_ground_CI[2])
rownames(femur_Cg_slope_size_glide_ground_table) <- ("femur_Cg_slope_size_glide_ground")
colnames(femur_Cg_slope_size_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_Cg_slope_size_glide_ground_table

# glide - tree
femur_Cg_slope_size_glide_tree_SE <- sqrt((sd(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)))

femur_Cg_slope_size_glide_tree_meandiff <- mean(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)

femur_Cg_slope_size_glide_tree_CI <- mean(femur_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_Cg_slope_size_glide_tree_SE  

femur_Cg_slope_size_glide_tree_table <- cbind(femur_Cg_slope_size_glide_tree_meandiff, femur_Cg_slope_size_glide_tree_CI[1], femur_Cg_slope_size_glide_tree_CI[2])
rownames(femur_Cg_slope_size_glide_tree_table) <- ("femur_Cg_slope_size_glide_tree")
colnames(femur_Cg_slope_size_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_Cg_slope_size_glide_tree_table

# ground - tree
femur_Cg_slope_size_ground_tree_SE <- sqrt((sd(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot)) + (sd(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)))

femur_Cg_slope_size_ground_tree_meandiff <- mean(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot)

femur_Cg_slope_size_ground_tree_CI <- mean(femur_Cg_pANCOVA_size_ecotype_ground_slope_boot) - mean(femur_Cg_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_Cg_slope_size_ground_tree_SE 

femur_Cg_slope_size_ground_tree_table <- cbind(femur_Cg_slope_size_ground_tree_meandiff, femur_Cg_slope_size_ground_tree_CI[1], femur_Cg_slope_size_ground_tree_CI[2])
rownames(femur_Cg_slope_size_ground_tree_table) <- ("femur_Cg_slope_size_ground_tree")
colnames(femur_Cg_slope_size_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_Cg_slope_size_ground_tree_table

#table
femur_Cg_slope_size_meandiff_table <- rbind(
  femur_Cg_slope_size_chip_glide_table,
  femur_Cg_slope_size_chip_ground_table,
  femur_Cg_slope_size_chip_tree_table,
  femur_Cg_slope_size_glide_ground_table,
  femur_Cg_slope_size_glide_tree_table,
  femur_Cg_slope_size_ground_tree_table)
round(femur_Cg_slope_size_meandiff_table, digits = 2)

##### ecotype model: femur_Cg ~ ecotype ####
femur_Cg_pANOVA_ecotype <- phylolm(femur_Cg~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_Cg_pANOVA_ecotype)

###### Make table of coefficients (intercept) ###### 
femur_Cg_pANOVA_ecotype_chip_int <- femur_Cg_pANOVA_ecotype$coefficients["(Intercept)"]
femur_Cg_pANOVA_ecotype_glide_int <- femur_Cg_pANOVA_ecotype$coefficients["(Intercept)"]+femur_Cg_pANOVA_ecotype$coefficients["ecotypegliding"]
femur_Cg_pANOVA_ecotype_ground_int <- femur_Cg_pANOVA_ecotype$coefficients["(Intercept)"]+femur_Cg_pANOVA_ecotype$coefficients["ecotypeground"]
femur_Cg_pANOVA_ecotype_tree_int <- femur_Cg_pANOVA_ecotype$coefficients["(Intercept)"]+femur_Cg_pANOVA_ecotype$coefficients["ecotypetree"]
femur_Cg_pANOVA_ecotype_chip_int_boot <- femur_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"]
femur_Cg_pANOVA_ecotype_glide_int_boot <- femur_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_Cg_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
femur_Cg_pANOVA_ecotype_ground_int_boot <- femur_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_Cg_pANOVA_ecotype$bootstrap[,"ecotypeground"]
femur_Cg_pANOVA_ecotype_tree_int_boot <- femur_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_Cg_pANOVA_ecotype$bootstrap[,"ecotypetree"]

femur_Cg_pANOVA_ecotype_int_cof <- rbind(
  femur_Cg_pANOVA_ecotype_chip_int,
  femur_Cg_pANOVA_ecotype_glide_int,
  femur_Cg_pANOVA_ecotype_ground_int,
  femur_Cg_pANOVA_ecotype_tree_int)
femur_Cg_pANOVA_ecotype_int_95 <- rbind(
  quantile(femur_Cg_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_Cg_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

femur_Cg_pANOVA_ecotype_cof <- data.frame(femur_Cg_pANOVA_ecotype_int_cof, femur_Cg_pANOVA_ecotype_int_95)
colnames(femur_Cg_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(femur_Cg_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(femur_Cg_pANOVA_ecotype_cof, digits = 3)

######  Plot femur_Cg ~ ecotype ###### 
(femur_Cg_vp_ecology <- ggplot(femur_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= femur_Cg, fill = ecotype)) + 
  scale_y_continuous(name="ln femur_Cg (mm)", limits = c(.45,.85), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(femur_Cg)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill = "grey") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color = "black", size = 13, angle = 0, hjust = .5, vjust = .5),
        axis.text.y=element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position= "none"))


######  femur_Cg_int_meandiff_table ###### 
# chip - glide
femur_Cg_int_chip_glide_SE <- sqrt((sd(femur_Cg_pANOVA_ecotype_chip_int_boot)^2/length(femur_Cg_pANOVA_ecotype_chip_int_boot)) + (sd(femur_Cg_pANOVA_ecotype_glide_int_boot)^2/length(femur_Cg_pANOVA_ecotype_glide_int_boot)))

femur_Cg_int_chip_glide_meandiff <- mean(femur_Cg_pANOVA_ecotype_chip_int_boot) - mean(femur_Cg_pANOVA_ecotype_glide_int_boot)

femur_Cg_int_chip_glide_CI <- mean(femur_Cg_pANOVA_ecotype_chip_int_boot) - mean(femur_Cg_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * femur_Cg_int_chip_glide_SE  # 95% confidence interval using z=1.96

femur_Cg_int_chip_glide_table <- cbind(femur_Cg_int_chip_glide_meandiff, femur_Cg_int_chip_glide_CI[1], femur_Cg_int_chip_glide_CI[2])
rownames(femur_Cg_int_chip_glide_table) <- ("femur_Cg_int_chip_glide")
colnames(femur_Cg_int_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_Cg_int_chip_glide_table

# chip - ground
femur_Cg_int_chip_ground_SE <- sqrt((sd(femur_Cg_pANOVA_ecotype_chip_int_boot)^2/length(femur_Cg_pANOVA_ecotype_chip_int_boot)) + (sd(femur_Cg_pANOVA_ecotype_ground_int_boot)^2/length(femur_Cg_pANOVA_ecotype_ground_int_boot)))

femur_Cg_int_chip_ground_meandiff <- mean(femur_Cg_pANOVA_ecotype_chip_int_boot) - mean(femur_Cg_pANOVA_ecotype_ground_int_boot)

femur_Cg_int_chip_ground_CI <- mean(femur_Cg_pANOVA_ecotype_chip_int_boot) - mean(femur_Cg_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_Cg_int_chip_ground_SE  # 95% confidence interval using z=1.96

femur_Cg_int_chip_ground_table <- cbind(femur_Cg_int_chip_ground_meandiff, femur_Cg_int_chip_ground_CI[1], femur_Cg_int_chip_ground_CI[2])
rownames(femur_Cg_int_chip_ground_table) <- ("femur_Cg_int_chip_ground")
colnames(femur_Cg_int_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_Cg_int_chip_ground_table

# chip - tree
femur_Cg_int_chip_tree_SE <- sqrt((sd(femur_Cg_pANOVA_ecotype_chip_int_boot)^2/length(femur_Cg_pANOVA_ecotype_chip_int_boot)) + (sd(femur_Cg_pANOVA_ecotype_tree_int_boot)^2/length(femur_Cg_pANOVA_ecotype_tree_int_boot)))

femur_Cg_int_chip_tree_meandiff <- mean(femur_Cg_pANOVA_ecotype_chip_int_boot) - mean(femur_Cg_pANOVA_ecotype_tree_int_boot)

femur_Cg_int_chip_tree_CI <- mean(femur_Cg_pANOVA_ecotype_chip_int_boot) - mean(femur_Cg_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_Cg_int_chip_tree_SE  #

femur_Cg_int_chip_tree_table <- cbind(femur_Cg_int_chip_tree_meandiff, femur_Cg_int_chip_tree_CI[1], femur_Cg_int_chip_tree_CI[2])
rownames(femur_Cg_int_chip_tree_table) <- ("femur_Cg_int_chip_tree")
colnames(femur_Cg_int_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_Cg_int_chip_tree_table

# glide - ground
femur_Cg_int_glide_ground_SE <- sqrt((sd(femur_Cg_pANOVA_ecotype_glide_int_boot)^2/length(femur_Cg_pANOVA_ecotype_glide_int_boot)) + (sd(femur_Cg_pANOVA_ecotype_ground_int_boot)^2/length(femur_Cg_pANOVA_ecotype_ground_int_boot)))

femur_Cg_int_glide_ground_meandiff <- mean(femur_Cg_pANOVA_ecotype_glide_int_boot) - mean(femur_Cg_pANOVA_ecotype_ground_int_boot)

femur_Cg_int_glide_ground_CI <- mean(femur_Cg_pANOVA_ecotype_glide_int_boot) - mean(femur_Cg_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_Cg_int_glide_ground_SE  

femur_Cg_int_glide_ground_table <- cbind(femur_Cg_int_glide_ground_meandiff, femur_Cg_int_glide_ground_CI[1], femur_Cg_int_glide_ground_CI[2])
rownames(femur_Cg_int_glide_ground_table) <- ("femur_Cg_int_glide_ground")
colnames(femur_Cg_int_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_Cg_int_glide_ground_table

# glide - tree
femur_Cg_int_glide_tree_SE <- sqrt((sd(femur_Cg_pANOVA_ecotype_glide_int_boot)^2/length(femur_Cg_pANOVA_ecotype_glide_int_boot)) + (sd(femur_Cg_pANOVA_ecotype_tree_int_boot)^2/length(femur_Cg_pANOVA_ecotype_tree_int_boot)))

femur_Cg_int_glide_tree_meandiff <- mean(femur_Cg_pANOVA_ecotype_glide_int_boot) - mean(femur_Cg_pANOVA_ecotype_tree_int_boot)

femur_Cg_int_glide_tree_CI <- mean(femur_Cg_pANOVA_ecotype_glide_int_boot) - mean(femur_Cg_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_Cg_int_glide_tree_SE  

femur_Cg_int_glide_tree_table <- cbind(femur_Cg_int_glide_tree_meandiff, femur_Cg_int_glide_tree_CI[1], femur_Cg_int_glide_tree_CI[2])
rownames(femur_Cg_int_glide_tree_table) <- ("femur_Cg_int_glide_tree")
colnames(femur_Cg_int_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_Cg_int_glide_tree_table

# ground - tree
femur_Cg_int_ground_tree_SE <- sqrt((sd(femur_Cg_pANOVA_ecotype_ground_int_boot)^2/length(femur_Cg_pANOVA_ecotype_ground_int_boot)) + (sd(femur_Cg_pANOVA_ecotype_tree_int_boot)^2/length(femur_Cg_pANOVA_ecotype_tree_int_boot)))

femur_Cg_int_ground_tree_meandiff <- mean(femur_Cg_pANOVA_ecotype_ground_int_boot) - mean(femur_Cg_pANOVA_ecotype_tree_int_boot)

femur_Cg_int_ground_tree_CI <- mean(femur_Cg_pANOVA_ecotype_ground_int_boot) - mean(femur_Cg_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_Cg_int_ground_tree_SE  # 95% confidence interval using z=1.96

femur_Cg_int_ground_tree_table <- cbind(femur_Cg_int_ground_tree_meandiff, femur_Cg_int_ground_tree_CI[1], femur_Cg_int_ground_tree_CI[2])
rownames(femur_Cg_int_ground_tree_table) <- ("femur_Cg_int_ground_tree")
colnames(femur_Cg_int_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_Cg_int_ground_tree_table

#table
femur_Cg_int_meandiff_table <- rbind(
  femur_Cg_int_chip_glide_table,
  femur_Cg_int_chip_ground_table,
  femur_Cg_int_chip_tree_table,
  femur_Cg_int_glide_ground_table,
  femur_Cg_int_glide_tree_table,
  femur_Cg_int_ground_tree_table)
round(femur_Cg_int_meandiff_table, digits = 2)


##### femur_Cg AIC #####
geiger::aicw(c(AIC(femur_Cg_pgls_size),
               AIC(femur_Cg_pANOVA_ecotype),
               AIC(femur_Cg_pANCOVA_size_ecotype),
               AIC(femur_Cg_pANCOVAnoint_size_ecotype),
               AIC(femur_Cg_pgls_null)))


#### femur_DE #####
##### null model: femur_DE ~ 1 ##### 
femur_DE_pgls_null <- phylolm(femur_DE~1, phy = tree_pr, model = "lambda", boot = 1000)


##### size model: femur_DE ~ size ##### 
femur_DE_pgls_size <- phylolm(femur_DE~lnfemur_size, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_DE_pgls_size)


##### size + ecotype model: femur_DE ~ size + ecotype ##### 
femur_DE_pANCOVAnoint_size_ecotype <- phylolm(femur_DE~lnfemur_size+ecotype, phy = tree_pr, model = "lambda", boot = 1000, lower.bound = 0)
summary(femur_DE_pANCOVAnoint_size_ecotype)


##### size * ecotype model: femur_DE ~ size * ecotype ##### 
femur_DE_pANCOVA_size_ecotype <- phylolm(femur_DE~lnfemur_size*ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_DE_pANCOVA_size_ecotype)

###### Make table of coefficients (intercept and slope) ###### 
femur_DE_pgls_size_all_int <- femur_DE_pgls_size$coefficients["(Intercept)"]
femur_DE_pgls_size_all_int_boot <- femur_DE_pgls_size$bootstrap[,"(Intercept)"]

femur_DE_pANCOVA_size_ecotype_chip_int <- femur_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]
femur_DE_pANCOVA_size_ecotype_glide_int <- femur_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_DE_pANCOVA_size_ecotype$coefficients["ecotypegliding"]
femur_DE_pANCOVA_size_ecotype_ground_int <- femur_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_DE_pANCOVA_size_ecotype$coefficients["ecotypeground"]
femur_DE_pANCOVA_size_ecotype_tree_int <- femur_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_DE_pANCOVA_size_ecotype$coefficients["ecotypetree"]
femur_DE_pANCOVA_size_ecotype_chip_int_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"]
femur_DE_pANCOVA_size_ecotype_glide_int_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_DE_pANCOVA_size_ecotype$bootstrap[,"ecotypegliding"]
femur_DE_pANCOVA_size_ecotype_ground_int_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_DE_pANCOVA_size_ecotype$bootstrap[,"ecotypeground"]
femur_DE_pANCOVA_size_ecotype_tree_int_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_DE_pANCOVA_size_ecotype$bootstrap[,"ecotypetree"]

femur_DE_pgls_size_all_slope <- femur_DE_pgls_size$coefficients["lnfemur_size"]
femur_DE_pgls_size_all_slope_boot <- femur_DE_pgls_size$bootstrap[,"lnfemur_size"]

femur_DE_pANCOVA_size_ecotype_chip_slope <- femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size"]
femur_DE_pANCOVA_size_ecotype_glide_slope <- femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypegliding"]
femur_DE_pANCOVA_size_ecotype_ground_slope <- femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypeground"]
femur_DE_pANCOVA_size_ecotype_tree_slope <- femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_DE_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypetree"]
femur_DE_pANCOVA_size_ecotype_chip_slope_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"]
femur_DE_pANCOVA_size_ecotype_glide_slope_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypegliding"]
femur_DE_pANCOVA_size_ecotype_ground_slope_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypeground"]
femur_DE_pANCOVA_size_ecotype_tree_slope_boot <- femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_DE_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypetree"]

femur_DE_pANCOVA_size_ecotype_int_cof <- rbind(
  femur_DE_pgls_size_all_int,
  femur_DE_pANCOVA_size_ecotype_chip_int,
  femur_DE_pANCOVA_size_ecotype_glide_int,
  femur_DE_pANCOVA_size_ecotype_ground_int,
  femur_DE_pANCOVA_size_ecotype_tree_int)
femur_DE_pANCOVA_size_ecotype_int_95 <- rbind(
  quantile(femur_DE_pgls_size_all_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_tree_int_boot, prob = c(0.025, 0.975)))
femur_DE_pANCOVA_size_ecotype_slope_cof <- rbind(
  femur_DE_pgls_size_all_slope,
  femur_DE_pANCOVA_size_ecotype_chip_slope,
  femur_DE_pANCOVA_size_ecotype_glide_slope,
  femur_DE_pANCOVA_size_ecotype_ground_slope,
  femur_DE_pANCOVA_size_ecotype_tree_slope)
femur_DE_pANCOVA_size_ecotype_slope_95 <- rbind(
  quantile(femur_DE_pgls_size_all_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_chip_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_glide_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_ground_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANCOVA_size_ecotype_tree_slope_boot, prob = c(0.025, 0.975)))

femur_DE_pANCOVA_size_ecotype_cof <- data.frame(femur_DE_pANCOVA_size_ecotype_int_cof, femur_DE_pANCOVA_size_ecotype_int_95, femur_DE_pANCOVA_size_ecotype_slope_cof, femur_DE_pANCOVA_size_ecotype_slope_95 )
colnames(femur_DE_pANCOVA_size_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95", "Slope", "Slope_L95", "Slope_U95")
rownames(femur_DE_pANCOVA_size_ecotype_cof) <- c("all", "chip", "gliding", "ground", "tree")
round(femur_DE_pANCOVA_size_ecotype_cof, digits = 3)

######  plot femur_DE ~ lnsize*ecology ###### 
plot(femur_DE ~ lnfemur_size, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)
abline(a = femur_DE_pgls_size_all_int, b = femur_DE_pgls_size_all_slope, col = "black", lty = 2)
abline(a = femur_DE_pANCOVA_size_ecotype_chip_int, b = femur_DE_pANCOVA_size_ecotype_chip_slope, col = "#073b4c", lty = 2)
abline(a = femur_DE_pANCOVA_size_ecotype_glide_int, b = femur_DE_pANCOVA_size_ecotype_glide_slope, col = "#118ab2", lty = 2)
abline(a = femur_DE_pANCOVA_size_ecotype_ground_int, b = femur_DE_pANCOVA_size_ecotype_ground_slope, col = "#ffd166", lty = 2)
abline(a = femur_DE_pANCOVA_size_ecotype_tree_int, b = femur_DE_pANCOVA_size_ecotype_tree_slope, col = "#06d6a0")


######femur_DE_slope_size_meandiff_table ###### 
# chip - glide
femur_DE_slope_size_chip_glide_SE <- sqrt((sd(femur_DE_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)))

femur_DE_slope_size_chip_glide_meandiff <- mean(femur_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)

femur_DE_slope_size_chip_glide_CI <- mean(femur_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_glide_slope_boot) + c(-1, 1) * 1.96 * femur_DE_slope_size_chip_glide_SE  # 95% confidence interval using z=1.96

femur_DE_slope_size_chip_glide_table <- cbind(femur_DE_slope_size_chip_glide_meandiff, femur_DE_slope_size_chip_glide_CI[1], femur_DE_slope_size_chip_glide_CI[2])
rownames(femur_DE_slope_size_chip_glide_table) <- ("femur_DE_slope_size_chip_glide")
colnames(femur_DE_slope_size_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_DE_slope_size_chip_glide_table

# chip - ground
femur_DE_slope_size_chip_ground_SE <- sqrt((sd(femur_DE_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)))

femur_DE_slope_size_chip_ground_meandiff <- mean(femur_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)

femur_DE_slope_size_chip_ground_CI <- mean(femur_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * femur_DE_slope_size_chip_ground_SE  # 95% confidence interval using z=1.96

femur_DE_slope_size_chip_ground_table <- cbind(femur_DE_slope_size_chip_ground_meandiff, femur_DE_slope_size_chip_ground_CI[1], femur_DE_slope_size_chip_ground_CI[2])
rownames(femur_DE_slope_size_chip_ground_table) <- ("femur_DE_slope_size_chip_ground")
colnames(femur_DE_slope_size_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_DE_slope_size_chip_ground_table

# chip - tree
femur_DE_slope_size_chip_tree_SE <- sqrt((sd(femur_DE_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)))

femur_DE_slope_size_chip_tree_meandiff <- mean(femur_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)

femur_DE_slope_size_chip_tree_CI <- mean(femur_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_DE_slope_size_chip_tree_SE  # 95% confidence interval using z=1.96

femur_DE_slope_size_chip_tree_table <- cbind(femur_DE_slope_size_chip_tree_meandiff, femur_DE_slope_size_chip_tree_CI[1], femur_DE_slope_size_chip_tree_CI[2])
rownames(femur_DE_slope_size_chip_tree_table) <- ("femur_DE_slope_size_chip_tree")
colnames(femur_DE_slope_size_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_DE_slope_size_chip_tree_table

# glide - ground
femur_DE_slope_size_glide_ground_SE <- sqrt((sd(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)))

femur_DE_slope_size_glide_ground_meandiff <- mean(femur_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)

femur_DE_slope_size_glide_ground_CI <- mean(femur_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * femur_DE_slope_size_glide_ground_SE  

femur_DE_slope_size_glide_ground_table <- cbind(femur_DE_slope_size_glide_ground_meandiff, femur_DE_slope_size_glide_ground_CI[1], femur_DE_slope_size_glide_ground_CI[2])
rownames(femur_DE_slope_size_glide_ground_table) <- ("femur_DE_slope_size_glide_ground")
colnames(femur_DE_slope_size_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_DE_slope_size_glide_ground_table

# glide - tree
femur_DE_slope_size_glide_tree_SE <- sqrt((sd(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)))

femur_DE_slope_size_glide_tree_meandiff <- mean(femur_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)

femur_DE_slope_size_glide_tree_CI <- mean(femur_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_DE_slope_size_glide_tree_SE  

femur_DE_slope_size_glide_tree_table <- cbind(femur_DE_slope_size_glide_tree_meandiff, femur_DE_slope_size_glide_tree_CI[1], femur_DE_slope_size_glide_tree_CI[2])
rownames(femur_DE_slope_size_glide_tree_table) <- ("femur_DE_slope_size_glide_tree")
colnames(femur_DE_slope_size_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_DE_slope_size_glide_tree_table

# ground - tree
femur_DE_slope_size_ground_tree_SE <- sqrt((sd(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_ground_slope_boot)) + (sd(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)))

femur_DE_slope_size_ground_tree_meandiff <- mean(femur_DE_pANCOVA_size_ecotype_ground_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_tree_slope_boot)

femur_DE_slope_size_ground_tree_CI <- mean(femur_DE_pANCOVA_size_ecotype_ground_slope_boot) - mean(femur_DE_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_DE_slope_size_ground_tree_SE 

femur_DE_slope_size_ground_tree_table <- cbind(femur_DE_slope_size_ground_tree_meandiff, femur_DE_slope_size_ground_tree_CI[1], femur_DE_slope_size_ground_tree_CI[2])
rownames(femur_DE_slope_size_ground_tree_table) <- ("femur_DE_slope_size_ground_tree")
colnames(femur_DE_slope_size_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_DE_slope_size_ground_tree_table

#table
femur_DE_slope_size_meandiff_table <- rbind(
  femur_DE_slope_size_chip_glide_table,
  femur_DE_slope_size_chip_ground_table,
  femur_DE_slope_size_chip_tree_table,
  femur_DE_slope_size_glide_ground_table,
  femur_DE_slope_size_glide_tree_table,
  femur_DE_slope_size_ground_tree_table)
round(femur_DE_slope_size_meandiff_table, digits = 2)

##### ecotype model: femur_DE ~ ecotype #####
femur_DE_pANOVA_ecotype <- phylolm(femur_DE~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_DE_pANOVA_ecotype)

# pANOVA
femur_DE_pANOVA_ecotype_BM_rrpp <- lm.rrpp(femur_DE~ecotype, iter = 999, Cov = vcv.phylo(tree_pr), SS.type = "II")
anova(femur_DE_pANOVA_ecotype_BM_rrpp)
femur_DE_pANOVA_ecotype_BM_rrpp_PW <- pairwise(femur_DE_pANOVA_ecotype_BM_rrpp, groups = ecotype )
summary(femur_DE_pANOVA_ecotype_BM_rrpp_PW, show.vectors = TRUE)

#anova
femur_DE_ANOVA_rrpp <- lm.rrpp(femur_DE~ecotype, iter = 999, SS.type = "II")
anova(femur_DE_ANOVA_rrpp)
femur_DE_ANOVA_rrpp_PW <- pairwise(femur_DE_ANOVA_rrpp, groups = ecotype )
summary(femur_DE_ANOVA_rrpp_PW, show.vectors = TRUE)

######  Make table of coefficients (intercept) ######  
femur_DE_pANOVA_ecotype_chip_int <- femur_DE_pANOVA_ecotype$coefficients["(Intercept)"]
femur_DE_pANOVA_ecotype_glide_int <- femur_DE_pANOVA_ecotype$coefficients["(Intercept)"]+femur_DE_pANOVA_ecotype$coefficients["ecotypegliding"]
femur_DE_pANOVA_ecotype_ground_int <- femur_DE_pANOVA_ecotype$coefficients["(Intercept)"]+femur_DE_pANOVA_ecotype$coefficients["ecotypeground"]
femur_DE_pANOVA_ecotype_tree_int <- femur_DE_pANOVA_ecotype$coefficients["(Intercept)"]+femur_DE_pANOVA_ecotype$coefficients["ecotypetree"]
femur_DE_pANOVA_ecotype_chip_int_boot <- femur_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"]
femur_DE_pANOVA_ecotype_glide_int_boot <- femur_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_DE_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
femur_DE_pANOVA_ecotype_ground_int_boot <- femur_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_DE_pANOVA_ecotype$bootstrap[,"ecotypeground"]
femur_DE_pANOVA_ecotype_tree_int_boot <- femur_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_DE_pANOVA_ecotype$bootstrap[,"ecotypetree"]

femur_DE_pANOVA_ecotype_int_cof <- rbind(
  femur_DE_pANOVA_ecotype_chip_int,
  femur_DE_pANOVA_ecotype_glide_int,
  femur_DE_pANOVA_ecotype_ground_int,
  femur_DE_pANOVA_ecotype_tree_int)
femur_DE_pANOVA_ecotype_int_95 <- rbind(
  quantile(femur_DE_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_DE_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

femur_DE_pANOVA_ecotype_cof <- data.frame(femur_DE_pANOVA_ecotype_int_cof, femur_DE_pANOVA_ecotype_int_95)
colnames(femur_DE_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(femur_DE_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(femur_DE_pANOVA_ecotype_cof, digits = 3)

######  Plot femur_DE ~ ecotype ######  
femur_DE_vp_ecology <- ggplot(femur_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= femur_DE, fill = ecotype)) + 
  scale_y_continuous(name="ln femur_DE (mm)", limits = c(11,25), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(femur_DE)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill = "grey") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color = "black", size = 13, angle = 0, hjust = .5, vjust = .5),
        axis.text.y=element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position= "none") 
femur_DE_vp_ecology


######  femur_DE_int_meandiff_table ######  
# chip - glide
femur_DE_int_chip_glide_SE <- sqrt((sd(femur_DE_pANOVA_ecotype_chip_int_boot)^2/length(femur_DE_pANOVA_ecotype_chip_int_boot)) + (sd(femur_DE_pANOVA_ecotype_glide_int_boot)^2/length(femur_DE_pANOVA_ecotype_glide_int_boot)))

femur_DE_int_chip_glide_meandiff <- mean(femur_DE_pANOVA_ecotype_chip_int_boot) - mean(femur_DE_pANOVA_ecotype_glide_int_boot)

femur_DE_int_chip_glide_CI <- mean(femur_DE_pANOVA_ecotype_chip_int_boot) - mean(femur_DE_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * femur_DE_int_chip_glide_SE  # 95% confidence interval using z=1.96

femur_DE_int_chip_glide_table <- cbind(femur_DE_int_chip_glide_meandiff, femur_DE_int_chip_glide_CI[1], femur_DE_int_chip_glide_CI[2])
rownames(femur_DE_int_chip_glide_table) <- ("femur_DE_int_chip_glide")
colnames(femur_DE_int_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_DE_int_chip_glide_table

# chip - ground
femur_DE_int_chip_ground_SE <- sqrt((sd(femur_DE_pANOVA_ecotype_chip_int_boot)^2/length(femur_DE_pANOVA_ecotype_chip_int_boot)) + (sd(femur_DE_pANOVA_ecotype_ground_int_boot)^2/length(femur_DE_pANOVA_ecotype_ground_int_boot)))

femur_DE_int_chip_ground_meandiff <- mean(femur_DE_pANOVA_ecotype_chip_int_boot) - mean(femur_DE_pANOVA_ecotype_ground_int_boot)

femur_DE_int_chip_ground_CI <- mean(femur_DE_pANOVA_ecotype_chip_int_boot) - mean(femur_DE_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_DE_int_chip_ground_SE  # 95% confidence interval using z=1.96

femur_DE_int_chip_ground_table <- cbind(femur_DE_int_chip_ground_meandiff, femur_DE_int_chip_ground_CI[1], femur_DE_int_chip_ground_CI[2])
rownames(femur_DE_int_chip_ground_table) <- ("femur_DE_int_chip_ground")
colnames(femur_DE_int_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_DE_int_chip_ground_table

# chip - tree
femur_DE_int_chip_tree_SE <- sqrt((sd(femur_DE_pANOVA_ecotype_chip_int_boot)^2/length(femur_DE_pANOVA_ecotype_chip_int_boot)) + (sd(femur_DE_pANOVA_ecotype_tree_int_boot)^2/length(femur_DE_pANOVA_ecotype_tree_int_boot)))

femur_DE_int_chip_tree_meandiff <- mean(femur_DE_pANOVA_ecotype_chip_int_boot) - mean(femur_DE_pANOVA_ecotype_tree_int_boot)

femur_DE_int_chip_tree_CI <- mean(femur_DE_pANOVA_ecotype_chip_int_boot) - mean(femur_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_DE_int_chip_tree_SE  # 95% confidence interval using z=1.96

femur_DE_int_chip_tree_table <- cbind(femur_DE_int_chip_tree_meandiff, femur_DE_int_chip_tree_CI[1], femur_DE_int_chip_tree_CI[2])
rownames(femur_DE_int_chip_tree_table) <- ("femur_DE_int_chip_tree")
colnames(femur_DE_int_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_DE_int_chip_tree_table

# glide - ground
femur_DE_int_glide_ground_SE <- sqrt((sd(femur_DE_pANOVA_ecotype_glide_int_boot)^2/length(femur_DE_pANOVA_ecotype_glide_int_boot)) + (sd(femur_DE_pANOVA_ecotype_ground_int_boot)^2/length(femur_DE_pANOVA_ecotype_ground_int_boot)))

femur_DE_int_glide_ground_meandiff <- mean(femur_DE_pANOVA_ecotype_glide_int_boot) - mean(femur_DE_pANOVA_ecotype_ground_int_boot)

femur_DE_int_glide_ground_CI <- mean(femur_DE_pANOVA_ecotype_glide_int_boot) - mean(femur_DE_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_DE_int_glide_ground_SE  # 95% confidence interval using z=1.96

femur_DE_int_glide_ground_table <- cbind(femur_DE_int_glide_ground_meandiff, femur_DE_int_glide_ground_CI[1], femur_DE_int_glide_ground_CI[2])
rownames(femur_DE_int_glide_ground_table) <- ("femur_DE_int_glide_ground")
colnames(femur_DE_int_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_DE_int_glide_ground_table

# glide - tree
femur_DE_int_glide_tree_SE <- sqrt((sd(femur_DE_pANOVA_ecotype_glide_int_boot)^2/length(femur_DE_pANOVA_ecotype_glide_int_boot)) + (sd(femur_DE_pANOVA_ecotype_tree_int_boot)^2/length(femur_DE_pANOVA_ecotype_tree_int_boot)))

femur_DE_int_glide_tree_meandiff <- mean(femur_DE_pANOVA_ecotype_glide_int_boot) - mean(femur_DE_pANOVA_ecotype_tree_int_boot)

femur_DE_int_glide_tree_CI <- mean(femur_DE_pANOVA_ecotype_glide_int_boot) - mean(femur_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_DE_int_glide_tree_SE  # 95% confidence interval using z=1.96

femur_DE_int_glide_tree_table <- cbind(femur_DE_int_glide_tree_meandiff, femur_DE_int_glide_tree_CI[1], femur_DE_int_glide_tree_CI[2])
rownames(femur_DE_int_glide_tree_table) <- ("femur_DE_int_glide_tree")
colnames(femur_DE_int_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_DE_int_glide_tree_table

# ground - tree
femur_DE_int_ground_tree_SE <- sqrt((sd(femur_DE_pANOVA_ecotype_ground_int_boot)^2/length(femur_DE_pANOVA_ecotype_ground_int_boot)) + (sd(femur_DE_pANOVA_ecotype_tree_int_boot)^2/length(femur_DE_pANOVA_ecotype_tree_int_boot)))

femur_DE_int_ground_tree_meandiff <- mean(femur_DE_pANOVA_ecotype_ground_int_boot) - mean(femur_DE_pANOVA_ecotype_tree_int_boot)

femur_DE_int_ground_tree_CI <- mean(femur_DE_pANOVA_ecotype_ground_int_boot) - mean(femur_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_DE_int_ground_tree_SE  # 95% confidence interval using z=1.96

femur_DE_int_ground_tree_table <- cbind(femur_DE_int_ground_tree_meandiff, femur_DE_int_ground_tree_CI[1], femur_DE_int_ground_tree_CI[2])
rownames(femur_DE_int_ground_tree_table) <- ("femur_DE_int_ground_tree")
colnames(femur_DE_int_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_DE_int_ground_tree_table

#table
femur_DE_int_meandiff_table <- rbind(
  femur_DE_int_chip_glide_table,
  femur_DE_int_chip_ground_table,
  femur_DE_int_chip_tree_table,
  femur_DE_int_glide_ground_table,
  femur_DE_int_glide_tree_table,
  femur_DE_int_ground_tree_table)
round(femur_DE_int_meandiff_table, digits = 2)


##### femur_DE AIC #####
geiger::aicw(c(AIC(femur_DE_pgls_size),
               AIC(femur_DE_pANOVA_ecotype),
               AIC(femur_DE_pANCOVA_size_ecotype),
               AIC(femur_DE_pANCOVAnoint_size_ecotype),
               AIC(femur_DE_pgls_null)))


#### femur_CSS #####
##### null model: femur_CSS ~ 1 ##### 
femur_CSS_pgls_null <- phylolm(femur_CSS~1, phy = tree_pr, model = "lambda", boot = 100)


##### size model: femur_CSS ~ size ##### 
femur_CSS_pgls_size <- phylolm(femur_CSS~lnfemur_size, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_CSS_pgls_size)


##### size+ecotype model: femur_CSS ~ size+ecotype ##### 
femur_CSS_pANCOVAnoint_size_ecotype <- phylolm(femur_CSS~lnfemur_size+ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_CSS_pANCOVAnoint_size_ecotype)


##### size*ecotype model: femur_CSS ~ size*ecotype ##### 
femur_CSS_pANCOVA_size_ecotype <- phylolm(femur_CSS~lnfemur_size*ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_CSS_pANCOVA_size_ecotype)


######  Make table of coefficients (intercept and slope) ######  
femur_CSS_pgls_size_all_int <- femur_CSS_pgls_size$coefficients["(Intercept)"]
femur_CSS_pgls_size_all_int_boot <- femur_CSS_pgls_size$bootstrap[,"(Intercept)"]

femur_CSS_pANCOVA_size_ecotype_chip_int <- femur_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]
femur_CSS_pANCOVA_size_ecotype_glide_int <- femur_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_CSS_pANCOVA_size_ecotype$coefficients["ecotypegliding"]
femur_CSS_pANCOVA_size_ecotype_ground_int <- femur_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_CSS_pANCOVA_size_ecotype$coefficients["ecotypeground"]
femur_CSS_pANCOVA_size_ecotype_tree_int <- femur_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]+femur_CSS_pANCOVA_size_ecotype$coefficients["ecotypetree"]
femur_CSS_pANCOVA_size_ecotype_chip_int_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"]
femur_CSS_pANCOVA_size_ecotype_glide_int_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_CSS_pANCOVA_size_ecotype$bootstrap[,"ecotypegliding"]
femur_CSS_pANCOVA_size_ecotype_ground_int_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_CSS_pANCOVA_size_ecotype$bootstrap[,"ecotypeground"]
femur_CSS_pANCOVA_size_ecotype_tree_int_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + femur_CSS_pANCOVA_size_ecotype$bootstrap[,"ecotypetree"]

femur_CSS_pgls_size_all_slope <- femur_CSS_pgls_size$coefficients["lnfemur_size"]
femur_CSS_pgls_size_all_slope_boot <- femur_CSS_pgls_size$bootstrap[,"lnfemur_size"]

femur_CSS_pANCOVA_size_ecotype_chip_slope <- femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size"]
femur_CSS_pANCOVA_size_ecotype_glide_slope <- femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypegliding"]
femur_CSS_pANCOVA_size_ecotype_ground_slope <- femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypeground"]
femur_CSS_pANCOVA_size_ecotype_tree_slope <- femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size"] + femur_CSS_pANCOVA_size_ecotype$coefficients["lnfemur_size:ecotypetree"]
femur_CSS_pANCOVA_size_ecotype_chip_slope_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"]
femur_CSS_pANCOVA_size_ecotype_glide_slope_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypegliding"]
femur_CSS_pANCOVA_size_ecotype_ground_slope_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypeground"]
femur_CSS_pANCOVA_size_ecotype_tree_slope_boot <- femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size"] + femur_CSS_pANCOVA_size_ecotype$bootstrap[,"lnfemur_size:ecotypetree"]

femur_CSS_pANCOVA_size_ecotype_int_cof <- rbind(
  femur_CSS_pgls_size_all_int,
  femur_CSS_pANCOVA_size_ecotype_chip_int,
  femur_CSS_pANCOVA_size_ecotype_glide_int,
  femur_CSS_pANCOVA_size_ecotype_ground_int,
  femur_CSS_pANCOVA_size_ecotype_tree_int)
femur_CSS_pANCOVA_size_ecotype_int_95 <- rbind(
  quantile(femur_CSS_pgls_size_all_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_tree_int_boot, prob = c(0.025, 0.975)))
femur_CSS_pANCOVA_size_ecotype_slope_cof <- rbind(
  femur_CSS_pgls_size_all_slope,
  femur_CSS_pANCOVA_size_ecotype_chip_slope,
  femur_CSS_pANCOVA_size_ecotype_glide_slope,
  femur_CSS_pANCOVA_size_ecotype_ground_slope,
  femur_CSS_pANCOVA_size_ecotype_tree_slope)
femur_CSS_pANCOVA_size_ecotype_slope_95 <- rbind(
  quantile(femur_CSS_pgls_size_all_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot, prob = c(0.025, 0.975)))

femur_CSS_pANCOVA_size_ecotype_cof <- data.frame(femur_CSS_pANCOVA_size_ecotype_int_cof, femur_CSS_pANCOVA_size_ecotype_int_95, femur_CSS_pANCOVA_size_ecotype_slope_cof, femur_CSS_pANCOVA_size_ecotype_slope_95 )
colnames(femur_CSS_pANCOVA_size_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95", "Slope", "Slope_L95", "Slope_U95")
rownames(femur_CSS_pANCOVA_size_ecotype_cof) <- c("all", "chip", "gliding", "ground", "tree")
round(femur_CSS_pANCOVA_size_ecotype_cof, digits = 3)

######  plot femur_CSS ~ lnsize*ecology ######  
plot(femur_CSS ~ lnfemur_size, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)
abline(a = femur_CSS_pgls_size_all_int, b = femur_CSS_pgls_size_all_slope, col = "black")
abline(a = femur_CSS_pANCOVA_size_ecotype_chip_int, b = femur_CSS_pANCOVA_size_ecotype_chip_slope, col = "#073b4c", lty = 2)
abline(a = femur_CSS_pANCOVA_size_ecotype_glide_int, b = femur_CSS_pANCOVA_size_ecotype_glide_slope, col = "#118ab2", lty = 2)
abline(a = femur_CSS_pANCOVA_size_ecotype_ground_int, b = femur_CSS_pANCOVA_size_ecotype_ground_slope, col = "#ffd166")
abline(a = femur_CSS_pANCOVA_size_ecotype_tree_int, b = femur_CSS_pANCOVA_size_ecotype_tree_slope, col = "#06d6a0", lty = 2)

######  femur_CSS_slope_size_meandiff_table######  
# chip - glide
femur_CSS_slope_size_chip_glide_SE <- sqrt((sd(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)))

femur_CSS_slope_size_chip_glide_meandiff <- mean(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)

femur_CSS_slope_size_chip_glide_CI <- mean(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot) + c(-1, 1) * 1.96 * femur_CSS_slope_size_chip_glide_SE  # 95% confidence slopeerval using z=1.96

femur_CSS_slope_size_chip_glide_table <- cbind(femur_CSS_slope_size_chip_glide_meandiff, femur_CSS_slope_size_chip_glide_CI[1], femur_CSS_slope_size_chip_glide_CI[2])
rownames(femur_CSS_slope_size_chip_glide_table) <- ("femur_CSS_slope_size_chip_glide")
colnames(femur_CSS_slope_size_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_CSS_slope_size_chip_glide_table

# chip - ground
femur_CSS_slope_size_chip_ground_SE <- sqrt((sd(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)))

femur_CSS_slope_size_chip_ground_meandiff <- mean(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)

femur_CSS_slope_size_chip_ground_CI <- mean(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * femur_CSS_slope_size_chip_ground_SE  # 95% confidence slopeerval using z=1.96

femur_CSS_slope_size_chip_ground_table <- cbind(femur_CSS_slope_size_chip_ground_meandiff, femur_CSS_slope_size_chip_ground_CI[1], femur_CSS_slope_size_chip_ground_CI[2])
rownames(femur_CSS_slope_size_chip_ground_table) <- ("femur_CSS_slope_size_chip_ground")
colnames(femur_CSS_slope_size_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_CSS_slope_size_chip_ground_table

# chip - tree
femur_CSS_slope_size_chip_tree_SE <- sqrt((sd(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)))

femur_CSS_slope_size_chip_tree_meandiff <- mean(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)

femur_CSS_slope_size_chip_tree_CI <- mean(femur_CSS_pANCOVA_size_ecotype_chip_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_CSS_slope_size_chip_tree_SE  # 95% confidence slopeerval using z=1.96

femur_CSS_slope_size_chip_tree_table <- cbind(femur_CSS_slope_size_chip_tree_meandiff, femur_CSS_slope_size_chip_tree_CI[1], femur_CSS_slope_size_chip_tree_CI[2])
rownames(femur_CSS_slope_size_chip_tree_table) <- ("femur_CSS_slope_size_chip_tree")
colnames(femur_CSS_slope_size_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_CSS_slope_size_chip_tree_table

# glide - ground
femur_CSS_slope_size_glide_ground_SE <- sqrt((sd(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)))

femur_CSS_slope_size_glide_ground_meandiff <- mean(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)

femur_CSS_slope_size_glide_ground_CI <- mean(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * femur_CSS_slope_size_glide_ground_SE  # 95% confidence slopeerval using z=1.96

femur_CSS_slope_size_glide_ground_table <- cbind(femur_CSS_slope_size_glide_ground_meandiff, femur_CSS_slope_size_glide_ground_CI[1], femur_CSS_slope_size_glide_ground_CI[2])
rownames(femur_CSS_slope_size_glide_ground_table) <- ("femur_CSS_slope_size_glide_ground")
colnames(femur_CSS_slope_size_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_CSS_slope_size_glide_ground_table

# glide - tree
femur_CSS_slope_size_glide_tree_SE <- sqrt((sd(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)))

femur_CSS_slope_size_glide_tree_meandiff <- mean(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)

femur_CSS_slope_size_glide_tree_CI <- mean(femur_CSS_pANCOVA_size_ecotype_glide_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_CSS_slope_size_glide_tree_SE  # 95% confidence slopeerval using z=1.96

femur_CSS_slope_size_glide_tree_table <- cbind(femur_CSS_slope_size_glide_tree_meandiff, femur_CSS_slope_size_glide_tree_CI[1], femur_CSS_slope_size_glide_tree_CI[2])
rownames(femur_CSS_slope_size_glide_tree_table) <- ("femur_CSS_slope_size_glide_tree")
colnames(femur_CSS_slope_size_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_CSS_slope_size_glide_tree_table

# ground - tree
femur_CSS_slope_size_ground_tree_SE <- sqrt((sd(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot)) + (sd(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)^2/length(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)))

femur_CSS_slope_size_ground_tree_meandiff <- mean(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot)

femur_CSS_slope_size_ground_tree_CI <- mean(femur_CSS_pANCOVA_size_ecotype_ground_slope_boot) - mean(femur_CSS_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * femur_CSS_slope_size_ground_tree_SE  # 95% confidence slopeerval using z=1.96

femur_CSS_slope_size_ground_tree_table <- cbind(femur_CSS_slope_size_ground_tree_meandiff, femur_CSS_slope_size_ground_tree_CI[1], femur_CSS_slope_size_ground_tree_CI[2])
rownames(femur_CSS_slope_size_ground_tree_table) <- ("femur_CSS_slope_size_ground_tree")
colnames(femur_CSS_slope_size_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_CSS_slope_size_ground_tree_table

#table
femur_CSS_slope_size_meandiff_table <- rbind(
  femur_CSS_slope_size_chip_glide_table,
  femur_CSS_slope_size_chip_ground_table,
  femur_CSS_slope_size_chip_tree_table,
  femur_CSS_slope_size_glide_ground_table,
  femur_CSS_slope_size_glide_tree_table,
  femur_CSS_slope_size_ground_tree_table)
round(femur_CSS_slope_size_meandiff_table, digits = 2)


##### ecotype model: femur_CSS ~ ecotype #####
femur_CSS_pANOVA_ecotype <- phylolm(femur_CSS~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_CSS_pANOVA_ecotype)

#pANOVA
femur_CSS_pANOVA_ecotype_BM_rrpp <- lm.rrpp(femur_CSS~ecotype, iter = 999, Cov = vcv.phylo(tree_pr), SS.type = "II")
anova(femur_CSS_pANOVA_ecotype_BM_rrpp)
femur_CSS_pANOVA_ecotype_BM_rrpp_PW <- pairwise(femur_CSS_pANOVA_ecotype_BM_rrpp, groups = ecotype )
summary(femur_CSS_pANOVA_ecotype_BM_rrpp_PW, show.vectors = TRUE)

#ANOVA
femur_CSS_ANOVA_rrpp <- lm.rrpp(femur_CSS~ecotype, iter = 999, SS.type = "II")
anova(femur_CSS_ANOVA_rrpp)
femur_CSS_ANOVA_rrpp_PW <- pairwise(femur_CSS_ANOVA_rrpp, groups = ecotype )
summary(femur_CSS_ANOVA_rrpp_PW, show.vectors = TRUE)

###### Make table of coefficients (intercept) ###### 
femur_CSS_pANOVA_ecotype_chip_int <- femur_CSS_pANOVA_ecotype$coefficients["(Intercept)"]
femur_CSS_pANOVA_ecotype_glide_int <- femur_CSS_pANOVA_ecotype$coefficients["(Intercept)"]+femur_CSS_pANOVA_ecotype$coefficients["ecotypegliding"]
femur_CSS_pANOVA_ecotype_ground_int <- femur_CSS_pANOVA_ecotype$coefficients["(Intercept)"]+femur_CSS_pANOVA_ecotype$coefficients["ecotypeground"]
femur_CSS_pANOVA_ecotype_tree_int <- femur_CSS_pANOVA_ecotype$coefficients["(Intercept)"]+femur_CSS_pANOVA_ecotype$coefficients["ecotypetree"]
femur_CSS_pANOVA_ecotype_chip_int_boot <- femur_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"]
femur_CSS_pANOVA_ecotype_glide_int_boot <- femur_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_CSS_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
femur_CSS_pANOVA_ecotype_ground_int_boot <- femur_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_CSS_pANOVA_ecotype$bootstrap[,"ecotypeground"]
femur_CSS_pANOVA_ecotype_tree_int_boot <- femur_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_CSS_pANOVA_ecotype$bootstrap[,"ecotypetree"]

femur_CSS_pANOVA_ecotype_int_cof <- rbind(
  femur_CSS_pANOVA_ecotype_chip_int,
  femur_CSS_pANOVA_ecotype_glide_int,
  femur_CSS_pANOVA_ecotype_ground_int,
  femur_CSS_pANOVA_ecotype_tree_int)
femur_CSS_pANOVA_ecotype_int_95 <- rbind(
  quantile(femur_CSS_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_CSS_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

femur_CSS_pANOVA_ecotype_cof <- data.frame(femur_CSS_pANOVA_ecotype_int_cof, femur_CSS_pANOVA_ecotype_int_95)
colnames(femur_CSS_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(femur_CSS_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(femur_CSS_pANOVA_ecotype_cof, digits = 3)

######  Plot femur_CSS ~ ecotype ###### 
femur_CSS_vp_ecology <- ggplot(femur_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= femur_CSS, fill = ecotype)) + 
  scale_y_continuous(name="ln femur_CSS (mm)", limits = c(1.05,2.6), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(femur_CSS)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill = "grey") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color = "black", size = 13, angle = 0, hjust = .5, vjust = .5),
        axis.text.y=element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position= "none") 
femur_CSS_vp_ecology


###### femur_CSS_int_meandiff_table###### 
# chip - glide
femur_CSS_int_chip_glide_SE <- sqrt((sd(femur_CSS_pANOVA_ecotype_chip_int_boot)^2/length(femur_CSS_pANOVA_ecotype_chip_int_boot)) + (sd(femur_CSS_pANOVA_ecotype_glide_int_boot)^2/length(femur_CSS_pANOVA_ecotype_glide_int_boot)))

femur_CSS_int_chip_glide_meandiff <- mean(femur_CSS_pANOVA_ecotype_chip_int_boot) - mean(femur_CSS_pANOVA_ecotype_glide_int_boot)

femur_CSS_int_chip_glide_CI <- mean(femur_CSS_pANOVA_ecotype_chip_int_boot) - mean(femur_CSS_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * femur_CSS_int_chip_glide_SE  # 95% confidence interval using z=1.96

femur_CSS_int_chip_glide_table <- cbind(femur_CSS_int_chip_glide_meandiff, femur_CSS_int_chip_glide_CI[1], femur_CSS_int_chip_glide_CI[2])
rownames(femur_CSS_int_chip_glide_table) <- ("femur_CSS_int_chip_glide")
colnames(femur_CSS_int_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_CSS_int_chip_glide_table

# chip - ground
femur_CSS_int_chip_ground_SE <- sqrt((sd(femur_CSS_pANOVA_ecotype_chip_int_boot)^2/length(femur_CSS_pANOVA_ecotype_chip_int_boot)) + (sd(femur_CSS_pANOVA_ecotype_ground_int_boot)^2/length(femur_CSS_pANOVA_ecotype_ground_int_boot)))

femur_CSS_int_chip_ground_meandiff <- mean(femur_CSS_pANOVA_ecotype_chip_int_boot) - mean(femur_CSS_pANOVA_ecotype_ground_int_boot)

femur_CSS_int_chip_ground_CI <- mean(femur_CSS_pANOVA_ecotype_chip_int_boot) - mean(femur_CSS_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_CSS_int_chip_ground_SE  # 95% confidence interval using z=1.96

femur_CSS_int_chip_ground_table <- cbind(femur_CSS_int_chip_ground_meandiff, femur_CSS_int_chip_ground_CI[1], femur_CSS_int_chip_ground_CI[2])
rownames(femur_CSS_int_chip_ground_table) <- ("femur_CSS_int_chip_ground")
colnames(femur_CSS_int_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_CSS_int_chip_ground_table

# chip - tree
femur_CSS_int_chip_tree_SE <- sqrt((sd(femur_CSS_pANOVA_ecotype_chip_int_boot)^2/length(femur_CSS_pANOVA_ecotype_chip_int_boot)) + (sd(femur_CSS_pANOVA_ecotype_tree_int_boot)^2/length(femur_CSS_pANOVA_ecotype_tree_int_boot)))

femur_CSS_int_chip_tree_meandiff <- mean(femur_CSS_pANOVA_ecotype_chip_int_boot) - mean(femur_CSS_pANOVA_ecotype_tree_int_boot)

femur_CSS_int_chip_tree_CI <- mean(femur_CSS_pANOVA_ecotype_chip_int_boot) - mean(femur_CSS_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_CSS_int_chip_tree_SE  

femur_CSS_int_chip_tree_table <- cbind(femur_CSS_int_chip_tree_meandiff, femur_CSS_int_chip_tree_CI[1], femur_CSS_int_chip_tree_CI[2])
rownames(femur_CSS_int_chip_tree_table) <- ("femur_CSS_int_chip_tree")
colnames(femur_CSS_int_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_CSS_int_chip_tree_table

# glide - ground
femur_CSS_int_glide_ground_SE <- sqrt((sd(femur_CSS_pANOVA_ecotype_glide_int_boot)^2/length(femur_CSS_pANOVA_ecotype_glide_int_boot)) + (sd(femur_CSS_pANOVA_ecotype_ground_int_boot)^2/length(femur_CSS_pANOVA_ecotype_ground_int_boot)))

femur_CSS_int_glide_ground_meandiff <- mean(femur_CSS_pANOVA_ecotype_glide_int_boot) - mean(femur_CSS_pANOVA_ecotype_ground_int_boot)

femur_CSS_int_glide_ground_CI <- mean(femur_CSS_pANOVA_ecotype_glide_int_boot) - mean(femur_CSS_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_CSS_int_glide_ground_SE

femur_CSS_int_glide_ground_table <- cbind(femur_CSS_int_glide_ground_meandiff, femur_CSS_int_glide_ground_CI[1], femur_CSS_int_glide_ground_CI[2])
rownames(femur_CSS_int_glide_ground_table) <- ("femur_CSS_int_glide_ground")
colnames(femur_CSS_int_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_CSS_int_glide_ground_table

# glide - tree
femur_CSS_int_glide_tree_SE <- sqrt((sd(femur_CSS_pANOVA_ecotype_glide_int_boot)^2/length(femur_CSS_pANOVA_ecotype_glide_int_boot)) + (sd(femur_CSS_pANOVA_ecotype_tree_int_boot)^2/length(femur_CSS_pANOVA_ecotype_tree_int_boot)))

femur_CSS_int_glide_tree_meandiff <- mean(femur_CSS_pANOVA_ecotype_glide_int_boot) - mean(femur_CSS_pANOVA_ecotype_tree_int_boot)

femur_CSS_int_glide_tree_CI <- mean(femur_CSS_pANOVA_ecotype_glide_int_boot) - mean(femur_CSS_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_CSS_int_glide_tree_SE  

femur_CSS_int_glide_tree_table <- cbind(femur_CSS_int_glide_tree_meandiff, femur_CSS_int_glide_tree_CI[1], femur_CSS_int_glide_tree_CI[2])
rownames(femur_CSS_int_glide_tree_table) <- ("femur_CSS_int_glide_tree")
colnames(femur_CSS_int_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_CSS_int_glide_tree_table

# ground - tree
femur_CSS_int_ground_tree_SE <- sqrt((sd(femur_CSS_pANOVA_ecotype_ground_int_boot)^2/length(femur_CSS_pANOVA_ecotype_ground_int_boot)) + (sd(femur_CSS_pANOVA_ecotype_tree_int_boot)^2/length(femur_CSS_pANOVA_ecotype_tree_int_boot)))

femur_CSS_int_ground_tree_meandiff <- mean(femur_CSS_pANOVA_ecotype_ground_int_boot) - mean(femur_CSS_pANOVA_ecotype_tree_int_boot)

femur_CSS_int_ground_tree_CI <- mean(femur_CSS_pANOVA_ecotype_ground_int_boot) - mean(femur_CSS_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_CSS_int_ground_tree_SE  

femur_CSS_int_ground_tree_table <- cbind(femur_CSS_int_ground_tree_meandiff, femur_CSS_int_ground_tree_CI[1], femur_CSS_int_ground_tree_CI[2])
rownames(femur_CSS_int_ground_tree_table) <- ("femur_CSS_int_ground_tree")
colnames(femur_CSS_int_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_CSS_int_ground_tree_table

#table
femur_CSS_int_meandiff_table <- rbind(
  femur_CSS_int_chip_glide_table,
  femur_CSS_int_chip_ground_table,
  femur_CSS_int_chip_tree_table,
  femur_CSS_int_glide_ground_table,
  femur_CSS_int_glide_tree_table,
  femur_CSS_int_ground_tree_table)
round(femur_CSS_int_meandiff_table, digits = 2)

##### femur_CSS AIC ####
geiger::aicw(c(AIC(femur_CSS_pgls_size),
               AIC(femur_CSS_pANOVA_ecotype),
               AIC(femur_CSS_pANCOVA_size_ecotype),
               AIC(femur_CSS_pANCOVAnoint_size_ecotype),
               AIC(femur_CSS_pgls_null)))


#### femur size ####
### femur_size ANCOVA ###
femur_size_pANOVA_ecotype <- phylolm(lnfemur_size~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(femur_size_pANOVA_ecotype)

femur_size_pANOVA_ecotype_chip_int <- femur_size_pANOVA_ecotype$coefficients["(Intercept)"]
femur_size_pANOVA_ecotype_glide_int <- femur_size_pANOVA_ecotype$coefficients["(Intercept)"]+femur_size_pANOVA_ecotype$coefficients["ecotypegliding"]
femur_size_pANOVA_ecotype_ground_int <- femur_size_pANOVA_ecotype$coefficients["(Intercept)"]+femur_size_pANOVA_ecotype$coefficients["ecotypeground"]
femur_size_pANOVA_ecotype_tree_int <- femur_size_pANOVA_ecotype$coefficients["(Intercept)"]+femur_size_pANOVA_ecotype$coefficients["ecotypetree"]
femur_size_pANOVA_ecotype_chip_int_boot <- femur_size_pANOVA_ecotype$bootstrap[,"(Intercept)"]
femur_size_pANOVA_ecotype_glide_int_boot <- femur_size_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_size_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
femur_size_pANOVA_ecotype_ground_int_boot <- femur_size_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_size_pANOVA_ecotype$bootstrap[,"ecotypeground"]
femur_size_pANOVA_ecotype_tree_int_boot <- femur_size_pANOVA_ecotype$bootstrap[,"(Intercept)"] + femur_size_pANOVA_ecotype$bootstrap[,"ecotypetree"]

femur_size_pANOVA_ecotype_int_cof <- rbind(
  femur_size_pANOVA_ecotype_chip_int,
  femur_size_pANOVA_ecotype_glide_int,
  femur_size_pANOVA_ecotype_ground_int,
  femur_size_pANOVA_ecotype_tree_int)
femur_size_pANOVA_ecotype_int_95 <- rbind(
  quantile(femur_size_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_size_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_size_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(femur_size_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

femur_size_pANOVA_ecotype_cof <- data.frame(femur_size_pANOVA_ecotype_int_cof, femur_size_pANOVA_ecotype_int_95)
colnames(femur_size_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(femur_size_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(femur_size_pANOVA_ecotype_cof, digits = 1)

femur_size_vp_ecology <- ggplot(femur_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= lnfemur_size, fill = ecotype)) + 
  scale_y_continuous(name="ln femur_size (mm)", limits = c(4.2,6.1), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(lnfemur_size)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill = "grey") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x=element_text(color = "black", size = 13, angle = 0, hjust = .5, vjust = .5),
        axis.text.y=element_text(color = "black", size = 14),
        axis.title.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.background = element_blank(),
        legend.position= "none") 
femur_size_vp_ecology


## femur_size_int_meandiff_table ##
# chip - glide
femur_size_int_chip_glide_SE <- sqrt((sd(femur_size_pANOVA_ecotype_chip_int_boot)^2/length(femur_size_pANOVA_ecotype_chip_int_boot)) + (sd(femur_size_pANOVA_ecotype_glide_int_boot)^2/length(femur_size_pANOVA_ecotype_glide_int_boot)))

femur_size_int_chip_glide_meandiff <- mean(femur_size_pANOVA_ecotype_chip_int_boot) - mean(femur_size_pANOVA_ecotype_glide_int_boot)

femur_size_int_chip_glide_CI <- mean(femur_size_pANOVA_ecotype_chip_int_boot) - mean(femur_size_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * femur_size_int_chip_glide_SE  # 95% confidence interval using z=1.96

femur_size_int_chip_glide_table <- cbind(femur_size_int_chip_glide_meandiff, femur_size_int_chip_glide_CI[1], femur_size_int_chip_glide_CI[2])
rownames(femur_size_int_chip_glide_table) <- ("femur_size_int_chip_glide")
colnames(femur_size_int_chip_glide_table) <- c("meandiff", "L95", "U95")
femur_size_int_chip_glide_table

# chip - ground
femur_size_int_chip_ground_SE <- sqrt((sd(femur_size_pANOVA_ecotype_chip_int_boot)^2/length(femur_size_pANOVA_ecotype_chip_int_boot)) + (sd(femur_size_pANOVA_ecotype_ground_int_boot)^2/length(femur_size_pANOVA_ecotype_ground_int_boot)))

femur_size_int_chip_ground_meandiff <- mean(femur_size_pANOVA_ecotype_chip_int_boot) - mean(femur_size_pANOVA_ecotype_ground_int_boot)

femur_size_int_chip_ground_CI <- mean(femur_size_pANOVA_ecotype_chip_int_boot) - mean(femur_size_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_size_int_chip_ground_SE  # 95% confidence interval using z=1.96

femur_size_int_chip_ground_table <- cbind(femur_size_int_chip_ground_meandiff, femur_size_int_chip_ground_CI[1], femur_size_int_chip_ground_CI[2])
rownames(femur_size_int_chip_ground_table) <- ("femur_size_int_chip_ground")
colnames(femur_size_int_chip_ground_table) <- c("meandiff", "L95", "U95")
femur_size_int_chip_ground_table

# chip - tree
femur_size_int_chip_tree_SE <- sqrt((sd(femur_size_pANOVA_ecotype_chip_int_boot)^2/length(femur_size_pANOVA_ecotype_chip_int_boot)) + (sd(femur_size_pANOVA_ecotype_tree_int_boot)^2/length(femur_size_pANOVA_ecotype_tree_int_boot)))

femur_size_int_chip_tree_meandiff <- mean(femur_size_pANOVA_ecotype_chip_int_boot) - mean(femur_size_pANOVA_ecotype_tree_int_boot)

femur_size_int_chip_tree_CI <- mean(femur_size_pANOVA_ecotype_chip_int_boot) - mean(femur_size_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_size_int_chip_tree_SE  # 95% confidence interval using z=1.96

femur_size_int_chip_tree_table <- cbind(femur_size_int_chip_tree_meandiff, femur_size_int_chip_tree_CI[1], femur_size_int_chip_tree_CI[2])
rownames(femur_size_int_chip_tree_table) <- ("femur_size_int_chip_tree")
colnames(femur_size_int_chip_tree_table) <- c("meandiff", "L95", "U95")
femur_size_int_chip_tree_table

# glide - ground
femur_size_int_glide_ground_SE <- sqrt((sd(femur_size_pANOVA_ecotype_glide_int_boot)^2/length(femur_size_pANOVA_ecotype_glide_int_boot)) + (sd(femur_size_pANOVA_ecotype_ground_int_boot)^2/length(femur_size_pANOVA_ecotype_ground_int_boot)))

femur_size_int_glide_ground_meandiff <- mean(femur_size_pANOVA_ecotype_glide_int_boot) - mean(femur_size_pANOVA_ecotype_ground_int_boot)

femur_size_int_glide_ground_CI <- mean(femur_size_pANOVA_ecotype_glide_int_boot) - mean(femur_size_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * femur_size_int_glide_ground_SE  # 95% confidence interval using z=1.96

femur_size_int_glide_ground_table <- cbind(femur_size_int_glide_ground_meandiff, femur_size_int_glide_ground_CI[1], femur_size_int_glide_ground_CI[2])
rownames(femur_size_int_glide_ground_table) <- ("femur_size_int_glide_ground")
colnames(femur_size_int_glide_ground_table) <- c("meandiff", "L95", "U95")
femur_size_int_glide_ground_table

# glide - tree
femur_size_int_glide_tree_SE <- sqrt((sd(femur_size_pANOVA_ecotype_glide_int_boot)^2/length(femur_size_pANOVA_ecotype_glide_int_boot)) + (sd(femur_size_pANOVA_ecotype_tree_int_boot)^2/length(femur_size_pANOVA_ecotype_tree_int_boot)))

femur_size_int_glide_tree_meandiff <- mean(femur_size_pANOVA_ecotype_glide_int_boot) - mean(femur_size_pANOVA_ecotype_tree_int_boot)

femur_size_int_glide_tree_CI <- mean(femur_size_pANOVA_ecotype_glide_int_boot) - mean(femur_size_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_size_int_glide_tree_SE  # 95% confidence interval using z=1.96

femur_size_int_glide_tree_table <- cbind(femur_size_int_glide_tree_meandiff, femur_size_int_glide_tree_CI[1], femur_size_int_glide_tree_CI[2])
rownames(femur_size_int_glide_tree_table) <- ("femur_size_int_glide_tree")
colnames(femur_size_int_glide_tree_table) <- c("meandiff", "L95", "U95")
femur_size_int_glide_tree_table

# ground - tree
femur_size_int_ground_tree_SE <- sqrt((sd(femur_size_pANOVA_ecotype_ground_int_boot)^2/length(femur_size_pANOVA_ecotype_ground_int_boot)) + (sd(femur_size_pANOVA_ecotype_tree_int_boot)^2/length(femur_size_pANOVA_ecotype_tree_int_boot)))

femur_size_int_ground_tree_meandiff <- mean(femur_size_pANOVA_ecotype_ground_int_boot) - mean(femur_size_pANOVA_ecotype_tree_int_boot)

femur_size_int_ground_tree_CI <- mean(femur_size_pANOVA_ecotype_ground_int_boot) - mean(femur_size_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * femur_size_int_ground_tree_SE  # 95% confidence interval using z=1.96

femur_size_int_ground_tree_table <- cbind(femur_size_int_ground_tree_meandiff, femur_size_int_ground_tree_CI[1], femur_size_int_ground_tree_CI[2])
rownames(femur_size_int_ground_tree_table) <- ("femur_size_int_ground_tree")
colnames(femur_size_int_ground_tree_table) <- c("meandiff", "L95", "U95")
femur_size_int_ground_tree_table

#table
femur_size_int_meandiff_table <- rbind(
  femur_size_int_chip_glide_table,
  femur_size_int_chip_ground_table,
  femur_size_int_chip_tree_table,
  femur_size_int_glide_ground_table,
  femur_size_int_glide_tree_table,
  femur_size_int_ground_tree_table)
round(femur_size_int_meandiff_table, digits = 2)


#### external shape vs internal structure ####
##### Cg #####
femur_shape_Cg_pls <- two.b.pls(femur_shape, femur_Cg,  iter = 999)
summary(femur_shape_Cg_pls)
plot(femur_shape_Cg_pls, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)

#chipmunk
femur_shape_Cg_pls_chipmunk <- two.b.pls(femur_shape_chipmunk, femur_Cg_chipmunk,  iter = 999)
summary(femur_shape_Cg_pls_chipmunk)
plot(femur_shape_Cg_pls_chipmunk)

#gliding
femur_shape_Cg_pls_gliding <- two.b.pls(femur_shape_gliding, femur_Cg_gliding,  iter = 999)
plot(femur_shape_Cg_pls_gliding)
summary(femur_shape_Cg_pls_gliding)

#ground
femur_shape_Cg_pls_ground <- two.b.pls(femur_shape_ground, femur_Cg_ground,  iter = 999)
plot(femur_shape_Cg_pls_ground)
summary(femur_shape_Cg_pls_ground)

#tree
femur_shape_Cg_pls_tree <- two.b.pls(femur_shape_tree, femur_Cg_tree,  iter = 999)
plot(femur_shape_Cg_pls_tree)
summary(femur_shape_Cg_pls_tree)

femur_shape_Cg_pls_table <- rbind(c(femur_shape_Cg_pls$r.pls, femur_shape_Cg_pls$Z, femur_shape_Cg_pls$P.value),
                                  c(femur_shape_Cg_pls_chipmunk$r.pls, femur_shape_Cg_pls_chipmunk$Z, femur_shape_Cg_pls_chipmunk$P.value),
                                  c(femur_shape_Cg_pls_gliding$r.pls, femur_shape_Cg_pls_gliding$Z, femur_shape_Cg_pls_gliding$P.value),
                                  c(femur_shape_Cg_pls_ground$r.pls, femur_shape_Cg_pls_ground$Z, femur_shape_Cg_pls_ground$P.value),
                                  c(femur_shape_Cg_pls_tree$r.pls, femur_shape_Cg_pls_tree$Z, femur_shape_Cg_pls_tree$P.value))
rownames(femur_shape_Cg_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(femur_shape_Cg_pls_table) <- c("r", "Z", "P")
round(femur_shape_Cg_pls_table, digits = 3)


##### DE #####
femur_shape_DE_pls <- two.b.pls(femur_shape, femur_DE,  iter = 999)
summary(femur_shape_DE_pls)
plot(femur_shape_DE_pls, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)

#chipmunk
femur_shape_DE_pls_chipmunk <- two.b.pls(femur_shape_chipmunk, femur_DE_chipmunk,  iter = 999)
summary(femur_shape_DE_pls_chipmunk)
plot(femur_shape_DE_pls_chipmunk)

#gliding
femur_shape_DE_pls_gliding <- two.b.pls(femur_shape_gliding, femur_DE_gliding,  iter = 999)
plot(femur_shape_DE_pls_gliding)
summary(femur_shape_DE_pls_gliding)

#ground
femur_shape_DE_pls_ground <- two.b.pls(femur_shape_ground, femur_DE_ground,  iter = 999)
plot(femur_shape_DE_pls_ground)
summary(femur_shape_DE_pls_ground)

#tree
femur_shape_DE_pls_tree <- two.b.pls(femur_shape_tree, femur_DE_tree,  iter = 999)
plot(femur_shape_DE_pls_tree)
summary(femur_shape_DE_pls_tree)

femur_shape_DE_pls_table <- rbind(c(femur_shape_DE_pls$r.pls, femur_shape_DE_pls$Z, femur_shape_DE_pls$P.value),
                                  c(femur_shape_DE_pls_chipmunk$r.pls, femur_shape_DE_pls_chipmunk$Z, femur_shape_DE_pls_chipmunk$P.value),
                                  c(femur_shape_DE_pls_gliding$r.pls, femur_shape_DE_pls_gliding$Z, femur_shape_DE_pls_gliding$P.value),
                                  c(femur_shape_DE_pls_ground$r.pls, femur_shape_DE_pls_ground$Z, femur_shape_DE_pls_ground$P.value),
                                  c(femur_shape_DE_pls_tree$r.pls, femur_shape_DE_pls_tree$Z, femur_shape_DE_pls_tree$P.value))
rownames(femur_shape_DE_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(femur_shape_DE_pls_table) <- c("r", "Z", "P")
round(femur_shape_DE_pls_table, digits = 3)


##### CSS #####
femur_shape_CSS_pls <- two.b.pls(femur_shape, femur_CSS,  iter = 999)
summary(femur_shape_CSS_pls)
plot(femur_shape_CSS_pls, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)

#chipmunk
femur_shape_CSS_pls_chipmunk <- two.b.pls(femur_shape_chipmunk, femur_CSS_chipmunk,  iter = 999)
summary(femur_shape_CSS_pls_chipmunk)
plot(femur_shape_CSS_pls_chipmunk)

#gliding
femur_shape_CSS_pls_gliding <- two.b.pls(femur_shape_gliding, femur_CSS_gliding,  iter = 999)
plot(femur_shape_CSS_pls_gliding)
summary(femur_shape_CSS_pls_gliding)

#ground
femur_shape_CSS_pls_ground <- two.b.pls(femur_shape_ground, femur_CSS_ground,  iter = 999)
plot(femur_shape_CSS_pls_ground)
summary(femur_shape_CSS_pls_ground)

#tree
femur_shape_CSS_pls_tree <- two.b.pls(femur_shape_tree, femur_CSS_tree,  iter = 999)
plot(femur_shape_CSS_pls_tree)
summary(femur_shape_CSS_pls_tree)

femur_shape_CSS_pls_table <- rbind(c(femur_shape_CSS_pls$r.pls, femur_shape_CSS_pls$Z, femur_shape_CSS_pls$P.value),
                                   c(femur_shape_CSS_pls_chipmunk$r.pls, femur_shape_CSS_pls_chipmunk$Z, femur_shape_CSS_pls_chipmunk$P.value),
                                   c(femur_shape_CSS_pls_gliding$r.pls, femur_shape_CSS_pls_gliding$Z, femur_shape_CSS_pls_gliding$P.value),
                                   c(femur_shape_CSS_pls_ground$r.pls, femur_shape_CSS_pls_ground$Z, femur_shape_CSS_pls_ground$P.value),
                                   c(femur_shape_CSS_pls_tree$r.pls, femur_shape_CSS_pls_tree$Z, femur_shape_CSS_pls_tree$P.value))
rownames(femur_shape_CSS_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(femur_shape_CSS_pls_table) <- c("r", "Z", "P")
round(femur_shape_CSS_pls_table, digits = 3)

