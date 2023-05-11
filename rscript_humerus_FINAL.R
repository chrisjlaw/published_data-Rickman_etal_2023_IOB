library(phytools)
library(geiger)
library(phylolm)
library(ggplot2)
library(forcats)
library(geomorph)
library(dplyr)
library(mvMORPH)

setwd("~/Documents/Projects/Project_Squirrels/squirrel_limbs")

#### data ####
humerus_data <- read.csv("data/humerus_combined_data.csv")
rownames(humerus_data) <- humerus_data[,1]
tree <- read.nexus("data/squirrel_mcctree.nex")

#remove species in phylogeny that aren't present in dataset
tree_pr <- treedata(phy = tree, data = humerus_data, sort=TRUE)$phy

humerus_data <- humerus_data[tree_pr$tip.label,]

size <- (as.numeric(as.character(humerus_data$geomean))) #save size as vector
names(size) <- rownames(humerus_data) #give species names for size vector
lnsize <- log(size) #natural log transform size vector (for stats purposes)

ecotype <- as.factor(humerus_data$ecotype) #save ecology as factor
names(ecotype) <- rownames(humerus_data) #give species names for ecology factor

col_ecotype <- c("#073b4c", "#118ab2", "#ffd166", "#06d6a0") #set colors for ecological groups: blue = flying; orange = ground; green = tree

#internal structures data 
humerus_Cg <- (as.numeric(as.character(humerus_data$Cg))) #save size as vector
names(humerus_Cg) <- rownames(humerus_data) #give species names for size vector

humerus_DE <- (as.numeric(as.character(humerus_data$DE))) #save size as vector
names(humerus_DE) <- rownames(humerus_data) #give species names for size vector

humerus_CSS <- (as.numeric(as.character(humerus_data$CSS))) #save size as vector
names(humerus_CSS) <- rownames(humerus_data) #give species names for size vector

#external structures data
humerus_rawlandmarks <- humerus_data[,12:593]
humerus_gpa <- gpagen(arrayspecs(humerus_rawlandmarks, 194, 3))
humerus_shape <- humerus_gpa$coords
lnhumerus_size <- log(humerus_gpa$Csize)
humerus_shape_2D <- two.d.array(humerus_shape)
humerus_shape_data <- data.frame(ecotype, lnhumerus_size, humerus_shape_2D)

##### data individual ecotypes #####
#chipmunk data
humerus_data_chipmunk <- filter(humerus_data, ecotype == "chipmunk")
humerus_Cg_chipmunk <- humerus_data_chipmunk$Cg
names(humerus_Cg_chipmunk) <- rownames(humerus_data_chipmunk)
humerus_DE_chipmunk <- humerus_data_chipmunk$DE
names(humerus_DE_chipmunk) <- rownames(humerus_data_chipmunk)
humerus_CSS_chipmunk <- humerus_data_chipmunk$CSS
names(humerus_CSS_chipmunk) <- rownames(humerus_data_chipmunk)
humerus_shape_data_chipmunk <- filter(humerus_shape_data, ecotype == "chipmunk")
humerus_lnsize_chipmunk <- humerus_shape_data_chipmunk$lnhumerus_size
names(humerus_lnsize_chipmunk) <- rownames(humerus_shape_data_chipmunk) 
humerus_shape_chipmunk <- arrayspecs(humerus_shape_data_chipmunk[,3:584], 194,3 )
humerus_shape_chipmunk_2D <- two.d.array(humerus_shape_chipmunk)
lnhumerus_size_chipmunk <- humerus_shape_data_chipmunk$lnhumerus_size
tree_pr_chipmunk <- treedata(phy = tree, data = humerus_shape_data_chipmunk, sort=TRUE)$phy

#gliding data
humerus_data_gliding <- filter(humerus_data, ecotype == "gliding")
humerus_Cg_gliding <- humerus_data_gliding$Cg
names(humerus_Cg_gliding) <- rownames(humerus_data_gliding)
humerus_DE_gliding <- humerus_data_gliding$DE
names(humerus_DE_gliding) <- rownames(humerus_data_gliding)
humerus_CSS_gliding <- humerus_data_gliding$CSS
names(humerus_CSS_gliding) <- rownames(humerus_data_gliding)
humerus_shape_data_gliding <- filter(humerus_shape_data, ecotype == "gliding")
humerus_lnsize_gliding <- humerus_shape_data_gliding$lnhumerus_size
names(humerus_lnsize_gliding) <- rownames(humerus_shape_data_gliding) 
humerus_shape_gliding <- arrayspecs(humerus_shape_data_gliding[,3:584], 194,3 )
humerus_shape_gliding_2D <- two.d.array(humerus_shape_gliding)
lnhumerus_size_gliding <- humerus_shape_data_gliding$lnhumerus_size
tree_pr_gliding <- treedata(phy = tree, data = humerus_shape_data_gliding, sort=TRUE)$phy

#ground data
humerus_data_ground <- filter(humerus_data, ecotype == "ground")
humerus_Cg_ground <- humerus_data_ground$Cg
names(humerus_Cg_ground) <- rownames(humerus_data_ground)
humerus_DE_ground <- humerus_data_ground$DE
names(humerus_DE_ground) <- rownames(humerus_data_ground)
humerus_CSS_ground <- humerus_data_ground$CSS
names(humerus_CSS_ground) <- rownames(humerus_data_ground)
humerus_shape_data_ground <- filter(humerus_shape_data, ecotype == "ground")
humerus_lnsize_ground <- humerus_shape_data_ground$lnhumerus_size
names(humerus_lnsize_ground) <- rownames(humerus_shape_data_ground) 
humerus_shape_ground <- arrayspecs(humerus_shape_data_ground[,3:584], 194,3 )
humerus_shape_ground_2D <- two.d.array(humerus_shape_ground)
lnhumerus_size_ground <- humerus_shape_data_ground$lnhumerus_size
tree_pr_ground <- treedata(phy = tree, data = humerus_shape_data_ground, sort=TRUE)$phy

#tree data
humerus_data_tree <- filter(humerus_data, ecotype == "tree")
humerus_Cg_tree <- humerus_data_tree$Cg
names(humerus_Cg_tree) <- rownames(humerus_data_tree)
humerus_DE_tree <- humerus_data_tree$DE
names(humerus_DE_tree) <- rownames(humerus_data_tree)
humerus_CSS_tree <- humerus_data_tree$CSS
names(humerus_CSS_tree) <- rownames(humerus_data_tree)
humerus_shape_data_tree <- filter(humerus_shape_data, ecotype == "tree")
humerus_lnsize_tree <- humerus_shape_data_tree$lnhumerus_size
names(humerus_lnsize_tree) <- rownames(humerus_shape_data_tree) 
humerus_shape_tree <- arrayspecs(humerus_shape_data_tree[,3:584], 194,3 )
humerus_shape_tree_2D <- two.d.array(humerus_shape_tree)
lnhumerus_size_tree <- humerus_shape_data_tree$lnhumerus_size
tree_pr_tree <- treedata(phy = tree, data = humerus_shape_data_tree, sort=TRUE)$phy


#### humerus_shape #####
##### humeral shape PCA #####
phylogenetic PCA
humerus_shape_pca <- gm.prcomp(humerus_shape, phy = tree_pr, GLS = FALSE)
summary(humerus_shape_pca)
plot(humerus_shape_pca, phylo = TRUE, main = "phylo PCA")

col_ecotype_phylomorph <-c(col_ecotype[as.factor(ecotype)],rep("white",tree_pr$Nnode))
names(col_ecotype_phylomorph)<-1:(length(tree_pr$tip)+tree_pr$Nnode)
phylomorphospace(tree_pr, humerus_shape_pca$x[,1:2], control = list(col.node = col_ecotype_phylomorph), label = "off")

# humerus phy signal
humerus_physig <- physignal(humerus_shape, phy = tree_pr, iter = 999)


##### null model #####
humerus_shape_null_mv <- mvgls(humerus_shape_2D ~ 1, tree = tree_pr, method = "LOOCV", model = "lambda")
save(humerus_shape_null_mv, file= "humerus_shape_null_mv.Rdata")
humerus_shape_null_mv_EIC <- EIC(humerus_shape_null_mv, nboot =1000)
save(humerus_shape_null_mv_EIC, file= "humerus_shape_null_mv_EIC.Rdata")


##### size model: shape ~ size  #####
humerus_shape_allometry_mv <- mvgls(humerus_shape_2D ~ lnhumerus_size, tree = tree_pr, method = "LOOCV", model = "lambda")
save(humerus_shape_allometry_mv, file= "humerus_shape_allometry_mv.Rdata")
humerus_shape_allometry_mv_EIC <- EIC(humerus_shape_allometry_mv, nboot =1000)
save(humerus_shape_allometry_mv_EIC, file= "humerus_shape_allometry_mv_EIC.Rdata")
humerus_shape_allometry_mv_test <- manova.gls(humerus_shape_allometry_mv, type = "II", test = "Pillai")
save(humerus_shape_allometry_mv_test, file= "humerus_shape_allometry_mv_test.Rdata")


##### size*ecotype model: shape ~ ecotype*size  #####
humerus_ecotype_mvpancova <- mvgls(humerus_shape_2D ~ lnhumerus_size*ecotype, tree = tree_pr, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpancova, file= "humerus_ecotype_mvpancova.Rdata")
humerus_ecotype_mvpancova_EIC <- EIC(humerus_ecotype_mvpancova, nboot =1000)
save(humerus_ecotype_mvpancova_EIC, file= "humerus_ecotype_mvpancova_EIC.Rdata")
humerus_ecotype_mvpancova_test <- manova.gls(humerus_ecotype_mvpancova, type = "II", test = "Pillai")
save(humerus_ecotype_mvpancova_test, file= "humerus_ecotype_mvpancova_test.Rdata")


##### size+ecotype model: shape ~ ecotype+size  #####
humerus_ecotype_mvpancova_noint <- mvgls(humerus_shape_2D ~ lnhumerus_size+ecotype, tree = tree_pr, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpancova_noint, file= "humerus_ecotype_mvpancova_noint.Rdata")
humerus_ecotype_mvpancova_noint_EIC <- EIC(humerus_ecotype_mvpancova_noint, nboot =1000)
save(humerus_ecotype_mvpancova_noint_EIC, file= "humerus_ecotype_mvpancova_noint_EIC.Rdata")
humerus_ecotype_mvpancova_noint_test <- manova.gls(humerus_ecotype_mvpancova_noint, type = "II", test = "Pillai")
save(humerus_ecotype_mvpancova_noint_test, file= "humerus_ecotype_mvpancova_noint_test.Rdata")


##### ecotype model: shape ~ ecotype  #####
humerus_ecotype_mvpanova <- mvgls(humerus_shape_2D ~ ecotype, tree = tree_pr, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpanova, file= "humerus_ecotype_mvpanova.Rdata")
humerus_ecotype_mvpanova_EIC <- EIC(humerus_ecotype_mvpanova, nboot =1000)
save(humerus_ecotype_mvpanova_EIC, file= "humerus_ecotype_mvpanova_EIC.Rdata")
humerus_ecotype_mvpanova_test <- manova.gls(humerus_ecotype_mvpanova, test = "Pillai")
save(humerus_ecotype_mvpanova_test, file= "humerus_ecotype_mvpanova_test.Rdata")

#pANOVA
humerus_ecotype_panova <- procD.pgls(humerus_shape ~ ecotype, phy = tree_pr, iter = 999)
summary(humerus_ecotype_panova)
humerus_ecotype_pw <- pairwise(humerus_ecotype_panova, groups = ecotype)
summary(humerus_ecotype_pw)

#ANOVA
humerus_ecotype_anova <- procD.lm(humerus_shape ~ ecotype, iter = 999)
summary(humerus_ecotype_anova)
humerus_ecotype_pancova_pw <- pairwise(humerus_ecotype_anova, groups = ecotype)
summary(humerus_ecotype_pancova_pw)


##### humerus_shape AIC ####
geiger::aicw(c(humerus_shape_allometry_mv_EIC$EIC,
               humerus_ecotype_mvpanova_EIC$EIC,
               humerus_ecotype_mvpancova_EIC$EIC,
               humerus_ecotype_mvpancova_noint_EIC$EIC,
               humerus_shape_null_mv_EIC$EIC))


##### within ecotypes allometry #####
###### chipmunk ####
humerus_ecotype_mvpancova_chipmunk <- mvgls(humerus_shape_chipmunk_2D ~ lnhumerus_size_chipmunk, tree = tree_pr_chipmunk, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpancova_chipmunk, file= "humerus_ecotype_mvpancova_chipmunk.Rdata")
humerus_ecotype_mvpancova_test_chipmunk <- manova.gls(humerus_ecotype_mvpancova_chipmunk, type = "II", test = "Pillai")
save(humerus_ecotype_mvpancova_test_chipmunk, file= "humerus_ecotype_mvpancova_test_chipmunk.Rdata")

humerus_shape_allometry_chipmunk <- procD.pgls(humerus_shape_chipmunk ~ lnhumerus_size_chipmunk, phy = tree_pr_chipmunk,  iter = 999)
summary(humerus_shape_allometry_chipmunk)
humerus_allo_plot <- plotAllometry(humerus_shape_allometry_chipmunk, size = lnhumerus_size_chipmunk, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

humerus_shape_allometry_chipmunk_np <- procD.lm(humerus_shape_chipmunk ~ lnhumerus_size_chipmunk, iter = 999)
summary(humerus_shape_allometry_chipmunk_np)
humerus_allo_plot_np <- plotAllometry(humerus_shape_allometry_chipmunk_np, size = lnhumerus_size_chipmunk, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

###### gliding ####
humerus_shape_data_gliding <- filter(humerus_shape_data, ecotype == "gliding")
humerus_shape_gliding <- arrayspecs(humerus_shape_data_gliding[,3:584], 194,3 )
humerus_shape_gliding_2D <- two.d.array(humerus_shape_gliding)
lnhumerus_size_gliding <- humerus_shape_data_gliding$lnhumerus_size
tree_pr_gliding <- treedata(phy = tree, data = humerus_shape_data_gliding, sort=TRUE)$phy

humerus_ecotype_mvpancova_gliding <- mvgls(humerus_shape_gliding_2D ~ lnhumerus_size_gliding, tree = tree_pr_gliding, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpancova_gliding, file= "humerus_ecotype_mvpancova_gliding.Rdata")
humerus_ecotype_mvpancova_test_gliding <- manova.gls(humerus_ecotype_mvpancova_gliding, type = "II", test = "Pillai")
save(humerus_ecotype_mvpancova_test_gliding, file= "humerus_ecotype_mvpancova_test_gliding.Rdata")

humerus_shape_allometry_gliding <- procD.pgls(humerus_shape_gliding ~ lnhumerus_size_gliding, phy = tree_pr_gliding,  iter = 999)
summary(humerus_shape_allometry_gliding)
humerus_allo_plot <- plotAllometry(humerus_shape_allometry_gliding, size = lnhumerus_size_gliding, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

humerus_shape_allometry_gliding_np <- procD.lm(humerus_shape_gliding ~ lnhumerus_size_gliding, iter = 999)
summary(humerus_shape_allometry_gliding_np)
humerus_allo_plot_np <- plotAllometry(humerus_shape_allometry_gliding_np, size = lnhumerus_size_gliding, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

###### ground ####
humerus_shape_data_ground <- filter(humerus_shape_data, ecotype == "ground")
humerus_shape_ground <- arrayspecs(humerus_shape_data_ground[,3:584], 194,3 )
humerus_shape_ground_2D <- two.d.array(humerus_shape_ground)
lnhumerus_size_ground <- humerus_shape_data_ground$lnhumerus_size
tree_pr_ground <- treedata(phy = tree, data = humerus_shape_data_ground, sort=TRUE)$phy

humerus_ecotype_mvpancova_ground <- mvgls(humerus_shape_ground_2D ~ lnhumerus_size_ground, tree = tree_pr_ground, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpancova_ground, file= "humerus_ecotype_mvpancova_ground.Rdata")
humerus_ecotype_mvpancova_test_ground <- manova.gls(humerus_ecotype_mvpancova_ground, type = "II", test = "Pillai")
save(humerus_ecotype_mvpancova_test_ground, file= "humerus_ecotype_mvpancova_test_ground.Rdata")

humerus_shape_allometry_ground <- procD.pgls(humerus_shape_ground ~ lnhumerus_size_ground, phy = tree_pr_ground,  iter = 999)
summary(humerus_shape_allometry_ground)
humerus_allo_plot <- plotAllometry(humerus_shape_allometry_ground, size = lnhumerus_size_ground, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

humerus_shape_allometry_ground_np <- procD.lm(humerus_shape_ground ~ lnhumerus_size_ground, iter = 999)
summary(humerus_shape_allometry_ground_np)
humerus_allo_plot_np <- plotAllometry(humerus_shape_allometry_ground_np, size = lnhumerus_size_ground, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

###### tree ####
humerus_shape_data_tree <- filter(humerus_shape_data, ecotype == "tree")
humerus_shape_tree <- arrayspecs(humerus_shape_data_tree[,3:584], 194,3 )
humerus_shape_tree_2D <- two.d.array(humerus_shape_tree)
lnhumerus_size_tree <- humerus_shape_data_tree$lnhumerus_size
tree_pr_tree <- treedata(phy = tree, data = humerus_shape_data_tree, sort=TRUE)$phy

humerus_ecotype_mvpancova_tree <- mvgls(humerus_shape_tree_2D ~ lnhumerus_size_tree, tree = tree_pr_tree, method = "LOOCV", model = "lambda")
save(humerus_ecotype_mvpancova_tree, file= "humerus_ecotype_mvpancova_tree.Rdata")
humerus_ecotype_mvpancova_test_tree <- manova.gls(humerus_ecotype_mvpancova_tree, type = "II", test = "Pillai")
save(humerus_ecotype_mvpancova_test_tree, file= "humerus_ecotype_mvpancova_test_tree.Rdata")

humerus_shape_allometry_tree <- procD.pgls(humerus_shape_tree ~ lnhumerus_size_tree, phy = tree_pr_tree,  iter = 999)
summary(humerus_shape_allometry_tree)
humerus_allo_plot <- plotAllometry(humerus_shape_allometry_tree, size = lnhumerus_size_tree, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)

humerus_shape_allometry_tree_np <- procD.lm(humerus_shape_tree ~ lnhumerus_size_tree, iter = 999)
summary(humerus_shape_allometry_tree_np)
humerus_allo_plot_np <- plotAllometry(humerus_shape_allometry_tree_np, size = lnhumerus_size_tree, logsz = FALSE, method = "RegScore", cex = 1.5, las = 1)




#### humerus_Cg #####
##### null model: humerus_Cg ~ 1 ####
humerus_Cg_pgls_null <- phylolm(humerus_Cg~1, phy = tree_pr, model = "lambda", boot = 1000)


##### size model: humerus_Cg ~ size ####
humerus_Cg_pgls_size <- phylolm(humerus_Cg~lnhumerus_size, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_Cg_pgls_size)


##### size+ecotype model: humerus_Cg ~ size+ecotype ####
humerus_Cg_pANCOVAnoint_size_ecotype <- phylolm(humerus_Cg~lnhumerus_size+ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_Cg_pANCOVAnoint_size_ecotype)


##### size*ecotype model: humerus_Cg ~ size*ecotype ####
humerus_Cg_pANCOVA_size_ecotype <- phylolm(humerus_Cg~lnhumerus_size*ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_Cg_pANCOVA_size_ecotype)

######  Make table of coefficients (intercept and slope) ######  
humerus_Cg_pgls_size_all_int <- humerus_Cg_pgls_size$coefficients["(Intercept)"]
humerus_Cg_pgls_size_all_int_boot <- humerus_Cg_pgls_size$bootstrap[,"(Intercept)"]

humerus_Cg_pANCOVA_size_ecotype_chip_int <- humerus_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]
humerus_Cg_pANCOVA_size_ecotype_glide_int <- humerus_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_Cg_pANCOVA_size_ecotype$coefficients["ecotypegliding"]
humerus_Cg_pANCOVA_size_ecotype_ground_int <- humerus_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_Cg_pANCOVA_size_ecotype$coefficients["ecotypeground"]
humerus_Cg_pANCOVA_size_ecotype_tree_int <- humerus_Cg_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_Cg_pANCOVA_size_ecotype$coefficients["ecotypetree"]
humerus_Cg_pANCOVA_size_ecotype_chip_int_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"]
humerus_Cg_pANCOVA_size_ecotype_glide_int_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"ecotypegliding"]
humerus_Cg_pANCOVA_size_ecotype_ground_int_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"ecotypeground"]
humerus_Cg_pANCOVA_size_ecotype_tree_int_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"ecotypetree"]

humerus_Cg_pgls_size_all_slope <- humerus_Cg_pgls_size$coefficients["lnhumerus_size"]
humerus_Cg_pgls_size_all_slope_boot <- humerus_Cg_pgls_size$bootstrap[,"lnhumerus_size"]

humerus_Cg_pANCOVA_size_ecotype_chip_slope <- humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size"]
humerus_Cg_pANCOVA_size_ecotype_glide_slope <- humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypegliding"]
humerus_Cg_pANCOVA_size_ecotype_ground_slope <- humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypeground"]
humerus_Cg_pANCOVA_size_ecotype_tree_slope <- humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_Cg_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypetree"]
humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"]
humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypegliding"]
humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypeground"]
humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot <- humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_Cg_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypetree"]

humerus_Cg_pANCOVA_size_ecotype_int_cof <- rbind(
  humerus_Cg_pgls_size_all_int,
  humerus_Cg_pANCOVA_size_ecotype_chip_int,
  humerus_Cg_pANCOVA_size_ecotype_glide_int,
  humerus_Cg_pANCOVA_size_ecotype_ground_int,
  humerus_Cg_pANCOVA_size_ecotype_tree_int)
humerus_Cg_pANCOVA_size_ecotype_int_95 <- rbind(
  quantile(humerus_Cg_pgls_size_all_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_tree_int_boot, prob = c(0.025, 0.975)))
humerus_Cg_pANCOVA_size_ecotype_slope_cof <- rbind(
  humerus_Cg_pgls_size_all_slope,
  humerus_Cg_pANCOVA_size_ecotype_chip_slope,
  humerus_Cg_pANCOVA_size_ecotype_glide_slope,
  humerus_Cg_pANCOVA_size_ecotype_ground_slope,
  humerus_Cg_pANCOVA_size_ecotype_tree_slope)
humerus_Cg_pANCOVA_size_ecotype_slope_95 <- rbind(
  quantile(humerus_Cg_pgls_size_all_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot, prob = c(0.025, 0.975)))

humerus_Cg_pANCOVA_size_ecotype_cof <- data.frame(humerus_Cg_pANCOVA_size_ecotype_int_cof, humerus_Cg_pANCOVA_size_ecotype_int_95, humerus_Cg_pANCOVA_size_ecotype_slope_cof, humerus_Cg_pANCOVA_size_ecotype_slope_95 )
colnames(humerus_Cg_pANCOVA_size_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95", "Slope", "Slope_L95", "Slope_U95")
rownames(humerus_Cg_pANCOVA_size_ecotype_cof) <- c("all", "chip", "gliding", "ground", "tree")
round(humerus_Cg_pANCOVA_size_ecotype_cof, digits = 3)

######   plot humerus_Cg ~ lnsize*ecology ######  
plot(humerus_Cg ~ lnhumerus_size, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)
abline(a = humerus_Cg_pgls_size_all_int, b = humerus_Cg_pgls_size_all_slope, col = "black")
abline(a = humerus_Cg_pANCOVA_size_ecotype_chip_int, b = humerus_Cg_pANCOVA_size_ecotype_chip_slope, col = "#073b4c", lty = 2)
abline(a = humerus_Cg_pANCOVA_size_ecotype_glide_int, b = humerus_Cg_pANCOVA_size_ecotype_glide_slope, col = "#118ab2", lty = 2)
abline(a = humerus_Cg_pANCOVA_size_ecotype_ground_int, b = humerus_Cg_pANCOVA_size_ecotype_ground_slope, col = "#ffd166", lty = 2)
abline(a = humerus_Cg_pANCOVA_size_ecotype_tree_int, b = humerus_Cg_pANCOVA_size_ecotype_tree_slope, col = "#06d6a0")


###### slope differences ######  
# chip - glide
humerus_Cg_slope_size_chip_glide_SE <- sqrt((sd(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)))

humerus_Cg_slope_size_chip_glide_meandiff <- mean(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)

humerus_Cg_slope_size_chip_glide_CI <- mean(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot) + c(-1, 1) * 1.96 * humerus_Cg_slope_size_chip_glide_SE  # 95% confidence interval using z=1.96

humerus_Cg_slope_size_chip_glide_table <- cbind(humerus_Cg_slope_size_chip_glide_meandiff, humerus_Cg_slope_size_chip_glide_CI[1], humerus_Cg_slope_size_chip_glide_CI[2])
rownames(humerus_Cg_slope_size_chip_glide_table) <- ("humerus_Cg_slope_size_chip_glide")
colnames(humerus_Cg_slope_size_chip_glide_table) <- c("meandiff", "L95", "U95")
humerus_Cg_slope_size_chip_glide_table

# chip - ground
humerus_Cg_slope_size_chip_ground_SE <- sqrt((sd(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)))

humerus_Cg_slope_size_chip_ground_meandiff <- mean(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)

humerus_Cg_slope_size_chip_ground_CI <- mean(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * humerus_Cg_slope_size_chip_ground_SE  # 95% confidence interval using z=1.96

humerus_Cg_slope_size_chip_ground_table <- cbind(humerus_Cg_slope_size_chip_ground_meandiff, humerus_Cg_slope_size_chip_ground_CI[1], humerus_Cg_slope_size_chip_ground_CI[2])
rownames(humerus_Cg_slope_size_chip_ground_table) <- ("humerus_Cg_slope_size_chip_ground")
colnames(humerus_Cg_slope_size_chip_ground_table) <- c("meandiff", "L95", "U95")
humerus_Cg_slope_size_chip_ground_table

# chip - tree
humerus_Cg_slope_size_chip_tree_SE <- sqrt((sd(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)))

humerus_Cg_slope_size_chip_tree_meandiff <- mean(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)

humerus_Cg_slope_size_chip_tree_CI <- mean(humerus_Cg_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * humerus_Cg_slope_size_chip_tree_SE  # 95% confidence interval using z=1.96

humerus_Cg_slope_size_chip_tree_table <- cbind(humerus_Cg_slope_size_chip_tree_meandiff, humerus_Cg_slope_size_chip_tree_CI[1], humerus_Cg_slope_size_chip_tree_CI[2])
rownames(humerus_Cg_slope_size_chip_tree_table) <- ("humerus_Cg_slope_size_chip_tree")
colnames(humerus_Cg_slope_size_chip_tree_table) <- c("meandiff", "L95", "U95")
humerus_Cg_slope_size_chip_tree_table

# glide - ground
humerus_Cg_slope_size_glide_ground_SE <- sqrt((sd(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)))

humerus_Cg_slope_size_glide_ground_meandiff <- mean(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)

humerus_Cg_slope_size_glide_ground_CI <- mean(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * humerus_Cg_slope_size_glide_ground_SE  # 95% confidence interval using z=1.96

humerus_Cg_slope_size_glide_ground_table <- cbind(humerus_Cg_slope_size_glide_ground_meandiff, humerus_Cg_slope_size_glide_ground_CI[1], humerus_Cg_slope_size_glide_ground_CI[2])
rownames(humerus_Cg_slope_size_glide_ground_table) <- ("humerus_Cg_slope_size_glide_ground")
colnames(humerus_Cg_slope_size_glide_ground_table) <- c("meandiff", "L95", "U95")
humerus_Cg_slope_size_glide_ground_table

# glide - tree
humerus_Cg_slope_size_glide_tree_SE <- sqrt((sd(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)))

humerus_Cg_slope_size_glide_tree_meandiff <- mean(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)

humerus_Cg_slope_size_glide_tree_CI <- mean(humerus_Cg_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * humerus_Cg_slope_size_glide_tree_SE  # 95% confidence interval using z=1.96

humerus_Cg_slope_size_glide_tree_table <- cbind(humerus_Cg_slope_size_glide_tree_meandiff, humerus_Cg_slope_size_glide_tree_CI[1], humerus_Cg_slope_size_glide_tree_CI[2])
rownames(humerus_Cg_slope_size_glide_tree_table) <- ("humerus_Cg_slope_size_glide_tree")
colnames(humerus_Cg_slope_size_glide_tree_table) <- c("meandiff", "L95", "U95")
humerus_Cg_slope_size_glide_tree_table

# ground - tree
humerus_Cg_slope_size_ground_tree_SE <- sqrt((sd(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot)) + (sd(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)^2/length(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)))

humerus_Cg_slope_size_ground_tree_meandiff <- mean(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot)

humerus_Cg_slope_size_ground_tree_CI <- mean(humerus_Cg_pANCOVA_size_ecotype_ground_slope_boot) - mean(humerus_Cg_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * humerus_Cg_slope_size_ground_tree_SE  # 95% confidence interval using z=1.96

humerus_Cg_slope_size_ground_tree_table <- cbind(humerus_Cg_slope_size_ground_tree_meandiff, humerus_Cg_slope_size_ground_tree_CI[1], humerus_Cg_slope_size_ground_tree_CI[2])
rownames(humerus_Cg_slope_size_ground_tree_table) <- ("humerus_Cg_slope_size_ground_tree")
colnames(humerus_Cg_slope_size_ground_tree_table) <- c("meandiff", "L95", "U95")
humerus_Cg_slope_size_ground_tree_table

#table
humerus_Cg_slope_size_meandiff_table <- rbind(
  humerus_Cg_slope_size_chip_glide_table,
  humerus_Cg_slope_size_chip_ground_table,
  humerus_Cg_slope_size_chip_tree_table,
  humerus_Cg_slope_size_glide_ground_table,
  humerus_Cg_slope_size_glide_tree_table,
  humerus_Cg_slope_size_ground_tree_table)
round(humerus_Cg_slope_size_meandiff_table, digits = 2)


##### ecotype model: humerus_Cg ~ ecotype ####
humerus_Cg_pANOVA_ecotype <- phylolm(humerus_Cg~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_Cg_pANOVA_ecotype)

######  Make table of coefficients (intercept) ######  
humerus_Cg_pANOVA_ecotype_chip_int <- humerus_Cg_pANOVA_ecotype$coefficients["(Intercept)"]
humerus_Cg_pANOVA_ecotype_glide_int <- humerus_Cg_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_Cg_pANOVA_ecotype$coefficients["ecotypegliding"]
humerus_Cg_pANOVA_ecotype_ground_int <- humerus_Cg_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_Cg_pANOVA_ecotype$coefficients["ecotypeground"]
humerus_Cg_pANOVA_ecotype_tree_int <- humerus_Cg_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_Cg_pANOVA_ecotype$coefficients["ecotypetree"]
humerus_Cg_pANOVA_ecotype_chip_int_boot <- humerus_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"]
humerus_Cg_pANOVA_ecotype_glide_int_boot <- humerus_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_Cg_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
humerus_Cg_pANOVA_ecotype_ground_int_boot <- humerus_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_Cg_pANOVA_ecotype$bootstrap[,"ecotypeground"]
humerus_Cg_pANOVA_ecotype_tree_int_boot <- humerus_Cg_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_Cg_pANOVA_ecotype$bootstrap[,"ecotypetree"]

humerus_Cg_pANOVA_ecotype_int_cof <- rbind(
  humerus_Cg_pANOVA_ecotype_chip_int,
  humerus_Cg_pANOVA_ecotype_glide_int,
  humerus_Cg_pANOVA_ecotype_ground_int,
  humerus_Cg_pANOVA_ecotype_tree_int)
humerus_Cg_pANOVA_ecotype_int_95 <- rbind(
  quantile(humerus_Cg_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_Cg_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

humerus_Cg_pANOVA_ecotype_cof <- data.frame(humerus_Cg_pANOVA_ecotype_int_cof, humerus_Cg_pANOVA_ecotype_int_95)
colnames(humerus_Cg_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(humerus_Cg_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(humerus_Cg_pANOVA_ecotype_cof, digits = 3)

######   Plot humerus_Cg ~ ecotype ######  
humerus_Cg_vp_ecology <- ggplot(humerus_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= humerus_Cg, fill = ecotype)) + 
  scale_y_continuous(name="ln humerus_Cg (mm)", limits = c(.45,.85), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(humerus_Cg)) + 
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
humerus_Cg_vp_ecology


######  humerus_Cg_int_meandiff_table ######  
# chip - glide
humerus_Cg_int_chip_glide_SE <- sqrt((sd(humerus_Cg_pANOVA_ecotype_chip_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_Cg_pANOVA_ecotype_glide_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_glide_int_boot)))

humerus_Cg_int_chip_glide_meandiff <- mean(humerus_Cg_pANOVA_ecotype_chip_int_boot) - mean(humerus_Cg_pANOVA_ecotype_glide_int_boot)

humerus_Cg_int_chip_glide_CI <- mean(humerus_Cg_pANOVA_ecotype_chip_int_boot) - mean(humerus_Cg_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * humerus_Cg_int_chip_glide_SE  # 95% confidence interval using z=1.96

humerus_Cg_int_chip_glide_table <- cbind(humerus_Cg_int_chip_glide_meandiff, humerus_Cg_int_chip_glide_CI[1], humerus_Cg_int_chip_glide_CI[2])
rownames(humerus_Cg_int_chip_glide_table) <- ("humerus_Cg_int_chip_glide")
colnames(humerus_Cg_int_chip_glide_table) <- c("meandiff", "L95", "U95")
humerus_Cg_int_chip_glide_table

# chip - ground
humerus_Cg_int_chip_ground_SE <- sqrt((sd(humerus_Cg_pANOVA_ecotype_chip_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_Cg_pANOVA_ecotype_ground_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_ground_int_boot)))

humerus_Cg_int_chip_ground_meandiff <- mean(humerus_Cg_pANOVA_ecotype_chip_int_boot) - mean(humerus_Cg_pANOVA_ecotype_ground_int_boot)

humerus_Cg_int_chip_ground_CI <- mean(humerus_Cg_pANOVA_ecotype_chip_int_boot) - mean(humerus_Cg_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_Cg_int_chip_ground_SE  # 95% confidence interval using z=1.96

humerus_Cg_int_chip_ground_table <- cbind(humerus_Cg_int_chip_ground_meandiff, humerus_Cg_int_chip_ground_CI[1], humerus_Cg_int_chip_ground_CI[2])
rownames(humerus_Cg_int_chip_ground_table) <- ("humerus_Cg_int_chip_ground")
colnames(humerus_Cg_int_chip_ground_table) <- c("meandiff", "L95", "U95")
humerus_Cg_int_chip_ground_table

# chip - tree
humerus_Cg_int_chip_tree_SE <- sqrt((sd(humerus_Cg_pANOVA_ecotype_chip_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_Cg_pANOVA_ecotype_tree_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_tree_int_boot)))

humerus_Cg_int_chip_tree_meandiff <- mean(humerus_Cg_pANOVA_ecotype_chip_int_boot) - mean(humerus_Cg_pANOVA_ecotype_tree_int_boot)

humerus_Cg_int_chip_tree_CI <- mean(humerus_Cg_pANOVA_ecotype_chip_int_boot) - mean(humerus_Cg_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_Cg_int_chip_tree_SE  #

humerus_Cg_int_chip_tree_table <- cbind(humerus_Cg_int_chip_tree_meandiff, humerus_Cg_int_chip_tree_CI[1], humerus_Cg_int_chip_tree_CI[2])
rownames(humerus_Cg_int_chip_tree_table) <- ("humerus_Cg_int_chip_tree")
colnames(humerus_Cg_int_chip_tree_table) <- c("meandiff", "L95", "U95")
humerus_Cg_int_chip_tree_table

# glide - ground
humerus_Cg_int_glide_ground_SE <- sqrt((sd(humerus_Cg_pANOVA_ecotype_glide_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_Cg_pANOVA_ecotype_ground_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_ground_int_boot)))

humerus_Cg_int_glide_ground_meandiff <- mean(humerus_Cg_pANOVA_ecotype_glide_int_boot) - mean(humerus_Cg_pANOVA_ecotype_ground_int_boot)

humerus_Cg_int_glide_ground_CI <- mean(humerus_Cg_pANOVA_ecotype_glide_int_boot) - mean(humerus_Cg_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_Cg_int_glide_ground_SE  

humerus_Cg_int_glide_ground_table <- cbind(humerus_Cg_int_glide_ground_meandiff, humerus_Cg_int_glide_ground_CI[1], humerus_Cg_int_glide_ground_CI[2])
rownames(humerus_Cg_int_glide_ground_table) <- ("humerus_Cg_int_glide_ground")
colnames(humerus_Cg_int_glide_ground_table) <- c("meandiff", "L95", "U95")
humerus_Cg_int_glide_ground_table

# glide - tree
humerus_Cg_int_glide_tree_SE <- sqrt((sd(humerus_Cg_pANOVA_ecotype_glide_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_Cg_pANOVA_ecotype_tree_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_tree_int_boot)))

humerus_Cg_int_glide_tree_meandiff <- mean(humerus_Cg_pANOVA_ecotype_glide_int_boot) - mean(humerus_Cg_pANOVA_ecotype_tree_int_boot)

humerus_Cg_int_glide_tree_CI <- mean(humerus_Cg_pANOVA_ecotype_glide_int_boot) - mean(humerus_Cg_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_Cg_int_glide_tree_SE  

humerus_Cg_int_glide_tree_table <- cbind(humerus_Cg_int_glide_tree_meandiff, humerus_Cg_int_glide_tree_CI[1], humerus_Cg_int_glide_tree_CI[2])
rownames(humerus_Cg_int_glide_tree_table) <- ("humerus_Cg_int_glide_tree")
colnames(humerus_Cg_int_glide_tree_table) <- c("meandiff", "L95", "U95")
humerus_Cg_int_glide_tree_table

# ground - tree
humerus_Cg_int_ground_tree_SE <- sqrt((sd(humerus_Cg_pANOVA_ecotype_ground_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_ground_int_boot)) + (sd(humerus_Cg_pANOVA_ecotype_tree_int_boot)^2/length(humerus_Cg_pANOVA_ecotype_tree_int_boot)))

humerus_Cg_int_ground_tree_meandiff <- mean(humerus_Cg_pANOVA_ecotype_ground_int_boot) - mean(humerus_Cg_pANOVA_ecotype_tree_int_boot)

humerus_Cg_int_ground_tree_CI <- mean(humerus_Cg_pANOVA_ecotype_ground_int_boot) - mean(humerus_Cg_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_Cg_int_ground_tree_SE  # 95% confidence interval using z=1.96

humerus_Cg_int_ground_tree_table <- cbind(humerus_Cg_int_ground_tree_meandiff, humerus_Cg_int_ground_tree_CI[1], humerus_Cg_int_ground_tree_CI[2])
rownames(humerus_Cg_int_ground_tree_table) <- ("humerus_Cg_int_ground_tree")
colnames(humerus_Cg_int_ground_tree_table) <- c("meandiff", "L95", "U95")
humerus_Cg_int_ground_tree_table

#table
humerus_Cg_int_meandiff_table <- rbind(
  humerus_Cg_int_chip_glide_table,
  humerus_Cg_int_chip_ground_table,
  humerus_Cg_int_chip_tree_table,
  humerus_Cg_int_glide_ground_table,
  humerus_Cg_int_glide_tree_table,
  humerus_Cg_int_ground_tree_table)
round(humerus_Cg_int_meandiff_table, digits = 2)


##### humerus_Cg AIC ####
geiger::aicw(c(AIC(humerus_Cg_pgls_size),
               AIC(humerus_Cg_pANOVA_ecotype),
               AIC(humerus_Cg_pANCOVA_size_ecotype),
               AIC(humerus_Cg_pANCOVAnoint_size_ecotype),
               AIC(humerus_Cg_pgls_null)))



#### humerus_DE #####
##### null model: humerus_DE ~ 1 ##### 
humerus_DE_pgls_null <- phylolm(humerus_DE~1, phy = tree_pr, model = "lambda", boot = 1000)

##### size model: humerus_DE ~ size ##### 
humerus_DE_pgls_size <- phylolm(humerus_DE~lnhumerus_size, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_DE_pgls_size)

##### size+ecotype model: humerus_DE ~ size+ecotype ##### 
humerus_DE_pANCOVA_size_ecotype_noint <- phylolm(humerus_DE~lnhumerus_size+ecotype, phy = tree_pr, model = "lambda", boot = 1000, lower.bound = 0)
summary(humerus_DE_pANCOVA_size_ecotype_noint)

##### size*ecotype model: humerus_DE ~ size*ecotype ##### 
humerus_DE_pANCOVA_size_ecotype <- phylolm(humerus_DE~lnhumerus_size*ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_DE_pANCOVA_size_ecotype)

######  Make table of coefficients (intercept and slope) ######  
humerus_DE_pgls_size_all_int <- humerus_DE_pgls_size$coefficients["(Intercept)"]
humerus_DE_pgls_size_all_int_boot <- humerus_DE_pgls_size$bootstrap[,"(Intercept)"]

humerus_DE_pANCOVA_size_ecotype_chip_int <- humerus_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]
humerus_DE_pANCOVA_size_ecotype_glide_int <- humerus_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_DE_pANCOVA_size_ecotype$coefficients["ecotypegliding"]
humerus_DE_pANCOVA_size_ecotype_ground_int <- humerus_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_DE_pANCOVA_size_ecotype$coefficients["ecotypeground"]
humerus_DE_pANCOVA_size_ecotype_tree_int <- humerus_DE_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_DE_pANCOVA_size_ecotype$coefficients["ecotypetree"]
humerus_DE_pANCOVA_size_ecotype_chip_int_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"]
humerus_DE_pANCOVA_size_ecotype_glide_int_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_DE_pANCOVA_size_ecotype$bootstrap[,"ecotypegliding"]
humerus_DE_pANCOVA_size_ecotype_ground_int_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_DE_pANCOVA_size_ecotype$bootstrap[,"ecotypeground"]
humerus_DE_pANCOVA_size_ecotype_tree_int_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_DE_pANCOVA_size_ecotype$bootstrap[,"ecotypetree"]

humerus_DE_pgls_size_all_slope <- humerus_DE_pgls_size$coefficients["lnhumerus_size"]
humerus_DE_pgls_size_all_slope_boot <- humerus_DE_pgls_size$bootstrap[,"lnhumerus_size"]

humerus_DE_pANCOVA_size_ecotype_chip_slope <- humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size"]
humerus_DE_pANCOVA_size_ecotype_glide_slope <- humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypegliding"]
humerus_DE_pANCOVA_size_ecotype_ground_slope <- humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypeground"]
humerus_DE_pANCOVA_size_ecotype_tree_slope <- humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_DE_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypetree"]
humerus_DE_pANCOVA_size_ecotype_chip_slope_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"]
humerus_DE_pANCOVA_size_ecotype_glide_slope_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypegliding"]
humerus_DE_pANCOVA_size_ecotype_ground_slope_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypeground"]
humerus_DE_pANCOVA_size_ecotype_tree_slope_boot <- humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_DE_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypetree"]

humerus_DE_pANCOVA_size_ecotype_int_cof <- rbind(
  humerus_DE_pgls_size_all_int,
  humerus_DE_pANCOVA_size_ecotype_chip_int,
  humerus_DE_pANCOVA_size_ecotype_glide_int,
  humerus_DE_pANCOVA_size_ecotype_ground_int,
  humerus_DE_pANCOVA_size_ecotype_tree_int)
humerus_DE_pANCOVA_size_ecotype_int_95 <- rbind(
  quantile(humerus_DE_pgls_size_all_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_tree_int_boot, prob = c(0.025, 0.975)))
humerus_DE_pANCOVA_size_ecotype_slope_cof <- rbind(
  humerus_DE_pgls_size_all_slope,
  humerus_DE_pANCOVA_size_ecotype_chip_slope,
  humerus_DE_pANCOVA_size_ecotype_glide_slope,
  humerus_DE_pANCOVA_size_ecotype_ground_slope,
  humerus_DE_pANCOVA_size_ecotype_tree_slope)
humerus_DE_pANCOVA_size_ecotype_slope_95 <- rbind(
  quantile(humerus_DE_pgls_size_all_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot, prob = c(0.025, 0.975)))

humerus_DE_pANCOVA_size_ecotype_cof <- data.frame(humerus_DE_pANCOVA_size_ecotype_int_cof, humerus_DE_pANCOVA_size_ecotype_int_95, humerus_DE_pANCOVA_size_ecotype_slope_cof, humerus_DE_pANCOVA_size_ecotype_slope_95 )
colnames(humerus_DE_pANCOVA_size_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95", "Slope", "Slope_L95", "Slope_U95")
rownames(humerus_DE_pANCOVA_size_ecotype_cof) <- c("all", "chip", "gliding", "ground", "tree")
round(humerus_DE_pANCOVA_size_ecotype_cof, digits = 3)

######   plot humerus_DE ~ lnsize*ecology ######  
plot(humerus_DE ~ lnhumerus_size, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)
abline(a = humerus_DE_pgls_size_all_int, b = humerus_DE_pgls_size_all_slope, col = "black", lty = 2)
abline(a = humerus_DE_pANCOVA_size_ecotype_chip_int, b = humerus_DE_pANCOVA_size_ecotype_chip_slope, col = "#073b4c", lty = 2)
abline(a = humerus_DE_pANCOVA_size_ecotype_glide_int, b = humerus_DE_pANCOVA_size_ecotype_glide_slope, col = "#118ab2", lty = 2)
abline(a = humerus_DE_pANCOVA_size_ecotype_ground_int, b = humerus_DE_pANCOVA_size_ecotype_ground_slope, col = "#ffd166", lty = 2)
abline(a = humerus_DE_pANCOVA_size_ecotype_tree_int, b = humerus_DE_pANCOVA_size_ecotype_tree_slope, col = "#06d6a0")


######  humerus_DE_slope_size_meandiff_table ######  
# chip - glide
humerus_DE_slope_size_chip_glide_SE <- sqrt((sd(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)))

humerus_DE_slope_size_chip_glide_meandiff <- mean(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)

humerus_DE_slope_size_chip_glide_CI <- mean(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot) + c(-1, 1) * 1.96 * humerus_DE_slope_size_chip_glide_SE  # 95% confidence interval using z=1.96

humerus_DE_slope_size_chip_glide_table <- cbind(humerus_DE_slope_size_chip_glide_meandiff, humerus_DE_slope_size_chip_glide_CI[1], humerus_DE_slope_size_chip_glide_CI[2])
rownames(humerus_DE_slope_size_chip_glide_table) <- ("humerus_DE_slope_size_chip_glide")
colnames(humerus_DE_slope_size_chip_glide_table) <- c("meandiff", "L95", "U95")
humerus_DE_slope_size_chip_glide_table

# chip - ground
humerus_DE_slope_size_chip_ground_SE <- sqrt((sd(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)))

humerus_DE_slope_size_chip_ground_meandiff <- mean(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)

humerus_DE_slope_size_chip_ground_CI <- mean(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * humerus_DE_slope_size_chip_ground_SE  # 95% confidence interval using z=1.96

humerus_DE_slope_size_chip_ground_table <- cbind(humerus_DE_slope_size_chip_ground_meandiff, humerus_DE_slope_size_chip_ground_CI[1], humerus_DE_slope_size_chip_ground_CI[2])
rownames(humerus_DE_slope_size_chip_ground_table) <- ("humerus_DE_slope_size_chip_ground")
colnames(humerus_DE_slope_size_chip_ground_table) <- c("meandiff", "L95", "U95")
humerus_DE_slope_size_chip_ground_table

# chip - tree
humerus_DE_slope_size_chip_tree_SE <- sqrt((sd(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot)) + (sd(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)))

humerus_DE_slope_size_chip_tree_meandiff <- mean(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)

humerus_DE_slope_size_chip_tree_CI <- mean(humerus_DE_pANCOVA_size_ecotype_chip_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * humerus_DE_slope_size_chip_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_slope_size_chip_tree_table <- cbind(humerus_DE_slope_size_chip_tree_meandiff, humerus_DE_slope_size_chip_tree_CI[1], humerus_DE_slope_size_chip_tree_CI[2])
rownames(humerus_DE_slope_size_chip_tree_table) <- ("humerus_DE_slope_size_chip_tree")
colnames(humerus_DE_slope_size_chip_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_slope_size_chip_tree_table

# glide - ground
humerus_DE_slope_size_glide_ground_SE <- sqrt((sd(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)))

humerus_DE_slope_size_glide_ground_meandiff <- mean(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)

humerus_DE_slope_size_glide_ground_CI <- mean(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot) + c(-1, 1) * 1.96 * humerus_DE_slope_size_glide_ground_SE  # 95% confidence interval using z=1.96

humerus_DE_slope_size_glide_ground_table <- cbind(humerus_DE_slope_size_glide_ground_meandiff, humerus_DE_slope_size_glide_ground_CI[1], humerus_DE_slope_size_glide_ground_CI[2])
rownames(humerus_DE_slope_size_glide_ground_table) <- ("humerus_DE_slope_size_glide_ground")
colnames(humerus_DE_slope_size_glide_ground_table) <- c("meandiff", "L95", "U95")
humerus_DE_slope_size_glide_ground_table

# glide - tree
humerus_DE_slope_size_glide_tree_SE <- sqrt((sd(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot)) + (sd(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)))

humerus_DE_slope_size_glide_tree_meandiff <- mean(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)

humerus_DE_slope_size_glide_tree_CI <- mean(humerus_DE_pANCOVA_size_ecotype_glide_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * humerus_DE_slope_size_glide_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_slope_size_glide_tree_table <- cbind(humerus_DE_slope_size_glide_tree_meandiff, humerus_DE_slope_size_glide_tree_CI[1], humerus_DE_slope_size_glide_tree_CI[2])
rownames(humerus_DE_slope_size_glide_tree_table) <- ("humerus_DE_slope_size_glide_tree")
colnames(humerus_DE_slope_size_glide_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_slope_size_glide_tree_table

# ground - tree
humerus_DE_slope_size_ground_tree_SE <- sqrt((sd(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot)) + (sd(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)^2/length(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)))

humerus_DE_slope_size_ground_tree_meandiff <- mean(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot)

humerus_DE_slope_size_ground_tree_CI <- mean(humerus_DE_pANCOVA_size_ecotype_ground_slope_boot) - mean(humerus_DE_pANCOVA_size_ecotype_tree_slope_boot) + c(-1, 1) * 1.96 * humerus_DE_slope_size_ground_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_slope_size_ground_tree_table <- cbind(humerus_DE_slope_size_ground_tree_meandiff, humerus_DE_slope_size_ground_tree_CI[1], humerus_DE_slope_size_ground_tree_CI[2])
rownames(humerus_DE_slope_size_ground_tree_table) <- ("humerus_DE_slope_size_ground_tree")
colnames(humerus_DE_slope_size_ground_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_slope_size_ground_tree_table

#table
humerus_DE_slope_size_meandiff_table <- rbind(
  humerus_DE_slope_size_chip_glide_table,
  humerus_DE_slope_size_chip_ground_table,
  humerus_DE_slope_size_chip_tree_table,
  humerus_DE_slope_size_glide_ground_table,
  humerus_DE_slope_size_glide_tree_table,
  humerus_DE_slope_size_ground_tree_table)
round(humerus_DE_slope_size_meandiff_table, digits = 2)


##### ecotype model: humerus_DE ~ ecotype #####
humerus_DE_pANOVA_ecotype <- phylolm(humerus_DE~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_DE_pANOVA_ecotype)

#pANOVA
humerus_DE_pANOVA_ecotype_BM_rrpp <- lm.rrpp(humerus_DE~ecotype, iter = 999, Cov = vcv.phylo(tree_pr), SS.type = "II")
anova(humerus_DE_pANOVA_ecotype_BM_rrpp)
humerus_DE_pANOVA_ecotype_BM_rrpp_PW <- pairwise(humerus_DE_pANOVA_ecotype_BM_rrpp, groups = ecotype )
summary(humerus_DE_pANOVA_ecotype_BM_rrpp_PW, show.vectors = TRUE)

#ANOVA
humerus_DE_ANOVA_rrpp <- lm.rrpp(humerus_DE~ecotype, iter = 999, SS.type = "II")
anova(humerus_DE_ANOVA_rrpp)
humerus_DE_ANOVA_rrpp_PW <- pairwise(humerus_DE_ANOVA_rrpp, groups = ecotype )
summary(humerus_DE_ANOVA_rrpp_PW, show.vectors = TRUE)

###### Make table of coefficients (intercept) ##### 
humerus_DE_pANOVA_ecotype_chip_int <- humerus_DE_pANOVA_ecotype$coefficients["(Intercept)"]
humerus_DE_pANOVA_ecotype_glide_int <- humerus_DE_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_DE_pANOVA_ecotype$coefficients["ecotypegliding"]
humerus_DE_pANOVA_ecotype_ground_int <- humerus_DE_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_DE_pANOVA_ecotype$coefficients["ecotypeground"]
humerus_DE_pANOVA_ecotype_tree_int <- humerus_DE_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_DE_pANOVA_ecotype$coefficients["ecotypetree"]
humerus_DE_pANOVA_ecotype_chip_int_boot <- humerus_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"]
humerus_DE_pANOVA_ecotype_glide_int_boot <- humerus_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_DE_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
humerus_DE_pANOVA_ecotype_ground_int_boot <- humerus_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_DE_pANOVA_ecotype$bootstrap[,"ecotypeground"]
humerus_DE_pANOVA_ecotype_tree_int_boot <- humerus_DE_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_DE_pANOVA_ecotype$bootstrap[,"ecotypetree"]

humerus_DE_pANOVA_ecotype_int_cof <- rbind(
  humerus_DE_pANOVA_ecotype_chip_int,
  humerus_DE_pANOVA_ecotype_glide_int,
  humerus_DE_pANOVA_ecotype_ground_int,
  humerus_DE_pANOVA_ecotype_tree_int)
humerus_DE_pANOVA_ecotype_int_95 <- rbind(
  quantile(humerus_DE_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_DE_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

humerus_DE_pANOVA_ecotype_cof <- data.frame(humerus_DE_pANOVA_ecotype_int_cof, humerus_DE_pANOVA_ecotype_int_95)
colnames(humerus_DE_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(humerus_DE_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(humerus_DE_pANOVA_ecotype_cof, digits = 3)

##### Plot humerus_DE ~ ecotype ##### 
humerus_DE_vp_ecology <- ggplot(humerus_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= humerus_DE, fill = ecotype)) + 
  scale_y_continuous(name="ln humerus_DE (mm)", limits = c(11,25), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(humerus_DE)) + 
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
humerus_DE_vp_ecology


##### humerus_DE_int_meandiff_table##### 
# chip - glide
humerus_DE_int_chip_glide_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_chip_int_boot)^2/length(humerus_DE_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_glide_int_boot)^2/length(humerus_DE_pANOVA_ecotype_glide_int_boot)))

humerus_DE_int_chip_glide_meandiff <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_glide_int_boot)

humerus_DE_int_chip_glide_CI <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_chip_glide_SE  # 95% confidence interval using z=1.96

humerus_DE_int_chip_glide_table <- cbind(humerus_DE_int_chip_glide_meandiff, humerus_DE_int_chip_glide_CI[1], humerus_DE_int_chip_glide_CI[2])
rownames(humerus_DE_int_chip_glide_table) <- ("humerus_DE_int_chip_glide")
colnames(humerus_DE_int_chip_glide_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_chip_glide_table

# chip - ground
humerus_DE_int_chip_ground_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_chip_int_boot)^2/length(humerus_DE_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_ground_int_boot)^2/length(humerus_DE_pANOVA_ecotype_ground_int_boot)))

humerus_DE_int_chip_ground_meandiff <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot)

humerus_DE_int_chip_ground_CI <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_chip_ground_SE  # 95% confidence interval using z=1.96

humerus_DE_int_chip_ground_table <- cbind(humerus_DE_int_chip_ground_meandiff, humerus_DE_int_chip_ground_CI[1], humerus_DE_int_chip_ground_CI[2])
rownames(humerus_DE_int_chip_ground_table) <- ("humerus_DE_int_chip_ground")
colnames(humerus_DE_int_chip_ground_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_chip_ground_table

# chip - tree
humerus_DE_int_chip_tree_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_chip_int_boot)^2/length(humerus_DE_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_tree_int_boot)^2/length(humerus_DE_pANOVA_ecotype_tree_int_boot)))

humerus_DE_int_chip_tree_meandiff <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot)

humerus_DE_int_chip_tree_CI <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_chip_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_int_chip_tree_table <- cbind(humerus_DE_int_chip_tree_meandiff, humerus_DE_int_chip_tree_CI[1], humerus_DE_int_chip_tree_CI[2])
rownames(humerus_DE_int_chip_tree_table) <- ("humerus_DE_int_chip_tree")
colnames(humerus_DE_int_chip_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_chip_tree_table

# glide - ground
humerus_DE_int_glide_ground_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_glide_int_boot)^2/length(humerus_DE_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_ground_int_boot)^2/length(humerus_DE_pANOVA_ecotype_ground_int_boot)))

humerus_DE_int_glide_ground_meandiff <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot)

humerus_DE_int_glide_ground_CI <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_glide_ground_SE  # 95% confidence interval using z=1.96

humerus_DE_int_glide_ground_table <- cbind(humerus_DE_int_glide_ground_meandiff, humerus_DE_int_glide_ground_CI[1], humerus_DE_int_glide_ground_CI[2])
rownames(humerus_DE_int_glide_ground_table) <- ("humerus_DE_int_glide_ground")
colnames(humerus_DE_int_glide_ground_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_glide_ground_table

# glide - tree
humerus_DE_int_glide_tree_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_glide_int_boot)^2/length(humerus_DE_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_tree_int_boot)^2/length(humerus_DE_pANOVA_ecotype_tree_int_boot)))

humerus_DE_int_glide_tree_meandiff <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot)

humerus_DE_int_glide_tree_CI <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_glide_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_int_glide_tree_table <- cbind(humerus_DE_int_glide_tree_meandiff, humerus_DE_int_glide_tree_CI[1], humerus_DE_int_glide_tree_CI[2])
rownames(humerus_DE_int_glide_tree_table) <- ("humerus_DE_int_glide_tree")
colnames(humerus_DE_int_glide_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_glide_tree_table

# ground - tree
humerus_DE_int_ground_tree_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_ground_int_boot)^2/length(humerus_DE_pANOVA_ecotype_ground_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_tree_int_boot)^2/length(humerus_DE_pANOVA_ecotype_tree_int_boot)))

humerus_DE_int_ground_tree_meandiff <- mean(humerus_DE_pANOVA_ecotype_ground_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot)

humerus_DE_int_ground_tree_CI <- mean(humerus_DE_pANOVA_ecotype_ground_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_ground_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_int_ground_tree_table <- cbind(humerus_DE_int_ground_tree_meandiff, humerus_DE_int_ground_tree_CI[1], humerus_DE_int_ground_tree_CI[2])
rownames(humerus_DE_int_ground_tree_table) <- ("humerus_DE_int_ground_tree")
colnames(humerus_DE_int_ground_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_ground_tree_table

#table
humerus_DE_int_meandiff_table <- rbind(
  humerus_DE_int_chip_glide_table,
  humerus_DE_int_chip_ground_table,
  humerus_DE_int_chip_tree_table,
  humerus_DE_int_glide_ground_table,
  humerus_DE_int_glide_tree_table,
  humerus_DE_int_ground_tree_table)
round(humerus_DE_int_meandiff_table, digits = 2)


##### humerus_DE AIC ####
geiger::aicw(c(AIC(humerus_DE_pgls_size),
               AIC(size_humerus_DE_pANOVA_ecotype),
               AIC(humerus_DE_pANCOVA_size_ecotype),
               AIC(humerus_DE_pANCOVA_size_ecotype_noint),
               AIC(humerus_DE_pgls_null)))


#### humerus_CSS #####
##### null model: humerus_CSS ~1 ##### 
humerus_CSS_pgls_null <- phylolm(humerus_CSS~1, phy = tree_pr, model = "lambda", boot = 1000)


##### size model: humerus_CSS ~ size ##### 
humerus_CSS_pgls_size <- phylolm(humerus_CSS~lnhumerus_size, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_CSS_pgls_size)


##### size+ecotype model: humerus_CSS ~ size+ecotype ##### 
humerus_CSS_pANCOVAnoint_size_ecotype <- phylolm(humerus_CSS~lnhumerus_size+ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_CSS_pANCOVAnoint_size_ecotype)


##### size*ecotype model: humerus_CSS ~ size*ecotype ##### 
humerus_CSS_pANCOVA_size_ecotype <- phylolm(humerus_CSS~lnhumerus_size*ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_CSS_pANCOVA_size_ecotype)

###### Make table of coefficients (intercept and slope) ###### 
humerus_CSS_pgls_size_all_int <- humerus_CSS_pgls_size$coefficients["(Intercept)"]
humerus_CSS_pgls_size_all_int_boot <- humerus_CSS_pgls_size$bootstrap[,"(Intercept)"]

humerus_CSS_pANCOVA_size_ecotype_chip_int <- humerus_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]
humerus_CSS_pANCOVA_size_ecotype_glide_int <- humerus_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_CSS_pANCOVA_size_ecotype$coefficients["ecotypegliding"]
humerus_CSS_pANCOVA_size_ecotype_ground_int <- humerus_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_CSS_pANCOVA_size_ecotype$coefficients["ecotypeground"]
humerus_CSS_pANCOVA_size_ecotype_tree_int <- humerus_CSS_pANCOVA_size_ecotype$coefficients["(Intercept)"]+humerus_CSS_pANCOVA_size_ecotype$coefficients["ecotypetree"]
humerus_CSS_pANCOVA_size_ecotype_chip_int_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"]
humerus_CSS_pANCOVA_size_ecotype_glide_int_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"ecotypegliding"]
humerus_CSS_pANCOVA_size_ecotype_ground_int_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"ecotypeground"]
humerus_CSS_pANCOVA_size_ecotype_tree_int_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"(Intercept)"] + humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"ecotypetree"]

humerus_CSS_pgls_size_all_slope <- humerus_CSS_pgls_size$coefficients["lnhumerus_size"]
humerus_CSS_pgls_size_all_slope_boot <- humerus_CSS_pgls_size$bootstrap[,"lnhumerus_size"]

humerus_CSS_pANCOVA_size_ecotype_chip_slope <- humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size"]
humerus_CSS_pANCOVA_size_ecotype_glide_slope <- humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypegliding"]
humerus_CSS_pANCOVA_size_ecotype_ground_slope <- humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypeground"]
humerus_CSS_pANCOVA_size_ecotype_tree_slope <- humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size"] + humerus_CSS_pANCOVA_size_ecotype$coefficients["lnhumerus_size:ecotypetree"]
humerus_CSS_pANCOVA_size_ecotype_chip_slope_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"]
humerus_CSS_pANCOVA_size_ecotype_glide_slope_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypegliding"]
humerus_CSS_pANCOVA_size_ecotype_ground_slope_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypeground"]
humerus_CSS_pANCOVA_size_ecotype_tree_slope_boot <- humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size"] + humerus_CSS_pANCOVA_size_ecotype$bootstrap[,"lnhumerus_size:ecotypetree"]

humerus_CSS_pANCOVA_size_ecotype_int_cof <- rbind(
  humerus_CSS_pgls_size_all_int,
  humerus_CSS_pANCOVA_size_ecotype_chip_int,
  humerus_CSS_pANCOVA_size_ecotype_glide_int,
  humerus_CSS_pANCOVA_size_ecotype_ground_int,
  humerus_CSS_pANCOVA_size_ecotype_tree_int)
humerus_CSS_pANCOVA_size_ecotype_int_95 <- rbind(
  quantile(humerus_CSS_pgls_size_all_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_tree_int_boot, prob = c(0.025, 0.975)))
humerus_CSS_pANCOVA_size_ecotype_slope_cof <- rbind(
  humerus_CSS_pgls_size_all_slope,
  humerus_CSS_pANCOVA_size_ecotype_chip_slope,
  humerus_CSS_pANCOVA_size_ecotype_glide_slope,
  humerus_CSS_pANCOVA_size_ecotype_ground_slope,
  humerus_CSS_pANCOVA_size_ecotype_tree_slope)
humerus_CSS_pANCOVA_size_ecotype_slope_95 <- rbind(
  quantile(humerus_CSS_pgls_size_all_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_chip_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_glide_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_ground_slope_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANCOVA_size_ecotype_tree_slope_boot, prob = c(0.025, 0.975)))

humerus_CSS_pANCOVA_size_ecotype_cof <- data.frame(humerus_CSS_pANCOVA_size_ecotype_int_cof, humerus_CSS_pANCOVA_size_ecotype_int_95, humerus_CSS_pANCOVA_size_ecotype_slope_cof, humerus_CSS_pANCOVA_size_ecotype_slope_95 )
colnames(humerus_CSS_pANCOVA_size_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95", "Slope", "Slope_L95", "Slope_U95")
rownames(humerus_CSS_pANCOVA_size_ecotype_cof) <- c("all", "chip", "gliding", "ground", "tree")
round(humerus_CSS_pANCOVA_size_ecotype_cof, digits = 3)

######  plot humerus_CSS ~ lnsize*ecology ###### 
plot(humerus_CSS ~ lnhumerus_size, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)
abline(a = humerus_CSS_pgls_size_all_int, b = humerus_CSS_pgls_size_all_slope, col = "black")
abline(a = humerus_CSS_pANCOVA_size_ecotype_chip_int, b = humerus_CSS_pANCOVA_size_ecotype_chip_slope, col = "#073b4c", lty = 2)
abline(a = humerus_CSS_pANCOVA_size_ecotype_glide_int, b = humerus_CSS_pANCOVA_size_ecotype_glide_slope, col = "#118ab2", lty = 2)
abline(a = humerus_CSS_pANCOVA_size_ecotype_ground_int, b = humerus_CSS_pANCOVA_size_ecotype_ground_slope, col = "#ffd166")
abline(a = humerus_CSS_pANCOVA_size_ecotype_tree_int, b = humerus_CSS_pANCOVA_size_ecotype_tree_slope, col = "#06d6a0", lty = 2)


##### ecotype model: humerus_CSS ~ ecotype #####
humerus_CSS_pANOVA_ecotype <- phylolm(humerus_CSS~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_CSS_pANOVA_ecotype)

#pANOVA #
humerus_CSS_pANOVA_ecotype_BM_rrpp <- lm.rrpp(humerus_CSS~ecotype, iter = 999, Cov = vcv.phylo(tree_pr), SS.type = "II")
anova(humerus_CSS_pANOVA_ecotype_BM_rrpp)
humerus_CSS_pANOVA_ecotype_BM_rrpp_PW <- pairwise(humerus_CSS_pANOVA_ecotype_BM_rrpp, groups = ecotype )
summary(humerus_CSS_pANOVA_ecotype_BM_rrpp_PW, show.vectors = TRUE)

#ANOVA 
humerus_CSS_ANOVA_rrpp <- lm.rrpp(humerus_CSS~ecotype, iter = 999, SS.type = "II")
anova(humerus_CSS_ANOVA_rrpp)
humerus_CSS_ANOVA_rrpp_PW <- pairwise(humerus_CSS_ANOVA_rrpp, groups = ecotype )
summary(humerus_CSS_ANOVA_rrpp_PW, show.vectors = TRUE)

######  Make table of coefficients (intercept) ######  
humerus_CSS_pANOVA_ecotype_chip_int <- humerus_CSS_pANOVA_ecotype$coefficients["(Intercept)"]
humerus_CSS_pANOVA_ecotype_glide_int <- humerus_CSS_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_CSS_pANOVA_ecotype$coefficients["ecotypegliding"]
humerus_CSS_pANOVA_ecotype_ground_int <- humerus_CSS_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_CSS_pANOVA_ecotype$coefficients["ecotypeground"]
humerus_CSS_pANOVA_ecotype_tree_int <- humerus_CSS_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_CSS_pANOVA_ecotype$coefficients["ecotypetree"]
humerus_CSS_pANOVA_ecotype_chip_int_boot <- humerus_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"]
humerus_CSS_pANOVA_ecotype_glide_int_boot <- humerus_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_CSS_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
humerus_CSS_pANOVA_ecotype_ground_int_boot <- humerus_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_CSS_pANOVA_ecotype$bootstrap[,"ecotypeground"]
humerus_CSS_pANOVA_ecotype_tree_int_boot <- humerus_CSS_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_CSS_pANOVA_ecotype$bootstrap[,"ecotypetree"]

humerus_CSS_pANOVA_ecotype_int_cof <- rbind(
  humerus_CSS_pANOVA_ecotype_chip_int,
  humerus_CSS_pANOVA_ecotype_glide_int,
  humerus_CSS_pANOVA_ecotype_ground_int,
  humerus_CSS_pANOVA_ecotype_tree_int)
humerus_CSS_pANOVA_ecotype_int_95 <- rbind(
  quantile(humerus_CSS_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_CSS_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

humerus_CSS_pANOVA_ecotype_cof <- data.frame(humerus_CSS_pANOVA_ecotype_int_cof, humerus_CSS_pANOVA_ecotype_int_95)
colnames(humerus_CSS_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(humerus_CSS_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(humerus_CSS_pANOVA_ecotype_cof, digits = 3)

######   Plot humerus_CSS ~ ecotype ######  
humerus_CSS_vp_ecology <- ggplot(humerus_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= humerus_CSS, fill = ecotype)) + 
  scale_y_continuous(name="ln humerus_CSS (mm)", limits = c(1.05,2.6), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(humerus_CSS)) + 
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
humerus_CSS_vp_ecology


######  humerus_DE_int_meandiff_table ######  
# chip - glide
humerus_DE_int_chip_glide_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_chip_int_boot)^2/length(humerus_DE_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_glide_int_boot)^2/length(humerus_DE_pANOVA_ecotype_glide_int_boot)))

humerus_DE_int_chip_glide_meandiff <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_glide_int_boot)

humerus_DE_int_chip_glide_CI <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_chip_glide_SE  # 95% confidence interval using z=1.96

humerus_DE_int_chip_glide_table <- cbind(humerus_DE_int_chip_glide_meandiff, humerus_DE_int_chip_glide_CI[1], humerus_DE_int_chip_glide_CI[2])
rownames(humerus_DE_int_chip_glide_table) <- ("humerus_DE_int_chip_glide")
colnames(humerus_DE_int_chip_glide_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_chip_glide_table

# chip - ground
humerus_DE_int_chip_ground_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_chip_int_boot)^2/length(humerus_DE_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_ground_int_boot)^2/length(humerus_DE_pANOVA_ecotype_ground_int_boot)))

humerus_DE_int_chip_ground_meandiff <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot)

humerus_DE_int_chip_ground_CI <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_chip_ground_SE  # 95% confidence interval using z=1.96

humerus_DE_int_chip_ground_table <- cbind(humerus_DE_int_chip_ground_meandiff, humerus_DE_int_chip_ground_CI[1], humerus_DE_int_chip_ground_CI[2])
rownames(humerus_DE_int_chip_ground_table) <- ("humerus_DE_int_chip_ground")
colnames(humerus_DE_int_chip_ground_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_chip_ground_table

# chip - tree
humerus_DE_int_chip_tree_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_chip_int_boot)^2/length(humerus_DE_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_tree_int_boot)^2/length(humerus_DE_pANOVA_ecotype_tree_int_boot)))

humerus_DE_int_chip_tree_meandiff <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot)

humerus_DE_int_chip_tree_CI <- mean(humerus_DE_pANOVA_ecotype_chip_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_chip_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_int_chip_tree_table <- cbind(humerus_DE_int_chip_tree_meandiff, humerus_DE_int_chip_tree_CI[1], humerus_DE_int_chip_tree_CI[2])
rownames(humerus_DE_int_chip_tree_table) <- ("humerus_DE_int_chip_tree")
colnames(humerus_DE_int_chip_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_chip_tree_table

# glide - ground
humerus_DE_int_glide_ground_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_glide_int_boot)^2/length(humerus_DE_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_ground_int_boot)^2/length(humerus_DE_pANOVA_ecotype_ground_int_boot)))

humerus_DE_int_glide_ground_meandiff <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot)

humerus_DE_int_glide_ground_CI <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_glide_ground_SE  # 95% confidence interval using z=1.96

humerus_DE_int_glide_ground_table <- cbind(humerus_DE_int_glide_ground_meandiff, humerus_DE_int_glide_ground_CI[1], humerus_DE_int_glide_ground_CI[2])
rownames(humerus_DE_int_glide_ground_table) <- ("humerus_DE_int_glide_ground")
colnames(humerus_DE_int_glide_ground_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_glide_ground_table

# glide - tree
humerus_DE_int_glide_tree_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_glide_int_boot)^2/length(humerus_DE_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_tree_int_boot)^2/length(humerus_DE_pANOVA_ecotype_tree_int_boot)))

humerus_DE_int_glide_tree_meandiff <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot)

humerus_DE_int_glide_tree_CI <- mean(humerus_DE_pANOVA_ecotype_glide_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_glide_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_int_glide_tree_table <- cbind(humerus_DE_int_glide_tree_meandiff, humerus_DE_int_glide_tree_CI[1], humerus_DE_int_glide_tree_CI[2])
rownames(humerus_DE_int_glide_tree_table) <- ("humerus_DE_int_glide_tree")
colnames(humerus_DE_int_glide_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_glide_tree_table

# ground - tree
humerus_DE_int_ground_tree_SE <- sqrt((sd(humerus_DE_pANOVA_ecotype_ground_int_boot)^2/length(humerus_DE_pANOVA_ecotype_ground_int_boot)) + (sd(humerus_DE_pANOVA_ecotype_tree_int_boot)^2/length(humerus_DE_pANOVA_ecotype_tree_int_boot)))

humerus_DE_int_ground_tree_meandiff <- mean(humerus_DE_pANOVA_ecotype_ground_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot)

humerus_DE_int_ground_tree_CI <- mean(humerus_DE_pANOVA_ecotype_ground_int_boot) - mean(humerus_DE_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_DE_int_ground_tree_SE  # 95% confidence interval using z=1.96

humerus_DE_int_ground_tree_table <- cbind(humerus_DE_int_ground_tree_meandiff, humerus_DE_int_ground_tree_CI[1], humerus_DE_int_ground_tree_CI[2])
rownames(humerus_DE_int_ground_tree_table) <- ("humerus_DE_int_ground_tree")
colnames(humerus_DE_int_ground_tree_table) <- c("meandiff", "L95", "U95")
humerus_DE_int_ground_tree_table

#table
humerus_DE_int_meandiff_table <- rbind(
  humerus_DE_int_chip_glide_table,
  humerus_DE_int_chip_ground_table,
  humerus_DE_int_chip_tree_table,
  humerus_DE_int_glide_ground_table,
  humerus_DE_int_glide_tree_table,
  humerus_DE_int_ground_tree_table)
round(humerus_DE_int_meandiff_table, digits = 2)



humerus_internal <- cbind(humerus_Cg, humerus_DE, humerus_CSS)
humerus_internal_pca <- gm.prcomp(humerus_internal, phy = tree_pr, GLS = FALSE)
summary(humerus_internal_pca)
plot(humerus_internal_pca, phylo = TRUE, main = "phylo PCA")

col_ecotype_phylomorph <-c(col_ecotype[as.factor(ecotype)],rep("white",tree_pr$Nnode))
names(col_ecotype_phylomorph)<-1:(length(tree_pr$tip)+tree_pr$Nnode)
phylomorphospace(tree_pr, humerus_internal_pca$x[,1:2], control = list(col.node = col_ecotype_phylomorph), label = "off")

##### humerus_CSS AIC ####
geiger::aicw(c(AIC(humerus_CSS_pgls_size),
               AIC(humerus_CSS_pANOVA_ecotype),
               AIC(humerus_CSS_pANCOVA_size_ecotype),
               AIC(humerus_CSS_pANCOVAnoint_size_ecotype),
               AIC(humerus_CSS_pgls_null)))


#### humerus size ####
### humerus_size ANCOVA ###
humerus_size_pANOVA_ecotype <- phylolm(lnhumerus_size~ecotype, phy = tree_pr, model = "lambda", boot = 1000)
summary(humerus_size_pANOVA_ecotype)

humerus_size_pANOVA_ecotype_chip_int <- humerus_size_pANOVA_ecotype$coefficients["(Intercept)"]
humerus_size_pANOVA_ecotype_glide_int <- humerus_size_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_size_pANOVA_ecotype$coefficients["ecotypegliding"]
humerus_size_pANOVA_ecotype_ground_int <- humerus_size_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_size_pANOVA_ecotype$coefficients["ecotypeground"]
humerus_size_pANOVA_ecotype_tree_int <- humerus_size_pANOVA_ecotype$coefficients["(Intercept)"]+humerus_size_pANOVA_ecotype$coefficients["ecotypetree"]
humerus_size_pANOVA_ecotype_chip_int_boot <- humerus_size_pANOVA_ecotype$bootstrap[,"(Intercept)"]
humerus_size_pANOVA_ecotype_glide_int_boot <- humerus_size_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_size_pANOVA_ecotype$bootstrap[,"ecotypegliding"]
humerus_size_pANOVA_ecotype_ground_int_boot <- humerus_size_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_size_pANOVA_ecotype$bootstrap[,"ecotypeground"]
humerus_size_pANOVA_ecotype_tree_int_boot <- humerus_size_pANOVA_ecotype$bootstrap[,"(Intercept)"] + humerus_size_pANOVA_ecotype$bootstrap[,"ecotypetree"]

humerus_size_pANOVA_ecotype_int_cof <- rbind(
  humerus_size_pANOVA_ecotype_chip_int,
  humerus_size_pANOVA_ecotype_glide_int,
  humerus_size_pANOVA_ecotype_ground_int,
  humerus_size_pANOVA_ecotype_tree_int)
humerus_size_pANOVA_ecotype_int_95 <- rbind(
  quantile(humerus_size_pANOVA_ecotype_chip_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_size_pANOVA_ecotype_glide_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_size_pANOVA_ecotype_ground_int_boot, prob = c(0.025, 0.975)),
  quantile(humerus_size_pANOVA_ecotype_tree_int_boot, prob = c(0.025, 0.975)))

humerus_size_pANOVA_ecotype_cof <- data.frame(humerus_size_pANOVA_ecotype_int_cof, humerus_size_pANOVA_ecotype_int_95)
colnames(humerus_size_pANOVA_ecotype_cof) <- c("intercept", "Int_L95", "Int_U95")
rownames(humerus_size_pANOVA_ecotype_cof) <- c("chip", "gliding", "ground", "tree")
round(humerus_size_pANOVA_ecotype_cof, digits = 1)

humerus_size_vp_ecology <- ggplot(humerus_data, aes(x=fct_relevel(ecotype, "chip", "gliding", "ground", "tree"), y= lnhumerus_size, fill = ecotype)) + 
  scale_y_continuous(name="ln humerus_size (mm)", limits = c(4.2,6.1), breaks = scales::pretty_breaks(n = 8)) +
  xlab("") +
  geom_violin(trim=TRUE) + 
  scale_fill_manual(values= col_ecotype) + 
  geom_hline(yintercept=mean(lnhumerus_size)) + 
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
humerus_size_vp_ecology


##### humerus_size_int_meandiff_table ##### 
# chip - glide
humerus_size_int_chip_glide_SE <- sqrt((sd(humerus_size_pANOVA_ecotype_chip_int_boot)^2/length(humerus_size_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_size_pANOVA_ecotype_glide_int_boot)^2/length(humerus_size_pANOVA_ecotype_glide_int_boot)))

humerus_size_int_chip_glide_meandiff <- mean(humerus_size_pANOVA_ecotype_chip_int_boot) - mean(humerus_size_pANOVA_ecotype_glide_int_boot)

humerus_size_int_chip_glide_CI <- mean(humerus_size_pANOVA_ecotype_chip_int_boot) - mean(humerus_size_pANOVA_ecotype_glide_int_boot) + c(-1, 1) * 1.96 * humerus_size_int_chip_glide_SE  # 95% confidence interval using z=1.96

humerus_size_int_chip_glide_table <- cbind(humerus_size_int_chip_glide_meandiff, humerus_size_int_chip_glide_CI[1], humerus_size_int_chip_glide_CI[2])
rownames(humerus_size_int_chip_glide_table) <- ("humerus_size_int_chip_glide")
colnames(humerus_size_int_chip_glide_table) <- c("meandiff", "L95", "U95")
humerus_size_int_chip_glide_table

# chip - ground
humerus_size_int_chip_ground_SE <- sqrt((sd(humerus_size_pANOVA_ecotype_chip_int_boot)^2/length(humerus_size_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_size_pANOVA_ecotype_ground_int_boot)^2/length(humerus_size_pANOVA_ecotype_ground_int_boot)))

humerus_size_int_chip_ground_meandiff <- mean(humerus_size_pANOVA_ecotype_chip_int_boot) - mean(humerus_size_pANOVA_ecotype_ground_int_boot)

humerus_size_int_chip_ground_CI <- mean(humerus_size_pANOVA_ecotype_chip_int_boot) - mean(humerus_size_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_size_int_chip_ground_SE  # 95% confidence interval using z=1.96

humerus_size_int_chip_ground_table <- cbind(humerus_size_int_chip_ground_meandiff, humerus_size_int_chip_ground_CI[1], humerus_size_int_chip_ground_CI[2])
rownames(humerus_size_int_chip_ground_table) <- ("humerus_size_int_chip_ground")
colnames(humerus_size_int_chip_ground_table) <- c("meandiff", "L95", "U95")
humerus_size_int_chip_ground_table

# chip - tree
humerus_size_int_chip_tree_SE <- sqrt((sd(humerus_size_pANOVA_ecotype_chip_int_boot)^2/length(humerus_size_pANOVA_ecotype_chip_int_boot)) + (sd(humerus_size_pANOVA_ecotype_tree_int_boot)^2/length(humerus_size_pANOVA_ecotype_tree_int_boot)))

humerus_size_int_chip_tree_meandiff <- mean(humerus_size_pANOVA_ecotype_chip_int_boot) - mean(humerus_size_pANOVA_ecotype_tree_int_boot)

humerus_size_int_chip_tree_CI <- mean(humerus_size_pANOVA_ecotype_chip_int_boot) - mean(humerus_size_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_size_int_chip_tree_SE  # 95% confidence interval using z=1.96

humerus_size_int_chip_tree_table <- cbind(humerus_size_int_chip_tree_meandiff, humerus_size_int_chip_tree_CI[1], humerus_size_int_chip_tree_CI[2])
rownames(humerus_size_int_chip_tree_table) <- ("humerus_size_int_chip_tree")
colnames(humerus_size_int_chip_tree_table) <- c("meandiff", "L95", "U95")
humerus_size_int_chip_tree_table

# glide - ground
humerus_size_int_glide_ground_SE <- sqrt((sd(humerus_size_pANOVA_ecotype_glide_int_boot)^2/length(humerus_size_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_size_pANOVA_ecotype_ground_int_boot)^2/length(humerus_size_pANOVA_ecotype_ground_int_boot)))

humerus_size_int_glide_ground_meandiff <- mean(humerus_size_pANOVA_ecotype_glide_int_boot) - mean(humerus_size_pANOVA_ecotype_ground_int_boot)

humerus_size_int_glide_ground_CI <- mean(humerus_size_pANOVA_ecotype_glide_int_boot) - mean(humerus_size_pANOVA_ecotype_ground_int_boot) + c(-1, 1) * 1.96 * humerus_size_int_glide_ground_SE  # 95% confidence interval using z=1.96

humerus_size_int_glide_ground_table <- cbind(humerus_size_int_glide_ground_meandiff, humerus_size_int_glide_ground_CI[1], humerus_size_int_glide_ground_CI[2])
rownames(humerus_size_int_glide_ground_table) <- ("humerus_size_int_glide_ground")
colnames(humerus_size_int_glide_ground_table) <- c("meandiff", "L95", "U95")
humerus_size_int_glide_ground_table

# glide - tree
humerus_size_int_glide_tree_SE <- sqrt((sd(humerus_size_pANOVA_ecotype_glide_int_boot)^2/length(humerus_size_pANOVA_ecotype_glide_int_boot)) + (sd(humerus_size_pANOVA_ecotype_tree_int_boot)^2/length(humerus_size_pANOVA_ecotype_tree_int_boot)))

humerus_size_int_glide_tree_meandiff <- mean(humerus_size_pANOVA_ecotype_glide_int_boot) - mean(humerus_size_pANOVA_ecotype_tree_int_boot)

humerus_size_int_glide_tree_CI <- mean(humerus_size_pANOVA_ecotype_glide_int_boot) - mean(humerus_size_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_size_int_glide_tree_SE  # 95% confidence interval using z=1.96

humerus_size_int_glide_tree_table <- cbind(humerus_size_int_glide_tree_meandiff, humerus_size_int_glide_tree_CI[1], humerus_size_int_glide_tree_CI[2])
rownames(humerus_size_int_glide_tree_table) <- ("humerus_size_int_glide_tree")
colnames(humerus_size_int_glide_tree_table) <- c("meandiff", "L95", "U95")
humerus_size_int_glide_tree_table

# ground - tree
humerus_size_int_ground_tree_SE <- sqrt((sd(humerus_size_pANOVA_ecotype_ground_int_boot)^2/length(humerus_size_pANOVA_ecotype_ground_int_boot)) + (sd(humerus_size_pANOVA_ecotype_tree_int_boot)^2/length(humerus_size_pANOVA_ecotype_tree_int_boot)))

humerus_size_int_ground_tree_meandiff <- mean(humerus_size_pANOVA_ecotype_ground_int_boot) - mean(humerus_size_pANOVA_ecotype_tree_int_boot)

humerus_size_int_ground_tree_CI <- mean(humerus_size_pANOVA_ecotype_ground_int_boot) - mean(humerus_size_pANOVA_ecotype_tree_int_boot) + c(-1, 1) * 1.96 * humerus_size_int_ground_tree_SE  # 95% confidence interval using z=1.96

humerus_size_int_ground_tree_table <- cbind(humerus_size_int_ground_tree_meandiff, humerus_size_int_ground_tree_CI[1], humerus_size_int_ground_tree_CI[2])
rownames(humerus_size_int_ground_tree_table) <- ("humerus_size_int_ground_tree")
colnames(humerus_size_int_ground_tree_table) <- c("meandiff", "L95", "U95")
humerus_size_int_ground_tree_table

#table
humerus_size_int_meandiff_table <- rbind(
  humerus_size_int_chip_glide_table,
  humerus_size_int_chip_ground_table,
  humerus_size_int_chip_tree_table,
  humerus_size_int_glide_ground_table,
  humerus_size_int_glide_tree_table,
  humerus_size_int_ground_tree_table)
round(humerus_size_int_meandiff_table, digits = 2)



#### external shape vs internal structure ####
##### Cg #####
humerus_shape_Cg_pls <- two.b.pls(humerus_shape, humerus_Cg,  iter = 999)
summary(humerus_shape_Cg_pls)
plot(humerus_shape_Cg_pls, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)

#chipmunk
humerus_shape_Cg_pls_chipmunk <- two.b.pls(humerus_shape_chipmunk, humerus_Cg_chipmunk,  iter = 999)
summary(humerus_shape_Cg_pls_chipmunk)
plot(humerus_shape_Cg_pls_chipmunk)

#gliding
humerus_shape_Cg_pls_gliding <- two.b.pls(humerus_shape_gliding, humerus_Cg_gliding,  iter = 999)
plot(humerus_shape_Cg_pls_gliding)
summary(humerus_shape_Cg_pls_gliding)

#ground
humerus_shape_Cg_pls_ground <- two.b.pls(humerus_shape_ground, humerus_Cg_ground,  iter = 999)
plot(humerus_shape_Cg_pls_ground)
summary(humerus_shape_Cg_pls_ground)

#tree
humerus_shape_Cg_pls_tree <- two.b.pls(humerus_shape_tree, humerus_Cg_tree,  iter = 999)
plot(humerus_shape_Cg_pls_tree)
summary(humerus_shape_Cg_pls_tree)

humerus_shape_Cg_pls_table <- rbind(c(humerus_shape_Cg_pls$r.pls, humerus_shape_Cg_pls$Z, humerus_shape_Cg_pls$P.value),
                                    c(humerus_shape_Cg_pls_chipmunk$r.pls, humerus_shape_Cg_pls_chipmunk$Z, humerus_shape_Cg_pls_chipmunk$P.value),
                                    c(humerus_shape_Cg_pls_gliding$r.pls, humerus_shape_Cg_pls_gliding$Z, humerus_shape_Cg_pls_gliding$P.value),
                                    c(humerus_shape_Cg_pls_ground$r.pls, humerus_shape_Cg_pls_ground$Z, humerus_shape_Cg_pls_ground$P.value),
                                    c(humerus_shape_Cg_pls_tree$r.pls, humerus_shape_Cg_pls_tree$Z, humerus_shape_Cg_pls_tree$P.value))
rownames(humerus_shape_Cg_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(humerus_shape_Cg_pls_table) <- c("r", "Z", "P")
round(humerus_shape_Cg_pls_table, digits = 3)


##### DE #####
humerus_shape_DE_pls <- two.b.pls(humerus_shape, humerus_DE,  iter = 999)
summary(humerus_shape_DE_pls)
plot(humerus_shape_DE_pls, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)

#chipmunk
humerus_shape_DE_pls_chipmunk <- two.b.pls(humerus_shape_chipmunk, humerus_DE_chipmunk,  iter = 999)
summary(humerus_shape_DE_pls_chipmunk)
plot(humerus_shape_DE_pls_chipmunk)

#gliding
humerus_shape_DE_pls_gliding <- two.b.pls(humerus_shape_gliding, humerus_DE_gliding,  iter = 999)
plot(humerus_shape_DE_pls_gliding)
summary(humerus_shape_DE_pls_gliding)

#ground
humerus_shape_DE_pls_ground <- two.b.pls(humerus_shape_ground, humerus_DE_ground,  iter = 999)
plot(humerus_shape_DE_pls_ground)
summary(humerus_shape_DE_pls_ground)

#tree
humerus_shape_DE_pls_tree <- two.b.pls(humerus_shape_tree, humerus_DE_tree,  iter = 999)
plot(humerus_shape_DE_pls_tree)
summary(humerus_shape_DE_pls_tree)

humerus_shape_DE_pls_table <- rbind(c(humerus_shape_DE_pls$r.pls, humerus_shape_DE_pls$Z, humerus_shape_DE_pls$P.value),
                                    c(humerus_shape_DE_pls_chipmunk$r.pls, humerus_shape_DE_pls_chipmunk$Z, humerus_shape_DE_pls_chipmunk$P.value),
                                    c(humerus_shape_DE_pls_gliding$r.pls, humerus_shape_DE_pls_gliding$Z, humerus_shape_DE_pls_gliding$P.value),
                                    c(humerus_shape_DE_pls_ground$r.pls, humerus_shape_DE_pls_ground$Z, humerus_shape_DE_pls_ground$P.value),
                                    c(humerus_shape_DE_pls_tree$r.pls, humerus_shape_DE_pls_tree$Z, humerus_shape_DE_pls_tree$P.value))
rownames(humerus_shape_DE_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(humerus_shape_DE_pls_table) <- c("r", "Z", "P")
round(humerus_shape_DE_pls_table, digits = 3)


##### CSS #####
humerus_shape_CSS_pls <- two.b.pls(humerus_shape, humerus_CSS,  iter = 999)
summary(humerus_shape_CSS_pls)
plot(humerus_shape_CSS_pls, pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5)

#chipmunk
humerus_shape_CSS_pls_chipmunk <- two.b.pls(humerus_shape_chipmunk, humerus_CSS_chipmunk,  iter = 999)
summary(humerus_shape_CSS_pls_chipmunk)
plot(humerus_shape_CSS_pls_chipmunk)

#gliding
humerus_shape_CSS_pls_gliding <- two.b.pls(humerus_shape_gliding, humerus_CSS_gliding,  iter = 999)
plot(humerus_shape_CSS_pls_gliding)
summary(humerus_shape_CSS_pls_gliding)

#ground
humerus_shape_CSS_pls_ground <- two.b.pls(humerus_shape_ground, humerus_CSS_ground,  iter = 999)
plot(humerus_shape_CSS_pls_ground)
summary(humerus_shape_CSS_pls_ground)

#tree
humerus_shape_CSS_pls_tree <- two.b.pls(humerus_shape_tree, humerus_CSS_tree,  iter = 999)
plot(humerus_shape_CSS_pls_tree)
summary(humerus_shape_CSS_pls_tree)

humerus_shape_CSS_pls_table <- rbind(c(humerus_shape_CSS_pls$r.pls, humerus_shape_CSS_pls$Z, humerus_shape_CSS_pls$P.value),
                                     c(humerus_shape_CSS_pls_chipmunk$r.pls, humerus_shape_CSS_pls_chipmunk$Z, humerus_shape_CSS_pls_chipmunk$P.value),
                                     c(humerus_shape_CSS_pls_gliding$r.pls, humerus_shape_CSS_pls_gliding$Z, humerus_shape_CSS_pls_gliding$P.value),
                                     c(humerus_shape_CSS_pls_ground$r.pls, humerus_shape_CSS_pls_ground$Z, humerus_shape_CSS_pls_ground$P.value),
                                     c(humerus_shape_CSS_pls_tree$r.pls, humerus_shape_CSS_pls_tree$Z, humerus_shape_CSS_pls_tree$P.value))
rownames(humerus_shape_CSS_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(humerus_shape_CSS_pls_table) <- c("r", "Z", "P")
round(humerus_shape_CSS_pls_table, digits = 3)
