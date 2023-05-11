library(phytools)
library(geiger)
library(phylolm)
library(ggplot2)
library(forcats)
library(geomorph)
library(dplyr)
library(mvMORPH)

load("~/Documents/Projects/Project_Squirrels/squirrel_limbs/Rdata/humerus_data.RData")
load("~/Documents/Projects/Project_Squirrels/squirrel_limbs/Rdata/femur_data.RData")

#### shape - humerus vs femur ####
limb_shape_pls <- two.b.pls(humerus_shape, femur_shape, iter = 999)
plot(limb_shape_pls$XScores[,1], limb_shape_pls$YScores[,1], pch = c(24, 21, 22, 23)[factor(ecotype)], bg = col_ecotype[factor(ecotype)], las =1, cex = 1.75, cex.axis=1.5, xlab = "humeral shape", ylab = "femoral shape")
summary(limb_shape_pls)

limb_shape_pls_chipmunk <- two.b.pls(humerus_shape_chipmunk, femur_shape_chipmunk,  iter = 999)
plot(limb_shape_pls_chipmunk)
summary(limb_shape_pls)

limb_shape_pls_gliding <- two.b.pls(humerus_shape_gliding, femur_shape_gliding,  iter = 999)
plot(limb_shape_pls_gliding)
summary(limb_shape_pls)

limb_shape_pls_ground <- two.b.pls(humerus_shape_ground, femur_shape_ground, iter = 999)
plot(limb_shape_pls_ground)
summary(limb_shape_pls)

limb_shape_pls_tree <- two.b.pls(humerus_shape_tree, femur_shape_tree, iter = 999)
plot(limb_shape_pls_tree)
summary(limb_shape_pls_tree)

shape_pls_table <- rbind(c(limb_shape_pls$r.pls, limb_shape_pls$Z, limb_shape_pls$P.value),
                         c(limb_shape_pls_chipmunk$r.pls, limb_shape_pls_chipmunk$Z, limb_shape_pls_chipmunk$P.value),
                         c(limb_shape_pls_gliding$r.pls, limb_shape_pls_gliding$Z, limb_shape_pls_gliding$P.value),
                         c(limb_shape_pls_ground$r.pls, limb_shape_pls_ground$Z, limb_shape_pls_ground$P.value),
                         c(limb_shape_pls_tree$r.pls, limb_shape_pls_tree$Z, limb_shape_pls_tree$P.value))
rownames(shape_pls_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(shape_pls_table) <- c("r", "Z", "P")
round(shape_pls_table, digits = 3)

                         
#### Cg phy t test ####
Cg_phyt.test <- phyl.pairedttest(tree_pr, humerus_Cg, femur_Cg, lambda = 0, h0 = 0)
Cg_phyt.test

Cg_phyt.test_chipmunk <- phyl.pairedttest(tree_pr_chipmunk, humerus_chipmunk_Cg, femur_chipmunk_Cg, lambda = 0, h0 = 0)
Cg_phyt.test_chipmunk

Cg_phyt.test_gliding <- phyl.pairedttest(tree_pr_gliding, humerus_gliding_Cg, femur_gliding_Cg, lambda = 0, h0 = 0)
Cg_phyt.test_gliding

Cg_phyt.test_ground <- phyl.pairedttest(tree_pr_ground, humerus_ground_Cg, femur_ground_Cg, lambda = 0, h0 = 0)
Cg_phyt.test_ground

Cg_phyt.test_tree <- phyl.pairedttest(tree_pr_tree, humerus_tree_Cg, femur_tree_Cg, lambda = 0, h0 = 0)
Cg_phyt.test_tree

Cg_phyt.test_table <- rbind(c(mean(humerus_Cg), mean(femur_Cg), Cg_phyt.test$dbar, Cg_phyt.test$P.dbar, Cg_phyt.test$lambda),
  c(mean(humerus_chipmunk_Cg), mean(femur_chipmunk_Cg), Cg_phyt.test_chipmunk$dbar, Cg_phyt.test_chipmunk$P.dbar, Cg_phyt.test_chipmunk$lambda),
  c(mean(humerus_gliding_Cg), mean(femur_gliding_Cg), Cg_phyt.test_gliding$dbar, Cg_phyt.test_gliding$P.dbar, Cg_phyt.test_gliding$lambda),
  c(mean(humerus_ground_Cg), mean(femur_ground_Cg), Cg_phyt.test_ground$dbar, Cg_phyt.test_ground$P.dbar, Cg_phyt.test_ground$lambda),
  c(mean(humerus_tree_Cg), mean(femur_tree_Cg), Cg_phyt.test_tree$dbar, Cg_phyt.test_tree$P.dbar, Cg_phyt.test_tree$lambda)
    )
rownames(Cg_phyt.test_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(Cg_phyt.test_table) <- c("humerus Cg", "femur Cg", "phylogenetic difference", "P", "lambda")
round(Cg_phyt.test_table, digits = 3)


#### DE phy t test ####
DE_phyt.test <- phyl.pairedttest(tree_pr, humerus_DE, femur_DE, lambda = 0, h0 = 0)
DE_phyt.test

DE_phyt.test_chipmunk <- phyl.pairedttest(tree_pr_chipmunk, humerus_chipmunk_DE, femur_chipmunk_DE, lambda = 0, h0 = 0)
DE_phyt.test_chipmunk

DE_phyt.test_gliding <- phyl.pairedttest(tree_pr_gliding, humerus_gliding_DE, femur_gliding_DE, lambda = 0, h0 = 0)
DE_phyt.test_gliding

DE_phyt.test_ground <- phyl.pairedttest(tree_pr_ground, humerus_ground_DE, femur_ground_DE, lambda = 0, h0 = 0)
DE_phyt.test_ground

DE_phyt.test_tree <- phyl.pairedttest(tree_pr_tree, humerus_tree_DE, femur_tree_DE, lambda = 0, h0 = 0)
DE_phyt.test_tree

DE_phyt.test_table <- rbind(c(mean(humerus_DE), mean(femur_DE), DE_phyt.test$dbar, DE_phyt.test$P.dbar, DE_phyt.test$lambda),
                            c(mean(humerus_chipmunk_DE), mean(femur_chipmunk_DE), DE_phyt.test_chipmunk$dbar, DE_phyt.test_chipmunk$P.dbar, DE_phyt.test_chipmunk$lambda),
                            c(mean(humerus_gliding_DE), mean(femur_gliding_DE), DE_phyt.test_gliding$dbar, DE_phyt.test_gliding$P.dbar, DE_phyt.test_gliding$lambda),
                            c(mean(humerus_ground_DE), mean(femur_ground_DE), DE_phyt.test_ground$dbar, DE_phyt.test_ground$P.dbar, DE_phyt.test_ground$lambda),
                            c(mean(humerus_tree_DE), mean(femur_tree_DE), DE_phyt.test_tree$dbar, DE_phyt.test_tree$P.dbar, DE_phyt.test_tree$lambda)
)
rownames(DE_phyt.test_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(DE_phyt.test_table) <- c("humerus DE", "femur DE", "phylogenetic difference", "P", "lambda")
round(DE_phyt.test_table, digits = 3)


#### CSS phy t test ####
CSS_phyt.test <- phyl.pairedttest(tree_pr, humerus_CSS, femur_CSS, lambda = 0, h0 = 0)
CSS_phyt.test

CSS_phyt.test_chipmunk <- phyl.pairedttest(tree_pr_chipmunk, humerus_chipmunk_CSS, femur_chipmunk_CSS, lambda = 0, h0 = 0)
CSS_phyt.test_chipmunk

CSS_phyt.test_gliding <- phyl.pairedttest(tree_pr_gliding, humerus_gliding_CSS, femur_gliding_CSS, lambda = 0, h0 = 0)
CSS_phyt.test_gliding

CSS_phyt.test_ground <- phyl.pairedttest(tree_pr_ground, humerus_ground_CSS, femur_ground_CSS, lambda = 0, h0 = 0)
CSS_phyt.test_ground

CSS_phyt.test_tree <- phyl.pairedttest(tree_pr_tree, humerus_tree_CSS, femur_tree_CSS, lambda = 0, h0 = 0)
CSS_phyt.test_tree

CSS_phyt.test_table <- rbind(c(mean(humerus_CSS), mean(femur_CSS), CSS_phyt.test$dbar, CSS_phyt.test$P.dbar, CSS_phyt.test$lambda),
                             c(mean(humerus_chipmunk_CSS), mean(femur_chipmunk_CSS), CSS_phyt.test_chipmunk$dbar, CSS_phyt.test_chipmunk$P.dbar, CSS_phyt.test_chipmunk$lambda),
                             c(mean(humerus_gliding_CSS), mean(femur_gliding_CSS), CSS_phyt.test_gliding$dbar, CSS_phyt.test_gliding$P.dbar, CSS_phyt.test_gliding$lambda),
                             c(mean(humerus_ground_CSS), mean(femur_ground_CSS), CSS_phyt.test_ground$dbar, CSS_phyt.test_ground$P.dbar, CSS_phyt.test_ground$lambda),
                             c(mean(humerus_tree_CSS), mean(femur_tree_CSS), CSS_phyt.test_tree$dbar, CSS_phyt.test_tree$P.dbar, CSS_phyt.test_tree$lambda)
)
rownames(CSS_phyt.test_table) <- c("all", "chipmunk", "gliding", "ground", "tree")
colnames(CSS_phyt.test_table) <- c("humerus CSS", "femur CSS", "phylogenetic difference", "P", "lambda")
round(CSS_phyt.test_table, digits = 3)






