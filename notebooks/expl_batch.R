library(knitr)
library(xtable) # table
library(mixOmics)
library(sva) # ComBat
library(ggplot2) # PCA sample plot with density
library(gridExtra) # PCA sample plot with density
library(limma) # removeBatchEffect (LIMMA)
library(vegan) # RDA
library(AgiMicroRna) # RLE plot
library(cluster) # silhouette coefficient
library(variancePartition) # variance calculation
library(pvca) # PVCA
library(pheatmap) # heatmap
library(ruv) # RUVIII
library(lmerTest) # lmer
library(bapred)
library(FactoMineR)# FAbatch
source('/Users/basilews/phylo/helpers.R')
source('/Users/basilews/phylo/batc_panc/DATA/exp/helpers.R')

load(file = '/Users/basilews/phylo/batc_panc/DATA/exp/microbiome_datasets.RData')
sponge.trt
#Prefiltering
sponge.tss <- sponge.tss + 0.01
typeof(sponge.tss)
class(sponge.tss) 
#Centered log-ratio transformation Microbiome data are compostional and with different library sizes. Using standard statistical methods on such data may lead to spurious results and therefore the data must be further transformed. The CLR is the transformation of choice.
sponge.tss.clr <- logratio.transfo(sponge.tss, logratio = 'CLR')
#Batch effect detection
#Principal component analysis (PCA) with density plot per component
class(sponge.tss.clr) <- 'matrix' 
sponge.pca.before <- pca(sponge.tss.clr, ncomp = 3)
Scatter_Density(data = sponge.pca.before$variates$X, batch = sponge.batch, 
                trt = sponge.trt, expl.var = sponge.pca.before$explained_variance, 
                xlim = c(-4.5,5), ylim = c(-3,4), 
                batch.legend.title = 'Gel (batch)', 
                trt.legend.title = 'Tissue (trt)', 
                title = 'Before batch effect correction (Sponge)')
#Density plot and box plot
sponge.before.df <- data.frame(value = sponge.tss.clr[,9], batch = sponge.batch)

box_plot_fun(data = sponge.before.df, x = sponge.before.df$batch,
             y = sponge.before.df$value, title = 'OTU9 (Sponge)',
             batch.legend.title = 'Gel (batch)')
ggplot(sponge.before.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())

sponge.lm <- lm(sponge.tss.clr[,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm)

#RLE plots
sponge.batch_c <- sponge.batch[sponge.trt == 'C']
sponge.batch_e <- sponge.batch[sponge.trt == 'E'] 

sponge.before_c <- sponge.tss.clr[sponge.trt == 'C', ]
sponge.before_e <- sponge.tss.clr[sponge.trt == 'E', ] 


RleMicroRna2(object = t(sponge.before_c), batch = sponge.batch_c, 
             maintitle = 'Sponge (tissue: choanosome)')

RleMicroRna2(object = t(sponge.before_e), batch = sponge.batch_e, 
             maintitle = 'Sponge (tissue: ectosome)')
class(sponge.batch)

# Heatmap
# Sponge data
# scale on OTUs
sponge.tss.clr.scale <- scale(sponge.tss.clr, center = T, scale = T) 
# scale on samples
sponge.tss.clr.scale <- scale(t(sponge.tss.clr.scale), center = T, scale = T) 

sponge.anno_col <- data.frame(Batch = sponge.batch, Tissue = sponge.trt)
sponge.anno_metabo_colors <- list(Batch = c('1' = '#388ECC', '2' = '#F68B33'), 
                                  Tissue = c(C = '#F0E442', E = '#D55E00'))

class(sponge.batch)
class(sponge.trt)
typeof(sponge.batch)
typeof(sponge.trt)
sponge.trt
pheatmap(sponge.tss.clr.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - Scaled')
dev.off() 
#Accounting for batch effects
#Linear model and linear mixed model
sponge.trt_p <- apply(sponge.tss.clr, 2, FUN = function(x){
  res.lm <- lm(x ~ sponge.trt + sponge.batch)
  summary.res <- summary(res.lm)
  p <- summary.res$coefficients[2,4]
})

sponge.trt_adjp <- p.adjust(sponge.trt_p, method = 'fdr')

#SVA
# sponge data
sponge.mod <- model.matrix( ~ sponge.trt) # full model
sponge.mod0 <- model.matrix( ~ 1,data = sponge.trt) # null model
sponge.sva.n <- num.sv(dat = t(sponge.tss.clr), mod = sponge.mod)
sponge.sva <- sva(dat = t(sponge.tss.clr), mod = sponge.mod, 
                  mod0 = sponge.mod0, n.sv = sponge.sva.n)
sponge.mod.bat <- cbind(sponge.mod, sponge.sva$sv)
sponge.mod0.bat <- cbind(sponge.mod0, sponge.sva$sv)

sponge.sva.trt_p <- f.pvalue(t(sponge.tss.clr), sponge.mod.bat, sponge.mod0.bat)
sponge.sva.trt_adjp <- p.adjust(sponge.sva.trt_p, method='fdr')

#RUV2
#This analysis is not optimal as we use pseudo negative control variables, but we wish to provide this as a practical example
sponge.nc <- sponge.trt_adjp > 0.05
sponge.ruv2 <- RUV2(Y = sponge.tss.clr, X = sponge.trt,ctl = sponge.nc, k = 3) # k is subjective
sponge.ruv2.trt_p <- sponge.ruv2$p
sponge.ruv2.trt_adjp <- p.adjust(sponge.ruv2.trt_p, method="fdr")

#RUV4
sponge.k.obj <- getK(Y = sponge.tss.clr, X = sponge.trt, ctl = sponge.nc)
sponge.k <- sponge.k.obj$k
sponge.k <- ifelse(sponge.k !=0, sponge.k, 1)
sponge.ruv4 <- RUV4(Y = sponge.tss.clr, X = sponge.trt, ctl = sponge.nc, k = sponge.k) 
sponge.ruv4.trt_p <- sponge.ruv4$p
sponge.ruv4.trt_adjp <- p.adjust(sponge.ruv4.trt_p, method="fdr")


#Correcting for batch effects
#BMC (batch mean centering)
# Sponge data
sponge.b1 <- scale(sponge.tss.clr[sponge.batch == 1, ], center = TRUE, scale = FALSE)
sponge.b2 <- scale(sponge.tss.clr[sponge.batch == 2, ], center = TRUE, scale = FALSE)
sponge.bmc <- rbind(sponge.b1, sponge.b2)
sponge.bmc <- sponge.bmc[rownames(sponge.tss.clr), ]

#ComBat
sponge.combat <- t(ComBat(t(sponge.tss.clr), batch = sponge.batch, 
                          mod = sponge.mod, par.prior = F, prior.plots = F))

#removeBatchEffect
# Sponge data
sponge.limma <- t(removeBatchEffect(t(sponge.tss.clr), batch = sponge.batch, 
                                    design = sponge.mod))
#FAbatch
sponge.fabatch.obj <- fabatch(x = sponge.tss.clr, 
                              y = as.factor(as.numeric(sponge.trt)), 
                              batch = sponge.batch)
sponge.fabatch <- sponge.fabatch.obj$xadj

#Percentile normalisation
sponge.tss <- t(apply(sponge.tss, 1, function(x){x/sum(x)}))
sponge.percentile <- percentile_norm(data = sponge.tss, batch = sponge.batch, 
                                     trt = sponge.trt)
#SVD-based method
# sponge data
sponge.sd <- apply(sponge.tss.clr, 2, sd) # calculate standard deviation
sponge.mean <- apply(sponge.tss.clr, 2, mean) # calculate mean
sponge.X <- scale(sponge.tss.clr, center = T, scale = T) # center and scale

sponge.m <- crossprod(sponge.X) # generate a square matrix
sponge.m.svd <- svd(sponge.m) # SVD 
# extract 1st singular vectors
sponge.a1 <- sponge.m.svd$u[ ,1] 
sponge.b1 <- sponge.m.svd$v[ ,1]

# deflate component 1 from the data
sponge.t1 <- sponge.X %*% sponge.a1 / drop(sqrt(crossprod(sponge.a1)))
sponge.c1 <- crossprod(sponge.X, sponge.t1) / drop(crossprod(sponge.t1))
sponge.svd.defl.matrix1  <- sponge.X - sponge.t1 %*% t(sponge.c1)

# add back mean and standard deviation
sponge.svd <- sponge.svd.defl.matrix1
sponge.svd[1:nrow(sponge.svd), 1:ncol(sponge.svd)] = NA
for(i in 1:ncol(sponge.svd.defl.matrix1)){
  for(j in 1:nrow(sponge.svd.defl.matrix1)){
    sponge.svd[j,i] = sponge.svd.defl.matrix1[j,i]*sponge.sd[i] + sponge.mean[i]
  }
}

#Methods evaluation
#Principal component analysis (PCA) with density plot per component
sponge.pca.before <- pca(sponge.tss.clr, ncomp = 3)
sponge.pca.bmc <- pca(sponge.bmc, ncomp = 3)
sponge.pca.combat <- pca(sponge.combat, ncomp = 3)
sponge.pca.limma <- pca(sponge.limma, ncomp = 3)
sponge.pca.percentile <- pca(sponge.percentile, ncomp = 3)
sponge.pca.svd <- pca(sponge.svd, ncomp = 3)
sponge.pca.plot.before = Scatter_Density(data = sponge.pca.before$variates$X, batch = sponge.batch, 
                                         trt = sponge.trt, expl.var = sponge.pca.before$explained_variance, 
                                         xlim = c(-4.5,5), ylim = c(-3,4), 
                                         batch.legend.title = 'Gel (batch)', 
                                         trt.legend.title = 'Tissue (trt)', 
                                         title = 'Before batch effect correction (Sponge)')
sponge.pca.plot.bmc = Scatter_Density(data = sponge.pca.bmc$variates$X, batch = sponge.batch, 
                                      trt = sponge.trt, expl.var = sponge.pca.bmc$explained_variance, 
                                      xlim = c(-4.5,5), ylim = c(-3,4), 
                                      batch.legend.title = 'Gel (batch)', 
                                      trt.legend.title = 'Tissue (trt)', 
                                      title = 'bmc (Sponge)')
sponge.pca.plot.combat = Scatter_Density(data = sponge.pca.combat$variates$X, batch = sponge.batch, 
                                         trt = sponge.trt, expl.var = sponge.pca.combat$explained_variance, 
                                         xlim = c(-4.5,5), ylim = c(-3,4), 
                                         batch.legend.title = 'Gel (batch)', 
                                         trt.legend.title = 'Tissue (trt)', 
                                         title = 'combat (Sponge)')
sponge.pca.plot.limma = Scatter_Density(data = sponge.pca.limma$variates$X, batch = sponge.batch, 
                                        trt = sponge.trt, expl.var = sponge.pca.limma$explained_variance, 
                                        xlim = c(-4.5,5), ylim = c(-3,4), 
                                        batch.legend.title = 'Gel (batch)', 
                                        trt.legend.title = 'Tissue (trt)', 
                                        title = 'Limma (Sponge)')
grid.arrange(sponge.pca.plot.before, sponge.pca.plot.bmc, 
             sponge.pca.plot.combat, sponge.pca.plot.limma, ncol = 2)
sponge.pca.plot.percentile = Scatter_Density(data = sponge.pca.percentile$variates$X, batch = sponge.batch, 
                                             trt = sponge.trt, expl.var = sponge.pca.percentile$explained_variance, 
                                             xlim = c(-4.5,5), ylim = c(-3,4), 
                                             batch.legend.title = 'Gel (batch)', 
                                             trt.legend.title = 'Tissue (trt)', 
                                             title = 'Percentile (Sponge)')
sponge.pca.plot.svd = Scatter_Density(data = sponge.pca.svd$variates$X, batch = sponge.batch, 
                                      trt = sponge.trt, expl.var = sponge.pca.svd$explained_variance, 
                                      xlim = c(-4.5,5), ylim = c(-3,4), 
                                      batch.legend.title = 'Gel (batch)', 
                                      trt.legend.title = 'Tissue (trt)', 
                                      title = 'SVD (Sponge)')
grid.arrange(sponge.pca.plot.before, sponge.pca.plot.percentile, 
             sponge.pca.plot.svd, ncol = 2)

#sponge data
sponge.before.df <- data.frame(value = sponge.tss.clr[ ,9], batch = sponge.batch)
sponge.boxplot.before <- box_plot_fun(data = sponge.before.df, 
                                      x = sponge.before.df$batch,
                                      y = sponge.before.df$value, 
                                      title = 'OTU9 - before (Sponge)', 
                                      batch.legend.title = 'Gel (batch)')

sponge.bmc.df <- data.frame(value = sponge.bmc[ ,9], batch = sponge.batch)
sponge.boxplot.bmc <- box_plot_fun(data = sponge.bmc.df, 
                                   x = sponge.bmc.df$batch,
                                   y = sponge.bmc.df$value, 
                                   title = 'OTU9 - BMC (Sponge)', 
                                   batch.legend.title = 'Gel (batch)')


sponge.combat.df <- data.frame(value = sponge.combat[ ,9], batch = sponge.batch)
sponge.boxplot.combat <- box_plot_fun(data = sponge.combat.df, 
                                      x = sponge.combat.df$batch, 
                                      y = sponge.combat.df$value, 
                                      title = 'OTU9 - ComBat (Sponge)',
                                      batch.legend.title = 'Gel (batch)')


sponge.limma.df <- data.frame(value = sponge.limma[ ,9], batch = sponge.batch)
sponge.boxplot.limma <- box_plot_fun(data = sponge.limma.df, 
                                     x = sponge.limma.df$batch, 
                                     y = sponge.limma.df$value, 
                                     title = 'OTU9 - rBE(Sponge)', 
                                     batch.legend.title = 'Gel (batch)')

sponge.percentile.df <- data.frame(value = sponge.percentile[ ,9], batch = sponge.batch)
sponge.boxplot.percentile <- box_plot_fun(data = sponge.percentile.df, 
                                          x = sponge.percentile.df$batch, 
                                          y = sponge.percentile.df$value, 
                                          title = 'OTU9 - PN (Sponge)', 
                                          batch.legend.title = 'Gel (batch)')

sponge.svd.df <- data.frame(value = sponge.svd[ ,9], batch = sponge.batch)
sponge.boxplot.svd <- box_plot_fun(data = sponge.svd.df, 
                                   x = sponge.svd.df$batch, 
                                   y = sponge.svd.df$value, 
                                   title = 'OTU9 - SVD (Sponge)', 
                                   batch.legend.title = 'Gel (batch)')

grid.arrange(sponge.boxplot.before, sponge.boxplot.bmc, 
             sponge.boxplot.combat, sponge.boxplot.limma, ncol = 2)
grid.arrange(sponge.boxplot.before, sponge.boxplot.percentile, 
             sponge.boxplot.svd, ncol = 2)
# density plot
# before
sponge.dens.before <- ggplot(sponge.before.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 - before (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())

# BMC
sponge.dens.bmc <- ggplot(sponge.bmc.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 - BMC (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())


# ComBat
sponge.dens.combat <- ggplot(sponge.combat.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 - ComBat (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())


# removeBatchEffect 
sponge.dens.limma <- ggplot(sponge.limma.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 - rBE (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())


# percentile normal
sponge.dens.percentile <- ggplot(sponge.percentile.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 - PN (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())


# SVD
sponge.dens.svd <- ggplot(sponge.svd.df, aes(x = value, fill = batch)) + 
  geom_density(alpha = 0.5) + scale_fill_manual(values = color.mixo(1:10)) + 
  labs(title = 'OTU9 - SVD (Sponge)', x = 'Value', fill = 'Gel (batch)') + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5), 
                     panel.grid = element_blank())
grid.arrange(sponge.dens.before, sponge.dens.bmc, 
             sponge.dens.combat, sponge.dens.limma, ncol = 2)
grid.arrange(sponge.dens.before, sponge.dens.percentile, 
             sponge.dens.svd, ncol = 2)
# P-values
sponge.lm.before <- lm(sponge.tss.clr[ ,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm.before)
sponge.lm.bmc <- lm(sponge.bmc[ ,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm.bmc)
sponge.lm.combat <- lm(sponge.combat[ ,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm.combat)

sponge.lm.limma <- lm(sponge.limma[ ,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm.limma)

sponge.lm.percentile <- lm(sponge.percentile[ ,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm.percentile)

sponge.lm.svd <- lm(sponge.svd[ ,9] ~ sponge.trt + sponge.batch)
summary(sponge.lm.svd)
# RLE plots
# sponge data
# BMC
sponge.bmc_c <- sponge.bmc[sponge.trt == 'C', ]
sponge.bmc_e <- sponge.bmc[sponge.trt == 'E', ] 

# ComBat
sponge.combat_c <- sponge.combat[sponge.trt == 'C', ]
sponge.combat_e <- sponge.combat[sponge.trt == 'E', ] 

# rBE
sponge.limma_c <- sponge.limma[sponge.trt == 'C', ]
sponge.limma_e <- sponge.limma[sponge.trt == 'E', ] 

# PN
sponge.percentile_c <- sponge.percentile[sponge.trt == 'C', ]
sponge.percentile_e <- sponge.percentile[sponge.trt == 'E', ] 

# SVD
sponge.svd_c <- sponge.svd[sponge.trt == 'C', ]
sponge.svd_e <- sponge.svd[sponge.trt == 'E', ] 
dev.off() 
par(mfrow = c(2,3), mai = c(0.4,0.6,0.3,0.1))

RleMicroRna2(object = t(sponge.before_c), batch = sponge.batch_c, 
             maintitle = 'Sponge: before (choanosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.bmc_c), batch = sponge.batch_c, 
             maintitle = 'Sponge: BMC (choanosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.combat_c), batch = sponge.batch_c, 
             maintitle = 'Sponge: ComBat (choanosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.limma_c), batch = sponge.batch_c, 
             maintitle = 'Sponge: rBE (choanosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.percentile_c), batch = sponge.batch_c, 
             maintitle = 'Sponge: PN (choanosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.svd_c), batch = sponge.batch_c, 
             maintitle = 'Sponge: SVD (choanosome)', title.cex = 1)

par(mfrow = c(1,1))
par(mfrow = c(2,3), mai = c(0.4,0.6,0.3,0.1))

RleMicroRna2(object = t(sponge.before_e), batch = sponge.batch_e, 
             maintitle = 'Sponge: before (ectosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.bmc_e), batch = sponge.batch_e, 
             maintitle = 'Sponge: BMC (ectosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.combat_e), batch = sponge.batch_e, 
             maintitle = 'Sponge: ComBat (ectosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.limma_e), batch = sponge.batch_e, 
             maintitle = 'Sponge: rBE (ectosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.percentile_e), batch = sponge.batch_e, 
             maintitle = 'Sponge: PN (ectosome)', title.cex = 1)

RleMicroRna2(object = t(sponge.svd_e), batch = sponge.batch_e, 
             maintitle = 'Sponge: SVD (ectosome)', title.cex = 1)
#Heatmap
# Sponge data
# before 
sponge.tss.clr.scale <- scale(sponge.tss.clr, center = T, scale = T) 
# scale on OTUs
sponge.tss.clr.scale <- scale(t(sponge.tss.clr.scale), center = T, scale = T) 
# scale on samples

sponge.anno_col <- data.frame(Batch = sponge.batch, Tissue = sponge.trt)
sponge.anno_metabo_colors <- list(Batch = c('1' = '#388ECC', '2' = '#F68B33'), 
                                  Tissue = c(C = '#F0E442', E = '#D55E00'))

dev.off()
pheatmap(sponge.tss.clr.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - before')
# ComBat 
sponge.combat.scale <- scale(sponge.combat, center = T, scale = T) 
# scale on OTUs
sponge.combat.scale <- scale(t(sponge.combat.scale), center = T, scale = T) 
# scale on samples

pheatmap(sponge.combat.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - ComBat')
# BMC 
sponge.bmc.scale <- scale(sponge.bmc, center = T, scale = T) 
# scale on OTUs
sponge.bmc.scale <- scale(t(sponge.bmc.scale), center = T, scale = T) 
# scale on samples

pheatmap(sponge.bmc.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - BMC')


# removeBatchEffect
sponge.limma.scale <- scale(sponge.limma, center = T, scale = T) 
# scale on OTUs
sponge.limma.scale <- scale(t(sponge.limma.scale), center = T, scale = T) 
# scale on samples

pheatmap(sponge.limma.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - removeBatchEffect')

# percentile normalisation
sponge.percentile.scale <- scale(sponge.percentile, center = T, scale = T) 
# scale on OTUs
sponge.percentile.scale <- scale(t(sponge.percentile.scale), center = T, scale = T) 
# scale on samples

pheatmap(sponge.percentile.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - percentile norm')


# SVD
sponge.svd.scale <- scale(sponge.svd, center = T, scale = T) 
# scale on OTUs
sponge.svd.scale <- scale(t(sponge.svd.scale), center = T, scale = T) 
# scale on samples

pheatmap(sponge.svd.scale, 
         scale = 'none', 
         cluster_rows = F, 
         cluster_cols = T, 
         fontsize_row = 5, fontsize_col = 8,
         fontsize = 8,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         treeheight_row = 30,
         annotation_col = sponge.anno_col,
         annotation_colors = sponge.anno_metabo_colors,
         border_color = 'NA',
         main = 'Sponge data - SVD')
#Linear model per variable
# Sponge data
sponge.form <- ~ sponge.trt + sponge.batch
sponge.info <- as.data.frame(cbind(rownames(sponge.tss.clr), sponge.trt, sponge.batch))
rownames(sponge.info) <- rownames(sponge.tss.clr)

# before
sponge.varPart.before <- fitExtractVarPartModel(exprObj = t(sponge.tss.clr), 
                                                formula = sponge.form, 
                                                data = sponge.info)

# BMC
sponge.varPart.bmc <- fitExtractVarPartModel(exprObj = t(sponge.bmc), 
                                             formula = sponge.form, 
                                             data = sponge.info)

# combat
sponge.varPart.combat <- fitExtractVarPartModel(exprObj = t(sponge.combat), 
                                                formula = sponge.form, 
                                                data = sponge.info)

# removeBatchEffect
sponge.varPart.limma <- fitExtractVarPartModel(exprObj = t(sponge.limma), 
                                               formula = sponge.form, 
                                               data = sponge.info)

# percentile normalisation
sponge.varPart.percentile <- fitExtractVarPartModel(exprObj = t(sponge.percentile), 
                                                    formula = sponge.form, 
                                                    data = sponge.info)

# svd
sponge.varPart.svd <- fitExtractVarPartModel(exprObj = t(sponge.svd), 
                                             formula = sponge.form, 
                                             data = sponge.info)
# extract the variance of trt and batch
# before
sponge.varmat.before <- as.matrix(sponge.varPart.before[ ,1:2])
# BMC
sponge.varmat.bmc <- as.matrix(sponge.varPart.bmc[ ,1:2])
# ComBat
sponge.varmat.combat <- as.matrix(sponge.varPart.combat[ ,1:2])
# removeBatchEffect
sponge.varmat.limma <- as.matrix(sponge.varPart.limma[ ,1:2])
# percentile normalisation
sponge.varmat.percentile <- as.matrix(sponge.varPart.percentile[ ,1:2])
# SVD
sponge.varmat.svd <- as.matrix(sponge.varPart.svd[ ,1:2])

# merge results
sponge.variance <- c(as.vector(sponge.varmat.before), as.vector(sponge.varmat.bmc),
                     as.vector(sponge.varmat.combat), as.vector(sponge.varmat.limma),
                     as.vector(sponge.varmat.percentile), as.vector(sponge.varmat.svd))

# add batch, trt and methods info
sponge.variance <- cbind(variance = sponge.variance, 
                         Type = rep(c('Tissue', 'Batch'), each = ncol(sponge.tss.clr)),
                         method = rep(c('Before', 'BMC', 'ComBat', 'rBE', 
                                        'PN', 'SVD'), each = 2*ncol(sponge.tss.clr)))
# reorder levels  
sponge.variance <- as.data.frame(sponge.variance)
sponge.variance$method <- factor(sponge.variance$method, 
                                 levels = unique(sponge.variance$method))
sponge.variance$variance <- as.numeric(as.character(sponge.variance$variance))

ggplot(sponge.variance, aes(x = Type, y = variance, fill = Type)) + 
  geom_boxplot() + facet_grid(cols = vars(method)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        strip.text = element_text(size = 12), panel.grid = element_blank(), 
        axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) + 
  labs(x = 'Type', y = 'Proportion Variance', name = 'Type') + ylim(0,1)

#RDA

#The pRDA method is a multivariate method to assess globally the effect of batch and treatment (two separate models).
# Sponge data
sponge.data.design <- numeric()
sponge.data.design$group <- sponge.trt
sponge.data.design$batch <- sponge.batch

# before
# conditioning on a batch effect
sponge.rda.before1 <- rda(sponge.tss.clr ~ group + Condition(batch), 
                          data = sponge.data.design)
sponge.rda.before2 <- rda(sponge.tss.clr ~ batch + Condition(group), 
                          data = sponge.data.design)

# amount of variance
sponge.rda.bat_prop.before <- sponge.rda.before1$pCCA$tot.chi*100/sponge.rda.before1$tot.chi
sponge.rda.trt_prop.before <- sponge.rda.before2$pCCA$tot.chi*100/sponge.rda.before2$tot.chi

# BMC
# conditioning on a batch effect
sponge.rda.bmc1 <- rda(sponge.bmc ~ group + Condition(batch), 
                       data = sponge.data.design)
sponge.rda.bmc2 <- rda(sponge.bmc ~ batch + Condition(group), 
                       data = sponge.data.design)

# amount of variance
sponge.rda.bat_prop.bmc <- sponge.rda.bmc1$pCCA$tot.chi*100/sponge.rda.bmc1$tot.chi
sponge.rda.trt_prop.bmc <- sponge.rda.bmc2$pCCA$tot.chi*100/sponge.rda.bmc2$tot.chi


# combat
# conditioning on a batch effect
sponge.rda.combat1 <- rda(sponge.combat ~ group + Condition(batch), 
                          data = sponge.data.design)
sponge.rda.combat2 <- rda(sponge.combat ~ batch + Condition(group), 
                          data = sponge.data.design)

# amount of variance
sponge.rda.bat_prop.combat <- sponge.rda.combat1$pCCA$tot.chi*100/sponge.rda.combat1$tot.chi
sponge.rda.trt_prop.combat <- sponge.rda.combat2$pCCA$tot.chi*100/sponge.rda.combat2$tot.chi


# limma
# conditioning on a batch effect
sponge.rda.limma1 <- rda(sponge.limma ~ group + Condition(batch), 
                         data = sponge.data.design)
sponge.rda.limma2 <- rda(sponge.limma ~ batch + Condition(group), 
                         data = sponge.data.design)

# amount of variance
sponge.rda.bat_prop.limma <- sponge.rda.limma1$pCCA$tot.chi*100/sponge.rda.limma1$tot.chi
sponge.rda.trt_prop.limma <- sponge.rda.limma2$pCCA$tot.chi*100/sponge.rda.limma2$tot.chi


# percentile
# conditioning on a batch effect
sponge.rda.percentile1 <- rda(sponge.percentile ~ group + Condition(batch), 
                              data = sponge.data.design)
sponge.rda.percentile2 <- rda(sponge.percentile ~ batch + Condition(group), 
                              data = sponge.data.design)

# amount of variance
sponge.rda.bat_prop.percentile <- sponge.rda.percentile1$pCCA$tot.chi*100/sponge.rda.percentile1$tot.chi
sponge.rda.trt_prop.percentile <- sponge.rda.percentile2$pCCA$tot.chi*100/sponge.rda.percentile2$tot.chi


# SVD
# conditioning on a batch effect
sponge.rda.svd1 <- rda(sponge.svd ~ group + Condition(batch), 
                       data = sponge.data.design)
sponge.rda.svd2 <- rda(sponge.svd ~ batch + Condition(group), 
                       data = sponge.data.design)

# amount of variance
sponge.rda.bat_prop.svd <- sponge.rda.svd1$pCCA$tot.chi*100/sponge.rda.svd1$tot.chi
sponge.rda.trt_prop.svd <- sponge.rda.svd2$pCCA$tot.chi*100/sponge.rda.svd2$tot.chi

# proportion
sponge.rda.prop.before <- c(sponge.rda.bat_prop.before, 
                            sponge.rda.trt_prop.before)
sponge.rda.prop.bmc <- c(sponge.rda.bat_prop.bmc, 
                         sponge.rda.trt_prop.bmc)
sponge.rda.prop.combat <- c(sponge.rda.bat_prop.combat, 
                            sponge.rda.trt_prop.combat)
sponge.rda.prop.limma <- c(sponge.rda.bat_prop.limma, 
                           sponge.rda.trt_prop.limma)
sponge.rda.prop.percentile <- c(sponge.rda.bat_prop.percentile, 
                                sponge.rda.trt_prop.percentile)
sponge.rda.prop.svd <- c(sponge.rda.bat_prop.svd, 
                         sponge.rda.trt_prop.svd)

# merge results
sponge.rda.prop.val <- c(sponge.rda.prop.before, sponge.rda.prop.bmc, 
                         sponge.rda.prop.combat, sponge.rda.prop.limma, 
                         sponge.rda.prop.percentile, sponge.rda.prop.svd)

# add batch, trt and method info
sponge.rda.prop <- data.frame(prop = sponge.rda.prop.val, 
                              prop.r = round(sponge.rda.prop.val, 2), 
                              Method = rep(c('Before', 'BMC', 'ComBat', 
                                             'rBE', 'PN', 'SVD'), each = 2), 
                              Type = rep(c('Batch', 'Tissue'), 6))

# reorder levels
sponge.rda.prop$Method <- factor(sponge.rda.prop$Method, 
                                 levels = unique(sponge.rda.prop$Method))

ggplot(data = sponge.rda.prop, aes(x = Method, y = prop, fill = Type)) + 
  geom_bar(stat = "identity", position = 'dodge', colour = 'black') + 
  geom_text(data = sponge.rda.prop, aes(Method, prop + 2.5, label = prop.r), 
            position = position_dodge(width = 0.9), size = 3) + theme_bw() + 
  labs(y = "Variance explained (%)") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.grid = element_blank(), axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), legend.title = element_text(size = 15),
        legend.text = element_text(size = 12)) + ylim(0,100)


#Silhouette coefficient
# Sponge data
sponge.silh.before <- calc.sil(sponge.pca.before$variates$X, 
                               y1 = sponge.batch, y2 = sponge.trt, 
                               name.y1 = 'Batch', name.y2 = 'Tissue')
sponge.silh.bmc <- calc.sil(sponge.pca.bmc$variates$X, 
                            y1 = sponge.batch, y2 = sponge.trt, 
                            name.y1 = 'Batch', name.y2 = 'Tissue')
sponge.silh.combat <- calc.sil(sponge.pca.combat$variates$X, 
                               y1 = sponge.batch, y2 = sponge.trt, 
                               name.y1 = 'Batch', name.y2 = 'Tissue')
sponge.silh.limma <- calc.sil(sponge.pca.limma$variates$X, 
                              y1 = sponge.batch, y2 = sponge.trt, 
                              name.y1 = 'Batch', name.y2 = 'Tissue')
sponge.silh.percentile <- calc.sil(sponge.pca.percentile$variates$X, 
                                   y1 = sponge.batch, y2 = sponge.trt, 
                                   name.y1 = 'Batch', name.y2 = 'Tissue')
sponge.silh.svd <- calc.sil(sponge.pca.svd$variates$X, 
                            y1 = sponge.batch, y2 = sponge.trt, 
                            name.y1 = 'Batch', name.y2 = 'Tissue')


sponge.silh.plot <- rbind(sponge.silh.before, sponge.silh.bmc, sponge.silh.combat, 
                          sponge.silh.limma, sponge.silh.percentile, sponge.silh.svd)
sponge.silh.plot$method <- c(rep('Before', nrow(sponge.silh.before)), 
                             rep('BMC', nrow(sponge.silh.bmc)),
                             rep('ComBat', nrow(sponge.silh.combat)),
                             rep('rBE', nrow(sponge.silh.limma)),
                             rep('PN', nrow(sponge.silh.percentile)),
                             rep('SVD', nrow(sponge.silh.svd))
)
sponge.silh.plot$method <- factor(sponge.silh.plot$method, 
                                  levels = unique(sponge.silh.plot$method))
sponge.silh.plot$Cluster <- factor(sponge.silh.plot$Cluster, 
                                   levels = unique(sponge.silh.plot$Cluster))
sponge.silh.plot$Type <- factor(sponge.silh.plot$Type, 
                                levels = unique(sponge.silh.plot$Type))


ggplot(sponge.silh.plot, aes(x = Type, y = silh.coeff, color = Cluster, shape = Type)) + 
  geom_point() + facet_grid(cols = vars(method)) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        strip.text = element_text(size = 12), panel.grid = element_blank(), 
        axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) + 
  scale_color_manual(values = c('#388ECC','#F68B33','#F0E442','#D55E00')) + 
  labs(x = 'Type', y = 'Silhouette Coefficient', name = 'Type') 

