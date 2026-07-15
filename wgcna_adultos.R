library(WGCNA)
library(flashClust)
library(ggdendro)
library(ggplot2)
library(patchwork)

df_exp <- read.csv(file = '..//TC_adultos_promediadoXtiempoXadjusted_tpm_estimate.csv', stringsAsFactors = FALSE, header = TRUE)
head(df_exp)

rownames(df_exp) <- df_exp$cell_type
df_exp <- df_exp[, !(names(df_exp) == "cell_type")]


# ---- 1) Data cleaning ----
#identificación de outliers:
gsg <-goodSamplesGenes(df_exp)
summary(gsg)

# (...) By viewing the gsg list object you can see it contains 3 logical vectors 
# (good genes, good samples and allOK). If you want to see if the function 
# identified any possible outlier all you have to do is evaluate the allOK vector.

gsg$allOK
#Como da false hacemos lo siguiente:
# Para sacar GENES
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(df_exp)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(df_exp)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  expression.data <- df_exp[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

# Para sacar SAMPLES
sampleTree <- hclust(dist(expression.data), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
abline(h = 225000, col = "red");

cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 225000, minSize = 10) #returns numeric vector

# Las hojas que están por encima de la linea roja pueden ser potencialmente 
# outliers, por lo que se eliminan.

#Remove outlier
expression.data <- expression.data[cut.sampleTree==1, ]

# ---- 2) Network construction ----
# ELEGIR SOFTPOWER: demora unos minutos
spt <- pickSoftThreshold(expression.data) 

#Scale free topology 
# Acá tenemos que ver con qué softpower tenemos todos los valores por encima del cutoff
# En este caso el cutoff es 0.8 que es lo que recomiendan ellos, pero puede ser flexible. 
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.80,col="red")

# Mean conectivity
par(mar=c(1,1,1,1))
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")

# MATRIZ DE ADYACENCIA
softPower <- 2 #esto lo definimos por la gráfica de arriba.
adjacency <- adjacency(expression.data, power = softPower) #demora unos minutos

# ---- 3) Module construction ----
TOM <- TOMsimilarity(adjacency) #demora cerca de una hora
TOM.dissimilarity <- 1-TOM

#creating the dendrogram 
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")
#plotting the dendrogram
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)
# Modulos
Modules <- cutreeDynamic(dendro = geneTree, 
                         distM = TOM.dissimilarity, 
                         deepSplit = 2, 
                         pamRespectsDendro = FALSE, 
                         minClusterSize = 50)
ModuleColors <- labels2colors(Modules) #assigns each module number a color
head(ModuleColors)
table(ModuleColors) #returns the counts for each color (aka the number of genes within each module)

length(unique(ModuleColors))


#plots the gene dendrogram with the module colors
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Module Eigengene Identification

MElist <- moduleEigengenes(expression.data, colors = ModuleColors) 
ncol(MElist$eigengenes)
MEs <- MElist$eigengenes 
head(MEs)

# Aca vemos si nos quedó algun módulo con nas
length(unique(ModuleColors))
anyNA(MElist$eigengenes)
colSums(is.na(MElist$eigengenes))

# Nos quedó el módulo gray, asi que lo sacamos
length(unique(ModuleColors))
MElist$eigengenes <- MElist$eigengenes[, colSums(is.na(MElist$eigengenes)) == 0]
ncol(MElist$eigengenes)

# ahora si calculamos la mat.dis
ME.dissimilarity <- 1 - cor(MElist$eigengenes, use = "pairwise.complete.obs")

METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 

#Para graficar:
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(METree)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75

merge <- mergeCloseModules(expression.data, ModuleColors, cutHeight = .25)

# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs
length(unique(mergedMEs))
# Para ponerle nombre a los módulos
if (is.null(names(mergedColors))) {
names(mergedColors) <- colnames(expression.data)
}

plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")

# 1. Eliminar MEgrey del objeto mergedMEs
mergedMEs_clean <- mergedMEs[, !grepl("^MEgrey$", colnames(mergedMEs))]

# 2. Calcular matriz de disimilitud entre los eigengenes fusionados
dissME_merged <- 1 - cor(mergedMEs_clean, use = "pairwise.complete.obs")

# 3. Clustering jerárquico
treeME_merged <- hclust(as.dist(dissME_merged), method = "average")

# 2. Extraer nombres de los módulos fusionados
# Obtener los nombres reales de los módulos sin el prefijo "ME"
moduleNames <- substring(colnames(mergedMEs_clean), 3)

# Asignar colores directamente (esto es solo simbólico, da colores únicos)
mergedModuleColors <- labels2colors(factor(moduleNames))
# Graficar el dendrograma con los colores ordenados
plotDendroAndColors(
  dendro = treeME_merged,
  colors = mergedModuleColors[treeME_merged$order],
  groupLabels = "Merged module colors",
  hang = -1,
  addGuide = TRUE,
  guideHang = 0.05,
  dendroLabels = FALSE
)

par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
plot(treeME_merged)
abline(h=.25, col = "red")

# -----------------------------------
#           GUARDAMOS TODO 
# -----------------------------------
table(ModuleColors)
save(adjacency, TOM, mergedColors, file = "network_objects_adultos.RData")
save(MElist, file = "adul_MElist.RData")
save(treeME_merged,dissME_merged, file = "adul_dendro.RData")








# ---- 4) Ver cómo quedan los genes de la muestra de entrenamiento ----

load("network_objects_adultos.RData")

# 1. Cargar lista de genes desde el archivo CSV
genes <- read.csv("C:\\Users\\masoz\\Desktop\\TesisMaestria\\Papers\\Paper1\\entrenamiento_sin_paralogos\\genes_actualizado.csv", stringsAsFactors = FALSE)

# 2. Filtrar genes que se usan para entrenar
filas_seleccionadas <- genes[genes$Usar.para.entrenar == "si", ]

# 3. Obtener la lista de WBIDs
muestra_positiva <- filas_seleccionadas$WBID

# 4. Filtrar los que están en la red (es decir, en ModuleColors)
names(ModuleColors) <- colnames(expression.data)
genes_en_red <- intersect(muestra_positiva, names(ModuleColors))

# 5. Obtener módulo correspondiente
modulos_genes_positivos <- ModuleColors[genes_en_red]

# 6. Ver cantidad de genes por módulo
table(modulos_genes_positivos)

# 7. Generar un data.frame con el detalle
df_modulos <- data.frame(
  Gen = genes_en_red,
  Modulo = modulos_genes_positivos
)

# 8. Crear conteo total de genes por módulo
conteo_modulos <- table(ModuleColors)
df_conteo <- as.data.frame(conteo_modulos)
colnames(df_conteo) <- c("Modulo", "Cantidad")

# 9. Contar cuántos genes de interés hay por módulo
conteo_positivos <- table(modulos_genes_positivos)
df_positivos <- as.data.frame(conteo_positivos)
colnames(df_positivos) <- c("Modulo", "Positivos")

# 10. Eliminar columna duplicada si ya existe
df_conteo$Positivos <- NULL

# 11. Unir conteos
df_conteo <- merge(df_conteo, df_positivos, by = "Modulo", all.x = TRUE)
df_conteo$Positivos[is.na(df_conteo$Positivos)] <- 0

# 12. Obtener orden de módulos del dendrograma
MEs_clean <- MElist$eigengenes[, !grepl("^MEgrey$", colnames(MElist$eigengenes))]
dissME <- 1 - cor(MEs_clean, use = "pairwise.complete.obs")
treeME <- hclust(as.dist(dissME), method = "average")
orden_dendrograma <- gsub("^ME", "", colnames(MEs_clean)[treeME$order])

# 13. Reordenar df_conteo según dendrograma
df_conteo$Modulo <- factor(df_conteo$Modulo, levels = orden_dendrograma)
df_conteo <- df_conteo[order(df_conteo$Modulo), ]

# 14. Gráfico de barras
p_barras <- ggplot(df_conteo, aes(x = Modulo, y = Cantidad)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  geom_text(aes(label = Positivos), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(x = NULL, y = "Cantidad de genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

# 15. Dendrograma compatible con patchwork
dendro_data <- ggdendrogram(treeME, rotate = FALSE, theme_dendro = FALSE)

# 16. Combinar ambos gráficos
dendro_data / p_barras + plot_layout(heights = c(1, 2))



# ---- 5) Predictions ----
# 1. Leer lista de genes predichos desde archivo de texto
predicciones <- scan("C:\\Users\\masoz\\Desktop\\TesisMaestria\\Papers\\Paper1\\entrenamiento_con_paralogos\\intersected_genes_randomx10_w_patalogs.txt",
                     what = character(), sep = "\n", quiet = TRUE)

#  2. Filtrar los genes que están en la red 
names(ModuleColors) <- colnames(expression.data)
genes_en_red <- intersect(predicciones, names(ModuleColors))

# 3. Obtener el módulo correspondiente para cada gen positivo
modulos_genes_positivos <- ModuleColors[genes_en_red]

# 4. Generar data frame con detalle (opcional) 
df_modulos <- data.frame(
  Gen = genes_en_red,
  Modulo = modulos_genes_positivos
)

# 5. Conteo total de genes por módulo
conteo_modulos <- table(ModuleColors)
df_conteo <- as.data.frame(conteo_modulos)
colnames(df_conteo) <- c("Modulo", "Cantidad")

# 6. Conteo de predicciones por módulo
conteo_positivos <- table(modulos_genes_positivos)
df_positivos <- as.data.frame(conteo_positivos)
colnames(df_positivos) <- c("Modulo", "Positivos")

# 7. Limpiar y unir conteos
df_conteo$Positivos <- NULL
df_conteo <- merge(df_conteo, df_positivos, by = "Modulo", all.x = TRUE)
df_conteo$Positivos[is.na(df_conteo$Positivos)] <- 0

# 8. Dendrograma de eigengenes y orden 
MEs_clean <- MElist$eigengenes[, !grepl("^MEgrey$", colnames(MElist$eigengenes))]
dissME <- 1 - cor(MEs_clean, use = "pairwise.complete.obs")
treeME <- hclust(as.dist(dissME), method = "average")
orden_dendrograma <- gsub("^ME", "", colnames(MEs_clean)[treeME$order])

# 9. Reordenar módulos para la gráfica 
df_conteo$Modulo <- factor(df_conteo$Modulo, levels = orden_dendrograma)
df_conteo <- df_conteo[order(df_conteo$Modulo), ]

# 10. Gráfico de barras 
p_barras <- ggplot(df_conteo, aes(x = Modulo, y = Cantidad)) +
  geom_bar(stat = "identity", fill = "grey", color = "black") +
  geom_text(aes(label = Positivos), vjust = -0.5, size = 3.5) +
  theme_minimal() +
  labs(x = NULL, y = "Cantidad de genes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")

# 11. Dendrograma 
dendro_data <- ggdendrogram(treeME, rotate = FALSE, theme_dendro = FALSE)

# 12. Combinar ambos gráficos
dendro_data / p_barras + plot_layout(heights = c(1, 2))





library(ggplot2)
library(ggdendro)
library(patchwork)

# ------------------ INPUT ------------------
# Cargar lista de genes predichos
predicciones <- scan("C:\\Users\\masoz\\Desktop\\TesisMaestria\\Papers\\Paper1\\entrenamiento_con_paralogos\\intersected_genes_randomx10_w_patalogs.txt",
                     what = character(), sep = "\n", quiet = TRUE)

# Cargar lista de genes positivos (muestra de entrenamiento)
genes <- read.csv("C:\\Users\\masoz\\Desktop\\TesisMaestria\\Papers\\Paper1\\entrenamiento_sin_paralogos\\genes_actualizado.csv",
                  stringsAsFactors = FALSE)
positivos <- genes$WBID[genes$Usar.para.entrenar == "si"]

# Asignar nombres
names(ModuleColors) <- colnames(expression.data)
total_genes <- names(ModuleColors)

# -------------------- 2. Conteos --------------------
df <- data.frame(Modulo = names(table(ModuleColors)))
df$Cantidad <- as.integer(table(ModuleColors)[df$Modulo])

# Positivos
positivos_en_red <- intersect(positivos, total_genes)
mod_positivos <- ModuleColors[positivos_en_red]
df$Positivos <- as.integer(table(factor(mod_positivos, levels = df$Modulo)))

# Predichos
pred_en_red <- intersect(predicciones, total_genes)
mod_pred <- ModuleColors[pred_en_red]
df$Predichos <- as.integer(table(factor(mod_pred, levels = df$Modulo)))

# -------------------- 3. Enriquecimiento y e-value --------------------
total_pos <- length(positivos_en_red)
total_pred <- length(pred_en_red)
total_all <- length(total_genes)

# Fisher + fold enrichment
df$Pvalor <- NA
df$Evalue <- NA
df$FoldEnrich_Pred <- NA
df$FoldEnrich_Pos <- NA

for (i in seq_len(nrow(df))) {
  # Valores para predichos
  a <- df$Predichos[i]
  b <- df$Cantidad[i] - a
  c <- total_pred - a
  d <- (total_all - df$Cantidad[i]) - c
  m <- matrix(c(a, b, c, d), nrow = 2)
  p <- fisher.test(m, alternative = "greater")$p.value
  df$Pvalor[i] <- p
  df$FoldEnrich_Pred[i] <- (a / df$Cantidad[i]) / (total_pred / total_all)
  
  # Valores para positivos
  a2 <- df$Positivos[i]
  b2 <- df$Cantidad[i] - a2
  c2 <- total_pos - a2
  d2 <- (total_all - df$Cantidad[i]) - c2
  df$FoldEnrich_Pos[i] <- (a2 / df$Cantidad[i]) / (total_pos / total_all)
}

df$Evalue <- p.adjust(df$Pvalor, method = "BH")
df$FoldEnrich_Pred[is.nan(df$FoldEnrich_Pred)] <- 0
df$FoldEnrich_Pos[is.nan(df$FoldEnrich_Pos)] <- 0

# -------------------- 4. Significancia y etiquetas --------------------
df$ColorPred <- ifelse(df$Predichos > 0, df$FoldEnrich_Pred, NA)
df$ColorPos <- ifelse(df$Positivos > 0, df$FoldEnrich_Pos, NA)

df$SignifPred <- ifelse(df$Evalue < 0.05 & df$Predichos > 0, "*", "")
df$SignifPos <- ifelse(df$Evalue < 0.05 & df$Positivos > 0, "*", "")

# -------------------- 5. Orden por dendrograma --------------------
MEs_clean <- MElist$eigengenes[, !grepl("^MEgrey$", colnames(MElist$eigengenes))]
dissME <- 1 - cor(MEs_clean, use = "pairwise.complete.obs")
treeME <- hclust(as.dist(dissME), method = "average")
orden_dendro <- gsub("^ME", "", colnames(MEs_clean)[treeME$order])

df$Modulo <- factor(df$Modulo, levels = orden_dendro)
df <- df[order(df$Modulo), ]

# -------------------- 6. Gráficos --------------------

# POSITIVOS
p_positivos <- ggplot(df, aes(x = Modulo, y = Cantidad)) +
  geom_bar(stat = "identity", aes(fill = ColorPos)) +
  geom_text(aes(label = ifelse(Positivos > 0, Positivos, "")), vjust = -0.5, size = 3.5) +
  geom_text(aes(label = SignifPos), vjust = -2, size = 6) +
  scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = "grey") +
  labs(title = "Genes de entrenamiento", y = "Genes totales", fill = "Fold enrichment") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

# PREDICHOS
p_predichos <- ggplot(df, aes(x = Modulo, y = Cantidad)) +
  geom_bar(stat = "identity", aes(fill = ColorPred)) +
  geom_text(aes(label = ifelse(Predichos > 0, Predichos, "")), vjust = -0.5, size = 3.5) +
  geom_text(aes(label = SignifPred), vjust = -2, size = 6) +
  scale_fill_gradient(low = "lightsalmon", high = "firebrick", na.value = "grey") +
  labs(title = "Genes predichos", y = "Genes totales", fill = "Fold enrichment") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

# DENDROGRAMA LIMPIO CON EJE X (módulos)
ddata <- dendro_data(treeME, type = "rectangle")

# Extraer etiquetas
labels_ddata <- label(ddata)
labels_ddata$label <- gsub("^ME", "", labels_ddata$label)

p_dendro <- ggplot() +
  geom_segment(data = segment(ddata),
               aes(x = x, y = y, xend = xend, yend = yend)) +
  scale_x_continuous(breaks = labels_ddata$x,
                     labels = labels_ddata$label) +
  labs(x = "Módulo", y = "Dissimilarity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )


# -------------------- 7. Layout final --------------------
p_dendro / p_positivos / p_predichos + plot_layout(heights = c(0.1, 1.8, 1.8))

ggsave("Enriquecimiento_positivos_adultos.png", plot = p_positivos, width = 6, height = 4, dpi = 300)
ggsave("Enriquecimiento_predichos_adultos.png", plot = p_predichos, width = 6, height = 4, dpi = 300)
ggsave("dendograma_grafo_adultos.png", plot = p_dendro, width = 6, height = 4, dpi = 300)



# ------------- otra visualización: 

library(ggplot2)
library(ggdendro)
library(patchwork)

# -------------------- 1. Dendrograma vertical --------------------
ddata <- dendro_data(treeME, type = "rectangle")

# Extraer etiquetas (módulos) y su orden
labels_ddata <- label(ddata)
labels_ddata$label <- gsub("^ME", "", labels_ddata$label)
modulos_ordenados <- labels_ddata$label
df$Modulo <- factor(df$Modulo, levels = modulos_ordenados)
df <- df[order(df$Modulo), ]

# Dendrograma rotado (izquierda)
p_dendro_vert <- ggplot() +
  geom_segment(data = segment(ddata),
               aes(x = y, y = x, xend = yend, yend = xend)) +
  scale_y_continuous(breaks = labels_ddata$x,
                     labels = labels_ddata$label,
                     expand = expansion(add = 0.5)) +
  labs(x = "Disimilitud", y = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        panel.background = element_blank())

# -------------------- 2. Preparar datos --------------------

# Tamaños (número de genes)
df$SizePos <- ifelse(df$Positivos > 0, df$Positivos, NA)
df$SizePred <- ifelse(df$Predichos > 0, df$Predichos, NA)

# Color por enrichment
df$ColorPos <- ifelse(df$Positivos > 0, df$FoldEnrich_Pos, NA)
df$ColorPred <- ifelse(df$Predichos > 0, df$FoldEnrich_Pred, NA)

# Asteriscos
df$AstPos <- ifelse(df$Evalue < 0.05 & df$Positivos > 0, "*", "")
df$AstPred <- ifelse(df$Evalue < 0.05 & df$Predichos > 0, "*", "")

# -------------------- 3. Gráfica de puntos (entrenamiento) --------------------
p_pos_dots <- ggplot(df, aes(y = Modulo, x = 1)) +
  geom_point(aes(size = SizePos, color = ColorPos)) +
  geom_text(aes(label = AstPos), hjust = -1, size = 6) +
  scale_color_gradient(low = "lightblue", high = "darkblue", na.value = "grey") +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Genes de entrenamiento", x = NULL, y = NULL, color = "Fold enrichment", size = "Cantidad") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")

# -------------------- 4. Gráfica de puntos (predicciones) --------------------
p_pred_dots <- ggplot(df, aes(y = Modulo, x = 1)) +
  geom_point(aes(size = SizePred, color = ColorPred)) +
  geom_text(aes(label = AstPred), hjust = -1, size = 6) +
  scale_color_gradient(low = "lightsalmon", high = "firebrick", na.value = "grey") +
  scale_size_continuous(range = c(2, 10)) +
  labs(title = "Genes predichos", x = NULL, y = NULL, color = "Fold enrichment", size = "Cantidad") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right")

# -------------------- 5. Layout horizontal --------------------
(p_dendro_vert | p_pos_dots | p_pred_dots) + plot_layout(widths = c(1.5, 1, 1))


# ------------- circular:
library(circlize)
library(stats)  # asegúrate de tener treeME como objeto hclust

dend <- as.dendrogram(treeME)
labels(dend) <- gsub("^ME", "", labels(dend))  # remover "ME"
orden_dendro <- labels(dend)  # nuevo orden
# Filtrar y ordenar
df <- df[df$Modulo %in% orden_dendro, ]
df <- df[match(orden_dendro, df$Modulo), ]
modules <- df$Modulo

all(modules == labels(dend))  # debería ser TRUE









# Quitar el prefijo "ME" de cada nombre
orden_dendro <- labels(as.dendrogram(treeME))
orden_dendro <- gsub("^ME", "", orden_dendro)

# Asegurar que df$Modulo sea carácter
df$Modulo <- as.character(df$Modulo)

# Filtrar y reordenar
orden_dendro_filtrado <- orden_dendro[orden_dendro %in% df$Modulo]
df <- df[match(orden_dendro_filtrado, df$Modulo), ]
modules <- df$Modulo




# 2. Inicializar circos
circos.clear()
circos.par(gap.degree = 1, cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = modules,
                  xlim = cbind(rep(0, length(modules)), rep(1, length(modules))))

# 3. Track externo: nombres de módulos
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(0.5, 0.5, CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.7)
}, track.height = 0.06, bg.border = NA)

# 4. Track puntos predichos
max_pred <- max(df$Predichos, na.rm = TRUE)
col_fun_pred <- colorRamp2(c(0, max(df$FoldEnrich_Pred, na.rm = TRUE)),
                           c("lightcoral", "firebrick"))
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  idx <- CELL_META$sector.index
  i <- which(df$Modulo == idx)
  size <- df$Predichos[i]; color <- df$FoldEnrich_Pred[i]
  signif <- df$Evalue[i] < 0.05 & size > 0
  if (size > 0) {
    circos.points(0.5, 0.5, pch = 16,
                  col = col_fun_pred(color),
                  cex = size / max_pred * 2)
    if (signif) circos.text(0.5, 0.8, "*", cex = 1.2,
                            col = "black", facing = "inside")
  }
}, track.height = 0.12, bg.border = NA)

# 5. Track puntos entrenamiento
max_pos <- max(df$Positivos, na.rm = TRUE)
col_fun_pos <- colorRamp2(c(0, max(df$FoldEnrich_Pos, na.rm = TRUE)),
                          c("lightblue", "darkblue"))
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  idx <- CELL_META$sector.index
  i <- which(df$Modulo == idx)
  size <- df$Positivos[i]
  color <- df$FoldEnrich_Pos[i]
  signif <- df$Evalue[i] < 0.05 & size > 0
  
  if (size > 0) {
    circos.points(0.5, 0.5, pch = 16,
                  col = col_fun_pos(color),
                  cex = size / max_pos * 2)
    if (signif) circos.text(0.5, 0.8, "*", cex = 1.2,
                            col = "black", facing = "inside")
  }
}, track.height = 0.12, bg.border = NA)


# 6. Track barras de cantidad total
max_count <- max(df$Cantidad)
circos.trackPlotRegion(ylim = c(0, max_count), panel.fun = function(x, y) {
  idx <- CELL_META$sector.index
  i <- which(df$Modulo == idx)
  val <- df$Cantidad[i]
  circos.rect(0.2, 0, 0.8, val, col = "grey50", border = NA)
}, track.height = 0.1, bg.border = NA)

# 7. Track dendrograma real (interior)
max_height <- attr(dend, "height")
circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,
                       track.height = 0.2, panel.fun = function(x, y) {
                         circos.dendrogram(dend, max_height = max_height)
                       })





library(circlize)
library(dendextend)

# 1. Preparar módulos y df
modules <- as.character(df$Modulo)
df <- df[match(modules, df$Modulo), ]

# 2. Arreglar nombres del dendrograma
dend <- as.dendrogram(treeME)
labels(dend) <- gsub("^ME", "", labels(dend))  # eliminar "ME" para que coincida con df

# 3. Reordenar módulos y df según el orden del dendrograma
orden_dendro <- labels(dend)
modules <- orden_dendro
df <- df[match(modules, df$Modulo), ]

# 4. Inicializar circos
circos.clear()
circos.par(gap.degree = 1, cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = modules,
                  xlim = cbind(rep(0, length(modules)), rep(1, length(modules))))

# 5. Track externo: nombres de módulos
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  circos.text(0.5, 0.5, CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.7)
}, track.height = 0.06, bg.border = NA)

# 6. Track puntos predichos
max_pred <- max(df$Predichos, na.rm = TRUE)
col_fun_pred <- colorRamp2(c(0, max(df$FoldEnrich_Pred, na.rm = TRUE)),
                           c("lightcoral", "firebrick"))
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  idx <- CELL_META$sector.index
  i <- which(df$Modulo == idx)
  size <- df$Predichos[i]; color <- df$FoldEnrich_Pred[i]
  signif <- df$Evalue[i] < 0.05 & size > 0
  if (size > 0) {
    circos.points(0.5, 0.5, pch = 16,
                  col = col_fun_pred(color),
                  cex = size / max_pred * 2)
    if (signif) circos.text(0.5, 0.8, "*", cex = 1.2,
                            col = "black", facing = "inside")
  }
}, track.height = 0.12, bg.border = NA)

# 7. Track puntos positivos
max_pos <- max(df$Positivos, na.rm = TRUE)
col_fun_pos <- colorRamp2(c(0, max(df$FoldEnrich_Pos, na.rm = TRUE)),
                          c("lightblue", "darkblue"))
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  idx <- CELL_META$sector.index
  i <- which(df$Modulo == idx)
  size <- df$Positivos[i]
  color <- df$FoldEnrich_Pos[i]
  signif <- df$Evalue[i] < 0.05 & size > 0
  if (size > 0) {
    circos.points(0.5, 0.5, pch = 16,
                  col = col_fun_pos(color),
                  cex = size / max_pos * 2)
    if (signif) circos.text(0.5, 0.8, "*", cex = 1.2,
                            col = "black", facing = "inside")
  }
}, track.height = 0.12, bg.border = NA)

# 8. Track barras cantidad total
max_count <- max(df$Cantidad)
circos.trackPlotRegion(ylim = c(0, max_count), panel.fun = function(x, y) {
  idx <- CELL_META$sector.index
  i <- which(df$Modulo == idx)
  val <- df$Cantidad[i]
  circos.rect(0.2, 0, 0.8, val, col = "grey50", border = NA)
}, track.height = 0.1, bg.border = NA)

# 9. Track dendrograma (interior)
circos.trackPlotRegion(ylim = c(0, attr(dend, "height")),
                       panel.fun = function(x, y) {
                         circos.dendrogram(dend, max_height = attr(dend, "height"))
                       },
                       track.height = 0.2, bg.border = NA)

