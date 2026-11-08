##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################

# Abderrahman Bourahi Taleb 
# NP: 154380


# Comprobamos en qué carpeta está trabajando R 
getwd()



### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)

# Columnas separadas por tabuladores y primera fila = nombres de columnas
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "	")

head(datos)   # Primeras 6 filas 
tail(datos)   # Últimas 6 filas 
dim(datos)    # Dimensiones: 619 genes (filas) x 7 columnas 



##################################################################################################################################


### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)

# Guardamos las posiciones de las columnas de cada grupo para no repetir números por todo el script
# (la 7ª columna, 1W5-COV, es la tercera réplica de la condición)
cols_control   <- 2:4   # C2, C1, C3
cols_condicion <- 5:7   # 1W5, 1W6, 1W5-COV

datos$media_control   <- rowMeans(datos[, cols_control])         # Media por gen, control
datos$de_control      <- apply(datos[, cols_control], 1, sd)     # Desviación estándar por gen, control
datos$media_condicion <- rowMeans(datos[, cols_condicion])       # Media por gen, condición
datos$de_condicion    <- apply(datos[, cols_condicion], 1, sd)   # Desviación estándar por gen, condición

head(datos)   # Comprobamos que se han añadido las columnas
dim(datos)    # Ahora 619 filas x 11 columnas (7 originales + 4 nuevas)


##################################################################################################################################


### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).

# Definimos los dos colores UNA sola vez; se reutilizan en el boxplot y en el heatmap
color_control   <- "steelblue"
color_condicion <- "darkorange"

par(mfrow = c(1, 2))   # Ventana gráfica con 1 fila y 2 paneles

hist(datos$media_control, breaks = 30,
     main = "Control", col = color_control,
     xlab = "Expresión media (TPM)", ylab = "Número de genes")

hist(datos$media_condicion, breaks = 30,
     main = "Condición", col = color_condicion,
     xlab = "Expresión media (TPM)", ylab = "Número de genes")

par(mfrow = c(1, 1))   # Restauramos la ventana a un solo panel

# Con los datos originales casi todos los genes se agrupan en la primera barra (hay unos pocos genes
# con expresiones enormes, como ssrA) y no se aprecia la forma de la distribución.
# Transformamos con log2. Sumamos 1 porque hay medias = 0 y log2(0) = -Infinito.
datos$log_control   <- log2(datos$media_control + 1)
datos$log_condicion <- log2(datos$media_condicion + 1)

par(mfrow = c(1, 2))

hist(datos$log_control, breaks = 30,
     main = "Control (escala log2)", col = color_control,
     xlab = "log2(Media + 1)", ylab = "Número de genes")

hist(datos$log_condicion, breaks = 30,
     main = "Condición (escala log2)", col = color_condicion,
     xlab = "log2(Media + 1)", ylab = "Número de genes")

par(mfrow = c(1, 1))


##################################################################################################################################


### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)

# Usamos las medias en log2 (columnas creadas en el apartado 3) para que las cajas se vean bien
boxplot(datos$log_control, datos$log_condicion,
        names = c("Control", "Condición"),
        col = c(color_control, color_condicion),
        main = "Expresión media por gen: Control vs Condición",
        ylab = "log2(Media + 1)",
        ylim = c(0, 20))   # Dejamos hueco arriba para que la leyenda no tape las cajas

legend("topleft", legend = c("Control", "Condición"),
       fill = c(color_control, color_condicion))


##################################################################################################################################


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.

# Test de Shapiro-Wilk (H0: los datos siguen una distribución normal)
shapiro.test(datos$media_control)
shapiro.test(datos$media_condicion)

# En ambos casos p < 2.2e-16, así que rechazamos la normalidad: no se puede usar un t-test.
# Usamos un test no paramétrico: Wilcoxon Mann-Whitney, con muestras independientes (paired = FALSE).

datos$p_valor <- apply(datos[, 2:7], 1, function(fila) {
  wilcox.test(fila[1:3], fila[4:6], paired = FALSE)$p.value   # 1:3 = control, 4:6 = condición
})

head(datos)   # Comprobamos la nueva columna
summary(datos$p_valor)


##################################################################################################################################


### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)

genes_005 <- datos[datos$p_valor < 0.05, ]   # Genes con p-valor < 0,05
nrow(genes_005)                              # Nº de genes significativos
min(datos$p_valor)                           # El p-valor más bajo de toda la tabla

# RESULTADO: ningún gen tiene p < 0,05 (el mínimo es 0,0636).
#
# Justificación:
# Con sólo 3 réplicas por grupo, el test de Wilcoxon tiene muy pocas ordenaciones posibles
# (20 al ordenar las 6 muestras juntas), así que el p-valor más bajo que puede dar es 0,1.
# Por eso es imposible llegar a p < 0,05, aunque la diferencia entre grupos sea grande
# (en nuestros datos el mínimo es 0,0636, por los valores repetidos).
# Con más réplicas sí saldrían genes significativos. Como no es posible,
# usamos un umbral menos estricto (p <= 0,1) para tener genes con los que seguir trabajando.

genes_dif <- datos[datos$p_valor <= 0.1, ]   # Genes diferencialmente expresados con el umbral relajado
nrow(genes_dif)                              # Nº de genes con p <= 0,1

# Para saber el sentido del cambio calculamos el log2 Fold Change (condición frente a control).
datos$log2FC <- datos$log_condicion - datos$log_control    # > 0: más expresado en condición / < 0: menos

# Clasificamos cada gen según su p-valor y el signo del cambio:
datos$estado <- "no significativo"
datos$estado[datos$p_valor <= 0.1 & datos$log2FC > 0] <- "sobreexpresado"
datos$estado[datos$p_valor <= 0.1 & datos$log2FC < 0] <- "reprimido"

table(datos$estado)   # Nº de genes en cada categoría

# RESULTADO: de los 468 genes con p <= 0,1, 69 están sobreexpresados y 399 reprimidos
# (los 151 restantes no son significativos).


##################################################################################################################################


### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)

# Matriz numérica sólo con los genes diferencialmente expresados y las 6 réplicas,
# con el nombre del gen como nombre de fila
matriz <- as.matrix(genes_dif[, 2:7])
rownames(matriz) <- genes_dif$Gene.Name
colnames(matriz) <- c("Ctrl_C2", "Ctrl_C1", "Ctrl_C3", "Cond_1W5", "Cond_1W6", "Cond_1W5.COV")

# Dendrogramas: agrupamos genes y réplicas usando 1 - correlación como distancia
# (cor() trabaja por columnas, por eso para los genes hay que trasponer la matriz con t())
dend_genes <- as.dendrogram(hclust(as.dist(1 - cor(t(matriz)))))
dend_muestras <- as.dendrogram(hclust(as.dist(1 - cor(matriz))))

# Barra de colores sobre las columnas para distinguir grupo (mismos colores que en los apartados 3 y 4)
colores_muestras <- c(rep(color_control, 3), rep(color_condicion, 3))

# Paleta: azul = menos expresado que la media del gen, naranja = más expresado
paleta <- colorRampPalette(c(color_control, "white", color_condicion))(100)

# Heatmap básico.
heatmap(matriz, Rowv = dend_genes, Colv = dend_muestras,
        col = paleta, ColSideColors = colores_muestras,
        labRow = FALSE,        # Con 468 genes los nombres no caben, los ocultamos
        cexCol = 0.8, margins = c(8, 2))

# Heatmap con leyenda de color. Instalamos el paquete sólo si todavía no lo tenemos:
if (!requireNamespace("gplots", quietly = TRUE)) install.packages("gplots")
library(gplots)

heatmap.2(matriz, Rowv = dend_genes, Colv = dend_muestras,
          col = paleta, ColSideColors = colores_muestras,
          scale = "row",                           # Z-score por gen: (valor - media del gen) / DE del gen
          trace = "none", density.info = "none",   # Quitamos las líneas cian internas y el histograma de la leyenda
          key.title = NA, key.xlab = "Z-score",    # Título de la leyenda de color
          labRow = FALSE, cexCol = 0.8, margins = c(8, 2))
