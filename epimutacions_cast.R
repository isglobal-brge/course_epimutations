## Cargamos la librería
library(epimutacions)

## Cargamos un paquete que nos permite acceder a los datasets y cargamos
## los datos que usaremos para el diseño caso-control
library(ExperimentHub)
eh <- ExperimentHub()
query(eh, c("epimutacionsData"))
methy <- eh[["EH6690"]]

## Exploramos la estructura del GenomicRatioSet
methy
assays(methy)$Beta[1:5,1:5]  ## Betas (nivel de metilación)
pData(methy)  ## Fenotipos (edad, sexo, status)

## Cargamos los datos que usaremos para el diseño leave-one-out
data(GRset)
GRset
assays(GRset)$Beta[1:5,1:5]  ## Betas (nivel de metilación)
pData(GRset)  ## Fenotipos (edad, sexo, status)

## La normalizacion debe hacerse de manera independiente
## pero hemos implementado una normalizacion conjunta que se 
## puede hacer de la siguiente manera (costoso computacionalmente):
## NO CORRER 
## Dataset_procesado <- epi_preprocess(directorio,
##                                   panel_referencia,
##                                   pattern = "SampleSheet.csv")


## Esta función nos permite cambiar los parámetros de cualquiera de los métodos.
## En la user's guide hay una tabla de referencia con la explicación de los
## parámetros modificables. Con este comando vemos los valores actuales:
epi_parameters()


## Cambiamos el p-valor umbral que determina qué regiones son
## outliers en el método manova:
parameters <- epi_parameters(manova = list("pvalue_cutoff" = 0.01))
parameters


## Diseño caso-control
##    1. Separamos los casos de los controles
case_samples <- methy[,methy$status == "case"]
control_samples <- methy[,methy$status == "control"]


##    2. Corremos la función epimutations() con el método manova
epi_mvo <- epimutations(case_samples, 
                        control_samples, 
                        method = "manova")


## Diseño leave-one-out
##    1. Corremos la función epimutations_one_leave_out() con el método manova
epi_mvo_one_leave_out<- epimutations_one_leave_out(GRset,
                                                   method = 'manova')



## Accedemos a los resultados:

## Dimensiones de la tabla - 51 líneas(epimutaciones) x 16 columnas (características)
dim(epi_mvo)
## Vemos las características de las primeras 12 epimutaciones
head(epi_mvo, 12)

## Anotación de las epimutaciones, para añadir información extra
rst_mvo <- annotate_epimutations(epi_mvo)


## Graficamos una epimutación en concreto (en este caso, la 15)
plot_epimutations(epi_mvo[epi_mvo$epi_id == 'epi_manova_1',], methy)


## Graficamos la 1a epimutación con anotación de genes
plot_epimutations(as.data.frame(epi_mvo[epi_mvo$epi_id == 'epi_manova_1',]), 
                  methy = methy, 
                  genes_annot = TRUE)



## Graficamos la 1a epimutación con anotación de metilación y acetilación
## H3K4Me3, H3K27Ac y H3K27Me3
plot_epimutations(as.data.frame(epi_mvo[epi_mvo$epi_id == 'epi_manova_1',]), 
                  methy =  methy, 
                  regulation = TRUE)
