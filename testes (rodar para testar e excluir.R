##-------------------GEODATASETS------------------------------------------------
#
#
#
# 
================================================================================
##--------------------Obtendo dados do Geo--------------------------------------

## Chamando as librarys:

library(Biobase)
library(GEOquery)

--------------------------------------------------------------------------------
## Obtendo os dados                                 #process

GSE42023 <- getGEO("GSE42023", GSEMatrix = TRUE,AnnotGPL = TRUE)

GSE42023 <- GSE42023[[1]]

## Entendendo o tipo de dados baixado

class(GSE42023)
length(GSE42023)
names(GSE42023)

# Transformando o Dado em um Elist

GSE42023 <- GSE42023[[1]]



# Obtendo os dados Suplementares               #raw data

supp_data <- getGEOSuppFiles("GSE42023")


# Gerando um dataframe com informações dos dados suplementares

rownames(supp_data) <- basename(rownames(supp_data))

supp_data

# os dados brutos vem em formato tar, para vêr o caminho de Download no R use: 

tarArchive <- rownames(supp_data)[1]
tarArchive

# a partir do caminho fornecido, encontre os arquivos em seu computador.
# É importante conhecer o tipo de plataforma utilizada na obtenção dos dados 
# que deseja analizar. 
#Normalmente, em caso de ensaios ILLUMINA BEADCHIP, 
# Geo fornece 2 arquivos suplementares, raw.tar ( geralmente contém 
#dados do experimento)

# non-normalize.txt (fornece indices de expressão genica e pvalores de detecção)

--------------------------------------------------------------------------------

#--------------------Target file------------------------------------------------
# Conhecendo as variáveis  e montando o desenho experimental de nosso estudo

names(pData(GSE42023)) 

# Criando um dataframe e incluindo uma nova variável 

clin_data <- as.data.frame(GSE42023@phenoData@data)
clin_data$group <- "  "

#Alocando seus pacientes aos grupos que criou 

clin_data[clin_data$`tumor type:ch1`!= 'non-habits associated tongue tumor',41] <- "habits"
clin_data[clin_data$`tumor type:ch1`== 'non-habits associated tongue tumor',41] <- "non_habits";


# Note que 41 , representa a localizaçao da coluna "group" na tabela

# Por fim criamos um fator com os grupos e seu respectivos samples

targets<- as.data.frame(cbind(clin_data$description,clin_data$group))
colnames(targets) <- c('description', 'group')
row.names(targets) <- targets$X
targets$X <-NULL
write.table(targets,"targets",sep = "\t")

--------------------------------------------------------------------------------
##=====================ILLUMINA bead array =====================================

---------------Importando os dados brutos para o ambiente R--------------------

--------------------------Arquivo TXT-------------------------------------------

#  Defina o diretório onde o arquivo poderá ser encontrado
  
setwd("/home/daniela/GSE42023/TESTES")   

# Importando os seus arquivos .txt

# Note que,para o exemplo aqui utilizado (beadchip illumina) precisamos tratar 
#manualmente o arquivo .txt
# Criando 2 arquivos, 1 apenas com p_valores de detecção
#                     1 com os niveis de expressão

require(data.table)

#Dados de expressão (exclui manualmente as colunas p_valor)

non_norm_counts <- read.table("GSE42023_non_normalized_counts.csv",
                              sep = ',',header = TRUE)

non_norm_counts_1 <- non_norm_counts

# P_valor de detecção  (exclui manualmente as colunas de expressão)

non_norm_pvalue <- read.table("GSE42023_non_normalized_pvalue.csv",
                              sep = ',',header = TRUE)

#####################não rodei ainda 17:52#######################
# Defina a coluna ID_REF como rownames

row.names(non_norm_counts_1) <- non_norm_counts_1$Probe_Id

row.names(non_norm_pvalue) <- non_norm_pvalue$Probe_Id

# Excluindo a coluna ID_REF
non_norm_counts_1$Probe_Id <- NULL
non_norm_pvalue$Probe_Id <- NULL

-------------------------------------------------------------------------------
------------------------Arquivo BGX--------------------------------------------
# Lendo o arquivo bgx que se encontra dentro da pasta RAW.TAR
  
library(limma)  

library(illuminaio)


bgxfile <- '/home/daniela/GSE42023/TESTES/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx'

annot_all <- illuminaio::readBGX(bgxfile)
annot<- illuminaio::readBGX(bgxfile)$probes

# se quiser personalizar as colunas   de interesse, use o comando:
annot<-annot[,which(colnames(annot) 
           %in% c('Source','Symbol','Transcript','ILMN_Gene','RefSeq_ID',
        'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
        'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
        'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
                                           'Ontology_Function', 'Synonyms'))]



##-----Alinhando  os dados brutos com o arquivo de Expressão--------------------

## BGX- Expression List


annot. <- annot[which(annot$Probe_Id %in% non_norm_counts$Probe_Id),]

annot. <- merge(non_norm_counts,annot,by.x = "Probe_Id", all.x = TRUE,all.y = FALSE)


##--------------------update the target info------------------------------------

targetinfo <- targets
rownames(targetinfo) <- targetinfo$description

annot.x <- non_norm_counts[which(colnames(non_norm_counts) %in% rownames(targetinfo)),]
annot.x <- non_norm_counts_1[,match(rownames(targets),colnames(non_norm_counts_1)]

if (!all(colnames(non_norm_counts_1) == rownames(targets)))
  stop('Target info is not aligned to expression data - they must be in the same order')

-------------------------------------------------------------------------------
  
# create a custom EListRaw object

project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- targets
project@.Data[[3]] <- annot.
project@.Data[[4]] <- non_norm_counts_1
project@.Data[[5]] <- NULL
project$E <- non_norm_counts_1
project$targets <- targetinfo
project$genes <- annot.
project$other$Detection <- non_norm_pvalue

# perform background correction on the fluorescent intensities


project.bgcorrect.norm <- neqc(project, offset = 16)

# filter out control probes, those with no symbol, and those that failed

annot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(project.bgcorrect.norm) %in% annot$Probe_Id),]
annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot
Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$Detection <= 0.05) >= 3
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)


# remove annotation columns we no longer need
project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Probe_Id',
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Chromosome','Probe_Coordinates','Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Entrez_Gene_ID','Symbol')]
head(project.bgcorrect.norm.filt$genes)

# summarise across genes by mean
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt, ID = project.bgcorrect.norm.filt$genes$Symbol)
dim(project.bgcorrect.norm.filt.mean)






