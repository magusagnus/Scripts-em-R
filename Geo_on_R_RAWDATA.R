###############################################################################
#                 Script Original por Daniela Paola                           #
#                 Para checar o trabalho dela acesse o git:                   #
#                 https://github.com/d-paola                                  #
###############################################################################

###---------------Trabalhando com dados do Geo Datasets no R--------------------
                                                                27/11/2020
#Este é um script em fase de escrita e adaptação para trabalhar com dados Brutos, 
#obtidos na plataforma NCBI GEO.
================================================================================
#================Preparando as bibliotecas que iremos usar=====================

#Instalando pacotes (se necessário):
  
#Instalando o Bioconductor manager:

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  
  install.packages("BiocManager")
  
}
# Instalando os pacotes necessários do Bioconductor:
                                  #Descomente se for a sua primeira vez# 

#BiocManager::install("BioBase")
#BiocManager::install("GEOquery")

================================================================================
##----------------------Obtendo os Dados a serem analisados---------------------

# Passo 1: chame as librarys:

library(Biobase)
library(GEOquery)

# Em seguida, para obter os dados , use o seguinte comando:

GSE42023 <- getGEO("GSE42023", GSEMatrix = TRUE,AnnotGPL = TRUE)

# Este comando vai gerar um arquivo tipo Large list, contendo todas as informa-
# cões e dados inerentes ao GSE estudado.

# GSE42023 Refere-se ao banco de "Genome wide expression profiling of anterior
#tongue cancer with no history of tobacco and alcohol use"


#Para prosseguir com nossas analises e conseguir extrair os dados de expressão,
#será necessario transformar o arquivo large.list, em um arquivo expression list,

#Usaremos portanto, o seguinte comando:

GSE42023 <- GSE42023[[1]]

GSE42023

# Esses dados, vêm  pré-processados, pelos próprios autores. Entretanto, para 
# este protocolo, queremos obter os dados brutos, para podermos tratar e 
# prosseguir com nossas análises

# Os arquivos de dados brutos, normalmente estão alocados nos dados suplementares:
# Para obtê-los:

supp_data <- getGEOSuppFiles("GSE42023")

# Gerando um dataframe com informações dos dados suplementares

rownames(supp_data) <- basename(rownames(supp_data))

supp_data

# os dados brutos vem em formato tar, para vêr o caminho de Download no R use: 

tarArchive <- rownames(supp_data)[1]
tarArchive

# A partir do caminho fornecido, encontre os arquivos em seu computador.

# É importante conhecer o tipo de plataforma utilizada na obtenção dos dados 
# que deseja analizar. 

# Normalmente, em caso de ensaios ILLUMINA BEADCHIP, 

# GEO fornece 2 arquivos suplementares:

#raw.tar ( geralmente contém dados do experimento)
#non-normalize.txt (fornece indices de expressão genica e pvalores de detecção)

================================================================================

##------------- Obtendo dados clinicos, criando e categorizando grupos----------

# Conhecendo as variáveis  e montando o desenho experimental de nosso estudo:

names(pData(GSE42023)) 

# Criando um dataframe e incluindo uma nova variável 

clin_data <- as.data.frame(GSE42023@phenoData@data)
clin_data$group <- "  "

#Alocando seus pacientes aos grupos que criou: 

clin_data[clin_data$`tumor type:ch1`== 'non-habits associated tongue tumor',41] <- "non_habits"
clin_data[clin_data$`tumor type:ch1`!= 'non-habits associated tongue tumor',41] <- "habits";


# Note que 41 , representa a localizaçao da coluna "group" na tabela

targets<- as.data.frame(cbind(clin_data$description,clin_data$group))
colnames(targets) <- c('description', 'group')
row.names(targets) <- targets$description
targets$description <-NULL
write.table(targets,"targets",sep = "\t")

================================================================================
# Para prosseguir com as análises, é preciso entender qual o perfil dos dados
# iremos trabalhar (.TXT , .CEL)
# Na propria pagina do GEO, voce conseque todas estas informacoes:
  
================================================================================
##                  Pré -processamento de Dados
##=====================ILLUMINA bead array =====================================
# Normalmente, os ensaios Illumina, alocados no Geo, fornecem arquivos padrão:

#Arquivo .txt: geralmente contem os dados de expressao bruta, apresentando 
#genes nas linhas, contagens brutas e pvalores de detecção  (de cada amostra)
#nas colunas


#Arquivo BGX: geralmente contem as informaçoes do ensaio, dados completos de 
# todas as sondas utilizadas ( probes e controles) ,alem de dados de anotação.

##-------------Importando os dados brutos para o ambiente R---------------------
--------------------------Arquivo TXT-------------------------------------------
  
#  Defina o diretório onde o arquivo poderá ser encontrado
  
setwd("/home/daniela/GSE42023/TESTES")   

# Importando os seus arquivos .txt

# Note que,para o exemplo aqui utilizado (beadchip illumina) precisamos tratar 
# manualmente o arquivo .txt
# Criando 2 arquivos, 1 apenas com p_valores de detecção
#                     1 com os niveis de expressão

require(data.table)

#Dados de expressão (exclui manualmente as colunas p_valor)

non_norm_counts <- read.table("GSE42023_non_normalized_counts.csv",
                              sep = ',',header = TRUE)


# P_valor de detecção  (exclui manualmente as colunas de expressão)

non_norm_pvalue <- read.table("GSE42023_non_normalized_pvalue.csv",
                              sep = ',',header = TRUE)

# Defina a coluna ID_REF como rownames

row.names(non_norm_counts) <- non_norm_counts$Probe_Id

row.names(non_norm_pvalue) <- non_norm_pvalue$Probe_Id

# Excluindo a coluna ID_REF
non_norm_counts$Probe_Id <- NULL
non_norm_pvalue$Probe_Id <- NULL

------------------------Arquivo BGX--------------------------------------------
# Lendo o arquivo bgx que se encontra dentro da pasta RAW.TAR
  
library(limma)  

library(illuminaio)


bgxfile <-'/home/daniela/GSE42023/TESTES/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx'


annot_all <- illuminaio::readBGX(bgxfile)
# Usando este comando, geramos um arquivo do tipo Elist, que nos fornece dados 
# de probes e controles de maneira separada.

annot<- illuminaio::readBGX(bgxfile)$probes
# Usando este comando chamamos apenas a planilha contendo os probes.

# Se quiser personalizar as colunas   de interesse, use o comando:
annot<-annot[,which(colnames(annot) 
            %in% c('Source','Symbol','Transcript','ILMN_Gene','RefSeq_ID',
        'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
        'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
        'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
        'Ontology_Function', 'Synonyms'))]
--------------------------------------------------------------------------------

##-----Alinhando  os dados brutos com o arquivo de Expressão--------------------

# Nesta etapa vamos alinhar  e unir os arquivo bgx e  txt:

annot <- annot[which(annot$Probe_Id %in% non_norm_counts$Probe_Id),]
annot <- merge(non_norm_counts,annot,by.x = "Probe_Id", all.x = TRUE,all.y = FALSE)

##--------------------Alinhando os arquivos de design---------------------------

#Nesta etapa, iremos alinhar o arquivo target e o arquivo de expressão, essa
# etapa é extremamente necessaria para a fase de analise de dados, especialmente 
# para a analise DEA.
# As amostras em Target precisam estar perfeitamente alinhadas com o arquivo de
# expressao.

targetinfo <- read.table("/home/daniela/GSE42023/TESTES/targets",sep = "\t", header = TRUE)
rownames(targetinfo) <- targetinfo
targetinfo$X <- NULL

targetinfo <- targets

annot.x <- non_norm_counts[which(colnames(non_norm_counts) %in% rownames(targetinfo)),]
if (!all(colnames(annot.x) == rownames(targetinfo)))
  stop('Target info is not aligned to expression data - they must be in the same order')
================================================================================
##---------Criando um Elist RAW (de dados brutos), personalizado----------------
# Este Elist, permite que todos os meus dados sejam alocados em apenas um objeto
# permitindo-me realizar as análises seguintes utilizando dados de apenas um lugar.

project <- new('EListRaw')
project@.Data[[1]] <- 'illumina'
project@.Data[[2]] <- targetinfo
#project@.Data[[3]] <- annot (inclua seu arquivo de anotação)
project@.Data[[3]] <- NULL
project@.Data[[4]] <- non_norm_counts # inclua aqui seu arquivo de expressão bruta
project@.Data[[5]] <- NULL
project$E <- non_norm_counts # arquivo de expressao bruta
project$targets <- targetinfo
#project$genes <- annot
project$genes <- NULL
project$other$Detection <- non_norm_pvalue #arquivo de pvalores

================================================================================
##----------Correção de fundo e  Normalização----------------------------------

# Para BeadArrays, a correção e normalização do plano de fundo são tratadas 
#por uma única função: neqc ()

# Realizar a  correção de fundo nas intensidades fluorescentes 'normexp' é 
# benéfico porque não resulta em valores negativos, o que significa que nenhum 
#dado é perdido para obter o deslocamento ideal (offset), consulte a resposta
#de Gordon Smyth, aqui: https://stat.ethz.ch/pipermail/bioconductor/2006-April/012554.html

#Normalize os dados com o método 'quantil', para ser consistente com "RMA"
#para matrizes Affymetrix


#As funções neqc executam a correção de fundo de normexp usando controles 
#negativos e, em seguida, normaliza os quantis e, finalmente, faz as transformações
#de LOG2. Ele também remove automaticamente as sondas de  controle, deixando 
# apenas as sondas regulares .


project.bgcorrect.norm <- neqc(project, offset = 16)


================================================================================
##------- Filtrando Sondas de controle, sem simbolo e as que falharam----------
  
#Filtrando as sondas que não são expressas. Mantemos sondas que são expressas em 
# pelo menos três matrizes de acordo com valores de p de detecção de 5%:

annot <- annot[which(annot$Probe_Id %in% rownames(project.bgcorrect.norm)),]
project.bgcorrect.norm <- project.bgcorrect.norm[which(rownames(
                            project.bgcorrect.norm) %in% annot$Probe_Id),]

annot <- annot[match(rownames(project.bgcorrect.norm), annot$Probe_Id),]
project.bgcorrect.norm@.Data[[3]] <- annot
project.bgcorrect.norm$genes <- annot
Control <- project.bgcorrect.norm$genes$Source=="ILMN_Controls"
NoSymbol <- project.bgcorrect.norm$genes$Symbol == ""
isexpr <- rowSums(project.bgcorrect.norm$other$Detection <= 0.05) >= 3
project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol
                                                                  & isexpr, ]
dim(project.bgcorrect.norm)

dim(project.bgcorrect.norm.filt)

##--------------- Removendo colunas anotadas que não usaremos------------------- 

project.bgcorrect.norm.filt$genes <- project.bgcorrect.norm.filt$genes[,c(
  'Probe_Id',
  'Definition','Ontology_Component','Ontology_Process','Ontology_Function',
  'Chromosome','Probe_Coordinates','Cytoband','Probe_Chr_Orientation',
  'RefSeq_ID','Entrez_Gene_ID','Symbol')]


##-------- Sumarizando os genes pela média , excluindo réplicas-----------------


project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt
                                ,ID = project.bgcorrect.norm.filt$genes$Symbol)

dim(project.bgcorrect.norm.filt.mean)


===============Processamento de Dados-------------------------------------------
#Differential Expression bethwen groups

##protocolo usersguide
  

targets <- factor(targets$group)
design <- model.matrix(~targets)
colnames(design) <- levels(targets)
fit <- lmFit(project.bgcorrect.norm.filt.mean$E, design)
contrasts<-makeContrasts(habits-non_habits,levels = design)
fit2 <- contrasts.fit(fit,contrasts)
fit2 <-eBayes(fit2,trend = TRUE)
summary(decideTests(fit2,method = "global"))

topTable(fit2,coef = 1)
levels(targets)

results <- decideTests(fit2, lfc = 1)
vennDiagram(results,include = c("up","down"))
==============================================================================
  

