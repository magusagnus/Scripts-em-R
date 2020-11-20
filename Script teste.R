#============== 18/11/2020   ==================================================

#teste

#library

library(Biobase)
library(GEOquery)
library(edgeR)
library(limma)

--------------------------------------------------------------------------------
#obtendo dadoos do Geo


GSE42023 <- getGEO("GSE42023", GSEMatrix = TRUE,AnnotGPL = TRUE)
GSE42023
  

class(GSE42023)
length(GSE42023)
names(GSE42023)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-------------------------getGEOsupp---------------------------------------------

# Agora vamos obter os arquivos suplementares, que são arquivos com dados brutos

supp_data <- getGEOSuppFiles("GSE42023")

# gerando um dataframe com informações dos dados suplementares

rownames(supp_data) <- basename(rownames(supp_data))

supp_data

# normalmente, os dados brutos vem em formato tar, para abri-los no R use: 

tarArchive <- rownames(supp_data)[1]
tarArchive
--------------------------------------------------------------------------------
# IMPORTANDO OS DADOS PARA O AMBIENTE R

  # defina o diretório de onde deseja extrair o arquivo

  setwd("/home/daniela/Área de Trabalho/R_GIT_HUB/Scripts-em-R/GSE42023/testes")   

# extraindo o arquivo

require(data.table)

#Apenas dados de expressão (exclui manualmente as colunas p_valor)

non_normalized_counts <- read.table("GSE42023_non_normalized (scan).csv",
                                    sep = ',',header = TRUE)

non_normalized_counts <- as.matrix(non_normalized_counts)  

show(non_normalized_data) 


# p_valor de detecção  (exclui manualmente as colunas de expressão)

non_normalized_pvalue <- read.table("GSE42023_non_normalized (p_value).csv",
                                  sep = ',',header = TRUE)

non_normalized_pvalue <- as.matrix(non_normalized_pvalue)
# é preciso setar a coluna ID_REF como rownames

row.names(non_normalized_counts) <- non_normalized_counts$X

row.names(non_normalized_pvalue) <- non_normalized_pvalue$X

# em seguida, excluir a coluna ID_REF
non_normalized_counts$X <- NULL
non_normalized_pvalue$X <- NULL


--------------------------------------------------------------------------------
  
  #----- criando a matriz target ( dados de design e grupos)-----------------
  

# Dando sequencia nas analises preliminares, precisamos conhecer os dados 
#"clinicos" e informacoes necessarias para definirmos os grupos.
# normalmente, estes dados sao encontrados em phenoData, aqui usado como "pData"

exp_Data <- GSE42023[[1]]

exp_Data


#conhecendo as variaveis de dados existentes e criando um DF:

names(pData(exp_Data))

clin_data <- as.data.frame(exp_Data@phenoData@data)


#criando a variável "group", para alocar informacoes e pacientes:

clin_data$group <- "  "

#alocando os pacientes e dados, para "group"

# Em nosso estudo, definimos dois grupos : "habits"/ "non-habits", para serem
#alocados nestes grupos,  a variável "tumor type" precisa atender aos requisitos:

# "tumor type" igual(==) a non-habits associated tongue tumor-> group "non-habits"

# "tumor type"  diferente (!=) de "non-habits associated..."-> group "habits"

clin_data[clin_data$`tumor type:ch1`== 'non-habits associated tongue tumor',41] <- "non_habits"

clin_data[clin_data$`tumor type:ch1`!= 'non-habits associated tongue tumor',41] <- "habits";

# Note que 41 , representa a localizaçao da coluna "group" na tabela


#--------------------Criando o data frame Group---------------------------------


# Agora, para facilitar a visualizacao as variáveis de interesse, geramos um Df
#"Group" contendo apenas o codigo da amostra e o grupo

group <- as.data.frame(cbind(clin_data$description,clin_data$group))

colnames(group) <- c('description', 'group_definition')

#agora, vamos fatorizar o Df para gerar as categorias:

group <- paste(group$group_definition, sep=",")
group <- factor(group)
table(group)
group  
--------------------------------------------------------------------------------  
#lendo o arquivo bgx de dados brutos

library(illuminaio)
  
my_idat_files <- paste("/home/daniela/Área de Trabalho/R_GIT_HUB/Scripts-em-R/GSE42023/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx")

bgxfile <-my_idat_files

bgx <-readBGX(bgxfile)

--------------------------------------------------------------------------------
# TRANSFORMANDO OS ARQUIVOS BRUTOS EM Elist
 
library(limma)

raw_data <- " "
raw_data <- new("EListRaw",raw_data)
raw_data@.Data [[1]] <- 'illumina'
raw_data@.Data [[2]] <- clin_data
raw_data@.Data [[3]]  <-bgx
raw_data@.Data [[4]]  <- non_normalized_counts
raw_data@.Data [[5]]  <-NULL

raw_data$E <- non_normalized_counts
raw_data$targets <- cbind(clin_data$description,clin_data$group)
raw_data$genes <-bgx
raw_data$other$Detection <- non_normalized_pvalue

--------------------------------------------------------------------------------
  

#20:46 TESTE2
--------------------------------------------------------------------------------
#Background correction / normalisation

library(limma)

raw_data.bgcorrect.norm <-neqc(raw_data, offset = 16) # recomendado pelo professor travis gordon


raw_data.bgcorrect.norm <-neqc(raw_data, offset = 16)

--------------------------------------------------------------------------------
#############3 não deu certo=================================================

#filter out control probes, those with no symbol, and those that failed


Control <- raw_data.bgcorrect.norm$genes$controls=="ILMN_Controls"

NoSymbol <- raw_data.bgcorrect.norm$genes$probes == " "

isexpr <- rowSums(non_normalized_pvalue <= 0.05) >= 3

raw_data.bgcorrect.norm.filt <-raw_data.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]

project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)

dim(teste)


table(Control)
table(NoSymbol)
table(isexpr)

=================================================================================
summarise across genes by mean
# ID is used to identify the replicates
project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                            ID = project.bgcorrect.norm.filt$genes$Symbol)
dim(project.bgcorrect.norm.filt.mean)
