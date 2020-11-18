##### Usando Geo datasets on R


################################################
#  Script Original por Daniela Paola           #
#  Para checar o trabalho dela acesse o git:   #
#  https://github.com/d-paola                  #
################################################

#Que comecem os jogos.


#----------------Preparando as bibliotecas que iremos usar---------------------

#Instalando o Bioconductor manager se necessário

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages("BiocManager")

}
  
# Instalando os pacotes necessários do Bioconductor, Descomente se for a sua primeira vez 

#BiocManager::install("BioBase")
#BiocManager::install("edgeR")
#BiocManager::install("GEOquery")


library(Biobase)
library(GEOquery)
library(edgeR)

##----------------- Passo 1: obter os dados genicos-----------------------------

# Este comando vai gerar um arquivo tipo Large list, contendo todas as informa-
# cões e dados inerentes ao GSE estudado.

# GSE42023 Refere-se ao banco de "Genome wide expression profiling of anterior
#tongue cancer with no history of tobacco and alcohol use"

GSE42023 <- getGEO("GSE42023")

class(GSE42023)
length(GSE42023)
names(GSE42023)

# Para prosseguir com nossas analises e conseguir extrair os dados de expressão,
#será necessario transformar o arquivo large.list, em um arquivo expression list,
#usaremos portanto, o seguinte comando

exp_Data <- GSE42023[[1]]

exp_Data

##----- --------------Passo 2 obter os dados de expressão-----------------------


# este comando vai gerar  um data frame com os dados de expressao genica,
# sendo as linhas, os genes e as colunas os samples

# expression_data <- as.data.frame(exprs(exp_Data))  salva como Data Frame

expression_data <- (exprs(exp_Data))               # salva como Matriz

dim(expression_data)


##----- Passo 3: Obter os dados clinicos,criar e categorizar grupos-------------


# Dando sequencia nas analises preliminares, precisamos conhecer os dados 
#"clinicos" e informacoes necessarias para definirmos os grupos.
# normalmente, estes dados sao encontrados em phenoData, aqui usado como "pData"

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

group <- as.data.frame(cbind(clin_data$geo_accession,clin_data$group))

colnames(group) <- c('geo_accession', 'group_definition')

#agora, vamos fatorizar o Df para gerar as categorias:

group <- paste(group$group_definition, sep=",")
group <- factor(group)
table(group)
group 
--------------------------------------------------------------------------------
##---------------------------12/11/2020-----------------------------------------
# 12:00

##--------------------------getGEOsupp------------------------------------------

# Agora vamos obter os arquivos suplementares, que são arquivos com dados brutos

supp_data <- getGEOSuppFiles("GSE42023")

# gerando um dataframe com informações dos dados suplementares

rownames(supp_data) <- basename(rownames(supp_data))

supp_data

# normalmente, os dados brutos vem em formato tar, para abri-los no R use: 

tarArchive <- rownames(supp_data)[1]
tarArchive

# seguindo o caminho de download, vamos descompactar as pastas e importar os 
#dados para o ambiente R

# converti o arquivo bgx para csv e importei o csv, tem os mesmos dados 
#feature data

#importei o arquivo nao normalizado


-------------------------------------------------------------------------------
#Após importar os arquivos de dados brutos, agora vou proceguir com as análises
  
#Passo 1: criar o df "raw_data" e mover os dados importados pra ela, alem disso
# farei também uma edição na tabela, extraindo apenas os dados brutos.
  
raw_data <-  GPL6883_HumanRef.8_V3_0_R0_11282963_A
  
raw_data_edited <- as.data.frame (raw_data [9:25219,1:28])

colnames(raw_data_edited) <- c("Species","Source","Search_Key","Transcript",
"ILMN_Gene","Source_Reference_ID","RefSeq_ID","Unigene_ID","Entrez_Gene_ID",
"GI","Accession","Symbol","Protein_Product","Probe_Id","Array_Address_Id",
"Probe_Type","Probe_Start","Probe_Sequence","Chromosome","Probe_Chr_Orientation",
"Probe_Coordinates","Cytoband","Definition","Ontology_Component","Ontology_Process",
"Ontology_Function","Synonyms","Obsolete_Probe")
  
--------------------------------------------------------------------------------
 # Agora, usando o outro arquivo de expressão, com dados brutos, vamos gerar
# um Df:
  
exprs_non_normalized <- GSE42023_non_normalized

exprs_non_normalized <- as.matrix(exprs_non_normalized)
exprs_non_normalized <- 
  
  
write.table(exprs_non_normalized,"expressao_dadosbrutos",sep = "\t")

write.table(raw_data_edited,"raw_data",sep = "\t")
---------------------------------------------------------------------------------


library(limma)  
--------------------------------------------------------------------------------
#           16/11/2020  

  
# tentando normalizar os dados brutos
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  














  
  





























##----- Passo 4: Anotar os Dados






##---- Passo 5:

##-----Passo 6:
