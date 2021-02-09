###############################################################################
#                 Script Original por Daniela Paola                           #
#                 Para checar o trabalho dela acesse o git:                   #
#                 https://github.com/d-paola                                  #
###############################################################################

###---------------Trabalhando com dados do Geo Datasets no R--------------------
27/11/2020
Este é um script em fase de escrita e adaptação para trabalhar com dados pré-
processados, obtidos na plataforma NCBI GEO.
================================================================================
  ##================Preparando as bibliotecas que iremos usar=====================

Instalando pacotes (se necessário):
  
  #Instalando o Bioconductor manager:
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    
    install.packages("BiocManager")
    
  }
# Instalando os pacotes necessários do Bioconductor:
#Descomente se for a sua primeira vez# 

#BiocManager::install("BioBase")
#BiocManager::install("GEOquery")

#================================================================================
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


# Para prosseguir com nossas analises e conseguir extrair os dados de expressão,
#será necessario transformar o arquivo large.list, em um arquivo expression list,
#usaremos portanto, o seguinte comando

GSE42023 <- GSE42023[[1]]

GSE42023

##----- --------------Passo 2 obter os dados de expressão-----------------------


# este comando vai gerar  um data frame com os dados de expressao genica,
# sendo as linhas, os genes e as colunas os samples

# expression_data <- as.data.frame(exprs(exp_Data))  salva como Data Frame

expression_data <- (exprs(GSE42023))               # salva como Matriz

dim(expression_data)


##----- Passo 3: Obter os dados clinicos,criar e categorizar grupos-------------


# Dando sequencia nas analises preliminares, precisamos conhecer os dados 
#"clinicos" e informacoes necessarias para definirmos os grupos.
# normalmente, estes dados sao encontrados em phenoData, aqui usado como "pData"

#conhecendo as variaveis de dados existentes e criando um DF:

names(pData(GSE42023))

clin_data <- as.data.frame(GSE42023@phenoData@data)


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


#--------------------Criando o data frame Targets---------------------------------


# Agora, para facilitar a visualizacao as variáveis de interesse, geramos um Df
#"Group" contendo apenas o codigo da amostra e o grupo

targets<- as.data.frame(cbind(clin_data$description,clin_data$group))
colnames(targets) <- c('description', 'group')
row.names(targets) <- targets$description
targets$description <-NULL
write.table(targets,"targets",sep = "\t")

================================================================================
  
