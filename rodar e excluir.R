##=====================ILLUMINA bead array =====================================
---------------Importando os dados brutos para o ambiente R--------------------
  
-------------------------Arquivo TXT-------------------------------------------
  #  defina o diretório onde o arquivo poderá ser encontrado
  
setwd("/home/daniela/GSE42023/TESTES")   

# Importando os seus arquivos .txt

# Note que,para o exemplo aqui utilizado (beadchip illumina) precisamos tratar 
#manualmente o arquivo .txt
# Criando 2 arquivos, 1 apenas com p_valores de detecção
#                     1 com os niveis de expressão

require(data.table)

#Dados de expressão (exclui manualmente as colunas p_valor)

non_norm_counts <- read.table("GSE42023_non_norm_counts.csv",
                              sep = ',',header = TRUE)

non_norm_counts <- as.matrix(non_norm_counts)  



# P_valor de detecção  (exclui manualmente as colunas de expressão)

non_norm_pvalue <- read.table("GSE42023_non_norm_pvalue.csv",
                              sep = ',',header = TRUE)

non_norm_pvalue <- as.matrix(non_norm_pvalue)
-------------------------------------------------------------------------------
  fazer depois de importar o bgx

# Defina a coluna ID_REF como rownames

row.names(non_norm_counts) <- non_norm_counts$Probe_Id

row.names(non_norm_pvalue) <- non_norm_pvalue$Probe_Id

# Excluindo a coluna ID_REF
non_norm_counts$Probe_Id <- NULL
non_norm_pvalue$Probe_Id <- NULL
-------------------------------------------------------------------------------
probes <- non_norm_counts$Probe_Id
probes <-as.data.frame(probes)
colnames(probes)<- c("Probe_Id")
------------------------Arquivo BGX--------------------------------------------
# Lendo o arquivo bgx que se encontra dentro da pasta RAW.TAR
  
library(limma)  

library(illuminaio)


bgxfile <- '/home/daniela/GSE42023/TESTES/GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx'

annot<- illuminaio::readBGX(bgxfile)$probes

annot<-annot[,which(colnames(annot) %in% c('Source','Symbol','Transcript','ILMN_Gene','RefSeq_ID',
                                           'Entrez_Gene_ID','Symbol','Protein_Product','Probe_Id','Probe_Type',
                                           'Probe_Start','Chromosome','Probe_Chr_Orientation','Probe_Coordinates',
                                           'Cytoband', 'Definition', 'Ontology_Component', 'Ontology_Process',
                                           'Ontology_Function', 'Synonyms'))]



#Alinhando o arquivo bgx 


annot <-merge(non_norm_counts,annot,by.x = "Probe_Id",all.x = TRUE, all.y = FALSE)




# update the target info

targetfile <-'/home/daniela/GSE42023/TESTES/targets'
targetinfo <-readTargets(targetfile,sep = '\t')
rownames(targetinfo) <- targetinfo$description
targetinfo$description <- NULL


annot_test <- match(rownames(targetinfo,colnames(annot)))
annot_test <- annot[match(rownames(targetinfo),colnames(annot))]
  
  if (!all(row.names(targetinfo)== colnames(annot)))
    stop('Target info is not aligned to expression data')

-------------------------------------------------------------------------------
 


--------------------------------------------------------------------------------
##------------------------------------------------------------------------------
## Alinhando os dados em BGX aos arquivos de expressã
#23/11/2020 comando merge (só consegue alinhar dois dataframes por vez )

# Passo 1 os arquivos a serem alinhados precisam estar em formato DF

probes <- bgx$probes
controls <- bgx$controls
counts <- non_norm_counts

#passo 2 unindo os dataframes por colunas

#--> unindo Counts e probes
raw_int <- merge(counts,probes,by.x = "Probe_Id",all.x = TRUE,all.y = FALSE)

#--> unindo Counts e controls   ### não deu certo (pq nao tem genes iguais)
raw_int<- merge(raw_int,controls, by.x = "Probe_Id",all.x = TRUE,all.y = FALSE)

--------------------------------------------------------------------------------
  #preparando os dados:
  # Defina a coluna ID_REF como rownames
  
  row.names(raw_int) <- raw_int$Probe_Id

# Excluindo a coluna ID_REF

raw_int$Probe_Id <- NULL
-------------------------------------------------------------------------------
---------------Conhecendo os dados e definindo extratificação da amostra-------
  
#Conhecendo as variáveis 
  
names(pData(exp_Data)) 

# Criando um dataframe e incluindo uma nova variável 

clin_data <- as.data.frame(exp_Data@phenoData@data)
clin_data$group <- "  "

#Alocando seus pacientes aos grupos que criou 

clin_data[clin_data$`tumor type:ch1`!= 'non-habits associated tongue tumor',41] <- "habits"
clin_data[clin_data$`tumor type:ch1`== 'non-habits associated tongue tumor',41] <- "non_habits";


# Note que 41 , representa a localizaçao da coluna "group" na tabela

# Por fim criamos um fator com os grupos e seu respectivos samples

group <- as.data.frame(cbind(clin_data$description,clin_data$group))

colnames(group) <- c('description', 'group_definition')


group <- paste(group$group_definition, sep=",")
group <- factor(group)
table(group)
levels(group)
group  
--------------------------------------------------------------------------------
  # TRANSFORMANDO OS ARQUIVOS BRUTOS EM Elist
  
  library(limma)

raw_data <- " "
raw_data <- new("EListRaw",raw_data)
raw_data@.Data [[1]] <- 'illumina'
raw_data@.Data [[2]] <- clin_data
raw_data@.Data [[3]]  <-raw_int
raw_data@.Data [[4]]  <- non_norm_counts
raw_data@.Data [[5]]  <-NULL

raw_data$E <- non_norm_counts
raw_data$targets <- cbind(clin_data$description,clin_data$group)
colnames(raw_data$targets) <- c('description', 'group_definition')
row.names(raw_data$targets) <-targets$description
write.table(raw_data$targets,"targets",sep = "\t")
raw_data$genes <-bgx
raw_data$other$Detection <- non_norm_pvalue
------------------------------------------------------------------------------------
  #Background correction / normalisation
  
  library(limma)

raw_data.bgcorrect.norm <-neqc(raw_data, offset = 16) # recomendado pelo professor travis gordon

dim(raw_data.bgcorrect.norm)

boxplot(log2(raw_data.bgcorrect.norm$E),range=0,ylab="log2 intensity")




As funções neqc executam a correção de fundo de normexp usando controles 
negativos e, em seguida, normaliza os quantis e, finalmente, 
log 2 transformações [36]. Ele também remove automaticamente as sondas de 
controle, deixando apenas as sondas regulares em y 
-------------------------------------------------------------------------------=
  #Filtrando as sondas que não são expressas. Mantemos sondas que são expressas em 
  #pelo menos três matrizes de acordo com valores de p de detecção de 5%:
  
  expressed <-rowSums(raw_data.bgcorrect.norm$other$Detection <0.05) >= 3

raw_data.bgcorrect.norm <- raw_data.bgcorrect.norm[expressed,]


dim(raw_data.bgcorrect.norm)

--------------------------------------------------------------------------------
  ############# não deu certo====================================================

#filter out control probes, those with no symbol, and those that failed


Control <- raw_data.bgcorrect.norm$genes$controls=="ILMN_Controls"

NoSymbol <- raw_data.bgcorrect.norm$genes$probes == " "

isexpr <- rowSums(non_norm_pvalue <= 0.05)

raw_data.bgcorrect.norm.filt <-raw_data.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]

project.bgcorrect.norm.filt <- project.bgcorrect.norm[!Control & !NoSymbol & isexpr, ]
dim(project.bgcorrect.norm)
dim(project.bgcorrect.norm.filt)


table(Control)
table(NoSymbol)
table(isexpr)

=================================================================================
  
  # ID is used to identify the replicates
  
  project.bgcorrect.norm.filt.mean <- avereps(project.bgcorrect.norm.filt,
                                              ID = project.bgcorrect.norm.filt$genes$Symbol)
dim(project.bgcorrect.norm.filt.mean)

================================================================================
  
  
  
  #arquivo targets
  
  targets <- cbind(clin_data$description,clin_data$group)
colnames(targets) <-c('description','group')
targets <-as.data.frame(targets)
row.names(targets) <- targets$description
targets$description <-NULL

write.table(targets,"targets",sep = "\t")


