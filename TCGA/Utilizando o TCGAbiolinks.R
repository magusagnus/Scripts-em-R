#####       Utilizando o TCGAbiolinks ######

# No caso desta ferramenta, faremos todas as análises através de linhas de 
# comando, por isto é de suma importancia estabelecer as etapas de todas as 
# análises


#Para usar a ferramenta, é necessário chama-la através dos comandos

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(plyr)
library(limma)
library(biomaRt)

#Quando algum pacote precisar ser instalado, voce pode usar este comando:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")



#Extraindo planilha de dados do TCGA


# Observe que, onde existe "TCGA-HNSC", voce deverá alterar para a sigla 
#relacionada ao cancer de sua escolha

#Ao rodar este comando, observe na aba Enviroment, será gerado um data frame
#entitulado "clinical", nele voce terá acesso às informacoes clinicas dos 
# dos pacientes acometidos pelo tipo de cancer avaliado

# Para fazer o download da planilha importada, use o comando:

write.table(clinical, "DADOS_CLINICOS_TCGA_HNSC.txt", sep="\t")

# Automaticamente, será criado um arquivo txt, que será salvo (geralmente no 
# diretório documentos de seu pc), abra este arquivo e substitua todos os "." 
#por ","

#Para  visualizar os dados possiveis de serem exportatos, use esse comando:

query <- GDCquery(project = "TCGA-HNSC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)


#Em seguida, após o devido carregamento do pacote e dados, chame o arquivo que 
#deseja obter

query <- GDCquery(project = "TCGA-HNSC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type = "clinical_patient_hnsc")
GDCdownload(query)

#Para gerar uma tabela plotavel no R, use o seguinte comando:

clinical.BCRtab.radiation <- GDCprepare(query)

clinical.BCRtab.all$clinical_patient_hnsc  %>% 
  head  %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))


# Finalmente, para escrever esta tabela em formato editavel e txt, use o comando:

write.table(clinical.BCRtab.all$clinical_patient_hnsc, 
            "DADOS_CLINICOS_TCGA_HNSC.txt", sep="\t")


# nota: Assim que o R rodar este comando, um arquivo txt será salvo em 
# seus arquivos (NORMALMENTE NA PASTA DOCUMENTOS)


# é importante salientar, que sempre que usar o comando write.table ,
#é importante mudar o titulo do arquivo que deseja salvar, pois o r sobrescreve 
#os arquivos, e voce pode acabar perdendo dados anteriores


# Referencia: (https://bioconductor.org/packages/release/bioc/vignettes/
#               TCGAbiolinks/inst/doc/clinical.html#Useful_information )



####### Obtendo os dados de expressão gênica ###########

#Agora que você teve acesso aos dados clinicos e definiu quais casos queria 
#analizar, vamos aos proximos passos:



#listSamples <- c(#coloque aqui o barcode dos pacientes que deseja estudar)
  

listSamples <- c("TCGA-BA-4077", "TCGA-BA-6871", "TCGA-BA-A4IG", "TCGA-BB-4225",
                 "TCGA-BB-4228", "TCGA-BB-7861", "TCGA-BB-7871", "TCGA-CN-A498",
                 "TCGA-CN-A6UY", "TCGA-CN-A6V6", "TCGA-CR-5250", "TCGA-CR-6472",
                 "TCGA-CR-6477", "TCGA-CV-5439", "TCGA-CV-6943", "TCGA-CV-6950",
                 "TCGA-CV-7406", "TCGA-DQ-7591", "TCGA-DQ-7593", "TCGA-DQ-7594",
                 "TCGA-F7-A61V", "TCGA-F7-A620", "TCGA-HD-8224", "TCGA-HD-8314",
                 "TCGA-MZ-A6I9", "TCGA-MZ-A7D7", "TCGA-P3-A5QE", "TCGA-BA-5151",
                 "TCGA-CN-4726", "TCGA-CN-4731", "TCGA-CN-4734", "TCGA-CN-A63V",
                 "TCGA-CQ-5334", "TCGA-CQ-6220", "TCGA-CQ-7064", "TCGA-CQ-A4C6",
                 "TCGA-CQ-A4CG", "TCGA-CQ-A4CI", "TCGA-CV-6940", "TCGA-CV-A464",
                 "TCGA-D6-A6EN", "TCGA-DQ-7588", "TCGA-F7-A624", "TCGA-H7-8501",
                 "TCGA-HD-A4C1", "TCGA-IQ-7631", "TCGA-P3-A6T2", "TCGA-QK-A6IG",
                 "TCGA-UF-A7JA", "TCGA-UF-A7JD", "TCGA-BA-5149", "TCGA-BA-5556",
                 "TCGA-BA-6872", "TCGA-BA-A6D8", "TCGA-BA-A6DD", "TCGA-BA-A6DF",
                 "TCGA-BB-8601", "TCGA-CN-4730", "TCGA-CN-5358", "TCGA-CN-5359",
                 "TCGA-CN-5364", "TCGA-CN-5373", "TCGA-CN-6016", "TCGA-CN-6995",
                 "TCGA-CN-A642", "TCGA-CQ-5324", "TCGA-CQ-5332", "TCGA-CQ-6218",
                 "TCGA-CQ-6228", "TCGA-CQ-7068", "TCGA-CQ-7072", "TCGA-CQ-A4C7",
                 "TCGA-CQ-A4C9", "TCGA-CQ-A4CB", "TCGA-CR-6491", "TCGA-CV-5436",
                 "TCGA-CV-6936", "TCGA-CV-6948", "TCGA-CV-6953", "TCGA-CV-7102",
                 "TCGA-CV-7235", "TCGA-CV-7407", "TCGA-CV-A45X", "TCGA-CV-A463",
                 "TCGA-CV-A6JD", "TCGA-CX-7086", "TCGA-CX-7219", "TCGA-CX-A4AQ",
                 "TCGA-D6-A6EO", "TCGA-F7-8489", "TCGA-HD-7832", "TCGA-HD-7917",
                 "TCGA-IQ-A61G", "TCGA-KU-A66T", "TCGA-MT-A67D", "TCGA-MT-A7BN",
                 "TCGA-P3-A6T0", "TCGA-P3-A6T4", "TCGA-P3-A6T7", "TCGA-P3-A6T8",
                 "TCGA-QK-A6II", "TCGA-QK-A6IJ", "TCGA-QK-A6VB", "TCGA-QK-A8Z7",
                 "TCGA-QK-A8Z9", "TCGA-T3-A92N", "TCGA-UF-A719", "TCGA-UF-A71A",
                 "TCGA-UF-A71E", "TCGA-UF-A7JC", "TCGA-UF-A7JO", "TCGA-UF-A7JT",
                 "TCGA-WA-A7GZ", "TCGA-BA-5558", "TCGA-CN-5369", "TCGA-CQ-5331",
                 "TCGA-CQ-7063", "TCGA-CR-6492", "TCGA-CV-5442", "TCGA-QK-A64Z",
                 "TCGA-BB-A5HU", "TCGA-BB-A5HZ", "TCGA-CN-4728", "TCGA-CN-4729",
                 "TCGA-CN-4740", "TCGA-CN-6018", "TCGA-CN-6020", "TCGA-CN-6994",
                 "TCGA-CQ-6227", "TCGA-CQ-7071", "TCGA-CQ-A4CD", "TCGA-CR-6471",
                 "TCGA-CR-6484", "TCGA-CR-7365", "TCGA-CR-7367", "TCGA-CR-7368",
                 "TCGA-CR-7369", "TCGA-CR-7373", "TCGA-CR-7376", "TCGA-CR-7377",
                 "TCGA-CR-7379", "TCGA-CR-7380", "TCGA-CR-7386", "TCGA-CR-7395",
                 "TCGA-CV-5966", "TCGA-CV-6937", "TCGA-CV-6938", "TCGA-CV-6942",
                 "TCGA-CV-6955", "TCGA-CV-6960", "TCGA-CV-7090", "TCGA-CV-7091",
                 "TCGA-CV-7095", "TCGA-CV-7097", "TCGA-CV-7099", "TCGA-CV-7100",
                 "TCGA-CV-7178", "TCGA-CV-7183", "TCGA-CV-7252", "TCGA-CV-7253",
                 "TCGA-CV-7254", "TCGA-CV-7263", "TCGA-CV-7409", "TCGA-CV-7411",
                 "TCGA-CV-7413", "TCGA-CV-7414", "TCGA-CV-7416", "TCGA-CV-7423",
                 "TCGA-CV-7425", "TCGA-CV-7427", "TCGA-CV-7428", "TCGA-CV-7429",
                 "TCGA-CV-7432", "TCGA-CV-7434", "TCGA-CV-7435", "TCGA-CV-7568",
                 "TCGA-CV-A45Q", "TCGA-CV-A45U", "TCGA-CV-A45V", "TCGA-CV-A6JE",
                 "TCGA-CV-A6JN", "TCGA-CV-A6JY", "TCGA-CV-A6JZ", "TCGA-CV-A6K2",
                 "TCGA-CX-7082", "TCGA-H7-7774", "TCGA-H7-8502", "TCGA-HD-A633",
                 "TCGA-HD-A6I0", "TCGA-HL-7533", "TCGA-MT-A67F", "TCGA-P3-A6T3",
                 "TCGA-RS-A6TO", "TCGA-4P-AA8J", "TCGA-BA-4074", "TCGA-BA-4075",
                 "TCGA-BA-5557", "TCGA-BA-6873", "TCGA-BA-7269", "TCGA-BA-A6DB",
                 "TCGA-BA-A6DE", "TCGA-BA-A6DG", "TCGA-BB-4224", "TCGA-BB-7863",
                 "TCGA-BB-7872", "TCGA-BB-A6UO", "TCGA-C9-A47Z", "TCGA-C9-A480",
                 "TCGA-CN-4725", "TCGA-CN-4733", "TCGA-CN-4736", "TCGA-CN-4737",
                 "TCGA-CN-4742", "TCGA-CN-5367", "TCGA-CN-5370", "TCGA-CN-6017",
                 "TCGA-CN-6019", "TCGA-CN-6024", "TCGA-CN-6996", "TCGA-CN-6998",
                 "TCGA-CN-A640", "TCGA-CQ-5325", "TCGA-CQ-5327", "TCGA-CQ-5329",
                 "TCGA-CQ-5330", "TCGA-CQ-5333", "TCGA-CQ-6219", "TCGA-CQ-6221",
                 "TCGA-CQ-6222", "TCGA-CQ-6224", "TCGA-CQ-6225", "TCGA-CQ-6229",
                 "TCGA-CQ-7065", "TCGA-CQ-7067", "TCGA-CQ-A4CA", "TCGA-CQ-A4CE",
                 "TCGA-CQ-A4CH", "TCGA-CR-6488", "TCGA-CR-6493", "TCGA-CR-7372",
                 "TCGA-CR-7382", "TCGA-CR-7390", "TCGA-CR-7391", "TCGA-CR-7392",
                 "TCGA-CR-7393", "TCGA-CR-7394", "TCGA-CR-7397", "TCGA-CR-7401",
                 "TCGA-CV-5970", "TCGA-CV-5971", "TCGA-CV-5973", "TCGA-CV-5976",
                 "TCGA-CV-5977", "TCGA-CV-5979", "TCGA-CV-6003", "TCGA-CV-6433",
                 "TCGA-CV-6436", "TCGA-CV-6441", "TCGA-CV-6933", "TCGA-CV-6934",
                 "TCGA-CV-6939", "TCGA-CV-6941", "TCGA-CV-6945", "TCGA-CV-6951",
                 "TCGA-CV-6952", "TCGA-CV-6954", "TCGA-CV-6956", "TCGA-CV-6959",
                 "TCGA-CV-6961", "TCGA-CV-7103", "TCGA-CV-7104", "TCGA-CV-7180",
                 "TCGA-CV-7236", "TCGA-CV-7238", "TCGA-CV-7243", "TCGA-CV-7255",
                 "TCGA-CV-7438", "TCGA-CV-7446", "TCGA-CV-A45P", "TCGA-CV-A45R",
                 "TCGA-CV-A45T", "TCGA-CV-A465", "TCGA-CV-A6JO", "TCGA-CV-A6JT",
                 "TCGA-CV-A6JU", "TCGA-CV-A6K0", "TCGA-CX-7085", "TCGA-D6-6515",
                 "TCGA-D6-6823", "TCGA-D6-6825", "TCGA-D6-8569", "TCGA-D6-A4Z9",
                 "TCGA-D6-A4ZB", "TCGA-D6-A6EM", "TCGA-DQ-5624", "TCGA-DQ-5625",
                 "TCGA-DQ-5630", "TCGA-DQ-5631", "TCGA-DQ-7592", "TCGA-F7-A50G",
                 "TCGA-F7-A50J", "TCGA-F7-A61S", "TCGA-F7-A61W", "TCGA-H7-A6C4",
                 "TCGA-HD-7831", "TCGA-HD-8634", "TCGA-HD-8635", "TCGA-HD-A6HZ",
                 "TCGA-IQ-A61E", "TCGA-IQ-A61H", "TCGA-IQ-A61J", "TCGA-IQ-A61K",
                 "TCGA-IQ-A61L", "TCGA-IQ-A6SG", "TCGA-IQ-A6SH", "TCGA-KU-A6H8",
                 "TCGA-MT-A51X", "TCGA-MT-A67A", "TCGA-P3-A5QA", "TCGA-QK-A652",
                 "TCGA-QK-AA3K", "TCGA-T2-A6WX", "TCGA-T2-A6WZ", "TCGA-UF-A7JS",
                 "TCGA-UP-A6WW", "TCGA-WA-A7H4", "TCGA-BA-A4IF", "TCGA-BA-A4II",
                 "TCGA-BA-A6DL", "TCGA-BA-A8YP", "TCGA-HD-7753", "TCGA-IQ-7630",
                 "TCGA-IQ-A61I", "TCGA-IQ-A61O", "TCGA-QK-A8ZA")

query.exp <- GDCquery(project = "TCGA-HNSC", 
                      legacy = TRUE,
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type = "results",
                      barcode = listSamples,
                      experimental.strategy = "RNA-Seq",
                      sample.type = c("Primary Tumor","Solid Tissue Normal"))
GDCdownload(query.exp)

#Para exportar o arquivo de expressão .rda , use o comando abaixo)
#normalmente os arquivos exportados são alocados na pasta documentos )

HNSC.exp <- GDCprepare(query = query.exp, save = TRUE,
                       save.filename = "HNSC_selectedExp.rda")

#Com o objetivo de  organizar nossas análises, faremos alguns comandos 
#de organizacao de dados

# get subtype information ( subtipos de tumor)
dataSubt <- TCGAquery_subtype(tumor = "HNSC")

# get clinical data (para trazer os dados clinicos para dentro da análise)
dataClin <- GDCquery_clinic(project = "TCGA-HNSC","clinical") 


# Which samples are Primary Tumor ( para definir quais amostras sao tumor 
#primario)

dataSmTP <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"TP") 

# which samples are solid tissue normal (para definir quais amostras tec.normal)

dataSmNT <- TCGAquery_SampleTypes(getResults(query.exp,cols="cases"),"NT")


####### Para preparar e filtrar os dados######

#preparando os dados :

dataPrep <- TCGAanalyze_Preprocessing(object = HNSC.exp, cor.cut = 0.6)                      

#normalizando os dados:
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                
#filtrando os dados:
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)   



######Análise de Expressão diferencial#######

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

#para exportar a tabela .txt, use o comando abaixo

write.table(dataDEGs, "HNSC_selected-.txt", sep="\t")


# Uma excelente maneira de visualizar  gráficamente a expressão diferencial,
# é atraves do Volcano Plot, para gerar este gráfico, use o comando

TCGAVisualize_volcano(x = dataDEGs$logFC,
                      y = dataDEGs$FDR,
                      filename = "HNSCselected_volcanoexp.png",
                      x.cut = 6,
                      y.cut = 10^-5,
                      names = rownames(dataDEGs),
                      color = c("black","red","darkgreen"),
                      names.size = 2,
                      xlab = " Gene expression fold change (Log2)",
                      legend = "State",
                      title = "Volcano plot (CIMP-high vs CIMP-low)",
                      width = 10)

#NOTE QUE NESTE PARA GERAR O VOLCANO, ESTABELECEMOS COMO PARAMETRO X.CUT=6, 
# ESTE VALOR É ESTABELECIDO PELO PROPRIO PESQUISADOR, UMA VEZ QUE A TITULO DE 
# APRESENTAÇÃO, ELE MOSTRA OS DEGS(genes diferencialmente expressos) a partir do 
# valor de corte definido em "x.cut" ou seja:  x.cut=6 -> o volcano 
#mostra apenas os genes com logFC acima de 6, ou  abaixo de -6.
#isso deve ser explicitado em sua metodologia!!!
 

######Gene ontology #####


# Outra visualização gráfica de suma importancia é a ontologia gênica  a qual
#estes degs estão envolvidos, para nos auxiliar no aprofundamento, analise e 
#discussão dos resultados obtidos em nosso estudo

ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Normal Vs Tumor",
                                RegulonList = rownames(dataDEGs))  

TCGAvisualize_EAbarplot(tf = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = rownames(dataDEGs),
                        nBar = 20)
group1 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
group2 <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))





###### Sobrevida#####
# essa analise nos devolve os genes diretamente relacionados 
# com a sobrevida dos pacientes, entretanto, para esta análise é necessário 
#estratificar os grupos devido à quantidade de genes aceita

dataSurv <- TCGAanalyze_SurvivalKM(clinical_patient = dataClin,
                                   dataGE = dataFilt,
                                   Genelist = rownames(dataDEGs),
                                   Survresult = FALSE,
                                   ThreshTop = 0.67,
                                   ThreshDown = 0.33,
                                   p.cut = 0.05, group1, group2)
require(dnet)
org.Hs.string <- dRDataLoader(RData = "org.Hs.string")


