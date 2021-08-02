##Using RToolbox to retrieve data from Firehose 
library(RTCGAToolbox)

##Mostra as rundatas dispon?veis e ativas para serem escolhidas
lastRunDate <- getFirehoseRunningDates()

##Seleciona uma data para se fazer an?lise
lastAnalyseDate <- getFirehoseAnalyzeDates(1)

# get DNA methylation data, RNAseq2 and clinical data for GBM
gbm.data <- getFirehoseData(dataset = "GBM",
                            runDate = lastAnalyseDate, gistic2_Date = getFirehoseAnalyzeDates(1),
                            Methylation = TRUE, Clinic = TRUE, RNAseq2_Gene_Norm = TRUE,
                            fileSizeLimit = 10000)

#Pega os dados que foram baixados acima
gbm.mut <- getData(gbm.data,"Mutation")
gbm.clin <- getData(gbm.data,"clinical") 
gbm.gistic <- getData(gbm.data,"GISTIC")

