#####USANDO O TCGAbiolinksGUI#####


# O TCGAbiolinksGUI traz uma Interfacie Grafica para o usuário, 
# através dessa ferramenta, os dados são fornecidos de forma mais amigável, e 
#especialmente para aqueles que enfrentam dificuldades em utilizar scripts e 
#linhas de comando no R,  tráz consigo as facilidades de utilizar apenas 2 
#comandos básicos.

# Todavia, por ser uma ferramenta de cunho acadêmico e subsidiada sem apoio 
#financeiro, é susceptivel a alguns bugs e falta de suporte.

##### OBSERVAÇÃO (22/10/2020): por falta de mantenedores e suporte, essa ferramenta está 
#desatualizada, e poderá ser descontinuada.



#Abaixo, seguem alguns tópicos importantes:

#######instalando o biocmanager#######

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")



# instalando o TCGAbiolinksGUI: Use o seguinte comando

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinksGUI", dependencies = TRUE)



##### Utilizando o TCGAbiolinksGUI #####

#Após instalado, para utilizar as funcionalidades do TCGAbiolinksGui é preciso 
#chamar a biblioteca correspondente

library(TCGAbiolinksGUI)
TCGAbiolinksGUI()

#Assim que chamarmos este comando, será aberta em seu navegador de internet,
#uma página html contendo o painel com todas as analises a serem realizadas



# Os comandos para cada análise especifica estão atrelados a interface gráfica,
#assim, não serão necessarios linhas de comandos
