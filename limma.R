#      16/11/2020 ------------   13:45

# trabalhando com Dados brutos

##-------------------- Illumina Beadchip---------------------------------------

# Limma Package
##-------------------------------Dados de Microarray---------------------------------------------------------------------

illumina BeadChips

Os BeadChips de genoma completo da Illumina requerem tratamento especial.

As imagens Illumina são digitalizadas pelo software BeadScan e o software BeadStudio ou GenomeStudio da Illumina 
pode ser usado para exportar os resumos de perfis de sonda.

Os perfis de resumo da sonda são arquivos delimitados por tabulação contendo os dados de intensidade.

Normalmente, todos os arrays processados de uma vez são gravados em um único arquivo, com várias colunas correspondentes a cada array.

Recomendamos que as intensidades sejam exportadas do GenomeStudio.sem correção de fundo ou normalização, pois essas etapas de pré-processamento podem ser melhor realizadas por
funções de limma.O GenomeStudio também pode ser solicitado a exportar perfis para as sondas de controle, e nós recomendamos que isso seja feito também.

Os arquivos Illumina diferem de outras plataformas porque cada arquivo de saída de imagem contém dados de várias matrizes e em que as intensidades para as sondas de controle são gravadas em um arquivo separado das sondas regulares.

Existem outros recursos desses arquivos que podem ser usados opcionalmente para pré-processamento e filtragem.Os arquivos de resumo do probe Illumina podem ser lidos pela função read.ilmn.


> x <- read.ilmn("probe profile.txt", ctrlfiles="control probe profile.txt")

onde 'probe profile.txt' é o nome do arquivo que contém perfis  das sondas principais. 

'control probe profile.txt' é o nome do arquivo que contém perfis para sondas de controle.


Se houver vários perfis de resumo de sonda a serem lidos e as amostras forem resumidas em um alvo, então a 
função read.ilmn.targets pode ser usada.



Ler os perfis da sonda de controle é opcional, mas recomendado.Se os perfis da sonda de controle são disponíveis, 
então os dados do Illumina podem ser corrigidos favoravelmente em segundo plano e normalizados usando a função  neqc ou nec.

Caso contrário, os dados do Illumina são corrigidos em segundo plano e normalizados como para outras plataformas de canal único.



##-------------------Normalization---------------------------------------------------------
            ######################
            #    Within Arrays   #
            #####################

Os métodos podem ser amplamente classificados em métodos que normalizam os valores M para cada matriz separadamente (normalização dentro da matriz)
e métodos que normalizam intensidades ou relações logarítmicas para serem comparáveis entre matrizes (entre matrizes)
Esta seção discute principalmente a normalização dentro da matriz, que geralmente é necessário para a análise de razão de log tradicional de dados de duas cores.

Normalização de loess da ponta de impressão é o método de normalização padrão e pode ser realizada por
> MA <- normalizeWithinArrays(RG)

Existem alguns casos notáveis em que isso não é apropriado. Por exemplo, os arrays Agilent não têm
grupos print-tip, então deve-se usar a normalização de loess global em vez disso:

> MA <- normalizeWithinArrays(RG, method="loess")  
  
O loess da ponta de impressão também não é confiável para matrizes pequenas com menos de, digamos, 150 pontos por grupo de ponta de impressão.
Mesmo matrizes maiores podem ter grupos de ponta de impressão particulares que são muito pequenos para normalização de loess de ponta de impressão se o número de pontos com valores M não ausentes for pequeno para uma ou mais das pontas de impressão
grupos. Nestes casos, deve-se usar normalização global "loess" ou então usar spline robusta
normalização:

> MA <- normalizeWithinArrays(RG, method="robustspline")

A normalização de Loess assume que a maior parte das sondas na matriz não são diferencialmente
expresso. Não pressupõe que existam números iguais de genes regulados para cima e para baixo ou
essa expressão diferencial é simétrica em relação a zero, desde que o ajuste de loess seja implementado em um
forma robusta, mas é necessário que haja um corpo substancial de sondas que não mudam
níveis de expressão. Oshlack et al [23] mostram que a normalização de loess pode tolerar até cerca de 30%
expressão diferencial assimétrica ao mesmo tempo que dá bons resultados.
Esta suposição pode ser suspeita
para matrizes de boutique onde o número total de genes únicos na matriz é pequeno, digamos menos que
150, particularmente se esses genes foram selecionados para serem expressos especificamente em um dos RNA
fontes. Em tal situação, a melhor estratégia é incluir nos arrays uma série de
pontos de controle expressos, como uma série de titulação de pontos de pool de toda a biblioteca, e para usar o método de aumento de peso discutido abaixo [23]. Um pool de biblioteca inteira significa que se faz um pool de um
biblioteca de sondas e imprime manchas do pool em várias concentrações [52]. A biblioteca deve
ser suficientemente grande para que possamos ter certeza de que a média de todas as sondas não é diferencialmente
expresso. Quanto maior a biblioteca, melhor. Bons resultados foram obtidos com pools de biblioteca
com apenas 500 clones. Na ausência de tais pontos de controle, a normalização de matrizes boutique
requer aconselhamento especializado.

Quaisquer pesos de qualidade pontuais encontrados em RG serão usados na normalização por padrão. Isso significa
por exemplo, que pontos com peso zero (sinalizados) não irão influenciar a normalização de outros
pontos. O uso de pesos de qualidade pontuais não resultará, no entanto, na remoção de quaisquer manchas do
objeto de dados. Mesmo pontos com peso zero serão normalizados e aparecerão no objeto de saída,
tais pontos simplesmente não terão qualquer influência sobre os outros pontos. Se você não deseja a qualidade do local
pesos a serem usados na normalização, seu uso pode ser anulado usando:

> MA <- normalizeWithinArrays(RG, weights=NULL)

O objeto de saída MA ainda conterá quaisquer pesos de qualidade spot encontrados em RG, mas esses pesos são
não usado na etapa de normalização.

Freqüentemente, é útil fazer uso de pontos de controle para auxiliar no processo de normalização. Por exemplo,
se as matrizes contêm uma série de pontos que são conhecidos de antemão por serem expressos de forma não diferencial,
esses pontos podem receber mais peso no processo de normalização. Pontos que são conhecidos com antecedência
para ser expresso diferencialmente pode ser reduzido. Suponha, por exemplo, que controlStatus ()
tem sido usado para identificar pontos de pico que são expressos diferencialmente e uma série de titulação de
manchas da biblioteca inteira que não devem ser expressas diferencialmente. Então, pode-se usar

> w <- modifyWeights(RG$weights, RG$genes$Status, c("spikein","titration"), c(0,2))
> MA <- normalizeWithinArrays(RG, weights=w)

para dar peso zero aos pontos de pico e peso duplo aos pontos de titulação. Este processo é
automatizado pelo método de normalização "controle", por exemplo

> csi <- RG$genes$Status=="titration"
> MA <- normalizeWithinArrays(RG, method="control", controlspots=csi)

Em geral, csi é um vetor de índice que especifica os pontos de controle expressos de forma não diferencial [23].
A ideia de aumentar a ponderação dos pontos de titulação segue o mesmo espírito da normalização composta
método proposto por [52], mas é mais flexível e geralmente aplicável. O código acima assume que
RG já contém pesos de qualidade pontuais. Se não, pode-se usar

> w <- modifyWeights(array(1,dim(RG)), RG$genes$Status, c("spikein","titration"), c(0,2))
> MA <- normalizeWithinArrays(RG, weights=w)


Limma contém alguns métodos de normalização mais sofisticados. Em particular, alguns métodos de normalização entre matrizes são discutidos na Seção 6.3 deste guia.
=============================================================================================================================
----------------------------------------------------------------------------------------------------------------------
                    ######################################
                    #    
                        between-Array Normalization    #
                    #                                   #      
                    ######################################


Esta seção eplora alguns dos métodos disponíveis para normalização entre matrizes de duas cores.  

Um recurso que distingue a maioria desses métodos da normalização dentro da matriz é o enfoque 
nos valores individuais de intensidade de vermelho e verde, e não apenas nas relações logarítmicas.

Estes métodos podem, portantox, ser chamados de canais individuais ou métodos de normalização de canais separados.
A normalização de canal individual é normalmente um pré-requisito para métodos de análise de canal individual como o fornecido por 'lmscFit ()'.

> load("Apoa1.RData")


Uma questão importante a considerar antes de normalizar entre matrizes é como a correção de background foi realizada.
Para que esta normalização entre arranjos seja eficaz é importante evitar valores  em razões logarítmicas que podem surgir de intensidades corrigidas negativas ou zero.

A função backgroundCorrect () fornece várias opções úteis. Para os fins desta seção, os dados foram foi corrigido usando o método "mínimo":
  
  
> RG.b <- backgroundCorrect (RG, método = "mínimo")



plotDensities exibe densidades empíricas suavizadas para os canais verdes e
vermelhos individuais
em todas as matrizes. Sem qualquer normalização, há uma variação considerável entre os dois canais
e entre matrizes:
  
> plotDensities(RG.b)


Após a normalização de loess dos valores M para cada matriz, as distribuições de
vermelho e verde tornam-se essencialmente o mesmo para cada array, embora ainda
haja uma variação considerável entre os arrays:
  
> MA.p <-normalizeWithinArrays(RG.b)
> plotDensities(MA.p)


A normalização de Loess não afeta os valores A. Aplicando normalização de quantis 
aos valores A torna as distribuições essencialmente as mesmas em matrizes,
bem como em canais:
  
  > MA.pAq <- normalizeBetweenArrays(MA.p, method="Aquantile")
> plotDensities(MA.pAq)

Aplicar a normalização de quantis diretamente às intensidades individuais de 
vermelho e verde produz um resultado semelhante, mas é um pouco mais barulhento:
  
  > MA.q <- normalizeBetweenArrays(RG.b, method="quantile")
> plotDensities(MA.q, col="black")


Existem outros métodos de normalização entre matrizes não explorados aqui. 
Por exemplo normalizeBetweenArrays com o método = "vsn" fornece uma interface 
para os métodos de normalização de estabilização de variância do vsn
pacote.

======================================================================================================================
  
  Comparando Populações de Células Progenitoras Mamárias com Ilumina BeadChips

Os arquivos de dados de expressão são fornecidos como arquivos gzip e precisarão
ser descompactados antes que possam ser usados para uma análise de limma.
  
->The target RNA samples

passo 1 ->  Lendo o arquivo que descreve as amostras de rna        
  
library(limma)

targets <- readTargets("NOME DO ARQUIVO")
targets


-> The expression profiles

As amostras de RNA foram hibridizadas com dois BeadChips Illumina HumanWG-6 v 3. 
Cada conta-O chip pode acomodar seis amostras.
As imagens do BeadChips foram digitalizadas e resumidasusando o BeadStudio. 
Perfis de teste resumidos não normalizados foram exportados do BeadStudio
para arquivos de texto delimitados por tabulação


Perfis de expressao nao normalizados foram exportados do BeadStudio para arquivos .txt  

Arquivos separados foram escritos para probes regulares e para probes de controle. 
O arquivo de sonda profile.txt contém os perfis de expressão para testes 
regulares, projetados para interrogar os níveis de expressão degenes.
control probe profile.txt contém os perfis das sondas de controle, 
incluindo controle negativos.

Observe que o BeadStudio por padrão grava os perfis de todas as amostras nos 
mesmos dois arquivos  

Lemos nos perfis de expressão para sondas regulares e de controle, informando 
read.ilm que nós desejamos ler os valores p de detecção, bem como os valores de 
expressão:  

x<- read.ilmn(files="probe profile.txt",ctrlfiles="control probe profile.txt",
                   + other.columns="Detection")


Este comando lê um objeto EListRaw . Existem cerca de 750 sondas negativas e cerca
de 49.000 sondas regulares:
  
  
  O componente E contém o valor da expressão para cada sonda

> options(digits=3)
> head(x$E)


As intensidades variam de cerca de 5 a 14 na escala log 2 :
  

Quantas sondas são verdadeiramente expressas?
  
Os valores de detecção contêm valores p para testar se cada sonda é mais 
intensa do que o negativosondas de controle. Valores pequenos são evidências de 
que a sonda corresponde a um gene verdadeiramente expresso:
  

Podemos ir além disso e estimar a proporção geral das sondas regulares que
correspondempara transcrição expressa, usando o método de Shi et al

> pe <- propexpr(x)
> dim(pe) <- c(4,3)
> dimnames(pe) <- list(CellType=c("MS","Stroma","ML","LP"),Donor=c(1,2,3))
> pe

Normalization and filtering

Background correction and normalize:
 
 > y <- neqc(x)

As funções neqc executam a correção de fundo de normexp usando controles 
negativos e, em seguida, normaliza os quantis e, finalmente, 
log 2 transformações [36]. Ele também remove automaticamente as sondas de 
controle, deixando apenas as sondas regulares em y 

dim (y)


Filtrando as sondas que não são expressas. Mantemos sondas que são expressas em 
pelo menos três matrizes de acordo com valores de p de detecção de 5%:

> expressed <- rowSums(y$other$Detection < 0.05) >= 3
> y <- y[expressed,]
> dim(y)
[1] 24691 12



Differential expression between cell types


Agora procuramos genes expressos diferencialmente. Fazemos todas as comparações 
de pares possíveis entreos tipos de células epiteliais, permitindo a correlação 
dentro dos doadores:

> fit <- lmFit(y,design,block=targets$Donor,correlation=dupcor$consensus.correlation)
> contrasts <- makeContrasts(ML-MS, LP-MS, ML-LP, levels=design)
> fit2 <- contrasts.fit(fit, contrasts)
> fit2 <- eBayes(fit2, trend=TRUE)
> summary(decideTests(fit2, method="global"))


Dez principais sondas expressas diferencialmente entre ML e MS

topTable(fit2, coef=1)

Genes de assinatura para células progenitoras luminais

Agora encontramos genes expressos exclusivamente em células LP, em comparação 
com MS e ML. Reformamos o linearmodelo, tornando LP o tipo de célula de referência:

> ct <- relevel(ct, ref="LP")
> design <- model.matrix(~ct)
> fit <- lmFit(y,design,block=targets$Donor,correlation=dupcor$consensus.correlation)
> fit2 <- fit[,c("ctMS","ctML")]
> fit2 <- eBayes(fit2, trend=TRUE)

Então encontramos todos os genes que são regulados positivamente em LP vs MS e 
ML, usando uma alteração de 2 vezese 5% FDR:
  
> results <- decideTests(fit2, lfc=1)
> vennDiagram(results, include=c("up","down"))  
  
  
Existem 197 genes de assinatura positivos e 413 negativos. Para ver os 
principais genes de assinatura positivacom suas alterações de dobra
  
> LP.sig <- rowSums (resultados> 0) == 2
> topTable (fit2 [LP.sig,])

================================================================================
using Agilent

Data availability

Todos os arquivos dos arquivos foram baixados de 
http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-33005.
Os dados também estão disponíveis como série GEO GSE33005.

Usanfo R, os dados podem ser baixados para o seu espaço de trabalho usando:
  
> URL <- "https://www.ebi.ac.uk/arrayexpress/files/E-GEOD-33005"
> SDRF.file <- "E-GEOD-33005.sdrf.txt"
> Data.file <- "E-GEOD-33005.raw.1.zip"
> download.file(paste(URL,SDRF.file,sep="/"), SDRF.file)
> download.file(paste(URL,Data.file,sep="/"), Data.file)
> unzip(Data.file)

Reading the data

Primeiro, lemos o arquivo de amostra e formato de relacionamento de dados (SDRF).
Isso é equivalente ao que éconhecido como "arquivo de destino" no limma

> SDRF <- read.delim("E-GEOD-33005.sdrf.txt",check.names=FALSE,stringsAsFactors=FALSE)
> SDRF[,c("Array Data File","Characteristics[treatment]")]



Estamos interessados na coluna de tratamento:
  
> Treatment <- SDRF[,"Characteristics[treatment]"]

Definimos o controle de solução salina para ser o primeiro nível
do fator de tratamento

> levels <- c("10 ml/kg saline","2 ml/kg corn oil","5 ml/kg corn oil","10 ml/kg corn oil")
> Treatment <- factor(Treatment, levels=levels)

Em seguida, leia os dados de intensidade:
> x <- read.maimages (SDRF [, "Array Data File"],source = "agilent", 
                      green.only = TRUE, other.columns = "gIsWellAboveBG")


Observe que lemos na coluna extra gIsWellAboveBG , que registra se a intensidade 
decada ponto é considerado acima do nível de fundo dessa matriz. Esta coluna nos ajudará mais tardecom filtragem de sonda

Os dados têm 44.254 sondas e 19 matrizes:
  
>dim (x)

Gene annotation

We can use the annotation package for this type of Agilent array, RnAgilentDesign028282.db, to get
gene symbols and Entrez Gene Ids from the probe Ids:

> library(RnAgilentDesign028282.db)
> x$genes$EntrezID <- mapIds(RnAgilentDesign028282.db, x$genes$ProbeName,
                              keytype="PROBEID", column="ENTREZID")
> x$genes$Symbol <- mapIds(RnAgilentDesign028282.db, x$genes$ProbeName,
     keytype="PROBEID", column="SYMBOL")

> x$genes[201:205,]
Background correction and normalize
We use normexp background correction followed by quantile normalization:

> y <- backgroundCorrect(x, method="normexp")  
> y <- normalizeBetweenArrays(y, method="quantile")


Gene filtering
We will filter out control probes as indicated by the ControlType column:

  > Control <- y$genes$ControlType==1L


We will also filter out probes with no Entrez Gene Id or Symbol
> NoSymbol <- is.na(y$genes$Symbol)

Finally, we will filter probes that don’t appear to be expressed. We keep probes
that are above background on at least four arrays (because there are four 
 replicates of each treatment):
  
  IsExpr <- rowSums(y$other$gIsWellAboveBG > 0) >= 4
Now we select the probes to keep in a new data object yfilt:
  > yfilt <- y[!Control & !NoSymbol & IsExpr, ]
> dim(yfilt)
[1] 27191 19
  
Para ficar bem organizado, removemos colunas de anotação de que não precisamos mais:
  
> yfilt$genes <- yfilt$genes[,c("ProbeName","Symbol","EntrezID")]
> head(yfilt$genes)


Differential expression

Now we can find genes differentially expressed for the corn oil treatments compared to the saline
control:
  
  > design <- model.matrix(~Treatment)
> fit <- lmFit(yfilt,design)
> fit <- eBayes(fit,trend=TRUE,robust=TRUE)
> summary(decideTests(fit[,-1]))


Gene ontology analysis

> g <- goana(fit, coef=4, species="Rn", geneid="EntrezID")
> topGO(g,n=20,truncate="50")

==============================================================================================


  
-------------------------------------------------------------------------------------------------



Capitulo 5-> Analise de Qualidade

Uma etapa essencial na análise de quaisquer dados de microarray é verificar
a qualidade dos dados das matrizes. 

Para dados de matriz de duas cores, uma etapa essencial é visualizar os gráficos MA dos dados não normalizados
para cada matriz.

A função plotMD () produz gráficos para matrizes individuais.

A funcão plotMA3by2 () fornece uma maneira fácil de produzir gráficos MA para todos os arrays em um grande experimento. 







-------------------------------------------------------------------------------------------------------------------------





