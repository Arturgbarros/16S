# Workflow 16S por linha de comando

# Baixando arquivos 
* Copiar os arquivos do sequênciamento 16S (HTTPS)
```bash
git clone https://github.com/Tutugb/16S.git
```
* Copiar os arquivos do sequênciamento 16S (SSH)
```bash
git clone git@github.com:Tutugb/16S.git
```
# Criação das pastas
Deixa o workflow é uma das partes mais importantes, pois é com essa organização que não haverá confusão em quais arquivos são os certos em cada etapa por isso é essencial a criação de uma pasta para o armazenamento de cada grupo de arquivos
* Criação de pastas iniciais
```bash
mkdir 16S 16S/00.dadosbrutos 16S/01.qualidade 16S/02.trimmagem 16S/02.trimmagem/00.qualidade 16S/merge 16S/qiime
```
Os arquivos iniciais devem ser inseridos na pasta 00.dadosbrutos
# Controle de qualidade
Primeiramente usa-se o programa **fastqc** como forma de avaliar a qualidade do sequênciamento das amostras onde le os arquivos "fastq" ou "fq" e os converte para uma forma gráfica para melhor interpretação
* Criação do ambiente para análise
```bash
conda create -n qualidade
```
* Ativação do ambiente
```bash
conda activate qualidade
```
* Instalando **fastqc**
```bash
conda install -c bioconda fastqc
```
* Entrar na pasta onde estão os arquivos brutos
```bash
cd 00.dadosbrutos
```
* Executando programa em um arquivo
```bash
fastqc genoma1_1.fq -t 26 -o ../01.qualidade/
```
* Executando programa com **loop**
```bash
for i in {1..28}; do fastqc genoma${i}_1* genoma${i}_2* -o ../01.qualidade/ -t 26 ; done
```
Após a obtenção de todos os arquivos .zip e .html é possível abrir cada um dos 28 arquivos .html e avaliar a qualidade de cada um, porém como isso é um trabalho muito exaustivo e repetitivo, usa-se a ferramenta abaixo para otimmização
## Junção dos arquivos do controle de qualidade
Uma forma de juntar todos os arquivos resultantes do processo do fastqc é através do programa **multiqc** o qual junta todos os arquivos ".html" em um só mostrando todas as análises do fastqc de forma rezumida
*Instalando **multiqc**
pode instalar no ambiente "qualidade"
```bash
conda install -c bioconda multiqc
```
* Executando o programa estando na pasta 16S
```bash
multiqc 01.qualidade/
```
Assim tendo o arquivo "multiqc_data.html" só abrí-lo e avaliar qual onde deve ser feita a trimmagem e quais parâmetros usar para essa e se deve ser ou não retirado adaptadores utilizados no sequênciamento

# Trimmagem
Após avaliação dos parâmetros os quais são necessários serem retirados das sequências é realizada a trimmagem um processo que se consiste em retirar partes das sequências que estejam com a qualidade ruim ou de adapatadores os quais não são necessários para as análises futuras e podem até atrapalhar e para isso usa-se o programa **trimmomatic** e **trim_galore**
* Criação do ambiente para análise
```bash
conda create -n trimmagem
```
* Ativação do ambiente
```bash
conda activate trimmagem
```
* Instalando **trimmomatic**
```bash
conda install -c bioconda trimmomatic
```
* Instalando **trim_galore**
```bash
conda install -c bioconda trim-galore
```
* Aplicando trim_galore
  utilizado para retirar os adaptadores das sequências
```bash
  for i in {1..28}; do trim_galore --paired -o trimmagem/trimgalore dadosbrutos/genoma${i}_1* dadosbrutos/genoma${i}_2*; done
```
# Merge
*instalando **pear**
```bash
conda install -c bioconda pear
```
* Executando pear
```bash
for i in {1..28}; do pear -f dadosbrutos/genoma${i}_1* -r dadosbrutos/genoma${i}_2* -o merge/B${i} -j 20; done
```
# Qiime
* instalando qiime
  
[Tutorial](https://github.com/qiime2/qiime2)
* Executando qiime
```bash
export TMPDIR=/home/nanopore/Desktop/Artur_Barros/Amplicon/16S/Rawdata/merge
tmp=$(mktemp -d --tmpdir)
export TMPDIR=$tmp

qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path mapfile_qiime2.csv --input-format SingleEndFastqManifestPhred33V2 --output-path demux.qza

qiime demux summarize --i-data single-end-demux.qza --o-visualization demux.qzv

qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 25 --p-trunc-len 325 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza

qiime dada2 denoise-single --i-demultiplexed-seqs demux.qza --p-trim-left 20 --p-trunc-len 300 --o-representative-sequences rep-seqs-dada2.qza --o-table table-dada2.qza --o-denoising-stats stats-dada2.qza

qiime metadata tabulate --m-input-file stats-dada2.qza --o-visualization stats-dada2.qzv

mv rep-seqs-dada2.qza rep-seqs.qza

mv table-dada2.qza table.qza

qiime feature-table summarize --i-table table.qza --m-sample-metadata-file mapfile_2.txt --o-visualization table.qzv

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza

qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza

qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 98470 --m-metadata-file mapfile_2.txt --output-dir core-metrics-results

qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file mapfile_2.txt --m-metadata-column Treat --o-visualization core-metrics-results/unweighted-unifrac-Treat-significance.qzv --p-pairwise

qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 98470 --m-metadata-file mapfile_2.txt --o-visualization alpha-rarefaction.qzv

qiime feature-table rarefy --i-table table.qza --p-sampling-depth 98470 --o-rarefied-table rarefied_table.qza

qiime feature-table filter-features --i-table rarefied_table.qza --p-min-frequency 3 --o-filtered-table rarefied_filtered_table.qza

qiime feature-classifier classify-sklearn --i-classifier classifier-bac-trainning.qza --i-reads rep-seqs.qza --o-classification taxonomy132_97.qza

qiime metadata tabulate --m-input-file taxonomy132_97.qza --o-visualization taxonomy132_97.qzv

qiime taxa barplot --i-table rarefied_filtered_table.qza --i-taxonomy taxonomy132_97.qza --m-metadata-file mapfile_2.txt --o-visualization taxa-bar-rare-nondoublt.qzv

qiime tools export --input-path rarefied_filtered_table.qza --output-path exported_rarefied

qiime tools export --input-path taxonomy132_97.qza --output-path exported_rarefied
```

* Alterar na tabela taxonomy.tsv Feature ID por #OTUID,Taxon por taxonomy,Confidence por confidence
```bash
biom add-metadata -i exported_rarefied/feature-table.biom -o table-with-taxonomy.biom --observation-metadata-fp exported_rarefied/taxonomy.tsv --sc-separated taxonomy

biom convert -i table-with-taxonomy.biom -o table-with-taxonomy.txt --to-tsv --header-key taxonomy
```
# Análises Taxonômicas no R 
* Nessa etapa é onde ocorre a visualização de todos os dados processados e assim conseguir tirar conclusões sobre a questão taxônomica
### Importação dos pacotes necessários para as análises aqui relatadas
```R
library(ExpDes.pt)
library(tibble) 
library(reshape2)
library(ggplot2) 
library(microeco)
library(BiocManager)
library(phyloseq)
library(dplyr)
library(tidyverse)
library(magrittr)
library(viridis)
library(vegan)
library(phyloseq)
library(file2meco)
library(compositions)
```
### Definição de onde serão extraídos os dados
```R
setwd("C:/Users/artur/Documents/usp_projeto1/16s")
```
### Importação e modificação dos dataframes gerados das análises por linha de comando
```R
otu_table_16S <- read.delim("otu16_table.csv", row.names = 1)
taxonomy_table_16S <- read.delim("tax16_table.csv", row.names = 1)
sample_info_16S <- read_delim("mapfile16_2.txt")
sample_info_16S <- column_to_rownames(sample_info_16S, var = "sample-id")
sample_info_16S = rename(sample_info_16S, Group=Treat)
otu_table <- as.matrix(otu_table_16S)
taxonomy_table <- as.matrix(taxonomy_table_16S)
otu_table = otu_table(otu_table, taxa_are_rows = TRUE)
taxonomy_table = tax_table(taxonomy_table)
sample_info = sample_data(sample_info_16S)
phylo <- phyloseq(otu_table, taxonomy_table, sample_info)
dataset <- phyloseq2meco(phylo)
```
### Gráficos referentes a abundância 
**OBS:** Os gráficos aqui presentes são os que mais se usam nos artigos científicos, mas o pacote tem inúmeros de gráficos distintos, isso é uma escolha muito pessoal
#### Análises de abundância
**Explicação:** As análise de abundância visam verificar qual organismos estão mais presentes nos tratamentos usados e analisar o quanto essa muda de tratamento para tratamento
* Abundância relativa de filo em gráfico de barra
```R
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 20))
```
![image](https://github.com/Arturgbarros/16S/assets/125391314/961a14b0-d458-43be-afbe-1e3ce84dce9f)


*  Abundância relativa de gênero em gráfico de heatmap
```R
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
```
![image](https://github.com/Arturgbarros/16S/assets/125391314/5dc9ab44-3032-4280-ab48-53041df8034c)


#### Análises de diversidade e riqueza
**Explicação riqueza:** Análises de diversidade visam analisar qual dos tratamentos tem um número maior de diferentes, porém é uma análise onde apenas a quantidade de organismos distintos importa

**Explicação diversidade:** Análises de diversidade visam analisar o quão bem distribuida está a riqueza 
* Análise de diversidade e riqueza por boxplot com anova
```R
t1 <- trans_alpha$new(dataset = dataset, group = "Group")#cálculo da diversidade
t1$cal_diff(method = "anova") #cálculo da diversidade pelo método anova
t1$cal_diff(method = "KW") #cálculo da diversidade pelo método Kruskal-Wallis 
t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)#cálculo da diversidade pelo método Kruskal-Wallis e dunn
t1$cal_diff(method = "wilcox") #cálculo da diversidade pelo método Wilcoxon 
t1$cal_diff(method = "t.test") #cálculo da diversidade pelo método t de Student
t1$plot_alpha(measure = "Observed", y_increase = 0.3)#gráfico de riqueza
t1$plot_alpha(measure = "Chao1", y_increase = 0.3)#gráfico de riqueza
t1$plot_alpha(measure = "Shannon", y_increase = 0.1)#gráfico de diversidade
t1$plot_alpha(measure = "Simpson", y_increase = 0.3)#gráfico de diversidade
```
![image](https://github.com/Arturgbarros/16S/assets/125391314/7c7b45d8-0339-49d0-937c-530d5cb7e1c1)

![image](https://github.com/Arturgbarros/16S/assets/125391314/7fdccd52-6710-4547-ae7c-742497e2029f)

![image](https://github.com/Arturgbarros/16S/assets/125391314/9205c893-8c6c-4743-8a40-61c68cbf4a91)

![image](https://github.com/Arturgbarros/16S/assets/125391314/ebd7e66e-8f02-48c1-a23e-683de2c4864e)

* PCoa análise de intersecção de organismos presentes nos tratamentos
```R
dataset$cal_betadiv()#útil para cálculo
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
class(t1$res_ordination)
pcoa <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull"))#mais fácil de ver as diferenças
pcoa + theme_bw()
t1$cal_manova(manova_all = FALSE)
t1$res_manova #resultado estatísticos da PCoa
```
![image](https://github.com/Arturgbarros/16S/assets/125391314/16b50104-79ac-401e-be13-88b70db13acb)

* LDA score
```R
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL, p_adjust_method = "none")#caso não tenha nenhum valor p ajustado significativo
#t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL)#caso tenha valor p ajustado significativo
t1$plot_diff_bar(threshold = 1.75)
```
![image](https://github.com/Arturgbarros/16S/assets/125391314/a84d31da-a61a-42b7-9f0c-9a20e9e72149)



#### Network
**Explicação:** Network serve para avaliar 
```R
dataset <- phyloseq2meco(phylo)
###Bulk
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="Bulk")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8) 

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6) 
t1$cal_module(method = "cluster_fast_greedy") 
t1$save_network(filepath = "Bulk.gexf")

###Control
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="Control")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8) 

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy") 
t1$save_network(filepath = "Control.gexf")


###Bac_para
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="Bac_para")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8)

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy") 
t1$save_network(filepath = "Bac_para.gexf") 

###Bac_barba
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="Bac_barba")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8)

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6) 
t1$cal_module(method = "cluster_fast_greedy")
t1$save_network(filepath = "Bac_barba.gexf") 

###pae_jami
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="pae_jami")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8) 

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy") 
t1$save_network(filepath = "pae_jami.gexf") 


###pae_poly
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="pae_poly")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8) 

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6)
t1$cal_module(method = "cluster_fast_greedy") 
t1$save_network(filepath = "pae_poly.gexf")

###pool
group_second <- clone(dataset)
group_second$tax_table %<>% tidy_taxonomy
group_second$sample_table <- subset(group_second$sample_table, Group=="pool")
group_second$tidy_dataset()

set.seed(1234)
t1 <- trans_network$new(dataset = group_second, cor_method = "sparcc", 
                        use_sparcc_method = "SpiecEasi", filter_thres = 0.0005, 
                        nThreads = 8) 

t1$cal_network(COR_p_thres = 0.05,  COR_p_adjust = "none", COR_cut = 0.6) 
t1$cal_module(method = "cluster_fast_greedy") 
t1$save_network(filepath = "pool.gexf") 
```
# Extraindo função dos táxons determinados com Picrust2 e qiime2
* Ativando ambiente de instalação 
```bash
conda activate qiime2
```

* Instalação do **picrust2**
```bash
conda install -c bioconda picrust2
```
* Executando o programa
```bash
qiime picrust2 full-pipeline --i-table table.qza --i-seq rep-seqs.qza --p-threads 16 --output-dir picrust2
```
* Converção das tabelas 
```bash
qiime feature-table summarize --i-table q2-picrust2_output/pathway_abundance.qza --o-visualization q2-picrust2_output/pathway_abundance.qzv
qiime tools export --input-path q2-picrust2_output/pathway_abundance.qza --output-path pathabun_exported
biom convert -i pathabun_exported/feature-table.biom -o pathabun_exported/feature-table.biom.tsv --to-tsv
```
# Análises Funcionais no R 
* Nessa etapa é onde ocorre a visualização de todos os dados processados e assim conseguir tirar conclusões sobre a questão funcional dos táxons analisados anteriormente
### Importação dos pacotes necessários para as análises 
```R
library(ExpDes.pt)
library(tibble) 
library(reshape2)
library(ggplot2) 
library(microeco)
library(BiocManager)
library(phyloseq)
library(dplyr)
library(tidyverse)
library(magrittr)
library(viridis)
library(vegan)
library(phyloseq)
library(file2meco)
library(compositions)
```
### Seleção da pasta de trabalho
```R
setwd("C:/Users/artur/Documents/usp_projeto1/16s")
```
### Importação e tratamento dos dados funcionais 
```R
kegggenes <- read.delim("~/usp_projeto1/16s/kegggenes.tsv")
ko_table <- read.csv("~/usp_projeto1/16s/ko_table.tsv", sep="")
colnames(kegggenes)[1] = 'KO'#mudar nome
keggegenes_table = merge.data.frame(ko_table,kegggenes,all = TRUE)#juntar tabela
keggegenes_table <- na.omit(keggegenes_table)#remover valores vazios
keggegenes_table <- filter(keggegenes_table, LV1 != "09160 Human Diseases",
                           LV1 != "09150 Organismal Systems")#retirar valores desnecessários para a análise
write.table(keggegenes_table,file = 'kegg.picrust.tsv')#salvar arquivo
#fazendo OTU_table
otu_table= keggegenes_table [,-c(1:6)]#criando tabela de OTUs
otu_table$OTU_ID <- paste("OTU_", 1:nrow(otu_table), sep = "")#criando coluna OTU_ID e numerando
otu_table <- otu_table[c("OTU_ID", names(otu_table)[-which(names(otu_table) == "OTU_ID")])]#passando coluna para o inicio
rownames(otu_table) <- NULL#removendo valores nulos
otu_table <- column_to_rownames(otu_table, var = "OTU_ID")#tornanado OTU_ID linhas
otu_table <- round(otu_table)
#fazendo tax_table
tax_table = keggegenes_table [,c(1:6)]
tax_table$OTU_ID <- paste("OTU_", 1:nrow(tax_table), sep = "")#criando coluna OTU_ID e numerando
tax_table <- tax_table[c("OTU_ID", names(tax_table)[-which(names(tax_table) == "OTU_ID")])]
rownames(tax_table) <- NULL
tax_table <- column_to_rownames(tax_table, var = "OTU_ID")
#importando metadata
sample_info_16S <- read_delim("mapfile16_2.txt")
sample_info_16S <- column_to_rownames(sample_info_16S, var = "sample-id")
sample_info_16S = rename(sample_info_16S, Group=Treat)
#reeimportando dados
otu_table <- as.matrix(otu_table)
taxonomy_table <- as.matrix(tax_table)
otu_table = otu_table(otu_table, taxa_are_rows = TRUE)
taxonomy_table = tax_table(taxonomy_table)
sample_info = sample_data(sample_info_16S)
phylo <- phyloseq(otu_table, taxonomy_table, sample_info)
dataset <- phyloseq2meco(phylo)
```
### Plotagens gráficas
* Não será explicada cada uma pois as plotagens são as mesmas da taxônomia, porém nesse caso analísa-se **abundância**, **diversidade** , **riqueza** , **LDA score** e **network** das funções dos táxons encontrados.
* Porém as plotagens abaixo se tornam um pouco mais complexas devido a necessidade de mudança de legendas (rótulo e inclinação), cores, normalização do heatmap e ordenação das plotagens
```R
#Ordenando 
dataset$sample_table$Group %<>% factor(., levels = c('Bulk','Control','Bac_barba','Bac_para','Pae_jami','Pae_poly','pool'))#deve ser executado antes das plotagens

#Definindo cores fixas
color_val = c("#080707",
             "#59008a",
             "#217bc9",
             "#00a6bc",
             "#60c99f",
             "#77e723",
             "#f6e300")#deve ser citado em plotagens que permitem mudança de cor 

#Alterando nomes
names = c('Bac_barba'='T1','Bac_para'='T2','Pae_poly'='T4','Pae_jami'='T3','pool'='SynCom')#deve ser citado em cada plotagem

#virar 45°
theme(axis.text.x = element_text(angle = 45, hjust = 1))#deve ser citado em cada plotagem
```
* Plotagens
```R
t1 <- trans_abund$new(dataset = dataset, taxrank = "LV3", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)+scale_x_discrete(labels = names)

t1 <- trans_abund$new(dataset = dataset, taxrank = "LV2", ntaxa = 15, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE, color_values = turbo(15))
g1 + theme_classic() + theme(axis.title.y = element_text(size = 20))+theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(labels = names)

t1 <- trans_alpha$new(dataset = dataset, group = "Group")
head(t1$data_stat)
t1$cal_diff(method = "KW")
head(t1$res_diff)
t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
head(t1$res_diff)
t1$cal_diff(method = "wilcox")
head(t1$res_diff)
t1$cal_diff(method = "t.test")
head(t1$res_diff)
t1$cal_diff(method = "anova")
head(t1$res_diff)

#análise de diversidade e riqueza por boxplot 
t1$cal_diff(method = "anova") 
t1 = t1$plot_alpha(measure = "Observed", y_increase = 0.3, color = color_val) + theme(axis.text.x = element_text(angle = 45, hjust = 1)+scale_x_discrete(labels = names))#gráfico de riqueza
t1 + scale_x_discrete(labels=names)
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova")
t1 = t1$plot_alpha(measure = "Simpson", y_increase = 0.3, color = color_val)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)+scale_x_discrete(labels=names))#gráfico de diversidade
t1 + scale_x_discrete(labels=names)
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova") 
t1= t1$plot_alpha(measure = "InvSimpson", y_increase = 0.3, color = color_val)+ theme(axis.text.x = element_text(angle = 45, hjust = 1)+scale_x_discrete(labels=names))#gráfico de diversidade
t1 + scale_x_discrete(labels=names)
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova") 
t1=t1$plot_alpha(measure = "Shannon", y_increase = 0.1, color = color_val)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(labels=names)#gráfico de diversidade
t1 + scale_x_discrete(labels=names)
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova") 
t1=t1$plot_alpha(measure = "Chao1", add_sig_text_size = 6, color = color_val)+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+scale_x_discrete(labels=names)#gráfico de riqueza
t1 + scale_x_discrete(labels=names)
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova") 



dataset$cal_betadiv()
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
pcoa <- t1$plot_ordination(plot_color = "Group", plot_shape = "Group",color = color_val, plot_type = c("point", "chull", "centroid"))#mais fácil de ver as diferenças
pcoa + theme_bw()
t1$cal_group_distance()
t1$cal_manova(manova_all = TRUE)
t1$res_manova
t1$cal_manova(manova_all = FALSE)
t1$res_manova
```
# Network
* Nessa etapa é preparada e visualizada a network com o software **Gephi**
#### Instalação do software
```python
https://gephi.org/
```
#### **Abre** o arquivo gerado pela as análises do R
![image](https://github.com/Tutugb/16S/assets/125391314/0111fee3-dd6b-402e-96eb-2b4c810d84fb)
#### **Escolher** o tipo de distribuição dos dados e **executar** 
![image](https://github.com/Tutugb/16S/assets/125391314/1b412687-ad39-4701-a03f-32330928ad12)
#### Reduzir a largura das arestas
![image](https://github.com/Tutugb/16S/assets/125391314/a1bd7a18-fe84-44b0-91ec-f4854ffd6b3a)
#### Alterar aparência dos nós referente a abundância relativa
![image](https://github.com/Tutugb/16S/assets/125391314/3a2627ce-5a2d-4e83-bef4-f766cd37cba4)
#### Escolha da palheta dos nós referente a abundância relativa
![image](https://github.com/Tutugb/16S/assets/125391314/fccba525-ef12-4201-9abe-341b74460810)
#### OBS: caso a visualização das arestas esteja bom, não precisa aplicar filtro, pular ao passo
#### Excluindo arestas em excesso
###### No laborátório de dados organizar dados de forma crescente pela coluna **Weigth**
![image](https://github.com/Tutugb/16S/assets/125391314/58aa6b01-c665-4529-95a0-ccf4ed55528d)
###### Excluir linhas conforme um filtro pessoal estabelecido
###### Colorir linhas da coluna label classificadas como "+" de cinza e linhas classificadas como "-" de vermelho
![image](https://github.com/Tutugb/16S/assets/125391314/b76a0ea6-26b5-4d34-91a4-c7af558f3b23)
#### Exeutar parãmetros estatítisticos
![image](https://github.com/Tutugb/16S/assets/125391314/70d31ff5-51f4-4cbd-b6cc-273e0c9594be)
#### Arrumar a visualização dos dados 

**OBS:** Parâmetros pessoais

![image](https://github.com/Tutugb/16S/assets/125391314/c321ca9c-249b-4275-87a4-94bacebf6b67)
#### Salvar o network
**OBS:** Salvar com "salvar como" para assim criar outro arquivo com a indentificação da network modificada
![image](https://github.com/Tutugb/16S/assets/125391314/32c86ba9-eaa9-492f-9373-bb535f4dc4d8)




