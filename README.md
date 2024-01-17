# Workflow 16S por linha de comando
#### sumário
# 0. Baixando arquivos 
* Copiar os arquivos do sequênciamento 16S (HTTPS)
```bash
git clone https://github.com/Tutugb/16S.git
```
* Copiar os arquivos do sequênciamento 16S (SSH)
```bash
git clone git@github.com:Tutugb/16S.git
```
# 1. Criação das pastas
Deixa o workflow é uma das partes mais importantes, pois é com essa organização que não haverá confusão em quais arquivos são os certos em cada etapa por isso é essencial a criação de uma pasta para o armazenamento de cada grupo de arquivos
* Criação de pastas iniciais
```bash
mkdir 16S 16S/00.dadosbrutos 16S/01.qualidade 16S/02.trimmagem 16S/02.trimmagem/00.qualidade 16S/merge 16S/qiime
```
Os arquivos iniciais devem ser inseridos na pasta 00.dadosbrutos
# 2. Controle de qualidade
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
## 2.1. Junção dos arquivos do controle de qualidade
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

# 3. Trimmagem
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
# Análises Taxonômicas no R 
* Nessa etapa é onde ocorre a visualização de todos os dados processados e assim conseguir tirar conclusões sobre a questão taxônomica
### 1. Importação dos pacotes necessários para as análises aqui relatadas
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
### 2. Definição de onde serão extraídos os dados
```R
setwd("C:/Users/artur/Documents/usp_projeto1/16s")
```
### 3. Importação e modificação dos dataframes gerados das análises por linha de comando
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
### 4. Gráficos referentes a abundância 
**OBS:** Os gráficos aqui presentes são os que mais se usam nos artigos científicos, mas o pacote tem inúmeros de gráficos distintos, isso é uma escolha muito pessoal
#### 4.1. Análises de abundância
**Explicação:** As análise de abundância visam verificar qual organismos estão mais presentes nos tratamentos usados e analisar o quanto essa muda de tratamento para tratamento
* Abundância relativa de filo em gráfico de barra
```R
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 20))
```
![image](https://github.com/Tutugb/16S/assets/125391314/9dd3fd78-e8ed-42db-a26b-1346db9711c6)

*  Abundância relativa de gênero em gráfico de heatmap
```R
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
```
![image](https://github.com/Tutugb/16S/assets/125391314/97710733-b4c9-44b9-a211-5f5b8007bd43)

#### 4.1. Análises de diversidade e riqueza
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
![image](https://github.com/Tutugb/16S/assets/125391314/9e1da583-7756-4a23-b5d0-6bfd853b6d7e)
![image](https://github.com/Tutugb/16S/assets/125391314/1dcf9fdc-c014-4bdf-8bcb-1f17a677669d)
![image](https://github.com/Tutugb/16S/assets/125391314/4979b437-ebd2-4708-9cea-325d370235b7)
![image](https://github.com/Tutugb/16S/assets/125391314/40f93336-4c07-4676-acd7-561b6e63b170)

* PCoa análise de intersecção de organismos presentes nos tratamentos
```R
dataset$cal_betadiv()#útil para cálculo
t1 <- trans_beta$new(dataset = dataset, group = "Group", measure = "bray")
t1$cal_ordination(ordination = "PCoA")
class(t1$res_ordination)
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull"))
t1$cal_manova(manova_all = FALSE)
t1$res_manova #resultado estatísticos da PCoa
```
![image](https://github.com/Tutugb/16S/assets/125391314/72d6caa6-4f79-4dc6-b79e-5ac52f5d7c6c)

* LDA score
```R
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL, p_adjust_method = "none")#caso não tenha nenhum valor p ajustado significativo
#t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL)#caso tenha valor p ajustado significativo
t1$plot_diff_bar(threshold = 1.75)
```
![image](https://github.com/Tutugb/16S/assets/125391314/1ddeca44-36b0-4398-aae5-f03e4c083408)

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
*Ativando ambiente de instalação 
```bash
conda activate qiime2
```

*Instalação do **picrust2**
```bash
conda install -c bioconda picrust2
```
*Executando o programa
```bash
qiime picrust2 full-pipeline --i-table table.qza --i-seq rep-seqs.qza --p-threads 16 --output-dir picrust2
```
*converção das tabelas 
```bash
qiime feature-table summarize --i-table q2-picrust2_output/pathway_abundance.qza --o-visualization q2-picrust2_output/pathway_abundance.qzv
qiime tools export --input-path q2-picrust2_output/pathway_abundance.qza --output-path pathabun_exported
biom convert -i pathabun_exported/feature-table.biom -o pathabun_exported/feature-table.biom.tsv --to-tsv
```
# Análises Funcionais no R 
