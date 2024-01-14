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
