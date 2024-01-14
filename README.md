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
# 1. Controle de qualidade
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
* Executando programa em um arquivo
```bash
fastqc genoma1_1.fq
```
* Executando programa com **loop**
```bash
for i in {1..28}; do fastqc genoma${i}_1* genoma${i}_2*; done
```
Após a obtenção de todos os arquivos .zip e .html é possível abrir cada um dos 28 arquivos .html e avaliar a qualidade de cada um, porém como isso é um trabalho muito exaustivo e repetitivo, usa-se a ferramenta abaixo para otimmização
## 1.1 
