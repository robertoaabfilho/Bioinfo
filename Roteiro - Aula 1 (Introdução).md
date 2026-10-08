# Material de Apoio — Aula 1
## Introdução ao Ambiente Computacional para Bioinformática

**Disciplina:** Introdução à Bioinformática — PPGED/IME
**Carga horária:** 4 horas (síncrona) + atividade assíncrona
**Para quem é este material:** alunos de pós-graduação em Engenharia de Defesa, **sem** experiência prévia em programação ou em bioinformática.

---

## Como usar este material

Este documento foi escrito para você consultar **antes, durante e depois** da aula:

- **Antes:** leia a Parte 1 (teoria) com calma. Você não precisa decorar nada — a ideia é chegar à aula já conhecendo as palavras.
- **Durante:** use a Parte 2 (prática) como roteiro. Cada passo tem o comando exato e o resultado esperado.
- **Depois:** use os exercícios, o checklist, a seção de problemas comuns e o glossário.

> **Dica de ouro:** se algo der erro, **não se desespere**. Erro é parte normal do trabalho com computação. A seção "Problemas comuns" (no fim) cobre os erros que mais aparecem, e há uma orientação de como pedir ajuda a uma IA (Claude, Gemini etc.).

### Mapa da aula

| Bloco | Tema | Tempo |
|---|---|---|
| 1 | Por que bioinformática precisa de computador | ~40 min |
| 2 | Terminal e linha de comando | ~1 h |
| 3 | Ambiente de trabalho: Conda, Python, R e Jupyter | ~50 min |
| 4 | Python e R aplicados à bioinformática | ~40 min |
| 5 | Repositórios públicos (NCBI, ENA, Ensembl) | ~30 min |
| 6 | Automatizando buscas: APIs do NCBI e do PubChem | ~30 min |

---

# PARTE 1 — TEORIA

## 1. Por que bioinformática precisa de computador?

### 1.1 O que é bioinformática

**Bioinformática** é a área que usa **computação, estatística e matemática** para armazenar, organizar e interpretar dados biológicos. Uma forma útil de pensar: é *engenharia de dados aplicada à biologia*.

Se você vem da Engenharia de Defesa, já conhece boa parte da lógica: processamento de sinais, análise de grandes volumes de dados, automação de rotinas. A novidade é o **tipo de dado** — sequências de DNA, RNA, proteínas, estruturas moleculares.

### 1.2 O escalão do problema

Um único genoma humano tem cerca de **3 bilhões de "letras"** (bases). Um sequenciador moderno pode gerar **gigabytes de dados por dia**. Nenhum ser humano consegue ler isso manualmente — por isso o computador não é opcional, é a ferramenta central.

### 1.3 Aplicações que interessam à Defesa

- **Identificação de patógenos:** descobrir rapidamente qual microrganismo está presente em uma amostra comparando sua sequência com bancos de dados.
- **Biodefesa:** monitorar toxinas, agentes biológicos e sua evolução.
- **Análise forense:** comparar amostras biológicas.
- **Química de proteínas e macromoléculas:** projetar ou avaliar moléculas (por exemplo, inibidores, antídotos, biossensores).

### 1.4 Conceitos biológicos mínimos

Você não precisa ser biólogo, mas precisa destas palavras:

| Termo | O que é, em linguagem simples |
|---|---|
| **DNA** | Molécula que guarda a "receita" do organismo. Escrita com 4 letras: **A, T, C, G** (bases). |
| **RNA** | Cópia de trabalho do DNA. Troca a letra T por **U**. |
| **Proteína** | Molécula que executa funções na célula. É uma cadeia de **aminoácidos** (20 tipos, cada um representado por uma letra, ex.: M, A, K). |
| **Gene** | Trecho do DNA que contém a instrução para fazer uma proteína (ou RNA funcional). |
| **Genoma** | Todo o DNA de um organismo. |
| **Sequência** | A sequência de letras de um DNA, RNA ou proteína — é o dado básico da bioinformática. |
| **Dogma central** | O fluxo da informação: **DNA → RNA → proteína**. |

### 1.5 Como uma sequência é guardada: o formato FASTA

O formato mais simples e universal para sequências é o **FASTA**. Um arquivo FASTA tem:

1. Uma linha de cabeçalho que começa com `>` (identificação da sequência).
2. Uma ou mais linhas com a sequência em si.

Exemplo (inventado, apenas para ilustrar):

```
>exemplo_gene_01 Descrição da sequência
ATGGCGTACGTTAGCTAGCTAGGCTAACGTTAGC
TAGCTAGCTAGGATCGATCGTACGATCGATCGAT
```

Você vai ver arquivos FASTA reais na prática (Passo 7).

---

## 2. Terminal e linha de comando

### 2.1 O que é

O **terminal** (ou *shell*) é uma janela onde você **digita comandos** para o computador, em vez de clicar em ícones. Parece antiquado, mas é **mais rápido, mais preciso e automatizável** — e praticamente toda ferramenta séria de bioinformática (BLAST, SAMtools, GROMACS etc.) funciona assim.

### 2.2 Os cinco comandos que você precisa hoje

| Comando | Para que serve | Exemplo |
|---|---|---|
| `pwd` | Mostra **em qual pasta** você está | `pwd` |
| `ls` | **Lista** os arquivos da pasta atual | `ls` |
| `cd` | **Entra** em uma pasta | `cd bioinfo_aula1` |
| `mkdir` | **Cria** uma pasta | `mkdir bioinfo_aula1` |
| `cp` / `mv` | **Copia** / **move (ou renomeia)** arquivos | `cp a.txt b.txt` |

Dois comandos extras úteis para "olhar" arquivos (Linux/macOS/WSL):

| Comando | Para que serve |
|---|---|
| `head arquivo` | Mostra as **primeiras linhas** de um arquivo |
| `wc -l arquivo` | Conta quantas **linhas** o arquivo tem |

> **Objetivo de hoje:** perder o medo do terminal, não virar especialista em Linux.

> **Atenção (Windows sem WSL):** o Prompt de Comando/PowerShell do Windows tem comandos parecidos, mas não idênticos (por exemplo, `dir` em vez de `ls`, `type` em vez de `head`). Para seguir o curso de forma consistente, recomendamos usar o **Anaconda Prompt** ou o **WSL** (explicado a seguir).

---

## 3. Seu ambiente de trabalho

Esta é a parte mais importante da aula: um ambiente bem montado evita 90% dos problemas futuros.

### 3.1 O problema que o Conda resolve

Bioinformática usa **dezenas de programas**, e cada um depende de **versões específicas** de outras bibliotecas. Instalar tudo "junto" quase sempre gera conflito: uma ferramenta atualiza uma biblioteca e outra para de funcionar.

O **Conda** é um **gerenciador de pacotes e de ambientes**. Pense assim:

> Um *ambiente* Conda é uma **gaveta**. Dentro dela ficam as ferramentas de um projeto. Outra gaveta pode ter versões diferentes, sem conflito.

### 3.2 Anaconda ou Miniconda?

| | **Anaconda** | **Miniconda** |
|---|---|---|
| O que é | Conda + Python + centenas de pacotes pré-instalados (incluindo Jupyter, Pandas, NumPy) | Apenas Conda + Python básico |
| Tamanho | Grande (vários GB) | Pequeno |
| Jupyter já vem? | **Sim** | **Não** — precisa instalar |
| Quando usar | macOS/Linux/Windows sem WSL, se quiser tudo pronto | Dentro do WSL ou quando se quer um ambiente enxuto |

### 3.3 Python

Linguagem de programação de uso geral, muito popular em ciência de dados. Em bioinformática, é usada para **manipular sequências, acessar bancos de dados, automatizar análises**. É relativamente fácil de ler — por isso é boa para iniciantes.

### 3.4 R

Linguagem criada para **estatística e gráficos**. Em bioinformática é muito usada em **expressão gênica** (por meio do **Bioconductor**, que veremos adiante). O R **não vem** no Anaconda nem no Miniconda — é instalado à parte.

### 3.5 Jupyter Notebook

Um **caderno de laboratório digital**: você mistura **texto, código e resultados** (tabelas, gráficos) no mesmo documento, em blocos chamados **células**. Executa-se uma célula por vez (`Shift + Enter`) e vê-se o resultado logo abaixo. É ideal para explorar dados e documentar o que foi feito.

### 3.6 Google Colab (plano B)

É um Jupyter que roda **no navegador**, na nuvem do Google, sem instalar nada. Se seu computador estiver dando muito trabalho, use-o para não ficar para trás. Link: https://colab.research.google.com

### 3.7 VS Code

O **Visual Studio Code** é um **editor de código** gratuito da Microsoft. Nele você cria arquivos, escreve código, abre o terminal integrado e executa scripts. Vamos usá-lo para rodar os scripts de API.

### 3.8 WSL — para quem usa Windows

Muitas ferramentas do curso (BWA, SAMtools, BLAST em linha de comando, GROMACS) **só rodam em Linux**. O **WSL (Windows Subsystem for Linux)** instala um "Linux dentro do Windows", sem precisar de máquina virtual nem de dual boot.

> **Importante:** o WSL é necessário para as **ferramentas de Linux das próximas unidades**. Para **esta aula** — incluindo os scripts de API (Passos 8 e 9) — você pode trabalhar **direto no Windows com o VS Code**, sem WSL.

---

## 4. Python e R aplicados à bioinformática

Não vamos ensinar a programar hoje; vamos apenas apresentar as ferramentas que reaparecerão no curso.

### 4.1 No Python

| Biblioteca | O que faz | Analogia |
|---|---|---|
| **Biopython** | Lê, escreve e manipula sequências biológicas; acessa o NCBI | O "canivete suíço" da bioinformática em Python |
| **NumPy** | Cálculos numéricos eficientes com vetores e matrizes | Calculadora científica turbinada |
| **Pandas** | Tabelas de dados (linhas e colunas), filtros, estatísticas | "Excel programável" e muito mais rápido |
| **Requests** | Faz pedidos a sites e APIs pela internet | Seu "navegador" dentro do código |

### 4.2 No R

- **R base / r-essentials:** o núcleo da linguagem e os pacotes mais usados em análise de dados.
- **Bioconductor:** o grande repositório de pacotes de bioinformática para R. Será essencial na **Unidade 4 (Transcriptômica)**, com ferramentas como DESeq2 e edgeR.

### 4.3 Python *versus* R

Não competem — **se complementam**. Python tende a ser mais forte em automação e manipulação geral de dados; R, em estatística e visualização de dados de expressão gênica. Em um projeto real, é comum usar os dois.

---

## 5. Onde os dados biológicos moram

### 5.1 NCBI

**National Center for Biotechnology Information** (EUA). É o maior repositório público de dados biológicos. Dentro dele:

- **GenBank:** banco de sequências de DNA/RNA submetidas por pesquisadores do mundo todo.
- **RefSeq:** coleção de sequências **curadas** (revisadas), de referência.
- **PubMed:** banco de artigos científicos da área biomédica.
- **Protein, Gene, Taxonomy, SRA** e outros.

Cada registro tem um **número de acesso** (*accession number*), por exemplo `NC_001988.2` — é o "CPF" daquela sequência, único e estável. Você vai guardar um hoje.

### 5.2 ENA

**European Nucleotide Archive** (Europa). É o equivalente europeu, forte em **dados brutos de sequenciamento**. Os dados são sincronizados entre NCBI, ENA e o banco japonês (DDBJ).

### 5.3 Ensembl

Banco focado em **genomas anotados**, ou seja, genomas que já vêm com a informação de "o que cada trecho do DNA faz" (genes, regiões regulatórias etc.).

### 5.4 PubChem

Banco público de **moléculas pequenas** (compostos químicos, fármacos, ligantes). Cada composto tem um **CID** (*Compound ID*). Será muito útil na **Unidade 7 (docking molecular)**, quando usaremos moléculas como ligantes.

---

## 6. APIs — automatizando a busca

### 6.1 O que é uma API

**API** (*Application Programming Interface*) é uma forma **padronizada** de um programa pedir informações a outro programa.

> **Analogia do garçom:** você não entra na cozinha para pegar a comida. Você faz um pedido pelo **cardápio** (a API); o garçom leva à cozinha (o banco de dados) e devolve **exatamente o que você pediu**.

Diferença para a busca manual:

| Busca manual (site) | API (código) |
|---|---|
| Você clica, digita, baixa | O script faz tudo sozinho |
| Bom para 1 sequência | Bom para 1 ou 10.000 |
| Difícil de repetir exatamente | **Reprodutível**: o mesmo código dá o mesmo resultado |

### 6.2 Como funciona por baixo (sem mistério)

1. Seu código envia uma **requisição** (um endereço web com parâmetros) ao servidor.
2. O servidor responde com dados, geralmente em **JSON** (um formato de texto organizado em "chave: valor") ou em **FASTA**.
3. Seu código lê a resposta e usa o que precisa.

Um **código de status** acompanha a resposta: `200` = deu certo; `404` = não encontrado; `503` = servidor ocupado/indisponível.

### 6.3 API do NCBI: Entrez + Biopython

- **Entrez** é o sistema de busca e recuperação de dados do NCBI.
- Você poderia montar as URLs "na mão", mas o módulo **`Bio.Entrez`** do Biopython faz isso por você.
- O NCBI exige que você se **identifique com um e-mail** (`Entrez.email`) para poder contatá-lo em caso de uso abusivo.
- Uma **API key** gratuita aumenta o limite de requisições (de 3 para 10 por segundo). É opcional, mas recomendada se você fizer muitas buscas.

### 6.4 API do PubChem: PUG-REST

- **PUG-REST** é a interface do PubChem baseada em **URLs**.
- Você monta um endereço no padrão `.../compound/name/NOME/property/PROPRIEDADES/JSON` e recebe os dados.
- Exemplo de propriedades: fórmula molecular, peso molecular, SMILES (um texto que descreve a estrutura química).

> **Boas práticas:** não envie centenas de requisições seguidas sem pausa; os servidores são públicos e compartilhados. Por isso nosso código do PubChem inclui **tentativas com espera** (*retry com backoff*).

---

## 7. Reprodutibilidade — a ideia que une tudo

Hoje você fará a mesma tarefa de duas formas: **manual** (site) e **automatizada** (API). O ganho central não é só velocidade — é **reprodutibilidade**: qualquer pessoa, rodando o seu código, obtém o mesmo resultado. Isso é o que transforma uma consulta pontual em um **pipeline de análise**, conceito que reaparece em todas as unidades do curso.

---

# PARTE 2 — PRÁTICA

> Siga os passos **na ordem**. Se travar, anote em qual passo parou e consulte "Problemas comuns" ao final. Todos os comandos devem ser digitados **exatamente** como aparecem (copie e cole quando possível).

## Visão geral do caminho

| Passo | O que você faz |
|---|---|
| 1 | Instala o Anaconda (ou usa o Colab) |
| 1.1 | *(Windows, para as próximas unidades)* WSL + Miniconda via VS Code |
| 2 | Confirma o Python e cria o ambiente `bioinfo` |
| 3 | Instala o R |
| 4 | Instala as bibliotecas |
| 5 | Primeiros comandos no terminal |
| 6 | Abre o Jupyter Notebook |
| 7 | Busca manual no NCBI |
| 8 | Busca via API do NCBI (script no VS Code) |
| 9 | Busca via API do PubChem (script no VS Code) |

---

## Passo 1 — Instale o Anaconda

1. Acesse https://www.anaconda.com/download
2. Baixe a versão do seu sistema (Windows, macOS ou Linux).
3. Execute o instalador com as opções padrão.
4. Abra o **Anaconda Prompt** (Windows) ou o **Terminal** (macOS/Linux) e confirme:

```bash
conda --version
```

**Resultado esperado:** um número de versão, como `conda 24.x.x`.

> **Plano B sem instalação:** use o Google Colab (https://colab.research.google.com). Você só precisa de uma conta Google.

---

## Passo 1.1 — Windows: WSL + Miniconda pelo VS Code

> **Quando fazer:** o WSL é para as **ferramentas Linux das unidades seguintes**. Se você está com pouco tempo, pode **pular este passo agora** e fazê-lo depois da aula; os Passos 8 e 9 funcionam sem ele. Usuários de macOS/Linux não precisam deste passo.

### 1.1.1 Instale o VS Code
Baixe em https://code.visualstudio.com/ e instale com as opções padrão.

### 1.1.2 Instale a extensão WSL
No VS Code: ícone de **Extensões** (barra lateral) → pesquise **WSL** (Microsoft) → **Install**. Instale também **Remote Development**, se não vier junto.

### 1.1.3 Instale o WSL com Ubuntu
No VS Code: **Terminal → New Terminal** (confirme que é **PowerShell**) e rode:

```powershell
wsl --install
```

Reinicie o computador quando solicitado. Depois, o Ubuntu abrirá e pedirá para você criar **usuário e senha** (a senha não aparece enquanto você digita — é normal).

> **Erro de virtualização?** É preciso habilitar **Intel VT-x** ou **AMD-V** na BIOS. Procure "habilitar virtualização BIOS" + modelo do seu notebook, ou leve a dúvida à aula.

### 1.1.4 Conecte o VS Code ao Ubuntu
No canto inferior esquerdo, clique no ícone azul `><` → **Connect to WSL**. Quando aparecer **WSL: Ubuntu** no canto inferior esquerdo, o terminal integrado (`Ctrl + '`) passa a ser o **bash do Ubuntu**.

### 1.1.5 Instale o Miniconda dentro do Ubuntu
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
Pressione `Enter` para avançar nos termos, digite `yes` para aceitar a licença e `yes` para inicializar o Miniconda. **Feche e reabra o terminal.**

### 1.1.6 Confirme
```bash
conda --version
```

> **Lembrete:** Miniconda **não** traz Jupyter nem bibliotecas — instalaremos no Passo 4.

---

## Passo 2 — Confirme o Python e crie o ambiente `bioinfo`

### 2.1 Verifique o Python
```bash
python --version
```
**Resultado esperado:** algo como `Python 3.11.x` (ou 3.12).

### 2.2 Crie o ambiente dedicado ao curso
```bash
conda create -n bioinfo python=3.11 -y
conda activate bioinfo
```

> **Regra prática:** sempre que abrir um novo terminal para trabalhar no curso, rode `conda activate bioinfo` antes. O nome do ambiente aparece entre parênteses no início da linha: `(bioinfo)`.

---

## Passo 3 — Instale o R

Com o ambiente `bioinfo` ativo:

```bash
conda install -c conda-forge r-base r-essentials -y
```

Confirme:
```bash
R --version
```

**O que cada pacote faz:** `r-base` é o núcleo da linguagem; `r-essentials` traz os pacotes mais usados de análise de dados.

> **Alternativa:** instalar pelo site oficial (https://cran.r-project.org/) e, opcionalmente, o **RStudio** (https://posit.co/download/rstudio-desktop/).

**Para desinstalar o R depois:**
```bash
conda activate bioinfo
conda remove r-base r-essentials
```
Para apagar o ambiente inteiro e recomeçar:
```bash
conda deactivate
conda env remove -n bioinfo
```

---

## Passo 4 — Instale as bibliotecas

### 4.1 Python
Com `bioinfo` ativo:

```bash
python -m pip install jupyter biopython requests pandas numpy
```

> **Por que `python -m pip` e não só `pip`?** Se o computador tem mais de um Python instalado, `pip` sozinho pode instalar no Python "errado". `python -m pip` garante que a instalação vai para o **mesmo Python** que o comando `python` usa.

Teste:
```bash
python -c "import Bio, requests, pandas, numpy; print('Bibliotecas OK')"
```

### 4.2 R (Bioconductor)
Digite `R` no terminal para abrir o console do R e execute:

```r
install.packages("BiocManager")
BiocManager::install()
```
Pode demorar vários minutos na primeira vez. Para sair: `q()` e responda `n`.

---

## Passo 5 — Primeiros comandos no terminal

Abra o terminal (Anaconda Prompt, terminal do VS Code, ou Terminal) com o ambiente `bioinfo` ativo e rode, um por vez:

```bash
pwd
ls
mkdir bioinfo_aula1
cd bioinfo_aula1
pwd
```

**Resultado esperado:** o último `pwd` mostra um caminho terminando em `bioinfo_aula1`.

> **Windows sem WSL:** `pwd` e `mkdir` funcionam no Anaconda Prompt/PowerShell; no lugar de `ls`, use `dir`.

---

## Passo 6 — Abra o Jupyter Notebook

Dentro da pasta `bioinfo_aula1`, com `bioinfo` ativo:

```bash
jupyter notebook
```

Uma aba abrirá no navegador. Clique em **New → Python 3 (ipykernel)**. Em uma célula, digite:

```python
print("Olá, bioinformática!")
```
e pressione `Shift + Enter`.

**Resultado esperado:** a frase aparece abaixo da célula.

Para encerrar o Jupyter: volte ao terminal e pressione `Ctrl + C` (duas vezes, se pedir confirmação).

---

## Passo 7 — Busca manual no NCBI (sem código)

1. Acesse https://www.ncbi.nlm.nih.gov/
2. Na caixa de busca, escolha a base **Nucleotide**.
3. Pesquise por algo de seu interesse. Sugestão ligada à biodefesa: `Bacillus anthracis toxin`.
4. Abra um resultado e **anote o número de acesso** (*accession*), por exemplo `NC_001988.2`.
5. Clique em **Send to → File → Format: FASTA → Create File** para baixar.
6. Abra o arquivo baixado em um editor de texto simples e observe: linha de cabeçalho com `>` e depois a sequência.

> **Anote o accession number.** Você o usará no Passo 8.

**Opcional — olhando o arquivo pelo terminal** (Linux/macOS/WSL): coloque o arquivo na pasta `bioinfo_aula1` e rode
```bash
head sequence.fasta
wc -l sequence.fasta
```
(substitua `sequence.fasta` pelo nome do seu arquivo).

---

## Passo 8 — A mesma busca, agora via API do NCBI (script no VS Code)

> Aqui você trabalha **direto no VS Code, no Windows ou onde estiver — não é preciso WSL**.

**1. Abra o VS Code** e crie um arquivo novo: **File → New Text File**.

**2. Salve** como `busca_ncbi.py` (`Ctrl + S`), em uma pasta de sua escolha (ex.: `Documentos\Bioinfo`).

**3. Garanta que o Biopython está instalado no Python certo.** Abra o terminal integrado (**Terminal → New Terminal**) e rode:

```bash
python -m pip install biopython requests
```

**4. Cole o código abaixo** em `busca_ncbi.py`, trocando o e-mail e o *accession* pelos seus:

```python
from Bio import Entrez, SeqIO

# Identificação exigida pelo NCBI (pode ser seu e-mail pessoal ou institucional)
Entrez.email = "seu_email@exemplo.com"

# Troque pelo número de acesso que você anotou no Passo 7
accession = "NC_001988.2"

# Faz o pedido ao NCBI (o "garçom" vai até o banco de dados)
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")

# Lê o resultado como um registro de sequência
registro = SeqIO.read(handle, "fasta")
handle.close()

# Mostra as informações básicas
print("ID:", registro.id)
print("Descrição:", registro.description)
print("Tamanho da sequência:", len(registro.seq), "bases")
print("Primeiras 100 bases:", registro.seq[:100])
```

**5. Salve** (`Ctrl + S`) e **execute** no terminal integrado:

```bash
python busca_ncbi.py
```

**Resultado esperado:** ID, descrição, tamanho e as primeiras 100 bases da sequência — o mesmo que você baixou manualmente, agora em segundos.

### Entendendo o código, linha por linha

| Trecho | O que faz |
|---|---|
| `from Bio import Entrez, SeqIO` | Traz do Biopython duas "ferramentas": `Entrez` (acessa o NCBI) e `SeqIO` (lê sequências) |
| `Entrez.email = ...` | Identifica você para o NCBI |
| `accession = ...` | Guarda o número de acesso em uma variável |
| `Entrez.efetch(...)` | Pede ao NCBI o registro (`db`=base, `id`=acesso, `rettype`=formato) |
| `SeqIO.read(handle, "fasta")` | Interpreta a resposta como FASTA |
| `print(...)` | Mostra informações na tela |

> **API key (opcional):** se for fazer muitas buscas, crie uma chave gratuita em https://www.ncbi.nlm.nih.gov/account/ e adicione `Entrez.api_key = "sua_chave"` logo abaixo da linha do e-mail.

---

## Passo 9 — Busque um composto no PubChem (script no VS Code)

**1.** Crie outro arquivo (**File → New Text File**) e salve como `busca_pubchem.py`.

**2.** Cole o código abaixo (você pode trocar o composto, ex.: `"penicillin"`, `"caffeine"`; use o **nome em inglês**):

```python
import requests
import time

# Nome do composto (em inglês)
composto = "ciprofloxacin"

# Monta o endereço da API PUG-REST do PubChem
url = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    f"{composto}/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"
)

tentativas = 5
espera = 10  # segundos; dobra a cada nova tentativa

for tentativa in range(1, tentativas + 1):
    resposta = requests.get(url)
    dados = resposta.json()

    if "PropertyTable" in dados:
        propriedades = dados["PropertyTable"]["Properties"][0]
        print("CID (Compound ID):", propriedades["CID"])
        print("Fórmula molecular:", propriedades["MolecularFormula"])
        print("Peso molecular:", propriedades["MolecularWeight"])
        print("SMILES:", propriedades.get("ConnectivitySMILES", "não disponível"))
        break
    else:
        print(f"Tentativa {tentativa} falhou (status {resposta.status_code}):",
              dados.get("Fault", dados))
        if tentativa < tentativas:
            print(f"Aguardando {espera}s antes de tentar novamente...")
            time.sleep(espera)
            espera *= 2
```

**3.** Salve e execute:

```bash
python busca_pubchem.py
```

**Resultado esperado (para ciprofloxacin):**
```
CID (Compound ID): 2764
Fórmula molecular: C17H18FN3O3
Peso molecular: 331.34
SMILES: ...
```

### Entendendo o código

| Trecho | O que faz |
|---|---|
| `requests.get(url)` | Envia o pedido à API |
| `.json()` | Converte a resposta de texto JSON em um dicionário Python |
| `if "PropertyTable" in dados` | Verifica se a resposta veio no formato esperado (sem erro) |
| `for ... range(...)` + `time.sleep` | Tenta de novo, esperando mais a cada vez, se o servidor estiver ocupado |
| `.get("ConnectivitySMILES", ...)` | Pega o SMILES; se a chave não existir, devolve "não disponível" em vez de dar erro |

> **Nota sobre o SMILES:** o PubChem renomeou internamente essa propriedade; ao pedir `CanonicalSMILES` na URL, a resposta traz o dado sob a chave `ConnectivitySMILES`. Por isso o código usa essa chave.

> **Anote o CID.** Ele será reutilizado quando formos ao docking molecular.

---

# EXERCÍCIOS PARA FIXAR

1. **Troque o alvo (NCBI).** Escolha outro gene ou proteína de interesse na sua pesquisa, anote o accession e rode o `busca_ncbi.py` novamente. Qual o tamanho da sequência?
2. **Troque o composto (PubChem).** Busque outros dois compostos (por exemplo, um antibiótico e um solvente). Compare fórmulas e pesos moleculares.
3. **Banco de proteínas.** No script do NCBI, troque `db="nucleotide"` por `db="protein"` e use o accession de uma proteína. O que muda na saída?
4. **Para quem quiser ir além.** Rode o script do NCBI para **três** accessions diferentes, usando uma lista e um `for`. Dica: coloque os accessions em uma lista `["NC_001988.2", "...", "..."]` e repita o bloco de busca dentro do `for`.
5. **Reflexão (3–5 linhas).** Em que situação do *seu* trabalho de pesquisa a busca automatizada via API traria mais ganho do que a busca manual?

---

# CHECKLIST FINAL

- [ ] Anaconda (ou Colab) funcionando — `conda --version` retorna uma versão
- [ ] *(Windows, opcional agora)* VS Code + WSL + Miniconda instalados
- [ ] Ambiente `bioinfo` criado e ativado
- [ ] Python confirmado (`python --version`)
- [ ] R instalado (`R --version`)
- [ ] Bibliotecas instaladas (`Bibliotecas OK` no teste do Passo 4.1)
- [ ] Bioconductor instalado no R
- [ ] Naveguei entre pastas pelo terminal
- [ ] Abri o Jupyter e rodei uma célula
- [ ] Encontrei e baixei uma sequência no site do NCBI
- [ ] Repeti a busca com `busca_ncbi.py`
- [ ] Busquei um composto com `busca_pubchem.py`
- [ ] Anotei o *accession number* e o **CID**

---

# PROBLEMAS COMUNS

**`conda` não é reconhecido.**
Windows: abra o **Anaconda Prompt**, não o Prompt de Comando comum. Linux/WSL: feche e reabra o terminal após instalar o Miniconda.

**`pip install` deu erro de permissão.**
Use `python -m pip install --user nome_do_pacote` ou confirme que o ambiente `bioinfo` está ativo.

**`ModuleNotFoundError: No module named 'Bio'` ao rodar o Passo 8.**
O terminal está usando um Python **diferente** daquele onde o Biopython foi instalado. Solução:
```bash
python -m pip install biopython requests
```
no **mesmo terminal** em que você roda `python busca_ncbi.py`. Se persistir, confira com `where python` (Windows) ou `which python` (Linux/macOS/WSL) qual Python está em uso e com `python -m pip show biopython` se o pacote está instalado nele.

**`KeyError: 'CanonicalSMILES'` no Passo 9.**
Use `propriedades.get("ConnectivitySMILES", "não disponível")`, como no código acima. O PubChem renomeou a chave na resposta.

**Erro `503` / `PUGREST.ServerBusy` no PubChem.**
O servidor está sobrecarregado (comum quando muitos alunos consultam ao mesmo tempo). O script já tenta novamente com espera crescente. Se persistir, teste a URL direto no navegador: se o erro aparecer lá, é instabilidade do PubChem. Espere alguns minutos ou continue o roteiro e volte depois.

**Erro `HTTPError` / composto não encontrado.**
Use o nome **em inglês** e confira a grafia.

**O Jupyter não abre no navegador.**
Copie o endereço que aparece no terminal (algo como `http://localhost:8888/...`) e cole manualmente no navegador.

**`R` não é reconhecido.**
Ative o ambiente em que instalou o R: `conda activate bioinfo`.

**`wsl --install` falha por virtualização.**
Habilite **Intel VT-x** ou **AMD-V** na BIOS.

**Tenho mais de uma versão de Python — como saber qual está em uso?**
```bash
python --version
where python     # Windows (PowerShell)
which python     # Linux/macOS/WSL
conda env list   # lista os ambientes; o ativo tem *
```
Para trocar de ambiente: `conda activate nome_do_ambiente`. Para voltar ao padrão: `conda deactivate`. Para criar um ambiente com outra versão: `conda create -n bioinfo python=3.10 -y`.

**Qual e-mail uso em `Entrez.email`?**
Pode ser o institucional (IME) ou pessoal. Serve só para o NCBI contatá-lo em caso de uso excessivo.

### Como usar uma IA (Claude, Gemini etc.) para resolver erros

Usar IA para depurar é prática comum e legítima. Para obter boas respostas:

1. **Copie a mensagem de erro completa**, não um resumo.
2. **Informe seu ambiente:** sistema operacional (Windows/WSL/macOS/Linux), se usa Anaconda ou Miniconda, e o ambiente Conda ativo.
3. **Cole o comando ou o trecho de código** que você rodou antes do erro.
4. **Entenda antes de executar** comandos que apagam coisas (`rm`, `conda remove`, `conda env remove`).
5. Se persistir, **leve o print e o texto do erro para a aula** — cada máquina tem particularidades (permissões, antivírus, versão do Windows).

---

# GLOSSÁRIO

| Termo | Significado |
|---|---|
| **Accession number** | Identificador único de um registro em um banco de dados (ex.: `NC_001988.2`) |
| **API** | Interface padronizada para um programa pedir dados a outro |
| **Base (nucleotídeo)** | "Letra" do DNA/RNA: A, T (U no RNA), C, G |
| **Biopython** | Biblioteca Python para bioinformática |
| **Bioconductor** | Repositório de pacotes de bioinformática para R |
| **CID** | *Compound ID* — identificador de um composto no PubChem |
| **Conda** | Gerenciador de pacotes e ambientes |
| **Célula (Jupyter)** | Bloco de código ou texto de um notebook |
| **Entrez** | Sistema de busca e recuperação de dados do NCBI |
| **FASTA** | Formato de texto para sequências (cabeçalho com `>` + sequência) |
| **Genoma** | Todo o DNA de um organismo |
| **JSON** | Formato de texto organizado em pares chave–valor |
| **Jupyter** | Ambiente de notebooks (texto + código + resultados) |
| **Kernel** | Motor que executa o código de um notebook |
| **Pipeline** | Sequência automatizada e reprodutível de etapas de análise |
| **PUG-REST** | API do PubChem baseada em URLs |
| **Requisição** | Pedido enviado a um servidor |
| **Shell / Terminal** | Interface de linha de comando |
| **SMILES** | Texto que descreve a estrutura de uma molécula |
| **WSL** | *Windows Subsystem for Linux* — Linux dentro do Windows |

---

# REFERÊNCIAS E LEITURAS

- Lesk, A. M. *Introduction to Bioinformatics*. Oxford University Press.
- Pevsner, J. *Bioinformatics and Functional Genomics*. Wiley-Blackwell.
- Mariano, D. *Python para Bioinformática*. Novatec, 2025.
- Documentação oficial: Biopython (https://biopython.org), NCBI Entrez Programming Utilities (https://www.ncbi.nlm.nih.gov/books/NBK25501/), PubChem PUG-REST (https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest), Conda (https://docs.conda.io), VS Code (https://code.visualstudio.com/docs).

---

*Traga o accession number e o CID que você encontrou para a próxima aula — vamos retomá-los nas unidades seguintes.*
