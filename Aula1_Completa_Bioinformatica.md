# Aula 1 – Introdução ao Ambiente Computacional para Bioinformática

**Disciplina:** Introdução à Bioinformática – PPGED/IME
**Carga horária:** 4 horas

Bem-vindo(a) à primeira aula da disciplina! Você não precisa ter nenhuma experiência prévia com programação ou bioinformática — vamos construir esse conhecimento juntos, passo a passo. Ao final desta aula, você terá o ambiente de trabalho pronto e terá feito sua primeira busca de dados biológicos, tanto manualmente quanto de forma automatizada.

Este documento está dividido em duas partes:
- **Parte A** — o conteúdo teórico da aula, explicando o "porquê" de cada ferramenta.
- **Parte B** — o roteiro prático, com todos os comandos que você vai executar no seu computador.

---

## PARTE A — O QUE VOCÊ VAI APRENDER NESTA AULA

### 1. Por que bioinformática precisa de computador?

Um único genoma humano tem cerca de 3 bilhões de "letras" (bases A, T, C, G) — impossível de analisar manualmente. Pense nisso como **big data biológico**: cada sequenciador moderno gera gigabytes de dados por dia, e bioinformática é, na prática, engenharia de dados aplicada à biologia.

Como você é aluno de Engenharia de Defesa, vale conectar isso a aplicações que já fazem sentido para você: identificação de patógenos, biodefesa, análise forense de amostras biológicas. A lógica computacional que você já domina (processamento de sinais, análise de dados) se aplica diretamente aqui.

### 2. Terminal e linha de comando

Esse costuma ser o ponto que mais assusta quem está começando — mas você vai perceber que é mais simples do que parece.

- **O que é o Terminal/Shell:** é uma forma de "conversar" com o computador digitando comandos, em vez de clicar em ícones. É mais rápido, e é como praticamente toda ferramenta séria de bioinformática funciona.
- Você vai aprender só o essencial: navegar entre pastas (`cd`), listar arquivos (`ls`), criar pastas (`mkdir`), copiar/mover arquivos (`cp`, `mv`). Nada além disso por enquanto — o objetivo é você perder o medo, não virar especialista em Linux.
- Na parte prática, você vai baixar um arquivo de sequência de DNA e "olhar" seu conteúdo pelo terminal, para sentir na pele que não é bicho de sete cabeças.

### 3. Seu ambiente de trabalho: Conda, Python, R e Jupyter

- **Anaconda/Miniconda (Conda):** pense nisso como um organizador de gavetas. A bioinformática usa dezenas de programas diferentes, cada um com suas próprias dependências (versões de bibliotecas). O Conda cria "ambientes" isolados para você, evitando que uma ferramenta quebre a outra.
- **Python e R:** as duas linguagens que você vai usar ao longo do curso. Nenhuma delas "vem pronta sozinha" — é o Conda que vai instalar e organizar as duas para você, cada uma em sua gaveta.
- **Jupyter Notebook:** seu caderno de laboratório digital. Nele você escreve texto explicativo, código, e vê os resultados (gráficos, tabelas) tudo no mesmo lugar, célula por célula. Vai ser sua ferramenta principal no curso.
- **Google Colab:** uma alternativa ao Jupyter que roda direto no navegador, sem instalação. Se tiver dificuldade para configurar o ambiente local, use-o como plano B.
- **Se você usa Windows:** várias ferramentas do curso (BWA, SAMtools, BLAST em linha de comando, GROMACS) só rodam em Linux. Por isso, você vai instalar o **WSL + Miniconda via VS Code** — o passo a passo está na Parte B, Passo 1.1.

### 4. Python e R aplicados à bioinformática

Você ainda não vai programar de fato nesta aula — a ideia aqui é só você situar as ferramentas que vai encontrar ao longo do curso:

- **Python com Biopython:** biblioteca feita especificamente para você manipular sequências de DNA/proteína e ler arquivos de bancos de dados biológicos. Vai ser seu "canivete suíço" no curso.
- **NumPy e Pandas:** bibliotecas para você lidar com números e tabelas de forma eficiente — como um "Excel programável e muito mais rápido".
- **R e Bioconductor:** outro ambiente que você vai usar bastante, principalmente para análise estatística de dados de expressão gênica (você verá isso com mais profundidade na Unidade 4).

Guarde esta ideia: Python e R não competem entre si — eles são complementares, cada um mais forte em partes diferentes do seu fluxo de trabalho.

### 5. Onde encontrar dados biológicos públicos

Você vai conhecer os principais repositórios onde os dados "moram":

- **NCBI (National Center for Biotechnology Information):** o maior repositório público de dados biológicos dos EUA — sequências, artigos (PubMed), estruturas de proteínas, etc.
- **ENA (European Nucleotide Archive):** o equivalente europeu, focado em dados brutos de sequenciamento.
- **Ensembl:** banco especializado em genomas anotados (já com a informação de "o que cada trecho do DNA faz").

Na parte prática, você vai buscar uma sequência real de um patógeno ou proteína de interesse diretamente no site do NCBI e baixar o arquivo manualmente.

### 6. Automatizando sua busca: APIs do NCBI e do PubChem

**O que é uma API (sem jargão):** pense na API como um garçom. Em vez de você entrar na cozinha (o banco de dados) para pegar a informação, você faz um pedido em um "cardápio" padronizado, e o sistema devolve exatamente o que você pediu — de forma automática e repetível. Isso é bem diferente da busca manual que você fez no tópico anterior: com a API, um script seu pode baixar, filtrar e organizar centenas de sequências sem que você precise clicar em nada — algo essencial quando seu trabalho final exigir análises em escala.

**API do NCBI (Entrez / Biopython)**
- O NCBI disponibiliza um sistema chamado **Entrez**, acessível via API, para você buscar e baixar dados (sequências, artigos, estruturas) diretamente por código.
- Em vez de montar a requisição "na unha", você vai usar o módulo `Bio.Entrez` do **Biopython** — isso simplifica bastante a sintaxe, mesmo se você nunca programou antes.
- Na parte prática, você vai buscar o mesmo gene/proteína que já encontrou manualmente, agora via código, e baixar automaticamente em formato FASTA — dá para comparar visualmente quanto tempo e cliques você economiza.
- Ponto de atenção: o NCBI pede sua identificação (parâmetro `email`) e recomenda uma **API key** gratuita para evitar bloqueios por excesso de requisições — pense nisso como um "crachá de identificação": não é obrigatório, mas evita dor de cabeça.

**API do PubChem**
- **PubChem:** banco público de moléculas pequenas (compostos químicos, fármacos, ligantes) — muito relevante se você trabalha com química de proteínas e macromoléculas, e você vai retomar isso na Unidade 6 (docking).
- **PUG-REST:** forma simplificada de acessar o PubChem via URL — você monta um endereço web seguindo um padrão, e o sistema devolve os dados (estrutura, peso molecular, fórmula) em texto/JSON.
- Na parte prática, você vai buscar um composto conhecido pelo nome e recuperar automaticamente sua estrutura química e identificador (**CID** — Compound ID) — isso vai preparar o terreno para quando você usar esse mesmo composto como ligante no docking molecular (Unidade 6).

### Fechamento do conteúdo teórico

Hoje você vai buscar dados manualmente (NCBI/ENA/Ensembl) e depois vai automatizar essa mesma busca via API — essa automação é o que separa uma consulta pontual de um **pipeline de análise reprodutível**, um conceito que você vai encontrar em todas as unidades seguintes do curso.

---

## PARTE B — SEU ROTEIRO PRÁTICO PASSO A PASSO

**Objetivo:** você vai instalar seu ambiente de trabalho, fazer sua primeira busca manual em bancos de dados públicos e, em seguida, repetir essa busca de forma automatizada via API.

> Você não precisa de nenhum conhecimento prévio de programação. Siga os passos na ordem em que aparecem. Se travar em algum ponto, anote onde parou e traga a dúvida para a próxima aula síncrona.

### Passo 1 — Instale o Anaconda (recomendado para macOS/Linux ou Windows sem WSL)

1. Acesse o site oficial: https://www.anaconda.com/download
2. Baixe a versão correspondente ao seu sistema operacional (Windows, macOS ou Linux).
3. Execute o instalador, mantendo as opções padrão sugeridas.
4. Ao final, abra o **Anaconda Navigator** (ou, no terminal, digite `conda --version` para confirmar que a instalação funcionou).

> **O que é o Conda?** Pense nele como um organizador de gavetas: cada "ambiente" que você criar é uma gaveta separada, com suas próprias ferramentas, sem bagunçar as outras.

**Alternativa sem instalação: Google Colab**

Se você preferir não instalar nada agora, acesse https://colab.research.google.com — é um Jupyter Notebook que roda no navegador, gratuito, bastando você ter uma conta Google.

---

### Passo 1.1 — Se você usa Windows: instale o WSL + Miniconda pelo VS Code (recomendado)

Várias ferramentas que você vai usar ao longo do curso (ex.: BWA, SAMtools, BLAST em linha de comando, GROMACS) foram feitas para rodar em **Linux** e não funcionam nativamente no Windows. A forma mais confiável de você acompanhar o curso sem dor de cabeça é instalar o **WSL (Windows Subsystem for Linux)** — um "Linux dentro do Windows" — pelo **VS Code**, e dentro dele instalar o **Miniconda** (versão enxuta do Anaconda, só com o gerenciador de pacotes e Python básico).

> Se você usa macOS ou Linux nativamente, pode pular esta seção e seguir direto para o Passo 2.

**1.1.1 — Instale o VS Code**

1. Baixe e instale o **Visual Studio Code** em https://code.visualstudio.com/ (versão para Windows, opções padrão).
2. Abra o VS Code após a instalação.

> **O que é o VS Code?** Um editor de código muito usado no mercado, que também funciona como um "painel de controle" para o WSL — você vai abrir o terminal Linux direto de dentro dele, sem precisar alternar entre vários programas.

**1.1.2 — Instale a extensão WSL no VS Code**

1. No VS Code, clique no ícone de **Extensões** na barra lateral esquerda (ícone de quadrados).
2. Pesquise por **WSL** (extensão oficial da Microsoft, ícone verde/azul) e clique em **Install**.
3. Pesquise também por **Remote Development** (pacote de extensões que inclui a WSL) e instale, caso ainda não tenha vindo instalado junto.

**1.1.3 — Instale o WSL com Ubuntu pelo terminal integrado do VS Code**

1. No VS Code, abra o terminal integrado pelo menu **Terminal → New Terminal** (ou atalho `Ctrl + ' `).
2. Confirme que o terminal aberto é o **PowerShell** (aparece indicado no canto superior direito do painel de terminal).
3. Execute o comando:

```powershell
wsl --install
```

4. Isso instala o WSL e, por padrão, a distribuição **Ubuntu**. Reinicie seu computador quando solicitado.
5. Após reiniciar, o Ubuntu vai abrir automaticamente (pode ser pelo próprio VS Code ou pelo menu Iniciar) e vai pedir para você criar um **usuário e senha** (pode ser qualquer nome de usuário; a senha não aparece na tela enquanto você digita — isso é normal).

> **Se der erro sobre virtualização:** isso costuma exigir uma opção na BIOS do seu computador (ex.: "Intel VT-x" ou "AMD-V"). Se não souber habilitar, procure o modelo do seu notebook + "habilitar virtualização BIOS" ou traga a dúvida para a aula.

**1.1.4 — Conecte o VS Code ao Ubuntu (WSL)**

1. No canto inferior esquerdo do VS Code, clique no ícone azul de **conexão remota** (dois colchetes `><`).
2. Selecione a opção **Connect to WSL** (ou **Conectar-se ao WSL**).
3. Uma nova janela do VS Code vai abrir, já conectada ao Ubuntu — repare que o canto inferior esquerdo passa a mostrar **WSL: Ubuntu**.
4. Abra o terminal integrado novamente (`Ctrl + ' `) nesta nova janela: agora o terminal é o **bash do Ubuntu**, não mais o PowerShell. É este terminal que você vai usar daqui em diante.

**1.1.5 — Instale o Miniconda dentro do Ubuntu (WSL)**

Dentro do terminal Ubuntu, execute os comandos abaixo, um de cada vez:

```bash
# Baixa o instalador do Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Executa o instalador
bash Miniconda3-latest-Linux-x86_64.sh
```

- Pressione `Enter` para avançar pelos termos, digite `yes` para aceitar a licença, e `yes` novamente quando for perguntado se deseja inicializar o Miniconda automaticamente.
- Feche e reabra o terminal Ubuntu para que as mudanças tenham efeito.

**1.1.6 — Confirme a instalação**

```bash
conda --version
```

Se aparecer um número de versão (ex.: `conda 24.x.x`), está tudo certo.

> **Por que Miniconda e não Anaconda completo dentro do WSL?** O Miniconda instala só o gerenciador de pacotes, sem dezenas de programas que você ainda não vai usar — deixa seu ambiente mais leve, o que é preferível dentro do WSL. Você vai instalar Python, R, Jupyter e as demais bibliotecas nos próximos passos (2, 3 e 4).

**Como isso muda o resto do seu roteiro:** daqui em diante, sempre que o roteiro mencionar "terminal", abra o **VS Code conectado ao WSL** (janela com "WSL: Ubuntu" no canto inferior esquerdo) e use o terminal integrado (`Ctrl + ' `), que roda bash do Ubuntu — não o Anaconda Prompt nem o PowerShell do Windows.

> **Dica:** o Jupyter Notebook que você abrir dentro do WSL também abre no navegador do Windows normalmente — basta copiar o link que aparece no terminal (algo como `http://localhost:8888/...`) e colar no navegador, como você verá no Passo 6. O VS Code também tem uma extensão própria para abrir notebooks `.ipynb` direto nele, caso prefira não usar o navegador.

---

### Passo 2 — Instale (ou confirme) o Python

> **Se você instalou o Anaconda completo (Passo 1):** o Python já vem incluso — pode pular direto para o item 3 abaixo.
> **Se você instalou o Miniconda (Passo 1.1, Windows/WSL):** o Python básico também já vem junto com o Miniconda, mas vale a pena confirmar e, opcionalmente, organizar um ambiente dedicado.

1. **Verifique se o Python já está disponível:**

```bash
python --version
```

Se aparecer algo como `Python 3.11.x`, está tudo certo — você não precisa instalar nada além disso.

2. **(Recomendado) Crie um ambiente dedicado para o curso.** Em vez de instalar tudo no ambiente `base` do Conda, é boa prática você criar uma "gaveta" só para esta disciplina:

```bash
conda create -n bioinfo python=3.11 -y
conda activate bioinfo
```

> A partir de agora, sempre que você abrir um novo terminal para trabalhar no curso, rode `conda activate bioinfo` primeiro — isso garante que os pacotes que você instalar não se misturem com outros projetos seus.

3. **Se o comando `python` não for reconhecido** (raro, mas pode acontecer em instalações manuais fora do Conda): baixe o instalador oficial em https://www.python.org/downloads/ e, no instalador do Windows, marque a opção **"Add python.exe to PATH"** antes de concluir.

---

### Passo 3 — Instale o R

O R **não vem incluído automaticamente** nem no Anaconda completo nem no Miniconda — você precisa instalá-lo separadamente, pelo próprio Conda ou pelo instalador oficial do R.

**Opção A — Instalar via Conda (recomendado, no mesmo ambiente do Python)**

Com o ambiente `bioinfo` já ativado (Passo 2):

```bash
conda install -c conda-forge r-base r-essentials -y
```

- `r-base`: o núcleo da linguagem R.
- `r-essentials`: um pacote "combo" com as bibliotecas mais usadas em análise de dados em R (equivalente ao que Pandas/NumPy são para o Python).

Confirme a instalação:

```bash
R --version
```

**Opção B — Instalar o R "puro" pelo site oficial (alternativa, fora do Conda)**

1. Acesse https://cran.r-project.org/ e baixe a versão para seu sistema operacional.
2. (Opcional, mas recomendado) Instale também o **RStudio** (https://posit.co/download/rstudio-desktop/), uma interface gráfica dedicada ao R, mais amigável que o terminal puro.

> **Por que preferir a Opção A (Conda)?** Ela mantém Python e R organizados dentro do mesmo ambiente `bioinfo`, evitando conflitos de versão — é a forma mais alinhada com o restante do seu roteiro, que usa Conda em todas as etapas.

**Se quiser desinstalar o R depois:**

```bash
conda activate bioinfo
conda remove r-base r-essentials
```

Se preferir apagar o ambiente inteiro (Python, R e bibliotecas) e recomeçar do zero:

```bash
conda deactivate
conda env remove -n bioinfo
```

---

### Passo 4 — Instale as bibliotecas (Python e R)

Com o ambiente `bioinfo` ativado (`conda activate bioinfo`), instale as bibliotecas que você vai usar ao longo do curso.

**4.1 — Bibliotecas Python**

```bash
pip install jupyter biopython requests pandas numpy
```

O que cada uma faz por você:
- **jupyter:** o caderno digital onde você vai escrever e rodar seu código (precisa instalar manualmente com Miniconda; já vem pronto só no Anaconda completo).
- **biopython:** acessa e manipula dados biológicos (sequências, buscas no NCBI).
- **requests:** faz pedidos a APIs pela internet (você vai usar no PubChem).
- **pandas / numpy:** organizam e calculam dados em tabelas e números, como uma "planilha programável".

Confirme que tudo foi instalado sem erros rodando:

```bash
python -c "import Bio, requests, pandas, numpy; print('Bibliotecas OK')"
```

**4.2 — Bibliotecas R (Bioconductor)**

O **Bioconductor** é o principal repositório de pacotes de bioinformática para R (você vai usá-lo com mais profundidade na Unidade 4 — Transcriptômica). Para instalá-lo, abra o R (digite `R` no terminal) e execute dentro do console do R:

```r
install.packages("BiocManager")
BiocManager::install()
```

> Isso pode demorar alguns minutos na primeira vez — é normal. Para sair do console do R depois, digite `q()` e responda `n` quando for perguntado se deseja salvar a sessão.

---

### Passo 5 — Teste seus primeiros comandos no Terminal

Abra o **Anaconda Prompt** (Windows sem WSL) ou o terminal do **VS Code conectado ao WSL** (Windows com WSL) ou o **Terminal** (macOS/Linux) — lembre-se de rodar `conda activate bioinfo` se você criou o ambiente dedicado no Passo 2 — e teste os comandos abaixo, um de cada vez:

```bash
# Mostra em qual pasta você está
pwd

# Lista os arquivos da pasta atual
ls

# Cria uma pasta nova chamada "bioinfo_aula1"
mkdir bioinfo_aula1

# Entra na pasta criada
cd bioinfo_aula1

# Confirma que está dentro dela
pwd
```

> **Dica:** se algum comando der erro "command not found" no Windows, use o **Anaconda Prompt** ou o terminal do VS Code conectado ao WSL, não o Prompt de Comando padrão do Windows.

---

### Passo 6 — Abra o Jupyter Notebook

No terminal, dentro da pasta `bioinfo_aula1` (e com o ambiente `bioinfo` ativado, se for o seu caso), digite:

```bash
jupyter notebook
```

Isso vai abrir uma aba no seu navegador. Clique em **New → Python 3** para criar seu primeiro notebook.

> **O que é o Jupyter?** Seu caderno digital: cada "célula" pode conter texto explicativo ou código. Você executa uma célula de cada vez com `Shift + Enter` e vê o resultado logo abaixo — ideal para você ir testando aos poucos.

---

### Passo 7 — Faça sua primeira busca manual no NCBI (sem código)

1. Acesse https://www.ncbi.nlm.nih.gov/
2. Na barra de busca, selecione a base **Nucleotide** (ou **Protein**, se preferir).
3. Pesquise por um termo de interesse — sugestão: `Bacillus anthracis toxin` (relevante para biodefesa) ou qualquer gene/proteína do seu interesse de pesquisa.
4. Abra o primeiro resultado relevante e observe:
   - O **número de acesso** (accession number), ex.: `NC_001988.2`
   - A opção **Send to → File → FASTA Format** para baixar a sequência.
5. Baixe o arquivo e abra-o em um editor de texto simples para ver como uma sequência FASTA é organizada (cabeçalho com `>` seguido das bases/aminoácidos).

> Guarde o número de acesso que você encontrou — você vai usá-lo no Passo 8, agora via código.

---

### Passo 8 — Repita a busca, agora via API (NCBI Entrez + Biopython)

**1. Abra o VS Code** e crie um arquivo novo: menu **File → New Text File**.

**2. Salve o arquivo** com o nome `busca_ncbi.py` (`Ctrl + S`, escolha a pasta de sua preferência).

**3. Garanta que o Biopython está instalado no Python certo.** Abra o terminal integrado (menu **Terminal → New Terminal**) e rode:

```bash
python -m pip install biopython requests
```

> **Por que `python -m pip install` e não só `pip install`?** Se o seu computador tem mais de uma instalação de Python (comum em quem já usava Python antes do curso, ou tem Anaconda + Python do Windows ao mesmo tempo), rodar `pip install` sozinho pode instalar o pacote em um Python diferente daquele que o VS Code usa para executar o script — causando o erro `ModuleNotFoundError: No module named 'Bio'`. O comando `python -m pip install` garante que a instalação vai exatamente para o Python que o comando `python` está apontando, evitando essa confusão.

**4. Cole o código abaixo** no arquivo `busca_ncbi.py`, substituindo o e-mail e o número de acesso pelos seus:

```python
from Bio import Entrez, SeqIO

# Identificação obrigatória exigida pelo NCBI (não precisa ser um e-mail "oficial")
Entrez.email = "seu_email@exemplo.com"

# Substitua pelo número de acesso que você encontrou no Passo 7
accession = "NC_001988.2"

# Faz o pedido (o "garçom" leva o pedido até o banco de dados do NCBI)
handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")

# Lê o resultado
registro = SeqIO.read(handle, "fasta")
handle.close()

# Mostra as informações básicas
print("ID:", registro.id)
print("Descrição:", registro.description)
print("Tamanho da sequência:", len(registro.seq), "bases")
print("Primeiras 100 bases:", registro.seq[:100])
```

**5. Salve novamente** (`Ctrl + S`) e **execute** no mesmo terminal:

```bash
python busca_ncbi.py
```

Compare: o mesmo resultado que você baixou manualmente no Passo 7 agora foi obtido automaticamente, em segundos, por código.

> **Sobre a API Key do NCBI (opcional, mas recomendado):** se você for fazer várias buscas seguidas, crie uma chave gratuita em https://www.ncbi.nlm.nih.gov/account/ e adicione `Entrez.api_key = "sua_chave_aqui"` logo abaixo da linha do e-mail. Isso evita que você seja bloqueado por excesso de requisições.

---

### Passo 9 — Busque um composto químico no PubChem (via API)

**1. Crie outro arquivo** no VS Code (**File → New Text File**) e salve como `busca_pubchem.py`.

**2. Cole o código abaixo**, substituindo o nome do composto se quiser:

```python
import requests
import time

# Nome do composto de interesse (pode trocar por outro, ex.: "penicillin", "caffeine")
composto = "ciprofloxacin"

# Monta o endereço da API PUG-REST do PubChem
url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{composto}/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON"

tentativas = 5
espera = 10  # segundos, vai dobrando a cada nova tentativa

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
        print(f"Tentativa {tentativa} falhou (status {resposta.status_code}):", dados.get("Fault", dados))
        if tentativa < tentativas:
            print(f"Aguardando {espera}s antes de tentar novamente...")
            time.sleep(espera)
            espera *= 2
```

**3. Salve** (`Ctrl + S`) e **execute** pelo terminal integrado:

```bash
python busca_pubchem.py
```

Você vai ver informações estruturais do composto retornadas diretamente do PubChem, sem precisar abrir o site manualmente. Como o código já tenta novamente automaticamente em caso de instabilidade do servidor (erro `503 ServerBusy`), você vai ver essas tentativas impressas uma a uma no terminal.

> Guarde o **CID** que você obteve — ele será reutilizado na Unidade 6 (Bioinformática Estrutural), quando esse mesmo composto poderá ser usado como ligante em um experimento de docking molecular.

> **Se mesmo com as tentativas o erro `503` persistir:** teste a mesma URL direto no navegador (substituindo `{composto}` pelo nome do composto). Se o erro aparecer lá também, é uma instabilidade real do lado do PubChem — pause e tente novamente mais tarde, ou continue o roteiro e volte a este passo depois; nada nas etapas seguintes depende dele.

---

## Seu checklist final

- [ ] Anaconda completo (ou Google Colab) instalado e funcionando — **ou** Miniconda + WSL (Windows)
- [ ] (Windows) VS Code + WSL + Miniconda instalados e conectados
- [ ] Ambiente dedicado `bioinfo` criado e ativado (opcional, mas recomendado)
- [ ] Python confirmado (`python --version`)
- [ ] R instalado, via Conda ou instalador oficial (`R --version`)
- [ ] Bibliotecas Python instaladas (jupyter, biopython, requests, pandas, numpy)
- [ ] Bioconductor instalado no R (`BiocManager`)
- [ ] Você conseguiu navegar entre pastas pelo terminal
- [ ] Jupyter Notebook aberto e testado
- [ ] Você encontrou uma sequência manualmente no site do NCBI
- [ ] Você repetiu a mesma busca via `busca_ncbi.py` (Bio.Entrez)
- [ ] Você buscou um composto via `busca_pubchem.py` (PubChem)
- [ ] Você anotou o accession number e o CID que obteve, para usar em unidades futuras

---

## Dúvidas comuns

**"pip install" deu erro de permissão.**
Tente `pip install --user biopython requests` ou confirme que o ambiente `bioinfo` está ativado (`conda activate bioinfo`).

**Ao rodar o Passo 8, aparece `ModuleNotFoundError: No module named 'Bio'`.**
Isso significa que o terminal está rodando o script com um **Python diferente** daquele onde você instalou o Biopython (comum em computadores com mais de uma instalação de Python — ex.: Anaconda + Python do Windows). O Passo 8 já orienta instalar com `python -m pip install biopython requests`, o que resolve isso na maioria dos casos, pois garante que o pacote vai para o mesmo Python usado ao rodar o script. Se ainda assim o erro persistir, confirme rodando `python -m pip show biopython` — se aparecer "Package(s) not found", repita a instalação exatamente no mesmo terminal onde você vai rodar o script.

**O Jupyter não abre no navegador automaticamente.**
Copie o endereço (algo como `http://localhost:8888/...`) que aparece no terminal e cole manualmente no seu navegador.

**Erro "HTTPError" na busca do PubChem.**
Verifique se o nome do composto está em inglês e escrito corretamente (a API do PubChem não busca por nomes em português).

**Aparece `KeyError: 'CanonicalSMILES'`, mesmo com CID, fórmula e peso molecular retornados corretamente.**
O PubChem renomeou essa propriedade recentemente: mesmo pedindo `CanonicalSMILES` na URL, a resposta da API agora traz esse dado sob a chave `ConnectivitySMILES`. O código do roteiro já foi atualizado para usar `propriedades.get("ConnectivitySMILES", "não disponível")` — se você copiou o código antes dessa correção, basta trocar `propriedades["CanonicalSMILES"]` por `propriedades.get("ConnectivitySMILES", "não disponível")` no seu arquivo.

**A busca no PubChem dá erro `503` (`ServerBusy` / "Too many requests or server too busy"), mesmo com as tentativas automáticas do código.**
Esse erro é do lado do servidor do PubChem, não do seu código. Para confirmar, cole a mesma URL usada na célula diretamente no seu navegador (substituindo `{composto}` pelo nome do composto usado). Se o erro aparecer lá também, é uma instabilidade real do serviço — não há nada a corrigir no seu ambiente. Nesse caso: espere alguns minutos e tente de novo, ou continue o roteiro pelos passos seguintes e volte a este depois (nenhuma etapa posterior depende dele).

**Não sei qual e-mail usar no Entrez.email.**
Pode ser seu e-mail institucional do IME ou pessoal — serve apenas para o NCBI entrar em contato caso você faça uso excessivo indevido.

**"wsl --install" deu erro sobre virtualização.**
Verifique nas configurações da BIOS do seu computador se a virtualização (Intel VT-x ou AMD-V) está habilitada.

**O comando `R` não é reconhecido depois que instalei via Conda.**
Confirme que o ambiente `bioinfo` está ativado (`conda activate bioinfo`) — o R instalado por Conda só fica disponível dentro do ambiente onde foi instalado.

**Meu computador tem mais de uma versão de Python instalada — como sei qual está sendo usada, e como troco?**
Isso é comum quando já existia um Python no seu sistema antes de você instalar o Anaconda/Miniconda (ex.: um Python que veio com o Windows, ou um `python3` padrão do Ubuntu). Para verificar qual versão está ativa no momento:

```bash
which python      # mostra o caminho do Python em uso (Linux/macOS/WSL)
where python       # equivalente no PowerShell do Windows
python --version   # mostra o número da versão
```

Para trocar de versão, existem duas situações:

- **Se o problema é qual ambiente Conda está ativo:** rode `conda env list` para ver todos os ambientes existentes e qual está marcado com `*` (o ativo). Troque com `conda activate nome_do_ambiente` (ex.: `conda activate bioinfo`) ou volte ao padrão do sistema com `conda deactivate`.
- **Se você quiser uma versão específica de Python dentro do mesmo ambiente:** recrie o ambiente já especificando a versão desejada, por exemplo:

```bash
conda create -n bioinfo python=3.10 -y
```

> **Dica:** para evitar confusão, sempre confira `python --version` logo depois de ativar um ambiente, antes de rodar qualquer código do curso — é o jeito mais rápido de você confirmar que está na versão certa.

**Travei em um erro que não está nesta lista — posso usar uma IA (Claude, Gemini, etc.) para resolver?**
Sim, e é uma prática legítima e cada vez mais comum em bioinformática. Algumas dicas para você tornar isso produtivo:

- **Copie a mensagem de erro completa**, não só um resumo — ferramentas como Claude ou Gemini diagnosticam muito melhor com o texto exato do erro (a última linha do erro costuma ser a mais informativa, mas o contexto acima ajuda).
- **Informe o ambiente que você está usando**: sistema operacional (Windows/WSL, macOS, Linux), se você está com Anaconda ou Miniconda, e o nome do ambiente Conda ativo — isso evita sugestões genéricas que não se aplicam ao seu caso.
- **Cole o comando que você rodou** logo antes do erro aparecer.
- Trate a resposta da IA como um **ponto de partida, não a palavra final** — principalmente para comandos que apagam arquivos/pacotes (`rm`, `conda remove`, `env remove`): entenda o que o comando faz antes de rodar.
- Se o erro persistir mesmo depois de você tentar a sugestão, traga o print/trecho para a aula síncrona — é normal que o ambiente de cada aluno tenha pequenas particularidades (versão do Windows, permissões, antivírus, etc.) que só se resolvem olhando o seu caso específico.

---

*Traga seu notebook preenchido para a próxima aula síncrona — vamos retomar os accession numbers e CIDs que você encontrou aqui nas unidades seguintes do curso.*
