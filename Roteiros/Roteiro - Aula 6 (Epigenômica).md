# Aula 6 — Epigenômica
### Metilação do DNA, modificações de histonas e acessibilidade da cromatina

**Disciplina:** Introdução à Bioinformática — PPGED / IME
**Duração:** 3 horas (online, síncrona)
**Formato:** teoria → intervalo → prática guiada → conexão e fechamento

---

## Antes da aula (5 minutos de preparo)

Para não perder tempo na prática, deixe tudo pronto antes de entrar:

1. **Conta Google ativa** — vamos usar o **Google Colab** (explicado abaixo).
2. **Navegador atualizado** (Chrome, Firefox ou Edge).
3. Abra estas duas abas e deixe-as carregadas:
   - Google Colab: <https://colab.research.google.com>
   - IGV-Web: <https://igv.org/app>
4. **Não é preciso instalar nada no seu computador.** Tudo roda na nuvem.

> **O que é o Google Colab?** É um "caderno de laboratório digital" que roda no navegador. Cada caderno (chamado *notebook*) é dividido em **células**: você escreve (ou cola) um trecho de código numa célula e clica no botão ▶ à esquerda dela para executá-lo. O código roda num computador do Google, não no seu. O resultado aparece logo abaixo da célula. Você não precisa saber programar para esta aula: todo o código está pronto e comentado; sua tarefa é **executar, observar e interpretar**.

---

## Objetivos da aula

Ao final, você deverá ser capaz de:

- Explicar, com suas palavras, o que é epigenética e por que ela importa.
- Descrever as três camadas epigenômicas estudadas hoje: **metilação do DNA**, **modificações de histonas** e **acessibilidade da cromatina**.
- Saber qual técnica de laboratório mede cada camada (**sequenciamento com bissulfito**, **ChIP-seq**, **ATAC-seq**) e qual ferramenta de bioinformática processa cada uma (**Bismark**, **methylKit**, **MACS**).
- Explorar dados epigenômicos públicos (ENCODE) num navegador genômico (**IGV** / **UCSC**).
- Executar uma chamada de picos de ChIP-seq e uma análise de metilação diferencial com dados públicos.
- Relacionar os resultados com a regulação da expressão gênica e com a química de proteínas (ponte com a Aula 5 — Proteômica).

---

## Mapa da aula

| Horário | Momento | Atividade |
|---|---|---|
| 0:00 – 0:05 | Abertura | Objetivos e a pergunta-guia |
| 0:05 – 1:20 | **Momento 1 — Teoria** | Conceitos, técnicas e ferramentas |
| 1:20 – 1:30 | **Intervalo** | 10 minutos |
| 1:30 – 2:35 | **Momento 2 — Prática** | Navegador genômico, ChIP-seq com MACS3, metilação diferencial |
| 2:35 – 3:00 | **Momento 3 — Conexão e fechamento** | Integração, aplicações em defesa, revisão e tarefas |

**Pergunta-guia da aula:**
> *Todas as células do seu corpo têm (praticamente) o mesmo DNA. Como, então, um neurônio "sabe" que é neurônio e uma célula do fígado "sabe" que é fígado?*

---

# MOMENTO 1 — TEORIA (0:05 – 1:20)

## Bloco 1 — O que é epigenética? (0:05 – 0:20)

### 1.1 Uma analogia para engenheiros

Pense no **genoma** como o **hardware e o firmware** de um equipamento: é o mesmo em todas as unidades que saem da fábrica. Pense no **epigenoma** como o **arquivo de configuração**: ele define quais funções ficam ligadas, quais ficam desligadas e com que intensidade, *sem alterar o firmware*.

- Uma célula do fígado e um neurônio têm o mesmo "firmware" (DNA).
- O que muda é a "configuração" (epigenoma): quais genes podem ser lidos e quais ficam trancados.

**Definição formal:** epigenética é o estudo de alterações **herdáveis** (de uma célula para suas células-filhas, e às vezes entre gerações) na atividade dos genes que **não envolvem mudança na sequência de DNA**.

**Epigenômica** é a mesma ideia em escala de genoma inteiro: em vez de olhar um gene por vez, medimos as marcas epigenéticas em milhões de posições simultaneamente — e para isso precisamos de bioinformática.

### 1.2 Revisão relâmpago: como o DNA cabe no núcleo

- Cada célula humana tem cerca de **2 metros de DNA** dentro de um núcleo com poucos micrômetros.
- O DNA se enrola em "carretéis" de proteínas chamadas **histonas**. Cada carretel tem 8 histonas (duas cópias de H2A, H2B, H3 e H4) e cerca de **147 pares de bases** de DNA enrolados em volta. Esse conjunto é o **nucleossomo**.
- O DNA + histonas + outras proteínas formam a **cromatina**.
- A cromatina pode estar:
  - **aberta (eucromatina)**: frouxa, acessível — genes podem ser lidos;
  - **fechada (heterocromatina)**: compactada — genes ficam silenciados.

### 1.3 As três camadas que estudaremos hoje

| Camada | O que é | Efeito típico | Técnica principal | Ferramenta de bioinformática |
|---|---|---|---|---|
| **Metilação do DNA** | Grupo metil (–CH₃) ligado à citosina | Em promotores, costuma **silenciar** genes | Sequenciamento com bissulfito (WGBS, RRBS) | **Bismark**, **methylKit** |
| **Modificações de histonas** | Marcas químicas nas "caudas" das histonas | Depende da marca: ativa **ou** reprime | ChIP-seq | **MACS2/MACS3** |
| **Acessibilidade da cromatina** | Quão "aberto" está o DNA | Aberto = região regulatória potencialmente ativa | ATAC-seq | **MACS2/MACS3** (e outras) |

### 1.4 Quem escreve, quem lê, quem apaga

As marcas epigenéticas são colocadas e removidas por **enzimas** — ou seja, são processos de **química de proteínas**:

- **Escritores (writers):** colocam a marca. Ex.: DNA metiltransferases (DNMT1, DNMT3A, DNMT3B); histona acetiltransferases (HATs); histona metiltransferases (ex.: EZH2).
- **Apagadores (erasers):** removem a marca. Ex.: enzimas TET (iniciam a remoção da metilação do DNA); histona desacetilases (HDACs); histona desmetilases (KDMs).
- **Leitores (readers):** reconhecem a marca e recrutam outras proteínas. Ex.: proteínas com bromodomínio leem acetilação; cromodomínios leem metilação de lisina.

> **Por que isso importa?** Várias drogas atuam exatamente nessas enzimas (ex.: inibidores de DNMT e de HDAC usados em alguns cânceres hematológicos). Entender o epigenoma é entender um alvo terapêutico.

---

## Bloco 2 — Metilação do DNA (0:20 – 0:40)

### 2.1 O que é

- Em mamíferos, a metilação ocorre principalmente na **citosina (C) seguida de guanina (G)** — o chamado **dinucleotídeo CpG** (o "p" é o fosfato entre elas).
- A citosina vira **5-metilcitosina (5mC)**.
- O genoma humano tem cerca de **28 milhões de CpGs**; a maioria está metilada.
- Exceção importante: as **ilhas CpG**, regiões ricas em CpG frequentemente localizadas em **promotores** de genes. Em geral elas ficam **não metiladas** quando o gene está ativo. Quando uma ilha CpG de promotor fica **metilada**, o gene tende a ser **silenciado**.

### 2.2 Como medir: o truque do bissulfito

O sequenciador não "enxerga" o grupo metil diretamente (na tecnologia Illumina). A solução é química:

1. Trata-se o DNA com **bissulfito de sódio**.
2. **Citosinas NÃO metiladas** são convertidas em **uracila (U)**, que no sequenciamento é lida como **timina (T)**.
3. **Citosinas metiladas (5mC)** ficam **protegidas** e continuam sendo lidas como **C**.

```
Referência:          A C G T T C G A C G
Após bissulfito:     A C G T T T G A C G
                       ↑       ↑     ↑
                    C mantido  C→T  C mantido
                   (metilado) (não metilado) (metilado)
```

**Consequência para a bioinformática:** comparando cada leitura com o genoma de referência, onde a referência tem C e a leitura tem C → metilado; onde a leitura tem T → não metilado.

O **nível de metilação** de um CpG é calculado assim:

> **% metilação = leituras com C ÷ (leituras com C + leituras com T) × 100**

Exemplo: 30 leituras cobrem um CpG; 24 mostram C e 6 mostram T → 80% metilado.

**Limitação importante:** o bissulfito não distingue 5mC de **5hmC** (5-hidroximetilcitosina, um intermediário da desmetilação).

### 2.3 Variações da técnica

| Técnica | Cobertura | Custo | Quando usar |
|---|---|---|---|
| **WGBS** (*Whole-Genome Bisulfite Sequencing*) | Genoma inteiro (~28 milhões de CpGs) | Alto | Estudo completo, padrão-ouro |
| **RRBS** (*Reduced Representation Bisulfite Sequencing*) | Fração enriquecida em ilhas CpG (alguns milhões de CpGs) | Médio | Foco em promotores, muitas amostras |
| **Microarranjos** (ex.: Illumina EPIC) | Centenas de milhares de CpGs pré-selecionados | Baixo | Estudos populacionais com milhares de pessoas |
| **EM-seq** (conversão enzimática) | Genoma inteiro | Alto | Alternativa ao bissulfito que danifica menos o DNA |
| **Nanopore / PacBio** | Leitura direta da base modificada | Variável | Leituras longas, sem conversão química |

### 2.4 Ferramenta: **Bismark**

- **O que é:** programa de linha de comando (gratuito, código aberto, Babraham Institute, Reino Unido) que **alinha leituras tratadas com bissulfito** ao genoma de referência e **determina o estado de metilação** de cada citosina.
- **Qual problema resolve:** um alinhador comum falharia, porque as leituras convertidas têm muitos T onde a referência tem C — pareceriam "erros".
- **Como funciona (ideia):** o Bismark cria versões convertidas *in silico* tanto das leituras quanto do genoma (todos os C → T, e todos os G → A para a fita complementar). Alinha as versões convertidas usando o **Bowtie2** (alinhador visto na Aula 2) e, depois, volta à sequência original para decidir, posição a posição, se o C estava metilado.
- **Principais saídas:**
  - arquivo de alinhamento (BAM);
  - relatório de metilação por citosina (contagens de C e T em cada posição);
  - relatórios-resumo, incluindo a **taxa de conversão do bissulfito** (controle de qualidade: deve ser > 99%).
- **Fluxo típico (apenas para referência — não executaremos hoje, pois exige o genoma completo e horas de processamento):**

```bash
# 1. Preparar o genoma de referência (feito uma única vez)
bismark_genome_preparation pasta_do_genoma/

# 2. Alinhar as leituras tratadas com bissulfito
bismark --genome pasta_do_genoma/ amostra.fastq.gz

# 3. Extrair o estado de metilação de cada citosina
bismark_methylation_extractor --bedGraph --gzip amostra_bismark_bt2.bam
```

> Lembre: antes do Bismark, os dados passam pelo mesmo controle de qualidade da Aula 2 (**FastQC**, **Trimmomatic/Cutadapt**).

### 2.5 Ferramenta: **methylKit**

- **O que é:** pacote para a linguagem **R**, distribuído pelo **Bioconductor** (repositório de pacotes de bioinformática em R visto na Aula 1).
- **O que faz:** recebe as tabelas de metilação (por exemplo, as geradas pelo Bismark) e realiza a **análise estatística**:
  - estatísticas descritivas (distribuição de % de metilação e de cobertura);
  - **filtro por cobertura** (descarta CpGs com poucas leituras, pouco confiáveis);
  - união das amostras nos CpGs em comum;
  - correlação, agrupamento e PCA entre amostras;
  - identificação de **CpGs diferencialmente metilados (DMCs)** e **regiões diferencialmente metiladas (DMRs)** entre grupos (ex.: tratado × controle).
- **Estatística:** por padrão usa **regressão logística**; quando há só uma amostra por grupo, usa o **teste exato de Fisher**. Corrige os valores-p para múltiplos testes (valor-q).
- **Critério usual de "diferencialmente metilado":** diferença de metilação ≥ 25 pontos percentuais **e** valor-q < 0,01.

> **Por que o valor-q?** Testamos milhares de CpGs ao mesmo tempo. Com valor-p < 0,05, esperaríamos 5% de falsos positivos — em 1 milhão de CpGs, seriam 50 mil "achados" falsos! O valor-q controla a **taxa de falsas descobertas (FDR)**. Mesma lógica do DESeq2 na Aula 4.

---

## Bloco 3 — Modificações de histonas e ChIP-seq (0:40 – 1:00)

### 3.1 O "código das histonas"

As histonas têm "caudas" (extremidades) que saem do nucleossomo. Os aminoácidos dessas caudas recebem modificações químicas pós-traducionais (**PTMs**, as mesmas da Aula 5). Os alvos principais são **lisinas (K)** e **argininas (R)**, mas também **serinas (S)**, **treoninas (T)**, **tirosinas (Y)**, **histidinas (H)**, **glutaminas (Q)**, **cisteínas (C)** e **prolinas (P)**. As modificações também ocorrem no "miolo" globular da histona, não só nas caudas.

Hoje conhecemos **mais de 20 tipos químicos** de modificação, em **centenas de sítios**. A lista continua crescendo: novos tipos foram descobertos quase todo ano desde 2011, e não existe um catálogo "fechado". O **Apêndice C** traz o catálogo completo do que já é conhecido, para consulta. Na aula, vamos nos concentrar nas marcas mais usadas em ChIP-seq.

**Como ler o nome de uma marca:** `H3K27me3` =

- **H3** → histona H3
- **K27** → lisina (K) na posição 27
- **me3** → trimetilação (três grupos metil)

**Marcas mais estudadas:**

| Marca | Onde aparece | Significado |
|---|---|---|
| **H3K4me3** | Promotores | Gene **ativo** ou pronto para ser ativado |
| **H3K27ac** | Promotores e *enhancers* (amplificadores) | Região regulatória **ativa** |
| **H3K4me1** | *Enhancers* | *Enhancer* ativo ou "em espera" |
| **H3K36me3** | Corpo dos genes | Gene sendo **transcrito** |
| **H3K27me3** | Promotores e grandes domínios | **Repressão** (via complexo Polycomb, enzima EZH2) |
| **H3K9me3** | Heterocromatina | **Repressão** estável |

**As grandes famílias de modificações de histonas** (detalhes de cada uma no Apêndice C):

| Família | Exemplos | Ideia central |
|---|---|---|
| **Pequenos grupos químicos** | Acetilação, metilação (de lisina e de arginina), fosforilação | Mudam a carga ou o "formato" do aminoácido e criam pontos de encaixe para proteínas leitoras |
| **Acilações "do metabolismo"** | Crotonilação, butirilação, succinilação, lactilação, β-hidroxibutirilação… | Moléculas do metabolismo viram marcas: o estado nutricional da célula fica "escrito" na cromatina |
| **Proteínas inteiras** | Ubiquitinação, SUMOilação | Uma pequena proteína é pendurada na histona |
| **Grupos grandes e incomuns** | ADP-ribosilação, O-GlcNAcilação, serotonilação, dopaminilação | Sinais de dano ao DNA, de nutrientes, de neurotransmissores |
| **Alterações do próprio aminoácido** | Citrulinação, isomerização de prolina, corte da cauda (*clipping*) | A cadeia lateral ou a própria cauda é transformada |
| **Não enzimáticas** | Glicação, formilação | "Acidentes" químicos causados por metabólitos reativos |

> **Metilação de arginina tem uma pegadinha:** a arginina pode receber dois grupos metil no **mesmo** nitrogênio (dimetilação **assimétrica**, `me2a`) ou um em **cada** nitrogênio (dimetilação **simétrica**, `me2s`). Na mesma posição, as duas formas podem ter efeitos opostos: `H4R3me2a` costuma **ativar** e `H4R3me2s` costuma **reprimir**.

### 3.2 Ponte com a Aula 5 (Proteômica) e com a química de proteínas

As marcas de histona são PTMs em proteínas e também podem ser estudadas por **espectrometria de massas**. Um exemplo que mostra por que a precisão do espectrômetro importa:

| Modificação | Variação de massa |
|---|---|
| Metilação (me1) | +14,016 Da |
| **Trimetilação (me3)** | **+42,047 Da** |
| **Acetilação (ac)** | **+42,011 Da** |

Trimetilação e acetilação diferem em apenas **0,036 Da**! Só instrumentos de alta resolução (como os Orbitrap usados com o MaxQuant) conseguem separá-las com segurança. Os anticorpos usados no ChIP-seq, por outro lado, distinguem as marcas pela **forma** do epítopo, não pela massa.

### 3.3 Técnica: **ChIP-seq** (*Chromatin ImmunoPrecipitation + sequencing*)

**Pergunta que o ChIP-seq responde:** *em quais posições do genoma está a proteína X (ou a marca de histona Y)?*

Passo a passo:

1. **Fixar:** formaldeído "cola" (liga covalentemente) as proteínas ao DNA onde elas estão.
2. **Fragmentar:** a cromatina é quebrada em pedaços de ~200–500 pb (por ultrassom ou enzima).
3. **Imunoprecipitar:** um **anticorpo** específico (ex.: anti-H3K27ac ou anti-CTCF) "pesca" só os fragmentos que contêm o alvo.
4. **Liberar e purificar** o DNA pescado.
5. **Sequenciar** esses fragmentos.
6. **Bioinformática:** alinhar as leituras ao genoma (**Bowtie2** ou **BWA**, Aula 2) e encontrar onde elas se acumulam → esses acúmulos são os **picos**.

**O controle (input):** uma amostra processada da mesma forma, **sem** o anticorpo. Ele mostra o "ruído de fundo" (regiões que naturalmente fragmentam mais, regiões repetitivas). Um pico verdadeiro deve ser **significativamente maior que o controle**.

**Variações modernas:** **CUT&RUN** e **CUT&Tag** fazem o mesmo papel com menos células e menos ruído.

### 3.4 Ferramenta: **MACS** (MACS2 / MACS3)

- **Nome:** *Model-based Analysis of ChIP-Seq*. É a ferramenta mais usada no mundo para **chamada de picos** (*peak calling*). Gratuita, código aberto, em Python.
- **MACS2 × MACS3:** o MACS3 é o sucessor do MACS2, mantido pela mesma equipe. O comando principal (`callpeak`) e os arquivos de saída são os mesmos. **Usaremos o MACS3 na prática.**
- **Como funciona (ideia):**
  1. **Estima o tamanho do fragmento:** leituras de uma fita (+) e da outra (−) se acumulam em lados opostos do sítio de ligação; a distância entre esses dois "montes" revela o tamanho médio do fragmento. (Se não houver sinal suficiente, informamos o tamanho manualmente — veremos isso na prática!)
  2. **Modela o fundo** com uma distribuição de **Poisson** (modelo estatístico para contagens de eventos raros), ajustada **localmente** a partir do controle.
  3. **Testa cada região:** "a quantidade de leituras aqui é maior do que o esperado pelo acaso?" → gera valor-p e valor-q.
- **Dois modos:**
  - **narrow** (padrão): picos estreitos — fatores de transcrição, H3K4me3, H3K27ac;
  - **broad** (`--broad`): domínios largos — H3K27me3, H3K9me3, H3K36me3.
- **Principal saída: arquivo `.narrowPeak`**, uma tabela com uma linha por pico:

| Coluna | Conteúdo |
|---|---|
| 1–3 | Cromossomo, início, fim do pico |
| 4 | Nome do pico |
| 5 | Pontuação (para exibição) |
| 6 | Fita (geralmente ".") |
| 7 | **Enriquecimento** (*fold enrichment*): quantas vezes o sinal é maior que o fundo |
| 8 | −log₁₀(valor-p) |
| 9 | **−log₁₀(valor-q)** — quanto maior, mais confiável |
| 10 | Posição do **cume** (*summit*), contada a partir do início do pico |

> **Lendo −log₁₀:** valor-q = 0,01 → 2; valor-q = 0,0001 → 4. **Quanto maior o número, mais significativo.**

---

## Bloco 4 — Acessibilidade da cromatina: ATAC-seq (1:00 – 1:10)

### 4.1 A ideia

Regiões regulatórias ativas (promotores, *enhancers*) precisam estar **abertas** para que fatores de transcrição se liguem. Medir onde a cromatina está aberta revela **onde está a "maquinaria de controle" de cada tipo celular**.

### 4.2 Técnica: **ATAC-seq** (*Assay for Transposase-Accessible Chromatin*)

- Usa uma enzima chamada **transposase Tn5**, carregada com adaptadores de sequenciamento.
- A Tn5 "corta e cola" esses adaptadores **preferencialmente no DNA aberto** (o DNA enrolado nos nucleossomos fica protegido).
- Sequenciando os fragmentos, vemos picos exatamente onde a cromatina está acessível.

**Vantagens:** rápida, não precisa de anticorpo, funciona com poucas células (inclusive célula única: *scATAC-seq*).

**Assinatura de qualidade:** a distribuição de tamanhos dos fragmentos mostra "degraus": fragmentos curtos (< 100 pb, regiões livres de nucleossomos) e picos em ~200 pb, ~400 pb (um, dois nucleossomos). Esse padrão em escada indica um bom experimento.

**Análise:** mesmo início da ChIP-seq (QC → alinhamento) e, para chamar os picos, o **MACS** também é muito usado. Antecessores: **DNase-seq** e **FAIRE-seq**.

### 4.3 Juntando as três camadas

| Estado do gene | Metilação do promotor | Histonas | Cromatina (ATAC) | Expressão (RNA-seq, Aula 4) |
|---|---|---|---|---|
| **Ativo** | Baixa | H3K4me3, H3K27ac | Aberta (pico) | Alta |
| **Reprimido (Polycomb)** | Baixa a variável | H3K27me3 | Fechada | Baixa |
| **Silenciado estável** | **Alta** | H3K9me3 | Fechada | Nenhuma |

---

## Bloco 5 — Bancos de dados e visualização (1:10 – 1:20)

### 5.1 Bancos de dados públicos

- **ENCODE** (*Encyclopedia of DNA Elements*, <https://www.encodeproject.org>): consórcio internacional que gerou milhares de experimentos (ChIP-seq, ATAC-seq, DNase-seq, WGBS, RNA-seq) em centenas de tipos celulares humanos e de camundongo. Dados já processados e com controle de qualidade padronizado.
- **Roadmap Epigenomics** (<https://egg2.wustl.edu/roadmap/web_portal/>): mapeou mais de uma centena de **epigenomas de referência** de tecidos e células humanas (inclusive tecidos primários e células-tronco). Ficou famoso pelos "estados de cromatina" — o genoma classificado em categorias como "promotor ativo", "enhancer", "reprimido por Polycomb".
- **GEO** (*Gene Expression Omnibus*, NCBI, <https://www.ncbi.nlm.nih.gov/geo/>): repositório onde pesquisadores depositam dados de publicações. Tem dados brutos e processados de praticamente qualquer estudo epigenômico publicado. (Os dados brutos ficam no **SRA**, visto na Aula 2.)

### 5.2 Ferramentas de visualização

- **IGV** (*Integrative Genomics Viewer*, Broad Institute): visualizador de genomas gratuito. Mostra o genoma como uma "régua" horizontal e empilha embaixo **trilhas (*tracks*)**: genes, leituras alinhadas, sinais de ChIP-seq, picos, metilação. Existe a versão para instalar (*IGV Desktop*) e a **versão web (IGV-Web)**, que usaremos — roda no navegador e já oferece acesso direto a trilhas do ENCODE.
- **UCSC Genome Browser** (Universidade da Califórnia em Santa Cruz, <https://genome.ucsc.edu>): navegador genômico online com centenas de trilhas prontas (genes, conservação, regulação, ENCODE, metilação). Excelente para consultas rápidas.

**Formatos de arquivo que aparecerão:**

| Formato | O que contém | Exemplo de uso |
|---|---|---|
| **BED / narrowPeak** | Lista de regiões (cromossomo, início, fim, ...) | Picos do MACS |
| **bedGraph / bigWig** | Um valor contínuo ao longo do genoma | Sinal de ChIP-seq, % de metilação |
| **BAM** | Leituras alinhadas | Saída do Bowtie2 / Bismark |

> **Atenção às versões do genoma:** o genoma humano de referência tem versões (**hg19** = GRCh37; **hg38** = GRCh38). As coordenadas mudam entre versões! Sempre confira se seus dados e o navegador estão na mesma versão.

**Antes do intervalo:** abra o Colab e crie um novo notebook (*Arquivo → Novo notebook*).

---

# INTERVALO (1:20 – 1:30)

☕ 10 minutos.

---

# MOMENTO 2 — PRÁTICA (1:30 – 2:35)

Três atividades, da mais visual para a mais analítica:

| Atividade | Tempo | Ferramentas | Dados públicos |
|---|---|---|---|
| **P1** — Explorando o epigenoma no navegador | 15 min | IGV-Web (ou UCSC) | ENCODE (células GM12878) |
| **P2** — Encontrando picos de ChIP-seq | 25 min | Google Colab + MACS3 | ChIP-seq de CTCF, cromossomo 22 (dados de teste do MACS) |
| **P3** — Metilação diferencial | 25 min | Google Colab + Python (pandas, SciPy) | Dados de exemplo do methylKit |

Todas as ferramentas são **gratuitas e de código aberto**.

---

## P1 — Explorando o epigenoma no navegador (1:30 – 1:45)

**Objetivo:** ver, com seus próprios olhos, a "assinatura" epigenética de um gene ativo e de um gene reprimido.

**Contexto:** usaremos dados do ENCODE para a linhagem celular **GM12878** (células linfoblastoides humanas, derivadas de linfócitos B). É uma das linhagens mais estudadas do mundo.

### Passo a passo (IGV-Web)

1. Acesse <https://igv.org/app>.
2. No menu **Genome**, selecione **Human (GRCh37/hg19)**.
   *(Usaremos hg19 porque os dados da atividade P2 estão nessa versão.)*
3. No menu **Tracks**, escolha **ENCODE Signals** (o nome do submenu pode variar ligeiramente). Abre-se uma tabela pesquisável.
4. Na caixa de busca, digite `GM12878` e selecione as trilhas de sinal (*signal*) para:
   - **H3K4me3** (promotor ativo)
   - **H3K27ac** (região regulatória ativa)
   - **H3K27me3** (repressão Polycomb)
   - **CTCF** (proteína "isoladora", usaremos na P2)
   - Se disponível: **DNase** (cromatina aberta)
5. Clique em **OK**. As trilhas aparecem empilhadas.
6. Na caixa de busca de posição (no topo), digite **`GAPDH`** e pressione Enter.
   - *GAPDH* é um gene "de manutenção" (*housekeeping*), ativo em praticamente todas as células.
7. Observe e anote: em quais trilhas há sinal alto no início do gene?
8. Agora digite **`HOXD13`** (ou `HOXD10`).
   - Os genes *HOX* controlam o desenvolvimento embrionário e ficam reprimidos em células sanguíneas adultas.
9. Observe e anote: o padrão mudou? Qual trilha "acendeu" agora?

### Plano B (se o IGV-Web não carregar as trilhas)

1. Acesse <https://genome.ucsc.edu> → **Genome Browser** → selecione **Human, hg19**.
2. Digite `GAPDH` e clique em **go**.
3. Abaixo da imagem, na seção **Regulation**, abra **ENCODE Regulation** e ative as trilhas de histonas (*Layered H3K4Me3*, *Layered H3K27Ac*) e de DNase.
4. Repita para `HOXD13`.

### Registre suas observações

| Trilha | GAPDH (sinal alto/baixo) | HOXD13 (sinal alto/baixo) |
|---|---|---|
| H3K4me3 | | |
| H3K27ac | | |
| H3K27me3 | | |
| DNase / cromatina aberta | | |

> **O que se espera:** em *GAPDH*, H3K4me3 e H3K27ac altos e cromatina aberta no promotor; H3K27me3 baixo. No grupo *HOXD*, um grande "tapete" de H3K27me3 e pouco sinal ativo — repressão por Polycomb. **Confirme ou refute com o que você viu.**

**Não feche o IGV-Web** — vamos voltar a ele na P2.

---

## P2 — Encontrando picos de ChIP-seq com MACS3 (1:45 – 2:10)

**Objetivo:** executar uma chamada de picos real e interpretar a saída.

**Os dados:** leituras de ChIP-seq da proteína **CTCF** e do respectivo controle (*input*), restritas a um trecho do **cromossomo 22 humano (hg19)**, com cerca de 50 mil leituras cada. São os dados de teste oficiais distribuídos pelo projeto MACS (derivados de dados públicos do ENCODE), pequenos o bastante para rodar em segundos.

> **Quem é a CTCF?** Uma proteína que se liga ao DNA numa sequência específica e funciona como "isolador": organiza o genoma em alças e domínios, impedindo que um *enhancer* ative o gene errado. Por ligar-se a um ponto preciso, gera **picos estreitos e nítidos** — ideal para aprender.

**Formato dos arquivos:** BED (uma linha por leitura alinhada: cromossomo, início, fim, ..., fita). Ou seja, o alinhamento (Bowtie2/BWA) já foi feito; começamos a partir dele.

No Colab, crie uma célula para cada bloco abaixo (botão **+ Código**), cole o conteúdo e execute com ▶. Linhas que começam com `#` são **comentários** — explicações para humanos, ignoradas pelo computador. Linhas que começam com `!` são **comandos de terminal** (Aula 1) executados de dentro do Colab.

### Célula 1 — Instalar o MACS3

```python
# Instala o MACS3 no computador do Google (leva ~1 minuto)
!pip install -q macs3
# Confere a versão instalada
!macs3 --version
```

### Célula 2 — Baixar os dados

```python
# Baixa as leituras de ChIP-seq (tratamento) e do controle (input)
url = "https://raw.githubusercontent.com/macs3-project/MACS/master/test/"
!wget -q {url}CTCF_SE_ChIP_chr22_50k.bed.gz
!wget -q {url}CTCF_SE_CTRL_chr22_50k.bed.gz
# Mostra as 5 primeiras leituras (cromossomo, início, fim, ..., fita)
!zcat CTCF_SE_ChIP_chr22_50k.bed.gz | head -5
# Conta quantas leituras há no arquivo
!zcat CTCF_SE_ChIP_chr22_50k.bed.gz | wc -l
```

**Pergunte-se:** quantas leituras há? Em que cromossomo estão? *(No teste desta aula: 49.622 leituras, todas no chr22.)*

### Célula 3 — Primeira tentativa de chamada de picos

```python
# -t: arquivo de tratamento (ChIP)    -c: arquivo de controle (input)
# -f BED: formato dos arquivos         -g hs: tamanho do genoma humano (homo sapiens)
# -n ctcf: prefixo dos arquivos de saída   -q 0.05: limite de valor-q
!macs3 callpeak -t CTCF_SE_ChIP_chr22_50k.bed.gz -c CTCF_SE_CTRL_chr22_50k.bed.gz \
    -f BED -g hs -n ctcf --outdir resultado -q 0.05
```

**Vai aparecer um AVISO (WARNING) e nenhum arquivo de picos será gerado!** Leia a mensagem: o MACS3 não encontrou pares de picos suficientes (ele precisa de pelo menos 100) para **estimar sozinho o tamanho do fragmento** (passo 1 do funcionamento, Bloco 3.4). Com apenas 50 mil leituras, o sinal é fraco demais para essa estimativa. **Isso é normal com poucos dados** — e a própria mensagem sugere a solução: usar `--nomodel` e `--extsize`.

> **Lição para a vida real:** mensagens de aviso de programas de bioinformática quase sempre dizem o que está errado e como corrigir. Leia-as antes de pedir ajuda.

### Célula 4 — Chamada de picos informando o tamanho do fragmento

```python
# --nomodel: não tente estimar o tamanho do fragmento
# --extsize 200: considere que cada fragmento tem ~200 pb (valor típico de ChIP-seq)
!macs3 callpeak -t CTCF_SE_ChIP_chr22_50k.bed.gz -c CTCF_SE_CTRL_chr22_50k.bed.gz \
    -f BED -g hs -n ctcf --outdir resultado -q 0.05 --nomodel --extsize 200
# Lista os arquivos gerados
!ls -lh resultado/
```

**Arquivos gerados:**

- `ctcf_peaks.narrowPeak` → **a lista de picos** (o mais importante);
- `ctcf_summits.bed` → a posição exata do cume de cada pico;
- `ctcf_peaks.xls` → a mesma lista, em formato de tabela legível (apesar da extensão, é texto).

### Célula 5 — Ler e explorar os picos

```python
import pandas as pd   # pandas: biblioteca de tabelas (Aula 1)

colunas = ["cromossomo", "inicio", "fim", "nome", "pontuacao", "fita",
           "enriquecimento", "menos_log10_p", "menos_log10_q", "cume"]
picos = pd.read_csv("resultado/ctcf_peaks.narrowPeak", sep="\t", names=colunas)
picos["largura"] = picos["fim"] - picos["inicio"]

print("Número de picos encontrados:", len(picos))
print("Largura média dos picos (pb):", round(picos["largura"].mean()))
# Os 5 picos mais significativos (maior -log10 do valor-q)
picos.sort_values("menos_log10_q", ascending=False).head(5)
```

> **Resultado obtido no teste desta aula (MACS3 3.0.5):** **802 picos**, com largura média de ~390 pb. O pico mais forte fica em **chr22:37.252.291–37.252.808**, com enriquecimento de ~58 vezes sobre o fundo.

### Célula 6 — Visualizar

```python
import matplotlib.pyplot as plt   # matplotlib: biblioteca de gráficos

fig, eixos = plt.subplots(1, 2, figsize=(11, 4))
eixos[0].hist(picos["largura"], bins=40, color="steelblue")
eixos[0].set_title("Largura dos picos")
eixos[0].set_xlabel("Largura (pb)"); eixos[0].set_ylabel("Número de picos")

eixos[1].scatter(picos["enriquecimento"], picos["menos_log10_q"], s=8, color="darkred")
eixos[1].set_title("Enriquecimento × significância")
eixos[1].set_xlabel("Enriquecimento sobre o controle (vezes)")
eixos[1].set_ylabel("-log10(valor-q)")
plt.tight_layout(); plt.show()
```

**Interprete:** os picos são estreitos (centenas de pb) ou largos (milhares)? Isso combina com o que se espera de um fator como a CTCF? Picos com maior enriquecimento são também os mais significativos?

### Célula 7 — Baixar o arquivo de picos e ver no IGV-Web

```python
from google.colab import files
files.download("resultado/ctcf_peaks.narrowPeak")
```

1. Volte à aba do **IGV-Web** (genoma **hg19**).
2. Menu **Tracks → Local File** → selecione `ctcf_peaks.narrowPeak`.
3. Na caixa de posição, digite `chr22:37,250,000-37,256,000` (região do pico mais forte).
4. **Compare** os seus picos com a trilha de CTCF do ENCODE carregada na P1. Eles coincidem?

> **Esse é o momento-chave da prática:** você acabou de reproduzir, com uma ferramenta de código aberto, um resultado de um consórcio internacional.

---

## P3 — Metilação diferencial (2:10 – 2:35)

**Objetivo:** reproduzir a lógica do **methylKit**: calcular % de metilação, filtrar, comparar grupos e encontrar CpGs diferencialmente metilados (DMCs).

**Os dados:** os arquivos de exemplo oficiais do methylKit — 4 amostras (2 "teste" e 2 "controle") de sequenciamento com bissulfito de um trecho do **cromossomo 21 humano** (versão hg18). Cada linha é um CpG, com: posição, **cobertura** (nº de leituras), **freqC** (% de C = metilado) e **freqT** (% de T = não metilado). É exatamente o tipo de tabela que sai do Bismark.

> **Por que Python e não o methylKit (R) diretamente?** Instalar o methylKit no Colab exige compilar dezenas de pacotes do Bioconductor e pode levar mais de 15 minutos. Para caber na aula, faremos **os mesmos passos** com bibliotecas Python leves (**pandas** para tabelas, **SciPy** e **statsmodels** para estatística — todas já vêm instaladas no Colab). O script equivalente em R/methylKit está no **Apêndice A**, para você executar em casa.

### Célula 8 — Carregar as quatro amostras

```python
import pandas as pd, numpy as np

base = "https://raw.githubusercontent.com/al2na/methylKit/master/inst/extdata/"
amostras = {"teste1": "test1.myCpG.txt",   "teste2": "test2.myCpG.txt",
            "controle1": "control1.myCpG.txt", "controle2": "control2.myCpG.txt"}

dados = {}
for nome, arquivo in amostras.items():
    df = pd.read_csv(base + arquivo, sep="\t")
    # Converte as porcentagens em contagens de leituras metiladas / não metiladas
    df["metilados"] = (df["coverage"] * df["freqC"] / 100).round().astype(int)
    df["nao_metilados"] = df["coverage"] - df["metilados"]
    dados[nome] = df
    print(f"{nome}: {len(df)} CpGs | cobertura mediana = {df.coverage.median():.0f} "
          f"| metilação média = {df.freqC.mean():.1f}%")

dados["teste1"].head()
```

### Célula 9 — Distribuição da metilação (equivalente a `getMethylationStats`)

```python
import matplotlib.pyplot as plt
fig, eixos = plt.subplots(1, 4, figsize=(14, 3), sharey=True)
for eixo, (nome, df) in zip(eixos, dados.items()):
    eixo.hist(df["freqC"], bins=20, color="seagreen")
    eixo.set_title(nome); eixo.set_xlabel("% metilação")
eixos[0].set_ylabel("Número de CpGs")
plt.tight_layout(); plt.show()
```

**Interprete:** o histograma tem forma de "U" (muitos CpGs perto de 0% e de 100%)? Esse padrão **bimodal** é típico: em cada célula, um CpG está metilado ou não; valores intermediários refletem mistura de células ou de alelos.

### Célula 10 — Filtrar e unir (equivalente a `filterByCoverage` e `unite`)

```python
# Mantém só CpGs com pelo menos 10 leituras (abaixo disso, a % é pouco confiável)
filtrados = {n: d[d["coverage"] >= 10] for n, d in dados.items()}

# Une as 4 amostras, mantendo só os CpGs presentes em todas
tabela = None
for n, d in filtrados.items():
    parte = d[["chrBase", "chr", "base", "metilados", "nao_metilados", "freqC"]].rename(
        columns={"metilados": f"met_{n}", "nao_metilados": f"naomet_{n}", "freqC": f"pct_{n}"})
    tabela = parte if tabela is None else tabela.merge(parte, on=["chrBase", "chr", "base"])

print("CpGs presentes nas 4 amostras:", len(tabela))
# Correlação entre amostras (equivalente a getCorrelation)
tabela[[c for c in tabela.columns if c.startswith("pct_")]].corr().round(2)
```

**Interprete a matriz de correlação:** as réplicas de um mesmo grupo (teste1 × teste2; controle1 × controle2) são mais parecidas entre si do que com o outro grupo? Isso é um bom sinal de qualidade experimental.

### Célula 11 — Teste estatístico (equivalente a `calculateDiffMeth`)

```python
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Soma as leituras das réplicas de cada grupo
t = tabela
t["met_teste"]   = t.met_teste1 + t.met_teste2
t["naomet_teste"] = t.naomet_teste1 + t.naomet_teste2
t["met_ctrl"]    = t.met_controle1 + t.met_controle2
t["naomet_ctrl"] = t.naomet_controle1 + t.naomet_controle2

t["pct_teste"] = 100 * t.met_teste / (t.met_teste + t.naomet_teste)
t["pct_ctrl"]  = 100 * t.met_ctrl  / (t.met_ctrl  + t.naomet_ctrl)
t["diferenca"] = t.pct_teste - t.pct_ctrl   # positiva = mais metilado no grupo teste

# Teste exato de Fisher em cada CpG: a proporção metilada difere entre grupos?
t["p"] = [fisher_exact([[a, b], [c, d]])[1] for a, b, c, d in
          zip(t.met_teste, t.naomet_teste, t.met_ctrl, t.naomet_ctrl)]
# Correção para múltiplos testes (Benjamini-Hochberg) -> valor-q
t["q"] = multipletests(t["p"], method="fdr_bh")[1]
print("Teste concluído para", len(t), "CpGs")
```

> **O teste exato de Fisher, em uma frase:** para cada CpG, monta-se uma tabela 2×2 (metilado/não metilado × teste/controle) e calcula-se a probabilidade de uma diferença tão grande surgir por acaso. *Nota:* o methylKit usa, por padrão, **regressão logística**, que considera a variação entre réplicas; nossa versão simplificada soma as réplicas. Os resultados são parecidos, mas não idênticos.

### Célula 12 — DMCs e gráfico de vulcão (equivalente a `getMethylDiff`)

```python
dmc = t[(t["diferenca"].abs() >= 25) & (t["q"] < 0.01)]
print("CpGs diferencialmente metilados (|diferença| ≥ 25% e q < 0,01):", len(dmc))
print("  Hipermetilados no grupo teste:", (dmc.diferenca > 0).sum())
print("  Hipometilados no grupo teste: ", (dmc.diferenca < 0).sum())

plt.figure(figsize=(6, 5))
cor = np.where(t.index.isin(dmc.index), "crimson", "lightgray")
plt.scatter(t["diferenca"], -np.log10(t["q"] + 1e-300), c=cor, s=8)
plt.axvline(25, ls="--", c="k"); plt.axvline(-25, ls="--", c="k")
plt.axhline(2, ls="--", c="k")
plt.xlabel("Diferença de metilação (teste − controle, %)")
plt.ylabel("-log10(valor-q)"); plt.title("Gráfico de vulcão — metilação diferencial")
plt.show()

dmc.sort_values("q")[["chr", "base", "pct_teste", "pct_ctrl", "diferenca", "q"]].head(10)
```

> **Resultado obtido no teste desta aula:** **963** CpGs em comum; **127 DMCs** (72 hipermetilados e 55 hipometilados no grupo teste). As réplicas de cada grupo tiveram correlação de 0,93–0,95, contra ~0,88 entre grupos.

**Interprete:** o que significa um ponto no canto superior direito do vulcão? E no canto superior esquerdo? Por que a maioria dos pontos fica cinza perto do centro?

---

# MOMENTO 3 — CONEXÃO ENTRE TEORIA E PRÁTICA E FECHAMENTO (2:35 – 3:00)

## 3.1 O que cada prática mostrou (2:35 – 2:45)

| Prática | Conceito teórico | O que você viu |
|---|---|---|
| **P1** | Marcas ativas × repressivas; código das histonas | *GAPDH* com H3K4me3/H3K27ac; *HOXD* coberto por H3K27me3 |
| **P2** | ChIP-seq, controle, modelo de Poisson, valor-q | ~800 picos estreitos de CTCF que coincidem com o ENCODE; a necessidade de ajustar o tamanho do fragmento quando há poucos dados |
| **P3** | Bissulfito, % de metilação, cobertura, FDR | Distribuição bimodal; réplicas concordantes; ~127 DMCs |

**O fluxo completo, do laboratório ao resultado:**

```
Amostra → técnica de laboratório → sequenciamento → FASTQ
   → controle de qualidade (FastQC, Trimmomatic)          [Aula 2]
   → alinhamento (Bowtie2/BWA; Bismark para bissulfito)   [Aula 2 + hoje]
   → quantificação (MACS para picos; methylKit p/ metilação) [hoje]
   → visualização (IGV, UCSC) e integração com RNA-seq     [hoje + Aula 4]
```

### Discussão em grupo (2 perguntas)

1. Se um gene tem **promotor muito metilado**, **H3K27me3 baixo** e **nenhum pico de ATAC-seq**, qual é o provável estado dele? Que dado da Aula 4 confirmaria a hipótese?
2. Na P2, se não tivéssemos o arquivo de **controle (input)**, que tipo de erro poderíamos cometer?

## 3.2 Por que isso interessa à Engenharia de Defesa? (2:45 – 2:52)

- **Ciência forense e identificação:** o padrão de metilação de certos CpGs muda de forma previsível com a idade. Os chamados **"relógios epigenéticos"** permitem **estimar a idade** de uma pessoa a partir de uma amostra biológica — informação útil quando não há correspondência em bancos de DNA.
- **Biomarcadores de exposição:** exposições ambientais e químicas deixam "cicatrizes" na metilação. O exemplo mais estabelecido é o **tabagismo**, que altera de forma consistente a metilação de genes como *AHRR*. Linha de pesquisa análoga busca marcadores de exposição a poluentes, solventes e radiação — relevante para **saúde ocupacional militar**.
- **Saúde do combatente:** há estudos em coortes militares investigando associações entre metilação e estresse pós-traumático e resiliência.
- **Biodefesa:** muitos patógenos manipulam a cromatina da célula hospedeira para driblar a resposta imune; mapear essas alterações ajuda a entender a infecção. Bactérias também têm padrões próprios de metilação do DNA, detectáveis diretamente por sequenciamento de leitura longa (Nanopore/PacBio), o que pode auxiliar na caracterização de linhagens.
- **Química de proteínas e macromoléculas:** as enzimas escritoras/apagadoras (DNMTs, HDACs, EZH2) são **alvos de fármacos** — tema que retomaremos na **Aula 7** (bioinformática estrutural e docking).

## 3.3 Revisão rápida (2:52 – 2:57)

Responda mentalmente (as respostas estão no final deste documento):

1. O bissulfito converte qual base em uracila: a citosina **metilada** ou a **não metilada**?
2. Qual marca de histona indica um **promotor ativo**: H3K4me3 ou H3K27me3?
3. Para que serve o **controle (input)** num ChIP-seq?
4. Na saída `.narrowPeak`, um valor de −log₁₀(q) igual a **5** é mais ou menos significativo que um valor igual a **2**?
5. Qual técnica usa a enzima **Tn5** para medir acessibilidade?

## 3.4 Tarefas e próxima aula (2:57 – 3:00)

**Para a lista de exercícios (entrega conforme cronograma da disciplina):**

1. **(P1 estendida)** No IGV-Web, escolha um gene relacionado à sua linha de pesquisa e descreva sua assinatura epigenética em GM12878 (H3K4me3, H3K27ac, H3K27me3, cromatina aberta). Inclua uma captura de tela.
2. **(P2 estendida)** Rode novamente o MACS3 com `-q 0.01` e com `--extsize 150`. Quantos picos você obtém em cada caso? Explique por que o número muda.
3. **(P3 estendida)** Refaça a análise exigindo cobertura mínima de 20 leituras. O número de DMCs aumentou ou diminuiu? Por quê?
4. **(Opcional, para quem quiser ir além)** Execute o script do **Apêndice A** (methylKit em R) e compare o número de DMCs com o obtido em Python.

**Ideia para o trabalho final (integrar duas técnicas):** escolha uma proteína do seu interesse, verifique no ENCODE/IGV o estado epigenético do gene que a codifica e combine com a análise estrutural da Aula 7 ou com dados de proteômica da Aula 5.

**Próxima aula (Aula 7):** Bioinformática estrutural — veremos em 3D algumas das proteínas de hoje (histonas, nucleossomo, DNMTs) com PyMOL/ChimeraX, e faremos docking com AutoDock Vina.

---

## Glossário

| Termo | Significado |
|---|---|
| **5mC** | 5-metilcitosina — citosina com grupo metil |
| **Bissulfito** | Reagente que converte C não metilada em U (lida como T) |
| **ChIP-seq** | Imunoprecipitação da cromatina + sequenciamento |
| **ATAC-seq** | Ensaio de acessibilidade da cromatina com transposase Tn5 |
| **Cobertura** | Número de leituras que cobrem uma posição |
| **CpG** | Citosina seguida de guanina no DNA |
| **DMC / DMR** | CpG / região diferencialmente metilada |
| **Enhancer** | Região do DNA que aumenta a transcrição de um gene, às vezes à distância |
| **FDR / valor-q** | Taxa de falsas descobertas; valor-p corrigido para múltiplos testes |
| **Histona** | Proteína em torno da qual o DNA se enrola |
| **Input** | Controle do ChIP-seq, sem anticorpo |
| **Nucleossomo** | Unidade de 8 histonas + ~147 pb de DNA |
| **Pico** | Região com acúmulo de leituras significativamente acima do fundo |
| **Polycomb** | Complexo de proteínas repressoras que deposita H3K27me3 |
| **PTM** | Modificação pós-traducional de proteína |
| **Trilha (track)** | Camada de dados exibida num navegador genômico |
| **Acilação** | Ligação de um grupo acila (derivado de um ácido, como acetil, crotonil ou succinil) a uma lisina |
| **Acil-CoA** | Forma "ativada" de um ácido no metabolismo; é a "tinta" usada pelas enzimas que acilam histonas |
| **Isobárico** | Duas modificações com a mesma massa (ou quase): um desafio para a espectrometria de massas |
| **PRMT** | Enzima que metila argininas (*Protein Arginine MethylTransferase*) |
| **Oncohistona** | Histona com mutação que causa câncer (ex.: H3K27M) |

---

## Apêndice A — Script equivalente em R com methylKit (para casa)

Execute no RStudio (instalado localmente) ou no Colab com ambiente R (*Ambiente de execução → Alterar tipo → R*). A instalação demora; faça-a antes de começar.

```r
# Instalação (uma única vez; pode levar de 10 a 20 minutos)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("methylKit")

library(methylKit)

# Arquivos de exemplo que vêm com o próprio pacote
arquivos <- list(
  system.file("extdata", "test1.myCpG.txt",    package = "methylKit"),
  system.file("extdata", "test2.myCpG.txt",    package = "methylKit"),
  system.file("extdata", "control1.myCpG.txt", package = "methylKit"),
  system.file("extdata", "control2.myCpG.txt", package = "methylKit"))

# Leitura: treatment = 1 (teste) ou 0 (controle)
obj <- methRead(arquivos,
                sample.id = list("teste1", "teste2", "controle1", "controle2"),
                assembly  = "hg18", treatment = c(1, 1, 0, 0), context = "CpG")

getMethylationStats(obj[[1]], plot = TRUE, both.strands = FALSE)  # distribuição
obj_f  <- filterByCoverage(obj, lo.count = 10, hi.perc = 99.9)     # filtro
unidos <- unite(obj_f, destrand = FALSE)                            # CpGs comuns
getCorrelation(unidos, plot = TRUE)                                 # correlação
clusterSamples(unidos, dist = "correlation", method = "ward", plot = TRUE)
PCASamples(unidos)                                                  # PCA

dif   <- calculateDiffMeth(unidos)                                  # teste
dmcs  <- getMethylDiff(dif, difference = 25, qvalue = 0.01)         # todos os DMCs
hiper <- getMethylDiff(dif, difference = 25, qvalue = 0.01, type = "hyper")
hipo  <- getMethylDiff(dif, difference = 25, qvalue = 0.01, type = "hypo")
nrow(dmcs); nrow(hiper); nrow(hipo)
```

Documentação oficial: vinheta do methylKit no Bioconductor (<https://bioconductor.org/packages/methylKit>).

---

## Apêndice B — Referências e documentação

- Allis, C. D.; Jenuwein, T.; Reinberg, D. (Eds.). *Epigenetics*. Cold Spring Harbor Laboratory Press. (bibliografia da disciplina)
- Krueger, F.; Andrews, S. R. Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. *Bioinformatics*, 2011.
- Akalin, A. et al. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biology*, 2012.
- Zhang, Y. et al. Model-based Analysis of ChIP-Seq (MACS). *Genome Biology*, 2008.
- Buenrostro, J. D. et al. Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position (ATAC-seq). *Nature Methods*, 2013.
- ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. *Nature*, 2012.
- Roadmap Epigenomics Consortium. Integrative analysis of 111 reference human epigenomes. *Nature*, 2015.
- Documentação: Bismark (<https://github.com/FelixKrueger/Bismark>), MACS3 (<https://github.com/macs3-project/MACS>), IGV (<https://igv.org>), UCSC Genome Browser (<https://genome.ucsc.edu>), ENCODE (<https://www.encodeproject.org>).

---

## Apêndice C — Catálogo das modificações de histonas

> **Como usar este apêndice:** é material de **consulta**, não de memorização. Ele reúne todos os tipos de modificação de histonas descritos até hoje, organizados por família. Para a aula e a prova, o essencial é o Bloco 3. Use este catálogo quando encontrar uma marca desconhecida num artigo, num banco de dados ou no seu trabalho final.
>
> **Um aviso honesto:** não existe uma lista definitiva. Desde 2011, a espectrometria de massas de alta resolução vem revelando novos tipos quase todo ano. Algumas modificações abaixo estão bem estabelecidas, com enzimas e funções conhecidas; outras foram vistas em poucos estudos. A coluna **"Status"** indica isso.

### C.1 Como ler o nome de qualquer marca

O padrão é sempre **histona + aminoácido + posição + modificação**:

| Exemplo | Histona | Aminoácido e posição | Modificação |
|---|---|---|---|
| `H3K27me3` | H3 | Lisina 27 | Trimetilação |
| `H4R3me2s` | H4 | Arginina 3 | Dimetilação simétrica |
| `H3S10ph` | H3 | Serina 10 | Fosforilação |
| `H2AK119ub1` | H2A | Lisina 119 | Monoubiquitinação |
| `H3Q5ser` | H3 | Glutamina 5 | Serotonilação |
| `γH2AX` | Variante H2A.X | Serina 139 | Fosforilação (nome histórico) |

**Abreviações das modificações:** ac (acetil), me1/me2/me3 (mono/di/trimetil), me2a/me2s (dimetil assimétrica/simétrica), ph (fosfo), ub (ubiquitina), su (SUMO), ar (ADP-ribose), cit (citrulina), pr (propionil), bu (butiril), cr (crotonil), hib (2-hidroxi-isobutiril), bhb (β-hidroxibutiril), mal (malonil), succ (succinil), glu (glutaril), bz (benzoil), la (lactil), fo (formil), og (O-GlcNAc), ser (serotonil), dop (dopaminil), his (histaminil).

### C.2 Mapa por aminoácido: quem pode receber o quê

| Aminoácido | Modificações conhecidas |
|---|---|
| **Lisina (K)** | Acetilação; mono/di/trimetilação; ubiquitinação; SUMOilação; todas as acilações (propionil, butiril, crotonil, 2-hidroxi-isobutiril, β-hidroxibutiril, malonil, succinil, glutaril, benzoil, lactil, formil, isonicotinil); hidroxilação; glicação; biotinilação (controversa) |
| **Arginina (R)** | Monometilação; dimetilação assimétrica e simétrica; citrulinação (deiminação); ADP-ribosilação; glicação |
| **Serina (S)** | Fosforilação; O-GlcNAcilação; ADP-ribosilação; acetilação do grupo amino da ponta da proteína (quando a serina é o 1º aminoácido) |
| **Treonina (T)** | Fosforilação; O-GlcNAcilação |
| **Tirosina (Y)** | Fosforilação; hidroxilação |
| **Histidina (H)** | Fosforilação |
| **Glutamina (Q)** | Metilação; serotonilação; dopaminilação; histaminilação |
| **Glutamato (E) e aspartato (D)** | ADP-ribosilação |
| **Cisteína (C)** | Glutationilação (só a H3 tem cisteína entre as histonas canônicas) |
| **Prolina (P)** | Isomerização cis/trans (muda a forma, não a massa) |
| **Ponta N-terminal (grupo α-amino)** | Acetilação; metilação |

Por que a **lisina** domina a lista? A cadeia lateral dela termina num grupo amino (–NH₃⁺) reativo, com carga positiva, que aceita praticamente qualquer grupo acila ou metil. Neutralizar ou "enfeitar" essa carga muda como a histona segura o DNA e quais proteínas se encaixam nela.

### C.3 As modificações clássicas

As oito famílias mais estudadas, com enzimas e funções bem conhecidas.

#### 1. Acetilação de lisina (Kac)

| | |
|---|---|
| **O que é** | Grupo acetil (–COCH₃), vindo do **acetil-CoA**, ligado à lisina |
| **Efeito químico** | **Neutraliza a carga positiva** da lisina → a histona "segura" o DNA com menos força |
| **Escritores** | Acetiltransferases (**HATs/KATs**): **p300/CBP** (H3K27ac, H3K18ac), **GCN5/PCAF** (H3K9ac), família **MYST** (TIP60 e **MOF**, que faz H4K16ac) |
| **Apagadores** | **HDACs** (dependem de zinco) e **sirtuínas** (SIRT1–7, dependem de NAD⁺) |
| **Leitores** | **Bromodomínios** (ex.: BRD4), domínios **YEATS** e **DPF** |
| **Sítios importantes** | H3K9, K14, K18, **K27**, **K56** (no miolo; marca histonas recém-montadas na replicação); H4K5, K8, K12, **K16** (descompacta a fibra de cromatina) |
| **Função geral** | Quase sempre associada a **cromatina ativa** |
| **Fármacos** | Inibidores de HDAC já aprovados (vorinostate, romidepsina) para linfomas |
| **Δ massa** | **+42,0106 Da** |

#### 2. Metilação de lisina (Kme1, Kme2, Kme3)

| | |
|---|---|
| **O que é** | 1, 2 ou 3 grupos metil (–CH₃), vindos da **S-adenosilmetionina (SAM)** |
| **Efeito químico** | **Não muda a carga**; muda o tamanho e cria um ponto de encaixe para leitores |
| **Escritores** | Metiltransferases com **domínio SET**: MLL/KMT2 (H3K4), SUV39H1, G9a e SETDB1 (H3K9), **EZH2** (H3K27, no complexo PRC2), SETD2 (H3K36me3), NSD1/2 (H3K36me2), SET8 e SUV4-20H (H4K20). Exceção sem domínio SET: **DOT1L** (H3K79) |
| **Apagadores** | **LSD1/KDM1A** (só me1 e me2); desmetilases **JmjC** (removem até me3; dependem de ferro e de **α-cetoglutarato**) |
| **Leitores** | **Cromodomínio** (HP1 lê H3K9me3; CBX do Polycomb lê H3K27me3), **dedos PHD** (H3K4me3), **Tudor** (53BP1 lê H4K20me2), **PWWP** (H3K36me3) |
| **Função** | **Depende do sítio**: H3K4me3, H3K36me3 e H3K79me → ativas; H3K9me3, H3K27me3 e H4K20me3 → repressivas |
| **Fármacos** | Inibidor de EZH2 aprovado (**tazemetostate**); inibidores de DOT1L e LSD1 em estudo |
| **Δ massa** | me1 **+14,0157** · me2 **+28,0313** · me3 **+42,0470 Da** |

> **Ligação com o metabolismo e o câncer:** as desmetilases JmjC (e as enzimas TET do DNA) precisam de α-cetoglutarato. Em tumores com mutação no gene *IDH*, a célula produz um "falso" α-cetoglutarato, o **2-hidroxiglutarato**, que bloqueia essas enzimas. Resultado: histonas e DNA ficam hipermetilados.

#### 3. Metilação de arginina (Rme1, Rme2a, Rme2s)

| | |
|---|---|
| **O que é** | Até **2** grupos metil no grupo guanidino da arginina (nunca 3) |
| **Formas** | **me1** (monometil); **me2a** (2 metis no mesmo nitrogênio, assimétrica); **me2s** (1 metil em cada nitrogênio, simétrica) |
| **Escritores** | **PRMTs**: tipo I (**PRMT1**, **CARM1/PRMT4**) → me2a; tipo II (**PRMT5**) → me2s; tipo III (PRMT7) → só me1 |
| **Apagadores** | **Controversos**: há candidatas (como a JMJD6), mas nenhuma consolidada; a marca parece mais estável que a metilação de lisina |
| **Leitores** | Domínios **Tudor** (ex.: TDRD3 lê me2a) |
| **Exemplos** | H4R3me2a (PRMT1) → ativa; **H4R3me2s** e **H3R8me2s** (PRMT5) → reprimem; H3R17me2a (CARM1) → ativa; H3R2me2 → regula a leitura de H3K4me3 |
| **Fármacos** | Inibidores de PRMT5 em ensaios clínicos para tumores sem o gene *MTAP* |
| **Δ massa** | me1 **+14,0157** · me2a e me2s **+28,0313 Da** (idênticas! ver C.7) |

#### 4. Fosforilação (S, T, Y, H)

| | |
|---|---|
| **O que é** | Grupo fosfato (–PO₃²⁻), vindo do **ATP** |
| **Efeito químico** | Adiciona **carga negativa**: repele o DNA, também negativo |
| **Escritores** | Cinases: **Aurora B** (H3S10ph e H3S28ph, na divisão celular), **Haspin** (H3T3ph), MSK1/2 (H3S10ph na ativação de genes), **ATM, ATR e DNA-PK** (γH2AX), JAK2 (H3Y41ph), WSTF (H2A.X Y142ph) |
| **Apagadores** | Fosfatases (PP1, PP2A; EYA para a tirosina 142 da H2A.X) |
| **Leitores** | Proteínas **14-3-3** (H3S10ph); domínios **BRCT** de MDC1 (γH2AX); survivina (H3T3ph) |
| **Funções** | **Mitose** (condensação dos cromossomos); **resposta a dano no DNA**; ativação rápida de genes |
| **Destaque** | **γH2AX**: em minutos após uma quebra da dupla fita do DNA, milhares de H2A.X ao redor da quebra são fosforiladas e formam um "foco" visível ao microscópio. É usada em **biodosimetria** (estimar a dose de radiação recebida por uma pessoa) |
| **Histidina** | A fosfo-histidina na H4 foi descrita há décadas, mas é instável em meio ácido e se perde nos protocolos comuns; é pouco estudada |
| **Δ massa** | **+79,9663 Da** |

#### 5. Ubiquitinação de lisina (Kub)

| | |
|---|---|
| **O que é** | Ligação de uma **proteína inteira**, a **ubiquitina** (76 aminoácidos, ~8,5 kDa) |
| **Escritores** | Ligases E3: **RING1B** (do complexo **PRC1**, faz H2AK119ub1), **RNF20/RNF40** (H2BK120ub1), **RNF168** (H2AK13/K15ub, no dano ao DNA) |
| **Apagadores** | Desubiquitinases (DUBs): **BAP1** (H2AK119ub1), USP22 (H2Bub, no complexo SAGA) |
| **Leitores** | 53BP1 lê H2AK15ub junto com H4K20me2 no reparo do DNA |
| **Funções** | **H2AK119ub1 → repressão Polycomb**; **H2BK120ub1 → transcrição ativa** e pré-requisito para H3K4me3 e H3K79me |
| **Em espectrometria** | A tripsina corta a ubiquitina e deixa só **2 glicinas** presas à lisina (o "resto GG") |
| **Δ massa** | Resto GG **+114,0429 Da** |

#### 6. SUMOilação de lisina (Ksu)

| | |
|---|---|
| **O que é** | Ligação da **SUMO**, uma proteína "prima" da ubiquitina (~11 kDa) |
| **Enzimas** | UBC9 (a única E2 da via) liga; proteases **SENP** removem |
| **Função** | Em geral **repressiva** (ex.: H4K12su) |
| **Status** | Estabelecida, mas pouco abundante e difícil de medir |
| **Δ massa** | Depende do resto deixado após a digestão (fragmento grande da SUMO) |

#### 7. ADP-ribosilação (S, E, D, R, K)

| | |
|---|---|
| **O que é** | Ligação de **ADP-ribose**, vinda do **NAD⁺**; pode ser uma unidade (mono) ou uma cadeia ramificada (poli-ADP-ribose, **PAR**) |
| **Escritores** | **PARP1/PARP2**, que, junto com o cofator **HPF1**, modificam principalmente **serinas** (H3S10, H3S28, H2BS6) |
| **Apagadores** | **PARG** (corta a cadeia) e **ARH3** (remove a última unidade da serina) |
| **Leitores** | **Macrodomínios** (inclusive da variante macroH2A1.1) e domínios PBZ |
| **Função** | Sinal de **dano ao DNA**: em segundos, a cromatina ao redor da quebra se abre e as proteínas de reparo são recrutadas |
| **Fármacos** | **Inibidores de PARP** (olaparibe) são usados em tumores com mutação em *BRCA1/2* |
| **Δ massa** | **+541,0611 Da** por unidade |

#### 8. Citrulinação ou deiminação (Rcit)

| | |
|---|---|
| **O que é** | A arginina é convertida em **citrulina**, um aminoácido sem carga |
| **Escritor** | **PAD4** (*Peptidyl Arginine Deiminase 4*), dependente de cálcio; age em H3R2, R8, R17, R26 e H4R3 |
| **Apagador** | Nenhum conhecido: é considerada **irreversível** |
| **Função** | Perda maciça de carga positiva → **descondensação da cromatina**; essencial para as **NETs** (redes de DNA que neutrófilos expelem para capturar bactérias) |
| **Doenças** | Artrite reumatoide (o exame **anti-CCP** detecta anticorpos contra proteínas citrulinadas), trombose |
| **Curiosidade** | Já foi relatado que a PAD4 converte arginina **monometilada** em citrulina, o que seria uma via indireta de "apagar" a metilação de arginina |
| **Δ massa** | **+0,9840 Da** |

### C.4 Acilações de lisina: o metabolismo escrevendo na cromatina

Depois da acetilação, descobriu-se que a lisina aceita **muitos outros grupos acila**. Cada um vem de uma molécula diferente do metabolismo (um **acil-CoA**). Quando a célula muda de estado (jejum, hipóxia, dieta, inflamação), a disponibilidade dessas moléculas muda, e as marcas também.

Em geral, a mesma maquinaria da acetilação faz o trabalho:
- **Escritores:** principalmente **p300/CBP**; também KAT2A e TIP60.
- **Apagadores:** **HDAC1–3** e **sirtuínas**. A SIRT5 é especialista em grupos com carga negativa (succinil, malonil, glutaril); a SIRT3 remove crotonil; a SIRT2 remove benzoil.
- **Leitores:** domínios **YEATS** (AF9, ENL, YEATS2), que preferem crotonil à acetil; domínios **DPF**; alguns bromodomínios.

| Modificação | Sigla | Origem metabólica | Exemplo de sítio | Destaque | Status | Δ massa (Da) |
|---|---|---|---|---|---|---|
| **Propionilação** | Kpr | Propionil-CoA (ácidos graxos de cadeia ímpar, alguns aminoácidos) | H3K23pr | Química muito parecida com a acetilação | Estabelecida | +56,0262 |
| **Butirilação** | Kbu | Butiril-CoA; butirato produzido pela microbiota intestinal | H4K5bu, H4K8bu | Liga a dieta e a microbiota à cromatina | Estabelecida | +70,0419 |
| **Crotonilação** | Kcr | Crotonil-CoA (degradação de lisina e triptofano; oxidação de ácidos graxos) | H3K9cr, H3K18cr | Em promotores ativos; lida por domínios YEATS; descoberta junto com dezenas de novos sítios em 2011 | Estabelecida | +68,0262 |
| **2-hidroxi-isobutirilação** | Khib | 2-hidroxi-isobutiril-CoA | H4K8hib | Associada a genes ativos; abundante em células germinativas | Estabelecida | +86,0368 |
| **β-hidroxibutirilação** | Kbhb | **β-hidroxibutirato**, um corpo cetônico | H3K9bhb | **Aumenta no jejum** e na dieta cetogênica; ativa genes do metabolismo de jejum | Estabelecida | +86,0368 |
| **Malonilação** | Kmal | Malonil-CoA (síntese de ácidos graxos) | Vários | Carga negativa; removida pela SIRT5 | Menos estudada | +86,0004 |
| **Succinilação** | Ksucc | **Succinil-CoA** (ciclo de Krebs) | H3K79succ | Troca carga +1 por −1; desestabiliza o nucleossomo | Estabelecida | +100,0160 |
| **Glutarilação** | Kglu | Glutaril-CoA (degradação de lisina e triptofano) | H4K91glu | Carga negativa; afeta a montagem da cromatina | Menos estudada | +114,0317 |
| **Benzoilação** | Kbz | Benzoil-CoA, a partir do **benzoato** (conservante de alimentos) | Vários | Primeira acilação com anel aromático; removida pela SIRT2 | Menos estudada | +104,0262 |
| **Lactilação** | Kla | **Lactato**, produto da glicólise | H3K18la | Descoberta em 2019; muda macrófagos de "inflamatórios" para "reparadores"; abundante em tumores (efeito Warburg) | Estabelecida, em expansão | +72,0211 |
| **Formilação** | Kfo | **Não enzimática**: subproduto de dano oxidativo ao DNA | Vários (inclusive em lisinas que seriam acetiladas ou metiladas) | Um "acidente" químico que pode bloquear marcas normais | Estabelecida como lesão | +27,9949 |
| **Isonicotinilação** | Kinic | Derivada do fármaco **isoniazida** (antituberculose) | Vários | Mostra que **um medicamento** pode virar marca de histona | Recente, poucos estudos | +105,0215 |

### C.5 Outras modificações

| Modificação | Aminoácido e sítio | Enzimas | O que se sabe | Status | Δ massa (Da) |
|---|---|---|---|---|---|
| **O-GlcNAcilação** | S/T; ex.: **H2BS112** | Escritor **OGT**, apagador **OGA** | Sensor de nutrientes (vem da glicose); H2BS112og favorece H2BK120ub. É frágil na fragmentação comum de MS (exige métodos especiais, como ETD) | Estabelecida | +203,0794 |
| **Serotonilação** | Glutamina, **H3Q5ser** | **Transglutaminase 2 (TGM2)** | Descoberta em 2019; quando ao lado de H3K4me3, reforça a ligação do complexo de transcrição TFIID | Estabelecida, recente | +159,0684 |
| **Dopaminilação** | Glutamina, **H3Q5dop** | TGM2 | Descoberta em 2020 em neurônios de recompensa; ligada à recaída no uso de cocaína em modelos animais | Recente | +136,0524 |
| **Histaminilação** | Glutamina, **H3Q5his** | TGM2 | Descrita mais recentemente, em neurônios que produzem histamina; função ainda em estudo | Recente | +94,0531 |
| **Metilação de glutamina** | **H2AQ104me** | **Fibrilarina** | Específica do DNA ribossômico; regula a transcrição pela RNA polimerase I | Estabelecida | +14,0157 |
| **Acetilação N-terminal** | Grupo amino da ponta de H2A e H4 | **NatD (NAA40)** | Diferente da acetilação de lisina; afeta o metabolismo e o envelhecimento em leveduras | Estabelecida | +42,0106 |
| **Metilação N-terminal** | Ponta da variante CENP-A (e de outras proteínas da cromatina) | **NRMT1** | Na CENP-A, ajuda a manter o centrômero | Estabelecida | +14,0157 a +42,0470 |
| **Hidroxilação** | Lisina (5-hidroxilisina) e tirosina | JMJD6 proposta para a lisina | Relatada por MS; função pouco clara. Atenção: tem a mesma massa da **oxidação artificial** que ocorre no preparo da amostra | Controversa | +15,9949 |
| **Isomerização de prolina** | **H3P38** (leveduras) | Isomerase **Fpr4** | A prolina alterna entre as formas *cis* e *trans* e controla a metilação vizinha (H3K36). **Invisível em massa**: estudada por métodos estruturais | Estabelecida, nicho | 0 |
| **Corte da cauda (*clipping*)** | Cauda N-terminal da H3 | Catepsina L, MMP-9 e outras proteases | Remove de uma vez a cauda inteira **e todas as marcas dela**; ocorre na diferenciação celular | Estabelecida | (perda do trecho cortado) |
| **Glutationilação** | Cisteína **H3C110** (e C96 na H3.1) | Não enzimática, por estresse oxidativo | Afeta a estabilidade do nucleossomo; liga o estado redox da célula à cromatina | Recente | +305,0682 |
| **Glicação** | Lisina e arginina | **Não enzimática**, por açúcares e metilglioxal (subproduto da glicólise); a DJ-1 foi proposta como "desglicase" | Altera a arquitetura da cromatina; aumentada em diabetes e câncer | Recente | MG-H1 (R): +54,0106 · CEL (K): +72,0211 · frutosil-lisina: +162,0528 |
| **Biotinilação** | Lisina | Holocarboxilase sintetase (proposta) | Relatada nos anos 2000; estudos posteriores sugerem que é raríssima ou artefato | **Contestada** | +226,0776 |

### C.6 Variantes de histonas e oncohistonas

Não são modificações químicas, mas também mudam a informação epigenética: são **versões alternativas** das histonas, codificadas por outros genes, que substituem as canônicas em certos lugares do genoma.

| Variante | Substitui | Onde / função |
|---|---|---|
| **H3.3** | H3 | Genes ativos e regulados; depositada fora da replicação (pelas chaperonas HIRA e DAXX/ATRX) |
| **CENP-A** | H3 | **Define o centrômero** de cada cromossomo |
| **H2A.Z** | H2A | Promotores e *enhancers*; facilita a ativação |
| **H2A.X** | H2A | Espalhada pelo genoma; sua fosforilação vira **γH2AX** no dano ao DNA |
| **macroH2A** | H2A | Silenciamento (inclusive do cromossomo X inativo); tem macrodomínio que lê ADP-ribose |
| **H1** (histona de ligação) | — | Fica "por fora" do nucleossomo, no DNA de ligação; compacta a fibra. Também recebe acetilação, metilação, fosforilação e outras PTMs |

**Oncohistonas** são histonas com **mutações** que causam câncer. O exemplo clássico é a **H3K27M**: a lisina 27 da H3.3 vira metionina em gliomas pediátricos. A metionina "trava" a enzima EZH2 e derruba a H3K27me3 no genoma inteiro. Outros exemplos: H3K36M (condrossarcoma) e H3G34R/V (gliomas).

### C.7 Conversa entre marcas (*crosstalk*)

As marcas não agem sozinhas. Uma pode facilitar, bloquear ou competir com outra. É isso que se chama de **"código das histonas"**:

- **Competição pelo mesmo sítio:** a H3K27 está **acetilada** (ativa) **ou trimetilada** (reprimida), nunca as duas ao mesmo tempo.
- **Interruptor fosfo-metil:** na mitose, a **H3S10ph** (vizinha da H3K9me3) "empurra" a proteína HP1 para fora, mesmo com a H3K9me3 presente.
- **Dependência entre histonas diferentes:** a **H2BK120ub1** é necessária para que a H3K4 e a H3K79 sejam metiladas.
- **Leitura combinada:** a 53BP1 só se fixa num local de dano ao DNA quando lê **ao mesmo tempo** a H4K20me2 e a H2AK15ub.

### C.8 Tabela-resumo de massas e os pares que confundem a espectrometria (ponte com a Aula 5)

Na proteômica, uma modificação é identificada pela **variação de massa (Δ massa)** que ela causa no peptídeo. Vários pares têm massas **iguais** (isoméricas) ou **quase iguais** (isobáricas) e exigem cuidado:

| Par | Δ massa (Da) | Diferença | Como distinguir |
|---|---|---|---|
| Trimetilação × acetilação | 42,0470 × 42,0106 | 0,036 | Alta resolução (Orbitrap); íon diagnóstico da lisina acetilada (m/z 126,091) |
| Dimetilação × formilação | 28,0313 × 27,9949 | 0,036 | Alta resolução |
| **Arginina me2a × me2s** | 28,0313 × 28,0313 | **0** | Fragmentação: a assimétrica perde **dimetilamina** (45,058); a simétrica perde **metilamina** (31,042) |
| **β-hidroxibutiril × 2-hidroxi-isobutiril** | 86,0368 × 86,0368 | **0** | São isômeros: tempo de retenção na cromatografia e íons diagnósticos |
| Malonil × (Kbhb/Khib) | 86,0004 × 86,0368 | 0,036 | Alta resolução |
| **Lactil × carboxietil (glicação, CEL)** | 72,0211 × 72,0211 | **0** | Isômeros: padrões sintéticos e tempo de retenção |
| Glutaril × resto de ubiquitina (GG) | 114,0317 × 114,0429 | 0,011 | Altíssima resolução; íons de fragmentação |
| Citrulinação × desamidação (de N e Q) | 0,9840 × 0,9840 | **0** | Localizar o resíduo exato (R ou N/Q) nos espectros |
| Hidroxilação × oxidação artificial | 15,9949 × 15,9949 | **0** | Controles de preparo de amostra; localizar o resíduo |

**Outras massas úteis:** fosforilação +79,9663 · crotonil +68,0262 · propionil +56,0262 · butiril +70,0419 · succinil +100,0160 · benzoil +104,0262 · O-GlcNAc +203,0794 · ADP-ribose +541,0611.

### C.9 Como cada modificação é estudada

**1. Anticorpos (ChIP-seq, CUT&RUN, CUT&Tag):** dizem **onde no genoma** está uma marca. Limitações:
- só existe para marcas com anticorpo disponível, e muitas das novas ainda não têm anticorpos confiáveis;
- o anticorpo pode reconhecer marcas parecidas (me2 × me3; me2a × me2s). Antes de usar um anticorpo, vale consultar o **Histone Antibody Specificity Database** (<https://www.histoneantibodies.com>), que testa a especificidade de anticorpos comerciais.

**2. Espectrometria de massas:** diz **quais** marcas existem, **quanto** de cada e **quais combinações** aparecem na mesma molécula. Três estratégias:

| Estratégia | Como funciona | Vantagem | Limitação |
|---|---|---|---|
| **Bottom-up** | Corta as histonas em peptídeos pequenos (com tripsina) | Mais sensível e mais usada | Perde a informação de quais marcas estavam juntas na mesma cauda |
| **Middle-down** | Corta só a cauda inteira (com a enzima Glu-C) | Mostra combinações de marcas na mesma cauda | Análise mais complexa |
| **Top-down** | Analisa a histona inteira, sem cortar | Vê a "proteoforma" completa | Exige instrumentos e softwares avançados |

> **Detalhe de protocolo importante:** histonas têm tantas lisinas e argininas que a tripsina as corta em pedaços pequenos demais. Por isso, o protocolo padrão de bottom-up **propionila quimicamente** as lisinas livres (com anidrido propiônico) antes da digestão. A consequência: para estudar a **propionilação natural**, usa-se um reagente marcado com isótopos pesados, para distinguir a marca biológica da química.

**Ferramentas de software:**
- **EpiProfile:** software livre especializado em identificar e quantificar PTMs de histonas a partir de dados de bottom-up.
- **MaxQuant** (Aula 5): busca qualquer PTM listada como "modificação variável".
- **Skyline** (Aula 5): quantificação direcionada de peptídeos de histonas modificados.

**Bancos de dados:**
- **HIstome**: catálogo de modificações de histonas humanas e das enzimas que as escrevem e apagam.
- **HistoneDB**, do NCBI: sequências de todas as variantes de histonas.
- **UniProt** e **PhosphoSitePlus**: modificações anotadas sítio a sítio.

### C.10 Referências do catálogo

- Strahl, B. D.; Allis, C. D. The language of covalent histone modifications. *Nature*, 2000. (proposta do "código das histonas")
- Rogakou, E. P. et al. DNA double-stranded breaks induce histone H2AX phosphorylation on serine 139. *Journal of Biological Chemistry*, 1998. (γH2AX)
- Tan, M. et al. Identification of 67 histone marks and histone lysine crotonylation as a new type of histone modification. *Cell*, 2011.
- Tessarz, P. et al. Glutamine methylation in histone H2A is an RNA-polymerase-I-dedicated modification. *Nature*, 2014.
- Zhao, Y.; Garcia, B. A. Comprehensive catalog of currently documented histone modifications. *Cold Spring Harbor Perspectives in Biology*, 2015.
- Xie, Z. et al. Metabolic regulation of gene expression by histone lysine β-hydroxybutyrylation. *Molecular Cell*, 2016.
- Zhang, D. et al. Metabolic regulation of gene expression by histone lactylation. *Nature*, 2019.
- Farrelly, L. A. et al. Histone serotonylation is a permissive modification that enhances TFIID binding to H3K4me3. *Nature*, 2019.
- Lepack, A. E. et al. Dopaminylation of histone H3 in ventral tegmental area regulates cocaine seeking. *Science*, 2020.

---

### Respostas da revisão rápida

1. A **não metilada** (a metilada fica protegida e continua sendo lida como C).
2. **H3K4me3**.
3. Medir o **ruído de fundo**, para que só regiões realmente enriquecidas pelo anticorpo sejam chamadas como picos.
4. **Mais** significativo (5 corresponde a q = 0,00001; 2 corresponde a q = 0,01).
5. **ATAC-seq**.
