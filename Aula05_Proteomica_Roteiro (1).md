# Aula 5 — Proteômica: identificação e quantificação de proteínas por espectrometria de massas

**Disciplina:** Introdução à Bioinformática — PPGED / IME
**Duração:** 3 horas (online, síncrona)
**Formato:** teoria (65 min) → intervalo (10 min) → prática (70 min) → conexão teoria–prática e fechamento (30 min)

---

## Antes da aula (faça com antecedência — 15 minutos)

A prática acontece toda no navegador. Você **não** precisa instalar nada nem saber programar. Mas precisa de três coisas prontas:

1. **Conta no Galaxy Europe** — acesse <https://usegalaxy.eu>, clique em *Login or Register* e crie uma conta gratuita com seu e-mail. Confirme o e-mail recebido. *(Sem conta confirmada, a ferramenta estatística da prática não roda.)*
2. **Navegador atualizado** (Chrome, Firefox ou Edge) e, de preferência, **dois monitores ou a tela dividida**: de um lado a aula, do outro o Galaxy.
3. **Leitura rápida (opcional, mas recomendada):** o resumo do artigo de Aebersold & Mann, *Mass-spectrometry-based proteomics*, Nature (2003), indicado na bibliografia da disciplina.

> **Material complementar — caderno em Python (Google Colab):** a prática inteira também está disponível em código Python comentado, para quem quiser ver "por dentro" o que as ferramentas fazem.
> - **Abrir direto no Colab:** <https://colab.research.google.com/drive/1guMQyhhtXbJcfGYorPtBkc58t9UBWwJJ>
> - **Arquivo no Google Drive (para baixar o .ipynb):** <https://drive.google.com/file/d/1guMQyhhtXbJcfGYorPtBkc58t9UBWwJJ/view?usp=sharing>
>
> Ao abrir, use **Arquivo → Salvar uma cópia no Drive** para ter a sua própria versão editável. Não é preciso instalar nada: o Colab roda no navegador, com uma conta Google.

---

## Objetivos da aula

Ao final, você será capaz de:

- Explicar **por que** estudar proteínas diretamente, e não só genes;
- Descrever, etapa por etapa, o **fluxo "bottom-up"** de um experimento de proteômica por espectrometria de massas;
- Entender **como um computador identifica** um peptídeo a partir de um espectro, e o que significa **FDR de 1%**;
- Distinguir as principais estratégias de **quantificação** (sem marcação, SILAC, TMT, direcionada);
- Explicar o que é o **interatoma** e como uma lista de proteínas alteradas vira uma **rede de interações** (STRING);
- Saber **para que serve cada ferramenta** (MaxQuant/Andromeda, Perseus, PatternLab, Skyline, MSstats, Galaxy, UniProt, PRIDE, STRING);
- **Executar** uma análise estatística real com dados públicos e **interpretar** um *volcano plot*.

---

## Mapa da aula

| Horário | Bloco | O que acontece |
|---|---|---|
| 0:00 – 0:05 | Abertura | Objetivos e o "fio condutor" da aula |
| 0:05 – 1:10 | **Momento 1 — Teoria** | Conceitos, técnicas e ferramentas (6 blocos curtos) |
| 1:10 – 1:20 | **Intervalo** | 10 minutos |
| 1:20 – 2:30 | **Momento 2 — Prática** | Parte A (navegador): da proteína ao peptídeo. Parte B (Galaxy): estatística em dados reais de câncer de pele |
| 2:30 – 3:00 | **Momento 3 — Conexão e fechamento** | Ligando cada clique a um conceito, aplicações em defesa, conexão com o trabalho final |

> **Fio condutor da aula:** *"Como sair de um tubo com milhares de proteínas misturadas e chegar a uma lista confiável de quais proteínas estão lá e quais mudaram de quantidade?"* Guarde essa pergunta: cada conceito da teoria responde a um pedaço dela, e a prática percorre o caminho inteiro.

---

# MOMENTO 1 — Teoria (0:05 – 1:10)

## Bloco 1 — Por que proteômica? (10 min)

### Do genoma ao proteoma

Nas aulas anteriores trabalhamos com DNA e RNA. Lembre do **dogma central**: DNA → RNA → proteína. Uma analogia útil para engenheiros:

- O **genoma** é o *conjunto completo de plantas de projeto* de uma fábrica. É praticamente igual em todas as células do corpo.
- O **transcriptoma** (RNA) é o conjunto de *ordens de produção emitidas* naquele momento.
- O **proteoma** é o conjunto de *máquinas que efetivamente estão no chão de fábrica, funcionando*.

Por que não basta olhar para o RNA? Porque **a quantidade de RNA nem sempre prevê a quantidade de proteína**. Entre a ordem de produção e a máquina pronta há muitas etapas: velocidade de tradução, tempo de vida da proteína, degradação. Além disso, as proteínas sofrem **modificações pós-traducionais (PTMs)** — fosforilação, glicosilação, acetilação, entre outras — que funcionam como "interruptores" e não aparecem no DNA nem no RNA.

### Definições essenciais

- **Proteoma:** conjunto de todas as proteínas presentes em uma célula, tecido ou fluido, em um dado momento e condição. É **dinâmico**: muda com o tempo, com o estresse, com uma doença.
- **Proteômica:** o estudo em larga escala do proteoma — *quais* proteínas estão presentes, *quanto* de cada uma, *como* estão modificadas e *com quem* interagem.

### Por que isso interessa à Engenharia de Defesa?

- **Detecção e identificação de toxinas proteicas** (por exemplo, ricina e toxinas botulínicas) diretamente por espectrometria de massas, com alta especificidade;
- **Biomarcadores de exposição** a agentes químicos: alguns agentes neurotóxicos se ligam covalentemente a proteínas do sangue (como a butirilcolinesterase e a albumina), deixando um "registro" detectável;
- **Identificação de microrganismos** (biodetecção e microbiologia forense) pelo perfil de proteínas;
- **Biomateriais e proteínas de interesse tecnológico** (adesivos, fibras, enzimas de descontaminação).

---

## Bloco 2 — A espectrometria de massas: a "balança de moléculas" (11 min)

### A ideia central

Um **espectrômetro de massas** é, essencialmente, uma **balança de altíssima precisão para moléculas**. Mas há um detalhe: ele só consegue "pesar" moléculas que tenham **carga elétrica** e estejam em **fase gasosa**. Por isso, ele não mede a massa diretamente, e sim a razão **massa/carga**, escrita **m/z**.

> **Exemplo:** um peptídeo de massa 1148,6 Da com uma carga positiva (+1) aparece em m/z ≈ 1149,6. Se ele tiver duas cargas (+2), aparece em m/z ≈ 575,3 — praticamente metade. (Cada carga vem de um próton, de massa ≈ 1,007 Da.)

*Da (Dalton) é a unidade de massa molecular: 1 Da ≈ massa de um átomo de hidrogênio.*

### As três partes de todo espectrômetro

| Parte | Função | Exemplos | Analogia |
|---|---|---|---|
| **Fonte de ionização** | Transforma as moléculas em íons em fase gasosa, sem destruí-las | **ESI** (*electrospray*: o líquido é borrifado sob alta tensão) e **MALDI** (um laser atinge a amostra misturada a uma matriz) | A "porta de entrada" que dá carga às moléculas |
| **Analisador de massas** | Separa os íons pela razão m/z | **Quadrupolo** (filtro por campos elétricos oscilantes), **TOF** (*time-of-flight*: íons mais leves chegam primeiro, como numa corrida), **Orbitrap** (íons orbitam um eletrodo e a frequência da órbita revela o m/z) | A "pista de corrida" ou o "filtro" |
| **Detector** | Conta os íons que chegam em cada m/z | Multiplicadores de elétrons, detecção de corrente induzida | O "placar" |

O resultado de uma medida é um **espectro de massas**: um gráfico com **m/z no eixo horizontal** e **intensidade** (quantos íons chegaram) no eixo vertical. Cada pico é uma espécie molecular.

### Duas palavras que voltarão o tempo todo

- **MS1 (ou MS):** o espectro dos peptídeos **inteiros**. Diz *quais massas* estão chegando e *com que intensidade* (útil para **quantificar**).
- **MS2 (ou MS/MS, espectrometria em tandem):** o instrumento **escolhe** um íon do MS1, **quebra** esse íon em pedaços (por colisão com um gás inerte — técnicas chamadas **CID** e **HCD**) e mede os pedaços. Diz *qual é a sequência* do peptídeo (útil para **identificar**).

---

## Bloco 3 — O fluxo "bottom-up" (shotgun): da amostra ao espectro (12 min)

Proteínas inteiras são grandes, difíceis de ionizar e de fragmentar de forma interpretável. A estratégia dominante é **quebrá-las em peptídeos antes de medir** e depois "remontar" a informação no computador. Por isso se chama **bottom-up** ("de baixo para cima"). Também é chamada de **shotgun** ("espingarda"), porque se analisa a mistura inteira de uma vez.

### As etapas

1. **Extração das proteínas** da amostra (células, tecido, plasma, uma cultura bacteriana).
2. **Digestão enzimática com tripsina.** A tripsina é uma enzima que funciona como uma **tesoura com regra fixa**: corta logo **depois** dos aminoácidos **lisina (K)** e **arginina (R)**, **exceto** quando o próximo aminoácido é **prolina (P)**. Como a regra é previsível, o computador consegue simular o corte (é a **digestão *in silico***). Bônus: todo peptídeo tríptico termina em K ou R, que têm carga positiva — ótimo para a ionização.
3. **Separação por cromatografia líquida (LC).** A mistura de dezenas de milhares de peptídeos passa por uma coluna fina. Peptídeos mais hidrofóbicos ficam retidos por mais tempo. Assim eles **saem da coluna aos poucos**, ao longo de 1–2 horas. O momento em que cada peptídeo sai é o **tempo de retenção**. O conjunto LC + espectrômetro é chamado de **LC-MS/MS**.
4. **Ionização (ESI)** à medida que os peptídeos saem da coluna.
5. **Aquisição de dados**, alternando MS1 e MS2. Há duas estratégias principais:
   - **DDA** (*data-dependent acquisition*): a cada ciclo, o instrumento escolhe os ~10–20 picos **mais intensos** do MS1 e fragmenta cada um separadamente. Simples de interpretar, mas peptídeos pouco abundantes podem ser "esquecidos" em algumas amostras — o que gera **valores ausentes**.
   - **DIA** (*data-independent acquisition*): o instrumento fragmenta **todos os íons** dentro de janelas largas de m/z, sistematicamente. Mais completo e reprodutível, porém os espectros ficam misturados e exigem software mais sofisticado.
6. **Arquivos brutos** (*raw files*) — cada corrida gera um arquivo de centenas de MB a alguns GB, em formato do fabricante (por exemplo `.raw` da Thermo). O formato aberto padrão para troca é o **mzML**.

> **Para contraste — top-down:** analisa as proteínas **inteiras**, sem digestão. Preserva combinações de modificações na mesma molécula (as chamadas *proteoformas*), mas é tecnicamente mais difícil. Nesta aula ficamos no bottom-up, que é o padrão.

---

## Bloco 4 — Identificação: como o computador descobre "quem é quem" (11 min)

### O espectro MS2 é uma escada

Quando um peptídeo é fragmentado, ele quebra preferencialmente nas **ligações peptídicas**. Cada quebra gera dois pedaços:

- **íons b**: contêm o início (lado N-terminal) do peptídeo;
- **íons y**: contêm o final (lado C-terminal).

A **diferença de massa entre dois íons b consecutivos** (ou dois y consecutivos) é exatamente a **massa de um aminoácido**. Ou seja: o espectro é uma **escada**, e o tamanho de cada degrau revela qual aminoácido está ali.

**Exemplo real — peptídeo `LVNEVTEFAK` da albumina humana** (massas para íons de carga +1, em Da):

| Quebra | Fragmento b | m/z (b) | Fragmento y | m/z (y) |
|---|---|---|---|---|
| 1 | L | 114,09 | VNEVTEFAK | 1036,53 |
| 2 | LV | 213,16 | NEVTEFAK | 937,46 |
| 3 | LVN | 327,20 | EVTEFAK | 823,42 |
| 4 | LVNE | 456,25 | VTEFAK | 694,38 |
| 5 | LVNEV | 555,31 | TEFAK | 595,31 |
| 6 | LVNEVT | 656,36 | EFAK | 494,26 |
| 7 | LVNEVTE | 785,40 | FAK | 365,22 |
| 8 | LVNEVTEF | 932,47 | AK | 218,15 |
| 9 | LVNEVTEFA | 1003,51 | K | 147,11 |

Repare: de `L` (114,09) para `LV` (213,16), a diferença é 99,07 — a massa do resíduo de **valina (V)**. Este peptídeo, de massa monoisotópica 1148,61 Da, é observado com frequência em análises de plasma e soro.

### Busca em banco de dados (*database search*) — o método padrão

Ler a escada "à mão" (sequenciamento *de novo*) é possível, mas lento e sujeito a erros. O método padrão é **comparar com o que já se sabe**:

1. Pega-se um arquivo **FASTA** com as sequências de **todas as proteínas conhecidas** do organismo (por exemplo, o proteoma humano baixado do **UniProt**). *FASTA é um formato de texto simples: uma linha de cabeçalho começando com `>` e, abaixo, a sequência em letras.*
2. O programa faz a **digestão *in silico*** de todas essas proteínas com a regra da tripsina, gerando milhões de peptídeos teóricos.
3. Para cada espectro experimental, o programa seleciona os peptídeos teóricos com a **mesma massa** (dentro de uma tolerância de poucos ppm — partes por milhão), calcula seus espectros teóricos e **pontua** a semelhança.
4. O melhor par espectro–peptídeo é chamado de **PSM** (*peptide-spectrum match*).

O programa que faz essa comparação é chamado de **motor de busca** (*search engine*). Exemplos: **Andromeda** (dentro do MaxQuant), **Comet** (usado pelo PatternLab), **MSFragger** (dentro do FragPipe), **Mascot** (comercial).

### Como saber se não estamos nos enganando? A estratégia *target-decoy* e o FDR

Com milhões de comparações, alguns acertos acontecem **por puro acaso**. A solução elegante:

- Junto com o banco verdadeiro (**target**), o programa cria um banco **falso** (**decoy**), geralmente **invertendo** cada sequência de proteína. Nenhum peptídeo real deveria combinar com ele.
- Toda vez que um espectro "acerta" um peptídeo decoy, sabemos que aquilo é um **falso positivo**. Contando quantos decoys passam por um limiar de pontuação, **estimamos** quantos falsos positivos existem entre os targets.
- Ajusta-se o limiar para que a **FDR** (*false discovery rate*, taxa de falsas descobertas) fique em **1%**: dos peptídeos/proteínas aceitos, espera-se que no máximo 1 em cada 100 seja um erro.

> **Analogia para a área de defesa:** é como testar um radar inserindo **alvos fictícios conhecidos** no ambiente. Contando quantas vezes o radar "trava" nos alvos fictícios, você estima sua taxa de falsos alarmes e calibra a sensibilidade.

### Inferência de proteínas: o problema do quebra-cabeça

Medimos **peptídeos**, mas queremos falar de **proteínas**. Problema: um mesmo peptídeo pode existir em várias proteínas parecidas (isoformas, famílias de genes). Por isso os programas reportam **grupos de proteínas** (*protein groups*): conjuntos de proteínas que não podem ser distinguidas pelos peptídeos observados. Termos que você vai ver:

- **Peptídeo único** (*unique*): aparece em apenas uma proteína/grupo — é a evidência mais forte;
- **Peptídeo *razor***: compartilhado, mas atribuído (pelo princípio da parcimônia, a "navalha de Occam") ao grupo que já tem mais evidências.

---

## Bloco 5 — Quantificação: "quanto" de cada proteína (9 min)

Saber que uma proteína **está lá** é metade da história. Em geral queremos saber se ela **aumentou ou diminuiu** entre condições (doente × saudável, exposto × controle).

| Estratégia | Como funciona | Vantagens | Limitações |
|---|---|---|---|
| **Sem marcação (*label-free*, LFQ)** | Cada amostra é medida em uma corrida separada. Compara-se a **intensidade do sinal MS1** de cada peptídeo entre corridas (ou, de forma mais simples, o **número de espectros**, o *spectral counting*). O algoritmo **MaxLFQ** do MaxQuant é um padrão de referência. | Barato, sem limite de número de amostras | Mais sensível à variação entre corridas; mais valores ausentes |
| **Marcação metabólica (SILAC)** | Células crescem com aminoácidos "pesados" (isótopos estáveis ¹³C, ¹⁵N). As amostras são misturadas: o mesmo peptídeo aparece como par leve/pesado. | Muito preciso | Só para células em cultura; poucas condições (2–3) |
| **Marcação isobárica (TMT, iTRAQ)** | Peptídeos de cada amostra recebem uma "etiqueta" química de mesma massa total. Ao fragmentar, cada etiqueta libera um **íon repórter** de massa diferente. | Até 16–18 amostras numa só corrida | Custo dos reagentes; distorção de razões |
| **Direcionada (SRM / PRM)** | O instrumento é programado para monitorar **apenas uma lista de peptídeos pré-escolhidos**. | Sensibilidade e exatidão máximas — ideal para **validar** biomarcadores e para **detecção de alvos específicos** (como uma toxina) | Só mede o que você já sabe que procura |

### O que acontece com os números depois (a parte estatística)

1. **Transformação logarítmica (log₂):** intensidades variam por várias ordens de grandeza. No log₂, **dobrar** a quantidade equivale a somar **1**. Um *log₂ fold change* de +1 = o dobro; −1 = a metade; +0,58 ≈ 1,5 vez.
2. **Normalização:** corrige diferenças técnicas entre corridas (por exemplo, uma amostra com um pouco mais de material injetado). Abordagem comum: igualar as medianas de todas as corridas.
3. **Valores ausentes:** em DDA, muitos peptídeos não são medidos em todas as amostras — frequentemente porque estavam **abaixo do limite de detecção** (ausência "não aleatória"). Eles podem ser **imputados** (preenchidos com valores estimados baixos).
4. **Teste estatístico por proteína:** teste t ou modelos lineares, gerando um **valor-p** para cada proteína.
5. **Correção para testes múltiplos:** testando 2.000 proteínas com p < 0,05, esperaríamos ~100 "falsos achados" só por acaso. Por isso usamos o **valor-p ajustado** (método de Benjamini–Hochberg), que controla a FDR.
6. **Volcano plot:** gráfico que resume tudo — **eixo x = log₂ fold change** (tamanho do efeito) e **eixo y = −log(valor-p ajustado)** (confiança). Proteínas relevantes ficam nos "braços" superiores: mudaram muito **e** com confiança.

---

## Bloco 6 — As ferramentas, os bancos de dados e o interatoma (12 min)

Para cada ferramenta: **o que é**, **o que faz no fluxo**, **licença/custo** e **onde roda**.

### Identificação e quantificação (do arquivo bruto à tabela)

**MaxQuant + Andromeda**
- *O que é:* a plataforma mais usada no mundo para análise de dados de proteômica por espectrometria de massas, desenvolvida pelo grupo de Jürgen Cox (Instituto Max Planck de Bioquímica, Alemanha). O **Andromeda** é o motor de busca embutido nele.
- *O que faz:* recebe os arquivos brutos + o FASTA e devolve tabelas prontas: `evidence.txt` (cada peptídeo medido), `peptides.txt` e `proteinGroups.txt` (cada grupo de proteínas, com intensidades e valores LFQ). Faz identificação, controle de FDR, quantificação (LFQ, SILAC, TMT) e o *match between runs* (transfere identificações entre corridas pelo tempo de retenção, reduzindo valores ausentes).
- *Licença:* **gratuito**, mas **não é de código aberto** (exige cadastro para download).
- *Onde roda:* Windows e Linux. Também está disponível **dentro do Galaxy**, sem instalação.

**PatternLab for Proteomics** — *ferramenta nacional*
- *O que é:* um ambiente computacional completo para proteômica *shotgun*, desenvolvido no Brasil pelo grupo de Paulo C. Carvalho (**Fiocruz Paraná**), com raízes no doutorado do autor na COPPE/UFRJ. A versão atual, **PatternLab V**, foi publicada na *Nature Protocols* (2022).
- *O que faz:* cobre o fluxo todo — preparação do banco FASTA, busca (motor **Comet**), filtragem estatística dos resultados, quantificação (sem marcação e com marcação), análise diferencial e interpretação biológica (Gene Ontology), com interface gráfica.
- *Licença:* **gratuito** (<http://patternlabforproteomics.org>). Roda em Windows.

**Skyline**
- *O que é:* software do laboratório MacCoss (Universidade de Washington, EUA), referência mundial em **proteômica direcionada**.
- *O que faz:* ajuda a **planejar** quais peptídeos e fragmentos monitorar (SRM/PRM), e a **visualizar e quantificar** os cromatogramas de cada alvo. Também é muito usado para DIA. É a ferramenta típica para construir um **ensaio de detecção** de uma proteína-alvo.
- *Licença:* **gratuito e de código aberto**. Roda em Windows.

### Estatística e visualização (da tabela à resposta biológica)

**Perseus**
- *O que é:* o "parceiro" do MaxQuant para análise estatística, do mesmo grupo (Cox).
- *O que faz:* abre o `proteinGroups.txt` e, por **menus e cliques**, permite filtrar contaminantes, transformar em log₂, imputar valores ausentes, fazer testes t, *volcano plots*, *heatmaps* e agrupamentos.
- *Licença:* gratuito, não é de código aberto. Roda em Windows.

**MSstats** — *usado na prática de hoje*
- *O que é:* pacote estatístico de **código aberto** para a linguagem R (projeto Bioconductor), do grupo de Olga Vitek (Northeastern University, EUA).
- *O que faz:* recebe a saída do MaxQuant (ou de Skyline, DIA-NN, etc.) e uma descrição do desenho experimental, e realiza: log₂, normalização, imputação, resumo dos peptídeos em uma quantidade por proteína, e comparação entre grupos usando **modelos lineares de efeitos mistos**, que levam em conta réplicas biológicas e técnicas. Gera tabela de resultados e *volcano plot*.

**Galaxy** — *usado na prática de hoje*
- *O que é:* plataforma web de **código aberto** para análises de bioinformática **sem programação**. Servidores públicos e gratuitos, como o **usegalaxy.eu** (Europa), executam as ferramentas em computadores remotos.
- *O que faz:* oferece centenas de ferramentas (inclusive MaxQuant e MSstats) por formulários. Cada passo fica registrado num **histórico**, o que torna a análise **reprodutível** e compartilhável.
- *Material de apoio:* a **Galaxy Training Network** (<https://training.galaxyproject.org>) mantém tutoriais gratuitos de proteômica. A prática de hoje é uma versão guiada e encurtada de um deles.

**Expasy PeptideMass** — *usado na prática de hoje*
- *O que é:* ferramenta web do **SIB – Instituto Suíço de Bioinformática** (<https://web.expasy.org/peptide_mass/>).
- *O que faz:* "corta" uma proteína com a enzima escolhida e calcula a massa de cada peptídeo gerado — exatamente a digestão *in silico* que os motores de busca fazem por dentro.

### Bancos de dados

**UniProt** (<https://www.uniprot.org>)
- Banco de dados central de **sequências e funções de proteínas**. Duas seções: **Swiss-Prot** (revisada manualmente por curadores, mais confiável) e **TrEMBL** (anotada automaticamente). Cada proteína tem um código de acesso — por exemplo, **P02768** é a albumina sérica humana. É de onde se baixa o **FASTA** para a busca.

**ProteomeXchange e PRIDE**
- **ProteomeXchange** (<https://www.proteomexchange.org>) é o **consórcio** que padroniza a submissão pública de dados de proteômica. Cada conjunto recebe um identificador **PXD** (por exemplo, **PXD006914**, o dado de hoje). Muitas revistas exigem esse depósito.
- **PRIDE** (<https://www.ebi.ac.uk/pride>), no Instituto Europeu de Bioinformática (EMBL-EBI), é o maior repositório membro do consórcio: guarda os **arquivos brutos**, os resultados e os metadados.

### Interatoma: *com quem* as proteínas interagem (4 min)

Proteínas raramente trabalham sozinhas. Elas se juntam em **complexos** — o ribossomo, por exemplo, reúne dezenas de proteínas — e passam sinais umas às outras em **cadeias**, como numa rede de comunicações. O conjunto de **todas as interações entre as proteínas** de uma célula é o **interatoma**.

> **Analogia:** se o proteoma é a *lista de funcionários* da fábrica, o interatoma é o *organograma* — quem trabalha em equipe com quem, e quem passa ordens para quem.

**Como se mapeia o interatoma — e a espectrometria de massas está no centro.** Três das quatro técnicas principais terminam exatamente no fluxo *bottom-up* que acabamos de estudar:

| Técnica | Como funciona | Analogia |
|---|---|---|
| **AP-MS** (purificação por afinidade + espectrometria de massas) | Uma proteína de interesse (a **isca**) recebe uma "etiqueta" e é "pescada" da célula com uma resina ou anticorpo. Tudo o que vem preso a ela é digerido e identificado por LC-MS/MS | Pescar um peixe e ver o que veio preso no anzol |
| **Marcação por proximidade** (BioID, TurboID, APEX) | A isca carrega uma enzima que "carimba" quimicamente (com biotina) as proteínas **vizinhas**, dentro da célula viva. As carimbadas são purificadas e identificadas por MS | Um spray de tinta: quem estava perto sai marcado |
| **XL-MS** (ligação cruzada + MS) | Um reagente químico "grampeia" aminoácidos de proteínas que **se tocam**. O MS identifica os pares grampeados e revela até a região do contato | Grampear duas folhas que estavam encostadas |
| **Duplo-híbrido em levedura** (Y2H) | Teste genético, par a par, **sem** espectrometria: se duas proteínas se ligam, a levedura liga um gene repórter | Um interruptor que só liga se as duas peças se encaixam |

*Escala:* o projeto **BioPlex** (Harvard), usando AP-MS com milhares de iscas, mapeou cerca de **120 mil interações** entre ~14 mil proteínas humanas.

**Bancos de dados de interações** (todos gratuitos):

**STRING** (<https://string-db.org>)
- *O que é:* o banco de interações proteína–proteína mais usado, mantido por um consórcio europeu (SIB, Universidade de Copenhague, EMBL). Cobre milhares de organismos e tem site e API.
- *De onde vêm as interações:* experimentos publicados (AP-MS, Y2H…), bancos curados de complexos e vias, **mineração de texto** da literatura (lembre da Unidade 1), coexpressão e comparação entre genomas.
- *Escore de confiança (0 a 1):* combina essas evidências — 0,4 = média; 0,7 = alta; 0,9 = altíssima.
- *Atenção:* por padrão, uma linha no STRING significa **associação funcional** ("trabalham juntas"), **não necessariamente contato físico**. Há uma opção para ver só a rede física.
- *Uso típico:* colar uma lista de proteínas (por exemplo, as significativas do *volcano plot*) e ver a rede que elas formam, os *hubs* (proteínas com muitas conexões) e as funções enriquecidas.

**BioGRID** e **IntAct**
- Bancos **curados manualmente**: especialistas leem artigos e registram cada interação com o método experimental que a demonstrou. O IntAct é mantido no EMBL-EBI (o mesmo instituto do PRIDE). São mais restritos que o STRING, porém cada interação tem uma evidência experimental direta.

> **Coexpressão ≠ interação.** Na transcriptômica aparecem **redes de coexpressão**: genes cujo RNA sobe e desce junto. Isso sugere uma relação funcional, mas **não prova** que as proteínas se tocam. Interatoma, no sentido estrito, é contato entre proteínas.

**E para a defesa?** Mapear o **interatoma patógeno–hospedeiro** — quais proteínas humanas uma toxina ou uma proteína viral "sequestra" — aponta alvos para contramedidas (antídotos, vacinas, fármacos). E conhecer os parceiros de uma proteína ajuda a escolher alvos e interfaces para o *docking* da Unidade 7.

> **Resumo em uma frase:** *UniProt fornece as sequências; o espectrômetro gera os espectros; MaxQuant/Andromeda ou PatternLab transformam espectros em listas de proteínas com FDR controlada; Perseus ou MSstats dizem o que mudou; o STRING mostra com quem essas proteínas interagem; Skyline valida alvos específicos; PRIDE guarda tudo para que outros possam reproduzir.*

---

# INTERVALO (1:10 – 1:20)

Dez minutos. **Ao voltar, já esteja logado no <https://usegalaxy.eu>** e com uma aba aberta em <https://web.expasy.org/peptide_mass/>.

---

# MOMENTO 2 — Prática (1:20 – 2:30)

## O estudo de caso

Vamos usar dados reais e públicos do estudo de **Föll et al. (2018), *Matrix Biology*** (depósito **PRIDE PXD006914**), na versão preparada para ensino pela Galaxy Training Network.

- **Pergunta biológica:** o carcinoma espinocelular cutâneo (cSCC) é um câncer de pele comum, geralmente causado por radiação ultravioleta, quase sempre curável — mas os casos que **metastatizam** são muito graves. Pacientes com a doença genética **epidermólise bolhosa distrófica recessiva (RDEB)** — em que falta o colágeno VII e a pele forma bolhas facilmente — desenvolvem um cSCC muito agressivo. *O que diferencia, em nível de proteínas, esses dois tumores agressivos de origens diferentes?*
- **Amostras:** tecido tumoral de **19 pacientes** — **13** com cSCC metastático e **6** com cSCC associado à RDEB — preservado em formol e parafina (como em arquivos de patologia hospitalar).
- **Técnica:** proteômica *bottom-up*, **sem marcação (LFQ)**, aquisição **DDA**, processamento com **MaxQuant**.

Rodar o MaxQuant nos 19 arquivos brutos levaria **horas**. Por isso, a Galaxy Training Network disponibiliza os **resultados do MaxQuant já prontos** no repositório Zenodo — e é deles que vamos partir. Na Parte A, porém, você vai reproduzir à mão a lógica que o MaxQuant executa por dentro.

---

## Parte A — Da proteína ao peptídeo: pensando como o motor de busca (20 min)

**Ferramentas:** UniProt e Expasy PeptideMass (navegador).

### Passo A1 — Encontre a proteína no UniProt (5 min)

1. Acesse <https://www.uniprot.org> e busque por **P02768**.
2. Abra a entrada **ALBU_HUMAN — Albumin**. A albumina é a proteína mais abundante do plasma sanguíneo.
3. Role até a seção **Sequence**. Observe o tamanho (em aminoácidos) e a massa (em Da).
4. Ainda na entrada, procure a seção **PTM/Processing**: veja que a proteína possui um peptídeo-sinal e um propeptídeo que são removidos — a proteína madura é menor que a sequência "do gene". *(Mais um motivo para medir proteínas diretamente!)*

### Passo A2 — Digestão *in silico* no PeptideMass (10 min)

1. Acesse <https://web.expasy.org/peptide_mass/>.
2. No campo da sequência, digite o código **P02768** (a ferramenta busca a sequência sozinha).
3. Configure:
   - **Enzyme:** *Trypsin*
   - **Allow for … missed cleavages:** *0* (supõe que a tripsina cortou em todos os sítios possíveis)
   - **Masses:** *monoisotopic*, como **[M+H]⁺**
   - **Display peptides with a mass bigger than:** *500* Dalton
4. Clique em **Perform** e examine a tabela de peptídeos.

**Perguntas para responder (anote):**

- **A.1** Quantos peptídeos foram gerados com massa acima de 500 Da? Por que faz sentido descartar os peptídeos muito pequenos? *(Pista: um peptídeo de 2 ou 3 aminoácidos existe em centenas de proteínas diferentes.)*
- **A.2** Encontre o peptídeo **LVNEVTEFAK** na lista. Qual a massa [M+H]⁺ informada? Compare com o valor da tabela de íons b/y da teoria (≈ 1149,6).
- **A.3** Calcule o m/z desse peptídeo com **duas cargas**: m/z = (M + 2 × 1,00728) / 2, onde M = 1148,61 Da é a massa neutra. *(Resposta esperada ≈ 575,3.)* No espectrômetro, a maioria dos peptídeos trípticos aparece com carga +2 ou +3.
- **A.4** Mude *missed cleavages* para **1** e rode de novo. O que acontece com o número de peptídeos? Por que o motor de busca precisa considerar cortes "perdidos"? *(Pista: na prática a tripsina não é perfeita; e lembre a exceção da prolina.)*

### Passo A3 — Conectando ao PRIDE (5 min)

1. Abra <https://www.ebi.ac.uk/pride/archive/projects/PXD006914>.
2. Identifique: o instrumento usado, a espécie, o artigo associado e a lista de arquivos.
3. **Pergunta A.5:** qual o tamanho aproximado de cada arquivo `.raw`? Multiplique por 19 amostras. Isso ajuda a entender por que análises de proteômica pedem computação robusta (e por que hoje partimos de resultados pré-processados).

---

## Parte B — Estatística em dados reais no Galaxy (50 min)

**Ferramentas:** Galaxy Europe e MSstats.
**Meta:** descobrir quais proteínas diferem de forma significativa entre cSCC metastático e cSCC-RDEB.

> **Dica de ouro:** o MSstats leva alguns minutos para rodar no servidor público (e pode haver fila). Por isso, **inicie o MSstats o quanto antes** (passo B4) e use o tempo de espera para explorar o relatório de qualidade (B5).

### Passo B1 — Crie um histórico (2 min)

1. Logado em <https://usegalaxy.eu>, clique no ícone **+** (*Create new history*) no painel da direita.
2. Dê o nome **Aula5_Proteomica_cSCC**.

*No Galaxy, um **histórico** é como um caderno de laboratório: guarda cada arquivo e cada passo executado, na ordem.*

### Passo B2 — Importe os dados do Zenodo (5 min)

1. Clique em **Upload** (ícone de seta para cima, à esquerda) → **Paste/Fetch data**.
2. Cole os **seis** endereços abaixo (um por linha) e clique em **Start**, depois em **Close**:

```
https://zenodo.org/records/17177151/files/input_protein_database.fasta
https://zenodo.org/records/17177151/files/input_annotation_file.tabular
https://zenodo.org/records/17177151/files/input_comparison_matrix.tabular
https://zenodo.org/records/17177151/files/MaxQuant%20Evidence.tabular
https://zenodo.org/records/17177151/files/MaxQuant%20Protein%20Groups.tabular
https://zenodo.org/records/17177151/files/PTXQC%20report.pdf
```

3. Aguarde os arquivos ficarem **verdes** (cinza = na fila; amarelo = processando; verde = pronto; vermelho = erro).
4. Renomeie cada arquivo (ícone de **lápis** → campo *Name* → *Save*):

| Arquivo importado | Renomeie para | O que é |
|---|---|---|
| input_protein_database.fasta | `protein database` | Sequências de proteínas humanas usadas na busca |
| input_annotation_file.tabular | `annotation file` | **Desenho experimental**: qual arquivo bruto é de qual paciente e de qual grupo |
| input_comparison_matrix.tabular | `comparison matrix` | Diz ao MSstats **quais grupos comparar** |
| MaxQuant Evidence.tabular | `Evidence` | Saída do MaxQuant: cada peptídeo medido em cada corrida |
| MaxQuant Protein Groups.tabular | `Protein Groups` | Saída do MaxQuant: cada grupo de proteínas |
| PTXQC report.pdf | `PTXQC report` | Relatório automático de **controle de qualidade** do MaxQuant |

5. Clique no **olho** 👁 sobre `annotation file` e `comparison matrix`. Veja como um desenho experimental é descrito: colunas *Raw.file*, *Condition*, *BioReplicate*, *Run*, *IsotopeLabelType*. Na matriz de comparação, **1** e **−1** indicam os dois grupos comparados.

### Passo B3 — Prepare as tabelas (8 min)

O MaxQuant marca automaticamente certas proteínas como **"potenciais contaminantes"** (com um `+` numa coluna própria) — por exemplo, queratinas humanas, que costumam vir da pele de quem manipula a amostra. **Mas aqui a amostra É pele!** Excluir essas proteínas jogaria fora biologia real. Por isso vamos: (1) remover apenas as proteínas **não humanas** e (2) "desmarcar" os contaminantes humanos.

Use a caixa de busca de ferramentas (à esquerda) para encontrar cada ferramenta pelo nome.

1. **Select** (*Select lines that match an expression*)
   - *Select lines from:* `Protein Groups`
   - *the pattern:* `(HUMAN)|(Majority)`
   - Renomeie a saída: `select protein groups`
2. **Select**
   - *Select lines from:* `Evidence`
   - *the pattern:* `(HUMAN)|(Sequence)`
   - Renomeie a saída: `select evidence`

   *Por que "Majority" e "Sequence"? São palavras que aparecem na **linha de cabeçalho** de cada tabela — assim o cabeçalho é mantido.*

3. **Replace Text in a specific column**
   - *File to process:* `select protein groups`
   - *in column:* `Column: 138` — *Find pattern:* `+` — *Replace with:* deixe **vazio**
   - Renomeie: `protein groups input for MSstats`
4. **Replace Text in a specific column**
   - *File to process:* `select evidence`
   - *in column:* `Column: 50` — *Find pattern:* `+` — *Replace with:* **vazio**
   - Renomeie: `evidence input for MSstats`

> **Confira antes de rodar:** clique no 👁 de `Protein Groups` e confirme que a coluna 138 é *Potential contaminant* (e a 50 em `Evidence`). Se o número for outro na sua versão, use o número correto.

### Passo B4 — Rode o MSstats (5 min para configurar; a execução segue em segundo plano)

Busque **MSstats** e configure:

- *input source:* `MaxQuant`
  - *evidence.txt – feature-level data:* `evidence input for MSstats`
  - *proteinGroups.txt – protein-level data:* `protein groups input for MSstats`
  - *annotation file:* `annotation file`
  - *Select Protein ID in evidence.txt:* `Leading razor protein column`
- Em **MaxQtoMSstatsFormat Options** → *Remove the proteins which have only 1 peptide and charge:* `Yes`
- Em **dataProcess Options** → *Select outputs:* marque `MSstats log`, `MSstats ProteinLevelData`, `Sample Quantification Matrix Table`
  - *dataProcess Plot Options* → *Select visualization outputs:* `MSstats QCPlot`; *Select protein IDs to draw plots:* `allonly`
- *Compare Groups:* `Yes`
  - *Comparison Matrix:* `comparison matrix`
  - *Select outputs:* `MSstats ComparisonResult`
  - *Comparison Visualization Options* → *Select visualization outputs:* `MSstats VolcanoPlot`
    - *Advanced visualization parameters* → *Involve fold change cutoff:* `1.5`; *logarithm transformation of adjusted p-value:* `2`; *Display protein names in Volcano Plot:* `No`

Clique em **Run Tool**. Enquanto roda, siga para o B5.

### Passo B5 — Enquanto espera: o controle de qualidade (10 min)

Abra o `PTXQC report` (👁). O **PTXQC** gera automaticamente gráficos de qualidade de uma análise do MaxQuant. Folheie e responda:

- **B.1** Quantas linhas tem o `Protein Groups` original (o Galaxy mostra o número de linhas ao expandir o arquivo)? Cada linha, exceto o cabeçalho, é um grupo de proteínas identificado com FDR de 1%. *(Esperado: ≈ 2.500 grupos de proteínas — 2.512 no tutorial de referência.)*
- **B.2** Procure no relatório a página sobre **contaminantes**. Que fração da intensidade das amostras vem de "potenciais contaminantes"? Por que isso não significa, aqui, que as amostras estão sujas?
- **B.3** Procure o gráfico de **intensidades por amostra**. Alguma amostra parece muito diferente das demais? O que isso poderia indicar tecnicamente?

### Passo B6 — Interprete os resultados do MSstats (15 min)

Quando as saídas do MSstats ficarem verdes:

1. **MSstats log** — procure quantas proteínas entraram na análise estatística. *(Esperado: ≈ 1.700 — 1.693 no tutorial de referência. Por que menos que no B.1? Lembre do filtro de "apenas 1 peptídeo" e da remoção de não humanas.)*
2. **MSstats QCPlot** — *boxplots* das intensidades em log₂ por corrida, **depois da normalização**. As medianas devem estar alinhadas. Esse é o efeito da normalização que vimos no Bloco 5.
3. **MSstats VolcanoPlot** — o gráfico principal da aula:
   - **eixo x:** log₂ fold change da comparação (metastático − RDEB). Valores **positivos** = mais abundante no cSCC **metastático**; **negativos** = mais abundante no cSCC-**RDEB**;
   - **eixo y:** −log₂ do valor-p **ajustado**;
   - linhas tracejadas: valor-p ajustado = 0,05 (horizontal) e fold change = 1,5 (verticais, em ±0,58 no log₂);
   - pontos coloridos = **significativos e com efeito relevante**.
4. **MSstats ComparisonResult** — a tabela por trás do gráfico. Colunas importantes: **Protein**, **log2FC**, **pvalue**, **adj.pvalue**.

Para filtrar a tabela, use a ferramenta **Filter** (*Filter data on any column using simple expressions*):

- *Filter:* `MSstats ComparisonResult`; *With following condition:* `c8<0.05`; *Number of header lines to skip:* `1`

*(Antes, confirme com o 👁 que a coluna 8 é `adj.pvalue`; ajuste se necessário.)*

**Perguntas:**

- **B.4** Quantas proteínas têm valor-p ajustado < 0,05? *(Esperado: em torno de 200; a grande maioria mais abundante no cSCC metastático.)*
- **B.5** Por que filtramos pelo valor-p **ajustado** e não pelo valor-p comum?
- **B.6** Procure na tabela o **colágeno XIV** (código UniProt **Q05707**, gene *COL14A1*). Em que direção ele mudou? Isso faz sentido sabendo que a RDEB é causada pela **deficiência de colágeno VII**? *(No artigo original, o colágeno XIV — mais abundante no cSCC-RDEB, possivelmente compensando a falta de colágeno VII — foi confirmado por coloração com anticorpos. O mesmo ocorreu com dois achados no sentido oposto, mais abundantes no metastático: o **componente P do amiloide sérico**, P02743, e a proteína de reparo de DNA **XRCC6**, P12956. A proteômica gerou as hipóteses; outra técnica as validou.)*
- **B.7** Escolha **uma** proteína significativa que você não conheça, abra a entrada dela no UniProt e escreva **uma frase** sobre sua função.

> **Atenção a um artefato comum:** às vezes surgem pontos "grudados" nos extremos do *volcano plot* (fold change ±infinito). Normalmente são proteínas **medidas em apenas um dos grupos** — não são necessariamente os achados mais importantes, e merecem cautela.

### Passo B7 — (Opcional, se sobrar tempo ou em casa) Da lista à rede no STRING (≈ 10 min)

1. No Galaxy, use **Filter** sobre `MSstats ComparisonResult` com a condição `c8<0.05 and c3>0.58` (*header lines to skip:* `1`). Resultado: as proteínas significativamente **mais abundantes no metastático**.
2. Os códigos na coluna 1 vêm no formato `sp|P02743|SAMP_HUMAN`. Para deixar só o código UniProt, use **Replace Text in a specific column** na coluna 1 com duas substituições (*Insert Replacement*): `sp\|` → vazio e `\|.*` → vazio.
3. Abra o resultado (👁) e copie a coluna de códigos.
4. Acesse <https://string-db.org> → **Search** → **Multiple proteins** → cole a lista → *Organism:* **Homo sapiens** → **Search** → **Continue**.
5. Explore a rede: cada **círculo é uma proteína** e cada **linha, uma interação**; as cores das linhas indicam o tipo de evidência. Em **Settings**, escolha *minimum required interaction score* = **high confidence (0.700)** e clique em *Update*.
6. Na aba **Analysis**, veja as funções **enriquecidas** (processos biológicos, vias KEGG) — com FDR, como no *volcano plot*.

- **B.8** Que processos biológicos aparecem? Fazem sentido para um tumor de pele (pense em matriz extracelular, adesão, resposta imune)?
- **B.9** Quais proteínas são **hubs** (muitas conexões)? Por que elas seriam boas candidatas para investigar primeiro?

> Quem preferir pode fazer este passo em Python no caderno **Aula05_Proteomica_Pratica.ipynb** (seção B11 — [abrir no Colab](https://colab.research.google.com/drive/1guMQyhhtXbJcfGYorPtBkc58t9UBWwJJ)), que envia a lista ao STRING pela API e desenha a rede.

### Plano B (se o Galaxy estiver lento ou fora do ar)

O professor compartilhará na tela um histórico já executado. A resolução de referência pública da Galaxy Training Network também pode ser consultada no tutorial *"MaxQuant and MSstats for the analysis of label-free data"* (<https://gxy.io/GTN:T00219>), na opção *Answer Histories*.

---

# MOMENTO 3 — Conexão teoria ↔ prática e fechamento (2:30 – 3:00)

## 3.1 — Cada clique, um conceito (10 min)

| O que você fez na prática | Conceito da teoria por trás |
|---|---|
| Buscou P02768 no UniProt e viu o peptídeo-sinal removido | Por que medir proteínas e não só genes (Bloco 1) |
| Digeriu a albumina no PeptideMass | Regra da tripsina e digestão *in silico* (Blocos 3 e 4) |
| Calculou m/z com carga +2 | Ionização e a razão m/z (Bloco 2) |
| Testou *missed cleavages* | Tamanho do espaço de busca do motor de busca (Bloco 4) |
| Abriu o depósito PXD006914 no PRIDE | ProteomeXchange, arquivos brutos e reprodutibilidade (Blocos 3 e 6) |
| Partiu do `Protein Groups` e do `Evidence` do MaxQuant | Identificação com FDR 1% e inferência de grupos de proteínas (Bloco 4) |
| Decidiu **manter** contaminantes humanos | Conhecer a biologia da amostra muda a análise — o software não pensa por você |
| Informou o desenho experimental (annotation file) | Réplicas biológicas e modelos estatísticos (Bloco 5) |
| Observou medianas alinhadas no QCPlot | Normalização (Bloco 5) |
| Filtrou por valor-p **ajustado** | Testes múltiplos e FDR na quantificação (Bloco 5) |
| Interpretou o *volcano plot* e o colágeno XIV | Da lista estatística à hipótese biológica — que depois exige validação |
| (Opcional) Levou a lista ao STRING | Interatoma: proteínas agem em complexos e vias (Bloco 6) |

## 3.2 — E no mundo real? Como isso se conecta à defesa (7 min)

- **Do *volcano plot* ao ensaio de campo:** a proteômica de descoberta (o que fizemos) gera uma lista de candidatos. Depois, um ensaio **direcionado** (SRM/PRM, montado no **Skyline**) mede só esses alvos, com sensibilidade e exatidão muito maiores. É assim que se desenvolve, por exemplo, um teste confirmatório para uma toxina proteica ou para um marcador de exposição a um agente químico.
- **Mesmo fluxo, outra pergunta:** trocando o FASTA humano pelo proteoma de uma bactéria, o mesmo raciocínio identifica um microrganismo a partir de suas proteínas.
- **Do interatoma às contramedidas:** saber com quem uma toxina ou proteína viral interage no hospedeiro aponta onde intervir. Hoje demos o primeiro passo — transformar uma lista em rede no STRING.
- **Integração com outras aulas:** a lista de proteínas alteradas alimenta a **análise de vias** (KEGG, Reactome) e de **redes** (STRING) da Unidade 9, agora com mais profundidade; uma proteína de interesse pode ter sua estrutura analisada na Unidade 7 (PDB, AlphaFold) e sua história evolutiva na Unidade 8.

## 3.3 — Ideias para o trabalho final (5 min)

O artigo final deve integrar **pelo menos duas técnicas** da disciplina. Exemplos que partem desta aula:

- Reanalisar um conjunto **PXD** público ligado à sua linha de pesquisa (com Galaxy/MSstats ou Perseus) **+** análise de vias/redes das proteínas alteradas (KEGG/STRING);
- Escolher uma proteína-alvo, prever seus peptídeos trípticos (PeptideMass) e propor um ensaio direcionado **+** modelagem estrutural da proteína (AlphaFold/SWISS-MODEL) para discutir a acessibilidade desses peptídeos;
- Identificar proteínas diferenciais **+** buscar homólogos (BLAST) e construir uma filogenia de uma família de interesse.

## 3.4 — Fechamento (8 min)

**Cinco ideias para levar:**

1. A espectrometria de massas mede **m/z**; MS1 serve sobretudo para **quantificar**, MS2 para **identificar**.
2. No **bottom-up**, proteínas viram peptídeos (tripsina), são separados (LC), ionizados (ESI) e fragmentados (íons **b** e **y**).
3. A identificação compara espectros com um banco de sequências, e a estratégia **target-decoy** controla a **FDR em 1%**.
4. A quantificação exige **log₂, normalização, tratamento de valores ausentes e correção de testes múltiplos** — e um bom **desenho experimental** (réplicas biológicas!).
5. Ferramentas resolvem o cálculo; a **interpretação depende de conhecer a amostra e a pergunta** (lembre dos "contaminantes" que eram pele).

**Perguntas de autoavaliação (responda em 2–3 linhas cada):**

1. Por que a tripsina é a enzima preferida na proteômica *bottom-up*?
2. Explique a estratégia *target-decoy* para um colega de engenharia que nunca ouviu falar de proteômica.
3. Qual a diferença entre DDA e DIA, e qual delas tende a gerar mais valores ausentes?
4. Em um *volcano plot*, onde estão as proteínas que mudaram muito, mas sem significância estatística? E as significativas com efeito pequeno?
5. Você quer verificar se uma toxina específica está presente numa amostra ambiental. Qual estratégia de quantificação e qual ferramenta usaria? Por quê?
6. Qual a diferença entre uma **interação física** e uma **associação funcional** no STRING? E por que uma rede de coexpressão não é um interatoma?

---

## Atividade assíncrona (para a lista de exercícios)

1. **Complete a Parte B**, caso não tenha terminado, e salve o *volcano plot*.
2. **Separe as proteínas** significativas em dois grupos (mais abundantes no metastático: `c3>0.58`; mais abundantes no RDEB: `c3<-0.58`, usando a ferramenta **Filter** sobre o resultado filtrado por `c8<0.05`). Quantas há em cada grupo?
3. **Escolha 5 proteínas** de cada grupo e pesquise no UniProt sua função. Elas sugerem algum processo biológico em comum (por exemplo, matriz extracelular, resposta inflamatória)?
4. **Faça o Passo B7** (STRING) para os dois grupos e compare as redes e os processos enriquecidos. Teste também a opção de rede **física** (*Settings → network type → physical*): a rede encolhe? Por quê?
5. **Caderno em Python (opcional):** o arquivo **Aula05_Proteomica_Pratica.ipynb** refaz toda a prática em Python comentado no Google Colab, inclusive a rede do STRING. [Abrir no Colab](https://colab.research.google.com/drive/1guMQyhhtXbJcfGYorPtBkc58t9UBWwJJ) · [Arquivo no Google Drive](https://drive.google.com/file/d/1guMQyhhtXbJcfGYorPtBkc58t9UBWwJJ/view?usp=sharing). Use **Arquivo → Salvar uma cópia no Drive** antes de editar.
6. **Desafio opcional (para quem quiser um primeiro contato com Python):** abra o Google Colab (<https://colab.research.google.com>), crie um caderno novo e cole o código abaixo, célula a célula. Ele faz automaticamente o que você fez no PeptideMass. Não é preciso entender cada linha — leia os comentários (linhas que começam com `#`).

```python
# Célula 1 — instala a biblioteca pyteomics (ferramentas de proteômica em Python, código aberto)
!pip install -q pyteomics

# Célula 2 — baixa a sequência da albumina humana direto do UniProt
import requests
fasta = requests.get("https://rest.uniprot.org/uniprotkb/P02768.fasta").text
sequencia = "".join(fasta.splitlines()[1:])   # junta as linhas, ignorando o cabeçalho ">"
print("Tamanho da proteína:", len(sequencia), "aminoácidos")

# Célula 3 — digestão in silico com tripsina (sem cortes perdidos)
from pyteomics import parser, mass
peptideos = parser.cleave(sequencia, parser.expasy_rules["trypsin"], missed_cleavages=0)

# Célula 4 — calcula a massa e o m/z (carga +2) de cada peptídeo e mostra os maiores que 500 Da
tabela = []
for p in peptideos:
    m = mass.calculate_mass(sequence=p)              # massa monoisotópica neutra
    if m > 500:
        mz2 = mass.calculate_mass(sequence=p, charge=2)  # m/z com duas cargas
        tabela.append((p, round(m, 2), round(mz2, 2)))
tabela.sort(key=lambda linha: linha[1])
print(len(tabela), "peptídeos acima de 500 Da")
for p, m, mz2 in tabela[:15]:
    print(f"{p:<30} massa={m:>9}  m/z(2+)={mz2:>8}")
```

*Pergunta do desafio:* o número de peptídeos bate com o que você obteve no PeptideMass? Se não bater exatamente, pense em por quê (dica: a sequência do UniProt inclui o peptídeo-sinal e o propeptídeo).

---

## Glossário rápido

| Termo | Significado |
|---|---|
| **Da (Dalton)** | Unidade de massa molecular (≈ massa de um átomo de H) |
| **m/z** | Razão massa/carga — o que o espectrômetro de fato mede |
| **ESI / MALDI** | Técnicas de ionização (borrifo sob alta tensão / laser sobre matriz) |
| **MS1 / MS2** | Espectro dos peptídeos inteiros / espectro dos fragmentos |
| **LC-MS/MS** | Cromatografia líquida acoplada à espectrometria de massas em tandem |
| **Bottom-up / shotgun** | Estratégia que digere proteínas em peptídeos antes da medida |
| **DDA / DIA** | Aquisição dependente / independente de dados |
| **Íons b e y** | Fragmentos com o início (b) ou o fim (y) do peptídeo |
| **PSM** | Par espectro–peptídeo atribuído pelo motor de busca |
| **FASTA** | Formato de texto para sequências biológicas |
| **Target-decoy** | Uso de um banco falso para estimar falsos positivos |
| **FDR** | Taxa de falsas descobertas (proporção esperada de erros entre os aceitos) |
| **Protein group** | Conjunto de proteínas indistinguíveis pelos peptídeos observados |
| **LFQ** | Quantificação sem marcação (*label-free*) |
| **SILAC / TMT** | Quantificação com marcação metabólica / isobárica |
| **SRM / PRM** | Proteômica direcionada (monitora alvos pré-definidos) |
| **log₂ fold change** | Variação entre grupos em escala log₂ (+1 = dobro) |
| **Valor-p ajustado** | Valor-p corrigido para muitos testes simultâneos |
| **PTM** | Modificação pós-traducional |
| **PXD** | Identificador de um conjunto de dados no ProteomeXchange |
| **Interatoma** | Conjunto de todas as interações entre as proteínas de uma célula |
| **AP-MS** | Purificação por afinidade seguida de espectrometria de massas: identifica quem vem junto com uma proteína-isca |
| **BioID / APEX** | Marcação por proximidade: carimba as proteínas vizinhas da isca na célula viva |
| **XL-MS** | Ligação cruzada + MS: "grampeia" proteínas que se tocam |
| **Hub** | Proteína com muitas conexões numa rede |
| **STRING** | Banco de interações e associações funcionais entre proteínas |

---

## Referências e links

- Aebersold, R.; Mann, M. Mass-spectrometry-based proteomics. *Nature* 422, 198–207 (2003).
- Cox, J.; Mann, M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. *Nature Biotechnology* 26, 1367–1372 (2008).
- Tyanova, S.; Temu, T.; Cox, J. The MaxQuant computational platform for mass spectrometry-based shotgun proteomics. *Nature Protocols* 11, 2301–2319 (2016).
- Santos, M. D. M. *et al.* Simple, efficient and thorough shotgun proteomic analysis with PatternLab V. *Nature Protocols* (2022).
- Choi, M. *et al.* MSstats: an R package for statistical analysis of quantitative mass spectrometry-based proteomic experiments. *Bioinformatics* 30, 2524–2526 (2014).
- Käll, L.; Vitek, O. Computational mass spectrometry–based proteomics. *PLoS Computational Biology* 7, e1002277 (2011).
- Föll, M. C. *et al.* Identification of tissue damage, extracellular matrix remodeling and bacterial challenge as common mechanisms associated with high-risk cutaneous squamous cell carcinomas. *Matrix Biology* 66, 1–21 (2018).
- Szklarczyk, D. *et al.* The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research* 51, D638–D646 (2023).
- Huttlin, E. L. *et al.* Dual proteome-scale networks reveal cell-specific remodeling of the human interactome. *Cell* 184, 3022–3040 (2021).
- Föll, M.; Fahrner, M. *MaxQuant and MSstats for the analysis of label-free data* (Galaxy Training Materials) — <https://gxy.io/GTN:T00219>
- MaxQuant e Perseus: <https://www.maxquant.org>
- PatternLab for Proteomics: <http://patternlabforproteomics.org>
- Skyline: <https://skyline.ms>
- UniProt: <https://www.uniprot.org> · PRIDE: <https://www.ebi.ac.uk/pride> · ProteomeXchange: <https://www.proteomexchange.org>
- Expasy PeptideMass: <https://web.expasy.org/peptide_mass/>
- Caderno da prática em Python (Google Colab): <https://colab.research.google.com/drive/1guMQyhhtXbJcfGYorPtBkc58t9UBWwJJ>
- STRING: <https://string-db.org> · BioGRID: <https://thebiogrid.org> · IntAct: <https://www.ebi.ac.uk/intact>
