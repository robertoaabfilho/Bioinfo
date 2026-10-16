# Aula 8 — Filogenética e Evolução Molecular

**Reconstruindo a história evolutiva a partir de sequências · Alinhamento múltiplo · Árvores de máxima verossimilhança · Relógio molecular**

Introdução à Bioinformática — PPGED/IME
Duração: 3 horas (180 min) · Modalidade: online síncrona
Docentes: Prof. Mariana Quezado (IME) e [seu nome]

---

## Antes de começar: o que você precisa (5 min de preparo, faça na véspera)

Nada será instalado na sua máquina. Tudo roda no **navegador**, no Google Colab.

| Item | Para quê | Link |
|---|---|---|
| Conta Google | Rodar o Google Colab (Python na nuvem) | [colab.research.google.com](https://colab.research.google.com) |
| Navegador atualizado (Chrome/Edge/Firefox) | Colab e visualizadores online de árvore | — |
| O notebook `Aula08_Filogenetica_Pratica.ipynb` | O roteiro prático executável | fornecido pela disciplina |

> **Você não precisa saber programar.** Nas células de código do Colab você vai apenas clicar no ▶ e ler o que sai. O objetivo é **executar, observar e interpretar** — exatamente como nas aulas anteriores.

> **Se a rede da sua instituição bloquear downloads automáticos:** algumas redes corporativas e militares bloqueiam o acesso direto a repositórios. Se a primeira célula de download falhar, baixe os arquivos `sequences.fasta` e `metadata.tsv` manualmente pelo navegador (as URLs estão na célula) e arraste-os para o painel de arquivos do Colab (ícone de pasta na barra lateral esquerda). O resto do notebook funciona igual. **Teste isso na véspera, não durante a aula.**

---

## Mapa da aula

| Momento | Tempo | O que acontece |
|---|---|---|
| **1 — Teoria** | 00:00 – 01:15 (75 min) | Por que reconstruir filogenias; homologia e alinhamento; métodos de construção de árvore; modelos de substituição; suporte de ramo; relógio molecular |
| **Intervalo** | 01:15 – 01:25 (10 min) | Café |
| **2 — Prática** | 01:25 – 02:40 (75 min) | Hands-on com dados públicos e software livre: FASTA → alinhamento → árvore → datação |
| **3 — Integração e fechamento** | 02:40 – 03:00 (20 min) | Costura teoria↔prática, erros comuns, como usar isso no artigo final, exercícios |

**O fio condutor da aula** é uma pergunta concreta de **epidemiologia molecular**:

> **Quando e por onde o vírus Zika chegou às Américas?**

Entre 2015 e 2016, uma grande epidemia de Zika atingiu as Américas. A pergunta "*desde quando o vírus circulava e qual foi a rota de entrada?*" não foi respondida por notificações de casos — foi respondida **lendo e comparando os genomas do vírus**. Vamos refazer essa análise, do arquivo bruto de sequências até a data estimada da introdução, usando sempre dados públicos (GenBank) e ferramentas de código aberto.

---

# MOMENTO 1 — Aula teórica (75 min)

## 1.1 O que é uma filogenia e por que ela importa (10 min)

Nas aulas anteriores trabalhamos com **sequências** (DNA, RNA, proteínas) e com **estruturas**. A filogenética faz uma pergunta diferente sobre as mesmas sequências:

> **Qual é a história evolutiva que relaciona estas sequências entre si?**

O resultado dessa pergunta é uma **árvore filogenética**: um diagrama em que as pontas são as sequências que observamos hoje e os pontos de ramificação (nós internos) são **ancestrais comuns hipotéticos**, que existiram no passado mas não foram observados.

Analogia de engenharia: é como reconstruir a **árvore de versões** de um projeto a partir apenas dos arquivos finais, sem ter o histórico de commits. As semelhanças e diferenças entre as versões atuais permitem inferir quais derivaram de quais e em que ordem.

### O que a filogenética responde

1. **Quem é parente de quem?** (topologia da árvore)
2. **Quanto divergiram?** (comprimento dos ramos)
3. **Qual a direção do tempo?** (onde fica a raiz, o ancestral de todos)
4. **Quando** os ancestrais existiram? (datação por relógio molecular)
5. **Por onde** um organismo/patógeno se dispersou? (filogeografia)

### Onde isso aparece na prática

| Campo | Pergunta filogenética |
|---|---|
| **Epidemiologia molecular** | De onde veio um surto? Desde quando circula? Como se espalhou? |
| **Vigilância genômica** | Uma nova amostra pertence a qual linhagem? É uma variante nova? |
| **Evolução de proteínas** | Quais resíduos são conservados (essenciais) ao longo da evolução? |
| **Taxonomia e biodiversidade** | Duas amostras são a mesma espécie? Qual a relação entre elas? |
| **Identificação de origem** | Uma amostra desconhecida é mais próxima de qual referência conhecida? |

> **Escopo desta aula.** Trabalhamos a filogenética como ferramenta de **leitura do passado a partir de dados públicos** — rastrear a origem e a dispersão de um patógeno, identificar a linhagem de uma amostra, medir conservação evolutiva. É também esse o enquadramento esperado no seu artigo final.

### Nosso caso de hoje: a chegada do Zika às Américas

- O vírus Zika é um **flavivírus** (parente da dengue e da febre amarela), com genoma de **RNA de fita simples de ~10,8 mil bases**.
- Vírus de RNA acumulam mutações **rápido** — e é justamente essa velocidade que transforma o genoma num "relógio" e num "diário de viagem".
- A hipótese histórica, que vamos **testar**, é: o vírus saiu do Sudeste Asiático → passou pela Oceania (epidemia na Polinésia Francesa, 2013–2014) → chegou às Américas, onde circulou **silenciosamente** antes de ser detectado.

---

## 1.2 O pré-requisito de tudo: homologia e alinhamento (12 min)

### Homologia

Duas posições (colunas) ou dois caracteres são **homólogos** quando descendem de um **mesmo ancestral comum**. Toda a filogenética repousa sobre uma suposição: a de que estamos comparando **coisas homólogas** — a mesma posição do genoma em diferentes organismos, não posições quaisquer que por acaso têm a mesma letra.

### Por que alinhar primeiro

Os genomas que baixamos **não têm o mesmo comprimento** — faltam pedaços nas pontas, há inserções e deleções (*indels*) no meio. Logo, "a posição 2.500 do genoma A" **não** é necessariamente "a posição 2.500 do genoma B". Comparar as sequências diretamente seria comparar coisas não homólogas.

O **alinhamento múltiplo de sequências (MSA)** resolve isso inserindo lacunas (`-`, *gaps*) de modo que cada **coluna** passe a conter caracteres homólogos. Depois de alinhar:

- todas as sequências têm o **mesmo comprimento** (número de colunas);
- cada coluna representa **uma posição evolutiva** comparável em todos.

> **Mensagem-chave:** o alinhamento é a etapa mais crítica e mais subestimada da filogenética. **Um alinhamento ruim produz uma árvore ruim com aparência perfeita.** O lixo entra silenciosamente e sai como um número plausível.

### A ferramenta: MAFFT

| Ferramenta | O que é | Quando usar |
|---|---|---|
| **MAFFT** | Alinhador múltiplo rápido e preciso (*Multiple Alignment using Fast Fourier Transform*). Código aberto. **É o que usaremos** | Padrão para a maioria dos casos; escala bem |
| **MUSCLE** | Outro alinhador clássico, boa precisão | Alternativa ao MAFFT (visto na Aula 2) |
| **Clustal Omega** | Muito usado, bom para muitas sequências de proteína | Grandes conjuntos de proteínas |
| **MAFFT --auto** | Deixa o programa escolher a estratégia conforme o tamanho | O que usamos hoje: simples e confiável |

Para sequências muito divergentes, visualiza-se o alinhamento (ex.: **Jalview**, **pyMSAviz**) e às vezes se removem colunas mal alinhadas (*trimming*, com ferramentas como trimAl ou Gblocks) antes de construir a árvore.

---

## 1.3 Como se constrói uma árvore (18 min)

Existem quatro grandes famílias de métodos. Conhecê-las é importante para **escolher** e para **ler a metodologia de artigos**.

### A) Métodos de distância (Neighbor-Joining, UPGMA)

Calcula-se uma **matriz de distâncias** entre todos os pares (quantas posições diferem, corrigidas pelo modelo) e, a partir dela, monta-se a árvore.

- **Vantagem:** muito rápido; bom para uma primeira olhada ou para milhares de sequências.
- **Limitação:** joga toda a informação das sequências numa única distância por par; perde sinal. O **Neighbor-Joining (NJ)** é o mais usado dessa família.

### B) Máxima parcimônia (MP)

Procura a árvore que exige o **menor número total de mudanças** (mutações) para explicar os dados.

- **Vantagem:** intuitivo ("a explicação mais simples").
- **Limitação:** pode se enganar sistematicamente quando as taxas de evolução são muito desiguais entre ramos (fenômeno da *long-branch attraction*).

### C) Máxima verossimilhança (ML) — **o que usaremos**

Dado um **modelo** explícito de como os nucleotídeos mudam ao longo do tempo, a **verossimilhança** mede *o quão provável é observar as sequências que de fato temos*, para uma determinada árvore (topologia + comprimentos de ramo). O método procura a árvore que **maximiza** essa probabilidade.

- **Vantagem:** estatisticamente bem fundamentado; usa toda a informação das colunas; é o padrão atual junto com a inferência bayesiana.
- **Limitação:** exige escolher um modelo; é mais pesado computacionalmente (mas rápido para conjuntos pequenos como o nosso).
- **Ferramentas:** **IQ-TREE** (o que usaremos), **RAxML**, **PhyML**.

### D) Inferência bayesiana (BI)

Em vez de uma única "melhor" árvore, estima a **distribuição de probabilidade** das árvores dadas as sequências, por amostragem (MCMC). Entrega, para cada agrupamento, uma **probabilidade posterior**.

- **Vantagem:** incorpora incerteza de forma natural; integra modelos complexos (ex.: datação com relógio relaxado).
- **Limitação:** mais lento; exige cuidado com convergência.
- **Ferramentas:** **MrBayes**, **BEAST/BEAST2** (este último é o padrão para datação filogenética sofisticada).

> **Panorama das ferramentas de hoje**
>
> | Ferramenta | Papel | Licença |
> |---|---|---|
> | **MEGA** | Suíte gráfica (alinhamento + árvore NJ/MP/ML + testes), ótima para aprender clicando | Gratuita |
> | **IQ-TREE 2** | Inferência por ML com seleção automática de modelo e bootstrap ultrarrápido. Linha de comando | Código aberto (GPL) |
> | **Jalview** | Visualização e edição de alinhamentos | Gratuita |
> | **FigTree / iTOL** | Visualização e anotação de árvores (figuras de artigo) | Gratuita / web |
> | **TreeTime** | Relógio molecular e datação de árvores | Código aberto |
>
> **MEGA** é a porta de entrada clicável; **IQ-TREE** é o motor reprodutível que roda no notebook de hoje. Os dois chegam a árvores de ML — a diferença é que o fluxo em linha de comando é **documentável e repetível**, o que o seu artigo vai exigir.

---

## 1.4 Modelos de substituição: por que a árvore precisa de um (12 min)

Nem toda mudança de base é igualmente provável. Trocar uma purina por outra purina (A↔G) ou uma pirimidina por outra (C↔T) — as **transições** — é bioquimicamente muito mais frequente do que trocar purina por pirimidina (as **transversões**). Se ignorarmos isso, subestimamos o número real de mudanças em ramos longos (porque mudanças múltiplas na mesma posição "escondem" umas às outras).

Um **modelo de substituição** é um conjunto de regras (uma matriz de taxas) que descreve essas probabilidades. Do mais simples ao mais geral:

| Modelo | O que assume |
|---|---|
| **JC69** (Jukes-Cantor) | Todas as bases igualmente frequentes; todas as trocas com a mesma taxa (o mais ingênuo) |
| **K2P / K80** (Kimura) | Distingue transições de transversões |
| **HKY / TN93** | Distingue transições/transversões **e** permite frequências de base desiguais |
| **GTR** (*General Time Reversible*) | O mais geral: seis taxas de troca livres + frequências livres |

A esses adicionam-se **"temperos"** que descrevem a variação de taxa entre posições:

- **+I** (*invariant sites*): uma fração das colunas **nunca** muda (sítios estruturalmente essenciais). Vimos isso na prática — a maioria das colunas do Zika é invariável.
- **+G** (*gamma*): as demais posições evoluem a velocidades diferentes, descritas por uma distribuição gama.

### Como escolher? ModelFinder

Não se escolhe o modelo "no olho". O **ModelFinder** (embutido no IQ-TREE, acionado com `-m MFP`) testa dezenas de combinações e escolhe a melhor por um critério que **equilibra ajuste e complexidade**:

- **BIC** (*Bayesian Information Criterion*) e **AIC** penalizam modelos com parâmetros demais, evitando *overfitting*.

> Na prática de hoje, o ModelFinder escolhe um modelo da família **TN+F+I** para o Zika: Tamura–Nei (transições ≠ transversões), frequências de base empíricas (`+F`) e uma fração de sítios invariáveis (`+I`). Isso não é burocracia: um modelo pobre distorce os comprimentos de ramo e, portanto, **as datas** que estimaremos no fim.

---

## 1.5 Quão confiável é a árvore? Suporte de ramo (12 min)

Um programa **sempre** devolve uma árvore — inclusive de dados sem nenhum sinal. A pergunta honesta é: **quão sustentado por os dados está cada agrupamento?**

### Bootstrap

A ideia: se o sinal filogenético for real, ele deve aparecer mesmo que embaralhemos levemente os dados. O **bootstrap** reamostra as **colunas** do alinhamento (com reposição), reconstrói a árvore, e repete isso centenas ou milhares de vezes. Para cada agrupamento (ramo), conta-se **em que porcentagem das réplicas ele reaparece**.

- O **bootstrap ultrarrápido (UFBoot)** do IQ-TREE (`-B 1000`) faz 1000 réplicas em segundos.
- **Como ler** (regra prática para UFBoot):
  - **≥ 95%** → ramo bem sustentado, confiável.
  - 80–95% → sustentação moderada.
  - **< 80%** → fraco; trate o agrupamento com cautela.

> **Interpretação correta:** bootstrap alto significa que **os dados são internamente consistentes** com aquele agrupamento — **não** que ele seja uma verdade biológica absoluta. Ramos internos muito curtos (pouca evidência, como numa radiação explosiva) frequentemente têm suporte baixo, e isso é honesto: deve ser **relatado**, não escondido.

### Probabilidade posterior (inferência bayesiana)

Na abordagem bayesiana, o suporte é a **probabilidade posterior** do agrupamento (0 a 1). Valores ≥ 0,95 são considerados fortes. **Atenção:** probabilidade posterior e bootstrap **não são a mesma escala** e não devem ser comparados diretamente — a posterior costuma ser sistematicamente mais alta.

---

## 1.6 Enraizamento: onde fica o passado (6 min)

A árvore que o ML produz é **não enraizada**: mostra parentescos, mas não sabe, por si só, qual é a direção do tempo. Para ler a árvore como uma história, é preciso **enraizá-la** — definir o ancestral de todos.

Há duas estratégias:

1. **Grupo externo (*outgroup*)** — inclui-se na análise um organismo **sabidamente mais distante** de todos os demais (ex.: uma linhagem antiga, ou uma espécie irmã). A raiz fica no ramo que separa o grupo externo do resto. **É o método preferível** quando se tem um bom grupo externo. Na prática de hoje, usamos as cepas asiáticas (as mais basais e antigas) como grupo externo.
2. **Enraizamento no ponto médio (*midpoint*)** — na ausência de grupo externo, coloca-se a raiz no ponto médio do maior caminho entre duas pontas. É uma aproximação; supõe taxas de evolução parecidas entre as linhagens.

> **Por que isso importa tanto na datação:** a **escolha da raiz muda todas as datas dos nós internos**. Uma raiz mal colocada produz datas sem sentido (às vezes até "taxas negativas" no relógio molecular). É por isso que, no fim da prática, você vai testar o efeito de trocar a raiz — e é por isso que artigos sérios justificam explicitamente o enraizamento.

---

## 1.7 O relógio molecular: de substituições para anos (5 min)

Os comprimentos de ramo de uma árvore de ML estão em **substituições por sítio** — não em anos. Para responder "*desde quando?*", precisamos de uma ponte: o **relógio molecular**.

> **A hipótese do relógio:** se as mutações se acumulam a uma taxa aproximadamente constante no tempo, então **quanto mais recente** foi coletada uma amostra, **mais substituições** ela terá acumulado desde o ancestral. A inclinação dessa relação é a **taxa de evolução** (substituições por sítio por ano).

### Como se testa e se usa

- **Regressão raiz-às-pontas (*root-to-tip*):** para cada amostra, gráfico de "distância à raiz" × "data de coleta". Se os pontos sobem numa reta, há **sinal de relógio**; o **R²** mede quão clock-like é o conjunto (quanto mais perto de 1, melhor). A inclinação é a taxa.
- **Relógio estrito × relaxado:** o relógio estrito supõe a mesma taxa em todos os ramos; os **relógios relaxados** permitem que a taxa varie entre linhagens (mais realista; padrão no BEAST).
- **Datação de nós (tMRCA):** com o relógio calibrado pelas datas de coleta, estima-se a idade de cada nó interno — inclusive o **ancestral comum mais recente** (*time to the Most Recent Common Ancestor*) de qualquer clado de interesse.

> Vírus de RNA, como o Zika, têm taxas altas — tipicamente **~10⁻³ substituições por sítio por ano**. Isso significa que, mesmo em poucos anos, o genoma acumula mudanças suficientes para datar eventos com resolução de meses. É o que torna possível a pergunta de hoje.

**Ferramentas:** **TreeTime** (rápido, o que usamos — estima taxa, testa sinal de relógio e data a árvore) e **BEAST/BEAST2** (bayesiano, mais sofisticado, para artigos dedicados à datação).

---

# INTERVALO — 10 minutos

Abra o Colab e deixe o notebook `Aula08_Filogenetica_Pratica.ipynb` carregado (*Arquivo → Fazer upload de notebook*). Rode já a **primeira célula de instalação** — ela leva 1–2 minutos e assim estará pronta quando voltarmos.

---

# MOMENTO 2 — Aula prática (75 min)

**Todos os dados são públicos. Todo o software é livre.**

**Conjunto de dados:** 34 genomas completos do vírus Zika depositados no **GenBank**, com data e local de coleta (o mesmo conjunto público do tutorial do projeto **Nextstrain/zika**). Amostras de 2013 a 2016, vindas do Sudeste Asiático, Oceania e Américas.

**Ferramentas:** MAFFT (alinhamento), IQ-TREE 2 (árvore de ML + bootstrap), Biopython (leitura e manipulação), TreeTime (relógio molecular), pyMSAviz (visualização do alinhamento).

| Bloco | Tempo | O que faz |
|---|---|---|
| A — Os dados | 10 min | Carregar e inspecionar os 34 genomas e os metadados |
| B — Alinhamento | 15 min | Alinhar com MAFFT e visualizar o MSA |
| C — Árvore de ML | 35 min | ModelFinder + IQ-TREE + bootstrap; ler o modelo e o suporte |
| D — Enraizar e ler | 20 min* | Enraizar por grupo externo; testar monofilia do clado americano |
| E — Datar | 20 min* | Relógio molecular com TreeTime; estimar o tMRCA americano |

\* Os blocos D e E são mais curtos na execução; o tempo extra é para discussão guiada.

> **Todas as células de código do notebook foram executadas e verificadas** antes da aula (MAFFT 7.505, IQ-TREE 2.3.6, TreeTime 0.12.1, Biopython, pyMSAviz 0.5.0). Os números esperados estão indicados em cada bloco. Se o seu resultado divergir muito, algo deu errado — não siga adiante.

---

## Bloco A — Os dados: 34 genomas (10 min)

Baixamos dois arquivos: `sequences.fasta` (os genomas) e `metadata.tsv` (cepa, acesso GenBank, data, país, região).

**O que você vai ver e anotar no caderno de laboratório:**

- **34 genomas**, cada um com ~9.000–10.800 nucleotídeos. **Eles não têm o mesmo comprimento** — faltam pedaços nas pontas. Esse é o motivo pelo qual precisamos alinhar antes de comparar.
- Distribuição por região: Américas (do Sul e do Norte/Central), Oceania e Sudeste Asiático. Datas de **2013 a 2016**.

> **Pergunta de checagem:** por que precisamos das **datas** e dos **locais**, e não só das sequências? (Resposta: as datas calibram o relógio molecular; os locais permitem ler a árvore como uma história de dispersão.)

---

## Bloco B — Alinhamento múltiplo com MAFFT (15 min)

Rodamos `mafft --auto`. Em 34 genomas de vírus, leva **poucos segundos**.

**Valor esperado:** o alinhamento sai com **34 sequências × ~10.800 colunas** — agora todas do mesmo comprimento, porque o MAFFT preencheu as pontas ausentes com `-`.

Em seguida, visualizamos um trecho (colunas 3000–3120) com o pyMSAviz.

**O que observar:**
- Colunas em que todos têm a mesma letra são **posições conservadas**.
- As poucas colunas com letras diferentes são os **sítios variáveis** — são **eles** que carregam o sinal filogenético.
- Num vírus recente como o Zika, quase tudo é conservado: a árvore será construída a partir de um número relativamente pequeno de posições que diferem (no nosso conjunto, cerca de **190 sítios informativos para parcimônia**).

---

## Bloco C — A árvore de máxima verossimilhança com IQ-TREE (35 min)

Rodamos o IQ-TREE com três ingredientes:

- `-m MFP` → **ModelFinder** escolhe o modelo de substituição;
- `-B 1000` → **bootstrap ultrarrápido** com 1000 réplicas (suporte de ramo);
- `--seed 42` → **reprodutibilidade** (resultado repetível por outra pessoa).

Leva ~10–30 s. Em vez de abrir o relatório inteiro, extraímos as linhas-chave do arquivo `.iqtree`.

**Valores esperados:**
- **Modelo escolhido:** `TN+F+I` (Tamura–Nei + frequências empíricas + sítios invariáveis). O `+I` confirma, estatisticamente, o que vimos no Bloco B: a maioria das colunas não varia.
- Log-verossimilhança em torno de **−18.092** e comprimento total da árvore **~0,039 substituições/sítio** (ramos curtos — é um surto recente, pouca divergência acumulada).

> **Por que o modelo não é um detalhe:** usar um modelo pobre distorce os comprimentos de ramo e, portanto, as datas estimadas no Bloco E. O ModelFinder resolve isso objetivamente, pelo BIC.

---

## Bloco D — Lendo e enraizando a árvore (20 min)

A árvore de ML é **não enraizada**. Nós a lemos com o Biopython, desenhamos em **ASCII** (sem precisar de biblioteca gráfica) e a **enraizamos pelo grupo externo** — as cepas do **Sudeste Asiático** (2013), as mais basais.

**Como ler o desenho ASCII:**
- Cada ponta é um genoma (nome traz país/ano); os nós internos são ancestrais hipotéticos.
- O comprimento horizontal é proporcional ao número de substituições.
- Cepas no mesmo ramo compartilham um ancestral comum mais recente.

**O que observar:** as cepas do Sudeste Asiático e da Oceania ficam na base; **todas as 24 cepas das Américas formam um único grande grupo (um clado)**.

Depois **verificamos isso quantitativamente**: o ancestral comum das 24 cepas americanas contém **exatamente** essas 24 cepas (nenhuma intrusa) e tem **bootstrap de ~92%**.

**Interpretação:** um clado americano **monofilético e bem sustentado** é evidência forte de uma **introdução única** do vírus nas Américas, seguida de dispersão local — e não de várias entradas independentes. No artigo, escreveríamos: *"as cepas americanas formaram um clado monofilético bem sustentado (UFBoot ≈ 92%), compatível com uma introdução única do vírus no continente."*

> **Cuidado:** 92% é alto, mas **não** é certeza. Com poucos sítios variáveis, alguns ramos internos têm suporte menor. Relatar o valor — e não arredondá-lo para "100% certo" — é parte da honestidade metodológica.

---

## Bloco E — Datando a árvore: relógio molecular com TreeTime (20 min)

Transformamos a árvore (em substituições) numa árvore **no tempo**. O TreeTime faz duas coisas:

1. Testa o sinal de relógio com a **regressão raiz-às-pontas**.
2. Estima a data de cada nó, inclusive o **ancestral comum das cepas americanas** (tMRCA).

**Valores esperados:**
- **Taxa de evolução: ~1,2 × 10⁻³ substituições/sítio/ano** — exatamente a ordem de grandeza típica do Zika e de vírus de RNA em geral.
- **R² ≈ 0,75** — boa parte da variação na distância à raiz é explicada pela data de coleta: **há sinal de relógio suficiente** para datar.
- **tMRCA das cepas americanas: ~meados de 2014** (estimativa central ≈ maio de 2014; IC 95% no primeiro semestre de 2014).

**A resposta da prática:** o ancestral comum de todas as cepas americanas data de **~meados de 2014** — **mais de um ano antes** de a epidemia ser oficialmente reconhecida no Brasil, em 2015. Essa é a conclusão dos estudos publicados: **o vírus circulou de forma críptica (não detectada) por cerca de um ano** antes de ser notado. A combinação **genomas + datas + filogenética** enxergou uma circulação que a vigilância tradicional não viu em tempo real.

> **Honestidade metodológica:** 34 genomas são uma amostra pequena; o intervalo de confiança é largo e as datas dependem do modelo e do enraizamento. A mensagem robusta não é a data exata, e sim a **ordem de grandeza**: circulação silenciosa de ~1 ano. **Reportar o intervalo de confiança — e não só o número central — é parte da análise.**

Ao final, salvamos a árvore datada em `.nexus`, pronta para gerar uma figura de artigo no **FigTree** ou no **iTOL**.

---

# MOMENTO 3 — Integração e fechamento (20 min)

## 3.1 Costurando teoria e prática (8 min)

| Na teoria falamos de… | Na prática você viu… |
|---|---|
| Homologia posicional | Genomas de comprimentos diferentes — impossível comparar sem alinhar |
| Alinhamento múltiplo (MSA) | 34 × ~10.800 colunas; quase tudo conservado, ~190 sítios informativos |
| Máxima verossimilhança | IQ-TREE maximizando a verossimilhança sob um modelo explícito |
| Modelo de substituição | ModelFinder escolheu TN+F+I; o `+I` = a maioria das colunas não muda |
| Transições × transversões | A matriz de taxas do modelo, com A↔G e C↔T favorecidas |
| Suporte de ramo (bootstrap) | UFBoot de ~92% no clado americano — alto, mas não certeza |
| Enraizamento por grupo externo | A raiz posta nas cepas asiáticas basais |
| Monofilia | O clado americano contém exatamente as 24 cepas americanas |
| Relógio molecular | Taxa ~1,2×10⁻³ subst./sítio/ano; R² ≈ 0,75 = há sinal temporal |
| Regressão raiz-às-pontas | Pontos subindo: amostras mais recentes, mais distantes da raiz |
| Datação de nós (tMRCA) | Ancestral americano em ~meados de 2014: circulação silenciosa de ~1 ano |

## 3.2 Os sete erros mais comuns (5 min)

1. **Confiar num alinhamento ruim.** Nunca rode a árvore sem **olhar** o alinhamento. Lixo no MSA vira árvore errada de aparência perfeita.
2. **Esquecer de escolher o modelo** (usar JC69 por preguiça). Deixe o ModelFinder decidir.
3. **Reportar a árvore sem o suporte.** Uma topologia sem bootstrap/posterior é uma opinião, não um resultado.
4. **Interpretar bootstrap alto como verdade absoluta** — ou esconder os ramos de suporte baixo.
5. **Enraizar errado (ou não enraizar)** e mesmo assim falar de "ancestral" e de datas. A raiz define o tempo.
6. **Datar sem checar o sinal de relógio.** Se o R² da regressão raiz-às-pontas for ~0, as datas não têm valor — e isso precisa ser dito.
7. **Reportar o tMRCA como um número exato**, sem o intervalo de confiança, e **não registrar** versões de software, semente e modelo. *Métodos com reamostragem são estocásticos: sem `--seed` e sem o modelo anotados, ninguém — nem você — repete o seu resultado.*

## 3.3 O que levar para o artigo final (5 min)

O trabalho final exige **integrar pelo menos duas técnicas** da disciplina. Combinações naturais com esta aula:

- **Filogenia + Sequências (Aula 2):** buscar homólogos da sua proteína/gene com BLAST, alinhar e construir a árvore — identificando quais sequências são realmente parentes e quais são contaminação ou paralogia.
- **Filogenia + Estrutural (Aula 7):** mapear os resíduos **conservados** ao longo da árvore na estrutura 3D — eles tendem a formar o sítio ativo. Conservação evolutiva como pista funcional.
- **Filogenia + Transcriptômica (Aula 4) / Proteômica (Aula 5):** situar um gene/proteína de interesse em sua família evolutiva para discutir função por homologia.
- **Filogenia + Epigenômica (Aula 6):** comparar enzimas modificadoras da cromatina entre espécies e discutir conservação.
- **Filogeografia (como hoje):** rastrear origem e dispersão de um patógeno a partir de genomas públicos datados.

**Checklist de metodologia para o seu artigo:**

- [ ] Origem das sequências (GenBank/outro) + números de acesso + datas/locais
- [ ] Alinhador **e versão**; se houve *trimming*, qual ferramenta e critério
- [ ] Método de inferência (ML/BI) e **programa com versão**
- [ ] **Modelo de substituição** e como foi escolhido (ModelFinder/BIC)
- [ ] Tipo e número de réplicas de **suporte** (UFBoot/bootstrap/posterior)
- [ ] **Semente** (seed) e número de réplicas, para reprodutibilidade
- [ ] **Como a árvore foi enraizada** (grupo externo qual / ponto médio) e justificativa
- [ ] Se houve datação: ferramenta, sinal de relógio (R²), taxa estimada, e **tMRCA com intervalo de confiança**
- [ ] Figura da árvore com escala de ramos, valores de suporte e identificação dos clados

## 3.4 Exercícios (entrega até a próxima aula) — 2 min

1. **(Obrigatório)** Rode o IQ-TREE trocando apenas a semente (`--seed 7` e `--seed 100`). O modelo escolhido muda? O suporte do clado americano muda muito? O que isso diz sobre reportar um único valor de bootstrap?
2. **(Obrigatório)** Troque o grupo externo: enraíze pela cepa mais antiga da Polinésia Francesa (`1_0087_PF`). A data estimada do ancestral americano muda? Discuta como a **escolha da raiz** afeta todas as datas.
3. **(Obrigatório)** Refaça a regressão raiz-às-pontas removendo as 5 cepas mais recentes. O R² melhora ou piora? Por quê? O que isso sugere sobre o papel da amplitude temporal da amostragem?
4. **(Obrigatório, diagnóstico)** Construa, **de propósito**, uma árvore a partir das sequências **não alinhadas** (pule o Bloco B e alimente o IQ-TREE com `sequences.fasta`). Relate o que acontece e por quê. *Este exercício existe para você ver com os próprios olhos a consequência de pular o alinhamento.*
5. **(Desafio)** Pegue uma proteína/gene da **sua linha de pesquisa**, baixe de 8 a 15 sequências homólogas do GenBank (use o BLAST da Aula 2), e refaça todo o pipeline (MAFFT → IQ-TREE → árvore enraizada). Entregue a figura da árvore com os suportes e um parágrafo interpretando os clados.
6. **(Desafio)** Compare, no conjunto do Zika, a árvore de ML (IQ-TREE) com uma árvore de Neighbor-Joining (pode usar o MEGA, clicando). As topologias batem? Onde divergem e por quê?

---

## Glossário rápido

| Termo | Significado |
|---|---|
| **Alinhamento (MSA)** | Arranjo das sequências em colunas de posições homólogas, com lacunas (`-`) |
| **Bootstrap (UFBoot)** | Suporte de ramo por reamostragem das colunas; % de réplicas em que o ramo aparece |
| **Clado** | Um ancestral e **todos** os seus descendentes |
| **Grupo externo (*outgroup*)** | Linhagem sabidamente mais distante, usada para enraizar a árvore |
| **Homologia** | Semelhança por descendência de um ancestral comum |
| **Indel** | Inserção ou deleção (origem das lacunas no alinhamento) |
| **Máxima verossimilhança (ML)** | Método que busca a árvore que torna os dados mais prováveis sob um modelo |
| **ModelFinder** | Rotina que escolhe o modelo de substituição pelo BIC/AIC |
| **Modelo de substituição** | Matriz de taxas de troca entre bases (JC69, HKY, GTR…), com +I/+G |
| **Monofilético** | Grupo que contém um ancestral e todos os seus descendentes, e só eles |
| **Nó interno** | Ancestral comum hipotético (não observado) |
| **Relógio molecular** | Relação entre tempo e número de substituições acumuladas |
| **Root-to-tip** | Regressão distância-à-raiz × data; testa o sinal de relógio (R²) |
| **Taxa de evolução** | Substituições por sítio por ano (inclinação do relógio) |
| **tMRCA** | Tempo do ancestral comum mais recente de um grupo |
| **Topologia** | O padrão de ramificação da árvore (quem é parente de quem) |
| **Transição / transversão** | Troca purina↔purina ou pirim.↔pirim. / purina↔pirimidina |

---

## Recursos e links

**Bases de dados**

- GenBank / NCBI Nucleotide — https://www.ncbi.nlm.nih.gov/nuccore
- TimeTree (tempos de divergência entre espécies) — http://www.timetree.org
- Nextstrain (vigilância genômica em tempo real) — https://nextstrain.org

**Ferramentas (todas livres ou gratuitas para uso acadêmico)**

- MAFFT — https://mafft.cbrc.jp/alignment/software/
- MUSCLE — https://www.drive5.com/muscle/
- IQ-TREE 2 — http://www.iqtree.org
- MEGA — https://www.megasoftware.net
- RAxML — https://cme.h-its.org/exelixis/web/software/raxml/
- MrBayes — https://nbisweden.github.io/MrBayes/
- BEAST 2 — https://www.beast2.org
- Jalview — https://www.jalview.org
- FigTree — http://tree.bio.ed.ac.uk/software/figtree/
- iTOL — https://itol.embl.de
- TreeTime — https://github.com/neherlab/treetime
- Biopython — https://biopython.org
- pyMSAviz — https://github.com/moshi4/pyMSAviz

**Dados da prática**

- Tutorial público Nextstrain/zika — https://github.com/nextstrain/zika-tutorial

**Leituras recomendadas**

- Yang, Z.; Rannala, B. *Molecular phylogenetics: principles and practice*. Nature Reviews Genetics, 2012.
- Nguyen, L.-T. et al. *IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies*. Mol. Biol. Evol., 2015.
- Minh, B. Q. et al. *IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era*. Mol. Biol. Evol., 2020.
- Kalyaanamoorthy, S. et al. *ModelFinder: fast model selection for accurate phylogenetic estimates*. Nature Methods, 2017.
- Hoang, D. T. et al. *UFBoot2: Improving the ultrafast bootstrap approximation*. Mol. Biol. Evol., 2018.
- Katoh, K.; Standley, D. M. *MAFFT multiple sequence alignment software version 7*. Mol. Biol. Evol., 2013.
- Sagulenko, P.; Puller, V.; Neher, R. A. *TreeTime: Maximum-likelihood phylodynamic analysis*. Virus Evolution, 2018.
- Faria, N. R. et al. *Zika virus in the Americas: Early epidemiological and genetic findings*. Science, 2016. (origem genômica do surto americano)
- Metsky, H. C. et al. *Zika virus evolution and spread in the Americas*. Nature, 2017.

---

*Material da disciplina Introdução à Bioinformática — PPGED/IME. Dados e ferramentas de uso público e livre. Todas as células de código da prática foram executadas e verificadas antes da aula.*
