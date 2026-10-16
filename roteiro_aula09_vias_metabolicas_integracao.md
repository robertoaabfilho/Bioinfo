# Aula 9 — Vias Metabólicas e Integração de Dados Ômicos

**Da lista de genes à biologia: enriquecimento funcional, KEGG, Reactome, STRING e metabolômica**

Introdução à Bioinformática — PPGED/IME
Duração: 3 horas (180 min) · Modalidade: online síncrona
Docentes: Prof. Mariana Quezado (IME) e [seu nome]

---

## Antes de começar: o que você precisa (5 min de preparo, faça na véspera)

Nada será instalado na sua máquina. Tudo roda no **navegador**.

| Item | Para quê | Link |
|---|---|---|
| Conta Google | Rodar o Google Colab (Python na nuvem) | [colab.research.google.com](https://colab.research.google.com) |
| Notebook `Aula09_Vias_Integracao_Pratica.ipynb` | Roteiro da prática | Repositório da disciplina, pasta `unidade-09-vias-metabolicas/` |
| Navegador atualizado | Abrir KEGG, Reactome e STRING | — |

> **Você não precisa saber programar.** Nas células de código você apenas clica no ▶ e lê o resultado. O objetivo é **executar, observar e interpretar**.

> **Teste de rede (importante):** a prática conversa com quatro serviços na internet: `rest.kegg.jp`, `biit.cs.ut.ee` (g:Profiler), `reactome.org` e `string-db.org`. Redes corporativas e militares às vezes bloqueiam algum deles. Abra o notebook na véspera e rode a **célula P1**: ela diz "OK" ou "FALHOU" para cada serviço. Se algo falhar, use outra rede (doméstica ou celular roteado).

---

## Mapa da aula

| Momento | Tempo | O que acontece |
|---|---|---|
| **1 — Teoria** | 00:00 – 01:15 (75 min) | Vias como redes; metabolismo e metabolômica; bancos de vias; identificadores; análise de enriquecimento (ORA e GSEA); redes de interação; integração multi‑ômica |
| **Intervalo** | 01:15 – 01:25 (10 min) | Café |
| **2 — Prática** | 01:25 – 02:35 (70 min) | Resposta celular à hipóxia: enriquecimento "na mão" → g:Profiler → Reactome → KEGG com genes **e** metabólitos → rede STRING |
| **3 — Integração e fechamento** | 02:35 – 03:00 (25 min) | Costura teoria↔prática, armadilhas, como usar no artigo final, autoavaliação |

**O fio condutor da aula** é a **hipóxia** — a falta de oxigênio nos tecidos. É um problema operacional direto para a Defesa: pilotos em alta altitude, tropas de montanha, tripulações de submarino, mergulhadores de combate e militares usando máscaras e roupas de proteção QBRN por longos períodos. A célula responde à hipóxia com um "interruptor" molecular chamado **HIF‑1**, que reprograma o metabolismo inteiro. Vamos descobrir essa reprogramação **a partir dos dados**, como faria um pesquisador.

**Onde esta aula se encaixa na disciplina:** nas aulas 4 (transcriptômica), 5 (proteômica) e 6 (epigenômica), cada análise terminou com uma **lista** — genes ou proteínas que mudaram. A pergunta que ficou no ar foi: *"e agora, o que essa lista significa?"* Esta aula responde a essa pergunta e conecta as camadas ômicas entre si.

---

# MOMENTO 1 — Aula teórica (75 min)

## 1.1 O problema da lista (8 min)

Imagine que seu experimento de RNA‑seq (Aula 4) terminou com **350 genes** significativamente aumentados. Ler gene por gene, procurando cada um no Google, levaria dias — e você provavelmente se fixaria nos genes que já conhece.

A **análise funcional** (ou **análise de vias**) troca a pergunta *"o que faz cada gene?"* por *"que **processos** estão super‑representados nessa lista?"*. É como um analista de inteligência que, diante de centenas de interceptações, não lê uma a uma: ele pergunta *quais assuntos aparecem mais do que o normal*.

Três ideias guiam a aula:

1. **Genes e proteínas trabalham em equipe.** Uma enzima sozinha não "faz" glicólise; dez enzimas em sequência fazem.
2. **Equipes conhecidas estão catalogadas** em bancos de dados de vias, construídos a partir de décadas de literatura.
3. **A estatística diz se a coincidência é real.** Encontrar 3 genes da glicólise numa lista de 350 pode ser acaso; encontrar 13 numa lista de 36, não.

## 1.2 O que é uma via biológica (10 min)

Uma **via** (*pathway*) é uma **sequência organizada de eventos moleculares** que produz um resultado. Para um engenheiro, a melhor analogia é uma **linha de produção** ou uma **cadeia logística**: matéria‑prima entra, passa por estações, e sai um produto.

| Tipo de via | O que flui | Exemplo | Analogia |
|---|---|---|---|
| **Metabólica** | Matéria (moléculas pequenas) | Glicólise: glicose → piruvato → lactato | Refinaria |
| **De sinalização** | Informação (sinais químicos) | Via do HIF‑1: falta de O₂ → ativa genes | Cadeia de comando |
| **Regulatória (gênica)** | Controle de expressão | Fator de transcrição liga dezenas de genes | Ordem de operações |

### Vias como grafos

Toda via pode ser representada como um **grafo** (como em teoria de redes):

- **Nós** (pontos): metabólitos, genes, proteínas;
- **Arestas** (linhas): reações, interações, ativações/inibições;
- **Direção**: numa via metabólica, a seta indica o sentido da reação.

Isso permite usar ferramentas de engenharia de redes: **grau** (quantas conexões um nó tem), **hubs** (nós muito conectados — pontos críticos, como um roteador central), **caminhos** e **comunidades**.

### Nomenclatura que vai aparecer

- **Enzima**: proteína que catalisa (acelera) uma reação química.
- **Número EC** (*Enzyme Commission*): "CPF" de uma atividade enzimática, com 4 níveis. Ex.: **EC 1.1.1.27** = lactato desidrogenase. O 1º número diz a classe (1 = oxirredutase, 2 = transferase, 3 = hidrolase…). Útil em **química de proteínas**: proteínas diferentes podem ter o mesmo EC se fazem a mesma reação.
- **Substrato** e **produto**: o que entra e o que sai da reação.
- **Cofator**: molécula auxiliar necessária à reação (ex.: NAD⁺, O₂, 2‑oxoglutarato).

## 1.3 Metabolismo e metabolômica em 10 minutos (10 min)

O **metaboloma** é o conjunto de **moléculas pequenas** (geralmente < 1.500 Da) de uma amostra: açúcares, aminoácidos, ácidos orgânicos, lipídios, nucleotídeos. Na cascata ômica:

> **Genoma** (o que *pode* acontecer) → **Transcriptoma** (o que está *sendo pedido*) → **Proteoma** (as *máquinas* disponíveis) → **Metaboloma** (o que *de fato aconteceu*)

O metaboloma é a camada **mais próxima do fenótipo** — uma "fotografia bioquímica" do instante da coleta. Por isso é muito usado na busca de **biomarcadores** de exposição a agentes químicos, de estresse fisiológico, fadiga e hipóxia.

### Como se mede

| Plataforma | Princípio | Pontos fortes | Limitações |
|---|---|---|---|
| **Espectrometria de massas (MS)**, acoplada a cromatografia líquida (LC‑MS) ou gasosa (GC‑MS) | Separa as moléculas e mede a razão **m/z** (a mesma ideia da Aula 5) | Muito sensível; milhares de sinais | Identificar o composto é difícil |
| **Ressonância magnética nuclear (RMN)** | Mede o comportamento dos núcleos atômicos num campo magnético | Quantitativa, reprodutível, não destrói a amostra | Menos sensível |

Duas estratégias:

- **Metabolômica‑alvo** (*targeted*): mede com precisão uma lista **pré‑definida** de compostos (ex.: 20 metabólitos do metabolismo energético). Ideal para confirmar uma hipótese.
- **Metabolômica não‑alvo** (*untargeted*): tenta capturar **tudo** o que o equipamento enxerga e depois descobre o que mudou. Ideal para gerar hipóteses — mas muitos sinais ficam sem identificação.

### O gargalo: dar nome ao sinal

Um pico de massa 89,024 pode ser lactato… ou outra molécula de mesma massa. A comunidade usa **níveis de confiança na identificação** (iniciativa MSI — *Metabolomics Standards Initiative*): do nível 1 (confirmado com padrão puro) ao nível 4 (sinal desconhecido). **Sempre reporte o nível** num artigo.

### Bancos e ferramentas de metabolômica

| Recurso | O que é | Para que usar |
|---|---|---|
| **HMDB** (*Human Metabolome Database*, Canadá) | Enciclopédia de ~200 mil metabólitos humanos, com concentrações normais em sangue/urina e espectros | Saber o que é um metabólito e se o valor é "normal" |
| **ChEBI** (EMBL‑EBI) | Dicionário curado de entidades químicas de interesse biológico, com ontologia | Identificador químico padronizado |
| **PubChem** (NCBI) | Maior banco aberto de compostos químicos (> 100 milhões) — usado na Aula 7 para ligantes | Estrutura, propriedades, toxicologia |
| **KEGG Compound** | Compostos com número **C** (ex.: C00186 = L‑lactato), ligados aos mapas de vias | Colocar o metabólito no mapa |
| **MetaboAnalyst** (Universidade McGill, web, gratuito) | Plataforma completa: estatística, enriquecimento de metabólitos e análise conjunta genes+metabólitos | Análise sem programar |
| **MetaboLights** (EBI) e **Metabolomics Workbench** (NIH) | Repositórios públicos de dados brutos de metabolômica — o "PRIDE" da metabolômica | Achar dados públicos para reanálise |

## 1.4 Os bancos de vias (15 min)

Cada banco tem uma "filosofia". Conhecê‑la evita interpretações erradas.

### KEGG — *Kyoto Encyclopedia of Genes and Genomes*

- **Quem/onde:** Universidade de Kyoto (Japão), desde 1995.
- **O que tem:** ~550 **mapas de vias desenhados à mão** (metabolismo, processamento de informação genética, sinalização, doenças, fármacos), cobrindo milhares de organismos.
- **Como funciona:** cada mapa de referência (ex.: `map00010`, glicólise) é desenhado uma vez; para cada organismo, o KEGG "acende" as caixas das enzimas que aquele genoma possui (ex.: `hsa00010` para humano, `eco00010` para *E. coli*). Isso é feito pelos **KO** (*KEGG Orthology*): grupos de genes que fazem a mesma função em espécies diferentes.
- **Outras partes:** KEGG **Compound** (C…), **Reaction** (R…), **Enzyme** (EC), **Module** (blocos funcionais menores), **BRITE** (hierarquias funcionais), **Drug** e **Disease**.
- **Ferramentas:** **KEGG Mapper** (web; pinta seus genes/metabólitos nos mapas — usaremos o *Color*) e a **API REST** (`rest.kegg.jp`), que permite consultas automatizadas.
- **Licença:** site e API gratuitos para uso **acadêmico**; o download em massa (FTP) é pago. Por isso muitas ferramentas abertas usam versões antigas do KEGG — fique atento à **versão**.
- **Cuidado:** os mapas são **genéricos** — não dizem em que tecido ou condição cada reação ocorre.

### Reactome

- **Quem/onde:** consórcio público (OICR‑Canadá, EMBL‑EBI, NYU e outros); **aberto** (licença CC0/CC‑BY).
- **Filosofia:** tudo é descrito como **reação** (entradas → saídas, com catalisador e localização na célula), cada uma **revisada por especialistas** e ligada às publicações que a comprovam.
- **Organização:** hierárquica — ex.: *Metabolismo → Metabolismo de carboidratos → Glicólise*.
- **Ferramentas:** *Pathway Browser* (visualização interativa) e **Analysis Service** (envia a lista, recebe o enriquecimento e um *token* — um código que guarda a análise no servidor por alguns dias).
- **Foco:** humano (outras espécies por inferência de ortologia).

### Gene Ontology (GO)

- **O que é:** não é um banco de vias, é uma **ontologia** — um **vocabulário controlado e hierárquico** para descrever o que os genes fazem. Três domínios:
  - **Processo biológico (BP):** o "objetivo" maior (ex.: *resposta à hipóxia*, *processo glicolítico*);
  - **Função molecular (MF):** a atividade bioquímica (ex.: *atividade de lactato desidrogenase*);
  - **Componente celular (CC):** onde atua (ex.: *mitocôndria*).
- **Estrutura:** um grafo dirigido acíclico (DAG) — termos específicos são "filhos" de termos gerais. Consequência prática: **resultados redundantes** (o mesmo grupo de genes aparece em termos‑pai e termos‑filho).
- **Códigos de evidência:** cada anotação diz de onde veio (experimento, inferência por sequência, anotação eletrônica automática — *IEA*). Anotações IEA são as menos confiáveis.

### Outros que você vai encontrar

| Recurso | Característica |
|---|---|
| **WikiPathways** | Vias mantidas pela comunidade, no estilo Wikipédia; aberto |
| **MSigDB** (Broad Institute) | Coleção de **conjuntos de genes** (inclui os 50 *Hallmarks*, como *HALLMARK_HYPOXIA*), usada pelo GSEA |
| **BRENDA** | Maior banco de **dados enzimáticos** (Km, Vmax, pH ótimo, inibidores) — essencial em química de proteínas |
| **Rhea** | Reações bioquímicas curadas, ligadas ao UniProt e ao ChEBI |
| **UniProt** | Já visto: cada proteína traz função, EC, vias e referências cruzadas para KEGG/Reactome |

## 1.5 O problema dos identificadores (5 min)

O mesmo gene tem "RGs" diferentes em cada banco:

| Banco | Exemplo (lactato desidrogenase A) |
|---|---|
| Símbolo HGNC | `LDHA` |
| NCBI Gene (Entrez) | `3939` |
| Ensembl | `ENSG00000134333` |
| UniProt | `P00338` |
| KEGG | `hsa:3939` |
| Reactome | entidade `R-HSA-…` |

O mesmo vale para metabólitos (nome comum, KEGG C…, HMDB, ChEBI, PubChem CID, InChIKey).

**Mapeamento de identificadores** (*ID mapping*) é a tradução entre eles. Problemas comuns:

- **Nomes antigos**: `ERO1A` se chamava `ERO1L`; o Excel transforma `SEPT1` e `MARCH1` em datas (problema tão comum que o HGNC renomeou esses genes em 2020);
- **Um‑para‑muitos**: um gene pode ter várias isoformas proteicas (vários UniProt);
- **Perdas silenciosas**: genes que não são traduzidos simplesmente desaparecem da análise.

**Ferramentas de mapeamento:** **g:Convert** (dentro do g:Profiler), **UniProt ID Mapping** (web), **Ensembl BioMart** (web) e o comando `conv` da API do KEGG. **Regra de ouro: sempre conte quantos identificadores foram mapeados.**

## 1.6 Análise de enriquecimento (17 min)

### ORA — Análise de super‑representação (*Over‑Representation Analysis*)

É o método mais usado. A lógica é a de uma **urna**:

- **N** = total de genes do **universo** (*background*) — as bolinhas da urna;
- **K** = genes do universo que pertencem à via — as bolinhas **vermelhas**;
- **n** = tamanho da sua lista — quantas bolinhas você tirou;
- **k** = genes da lista que estão na via — quantas vermelhas saíram.

Se tirarmos ao acaso, esperamos **n·K/N** vermelhas. O **teste hipergeométrico** (equivalente ao teste exato de Fisher unilateral) calcula a probabilidade de sair **k ou mais** vermelhas por puro acaso — o **valor‑p**.

**Exemplo numérico:** universo de 8.000 genes; a glicólise tem 67; sua lista tem 36. Esperado por acaso: 36 × 67 / 8.000 ≈ **0,3 gene**. Se encontrar **13**, a probabilidade de isso ser acaso é astronomicamente pequena.

#### Três decisões que mudam o resultado

1. **O universo.** Deve conter **todos os genes que poderiam ter entrado na lista**. Num experimento de proteômica que detectou 3.000 proteínas, o universo são essas 3.000 — **não** o genoma inteiro. Usar o genoma inteiro infla artificialmente a significância de proteínas abundantes (que são as que o equipamento detecta).
2. **A correção para testes múltiplos.** Testar 400 vias a 5% gera ~20 falsos positivos esperados. Correções:
   - **Bonferroni:** multiplica o p pelo número de testes — muito rigorosa;
   - **Benjamini‑Hochberg (FDR):** controla a **proporção de falsas descobertas** entre os resultados declarados — o padrão na área (visto nas aulas 4 e 5);
   - **g:SCS:** a correção do g:Profiler, adaptada à hierarquia do GO.
3. **O ponto de corte da lista.** A ORA exige decidir quem entra (ex.: FDR < 0,05 e |log₂FC| > 1). Genes logo abaixo do corte são ignorados — informação perdida.

### GSEA — Análise de enriquecimento de conjuntos de genes (*Gene Set Enrichment Analysis*)

Resolve o problema do ponto de corte: em vez de uma lista, usa **todos os genes ordenados** do mais aumentado ao mais diminuído.

- Percorre a lista de cima para baixo; cada vez que encontra um gene da via, um **escore acumulado** sobe; quando não encontra, desce.
- Se os genes da via se concentram no **topo** (ou no fundo), o escore atinge um pico alto: a via está **coordenadamente** alterada, mesmo que cada gene mude pouco.
- Resultado: **NES** (escore de enriquecimento normalizado) e FDR.
- **Ferramentas:** **GSEA** (Broad Institute, programa desktop gratuito), **fgsea** e **clusterProfiler** (pacotes R/Bioconductor), **GSEApy** (Python).

| | ORA | GSEA |
|---|---|---|
| Entrada | Lista de genes (com corte) | Todos os genes com um valor (ex.: log₂FC) |
| Detecta | Vias com muitos genes muito alterados | Também mudanças pequenas e coordenadas |
| Simplicidade | Alta | Média |
| Exemplo de ferramenta | g:Profiler, Enrichr, DAVID, STRING | GSEA, fgsea, clusterProfiler |

Existem ainda métodos **topológicos** (ex.: SPIA), que consideram a **posição** do gene na via (um gene no início de uma cascata pesa mais). São mais sofisticados e menos usados.

### Ferramentas de enriquecimento (o que cada uma é)

| Ferramenta | Tipo | Comentário |
|---|---|---|
| **g:Profiler** (Univ. Tartu, Estônia) | Web + Python/R, aberto | Várias fontes de uma vez; atualizado com frequência; faz ID mapping sozinho. **Usaremos na prática.** |
| **Enrichr** (Mount Sinai, EUA) | Web + API | Centenas de bibliotecas de conjuntos de genes; muito visual |
| **DAVID** (NIH, EUA) | Web | Clássico e muito citado; agrupa termos redundantes |
| **clusterProfiler** | Pacote R/Bioconductor | Padrão em R para ORA e GSEA, com ótimos gráficos |
| **STRING** | Web + API | Além da rede, faz enriquecimento da lista enviada |
| **Cytoscape** + apps (**ClueGO**, **EnrichmentMap**, **stringApp**) | Programa desktop aberto | Visualiza resultados de enriquecimento como redes de termos |

### Armadilhas clássicas

- **Viés de conhecimento:** genes muito estudados têm muitas anotações e aparecem em tudo; genes pouco estudados nunca enriquecem nada.
- **Circularidade:** se sua lista veio de um conjunto curado (como a nossa, de alvos do HIF), é óbvio que a via do HIF vai aparecer.
- **Vias "de doença" não significam doença:** *"Pathways in cancer"* aparece em quase toda análise porque reúne genes de sinalização genéricos.
- **Enriquecimento ≠ ativação:** a ORA diz que a via está *envolvida*, não se está *ligada* ou *desligada*.
- **Versão importa:** bancos mudam; registre a versão e a data no método do artigo.

## 1.7 Redes de interação proteína‑proteína (5 min)

Já vimos o **STRING** na Aula 5. Relembrando o essencial:

- **O que é:** banco de **associações funcionais** entre proteínas (SIB‑Suíça, Univ. de Copenhague, EMBL), para milhares de organismos.
- **Evidências (canais):** experimentos, bancos curados, mineração de texto, coexpressão, vizinhança genômica, fusão gênica e co‑ocorrência filogenética.
- **Escore de confiança (0–1):** 0,4 média; 0,7 alta; 0,9 altíssima.
- **Cuidado:** associação funcional **≠** contato físico (há a opção de rede só física).

**Cytoscape** (consórcio aberto, desktop): o "AutoCAD" das redes biológicas. Importa redes do STRING (app *stringApp*), colore nós por expressão, detecta módulos e gera figuras de publicação. Em Python, o equivalente programático é o **networkx**, que usaremos na prática.

## 1.8 Integração multi‑ômica e aplicações (5 min)

**Integrar** é combinar camadas para que uma confirme, complemente ou explique a outra. Três níveis, do mais simples ao mais complexo:

1. **Mapeamento conjunto em vias:** pintar genes **e** metabólitos alterados no mesmo mapa (KEGG Mapper, Reactome). Se as duas camadas apontam para a mesma via, a evidência é mais forte. **Faremos isto na prática.**
2. **Análise conjunta de vias:** testar o enriquecimento usando genes e metabólitos juntos (ex.: *Joint Pathway Analysis* do **MetaboAnalyst**).
3. **Integração estatística:** modelos que procuram padrões comuns entre várias matrizes de dados — **MOFA** (fatores latentes; R/Python), **mixOmics/DIABLO** (R). Exigem as mesmas amostras medidas em todas as camadas.

**Cuidado:** correlação entre camadas **não** é causalidade; e camadas têm tempos diferentes (um transcrito muda em minutos, uma proteína em horas).

### Aplicações em química de proteínas e Defesa

- **Agentes neurotóxicos:** a acetilcolinesterase (AChE, Aula 7) está na via da **sinapse colinérgica** (KEGG `hsa04725`). A análise de vias de amostras expostas a organofosforados mostra que processos além da AChE são afetados (estresse oxidativo, metabolismo lipídico) — candidatos a biomarcadores e contramedidas.
- **Hipóxia e desempenho operacional:** perfis metabólicos (lactato, succinato) e genes do HIF como marcadores de aclimatação à altitude.
- **Biodefesa:** comparar as vias metabólicas presentes no genoma de um patógeno (via KEGG por organismo) com as do hospedeiro identifica **enzimas exclusivas do patógeno** — alvos para novos antimicrobianos.
- **Engenharia de proteínas e biocatálise:** números EC, BRENDA e Rhea orientam a escolha de enzimas para **degradar agentes químicos** (ex.: fosfotriesterases que hidrolisam organofosforados).

---

## ☕ INTERVALO (10 min) — 01:15 às 01:25

Aproveite para abrir o Colab e carregar o notebook `Aula09_Vias_Integracao_Pratica.ipynb`.

---

# MOMENTO 2 — Aula prática (70 min)

## O cenário

Você recebeu duas tabelas de um experimento com células expostas a **1% de oxigênio** (hipóxia) comparadas a células em oxigênio normal:

1. **36 genes** aumentados (alvos clássicos do HIF‑1, curados da literatura);
2. **7 metabólitos** medidos por espectrometria de massas, com o log₂ da mudança (valores ilustrativos, coerentes com a fisiologia).

**Missão:** descobrir, a partir dos dados, **que processos a célula reprogramou** e verificar se genes e metabólitos contam a mesma história.

## As ferramentas da prática

| Ferramenta | O que é | Como usamos |
|---|---|---|
| **Google Colab** | Ambiente gratuito do Google para rodar Python no navegador | Executar o notebook |
| **Python** | Linguagem de programação mais usada em ciência de dados | Linguagem do notebook |
| **pandas** | Biblioteca de tabelas ("Excel dentro do Python") | Organizar resultados |
| **scipy** | Biblioteca de estatística e matemática científica | Teste hipergeométrico |
| **matplotlib** | Biblioteca de gráficos | Dot plot e rede |
| **networkx** | Biblioteca de redes/grafos | Rede de proteínas, grau e hubs |
| **requests** | Biblioteca para conversar com APIs | Consultar KEGG, Reactome e STRING |
| **API do KEGG** | Acesso programático ao KEGG | Genes, vias, metabólitos e mapas |
| **g:Profiler** (`gprofiler-official`) | Serviço de enriquecimento + cliente Python | ORA com GO, KEGG, Reactome e WikiPathways |
| **Reactome Analysis Service** | API de análise do Reactome | ORA e link para o navegador gráfico |
| **KEGG Mapper — Color** | Ferramenta web do KEGG | Pintar genes e metabólitos no mapa |
| **API do STRING** | Acesso programático ao STRING | Baixar a rede de interações |

## Passo a passo

Os tempos são aproximados. **Leia a explicação antes de cada célula e responda às perguntas no próprio notebook** (dê duplo clique numa célula de texto para escrever).

### P1. Preparar o ambiente (5 min)
Rode as duas primeiras células. A primeira instala o g:Profiler e carrega as bibliotecas; a segunda testa os quatro serviços.
**Resultado esperado:** quatro linhas "OK".

### P2. Conhecer os dados (5 min)
Veja a lista de genes (agrupada por função apenas para sua orientação) e a tabela de metabólitos.
**Pergunta P2:** antes de rodar a análise, qual processo você aposta que vai aparecer?

### P3. Traduzir identificadores (5 min)
O notebook baixa do KEGG a lista de ~25 mil genes humanos e traduz nossos símbolos para identificadores `hsa:`. Observe que guardamos também os **sinônimos** (nomes antigos).
**Resultado esperado:** 36 de 36 traduzidos.

### P4. Enriquecimento "na mão" (15 min)
É o coração da aula. Você vai:
1. Baixar as tabelas **gene → via** e os nomes das vias;
2. Fazer **um** teste hipergeométrico, para a glicólise, vendo N, K, n e k explicitamente;
3. Repetir para **todas** as vias, aplicar Benjamini‑Hochberg e ver o *dot plot*.

**Resultado esperado:** *Glycolysis / Gluconeogenesis* e *HIF‑1 signaling pathway* no topo, com FDR minúsculo. Também aparecerão vias de metabolismo de carbono, frutose/manose e algumas de doenças.

### P5. Controle negativo (5 min)
Uma lista **aleatória** do mesmo tamanho passa pela mesma análise.
**Resultado esperado:** nenhuma (ou quase nenhuma) via com FDR < 5%, mas algumas com p < 0,05 sem correção — os "falsos alarmes".

### P6. g:Profiler (8 min)
A mesma lista vai ao g:Profiler, contra GO:BP, KEGG, Reactome e WikiPathways.
**Observe:** a concordância com a nossa análise manual, as pequenas diferenças (versões e universos diferentes) e a **redundância** dos termos do GO.

### P7. Reactome (7 min)
A lista vai ao Analysis Service do Reactome. **Abra o link** gerado e explore o *sunburst* e o diagrama da glicólise, que mostra **reações**, não só genes.

### P8. Integração genes + metabólitos (10 min)
1. O notebook descobre em quais vias estão os 7 metabólitos e cruza com o enriquecimento de genes;
2. Mostra o mapa da glicólise do KEGG;
3. Gera o bloco de cores para o **KEGG Mapper — Color**: genes da lista em **laranja**, metabólitos aumentados em **vermelho**, diminuídos em **azul**.

**Como usar o KEGG Mapper — Color:** acesse <https://www.kegg.jp/kegg/mapper/color.html>, cole o bloco impresso pelo notebook na caixa de texto, em *Search against* digite `hsa`, clique em **Exec** e escolha a via *Glycolysis / Gluconeogenesis*.

**Observe:** o ciclo de Krebs tem metabólitos alterados mas **não** aparece no enriquecimento de genes — uma camada vê o que a outra não vê.

### P9. Rede STRING (8 min)
A lista vai ao STRING (confiança ≥ 0,7); o `networkx` monta o grafo, calcula o grau e desenha a rede, com os genes da glicólise em laranja.
**Observe:** o aglomerado denso das enzimas glicolíticas e as proteínas isoladas.

### P10. Salvar (2 min)
Cinco arquivos CSV são salvos e baixados — prontos para o material suplementar de um artigo.

> **Se algo der errado:** (1) confira se rodou todas as células anteriores, na ordem; (2) leia a última linha da mensagem de erro; (3) se for erro de conexão, rode a célula de novo — os servidores públicos às vezes ficam lentos; (4) use *Ambiente de execução → Reiniciar e executar tudo*.

---

# MOMENTO 3 — Conexão entre teoria e prática e fechamento (25 min)

## 3.1 O que cada etapa da prática ilustrou (8 min)

| O que você fez | Conceito teórico | Lição |
|---|---|---|
| Traduzir símbolos para `hsa:` | Mapeamento de identificadores (1.5) | Conte sempre quantos foram mapeados; use sinônimos |
| `hypergeom.sf(k-1, N, K, n)` | ORA e a lógica da urna (1.6) | O valor‑p compara o observado com o esperado por acaso |
| Definir o `UNIVERSO` | Background (1.6) | O universo muda a significância |
| Benjamini‑Hochberg | Testes múltiplos (1.6) | Sem correção, até lista aleatória "dá resultado" |
| Lista aleatória | Controle negativo | Valide o método antes de confiar nele |
| HIF‑1 no topo | Circularidade (1.6) | Listas curadas reencontram as vias de onde vieram |
| g:Profiler × análise manual | Versões e escolhas (1.4, 1.6) | Registre ferramenta, versão, fontes, universo e correção |
| Termos GO redundantes | Hierarquia do GO (1.4) | Agrupe ou escolha termos representativos |
| Reactome por reações | Filosofias dos bancos (1.4) | KEGG = mapas genéricos; Reactome = reações curadas |
| Genes + metabólitos no mapa | Integração por mapeamento (1.8) | Duas camadas concordando = evidência mais forte |
| Krebs só nos metabólitos | Limites dos mapas (1.8) | Uma camada vê o que a outra não vê |
| Rede STRING e hubs | Redes (1.7) | Associação ≠ contato físico; isolado ≠ irrelevante |

## 3.2 A história biológica que os dados contaram (5 min)

Juntando tudo: sem oxigênio, a mitocôndria não consegue produzir energia com eficiência. O HIF‑1 é estabilizado (porque as enzimas **EGLN**, que precisam de O₂ e 2‑oxoglutarato para destruí‑lo, param de funcionar) e liga dezenas de genes. A célula:

1. **Aumenta a entrada de glicose** (transportadores SLC2A1/SLC2A3) — por isso a glicose cai;
2. **Acelera a glicólise** (HK2, PFKP, ALDOA, PGK1, ENO1, PKM…) — a via enriquecida nº 1;
3. **Desvia o piruvato para lactato** (LDHA) e **bloqueia a entrada na mitocôndria** (PDK1) — por isso o lactato sobe;
4. **Controla a acidez** (CA9) e **pede mais vasos sanguíneos** (VEGFA, ADM) e **mais hemácias** (EPO).

Esse é o mesmo princípio fisiológico por trás do **aumento do lactato no sangue** em esforço físico intenso e em altitude — e um exemplo de como uma análise de vias transforma uma lista em um **mecanismo testável**.

## 3.3 Armadilhas mais comuns em artigos (4 min)

- Usar o genoma inteiro como universo em dados de proteômica ou metabolômica;
- Relatar valores‑p sem correção;
- Não informar a versão do banco nem a data da análise;
- Tratar *"Pathways in cancer"* ou outras vias de doença como achado biológico;
- Interpretar uma aresta do STRING como interação física comprovada;
- Concluir "a via está ativada" a partir de ORA (que não tem direção).

## 3.4 Ideias para o trabalho final (3 min)

O artigo final deve integrar **pelo menos duas técnicas** da disciplina. Exemplos a partir desta aula:

- Reanalisar um conjunto público de proteômica (PRIDE, Aula 5) ou RNA‑seq (GEO, Aula 4) **+** enriquecimento de vias e rede STRING;
- Escolher uma enzima de interesse (ex.: uma fosfotriesterase ou a AChE), localizá‑la nas vias KEGG/Reactome **+** modelar sua estrutura (Aula 7) **+** filogenia da família (Aula 8);
- Reanalisar um estudo de metabolômica do **MetaboLights** ou **Metabolomics Workbench** **+** mapeamento conjunto com genes de um estudo de transcriptômica na mesma condição.

## 3.5 Fechamento (5 min)

**Cinco ideias para levar:**

1. Uma lista de genes só ganha sentido biológico quando é **projetada em vias, ontologias e redes**.
2. A **ORA** é um teste de urna: compare o **observado** com o **esperado** — e o esperado depende do **universo**.
3. **Sempre corrija para testes múltiplos** e **faça um controle negativo** quando possível.
4. Bancos têm **filosofias diferentes** (KEGG = mapas genéricos; Reactome = reações curadas; GO = vocabulário hierárquico; STRING = associações com evidência). Use mais de um e **registre as versões**.
5. **Integrar camadas ômicas** fortalece conclusões — e revela o que cada camada sozinha não enxerga.

**Perguntas de autoavaliação (responda em 2–3 linhas cada):**

1. Explique a análise de super‑representação usando a analogia da urna, para um colega que nunca estudou estatística.
2. Seu experimento de proteômica detectou 2.500 proteínas e 120 mudaram. Qual universo você usa na ORA? Por quê?
3. Qual a diferença essencial entre ORA e GSEA? Em que situação a GSEA detecta algo que a ORA perde?
4. Por que o termo *"Pathways in cancer"* aparece em tantas análises, mesmo sem câncer envolvido?
5. Dê um exemplo de como dados de metabolômica podem confirmar — ou contradizer — uma análise de vias feita com genes.
6. Num estudo sobre exposição a organofosforados, em que via KEGG você esperaria encontrar a AChE? Que outra camada ômica você mediria para confirmar o efeito, e por quê?

**Exercícios assíncronos:** os cinco exercícios no final do notebook fazem parte da lista da unidade.

---

## Glossário

| Termo | Significado |
|---|---|
| **API** | Interface que permite a um programa pedir dados a um serviço (o "garçom" da Aula 1) |
| **Background / universo** | Conjunto de todos os genes que poderiam ter entrado na lista |
| **Cofator** | Molécula auxiliar necessária a uma reação enzimática |
| **EC (número)** | Código da atividade de uma enzima (ex.: EC 1.1.1.27) |
| **Enriquecimento** | Presença de um grupo de genes na lista acima do esperado por acaso |
| **FDR** | Taxa de falsas descobertas; valor‑p ajustado para testes múltiplos |
| **Gene Ontology (GO)** | Vocabulário hierárquico que descreve processos, funções e localizações |
| **GSEA** | Enriquecimento usando a lista completa ordenada de genes |
| **HIF‑1** | Fator de transcrição ativado pela falta de oxigênio |
| **Hipóxia** | Oferta insuficiente de oxigênio aos tecidos |
| **Hub** | Nó muito conectado numa rede |
| **ID mapping** | Tradução de identificadores entre bancos |
| **KO (KEGG Orthology)** | Grupo de genes de espécies diferentes com a mesma função |
| **Metaboloma** | Conjunto de moléculas pequenas de uma amostra |
| **NES** | Escore de enriquecimento normalizado do GSEA |
| **Ontologia** | Vocabulário controlado com relações hierárquicas entre termos |
| **ORA** | Análise de super‑representação (teste da urna) |
| **Teste hipergeométrico** | Teste estatístico de sorteio sem reposição, base da ORA |
| **Token (Reactome)** | Código que guarda temporariamente sua análise no servidor |
| **Via (pathway)** | Sequência organizada de eventos moleculares com um resultado |

## Referências e links

**Bancos e ferramentas**
- KEGG — <https://www.kegg.jp> · API: <https://www.kegg.jp/kegg/rest/keggapi.html> · Mapper: <https://www.kegg.jp/kegg/mapper/>
- Reactome — <https://reactome.org>
- Gene Ontology — <https://geneontology.org>
- g:Profiler — <https://biit.cs.ut.ee/gprofiler>
- STRING — <https://string-db.org>
- WikiPathways — <https://www.wikipathways.org>
- MSigDB / GSEA — <https://www.gsea-msigdb.org>
- Enrichr — <https://maayanlab.cloud/Enrichr>
- DAVID — <https://davidbioinformatics.nih.gov>
- Cytoscape — <https://cytoscape.org>
- MetaboAnalyst — <https://www.metaboanalyst.ca>
- HMDB — <https://hmdb.ca> · ChEBI — <https://www.ebi.ac.uk/chebi> · PubChem — <https://pubchem.ncbi.nlm.nih.gov>
- BRENDA — <https://www.brenda-enzymes.org> · Rhea — <https://www.rhea-db.org>
- MetaboLights — <https://www.ebi.ac.uk/metabolights> · Metabolomics Workbench — <https://www.metabolomicsworkbench.org>

**Leituras**
- Kanehisa, M. et al. KEGG: biological systems database as a model of the real world. *Nucleic Acids Research* (edição anual de bancos de dados).
- Milacic, M. et al. The Reactome Pathway Knowledgebase. *Nucleic Acids Research* (edição anual de bancos de dados).
- Kolberg, L. et al. g:Profiler — interoperable web service for functional enrichment analysis and gene identifier mapping. *Nucleic Acids Research*, 2023.
- Szklarczyk, D. et al. The STRING database. *Nucleic Acids Research* (edição anual de bancos de dados).
- Subramanian, A. et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*, 2005.
- Khatri, P.; Sirota, M.; Butte, A. J. Ten years of pathway analysis: current approaches and outstanding challenges. *PLoS Computational Biology*, 2012.
- Timmons, J. A.; Szkop, K. J.; Gallagher, I. J. Multiple sources of bias confound functional enrichment analysis of global -omics data. *Genome Biology*, 2015.
- Semenza, G. L. Hypoxia-inducible factors in physiology and medicine. *Cell*, 2012.
- Pevsner, J. *Bioinformatics and Functional Genomics* (capítulos de análise funcional e vias).
