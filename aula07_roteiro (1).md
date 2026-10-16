# Aula 7 — Bioinformática Estrutural

**Visualização molecular · Modelagem por homologia · Docking · Dinâmica molecular**

Introdução à Bioinformática — PPGED/IME
Duração: 3 horas (180 min) · Modalidade: online síncrona
Docentes: Prof. Mariana Quezado (IME) e [seu nome]

> **Nota sobre a carga horária:** a ementa reserva 4 h para a Unidade 7. Este roteiro cobre 3 h de encontro síncrono; a hora restante está distribuída nos exercícios assíncronos da seção 3.4, que são parte obrigatória da unidade.

---

## Antes de começar: o que você precisa (5 min de preparo, faça na véspera)

Nada será instalado na sua máquina. Tudo roda no **navegador**.

| Item | Para quê | Link |
|---|---|---|
| Conta Google | Rodar o Google Colab (Python na nuvem) | [colab.research.google.com](https://colab.research.google.com) |
| Navegador atualizado (Chrome/Edge/Firefox) | Visualizador 3D do PDB | — |
| (Opcional) ChimeraX instalado | Visualização com mais recursos | [cgl.ucsf.edu/chimerax](https://www.cgl.ucsf.edu/chimerax/download.html) |

> **Você não precisa saber programar.** Nas células de código do Colab você vai apenas clicar no ▶ e ler o que sai. O objetivo é **executar, observar e interpretar** — exatamente como nas aulas anteriores.

> **O caderno da prática.** Todo o código dos Blocos B, C e D está no notebook **`Aula07_Estrutural_Pratica.ipynb`**, comentado linha a linha e já testado. Faça o *upload* dele no Colab (*Arquivo → Fazer upload de notebook*) e execute as células em ordem. Este roteiro explica **por que** cada passo existe; o notebook **executa** cada passo. Os trechos de código que aparecem adiante são os mesmos do notebook, reproduzidos aqui para leitura e para a prova.

> **Se a rede da sua instituição bloquear downloads automáticos:** algumas redes corporativas e militares bloqueiam `files.rcsb.org`. Se a célula de download falhar, baixe o arquivo `4EY7.pdb` manualmente pelo navegador (em [rcsb.org/structure/4EY7](https://www.rcsb.org/structure/4EY7) → *Download Files* → *PDB Format*) e arraste-o para o painel de arquivos do Colab (ícone de pasta na barra lateral esquerda). O resto do notebook funciona igual. **Teste isso na véspera, não durante a aula.**

---

## Mapa da aula

| Momento | Tempo | O que acontece |
|---|---|---|
| **1 — Teoria** | 00:00 – 01:15 (75 min) | Por que estrutura importa; o PDB; visualização; predição de estrutura (homologia e AlphaFold); docking; dinâmica molecular |
| **Intervalo** | 01:15 – 01:25 (10 min) | Café |
| **2 — Prática** | 01:25 – 02:40 (75 min) | Hands-on com dados públicos e software livre: PDB → modelo → docking → análise de trajetória |
| **3 — Integração e fechamento** | 02:40 – 03:00 (20 min) | Costura teoria↔prática, erros comuns, como usar isso no artigo final, exercícios |

**O fio condutor da aula** é um alvo de interesse direto para a Engenharia de Defesa: a **acetilcolinesterase (AChE)**, a enzima inibida por agentes neurotóxicos organofosforados (sarin, VX) e alvo de fármacos e antídotos. Usaremos sempre dados públicos e ferramentas abertas.

---

# MOMENTO 1 — Aula teórica (75 min)

## 1.1 Por que a estrutura 3D importa (10 min)

Nas aulas 2 a 6 trabalhamos com **sequências** (texto: ATGC…, ou a sequência de aminoácidos) e com **níveis de expressão** (números). Mas uma proteína não funciona como um texto: ela funciona como uma **peça mecânica tridimensional**.

A ideia central:

> **Sequência → Estrutura → Função.**
> A sequência de aminoácidos determina como a cadeia se dobra no espaço; o formato dobrado cria bolsões, fendas e superfícies; e é a geometria + a química desses bolsões que define o que a proteína faz e com quem ela interage.

Analogia de engenharia: saber a lista de peças de um motor (sequência) não permite prever o torque. É preciso saber **como as peças estão montadas** (estrutura).

### O que a bioinformática estrutural responde

1. **Onde** uma molécula pequena se liga a uma proteína (o "sítio ativo")?
2. **Quão forte** é essa ligação (afinidade)?
3. **Por que** uma mutação quebra a função de uma proteína?
4. Como a proteína **se move** — ela é rígida ou respira, abre e fecha?

### Por que isso interessa à Engenharia de Defesa

| Problema de defesa | Pergunta estrutural |
|---|---|
| Agentes neurotóxicos (sarin, VX, novichoks) | Como o organofosforado ataca a serina catalítica da acetilcolinesterase? |
| Antídotos e reativadores (pralidoxima, obidoxima, HI-6) | Que molécula consegue entrar na garganta estreita da enzima e remover o fosfato? |
| Profilaxia e *bioscavengers* | A butirilcolinesterase pode "sequestrar" o agente antes que ele atinja o alvo? |
| Detecção e biossensores | Que anticorpo ou aptâmero reconhece a toxina com especificidade? |
| Biodefesa (toxinas: ricina, botulínica) | Onde está o sítio catalítico e que inibidor o bloqueia? |
| Materiais e enzimas de descontaminação | Como melhorar por engenharia uma enzima que degrada organofosforados (fosfotriesterase)? |

> **Escopo desta aula.** Trabalhamos do lado **defensivo e terapêutico** do problema: entender o mecanismo de inibição para desenvolver antídotos, profilaxia, detecção e descontaminação. É também esse o enquadramento esperado no seu artigo final.

### Nosso alvo de hoje: a acetilcolinesterase (AChE)

- **Função fisiológica:** degrada o neurotransmissor acetilcolina na fenda sináptica — é o "freio" do impulso nervoso. Está entre as enzimas mais rápidas conhecidas.
- **Se for inibida:** a acetilcolina acumula, o músculo não relaxa → crise colinérgica, paralisia respiratória. É exatamente esse o mecanismo de letalidade dos agentes nervosos e dos inseticidas organofosforados.
- **Se for inibida "na medida certa":** é terapia para Alzheimer — é assim que funciona a **donepezila** (Aricept®).
- **Detalhe estrutural que faz toda a diferença:** o sítio ativo fica no fundo de uma **garganta estreita de ~20 Å** (*aromatic gorge*), forrada de aminoácidos aromáticos. Qualquer molécula candidata precisa caber e descer por esse túnel. Esse é um problema de **geometria** — logo, um problema de bioinformática estrutural.

---

## 1.2 De onde vêm as estruturas: o PDB (12 min)

### O Protein Data Bank (PDB)

O **PDB** ([rcsb.org](https://www.rcsb.org)) é o repositório público mundial de estruturas 3D de macromoléculas biológicas. Criado em 1971, hoje contém mais de 200 mil estruturas depositadas. É o equivalente estrutural do que o GenBank é para sequências (Aula 2).

Cada estrutura tem um **código de 4 caracteres** (ex.: `4EY7`, `1EVE`, `1AO6`). É esse código que você cita em um artigo.

### Como as estruturas são obtidas (e por que isso muda sua confiança nelas)

| Técnica | O que é | Limite prático |
|---|---|---|
| **Cristalografia de raios X** | A proteína é cristalizada; o padrão de difração dos raios X é convertido em um mapa de densidade eletrônica, onde os átomos são "encaixados" | ~85–90% do PDB. Exige cristalizar — difícil para proteínas de membrana. A proteína está "congelada" em uma conformação |
| **Criomicroscopia eletrônica (crio-EM)** | Milhares de imagens de partículas individuais congeladas são combinadas em um mapa 3D | Revolucionou complexos grandes (ribossomos, vírus). Resolução hoje chega a nível atômico |
| **RMN (NMR)** | Resolve a estrutura em solução; entrega um **conjunto** de modelos (ex.: 20 modelos) | Proteínas pequenas (< ~30 kDa). Vantagem: mostra flexibilidade |

### Três números que você SEMPRE deve olhar antes de usar uma estrutura

1. **Resolução (em Å)** — só existe para raios X/crio-EM. Quanto **menor**, melhor.
   - < 1,5 Å: excelente, dá para ver átomos individualmente
   - 1,5 – 2,5 Å: muito boa — padrão para docking
   - 2,5 – 3,5 Å: utilizável, mas cadeias laterais podem estar mal definidas
   - \> 3,5 Å: cuidado — não use para docking sem justificativa
2. **R-free** — mede o quanto o modelo concorda com dados que **não** foram usados no ajuste. Valor típico bom: < 0,25. É a "validação cruzada" da cristalografia.
3. **B-factor (fator de temperatura)** — um número por átomo que indica o quanto ele "treme" ou está mal definido. B alto = região flexível ou pouco confiável. É possível colorir a estrutura por B-factor no visualizador.

> **Armadilha clássica:** aluno baixa uma estrutura de 3,8 Å, faz docking, reporta afinidade com 3 casas decimais. O dado de entrada não sustenta essa precisão.

### Um quarto número, que quase todo mundo esquece: **quantas cópias da proteína há no arquivo?**

Um arquivo do PDB contém a **unidade assimétrica** do cristal — que pode ter 1, 2, 4 ou mais cópias da mesma proteína, identificadas por **cadeias** (*chains*): A, B, C…

A nossa estrutura de hoje, a **4EY7**, tem **duas cópias da AChE: cadeias A e B** (mais duas cadeias de glicano, C e D). Isso tem uma consequência prática imediata e silenciosa:

> **Se você extrair "o ligante donepezila" do arquivo sem filtrar a cadeia, você extrai DUAS moléculas de donepezila — uma na cadeia A e outra na cadeia B — e o "centro" que você calcular cairá no meio do caminho entre os dois sítios ativos, provavelmente dentro do solvente ou do miolo da proteína.** O docking roda, não dá erro nenhum, e o resultado é lixo.

Esse é o tipo de erro que não aparece como mensagem de erro, e sim como um número plausível e errado. **Na prática de hoje você vai filtrar explicitamente a cadeia A** — e vamos conferir o resultado.

### Formatos de arquivo que você vai encontrar

| Extensão | O que é |
|---|---|
| `.pdb` | Formato clássico, texto de colunas de largura fixa. Legível a olho nu. Limitado (máx. ~99.999 átomos) |
| `.cif` / `.mmCIF` | Formato moderno e **padrão oficial** desde 2014. Sem limite de tamanho, mais metadados |
| `.pdbqt` | Formato do AutoDock/Vina: é um PDB **+ cargas parciais (Q) + tipos de átomo (T) + ligações rotacionáveis**. É o que o docking consome |
| `.sdf` / `.mol2` | Formatos de moléculas pequenas (ligantes), vindos do PubChem/ChEMBL |

Anatomia mínima de uma linha de um arquivo `.pdb`:

```
ATOM    145  CA  SER A  203      12.345  67.890  10.111  1.00 18.52           C
 |       |   |   |  |   |          |       |       |      |     |             |
tipo  n°át. átomo res cad n°res     x       y       z   ocup. B-factor      elemento
```

- `ATOM` = átomo da proteína · `HETATM` = tudo que não é proteína/ácido nucleico (ligante, água, íon)
- `HOH` = água · `SER 203` = a serina catalítica da AChE humana
- **A posição das colunas é fixa.** O nome do resíduo está nas colunas 18–20 e a cadeia na coluna 22. Guarde isso: é por isso que, na prática, vamos filtrar o arquivo **por posição de coluna** e não por "procurar o texto em qualquer lugar da linha".

### Bases de dados complementares

- **PubChem** ([pubchem.ncbi.nlm.nih.gov](https://pubchem.ncbi.nlm.nih.gov)) — >100 milhões de compostos químicos. É daqui que tiramos o ligante (em `.sdf` ou como texto SMILES).
- **ChEMBL** ([ebi.ac.uk/chembl](https://www.ebi.ac.uk/chembl)) — compostos **com dados de bioatividade medidos** (IC50, Ki). Essencial para validar um docking contra a realidade experimental.
- **UniProt** (Aula 2) — a ficha funcional da proteína; traz links para PDB e AlphaFold. A AChE humana é a **P22303**.
- **AlphaFold DB** ([alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk)) — >200 milhões de estruturas **previstas** computacionalmente.

---

## 1.3 Visualização molecular (8 min)

Visualizar não é enfeite: é a etapa de **controle de qualidade** e de **formulação de hipóteses**. Erros grosseiros (cadeia faltando, ligante no lugar errado, duas cópias da proteína — exatamente o caso da 4EY7) aparecem em 10 segundos na tela e passariam despercebidos no arquivo de texto.

### As três ferramentas

| Ferramenta | O que é | Quando usar |
|---|---|---|
| **Mol\*** (MolStar) | Visualizador que roda **dentro do navegador**; é o visualizador oficial do site do RCSB PDB. Código aberto. Zero instalação | Inspeção rápida, aula, apresentação. **É o que usaremos na prática** |
| **PyMOL** | O mais usado em publicações. Tem versão de código aberto (*open-source PyMOL*) e versão comercial. Controlado por uma **linha de comando própria** → reprodutível e automatizável | Figuras de artigo, scripts repetíveis |
| **UCSF ChimeraX** | Sucessor do Chimera, desenvolvido na UCSF. Gratuito para uso acadêmico, excelente com crio-EM, mapas de densidade e estruturas enormes | Análises mais pesadas, vídeos, crio-EM |

### Vocabulário de representações (o que significa cada "estilo")

- **Cartoon / Ribbon** — mostra só o esqueleto e o enovelamento (hélices α como espirais, folhas β como setas). É o mapa geral.
- **Sticks (bastões)** — mostra as ligações químicas. Usado para o ligante e para os resíduos do sítio ativo.
- **Spheres / Surface (superfície)** — mostra o volume ocupado. **A superfície é indispensável para ver uma cavidade** — é assim que a "garganta" da AChE fica visível.
- **Coloração** — por cadeia, por tipo de estrutura secundária, por B-factor, por hidrofobicidade (apolar = o que atrai um fármaco gorduroso).

Três comandos de PyMOL que resolvem 80% dos casos (guarde para depois):

```
fetch 4EY7                                            # baixa a estrutura direto do PDB
show surface                                          # mostra a superfície molecular
select sitio, byres (polymer within 5 of resn E20)    # resíduos a 5 Å do ligante E20
```

---

## 1.4 Quando não existe estrutura: modelagem e predição (18 min)

Problema real: existem centenas de milhões de sequências de proteínas conhecidas, mas só ~200 mil estruturas experimentais. **Para a maioria das proteínas, não há estrutura resolvida.** O que fazer?

### A) Modelagem por homologia (ou modelagem comparativa)

**Premissa:** a estrutura é muito mais conservada pela evolução do que a sequência. Se duas proteínas têm sequências parecidas, elas se dobram de maneira parecida.

Os quatro passos:

1. **Buscar molde (*template*):** rodar BLAST (Aula 2) da sua sequência contra as sequências do PDB. Encontrar uma proteína homóloga **que já tenha estrutura**.
2. **Alinhar** sua sequência (alvo) com a do molde — é aqui que os erros nascem. Um alinhamento deslocado produz um modelo errado com aparência perfeita.
3. **Construir** o modelo: copiar as coordenadas do esqueleto do molde e reconstruir as alças (*loops*) e as cadeias laterais diferentes.
4. **Validar** o modelo (ver abaixo).

**Regra prática da identidade de sequência entre alvo e molde:**

| Identidade | Qualidade esperada |
|---|---|
| > 50% | Modelo muito bom, comparável a uma estrutura de média resolução |
| 30 – 50% | Modelo útil, com erros em alças e cadeias laterais |
| < 30% | "Zona crepuscular" (*twilight zone*) — alto risco de modelo errado |

**Ferramentas:**

- **SWISS-MODEL** ([swissmodel.expasy.org](https://swissmodel.expasy.org)) — servidor web do Instituto Suíço de Bioinformática. Você cola a sequência, ele faz tudo automaticamente e devolve o modelo com métricas de qualidade. Gratuito, sem instalação.
- **MODELLER** — programa clássico (linha de comando, Python), mais controle, licença acadêmica gratuita. É o motor por trás de muitos servidores.

**Métricas de qualidade do SWISS-MODEL:**

- **GMQE** (0 a 1): estimativa global de qualidade; quanto mais perto de 1, melhor.
- **QMEANDisCo** (0 a 1): qualidade local, resíduo a resíduo. Abaixo de 0,6 → trecho pouco confiável.
- **Gráfico de Ramachandran:** mostra se os ângulos da cadeia principal são geometricamente permitidos. Muitos resíduos em região proibida = modelo com problema físico.

### B) Predição por aprendizado profundo: AlphaFold

Em 2020-2021 o **AlphaFold 2** (DeepMind) mudou o campo: passou a prever estruturas com acurácia frequentemente comparável à experimental, **sem precisar de um molde óbvio**. Em 2024 o Prêmio Nobel de Química foi atribuído a trabalhos nessa área (David Baker; Demis Hassabis e John Jumper).

Como funciona, em uma frase acessível: o método monta um **alinhamento múltiplo de milhares de sequências homólogas** de várias espécies (um MSA, como o da Aula 2) e aprende que **pares de posições que mutam juntos ao longo da evolução provavelmente estão próximos no espaço 3D** — essas restrições de contato são então convertidas em coordenadas por uma rede neural.

**Como usar na prática:**

- **AlphaFold DB** — para a maioria das proteínas de organismos modelo (incluindo todo o proteoma humano), **o modelo já está pronto e é só baixar**. Instantâneo.
- **ColabFold** ([github.com/sokrypton/ColabFold](https://github.com/sokrypton/ColabFold)) — roda AlphaFold2 de graça no Google Colab, em minutos. É a porta de entrada prática para sequências novas (ex.: uma proteína que você sequenciou).
- **Boltz / Chai-1 / AlphaFold 3** — geração seguinte, prevê também complexos com ligantes, DNA e RNA. Boltz e Chai são **abertos**; o AlphaFold 3 tem uso restrito via servidor.

**As duas métricas que você precisa saber ler:**

| Métrica | O que mede | Como interpretar |
|---|---|---|
| **pLDDT** (0–100), por resíduo | Confiança **local** | > 90 muito alta · 70–90 boa (esqueleto confiável) · 50–70 baixa · **< 50 provavelmente região desordenada** (não é "erro", é proteína sem forma fixa) |
| **PAE** (Å), matriz entre pares | Confiança na **posição relativa** entre dois trechos | Blocos de PAE baixo = domínios rígidos bem posicionados; PAE alto entre dois domínios = o modelo não sabe a orientação relativa entre eles |

Nos modelos do AlphaFold, o campo do **B-factor do arquivo PDB é reaproveitado para guardar o pLDDT** — por isso colorir por B-factor mostra diretamente a confiança (azul = confiável, laranja = não). **Cuidado:** isso significa que, num modelo do AlphaFold, o "B-factor" **não** é mobilidade atômica. Se você aplicar a interpretação de cristalografia a um modelo previsto, vai concluir o oposto do correto (lá, valor **alto** = confiável).

**Limitações honestas (cobradas na avaliação):**

- Prevê **uma** conformação, normalmente o estado mais populado — não prevê o movimento nem estados alternativos.
- Não prevê, por si só, ligantes, íons, cofatores ou modificações pós-traducionais.
- Não prevê o efeito de uma mutação pontual na estabilidade de forma confiável.
- Sítios ativos previstos podem ter cadeias laterais em rotâmeros errados — o que **degrada o docking**.

> **Decisão prática:** existe estrutura experimental de boa resolução, com ligante parecido com o seu, do organismo certo? **Use a experimental.** Não existe? Use AlphaFold/homologia — e declare isso claramente na metodologia do seu artigo.

---

## 1.5 Docking molecular (18 min)

### O que é

**Docking** responde a duas perguntas sobre um par proteína–molécula pequena:

1. **Pose:** em que posição e orientação o ligante se encaixa no sítio?
2. **Score:** qual a estimativa de energia dessa interação (em kcal/mol, valores negativos = favorável)?

Analogia: é um problema de **encaixe de peças com tolerância e com química** — como encontrar a chave certa para uma fechadura, sabendo que tanto a chave quanto a fechadura são um pouco flexíveis.

### As duas engrenagens de qualquer programa de docking

1. **Algoritmo de busca** — explora o espaço de possibilidades: 3 translações + 3 rotações + 1 ângulo para cada **ligação rotacionável** do ligante. A donepezila, por exemplo, tem **6 ligações rotacionáveis** — ou seja, um espaço de busca de 12 dimensões. Por isso usa-se busca estocástica (Monte Carlo, algoritmos genéticos), **não** busca exaustiva.
2. **Função de pontuação (*scoring function*)** — estima a energia de cada pose somando termos: repulsão estérica, atração de van der Waals, ligações de hidrogênio, efeito hidrofóbico, eletrostática e penalidade por perda de entropia conformacional.

### Ferramentas abertas

| Ferramenta | Característica |
|---|---|
| **AutoDock Vina** | O mais usado do mundo. Rápido, código aberto (Apache 2.0), precisão boa. **É o que usaremos.** Versão atual: 1.2.x |
| **AutoDock4 / AutoDockTools (ADT)** | Geração anterior, mais configurável, mais lento. ADT é a interface gráfica clássica de preparação |
| **Smina** | *Fork* do Vina que permite customizar a função de pontuação e minimizar poses |
| **GNINA** | *Fork* do Smina com **rescoring por rede neural convolucional** — costuma melhorar a seleção de poses |
| **PyRx** | Interface gráfica que empacota o Vina para triagem virtual em lote |
| **Open Babel** | "Canivete suíço" de conversão entre formatos químicos; adiciona hidrogênios e calcula cargas |
| **Meeko / RDKit** | Bibliotecas Python para preparar ligantes corretamente (torções, cargas, protonação). **Meeko é a preparação de ligante recomendada hoje pelos próprios autores do Vina** — melhor que o Open Babel para esse papel |

### O fluxo de trabalho (decore esta sequência)

```
1. Receptor   : baixar do PDB → escolher UMA cadeia → remover águas e ligantes
                → adicionar hidrogênios → gerar .pdbqt (com cargas parciais)
2. Ligante    : baixar do PubChem (.sdf/SMILES) → gerar conformação 3D → protonação em pH 7,4
                → definir ligações rotacionáveis → gerar .pdbqt
3. Caixa      : definir center_x/y/z e size_x/y/z (a "caixa de busca")
4. Rodar Vina : exhaustiveness, num_modes, energy_range, seed
5. Validar    : redocking → RMSD
6. Analisar   : visualizar a pose, mapear as interações
```

### Os parâmetros que você vai ajustar

- **A caixa (*grid box*)** — define *onde* o programa procura.
  - **Docking dirigido (*site-specific*):** caixa pequena (~20–25 Å de lado) centrada no sítio ativo conhecido. É o padrão quando já se sabe onde o ligante liga. Mais rápido e mais confiável.
  - **Docking cego (*blind*):** caixa cobrindo a proteína inteira. Usado quando o sítio é desconhecido. Muito mais ruidoso.
  - **Como achar o centro:** se a estrutura cristalográfica já vem com um ligante no sítio, basta calcular o **centro geométrico desse ligante** — **de uma única cadeia** (ver 1.2). É o que faremos.
- **`exhaustiveness`** — quantas vezes a busca é repetida de forma independente. Padrão 8. Valores maiores (16, 32) dão resultados mais reprodutíveis e são mais lentos. **Docking é estocástico: rodar duas vezes pode dar valores ligeiramente diferentes.**
- **`num_modes`** — quantas poses retornar (padrão 9).
- **`energy_range`** — diferença máxima de energia, em kcal/mol, entre a melhor pose e as demais reportadas.
- **`seed`** — a semente do gerador de números aleatórios. Fixá-la é o que torna seu resultado **reproduzível por outra pessoa**.

### Validação: redocking e RMSD (a parte que todo mundo pula e o avaliador cobra)

**Redocking:** pegue a estrutura cristalográfica que já contém um ligante, **retire** esse ligante, e peça ao programa que o recoloque. Depois compare a pose prevista com a pose real.

A métrica é o **RMSD** (*Root Mean Square Deviation*, desvio quadrático médio) entre os átomos pesados:

- **RMSD < 2,0 Å → o protocolo reproduziu a realidade. Protocolo validado.**
- RMSD 2–3 Å → parcialmente correto (geralmente orientação certa, detalhes errados)
- RMSD > 3 Å → o protocolo não serve para esse alvo; revise preparação, caixa e protonação

Sem esse passo, um score de −9,5 kcal/mol não significa nada.

> **Atenção à armadilha lógica do redocking:** ele valida o **protocolo**, não o **score**. Um redocking com RMSD de 0,3 Å prova que seu pipeline encontra a pose certa quando ela existe; não prova que o número em kcal/mol corresponde à afinidade medida. São duas afirmações diferentes, e confundi-las é o erro mais comum em trabalhos de disciplina.

### Limitações — leia duas vezes

1. **Score ≠ afinidade experimental.** A correlação entre score do Vina e Ki/IC50 medido é modesta. Docking **classifica e prioriza** candidatos; não mede afinidade.
2. **O receptor é tratado como rígido** por padrão. Na realidade a proteína se adapta ao ligante (*induced fit*). Isso explica muitos falsos negativos.
3. **A água é ignorada** (modelo de solvente implícito), mas há casos em que uma molécula de água faz a ponte essencial entre ligante e proteína.
4. **Protonação importa muito.** Um nitrogênio protonado ou não em pH 7,4 muda completamente a eletrostática. A maior parte dos dockings ruins começa em uma preparação ruim, não no algoritmo.
5. **Entropia e dessolvatação** são tratadas de forma muito aproximada.

> Em resumo: docking é um **filtro barato e rápido** no topo de um funil. Confirmação vem de dinâmica molecular, de métodos de energia livre (MM/PBSA, FEP) e, no fim, do laboratório.

---

## 1.6 Dinâmica molecular (12 min)

### O que é

Docking entrega uma **fotografia**. A proteína real é um **filme**: ela vibra, respira, abre e fecha cavidades, em escalas de femtossegundos a milissegundos.

A **dinâmica molecular (MD)** calcula esse filme resolvendo as equações de Newton (F = ma) para cada átomo do sistema, em passos de tempo de ~2 femtossegundos (2 × 10⁻¹⁵ s), repetidamente, milhões de vezes.

### Campo de força (*force field*)

As forças não vêm da mecânica quântica (seria caro demais) e sim de um **campo de força**: um conjunto de equações clássicas simples (molas para ligações e ângulos, cossenos para torções, Lennard-Jones + Coulomb para interações não ligadas) com milhares de parâmetros ajustados a dados experimentais e quânticos.

Famílias principais: **AMBER**, **CHARMM**, **OPLS-AA**, **GROMOS**. Para moléculas pequenas (seu ligante!) usam-se extensões como GAFF (AMBER) ou CGenFF (CHARMM) — e **parametrizar um ligante novo é uma das maiores fontes de erro em MD**.

### O protocolo padrão de uma simulação

```
1. Preparar a estrutura : completar átomos faltantes, definir protonação
2. Solvatar             : colocar a proteína em uma caixa de água explícita
3. Neutralizar          : adicionar íons Na+/Cl- (concentração fisiológica ~0,15 M)
4. Minimizar energia    : eliminar colisões atômicas (descida de gradiente)
5. Equilibrar NVT       : estabilizar a temperatura (310 K ≈ 37 °C)
6. Equilibrar NPT       : estabilizar a pressão/densidade (1 bar)
7. Produção             : a simulação de verdade — de ns a µs
8. Análise              : RMSD, RMSF, Rg, ligações de H, SASA, PCA
```

### Ferramentas

| Ferramenta | Característica |
|---|---|
| **GROMACS** | Livre (LGPL), extremamente otimizado, padrão acadêmico. Linha de comando (`gmx`) |
| **OpenMM** | Biblioteca em Python, código aberto, excelente em GPU — fácil de *scriptar*, roda no Colab |
| **NAMD** | Gratuito para uso acadêmico, forte em sistemas muito grandes, integrado ao VMD |
| **AMBER** | Suíte consagrada; `AmberTools` é gratuito, o motor `pmemd` é pago |
| **MDAnalysis / MDTraj** | Bibliotecas Python para **analisar** a trajetória gerada |

### As análises que você vai reportar

| Métrica | O que mede | Como ler |
|---|---|---|
| **RMSD** vs. tempo | Quanto a estrutura se afastou da inicial | Deve **subir e estabilizar em um platô**. Platô = sistema equilibrado. Subida contínua = não convergiu, **ou** está ocorrendo uma transição conformacional real (ver a prática!) |
| **RMSF** por resíduo | Flexibilidade de cada resíduo | Picos = alças e terminais. Sítio ativo rígido é bom sinal |
| **Raio de giro (Rg)** | Compactação global | Estável = proteína mantém o enovelamento. Crescendo = está abrindo ou desenovelando |
| **Ligações de hidrogênio** ligante–proteína ao longo do tempo | Persistência do complexo | Contatos mantidos em >70% do tempo são os determinantes reais da ligação |
| **SASA** | Área exposta ao solvente | Queda = enterramento do ligante |

### O custo computacional (sem romantismo)

- Sistema típico proteína + água: 50 a 150 mil átomos.
- Em uma GPU moderna: ~100–500 ns por dia.
- Eventos biológicos relevantes (abertura de um canal, dissociação de um ligante) podem levar µs–ms → **ainda fora de alcance de rotina**.
- Por isso: **MD de 100 ns não "prova" afinidade**; ela testa a **estabilidade** da pose que o docking propôs. É assim que você deve escrever no artigo.

---

## 1.7 O pipeline integrado (5 min)

```
      sequência (UniProt)
             │
             ▼
   existe estrutura no PDB?  ──sim──►  baixar (checar resolução, R-free, nº de cadeias)
             │ não                              │
             ▼                                  │
   SWISS-MODEL / AlphaFold DB / ColabFold       │
   (checar GMQE / pLDDT / PAE)                  │
             └───────────────┬──────────────────┘
                             ▼
                  preparação do receptor
            (1 cadeia, H, cargas, remoção de água)
                             │
     ligante (PubChem) ──────┤
                             ▼
                    DOCKING (AutoDock Vina)
                             │
                   validação por REDOCKING (RMSD < 2 Å)
                             │
                             ▼
              DINÂMICA MOLECULAR (GROMACS/OpenMM)
                  a pose sobrevive 50–100 ns?
                             │
                             ▼
                 energia livre (MM/PBSA, FEP)
                             │
                             ▼
                     ENSAIO EXPERIMENTAL
```

**Mensagem central:** cada etapa é um **filtro** que reduz o número de candidatos e aumenta o custo por candidato. Bioinformática estrutural não substitui o experimento — ela decide **quais** experimentos valem o dinheiro.

---

# INTERVALO — 10 minutos

Abra o Colab e deixe um notebook novo pronto (*Arquivo → Novo notebook*).

---

# MOMENTO 2 — Aula prática (75 min)

**Todos os dados são públicos. Todo o software é livre.**

**Alvo:** acetilcolinesterase humana (AChE), PDB **4EY7** — complexada com **donepezila**, raios X, resolução **2,35 Å**, R-free **0,211**, *Homo sapiens*, **cadeias A e B** + glicanos (C, D), ligantes E20 (donepezila), NAG, EDO, NO3.
**Ligantes do exercício:** donepezila (PubChem CID 3152), galantamina (CID 9651), tacrina (CID 1935).

| Bloco | Tempo | Ferramenta |
|---|---|---|
| A — Explorar a estrutura | 12 min | RCSB PDB + Mol\* (navegador) |
| B — Modelo previsto | 15 min | AlphaFold DB + SWISS-MODEL |
| C — Docking e validação | 33 min | AutoDock Vina + Open Babel + RDKit (Colab) |
| D — Trajetória de MD | 15 min | MDAnalysis (Colab), trajetória pública pré-calculada |

> **Todas as células de código abaixo foram executadas e verificadas** antes da aula, com AutoDock Vina 1.2.5, Open Babel 3.1.0, RDKit e MDAnalysis. Os números esperados estão indicados em cada bloco. Se o seu resultado divergir muito, algo deu errado na preparação — não siga adiante.

---

## Bloco A — Explorando uma estrutura experimental (12 min)

**Só navegador. Sem código.**

1. Acesse [rcsb.org](https://www.rcsb.org) e busque por **4EY7**.
2. Na página de resumo, anote em seu caderno de laboratório:
   - Método experimental e **resolução**
   - **R-free**
   - Organismo
   - **Quantas cadeias de polímero existem** (seção *Macromolecules*) — este é o item novo e crítico
   - Quais ligantes (*Ligands*) aparecem — procure o código **E20** (donepezila)
3. Clique em **3D View / Structure** para abrir o **Mol\***.
4. Faça nesta ordem:
   - Gire a molécula (botão esquerdo) e dê zoom (*scroll*).
   - **Conte as cópias da proteína na tela.** Colora por cadeia (*Polymer Chain*) — você deve ver duas unidades da AChE, cada uma com sua própria donepezila.
   - Mude a representação para **Cartoon** — identifique hélices α e folhas β.
   - Ative **Surface** e procure a **abertura da garganta**. Use o modo de corte (*clip*) para "cortar" a proteína ao meio e enxergar o túnel.
   - Na lista de componentes, selecione **uma** das moléculas de E20, clique em *Focus* → a câmera vai até o sítio ativo.
   - Peça para mostrar as **interações** do ligante com os resíduos vizinhos.

**Perguntas de checagem (responda no chat):**

- A donepezila está na superfície ou enterrada?
- Quantos resíduos aromáticos (TRP, TYR, PHE) você vê ao redor dela?
- Localize **TRP86** (fundo da garganta) e **TRP286** (entrada, sítio periférico). Por que uma molécula que toca os dois extremos seria um inibidor potente?
- Olhando as duas cópias na tela: se você calculasse o ponto médio entre as duas donepezilas, esse ponto cairia dentro de um sítio ativo? **Guarde a resposta — ela explica o passo C4.**

---

## Bloco B — E quando não há estrutura experimental? (15 min)

### B1 — AlphaFold DB (7 min)

1. Acesse [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk) e busque **P22303** (AChE humana no UniProt).
2. Observe o modelo colorido por **pLDDT** (azul-escuro = alta confiança; laranja = baixa).
3. Identifique: o **núcleo globular** é azul; os **terminais e alças** são laranja. Por quê?
4. Abra o gráfico **PAE**. Identifique os blocos (domínios bem definidos).
5. Baixe o arquivo `.pdb` do modelo.

> Note que o AlphaFold entrega **um monômero**, sem glicanos, sem ligante, sem águas — ou seja, exatamente o que você teria de montar à mão a partir do cristal.

### B2 — Comparando previsão com realidade (8 min)

No Colab, sobrepomos o modelo do AlphaFold à estrutura cristalográfica 4EY7 e calculamos o RMSD do esqueleto. Mas há uma armadilha a resolver **antes**.

> **Não pareie os resíduos pelo número.** A AChE humana é sintetizada com um **peptídeo-sinal de 31 aminoácidos**, cortado antes de a proteína ficar pronta. O **AlphaFold** numera pela sequência do UniProt, que **inclui** o peptídeo-sinal; o **cristal** numera pela proteína **madura**, que não o tem. As duas numerações estão, portanto, **deslocadas em 31**. Parear pelo número compararia resíduos diferentes e produziria um RMSD enorme e sem sentido — mais um erro que não gera mensagem de erro.
>
> **A solução robusta:** alinhar as duas **sequências** (como na Aula 2) e parear os resíduos **pelo alinhamento**. De brinde, o código descobre e imprime o deslocamento.

```python
from Bio.PDB import PDBParser, Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.Align import PairwiseAligner

leitor = PDBParser(QUIET=True)

def carbonos_alfa(caminho, cadeia="A"):
    """Devolve (número do resíduo, letra do aminoácido, átomo CA) de cada resíduo."""
    estrutura = leitor.get_structure("x", caminho)
    return [(r.id[1], seq1(r.get_resname()), r["CA"])
            for r in estrutura[0][cadeia] if is_aa(r, standard=False) and "CA" in r]

cristal = carbonos_alfa("4EY7.pdb")          # cadeia A do cristal
modelo  = carbonos_alfa("AF-P22303.pdb")     # modelo do AlphaFold

# Alinha as SEQUÊNCIAS e pareia os resíduos pelo alinhamento
alinhador = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                            open_gap_score=-11, extend_gap_score=-1)
aln = alinhador.align("".join(x[1] for x in cristal), "".join(x[1] for x in modelo))[0]

pares = []
for (i0, i1), (j0, j1) in zip(aln.aligned[0], aln.aligned[1]):
    pares += [(cristal[i0 + k], modelo[j0 + k]) for k in range(i1 - i0)]

deslocamentos = [b[0] - a[0] for a, b in pares]
print("resíduos pareados:", len(pares))
print("deslocamento de numeração:", max(set(deslocamentos), key=deslocamentos.count))

sup = Superimposer()
sup.set_atoms([a[2] for a, _ in pares], [b[2] for _, b in pares])
print("RMSD do esqueleto (Å): %.2f" % sup.rms)
```

**Resultado esperado:** deslocamento de **+31** (o comprimento do peptídeo-sinal, conferível no UniProt) e RMSD do esqueleto (Cα) na faixa de **1 a 2 Å** — ou seja, o AlphaFold acertou o enovelamento global. Os resíduos mais divergentes serão alças e terminais, justamente os de pLDDT baixo.

> **No notebook da prática** há um passo adicional que vale a pena: o gráfico do **desvio resíduo a resíduo** contra o **pLDDT**, que mostra uma correlação **negativa** — o AlphaFold "sabe onde não sabe" — e uma tabela com a qualidade do modelo **especificamente nos resíduos do sítio ativo** (TRP86, SER203, HIS447…), que é o que de fato decide se dá para dockar nele.

**Discussão guiada:** se o enovelamento está certo, por que ainda assim é arriscado fazer docking diretamente no modelo do AlphaFold?

*Resposta a construir com a turma:* as cadeias laterais do sítio ativo podem estar em rotâmeros do estado *apo* (sem ligante), a garganta pode estar "fechada", não há a água nem o íon estrutural, e o modelo não tem o glicano. O enovelamento global está certo; os **detalhes de que o docking depende** é que são incertos. Um RMSD global de 1 Å conviver com um sítio ativo inutilizável não é contradição — é a regra.

### B3 — SWISS-MODEL (demonstração, 2 min)

O docente cola uma sequência sem estrutura no [swissmodel.expasy.org](https://swissmodel.expasy.org) e mostra ao vivo a seleção de molde, o **GMQE** e o **QMEANDisCo**. Você repetirá esse passo no exercício de casa com a proteína da sua linha de pesquisa.

---

## Bloco C — Docking com AutoDock Vina (33 min)

Tudo no Colab. **Execute as células em ordem.**

### C1 — Instalar as ferramentas (3 min)

```python
# Vina (binário oficial do GitHub), Open Babel e bibliotecas de química
!wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O vina
!chmod +x vina
!apt-get -qq install -y openbabel > /dev/null
!pip install -q rdkit
!./vina --version
!obabel -V
```

*O que cada um faz:* **Vina** é o motor de docking; **Open Babel** converte formatos, adiciona hidrogênios e calcula cargas; **RDKit** é a biblioteca de quimioinformática que gera a conformação 3D do ligante a partir do texto SMILES.

**Saída esperada:** `AutoDock Vina v1.2.5` e `Open Babel 3.x`.

### C2 — Baixar e preparar o receptor (6 min)

```python
# Baixar a estrutura do PDB
!wget -q https://files.rcsb.org/download/4EY7.pdb
!ls -lh 4EY7.pdb        # confira que tem ~1 MB; se tiver 0, a rede bloqueou (veja "Antes de começar")

# --- Receptor: SÓ a cadeia A, só a proteína (sem águas, sem glicano, sem ligante) ---
# Filtramos POR POSIÇÃO DE COLUNA: coluna 22 = identificador da cadeia.
!awk '/^ATOM/ && substr($0,22,1)=="A"' 4EY7.pdb > receptor_limpo.pdb

# --- Ligante cristalográfico: E20 da cadeia A APENAS (será a referência do redocking) ---
# colunas 18-20 = nome do resíduo ; coluna 22 = cadeia
!awk '/^HETATM/ && substr($0,18,3)=="E20" && substr($0,22,1)=="A"' 4EY7.pdb > ligante_cristal.pdb

# Conferência obrigatória
!echo "átomos do receptor:" $(grep -c . receptor_limpo.pdb)
!echo "átomos do ligante cristalográfico:" $(grep -c . ligante_cristal.pdb)
```

**Valores esperados:** receptor com alguns milhares de linhas; **ligante cristalográfico com 28 átomos** (a donepezila tem 28 átomos pesados: C₂₄H₂₉NO₃ → 24 C + 1 N + 3 O).

> **Se o ligante vier com ~56 átomos, você pegou as duas cópias** (cadeias A e B) e o centro da caixa sairá errado. Volte e confira o filtro de cadeia. Esta checagem de 5 segundos é a diferença entre um docking válido e um resultado sem sentido.

```python
# Adicionar hidrogênios em pH 7,4, calcular cargas parciais e gerar o PDBQT rígido
!obabel -ipdb receptor_limpo.pdb -opdbqt -O receptor.pdbqt -xr -p 7.4 --partialcharge gasteiger

# Conferência: o receptor precisa ter hidrogênios polares (tipo HD) e cargas não nulas
!echo "tipos de átomo:" && awk '/^ATOM/{print $NF}' receptor.pdbqt | sort | uniq -c
```

**Pare e entenda cada opção** (isto cai na prova):

| Opção | O que faz de verdade |
|---|---|
| `-xr` | Escreve o receptor como molécula **rígida** (sem árvore de torções). Obrigatório para receptor |
| `-p 7.4` | Adiciona hidrogênios conforme o estado de protonação esperado em **pH 7,4** |
| `--partialcharge gasteiger` | Calcula as **cargas parciais** de Gasteiger. **Sem esta opção o arquivo sai com todas as cargas = 0,000** e a eletrostática do docking simplesmente desaparece |
| (padrão, sem `-xh`) | O Open Babel **remove os hidrogênios apolares** e mantém os polares — exatamente o modelo de "átomos unidos" que o Vina espera |

> **Duas armadilhas de sintaxe que custam horas.**
> **(1)** No Open Babel 3, `-xh` **não** significa "adicione só hidrogênios polares" — significa *preserve hydrogens*, mantendo **todos** eles, inclusive os apolares que o Vina não quer.
> **(2)** A opção global `-d` **apaga todos os hidrogênios**. Combinada com `-p 7.4`, ela desfaz o trabalho que acabou de ser feito: o receptor sai **sem nenhum átomo HD**, e os nitrogênios doadores de ligação de hidrogênio são tipados como `NA` (aceptores). O docking roda sem erro e perde todas as ligações de hidrogênio doadas pela proteína.
> Em um teste controlado, o mesmo receptor preparado com `-xr -p 7.4` ficou com **509 hidrogênios polares (HD)**; preparado com `-xr -xh -d`, ficou com **zero**. O arquivo "funciona" nos dois casos — só um deles está certo.

Removemos as águas porque o Vina usa solvente implícito; removemos o ligante para não entupir o sítio; e removemos a cadeia B porque um receptor é **uma** proteína.

### C3 — Preparar o ligante (5 min)

```python
# Donepezila a partir do SMILES (PubChem CID 3152)
smiles = "COc1cc2c(cc1OC)C(=O)C(C2)CC1CCN(Cc2ccccc2)CC1"

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Lipinski

mol = Chem.MolFromSmiles(smiles)
print("fórmula:", rdMolDescriptors.CalcMolFormula(mol))          # esperado: C24H29NO3
print("massa molar: %.1f g/mol" % Descriptors.MolWt(mol))        # esperado: ~379,5
print("ligações rotacionáveis:", Lipinski.NumRotatableBonds(mol))# esperado: 6

mol = Chem.AddHs(mol)                        # hidrogênios explícitos, necessários para a geometria 3D
AllChem.EmbedMolecule(mol, randomSeed=42)    # gera coordenadas 3D
AllChem.MMFFOptimizeMolecule(mol)            # otimiza a geometria com o campo de força MMFF94
Chem.MolToMolFile(mol, "donepezila.sdf")

!obabel -isdf donepezila.sdf -opdbqt -O ligante.pdbqt --partialcharge gasteiger
!echo "ligações rotacionáveis marcadas no PDBQT:" $(grep -c BRANCH ligante.pdbqt)
```

**SMILES** é a notação em texto de uma molécula — você a copia direto da página do PubChem. O RDKit transforma esse texto em **coordenadas 3D**; o Open Babel adiciona cargas de Gasteiger e marca as **ligações rotacionáveis** que o Vina poderá girar.

> **Duas conferências, e nenhuma delas basta sozinha.**
>
> **1. A fórmula molecular.** Um SMILES colado errado produz uma molécula parecida que roda sem erro nenhum. `C24H29NO3`, massa 379,5 g/mol e 6 torções é a donepezila; qualquer outra coisa, pare. No notebook da prática, os SMILES são **pedidos ao PubChem pela API** e a fórmula devolvida é conferida contra a esperada — o que protege contra um CID digitado errado, que traria outra molécula em silêncio.
>
> **2. A geometria é possível?** Um SMILES pode ter a fórmula **correta** e descrever uma conectividade que **não fecha no espaço**: a fórmula bate e a molécula não existe. O teste é pedir ao RDKit que gere a conformação 3D — se `EmbedMolecule` falhar, o texto está errado, por mais convincente que pareça. *Isto não é hipotético: foi exatamente o erro encontrado ao preparar o notebook desta aula, com um SMILES de galantamina que passava na conferência da fórmula e falhava na da geometria.*

> **Na vida real, use Meeko.** O Open Babel é suficiente para a aula, mas a preparação de ligante recomendada pelos autores do Vina é a biblioteca **Meeko** (com RDKit), que trata torções, anéis macrocíclicos e estados de protonação com muito mais cuidado. É o que o script `convert_to_pdbqt.py` da disciplina usa (ver fim do Bloco C).

### C4 — Definir a caixa de busca (4 min)

```python
import numpy as np

# Lê as coordenadas x, y, z das colunas fixas do PDB (31-38, 39-46, 47-54)
coords = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                   for l in open("ligante_cristal.pdb") if l.startswith("HETATM")])

print("átomos lidos:", len(coords))        # DEVE ser 28 — se for ~56, você pegou as duas cadeias
cx, cy, cz = coords.mean(axis=0)
print("center_x = %.3f  center_y = %.3f  center_z = %.3f" % (cx, cy, cz))

# Diagnóstico: a maior dimensão do ligante indica o tamanho mínimo da caixa
print("extensão do ligante (Å):", np.round(coords.max(axis=0) - coords.min(axis=0), 1))
```

Usaremos uma caixa de **22 × 22 × 22 Å** centrada nesse ponto: grande o bastante para o ligante girar livremente (regra prática: a aresta deve ter ao menos a maior dimensão do ligante + ~8 Å), pequena o bastante para não desperdiçar busca.

### C5 — Rodar o docking (8 min)

```python
!./vina --receptor receptor.pdbqt --ligand ligante.pdbqt \
        --center_x {cx} --center_y {cy} --center_z {cz} \
        --size_x 22 --size_y 22 --size_z 22 \
        --exhaustiveness 16 --num_modes 9 --seed 42 \
        --out poses.pdbqt
```

**Tempo esperado:** cerca de **1 minuto** por ligante nesta configuração (medido em uma CPU de 2 núcleos, como a do Colab gratuito). Se levar muito mais, verifique se a caixa não ficou gigante.

Leia a tabela de saída:

```
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -12.4      0.000      0.000
   2       -11.1      1.873      2.994
   ...
```

- **mode 1** é sempre a melhor pose segundo o score.
- As colunas de RMSD medem a distância **entre as poses**, não em relação ao cristal. Um `rmsd u.b.` grande entre o modo 1 e o 2 significa que o programa encontrou duas soluções geometricamente bem diferentes — sinal de que vale olhar as duas no visualizador.
- `--seed 42` fixa o gerador aleatório → **reprodutibilidade**. Sem isso, duas execuções dão números ligeiramente diferentes.

### C6 — Validar: redocking e RMSD (5 min)

```python
# Extrair a melhor pose (modo 1) e comparar com a pose cristalográfica
!obabel -ipdbqt poses.pdbqt -opdb -O pose_melhor.pdb -f 1 -l 1
!obrms ligante_cristal.pdb pose_melhor.pdb
```

O `obrms` (parte do Open Babel) faz o **pareamento de átomos por grafo** antes de calcular o RMSD — ou seja, ele reconhece que o átomo nº 7 de um arquivo corresponde ao nº 12 do outro, o que um RMSD "ingênuo" erraria. Ele considera apenas os átomos presentes em ambos os arquivos.

**Critério:** RMSD **< 2,0 Å** → protocolo validado; podemos confiar nele para testar moléculas novas.

> **Validação do roteiro:** este mesmo pipeline (Open Babel → Vina 1.2.5 → `obrms`), aplicado a um caso-teste de referência com o ligante cristalográfico conhecido, reproduziu a pose experimental com **RMSD de 0,26 Å**. O protocolo funciona — se o seu redocking der > 2 Å, o problema está na sua preparação, muito provavelmente no centro da caixa (C4) ou nas cargas do receptor (C2).

### C7 — Triagem comparativa (2 min, se o tempo permitir)

Repita o docking para **galantamina** (CID 9651) e **tacrina** (CID 1935) e monte a tabela:

| Ligante | CID | Átomos pesados | Torções | Afinidade (kcal/mol) | *Ligand efficiency* | Comentário |
|---|---|---|---|---|---|---|
| Donepezila | 3152 | 28 | 6 | | | inibidor de alta afinidade, ocupa toda a garganta |
| Galantamina | 9651 | 21 | 1 | | | inibidor mais fraco, molécula rígida, de origem natural |
| Rivastigmina | 77991 | 18 | 4 | | | **carbamato** — a mesma família química dos inseticidas carbamatos |
| Tacrina | 1935 | 15 | 0 | | | pequena e plana, liga no fundo da garganta |

**Discussão:** a ordem do score reproduz a ordem experimental de potência? Onde ela falha e por quê?

*Pista para a discussão:* o score do Vina tem um viés conhecido a favor de moléculas **maiores**, porque mais átomos significam mais termos favoráveis de van der Waals. Comparar o score bruto de uma molécula de 28 átomos com uma de 15 é comparar coisas diferentes. É por isso que existem métricas normalizadas pelo tamanho (como a *ligand efficiency* = −ΔG/nº de átomos pesados).

> **Material de apoio da disciplina** — dois repositórios públicos, prontos para o trabalho final:
>
> - **`convert_to_pdbqt.py`** — [github.com/robertoaabfilho/to_pdbqt](https://github.com/robertoaabfilho/to_pdbqt). Converte pastas inteiras de `.sdf`, `.mol2`, `.pdb` e `.cif` para PDBQT, detectando automaticamente se o arquivo é ligante ou receptor. Usa **Meeko + RDKit** para ligantes (a preparação recomendada) e Open Babel/gemmi para receptores.
> - **`vina_docking.py`** — [github.com/robertoaabfilho/vina_docking](https://github.com/robertoaabfilho/vina_docking). Roda o Vina em lote no produto cartesiano **todos os ligantes × todas as proteínas**, com um único `config.txt` (formato INI, uma seção por proteína com centro, tamanho da caixa e parâmetros), retomada automática e saída consolidada em CSV.
>
> Use-os quando precisar testar dezenas ou centenas de moléculas — é exatamente o cenário do Exercício 4.

---

## Bloco D — Dinâmica molecular: analisando uma trajetória (15 min)

Uma simulação de verdade não cabe em 15 minutos de aula: o complexo AChE–donepezila solvatado tem ~100 mil átomos e 100 ns levariam horas de GPU. Então vamos **analisar uma trajetória pronta** — e aprender a ler os gráficos, que é a habilidade transferível.

Usaremos a trajetória de referência que acompanha a própria biblioteca **MDAnalysis**: a **adenilato quinase (AdK) de *E. coli*** (estrutura de partida 4AKE), ~1 ns, 98 quadros, 3.341 átomos. É um conjunto de dados público, citável e que já vem instalado — ninguém precisa baixar nada.

> **Por que a AdK e não a AChE?** Porque a AdK faz, em 1 ns, algo que a AChE não faria: ela **abre os seus domínios de forma visível**. É um dos exemplos didáticos mais usados do mundo justamente porque mostra, em poucos quadros, a diferença entre "o sistema não convergiu" e "está acontecendo uma transição conformacional real" — a interpretação mais difícil e mais cobrada em MD.

```python
!pip install -q MDAnalysis MDAnalysisTests

import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from MDAnalysis.tests.datafiles import PSF, DCD
import numpy as np, matplotlib.pyplot as plt

u = mda.Universe(PSF, DCD)
print("átomos:", len(u.atoms), "| resíduos:", len(u.residues), "| quadros:", len(u.trajectory))
print("tempo total: %.0f ps" % u.trajectory.totaltime)
```

**Saída esperada:** 3341 átomos · 214 resíduos · 98 quadros · 97 ps.

### D1 — RMSD × tempo

```python
R = rms.RMSD(u, u, select="backbone", ref_frame=0).run()
t, rmsd = R.results.rmsd[:, 1], R.results.rmsd[:, 2]

plt.figure(figsize=(7, 4))
plt.plot(t, rmsd, lw=1.8)
plt.xlabel("quadro"); plt.ylabel("RMSD do esqueleto (Å)")
plt.title("RMSD em relação ao quadro inicial"); plt.grid(alpha=.3); plt.show()
print("RMSD final: %.2f Å  (máximo %.2f Å)" % (rmsd[-1], rmsd.max()))
```

**Resultado esperado:** o RMSD **sobe continuamente** até ~**6,8 Å** e **não** atinge platô.

**Interprete — e resista à conclusão fácil.** Em um complexo proteína–ligante equilibrado, isso significaria "não convergiu, simulação inútil". Aqui, significa outra coisa: a proteína está executando a **transição do estado fechado para o aberto**, que é a sua função biológica. Como distinguir os dois casos? Com a métrica seguinte.

### D2 — Raio de giro: a proteína está abrindo ou desenovelando?

```python
prot = u.select_atoms("protein")
rg = np.array([prot.radius_of_gyration() for _ in u.trajectory])

plt.figure(figsize=(7, 4))
plt.plot(rg, lw=1.8, color="darkorange")
plt.xlabel("quadro"); plt.ylabel("raio de giro (Å)")
plt.title("Compactação global"); plt.grid(alpha=.3); plt.show()
print("Rg: %.2f Å  ->  %.2f Å" % (rg[0], rg[-1]))
```

**Resultado esperado:** o Rg cresce de ~**16,7 Å** para ~**19,6 Å**.

**Interprete:** a proteína ficou **mais aberta** (+3 Å de raio de giro) mas manteve a estrutura secundária — é abertura de domínios, não desenovelamento. Se fosse desenovelamento, o Rg dispararia e a estrutura secundária colapsaria. **A leitura correta de uma trajetória nunca vem de um gráfico só.**

### D3 — RMSF por resíduo: quem se move?

```python
align.AlignTraj(u, u, select="protein and name CA", in_memory=True).run()
ca = u.select_atoms("protein and name CA")
F = rms.RMSF(ca).run()

plt.figure(figsize=(9, 4))
plt.plot(ca.resids, F.results.rmsf, lw=1.5, color="seagreen")
plt.xlabel("número do resíduo"); plt.ylabel("RMSF (Å)")
plt.title("Flexibilidade por resíduo"); plt.grid(alpha=.3); plt.show()
print("RMSF média %.2f Å | máxima %.2f Å no resíduo %d"
      % (F.results.rmsf.mean(), F.results.rmsf.max(), ca.resids[F.results.rmsf.argmax()]))
```

**Resultado esperado:** RMSF média ~**1,9 Å**, com máximo de ~**5,7 Å** em torno do **resíduo 149**.

**Interprete:** os picos não estão distribuídos ao acaso — concentram-se nos domínios móveis (LID e NMP, as "tampas" da AdK), enquanto o núcleo permanece rígido. **É exatamente esse tipo de gráfico que você usaria para argumentar que o seu sítio ativo é estável enquanto a periferia respira.**

> **Nota de honestidade metodológica:** o passo `align.AlignTraj` é obrigatório antes do RMSF. Sem alinhar, você mede o movimento da proteína **como um corpo rígido pela caixa de simulação** (translação e rotação globais), não a flexibilidade interna. É um erro silencioso muito comum: os números saem grandes e plausíveis, e não querem dizer nada.

### Para rodar a sua própria simulação depois

A sequência mínima no GROMACS (referência, não precisa executar hoje):

```bash
gmx pdb2gmx  -f proteina.pdb -o proc.gro -water tip3p      # topologia + hidrogênios
gmx editconf -f proc.gro -o box.gro -c -d 1.0 -bt cubic    # caixa
gmx solvate  -cp box.gro -cs spc216.gro -o solv.gro -p topol.top
gmx grompp   -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
gmx genion   -s ions.tpr -o solv_ions.gro -p topol.top -neutral -conc 0.15
gmx grompp -f minim.mdp ... && gmx mdrun -v -deffnm em     # minimização
gmx grompp -f nvt.mdp   ... && gmx mdrun -deffnm nvt       # equilíbrio NVT
gmx grompp -f npt.mdp   ... && gmx mdrun -deffnm npt       # equilíbrio NPT
gmx grompp -f md.mdp    ... && gmx mdrun -deffnm md        # produção
```

Para um complexo com ligante, o passo extra (e mais delicado) é **parametrizar o ligante** — tipicamente com CGenFF (CHARMM) ou ACPYPE/GAFF (AMBER). É aí que mora a maior parte dos erros em MD de complexos.

---

# MOMENTO 3 — Integração e fechamento (20 min)

## 3.1 Costurando teoria e prática (8 min)

Revisitamos cada conceito teórico pelo que vimos na tela:

| Na teoria falamos de… | Na prática você viu… |
|---|---|
| Resolução e R-free | 4EY7 a 2,35 Å, R-free 0,211 — boa o bastante para docking |
| Unidade assimétrica e cadeias | A 4EY7 tem **duas** AChE; filtrar a cadeia A foi o que salvou o centro da caixa |
| Colunas de largura fixa do PDB | Por que filtramos por `substr($0,22,1)` e não por "procurar o texto" |
| B-factor | No cristal é mobilidade atômica; no modelo do AlphaFold o mesmo campo guarda o **pLDDT** — e a interpretação se inverte |
| A garganta aromática da AChE | A cavidade revelada pela representação de **superfície** com corte |
| Modelagem ≠ estrutura experimental | Enovelamento correto (RMSD ~1–2 Å), mas cadeias laterais do sítio incertas |
| Protonação e cargas | `-p 7.4` e `--partialcharge gasteiger`; sem eles, receptor sem doadores de H e com carga zero |
| Ligações rotacionáveis | Os 6 `BRANCH` que o Open Babel marcou no PDBQT e que o Vina girou durante a busca |
| Função de pontuação | A coluna "affinity" da tabela de saída — **uma estimativa, não uma medida** — e seu viés a favor de moléculas grandes |
| Validação | O redocking com RMSD < 2 Å; sem ele, o número não tem valor |
| Receptor rígido | A razão pela qual passamos para MD: deixar a proteína responder ao ligante |
| Convergência | O RMSD que **não** fez platô — e o Rg que mostrou que era abertura de domínio, não falha |
| Referencial de análise | O `AlignTraj` antes do RMSF: sem ele mede-se o movimento da caixa, não da proteína |

## 3.2 Os oito erros mais comuns (5 min)

1. Esquecer de **remover as águas** (ou removê-las quando uma delas é estrutural e essencial).
2. **Não filtrar a cadeia** — manter duas cópias da proteína no receptor, ou calcular o centro da caixa a partir de dois ligantes.
3. Não **adicionar hidrogênios** / ignorar o estado de protonação em pH 7,4 — ou anulá-los com uma opção contraditória.
4. Gerar o receptor **sem cargas parciais** (todas zero) e não perceber.
5. Caixa de busca **mal posicionada ou pequena demais** — o ligante não cabe e o score "explode".
6. Usar uma estrutura de **baixa resolução** ou do organismo errado (modelo de camundongo para discutir um fármaco humano).
7. **Não validar** o protocolo com redocking — ou validar e depois tratar o redocking como se ele tivesse validado o *score*.
8. Reportar o score do Vina **como se fosse afinidade medida** (ΔG experimental), e **não registrar** versões, sementes e parâmetros. *Docking é estocástico: sem `--seed` e sem `exhaustiveness` anotados, ninguém — nem você — repete seu resultado.*

## 3.3 O que levar para o artigo final (5 min)

O trabalho final exige **integrar pelo menos duas técnicas** da disciplina. Combinações naturais com esta aula:

- **Estrutural + Sequências (Aula 2):** buscar homólogos da sua proteína-alvo com BLAST, alinhar, identificar resíduos conservados e mostrar que eles formam o sítio ativo na estrutura 3D.
- **Estrutural + Transcriptômica (Aula 4):** um gene superexpresso nos seus dados de RNA-seq vira o alvo; você prevê a estrutura da proteína e faz docking de inibidores candidatos.
- **Estrutural + Proteômica (Aula 5):** uma modificação pós-traducional identificada por espectrometria de massas é mapeada na estrutura 3D — ela fica no sítio ativo, numa interface, ou numa alça exposta?
- **Estrutural + Epigenômica (Aula 6):** as enzimas escritoras e apagadoras de marcas epigenéticas (DNMT1, HDACs, EZH2) são alvos de fármacos com estrutura no PDB — feche o ciclo da aula passada fazendo docking em uma delas.
- **Estrutural + Filogenia (Aula 8):** comparar o sítio ativo entre espécies/variantes e discutir seletividade (por exemplo, por que um inseticida atinge o inseto e menos o mamífero).
- **Estrutural + Vias/Redes (Aula 9):** situar o alvo em uma via do KEGG ou em uma rede do STRING e justificar a escolha.

**Checklist de metodologia para o seu artigo:**

- [ ] Código PDB (ou fonte do modelo: AlphaFold DB / SWISS-MODEL) + resolução + R-free (ou pLDDT/GMQE)
- [ ] **Qual cadeia foi usada** e o que foi removido (águas, íons, glicanos, outras cadeias)
- [ ] Como o receptor foi preparado (hidrogênios, pH, método de cargas, software **e versão**)
- [ ] Origem e preparação do ligante (CID do PubChem, SMILES, gerador de conformação 3D)
- [ ] Centro e dimensões da caixa
- [ ] `exhaustiveness`, `num_modes`, `seed`, número de réplicas
- [ ] **Validação por redocking com RMSD reportado**
- [ ] Se houve MD: campo de força, modelo de água, íons, tempo de simulação, critério de convergência, e **se a trajetória foi alinhada** antes das análises
- [ ] Figuras com escala e legenda (nunca uma figura de PyMOL/Mol\* sem indicar os resíduos)

## 3.4 Exercícios (entrega até a Aula 8) — 2 min

1. **(Obrigatório)** Escolha uma proteína da sua linha de pesquisa. Verifique se há estrutura no PDB. Se houver, registre código, método, resolução, R-free **e número de cadeias**; se não houver, gere um modelo no SWISS-MODEL **e** busque no AlphaFold DB, e compare GMQE/QMEANDisCo com o pLDDT. Entregue meia página com a justificativa da estrutura que você usaria.
2. **(Obrigatório)** Refaça o docking do Bloco C com `exhaustiveness` 8, 16 e 32, três vezes cada (varie a `--seed`). Monte uma tabela com média e desvio-padrão da melhor afinidade. **Pergunta:** a partir de qual valor o resultado estabiliza? O que isso diz sobre reportar um único número com duas casas decimais?
3. **(Obrigatório)** Faça o docking dos demais ligantes da tabela do Bloco C7 e compare a ordem dos scores com os valores de IC50 publicados (busque no ChEMBL). A ordem bate? Escreva um parágrafo discutindo as causas de divergência — inclua o efeito do tamanho do ligante e calcule a *ligand efficiency* de cada um.
4. **(Obrigatório, diagnóstico)** Refaça **de propósito** o passo C2 sem filtrar a cadeia (`awk '/^HETATM/ && substr($0,18,3)=="E20"'`), recalcule o centro e rode o docking. Relate: quantos átomos o ligante "cristalográfico" passou a ter, para onde foi o centro da caixa, qual a afinidade reportada e qual o RMSD do redocking. **Este exercício existe para você ver com os próprios olhos um resultado errado que não dá erro.**
5. **(Desafio)** Use o `vina_docking.py` do repositório da disciplina para dockar um conjunto de 20 ligantes do PubChem contra a 4EY7 e contra a butirilcolinesterase (PDB **1P0I**). Discuta a seletividade entre as duas enzimas — o que isso significa para um *bioscavenger*?
6. **(Desafio)** No notebook do Bloco D, calcule o número de ligações de hidrogênio intramoleculares ao longo do tempo e verifique se a estrutura secundária se mantém durante a abertura dos domínios. Isso confirma ou refuta a interpretação de "abertura, não desenovelamento"?

---

## Glossário rápido

| Termo | Significado |
|---|---|
| **Å (ångström)** | 10⁻¹⁰ m. Uma ligação C–C mede ~1,5 Å |
| **Apo / Holo** | Estrutura sem ligante / com ligante |
| **B-factor** | Mobilidade/incerteza de cada átomo (nos modelos AlphaFold guarda o pLDDT, com interpretação invertida) |
| **Cadeia (*chain*)** | Identificador de cada polímero no arquivo PDB (A, B, C…). Um arquivo pode ter várias cópias da mesma proteína |
| **Campo de força** | Conjunto de equações e parâmetros clássicos que define as forças entre átomos na MD |
| **Docking** | Previsão da pose e do score de um ligante em um receptor |
| **Exhaustiveness** | Esforço de busca do Vina (nº de buscas independentes) |
| **GMQE / QMEANDisCo** | Métricas de qualidade de modelos do SWISS-MODEL |
| **Grid box** | Caixa onde o docking procura |
| **Induced fit** | Adaptação conformacional da proteína ao ligante |
| **Ligand efficiency** | −ΔG dividido pelo nº de átomos pesados; corrige o viés de tamanho do score |
| **MSA** | Alinhamento múltiplo de sequências |
| **PAE** | Erro esperado na posição relativa entre dois resíduos (AlphaFold) |
| **PDBQT** | Formato do AutoDock: PDB + cargas + tipos de átomo + torções |
| **pLDDT** | Confiança local do AlphaFold (0–100) |
| **Pose** | Posição + orientação + conformação do ligante no sítio |
| **RMSD** | Desvio quadrático médio entre duas estruturas (Å) |
| **RMSF** | Flutuação de cada resíduo ao longo de uma trajetória (exige trajetória alinhada) |
| **Rotâmero** | Uma das conformações possíveis de uma cadeia lateral |
| **SMILES** | Notação textual de uma molécula |
| **Scoring function** | Equação que estima a energia de uma pose |
| **Solvente implícito/explícito** | Água representada por um modelo contínuo / por moléculas reais |
| **Unidade assimétrica** | O conteúdo depositado no arquivo PDB — não necessariamente a unidade biológica |

---

## Recursos e links

**Bases de dados**

- RCSB PDB — https://www.rcsb.org
- PDBe (Europa) — https://www.ebi.ac.uk/pdbe
- AlphaFold DB — https://alphafold.ebi.ac.uk
- PubChem — https://pubchem.ncbi.nlm.nih.gov
- ChEMBL — https://www.ebi.ac.uk/chembl
- UniProt — https://www.uniprot.org

**Ferramentas (todas livres ou gratuitas para uso acadêmico)**

- Mol\* — https://molstar.org
- UCSF ChimeraX — https://www.cgl.ucsf.edu/chimerax
- PyMOL (open-source) — https://github.com/schrodinger/pymol-open-source
- SWISS-MODEL — https://swissmodel.expasy.org
- ColabFold — https://github.com/sokrypton/ColabFold
- AutoDock Vina — https://vina.scripps.edu · releases: https://github.com/ccsb-scripps/AutoDock-Vina/releases
- Smina / GNINA — https://github.com/gnina/gnina
- Open Babel — https://openbabel.org
- RDKit — https://www.rdkit.org · Meeko — https://github.com/forlilab/Meeko
- GROMACS — https://www.gromacs.org · OpenMM — https://openmm.org
- MDAnalysis — https://www.mdanalysis.org

**Scripts da disciplina**

- `convert_to_pdbqt.py` — https://github.com/robertoaabfilho/to_pdbqt
- `vina_docking.py` — https://github.com/robertoaabfilho/vina_docking

**Leituras recomendadas**

- Berman, H. M. et al. *The Protein Data Bank*. Nucleic Acids Research, 2000.
- Trott, O.; Olson, A. J. *AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading*. J. Comput. Chem., 2010.
- Eberhardt, J. et al. *AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings*. J. Chem. Inf. Model., 2021.
- Jumper, J. et al. *Highly accurate protein structure prediction with AlphaFold*. Nature, 2021.
- Waterhouse, A. et al. *SWISS-MODEL: homology modelling of protein structures and complexes*. Nucleic Acids Res., 2018.
- Abraham, M. J. et al. *GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers*. SoftwareX, 2015.
- Michaud-Agrawal, N. et al. *MDAnalysis: A toolkit for the analysis of molecular dynamics simulations*. J. Comput. Chem., 2011.
- Cheung, J. et al. *Structures of human acetylcholinesterase in complex with pharmacologically important ligands*. J. Med. Chem., 2012. (origem da estrutura 4EY7)
- O'Boyle, N. M. et al. *Open Babel: An open chemical toolbox*. J. Cheminformatics, 2011.

---

*Material da disciplina Introdução à Bioinformática — PPGED/IME. Dados e ferramentas de uso público e livre. Todas as células de código deste roteiro foram executadas e verificadas antes da aula.*
