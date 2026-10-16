# Aula 4: Transcriptômica

### RNA-seq e análise de expressão diferencial

**Introdução à Bioinformática – PPGED/IME**

---

## Como usar este material

Este roteiro acompanha a Aula 4. Ele tem duas partes:

- **Parte 1 – Conceitos:** o que é o transcriptoma, como ele é medido e quais ferramentas são usadas em cada etapa. Leia antes ou depois da aula, no seu ritmo.
- **Parte 2 – Guia da prática:** um passo a passo do notebook [Prática_de_Transcriptoma.ipynb](https://drive.google.com/file/d/1DSFiqZ7GWNl56d0P9gRzePEwxhqa0ePy/view?usp=sharing), explicando o que cada etapa faz e como interpretar o que aparece na tela.

Você **não precisa saber programar**. Na prática, o código já está pronto: seu trabalho é executar cada célula e entender o resultado. Ao final há um glossário e a lista de exercícios.

### O que você vai aprender

Ao final desta aula, você deve conseguir:

1. Explicar a diferença entre genoma e transcriptoma e o que o RNA-seq mede.
2. Descrever as etapas de um experimento de RNA-seq, da amostra até a lista de genes.
3. Dizer o que fazem as ferramentas HISAT2, STAR e Trinity, e quando usar cada uma.
4. Explicar por que comparar genes exige normalização e estatística.
5. Executar uma análise com o DESeq2 e interpretar seus resultados e gráficos.

---

# PARTE 1: CONCEITOS

## 1. Por que estudar o transcriptoma?

Comece com uma pergunta: todas as células do seu corpo têm o mesmo DNA. Uma célula do cérebro e uma do fígado carregam exatamente as mesmas instruções. Então por que elas são tão diferentes?

A resposta não está no DNA que elas **têm**, mas nas partes do DNA que elas estão **usando**.

### A analogia da biblioteca

Imagine uma grande fábrica com uma biblioteca de manuais: como montar cada peça, como consertar cada máquina, como reagir a um incêndio.

- A biblioteca completa é o **genoma**: todas as instruções do organismo, guardadas no DNA.
- Ninguém leva o livro original para a linha de produção. O que vai para o chão de fábrica são **fotocópias** dos manuais em uso naquele dia. O conjunto dessas cópias é o **transcriptoma**.

O neurônio e a célula do fígado têm a mesma biblioteca, mas fazem cópias de manuais diferentes. Mesma biblioteca, cópias diferentes, funções diferentes.

### O dogma central

Trocando a analogia pelos nomes reais:

```
DNA  ──(transcrição)──►  RNA mensageiro (mRNA)  ──(tradução)──►  Proteína
manual original           fotocópia                               peça fabricada
```

- Um trecho de DNA com uma instrução se chama **gene**.
- A célula copia o gene em uma molécula de **mRNA**. Esse processo se chama **transcrição** (daí o nome *transcriptoma*).
- A célula lê o mRNA e fabrica uma **proteína**, que executa o trabalho. Esse processo se chama **tradução**.

Medir o mRNA mostra quais instruções estão em uso, ou seja, **o que a célula está fazendo agora**.

### O transcriptoma muda o tempo todo

O genoma praticamente não muda ao longo da vida. O transcriptoma muda a cada instante.

Pense em cada gene como uma lâmpada com **regulador de intensidade**: pode estar apagada, fraca ou acesa no máximo. Quando a célula é estressada, infectada ou exposta a uma substância química, ela ajusta essas lâmpadas. Medir o transcriptoma é fotografar esse painel. Comparar duas fotografias, antes e depois de um estímulo, mostra como a célula reagiu.

### Por que isso interessa à Engenharia de Defesa?

| Problema | O que o transcriptoma revela |
|---|---|
| **Agentes químicos e toxinas** | Quais genes de proteção a célula ativa e com que intensidade. Ajuda a entender o mecanismo de ação e a buscar contramedidas. |
| **Patógenos** | Quais genes uma bactéria liga ao infectar um hospedeiro. Aponta alvos para medicamentos, vacinas e testes de diagnóstico. |
| **Radiação** | Genes de reparo do DNA que mudam com a exposição. Podem servir como **biomarcadores** para triagem após um incidente. |

A lógica é sempre a mesma: comparar quem foi exposto com quem não foi e descobrir o que mudou.

---

## 2. Do tubo de ensaio aos arquivos

Guarde este diagrama. Ele resume todo o experimento e vai aparecer ao longo da aula:

```
Amostra biológica → extração de RNA → conversão em cDNA → sequenciador
      → arquivos FASTQ (milhões de "leituras")
      → controle de qualidade (FastQC)
      → mapeamento OU montagem
      → tabela de contagens (gene × amostra)
      → análise estatística → genes diferencialmente expressos
```

### Etapa de laboratório: do RNA ao cDNA

Coleta-se a amostra (tecido, cultura de células ou de bactérias), rompem-se as células e separa-se o RNA. Há um problema: o RNA é **frágil**, e os sequenciadores mais usados foram feitos para ler **DNA**.

A solução é converter o RNA em DNA com uma enzima chamada **transcriptase reversa**, que faz o caminho inverso da transcrição. A cópia resultante se chama **cDNA** (DNA complementar).

> **Analogia:** é como se as fotocópias fossem feitas em papel térmico, que desbota em poucos dias. Para não perder a informação, digitalizamos cada uma. O conteúdo é o mesmo, mas o formato agora é estável.

O material também é cortado em fragmentos pequenos, do tamanho que a máquina consegue ler.

### O sequenciador e as leituras

O **sequenciador** lê a ordem das letras (A, C, G, T) de cada fragmento. Cada fragmento lido gera uma **leitura** (*read*), um trecho curto de 50 a 150 letras. Um gene tem milhares de letras, então cada leitura é só um pedacinho dele.

> **Analogia:** passe todas as fotocópias da fábrica numa picotadora de papel. Cada tira é uma leitura.

Uma única amostra gera **dezenas de milhões** de leituras. Por isso a análise exige computação.

### O arquivo FASTQ

As leituras chegam num arquivo de texto chamado **FASTQ**, já visto na Aula 2. Cada leitura ocupa 4 linhas:

```
@leitura_001                ← nome (etiqueta de identificação)
ACGTTGCAAGTCCGATTGCA        ← a sequência
+                           ← separador
IIIIIHHHGGFFF###IIII        ← qualidade de cada letra
```

A última linha é uma "nota de confiança" letra por letra. Qualidade 30, por exemplo, significa uma chance de erro em mil. Os símbolos `#` indicam letras de baixa confiança. O **FastQC** analisa essa informação antes de seguir adiante.

### O princípio da medição

Esta é a ideia central de todo o RNA-seq:

> Quanto mais ativo um gene, mais cópias de mRNA ele tem, e mais leituras dele aparecem.
> **Contar leituras é medir a atividade dos genes (expressão gênica).**

### Réplicas biológicas

Você não aprovaria um material com base num único ensaio de tração: cada corpo de prova dá um valor um pouco diferente. Com células é igual. Dois animais, duas culturas ou dois pacientes nunca têm exatamente a mesma atividade gênica.

Sem várias amostras por condição, não há como saber se uma diferença veio do tratamento ou do acaso. A recomendação é ter **pelo menos 3 réplicas biológicas por condição**.

> ⚠️ Réplica biológica é **outra** amostra (outro organismo, outra cultura). Sequenciar o mesmo tubo duas vezes não conta.

---

## 3. Dois caminhos: com ou sem genoma de referência

Uma leitura não traz etiqueta dizendo de que gene veio. Descobrir isso é como montar um quebra-cabeça gigante, e há duas situações:

- **Com a foto da caixa:** o organismo já tem o genoma sequenciado (humano, camundongo, *E. coli*). Esse genoma conhecido é o **genoma de referência**. Basta achar o lugar de cada peça na foto: é o **mapeamento**.
- **Sem a foto:** o organismo é pouco estudado (um microrganismo do solo, uma planta nativa da Amazônia). É preciso encaixar as peças umas nas outras: é a **montagem *de novo*** (do zero).

### Por que o RNA é um desafio: éxons e íntrons

Na maioria dos genes de plantas e animais, o gene no DNA não é um bloco contínuo. Ele tem trechos que vão para o mRNA, os **éxons**, intercalados com trechos removidos, os **íntrons**. Depois da transcrição, a célula corta os íntrons e cola os éxons. Esse processo se chama **splicing**.

```
Gene no DNA:   [éxon 1]───íntron───[éxon 2]───íntron───[éxon 3]
mRNA:          [éxon 1][éxon 2][éxon 3]
```

Uma leitura pode cair justamente na emenda entre dois éxons. Ao procurá-la no genoma, ela aparece "partida": metade de um lado do íntron, metade do outro. As ferramentas da Aula 2 (BWA e Bowtie2) não foram feitas para isso. O RNA exige alinhadores que **saltam íntrons**.

### Caminho A: mapeamento

**HISAT2**
- Programa de linha de comando que descobre **de que posição do genoma veio cada leitura**, inclusive quando ela atravessa um íntron.
- Usa um **índice** do genoma, como o índice remissivo de um livro técnico: em vez de folhear o livro inteiro, vai direto à página certa.
- Ponto forte: usa **pouca memória** (cerca de 8 GB para o genoma humano). Roda num computador comum.

**STAR**
- Faz a mesma tarefa que o HISAT2.
- É extremamente rápido e preciso, por isso é o padrão em grandes consórcios de pesquisa.
- Exige muita memória: cerca de **32 GB** para o genoma humano.

> **Regra prática:** computador modesto → HISAT2; servidor ou cluster → STAR.

**A saída** dos dois é um arquivo **SAM** (ou sua versão compactada, **BAM**), que diz: "a leitura X caiu na posição Y do cromossomo Z". Esses arquivos são manipulados com o **SAMtools** (Aula 2) e podem ser vistos no **IGV**, onde as leituras aparecem empilhadas sobre os éxons.

**Ensembl**
Portal público e gratuito do Instituto Europeu de Bioinformática, com genomas de centenas de organismos. Para RNA-seq, baixamos dele:
- o **genoma** (formato FASTA): "a foto da caixa", a sequência de cada cromossomo;
- a **anotação** (formato GTF): "o mapa dos genes", onde cada gene e cada éxon começa e termina.

Cada gene tem um código estável no Ensembl. Por exemplo, o gene humano *FKBP5* é o **ENSG00000096060**. Quando um gene aparece nos seus resultados, é no Ensembl que você descobre o que ele faz.

**Contagem**
Programas auxiliares (como *featureCounts* ou *HTSeq*) cruzam as posições das leituras com o mapa dos genes e contam quantas leituras caíram em cada gene. O resultado é a **tabela de contagens**: uma linha por gene, uma coluna por amostra. É o produto central do RNA-seq.

### Caminho B: montagem *de novo*

**Trinity**
- Reconstrói os transcritos a partir das próprias leituras, sem referência.
- Lembre da picotadora: como muitas cópias do mesmo manual foram picotadas, os trechos se sobrepõem. Encaixando as sobreposições, dá para reconstruir o texto sem nunca ter visto o original.
- Funciona em três etapas, com nomes inspirados na metamorfose de uma borboleta:
  - **Inchworm** (lagarta): junta leituras em sequências mais longas;
  - **Chrysalis** (casulo): agrupa as sequências que parecem vir do mesmo gene;
  - **Butterfly** (borboleta): resolve as variações e entrega os transcritos finais.
- A saída é um arquivo FASTA com transcritos ainda sem nome. Para identificá-los, compara-se cada um com bancos de dados usando o **BLAST** (Aula 2).
- É pesado: exige muita memória e pode levar horas ou dias. Roda em servidores ou na nuvem.

### Como escolher

| Situação | Caminho | Ferramenta |
|---|---|---|
| Organismo com genoma de referência, computador modesto | Mapeamento | HISAT2 |
| Organismo com genoma de referência, servidor ou cluster | Mapeamento | STAR |
| Organismo sem genoma de referência | Montagem *de novo* | Trinity + BLAST |

A primeira pergunta é sempre: **o meu organismo tem genoma de referência?** Confira no Ensembl ou no NCBI antes de começar.

---

## 4. Da contagem à descoberta: estatística

Com a tabela de contagens em mãos, queremos saber **quais genes mudaram entre as condições**. Isso se chama **expressão diferencial**. Mas não basta comparar os números diretamente, por três motivos.

**1. Profundidade de sequenciamento diferente.**
Se a amostra A gerou 40 milhões de leituras e a B gerou 20 milhões, todos os genes de A parecem "duas vezes mais expressos". É preciso **normalizar**.
> Analogia: comparar o consumo de combustível de dois veículos sem dividir pela distância percorrida.

**2. Variação biológica.**
Réplicas da mesma condição nunca dão o mesmo número. É preciso saber se a diferença observada é maior que esse "ruído" natural.

**3. Milhares de testes ao mesmo tempo.**
Testar 20 mil genes com o critério p < 0,05 gera cerca de mil "falsos positivos" só por acaso. Por isso se usa o **p-valor ajustado (padj)**, corrigido para comparações múltiplas.

### As medidas que você precisa saber ler

| Medida | O que significa | Como ler |
|---|---|---|
| **baseMean** | Expressão média do gene | Genes muito pouco expressos dão resultados frágeis |
| **log2FoldChange** | Quanto o gene mudou, em escala log de base 2 | +1 = dobrou; −1 = caiu pela metade; +3 = aumentou 8× |
| **lfcSE** | Incerteza (erro padrão) do log2FoldChange | Valor alto = estimativa pouco confiável |
| **padj** | Probabilidade de a diferença ser acaso, já corrigida | Abaixo de 0,05 = diferença confiável |

> **Para converter log2FoldChange em "vezes":** calcule 2 elevado ao valor. Exemplo: 2^4,57 ≈ 24 vezes.

### As ferramentas

**DESeq2**
Pacote da linguagem R, distribuído pelo **Bioconductor** (o repositório de pacotes de R para biologia). Recebe a tabela de contagens e a descrição das amostras, normaliza, estima o ruído de cada gene e testa cada um. É robusto com poucas réplicas e muito usado na literatura. É a ferramenta da nossa prática.

**edgeR**
Também é um pacote R/Bioconductor, com o mesmo modelo estatístico de base, mas com outra forma de normalizar (chamada TMM). Os dois são confiáveis e dão resultados parecidos. O importante é escolher um e informá-lo corretamente no seu artigo.

### Os gráficos

- **PCA:** cada ponto é uma **amostra**. Pontos próximos têm perfis parecidos. Serve para verificar se o experimento é consistente.
- **MA plot:** cada ponto é um **gene**. Mostra a expressão média (eixo X) contra a mudança (eixo Y).
- **Volcano plot:** cada ponto é um **gene**. Cruza o tamanho da mudança (eixo X) com a certeza estatística (eixo Y). Os genes interessantes ficam nos dois "braços" do vulcão.

---

# PARTE 2: GUIA DA PRÁTICA

## 5. O experimento que vamos analisar

Vamos usar o conjunto de dados **airway**, de um estudo real publicado por Himes e colaboradores na revista *PLoS One* em 2014 (dados públicos no GEO, código GSE52778).

- **As células:** músculo liso das vias aéreas, o músculo que envolve os brônquios. Na asma, ele se contrai e inflama, estreitando a passagem do ar.
- **O tratamento:** **dexametasona**, um corticoide sintético potente, usado no tratamento da asma. As células receberam 1 micromolar durante 18 horas.
- **O desenho:** células de **4 doadores**. As células de cada doador foram divididas em duas partes: uma controle e uma tratada.
- **A pergunta:** quais genes o corticoide liga e desliga nessas células?

**Como os dados foram gerados:** as leituras foram baixadas do SRA, mapeadas com o **STAR** no genoma humano usando as anotações do **Ensembl** e contadas por gene. É exatamente o fluxo da Parte 1. Nós começamos a partir da **tabela de contagens** pronta.

### Antes de começar

1. Abra o notebook da prática: **[Prática_de_Transcriptoma.ipynb](https://drive.google.com/file/d/1DSFiqZ7GWNl56d0P9gRzePEwxhqa0ePy/view?usp=sharing)**.
   - Se o Google Drive mostrar só uma pré-visualização, clique em **Abrir com → Google Colaboratory**. O Colab funciona no navegador, sem instalar nada.
   - Em seguida, clique em **Arquivo → Salvar uma cópia no Drive**. Assim você trabalha na sua própria cópia, sem alterar o arquivo original.
2. Mude a linguagem para R: *Ambiente de execução → Alterar tipo de ambiente de execução → R*.
3. Execute as células **em ordem, de cima para baixo** (Shift + Enter).
4. Se aparecer um erro de "objeto não encontrado", provavelmente uma célula anterior foi pulada. Use *Ambiente de execução → Executar anteriores*.

---

## 6. Passo a passo do notebook

### Seção 1: Preparação do ambiente

Instala os pacotes da aula: **airway** (os dados), **DESeq2** (a análise) e **ggplot2** (gráficos).

⏱️ Leva de 5 a 10 minutos. Mensagens em vermelho durante a instalação são normais.

### Seção 2: Conhecendo os dados

O objeto `airway` guarda três peças juntas, como uma planilha central com duas fichas presas nas bordas:

```
                    colData (ficha das amostras)
                    ┌──────────────────────────────┐
                    │ doador, tratamento, códigos  │
                    └──────────────────────────────┘
┌───────────────┐   ┌──────────────────────────────┐
│ rowData       │   │ assay: tabela de contagens   │
│ (ficha dos    │   │ 63.677 genes × 8 amostras    │
│  genes)       │   │                              │
└───────────────┘   └──────────────────────────────┘
```

**2.1 Ficha das amostras.** As colunas mais importantes:

| Coluna | Significado |
|---|---|
| `cell` | Doador das células (N61311, N052611, N080611, N061011) |
| `dex` | Tratamento: `untrt` = controle; `trt` = tratado com dexametasona |
| `avgLength` | Comprimento médio das leituras |
| `SampleName`, `Run` e outras | Códigos de catálogo no GEO e no SRA |

👀 **Repare:** cada doador aparece duas vezes, uma como controle e outra como tratado. É um **desenho pareado**, com **4 réplicas biológicas por condição**.

**2.2 Tabela de contagens.** Cada linha é um gene (código ENSG do Ensembl) e cada coluna é uma amostra. Cada número é quantas leituras caíram naquele gene, naquela amostra.

- Um gene com valores na casa das centenas está ativo.
- O gene ENSG00000000005 tem **zero em todas as amostras**: é o *TNMD*, típico de tendão, desligado nessas células.

⚠️ **Não compare esses números diretamente.** As amostras foram sequenciadas com profundidades diferentes.

**2.3 Tamanho do experimento.** O total de leituras por amostra varia de cerca de 15 a 31 milhões. Essa diferença de quase 2 vezes é o que a normalização vai corrigir. A contagem por tipo de gene mostra que nem todo "gene" codifica proteína: há também RNAs não codificantes.

### Seção 3: Desenho experimental

Três passos:
1. Definimos **"não tratado"** como referência. Assim, **log2FoldChange positivo = mais expresso no tratado**.
2. Informamos o desenho `~ cell + dex`. O DESeq2 **desconta a diferença natural entre os doadores** e isola o efeito do tratamento, como num ensaio de "antes e depois" no mesmo corpo de prova.
3. Removemos genes com menos de 10 leituras no total. Sobram **22.369 genes**.

### Seção 4: Controle de qualidade com PCA

Antes de procurar genes, verificamos se as amostras se comportam como esperado.

**O que é PCA.** Cada amostra é descrita pelos 500 genes que mais variam: é um ponto num espaço de 500 dimensões, impossível de visualizar. A PCA (Análise de Componentes Principais) encontra as direções em que as amostras mais diferem e projeta tudo num plano.
> Analogia: a sombra de um objeto 3D na parede. A PCA escolhe o ângulo de iluminação que deixa a sombra mais informativa.

**O que você deve ver:**
- **PC1 (eixo horizontal), 41% da variação:** separa **não tratados (esquerda)** de **tratados (direita)**. O tratamento é a maior fonte de variação do experimento.
- **PC2 (eixo vertical), 26% da variação:** separa os **doadores**. Cada pessoa mantém sua "assinatura" mesmo após o tratamento.
- Nenhuma amostra isolada fora do seu grupo: sinal de dados limpos.

**4.2 Quais genes formam os eixos?** Cada eixo é uma combinação dos 500 genes, cada um com um **peso**. Os 10 genes com maior peso no PC1 (o eixo do tratamento) são:

> ZBTB16, SPARCL1, GPX3, SAMHD1, **FKBP5**, KLF15, STEAP4, ANGPTL7, MAOA, **TSC22D3**

💡 **Repare em algo notável:** a PCA não sabia quem era tratado nem quais genes procurar. Mesmo assim, entre os genes que mais pesam estão o **FKBP5** e o **TSC22D3** (também chamado GILZ, "zíper de leucina induzido por glicocorticoide"), marcadores clássicos da ação dos corticoides. O método "redescobriu" sozinho a biologia conhecida.

⚠️ Os pesos indicam **candidatos**. A lista oficial de genes alterados vem do DESeq2.

**4.3 (opcional)** mostra que os 4 primeiros eixos explicam quase 97% da variação, e que PC2, PC3 e PC4 refletem principalmente as diferenças entre os 4 doadores.

### Seção 5: Análise de expressão diferencial

A função `DESeq` faz tudo de uma vez. As mensagens na tela mostram cada etapa:

| Mensagem | O que acontece |
|---|---|
| `estimating size factors` | **Normalização:** calcula um fator por amostra para corrigir a profundidade |
| `estimating dispersions` | **Ruído:** mede quanto cada gene varia entre réplicas |
| `fitting model and testing` | **Teste:** verifica, gene a gene, se o tratamento mudou a expressão além do ruído |

> **Como o DESeq2 lida com poucas réplicas:** estimar o ruído de um gene com só 4 amostras é instável. Mas há 20 mil genes. O DESeq2 usa o comportamento típico de todos eles para corrigir a estimativa de cada um. É como estimar a incerteza de um sensor usando o comportamento de milhares de sensores do mesmo tipo.

**O placar (`summary`):**

| Resultado | Valor |
|---|---|
| Genes analisados | 22.369 |
| Subiram com o tratamento (padj < 0,05) | 2.208 (9,9%) |
| Desceram com o tratamento (padj < 0,05) | 1.822 (8,1%) |
| Genes descartados por valores discrepantes | 0 |
| Genes com expressão baixa demais para o teste | 5.204 (23%) |

**Conclusão:** a dexametasona alterou de forma estatisticamente confiável cerca de **4.030 genes**, 18% dos analisados. É uma resposta ampla, coerente com o fato de os corticoides agirem diretamente sobre a expressão gênica.

**5.1 Fatores de normalização:** um número por amostra. Acima de 1 = amostra com mais leituras que a média. Variam de 0,67 a 1,40.

### Seção 6: Resultados, três perguntas diferentes

Uma lista de genes pode ser ordenada de várias formas, e cada uma responde a uma pergunta diferente.

**6.1 Quais mudanças são mais certas? (ordenado por padj)**

Os 10 primeiros: SPARCL1, CACNB2, DUSP1, SAMHD1, MAOA, GPX3, STEAP2, NEXN, MT2A, ADAMTS1. Todos **subiram**, entre cerca de 4 e 24 vezes, com p-valores extremamente pequenos.

⚠️ Um padj menor **não** significa que o gene é mais importante biologicamente. Significa apenas mais certeza de que a mudança é real.

**6.2 Quais genes mudaram mais? (ordenado por log2FoldChange)**

| Gene | baseMean | log2FoldChange | Aumento | lfcSE |
|---|---|---|---|---|
| ALOX15B | 67 | 9,51 | ≈ 730× | 1,05 |
| ZBTB16 | 385 | 7,35 | ≈ 160× | 0,54 |
| RP11-357D18.1 | 56 | 6,33 | ≈ 80× | 0,68 |
| STEAP4 | 286 | 5,21 | ≈ 37× | 0,49 |
| RP11-434D9.1 | 9 | 5,09 | ≈ 34× | 1,15 |
| ANGPTL7 | 255 | 5,08 | ≈ 34× | 0,76 |
| PRODH | 55 | 4,89 | ≈ 30× | 0,52 |
| LINC00890 | 7 | 4,88 | ≈ 30× | 1,19 |
| FAM107A | 160 | 4,84 | ≈ 29× | 0,33 |
| LGI3 | 21 | 4,76 | ≈ 27× | 0,80 |

👀 **Repare:**
- O **ZBTB16**, primeiro da PCA, aparece aqui com aumento de cerca de 160 vezes.
- Genes com **baseMean baixo** (7 a 9 leituras) e **lfcSE alto** (acima de 1) têm mudanças enormes, mas incertas. Sair de 1 para 30 leituras já é "30 vezes".
- O **FAM107A** mudou cerca de 29 vezes com incerteza muito baixa (0,33): é um resultado bem mais confiável.
- Nomes como **RP11-...** e **LINC...** indicam, em geral, RNAs que não viram proteína.

> **Analogia:** dois sensores registram "aumento de 30 vezes". Um saiu de 1 para 30 unidades, perto do limite de detecção. O outro saiu de 100 para 3.000, com medições estáveis. Você confia muito mais no segundo.

**6.3 Quais mudaram mais de forma confiável? (log2FoldChange encolhido)**

A função `lfcShrink` "encolhe" em direção a zero as mudanças incertas e preserva as bem medidas. Compare com a lista 6.2: genes com lfcSE alto devem cair na classificação.

**6.4 Exigindo uma mudança mínima.** Com `lfcThreshold = 1`, o DESeq2 testa se o gene mudou **mais do que o dobro**, e não apenas "se mudou". Compare o novo placar com o da Seção 5.

**Resumo:**

| Lista | Pergunta que responde | Cuidado |
|---|---|---|
| Ordenada por padj | Quais mudanças são mais certas? | Mistura mudanças grandes e pequenas |
| Ordenada por log2FoldChange | Quais genes mudaram mais? | Favorece genes pouco expressos |
| Ordenada por log2FoldChange encolhido | Quais mudaram mais de forma confiável? | Melhor para escolher candidatos |

Os melhores candidatos para um artigo são os que aparecem bem em várias listas: mudança grande, expressão razoável e alta certeza.

### Seção 7: Visualização

- **MA plot:** os pontos coloridos são os genes significativos. Acima da linha zero, subiram; abaixo, desceram.
- **Volcano plot:** em vermelho, genes com padj < 0,05 **e** mudança acima de 2 vezes. Quanto mais alto o ponto, maior a certeza; quanto mais distante do centro, maior a mudança. Os 10 genes mais significativos aparecem com o nome.

### Seção 8: Exportar os resultados

Gera a planilha `resultados_DESeq2.csv`, com todos os genes e seus nomes. Para baixar: clique no ícone de pasta na barra lateral do Colab, depois nos três pontinhos ao lado do arquivo → *Fazer download*. O arquivo abre no Excel.

---

## 7. Limitações do experimento

Pensar criticamente nos dados é parte da análise. Neste estudo:

- **Cultura de células:** o resultado mostra a resposta de células isoladas, não de um pulmão inteiro, com sistema imune e outros tecidos.
- **Uma dose e um tempo:** só 1 micromolar por 18 horas. Genes de resposta mais lenta ou dependentes da dose podem não aparecer.
- **Quatro doadores:** suficiente para a estatística, mas pequeno para generalizar para a população.
- **Genoma antigo:** os dados usam uma versão de 2014 do genoma humano (GRCh37, Ensembl 75). Alguns genes podem ter mudado de nome ou posição no site atual.

---

## 8. Para o seu trabalho final

O RNA-seq pode ser uma das técnicas do seu artigo. Dados públicos estão disponíveis em:

- **GEO** (Gene Expression Omnibus): estudos de expressão gênica, frequentemente com tabelas de contagens prontas;
- **SRA** (Sequence Read Archive): os arquivos brutos de sequenciamento.

Ideias de perguntas ligadas à defesa: como uma bactéria responde a um antimicrobiano? Como células humanas respondem a um composto tóxico? Que genes mudam após exposição à radiação?

Uma combinação natural é com a **proteômica** (Aula 5): o mRNA de um gene subiu, mas a proteína também subiu?

---

## 9. Exercícios

1. Na Seção 6.4, troque `lfcThreshold = 1` por `lfcThreshold = 2`. Quantos genes continuam significativos? O que isso diz sobre a resposta ao corticoide?
2. Compare as listas das Seções 4.2 (pesos do PC1), 6.1 e 6.3. Quais genes aparecem em mais de uma? Por quê?
3. Escolha um gene da Seção 6.3, pesquise o código ENSG em [ensembl.org](https://www.ensembl.org) e descreva sua função em 3 linhas.
4. **Desafio:** repita a análise com o pacote **edgeR** e compare os 10 primeiros genes com os do DESeq2.

---

## 10. Glossário

| Termo | Definição |
|---|---|
| **Anotação (GTF)** | Arquivo que informa onde começa e termina cada gene e cada éxon no genoma |
| **BAM / SAM** | Arquivo com a posição de cada leitura no genoma (BAM é a versão compactada) |
| **Biomarcador** | Sinal mensurável que indica um estado biológico, como exposição a um agente |
| **Bioconductor** | Repositório de pacotes da linguagem R para análise de dados biológicos |
| **cDNA** | DNA complementar, cópia estável do RNA feita pela transcriptase reversa |
| **Desenho pareado** | Experimento em que cada indivíduo fornece uma amostra para cada condição |
| **Dispersão** | Medida de quanto um gene varia entre réplicas da mesma condição |
| **Éxon** | Trecho do gene que permanece no mRNA |
| **Expressão diferencial** | Diferença de atividade de um gene entre duas condições |
| **FASTQ** | Arquivo de texto com as leituras e a qualidade de cada letra |
| **Gene** | Trecho de DNA com a instrução para produzir um RNA ou uma proteína |
| **Genoma** | Conjunto completo do DNA de um organismo |
| **Genoma de referência** | Genoma já sequenciado e publicado, usado como "foto da caixa" |
| **Índice** | Estrutura que permite buscar rapidamente sequências no genoma |
| **Íntron** | Trecho do gene removido do mRNA durante o splicing |
| **Leitura (read)** | Fragmento curto de sequência lido pelo sequenciador |
| **log2FoldChange** | Tamanho da mudança de um gene, em escala logarítmica de base 2 |
| **Mapeamento** | Localização de cada leitura no genoma de referência |
| **Montagem *de novo*** | Reconstrução dos transcritos sem genoma de referência |
| **mRNA** | RNA mensageiro, cópia de um gene usada para fabricar proteína |
| **Normalização** | Correção das diferenças de profundidade de sequenciamento entre amostras |
| **padj** | p-valor ajustado para o grande número de testes |
| **PCA** | Análise de Componentes Principais: resume muitas variáveis em poucos eixos |
| **Réplica biológica** | Amostra independente da mesma condição (outro organismo ou cultura) |
| **RNA-seq** | Técnica que sequencia o RNA para medir a atividade de milhares de genes |
| **Splicing** | Remoção dos íntrons e união dos éxons no RNA |
| **Tabela de contagens** | Tabela com o número de leituras por gene e por amostra |
| **Transcrição** | Cópia de um gene do DNA em RNA |
| **Transcriptoma** | Conjunto de RNAs presentes numa célula num dado momento |
| **Tradução** | Leitura do mRNA para fabricar uma proteína |

---

## 11. Para saber mais

- Himes BE et al. *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.* PLoS One, 2014. (O artigo que gerou os dados da prática.)
- Documentação oficial do **DESeq2** e do **edgeR** (site do Bioconductor).
- Documentação oficial do **HISAT2**, do **STAR** e do **Trinity**.
- Portal **Ensembl**: [ensembl.org](https://www.ensembl.org).
- Pevsner, J. *Bioinformatics and Functional Genomics.* Wiley-Blackwell. (Bibliografia da disciplina.)
