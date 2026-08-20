# Aula 2 — Qualidade e Processamento de Dados de NGS; Identificação e Comparação de Sequências Biológicas

**Disciplina:** Introdução à Bioinformática — PPGED/IME
**Carga horária da aula:** 3 horas (teoria + prática guiada)
**Público-alvo:** Pós-graduandos em Engenharia de Defesa, sem experiência prévia em programação ou biologia molecular aprofundada

---

## Ideia-guia da aula

Vamos usar uma metáfora que vai acompanhar toda a aula:

> **Um genoma sequenciado é como um livro despedaçado.** O sequenciador não lê o "livro" (o DNA) inteiro de uma vez — ele produz milhões de pedacinhos de texto (as *reads*), com erros de leitura aqui e ali. Nosso trabalho como bioinformatas é: (1) checar a qualidade de cada pedacinho, (2) jogar fora ou aparar o que está ruim, (3) remontar o quebra-cabeça comparando com um "livro de referência" já conhecido, e (4) buscar por trechos parecidos em bibliotecas mundiais de sequências.

Cada ferramenta desta aula resolve uma dessas quatro etapas.

---

## Parte 1 — TEORIA (1,5h)

### 1.1 Por que isso importa para Engenharia de Defesa 

Antes de entrar nas ferramentas, contextualizamos o "porquê":

- **Biovigilância e biodefesa:** identificar rapidamente um patógeno (vírus, bactéria) a partir de amostras ambientais ou clínicas.
- **Identificação forense biológica:** comparar material genético encontrado em campo com bancos de dados de referência.
- **Verificação de agentes biológicos:** distinguir uma cepa natural de uma manipulada, comparando sequências.
- **Análise de proteínas de interesse:** para quem trabalha com química de proteínas e macromoléculas, entender a sequência é o primeiro passo antes de qualquer análise estrutural (que veremos na Unidade 7).

### 1.2 O que é uma "read" e o formato FASTQ

- Sequenciadores modernos (NGS — *Next Generation Sequencing*) não leem uma molécula de DNA inteira. Eles quebram o material genético em milhões de fragmentos pequenos e leem cada um separadamente. Cada fragmento lido é uma **read**.
- O resultado é salvo em um arquivo de texto no formato **FASTQ**. Vamos abrir um arquivo real e explicar linha por linha:

```
@SRR000001.1 primeiro identificador da read
GATTACAGGACATTACAGGATTACA
+
IIIIIIIIIIIIII!!!!IIIIIII
```

  1. **Linha 1 (`@...`):** identificador único da read (qual sequenciador, qual "poço", qual leitura).
  2. **Linha 2:** a sequência de bases lida (A, T, C, G).
  3. **Linha 3 (`+`):** apenas um separador.
  4. **Linha 4:** a **qualidade** de cada base lida, codificada em caracteres (quanto mais para o fim do alfabeto ASCII, melhor a qualidade).

- **Phred Score:** cada caractere da linha 4 representa uma nota de confiança (0 a ~40) de que aquela base foi lida corretamente. Analogia: é como o "nível de certeza" de quem está ditando um texto por telefone com uma ligação ruidosa — quanto maior a nota, mais confiável a "letra" que foi ouvida.

#### Aprofundando: por que a escala vai só até ~40, e por que aparecem letras como "I"?

Essa dúvida costuma surgir na primeira vez que se olha um arquivo FASTQ de verdade, então vale a pena parar e explicar em duas partes.

**Parte 1 — a escala é logarítmica, não vai de 0 a 100.** O Phred Score não é um número arbitrário: ele é calculado por uma fórmula que relaciona a nota à probabilidade de erro na leitura daquela base:

$$Q = -10 \times \log_{10}(P_{erro})$$

Alguns valores de referência:

| Phred Score (Q) | Probabilidade de erro | Em palavras |
|---|---|---|
| Q10 | 1 em 10 | 90% de acerto |
| Q20 | 1 em 100 | 99% de acerto |
| Q30 | 1 em 1.000 | 99,9% de acerto |
| Q40 | 1 em 10.000 | 99,99% de acerto |

O "teto" perto de 40 não é uma limitação matemática (a escala poderia seguir subindo) — é uma limitação **física/tecnológica** dos sequenciadores atuais (como os da Illumina): eles não conseguem, na prática, garantir uma confiança maior que "1 erro a cada 10.000 bases" de forma sustentável. Reportar Q60, por exemplo, seria uma precisão que o instrumento não sustenta de verdade — então os fabricantes "saturam" a escala ali por volta de 40-41, como se dissessem "este é o nosso máximo de confiança possível".

**Parte 2 — por que aparecem letras (como "I") em vez de números?** Cada Phred Score é convertido em **um único caractere ASCII**, para ocupar menos espaço no arquivo. A regra de conversão mais usada hoje (padrão **Phred+33**, também chamado "Sanger") é simplesmente:

$$\text{código ASCII do caractere} = Q + 33$$

Ou seja:
- Q = 0 → ASCII 33 → o caractere **`!`**
- Q = 40 → ASCII 73 → o caractere **`I`**

No exemplo de FASTQ mostrado acima, a linha de qualidade foi escrita de propósito mostrando as duas pontas da escala:
```
IIIIIIIIIIIIII!!!!IIIIIII
```
Os `I` representam bases lidas com qualidade máxima (Q40, ótima confiança), e os `!` no meio representam bases de qualidade péssima (Q0) — ilustrando que uma mesma read pode ter trechos muito bons e um trecho ruim no meio. É exatamente esse tipo de padrão que o FastQC detecta nos seus gráficos, e que o Trimmomatic depois corta (Seção 1.4).

*Dica para a aula:* quando os alunos abrirem um arquivo FASTQ real no Galaxy, essa sequência de símbolos vai parecer estranha à primeira vista — reforçar que cada caractere ali é, na verdade, um número escondido representando a confiança do sequenciador naquela letra específica.

### 1.3 Controle de qualidade: FastQC e MultiQC

- **FastQC**
  O que é: um programa que faz um "raio-X" de um arquivo FASTQ inteiro e gera um relatório visual (HTML) com gráficos prontos: qualidade média por posição na read, conteúdo de GC, presença de adaptadores de sequenciamento, sequências duplicadas, entre outros.
  Por que importa: em vez de olhar milhões de linhas de texto, o pesquisador olha um punhado de gráficos e decide rapidamente se os dados estão bons o suficiente para prosseguir.
  Como se usa: aponta-se o arquivo FASTQ, clica em "rodar", e o relatório é gerado automaticamente — não é necessário programar.

- **MultiQC**
  O que é: uma ferramenta que agrega os relatórios de várias amostras (por exemplo, 50 relatórios do FastQC) em um único painel comparativo.
  Por que importa: em estudos de vigilância epidemiológica ou biodefesa, é comum processar dezenas ou centenas de amostras ao mesmo tempo — o MultiQC evita que o analista tenha que abrir um relatório de cada vez.

### 1.4 Limpeza dos dados: Trimmomatic e Cutadapt

- Depois do diagnóstico (FastQC), muitas vezes descobrimos que as pontas das reads têm qualidade baixa, ou que sobraram pedaços de "adaptadores" (sequências artificiais usadas no processo químico de sequenciamento, que não fazem parte do organismo).
- **Trimmomatic** e **Cutadapt**
  O que fazem: funcionam como uma "tesoura de aparar bordas" — cortam as partes de baixa qualidade nas extremidades das reads e removem sequências de adaptador contaminantes.
  Analogia: é como aparar as bordas tremidas ou fora de foco de uma fotografia, mantendo só a parte nítida e útil.
  Resultado: um novo arquivo FASTQ, mais "limpo", pronto para as etapas seguintes.
- Prática recomendada: rodar o FastQC de novo *depois* da limpeza, para confirmar visualmente que a qualidade melhorou (comparação "antes x depois").

### 1.5 Alinhamento contra um genoma de referência: BWA, Bowtie2 e SAMtools

- Depois de limpos, os pedacinhos (reads) precisam ser "encaixados" de volta em um genoma de referência conhecido — como resolver um quebra-cabeça enorme, comparando cada peça com a imagem da caixa.
- **BWA** (Burrows-Wheeler Aligner) e **Bowtie2**
  O que fazem: recebem as reads e um genoma de referência, e encontram a posição exata (ou mais provável) de onde cada read veio dentro desse genoma.
  Por que dois programas parecidos existem: ambos resolvem o mesmo problema com algoritmos ligeiramente diferentes; BWA é muito usado em genômica humana/clínica, Bowtie2 é muito usado em RNA-seq e metagenômica — na prática, os alunos vão ver que a lógica de uso é semelhante.
- **SAMtools**
  O que é: o "canivete suíço" para manipular os arquivos gerados pelo alinhamento (formatos **SAM** e **BAM** — basicamente uma tabela de "qual read caiu em qual posição do genoma"). Serve para ordenar, filtrar, indexar e converter esses arquivos, preparando-os para visualização ou análises posteriores.

### 1.6 Mudando de escala: busca e comparação de sequências com BLAST, Clustal Omega e MUSCLE

- Até agora falamos de alinhar *milhões* de reads pequenas contra *um* genoma de referência específico. Agora mudamos de problema: e se eu tiver **uma única sequência** (um gene, uma proteína) e quiser saber **quais sequências parecidas já existem no mundo**?
- **BLAST** (Basic Local Alignment Search Tool)
  O que é: a ferramenta mais usada em bioinformática para comparar uma sequência de interesse contra bancos de dados mundiais (GenBank, UniProt, RefSeq), retornando os "parentes" mais parecidos, com um escore de similaridade.
  Variantes: `blastn` (DNA contra DNA), `blastp` (proteína contra proteína), `blastx` (DNA traduzido contra proteínas).
  Como usar sem programar: existe uma versão totalmente online, no site do NCBI — basta colar a sequência e clicar em "buscar".
- **Clustal Omega** e **MUSCLE**
  O que fazem: enquanto o BLAST compara uma sequência contra um banco gigante, essas ferramentas alinham **várias sequências ao mesmo tempo**, lado a lado, destacando regiões conservadas (iguais) e variáveis entre elas — útil, por exemplo, para comparar a mesma proteína em espécies diferentes.

---

## Parte 2 — PRÁTICA GUIADA (1,5h, síncrona) + Atividade assíncrona

### Por que usar o Galaxy nesta turma

Como a turma não tem experiência prévia com linha de comando, toda a prática síncrona será conduzida na plataforma **Galaxy** (https://usegalaxy.org), uma interface web gratuita do tipo "clicar e escolher": o aluno sobe o arquivo, seleciona a ferramenta (FastQC, Trimmomatic, BWA, SAMtools) em um menu, e vê o resultado — sem escrever código. O terminal Linux e o Python/Biopython aparecerão apenas como demonstração opcional de observação, não como exigência.

### Roteiro passo a passo da prática síncrona (professor demonstra e alunos replicam)

**Passo 0 — Preparação**
1. Criar conta gratuita em https://usegalaxy.org (pode ser feito antes da aula).
2. Abrir uma nova "History" (histórico de trabalho) com um nome identificável, ex: `Aula2_Turma_PPGED`.

**Passo 1 — Obter um dataset público de exemplo**
1. Acessar o **SRA (Sequence Read Archive)**, o repositório mundial de dados brutos de sequenciamento, mantido pelo NCBI.
2. O professor já disponibiliza previamente o identificador de uma amostra pequena e didática (ver quadro "Amostras recomendadas" abaixo), para não gastar tempo de aula com downloads grandes.
3. Importar esse dataset diretamente para o Galaxy usando a ferramenta **"Faster Download and Extract Reads in FASTQ" (fastq-dump)** (link na tabela de materiais). Basta colar o *accession* desejado (ex.: `SRR957824` ou `SRR453566`) no campo de entrada da ferramenta e executar — o Galaxy baixa e converte os dados automaticamente para o formato FASTQ.

#### Amostras recomendadas (procarioto e eucarioto)

Os dados do SRA são sempre **brutos**, exatamente como saíram do sequenciador — ou seja, **não vêm trimados nem limpos**. Isso é bom para a aula: é justamente essa "sujeira" (queda de qualidade nas pontas, resíduo de adaptador) que o FastQC vai revelar e o Trimmomatic vai corrigir.

| | Procarioto | Eucarioto |
|---|---|---|
| **Organismo** | *Escherichia coli* O157 (EHEC) | Levedura (*Saccharomyces cerevisiae*) |
| **Accession (SRA)** | `SRR957824` | `SRR453566` |
| **Contexto biológico** | Cepa causadora de um surto real de doença gastrointestinal grave, investigado em Saint Louis (EUA, 2011) — ótimo gancho para discutir biovigilância e investigação de surtos | Organismo-modelo mais estudado da biologia depois da *E. coli*; genoma pequeno para um eucarioto (~12 Mb) e muito bem anotado no Ensembl |
| **Tamanho do genoma** | ~5 Mb (bactéria) | ~12 Mb |
| **Por que serve bem para a aula** | Genoma pequeno, processa rápido; dados conhecidamente com queda de qualidade nas pontas e resíduo de adaptador — ótimo para o exercício de FastQC → Trimmomatic | Introduz conceitos de eucarioto (núcleo, cromossomos lineares, alguns genes com íntrons) sem a complexidade de um genoma gigante como o humano |
| **Genoma de referência sugerido** | Disponível no NCBI Genome / Ensembl Bacteria | Disponível no Ensembl Fungi |
| **Uso futuro na disciplina** | Pode ser retomado na Unidade 3 (Genômica) | Pode ser retomado na Unidade 4 (Transcriptômica), já que a amostra é originalmente de RNA-seq |

*Observação:* como os arquivos completos podem ter alguns milhões de reads, o professor prepara previamente uma versão reduzida (subconjunto de algumas dezenas de milhares de reads) de cada amostra, para caber confortavelmente no tempo da aula síncrona.

**Passo 2 — Rodar o FastQC**
1. No menu de ferramentas do Galaxy, buscar por "FastQC".
2. Selecionar o arquivo FASTQ importado e executar.
3. Abrir o relatório HTML gerado e interpretar juntos, em turma:
   - O gráfico de qualidade por posição está caindo nas pontas?
   - Há sinal de adaptadores?
   - O conteúdo de GC está dentro do esperado para o organismo?

**Passo 3 — Limpeza com Trimmomatic**
1. Buscar a ferramenta "Trimmomatic" no Galaxy.
2. Configurar parâmetros básicos (corte por qualidade mínima, remoção de adaptadores) — o professor explica cada parâmetro em linguagem simples, sem jargão excessivo.
3. Rodar sobre o dataset original.

**Passo 4 — Comparar antes x depois**
1. Rodar o FastQC novamente sobre o arquivo já limpo pelo Trimmomatic.
2. Colocar os dois relatórios lado a lado e discutir as diferenças visuais.

**Passo 5 — Alinhamento com Bowtie2**
1. Buscar a ferramenta "Bowtie2" no Galaxy (ver observação sobre BWA x Bowtie2 abaixo).
2. Selecionar o genoma de referência correspondente à amostra usada (ver quadro "Genomas de referência" abaixo).
3. Rodar o alinhamento das reads limpas contra essa referência.
4. Abrir o arquivo BAM resultante e observar, de forma exploratória, o conceito de "cobertura" (quantas reads caíram em cada posição do genoma).

#### Genomas de referência recomendados

| | Procarioto (SRR957824) | Eucarioto (SRR453566) |
|---|---|---|
| **Organismo** | *Escherichia coli* O157:H7 str. Sakai | Levedura (*Saccharomyces cerevisiae*) |
| **Accession RefSeq (cromossomo)** | `NC_002695.1` | — |
| **Assembly completo** | `GCF_000008865.2` (ASM886v2) | `R64-1-1` |
| **Onde obter** | NCBI Genome / Ensembl Bacteria | Ensembl Fungi |
| **Observação** | A cepa Sakai é a referência-padrão mundial para *E. coli* O157:H7 e é a mesma linhagem do surto que originou a amostra — coerência total entre amostra e referência. Existem ainda dois plasmídeos associados (pO157 e pOsak1), mas não são necessários para este exercício introdutório | Assembly padrão e mais usado para esse organismo-modelo |

#### BWA ou Bowtie2 — por que usamos Bowtie2 nesta aula

Não existe um "vencedor" absoluto entre as duas ferramentas — a escolha depende do tipo de dado:

| Critério | BWA | Bowtie2 |
|---|---|---|
| Tipo de dado mais adequado | Reads mais longas, incluindo dados com indels maiores (ex.: PacBio, ou Illumina com reads mais longas) | Reads curtas, típicas de Illumina (nosso caso: 2×150bp) |
| Uso mais comum em RNA-seq | Menos usado diretamente | Muito usado (base de ferramentas como HISAT2/Tophat, revisitadas na Unidade 4) |
| Facilidade de configuração no Galaxy | Boa, mas com mais parâmetros técnicos | Interface um pouco mais simples para iniciantes |

Optamos pelo **Bowtie2** nesta aula porque: (1) nosso dataset é Illumina de reads curtas, cenário em que o Bowtie2 é considerado referência; (2) ele antecipa uma ferramenta-irmã (HISAT2) que os alunos vão rever na Unidade 4; (3) tende a exigir menos parâmetros técnicos no primeiro contato. Vale mencionar aos alunos, como nota conceitual, que na prática de mercado os dois são frequentemente usados de forma intercambiável para esse tipo de dado — a escolha muitas vezes reflete apenas o pipeline já estabelecido no laboratório.

**Passo 6 — Manipulação com SAMtools**
1. Usar a ferramenta "Samtools flagstat" no Galaxy para gerar um resumo estatístico simples do alinhamento (quantas reads alinharam, quantas não alinharam).
2. Interpretar juntos o resultado.

**Passo 7 — Busca no BLAST**
1. Acessar diretamente o site do NCBI BLAST (https://blast.ncbi.nlm.nih.gov).
2. Colar uma sequência de proteína ou gene de interesse (o professor sugere uma sequência ligada a alguma linha de pesquisa da turma, ex.: uma proteína usada em química de macromoléculas).
3. Rodar `blastp` e interpretar a tabela de resultados: percentual de identidade, valor de *E-value* (quanto menor, mais significativa a similaridade encontrada), organismo de origem do "parente" encontrado.

**Passo 8 — Alinhamento múltiplo com Clustal Omega**
1. Pegar 3-4 sequências homólogas (parecidas) encontradas no passo anterior via BLAST.
2. Colar no Clustal Omega online (https://www.ebi.ac.uk/jdispatcher/msa/clustalo).
3. Observar visualmente as regiões conservadas (coloridas de forma idêntica entre as sequências) e as regiões variáveis.

**Passo 9 (complementar/opcional) — Análise de arquivos Sanger (.ab1): qualidade Phred e montagem de contig**

Esta etapa é opcional e serve para conectar o conceito de Phred Score (Seção 1.2) com sua origem histórica no sequenciamento Sanger, além de introduzir a lógica de montagem de contigs em miniatura (2 leituras), antes da montagem em larga escala vista na Unidade 3.

1. **Obter os dados brutos:** baixar um par de arquivos `.ab1` (forward e reverse) do repositório público no Zenodo:
   👉 https://zenodo.org/records/7104640
   Esse conjunto de dados foi compartilhado abertamente pelos autores de um artigo científico sobre variantes genéticas associadas a distonia, e é o mesmo dataset usado no tutorial oficial do Galaxy Training Network para manejo de arquivos AB1.

2. **Analisar a qualidade (Phred) do cromatograma:** acessar a calculadora/visualizador do ConductScience e fazer upload de um dos arquivos `.ab1` baixados:
   👉 https://conductscience.com/tools/file-formats/ab1-chromatogram#tool-calculator
   Observar o eletroferograma (picos coloridos, um por base) e os valores de Phred Score exibidos por posição — os alunos devem identificar visualmente as regiões de baixa qualidade nas pontas da leitura (mesmo padrão discutido na Seção 1.2, agora aplicado a uma leitura Sanger longa em vez de reads curtas de NGS).

3. **Montar o contig (consenso):** com as duas leituras (forward + reverse) já identificadas, colar/enviar as sequências (formato FASTA) no CAP3 do PRABI-Doua:
   👉 https://doua.prabi.fr/cgi-bin/run_cap3
   O CAP3 vai alinhar as duas leituras, identificar a região de sobreposição e gerar a sequência consenso (contig) final.

4. **Discussão em turma:** comparar esse processo com o alinhamento em larga escala feito nos Passos 5-6 (Bowtie2/SAMtools) — a lógica de "sobrepor e resolver discordâncias usando qualidade" é a mesma, só que aqui em escala reduzida e visualmente mais fácil de acompanhar passo a passo.

---

## Atividade prática assíncrona (entre a Aula 2 e a Aula 3)

**Entrega:** resumo de 1 página (PDF ou Word), com prints de tela e interpretação em texto próprio.

**Roteiro para o aluno:**
1. Escolher (ou usar o mesmo da aula) um dataset pequeno do SRA.
2. Rodar o FastQC nesse dataset via Galaxy e escrever 3-4 linhas interpretando o relatório (qualidade geral, presença de adaptadores).
3. Escolher uma sequência de proteína ou gene relacionada à sua própria linha de pesquisa (química de proteínas/macromoléculas) e rodar uma busca BLAST online.
4. Registrar: qual foi o "parente" mais próximo encontrado, qual o percentual de identidade e o *E-value*, e uma frase de interpretação sobre o que isso significa biologicamente.
5. Enviar o resumo com prints de cada etapa.

**Critério de avaliação da atividade:** clareza da interpretação (não a perfeição técnica), demonstrando que o aluno entendeu o *significado* de cada resultado — não apenas que "rodou o programa".

---

## Materiais e links de apoio

| Ferramenta | Acesso | Necessita instalação? |
|---|---|---|
| Galaxy | https://usegalaxy.org | Não (web) |
| SRA (NCBI) | https://www.ncbi.nlm.nih.gov/sra | Não (web) |
| Ferramenta fastq-dump (download SRA no Galaxy) | https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fsra_tools%2Ffastq_dump%2F3.1.1%2Bgalaxy1&version=latest | Não (web) |
| FastQC / MultiQC | Via Galaxy — https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fdevteam%2Ffastqc%2Ffastqc%2F0.74%2Bgalaxy1&version=latest | Não |
| Trimmomatic / Cutadapt | Via Galaxy — https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fpjbriggs%2Ftrimmomatic%2Ftrimmomatic%2F0.39%2Bgalaxy2&version=latest | Não |
| BWA / Bowtie2 / SAMtools | Via Galaxy — Bowtie2: https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fdevteam%2Fbowtie2%2Fbowtie2%2F2.5.5%2Bgalaxy0&version=latest / SAMtools flagstat: https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fdevteam%2Fsamtools_flagstat%2Fsamtools_flagstat%2F2.0.8&version=latest | Não |
| BLAST | https://blast.ncbi.nlm.nih.gov | Não (web) |
| Clustal Omega | https://www.ebi.ac.uk/jdispatcher/msa/clustalo | Não (web) |
| Dataset de exemplo .ab1 (Zenodo) | https://zenodo.org/records/7104640 | Não (web) |
| Visualizador Phred/AB1 (ConductScience) | https://conductscience.com/tools/file-formats/ab1-chromatogram#tool-calculator | Não (web) |
| CAP3 — montagem de contig (PRABI-Doua) | https://doua.prabi.fr/cgi-bin/run_cap3 | Não (web) |

### Amostras de exemplo (SRA)

| Tipo | Organismo | Accession | Observação |
|---|---|---|---|
| Procarioto | *E. coli* O157 (EHEC) | SRR957824 | Dados brutos, com queda de qualidade nas pontas e resíduo de adaptador — ideal para o exercício de limpeza |
| Eucarioto | Levedura (*S. cerevisiae*) | SRR453566 | Dados brutos de RNA-seq; genoma pequeno e bem anotado (Ensembl Fungi) |

## Bibliografia de apoio para esta aula

- Pevsner, J. *Bioinformatics and Functional Genomics*. Wiley-Blackwell.
- Compeau, P.; Pevzner, P. *Bioinformatics Algorithms: An Active Learning Approach*.
- Mariano, D. *Python para Bioinformática: Fundamentos de Programação para Bioinformática e Biologia Computacional*. Novatec Editora, 2025.
- Documentação oficial: NCBI BLAST, SRA, FastQC, Trimmomatic.
