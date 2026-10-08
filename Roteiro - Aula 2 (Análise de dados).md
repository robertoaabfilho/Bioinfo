# Material de Apoio do Aluno — Aula 2
## Qualidade e Processamento de Dados de NGS · Identificação e Comparação de Sequências

**Disciplina:** Introdução à Bioinformática — PPGED/IME
**Unidade:** 2 · **Carga horária:** 4 h (2 h de teoria + 2 h de prática guiada) + atividade assíncrona
**Para quem é este material:** alunos de pós-graduação em Engenharia de Defesa, **sem** necessidade de saber programar e **sem** formação aprofundada em biologia.

---

## Como usar este material

- **Antes da aula:** leia a Parte 1 (teoria) e faça o checklist da Seção 0. Você chega mais confiante e aproveita melhor a prática.
- **Durante a aula:** deixe a Parte 2 (prática) aberta ao lado do Galaxy. Cada passo diz *o que fazer*, *o que você deve ver* e *o que fazer se der errado*.
- **Depois da aula:** use a Parte 3 para a atividade avaliativa e o Glossário (Parte 5) para tirar dúvidas de vocabulário.

> **Você não precisa escrever uma linha de código nesta aula.** Todas as ferramentas serão usadas por interface web (cliques e menus).

---

## 0. Checklist antes da aula

- [ ] Conta gratuita criada em **https://usegalaxy.org** (Menu *User → Register*). Use um e-mail que você acessa, pois o Galaxy envia confirmação.
- [ ] Navegador atualizado (Chrome, Firefox ou Edge) e conexão estável.
- [ ] Lista de ferramentas favoritas aberta (Seção 6 — Links).
- [ ] Uma sequência de proteína ou gene da **sua linha de pesquisa** salva em um arquivo de texto (vai ser usada no BLAST).
- [ ] Pasta no computador chamada `Aula2_Bioinfo` para guardar prints e arquivos baixados.

---

## Mapa da aula

A ideia que vai nos acompanhar do começo ao fim:

> **Um genoma sequenciado é como um livro despedaçado.**
> O sequenciador não lê o "livro" (o DNA) inteiro: ele produz **milhões de pedacinhos de texto** (as *reads*), com erros de leitura aqui e ali.

O nosso trabalho é um fluxo de quatro etapas:

| Etapa | Pergunta | Ferramentas |
|---|---|---|
| **1. Checar** | Os pedacinhos têm boa qualidade? | FastQC, MultiQC |
| **2. Aparar** | Como removo o que está ruim? | Trimmomatic, Cutadapt |
| **3. Remontar** | De onde no genoma cada pedacinho veio? | Bowtie2, BWA, SAMtools |
| **4. Buscar e comparar** | Que sequências parecidas já existem no mundo? | BLAST, Clustal Omega, MUSCLE |

```
FASTQ (dados brutos) → FastQC → Trimmomatic → FastQC de novo → Bowtie2 → SAMtools
                                                                       ↘
                         (outra pergunta, outra escala)   BLAST → Clustal Omega
```

### Objetivos de aprendizagem

Ao final desta aula você será capaz de:

1. Explicar o que é uma *read* e ler as quatro linhas de um registro FASTQ.
2. Interpretar o Phred Score e decodificar um caractere de qualidade.
3. Ler um relatório do FastQC e decidir se os dados precisam de limpeza.
4. Executar a limpeza com Trimmomatic e comprovar a melhora.
5. Alinhar reads a um genoma de referência com Bowtie2 e interpretar o `samtools flagstat`.
6. Explicar as colunas de um arquivo SAM, incluindo FLAG e CIGAR.
7. Fazer uma busca BLAST e interpretar identidade, cobertura e E-value.
8. Gerar e ler um alinhamento múltiplo com Clustal Omega.

---

# PARTE 1 — TEORIA

## 1.1 Por que isso importa para a Engenharia de Defesa

Antes das ferramentas, o "porquê":

- **Biovigilância e biodefesa:** identificar rapidamente um patógeno (vírus, bactéria) a partir de amostras ambientais ou clínicas.
- **Identificação forense biológica:** comparar material genético encontrado em campo com bancos de dados de referência.
- **Verificação de agentes biológicos:** comparar sequências para avaliar se uma cepa se afasta do que existe na natureza.
- **Química de proteínas e macromoléculas:** conhecer a sequência é o primeiro passo antes de qualquer análise estrutural (Unidade 7).

Em todos os casos, o ponto de partida são dados brutos de sequenciamento. **Se a qualidade dos dados for ruim e ninguém perceber, todas as conclusões seguintes ficam comprometidas.** Daí a ênfase em controle de qualidade nesta aula.

---

## 1.2 O que é uma *read* e o formato FASTQ

### Sequenciamento de nova geração (NGS)

Os sequenciadores modernos (NGS — *Next Generation Sequencing*) **não leem uma molécula de DNA de ponta a ponta**. Eles:

1. quebram o material genético em milhões de fragmentos pequenos;
2. leem cada fragmento separadamente;
3. entregam o resultado como uma enorme lista de textos curtos.

Cada fragmento lido é uma **read** (em português, "leitura").

### Reads *single-end* e *paired-end*

Um fragmento de DNA tem duas pontas. O sequenciador pode:

- ler **só uma ponta** → *single-end* (um arquivo);
- ler **as duas pontas** do mesmo fragmento → *paired-end* (dois arquivos).

Em *paired-end*, os dois arquivos são chamados de **R1** (também *forward*) e **R2** (também *reverse*). A *n*-ésima read de R1 e a *n*-ésima read de R2 vêm do **mesmo fragmento** e precisam permanecer sincronizadas durante toda a análise. Isso será importante na prática (Passo 5).

Analogia: é como fotografar a mesma página de um livro pelas duas margens. As duas fotos juntas dizem mais do que cada uma sozinha, inclusive sobre a distância entre elas.

### O formato FASTQ

O resultado fica em um arquivo de texto chamado **FASTQ**. Cada read ocupa **4 linhas**:

```
@SRR000001.1 identificador da read
GATTACAGGACATTACAGGATTACA
+
IIIIIIIIIIIIII!!!!IIIIIII
```

| Linha | Conteúdo | Para que serve |
|---|---|---|
| 1 (`@...`) | Identificador | Nome único da read (sequenciador, posição física, número da leitura) |
| 2 | Sequência | As bases lidas: A, T, C, G (e `N` quando o sequenciador não conseguiu decidir) |
| 3 (`+`) | Separador | Só marca o fim da sequência; às vezes repete o identificador |
| 4 | Qualidade | Uma "nota de confiança" para **cada** base da linha 2, em forma de caractere |

> A linha 4 tem **exatamente o mesmo número de caracteres** que a linha 2. O 1º caractere da linha 4 é a nota da 1ª base, e assim por diante.

### Phred Score: a nota de confiança de cada base

O **Phred Score (Q)** diz quão confiável é a leitura de uma base. Analogia: o nível de certeza de quem dita um texto por telefone numa ligação ruidosa. Quanto maior a nota, mais confiável a "letra" ouvida.

Ele é calculado a partir da **probabilidade de erro** daquela base:

$$Q = -10 \times \log_{10}(P_{erro})$$

| Phred (Q) | Probabilidade de erro | Em palavras |
|---|---|---|
| Q10 | 1 em 10 | 90% de acerto |
| Q20 | 1 em 100 | 99% de acerto |
| Q30 | 1 em 1.000 | 99,9% de acerto |
| Q40 | 1 em 10.000 | 99,99% de acerto |

Dois pontos para fixar:

- **A escala é logarítmica** (como decibéis ou a escala Richter). Cada +10 na nota significa 10 vezes menos erro. Por isso Q30 é **muito** melhor que Q20, não "só 10 pontos melhor".
- **Q30 é a referência de "boa qualidade"** mais usada na prática: 1 erro a cada mil bases.

**Por que a escala para perto de Q40?** Não é limite matemático (a fórmula permitiria Q50, Q60...). É limite **tecnológico**: os sequenciadores atuais não conseguem sustentar, de forma confiável, mais do que ~1 erro em 10.000 bases. Por isso os fabricantes limitam o valor máximo reportado, por volta de 40 a 41.

### Por que aparecem letras e símbolos como `I`, `!` e `?`?

Para economizar espaço, cada nota Q é gravada como **um único caractere**. A regra mais usada hoje (**Phred+33**, também chamada "Sanger") é:

$$\text{código ASCII do caractere} = Q + 33$$

| Caractere | Código ASCII | Q | Interpretação |
|---|---|---|---|
| `!` | 33 | 0 | Péssima |
| `#` | 35 | 2 | Muito baixa |
| `5` | 53 | 20 | Razoável (99%) |
| `?` | 63 | 30 | Boa (99,9%) |
| `I` | 73 | 40 | Excelente (máxima) |

No exemplo acima, `IIIIIIIIIIIIII!!!!IIIIIII` descreve uma read com trechos de qualidade máxima e um trecho péssimo no meio (`!!!!`).

> **Dúvida comum:** "Vi um `?` ao lado do `+`. Isso é um erro?" Não. O `+` está na linha 3 e o `?` na linha 4 (qualidade), uma embaixo da outra. O `?` é só mais um caractere da escala (Q30). **Qualquer símbolo estranho na linha 4 é apenas código de qualidade.** Você não precisa decorar a tabela: o FastQC faz a conversão em gráficos.

**Mini-exercício (2 min):** a linha de qualidade de uma read começa com `II5?#`. Qual é o Q de cada base?
*Resposta:* `I` = 40, `I` = 40, `5` = 20, `?` = 30, `#` = 2.

---

## 1.3 Controle de qualidade: FastQC e MultiQC

### FastQC

**O que é.** Um programa que faz um "raio-X" de um arquivo FASTQ e gera um relatório visual (HTML) com gráficos prontos. Em vez de ler milhões de linhas de texto, você olha alguns gráficos.

**Como se usa.** Você aponta o arquivo e clica em executar. Não precisa programar.

**Os módulos que mais importam:**

| Módulo do relatório | O que mostra | O que você procura |
|---|---|---|
| **Per base sequence quality** | Qualidade (Phred) em cada posição da read | Queda de qualidade no final da read é comum em Illumina; qualidade na zona vermelha (< Q20) em boa parte da read é um alerta |
| **Per sequence quality scores** | Distribuição da qualidade média por read | Um pico em Q30 ou mais é bom |
| **Per base sequence content** | Proporção de A, T, C, G em cada posição | Linhas aproximadamente paralelas; desvios nas primeiras posições são comuns |
| **GC content** | Distribuição do percentual de G+C | Curva próxima de um "sino" compatível com o organismo; picos extras podem sugerir contaminação |
| **Adapter content** | Presença de sequências de adaptador | Quanto mais perto de 0%, melhor |
| **Overrepresented sequences** | Sequências que se repetem demais | Pode indicar adaptadores ou contaminantes |
| **Sequence duplication levels** | Quanto do dado é duplicado | Interpretação depende do tipo de experimento |

**Os "semáforos" do FastQC.** Cada módulo recebe ✅ (passou), ⚠️ (aviso) ou ❌ (falhou). São **regras gerais**, não sentenças: um ❌ pode ser normal para certos experimentos (por exemplo, RNA-seq costuma "falhar" em conteúdo de bases nas primeiras posições). **O que vale é a sua interpretação**, não a cor.

### MultiQC

**O que é.** Uma ferramenta que reúne os relatórios de **várias** amostras (por exemplo, 50 relatórios do FastQC) em **um único painel comparativo**.

**Por que importa.** Em biovigilância é comum processar dezenas ou centenas de amostras. Abrir um relatório de cada vez é inviável. O MultiQC mostra de uma vez qual amostra "destoa" das demais.

---

## 1.4 Limpeza dos dados: Trimmomatic e Cutadapt

### O problema

Depois do diagnóstico, é comum encontrar:

1. **Pontas de baixa qualidade** nas reads (a qualidade tende a cair ao longo da leitura);
2. **Adaptadores**: sequências artificiais, adicionadas no preparo da amostra, que não fazem parte do organismo e podem aparecer no fim da read quando o fragmento é curto.

### A solução

**Trimmomatic** e **Cutadapt** funcionam como uma **tesoura de aparar bordas**: cortam as partes ruins e removem os adaptadores. Analogia: aparar as bordas tremidas de uma foto e ficar só com a parte nítida.

Operações principais do Trimmomatic (nomes em inglês, como aparecem na ferramenta):

| Operação | O que faz | Valor de partida comum |
|---|---|---|
| **ILLUMINACLIP** | Remove sequências de adaptador Illumina | Usar o arquivo de adaptadores correspondente ao kit |
| **SLIDINGWINDOW** | Percorre a read em "janelas" e corta quando a qualidade média da janela cai abaixo de um limite | Janela de 4 bases, qualidade média mínima 20 |
| **LEADING / TRAILING** | Remove bases de baixa qualidade no início/fim | Q mínimo de 3 |
| **MINLEN** | Descarta reads que ficaram curtas demais | 36 bases |

> Os valores acima são pontos de partida usuais na literatura. Na aula, o professor indicará os valores a usar.

### Saídas em *paired-end*

Ao aparar reads pareadas, pode acontecer de **uma** read do par ser descartada e a outra sobreviver. Por isso o Trimmomatic gera **quatro** arquivos:

| Arquivo | Conteúdo | Usar na etapa seguinte? |
|---|---|---|
| R1 paired | R1 cujo par (R2) também sobreviveu | ✅ Sim |
| R2 paired | R2 cujo par (R1) também sobreviveu | ✅ Sim |
| R1 unpaired | R1 "órfãs" | ❌ Não, neste exercício |
| R2 unpaired | R2 "órfãs" | ❌ Não, neste exercício |

**Boa prática:** rodar o FastQC **de novo** depois da limpeza e comparar *antes × depois*.

---

## 1.5 Alinhamento contra um genoma de referência

### A ideia

Depois de limpas, as reads precisam ser "encaixadas" de volta em um **genoma de referência**, ou seja, um genoma já montado e conhecido. Analogia: resolver um quebra-cabeça gigante comparando cada peça com a imagem da caixa.

**Importante, para não confundir:**

| | O que é | Exemplo | Formato |
|---|---|---|---|
| **Amostra (run)** | Os dados brutos que **você quer analisar**: milhões de fragmentos | `SRR957824` | FASTQ |
| **Genoma de referência** | O "livro completo" já montado, usado como gabarito | `NC_002695.1` | FASTA |

Uma conclusão importante: **raramente temos a referência exata da amostra**. Normalmente alinhamos contra o parente próximo mais bem documentado. O alinhamento revela **onde a amostra é igual e onde diverge** da referência (mutações, inserções, deleções), e muitas vezes é essa diferença que interessa.

### BWA e Bowtie2

São dois **alinhadores**: recebem reads e uma referência e descobrem a posição de origem de cada read.

| Critério | BWA (Burrows-Wheeler Aligner) | Bowtie2 |
|---|---|---|
| Dado mais adequado | Reads mais longas, indels maiores | Reads curtas de Illumina |
| Uso típico | Genômica humana/clínica | RNA-seq e metagenômica |
| Nesta aula | Citado como alternativa | **Usado na prática** |

Para este tipo de dado os dois são frequentemente intercambiáveis. A escolha costuma refletir o *pipeline* já estabelecido no laboratório.

### O formato SAM/BAM: o "resultado" do alinhamento

O alinhador devolve um arquivo que diz, **read por read**, onde ela caiu. Esse formato tem duas versões do mesmo conteúdo:

- **SAM** (*Sequence Alignment/Map*): texto legível por humanos;
- **BAM**: versão binária e comprimida do SAM (menor e mais rápida). É o que você terá no Galaxy.

**Estrutura de um SAM.** Linhas que começam com `@` são **cabeçalho** (por exemplo, `@SQ` lista as sequências de referência e seus tamanhos). As demais são **alinhamentos**: uma linha por read, com colunas separadas por tabulação.

Exemplo clássico (do artigo que definiu o formato):

```
@SQ SN:ref LN:45
r001  163  ref   7  30  8M2I4M1D3M  =  37   39  TTAGATAAAGGATACTA  *
r002    0  ref   9  30  3S6M1P1I4M   *   0    0  AAAAGATAAGGATA     *
r003    0  ref   9  30  5H6M         *   0    0  AGCTAA             *   NM:i:1
r004    0  ref  16  30  6M14N5M      *   0    0  ATAGCTTCAGC        *
r003   16  ref  29  30  6H5M         *   0    0  TAGGC              *   NM:i:0
r001   83  ref  37  30  9M          =   7  -39  CAGCGCCAT          *
```

As **11 colunas obrigatórias** (tomando a 1ª linha, `r001`, como exemplo):

| # | Campo | Valor | Significado |
|---|---|---|---|
| 1 | QNAME | `r001` | Nome da read (o mesmo da linha `@` do FASTQ) |
| 2 | FLAG | `163` | Código numérico que empacota várias informações "sim/não" (veja abaixo) |
| 3 | RNAME | `ref` | Em qual sequência de referência (cromossomo) a read caiu |
| 4 | POS | `7` | Posição inicial do alinhamento (começa em 1) |
| 5 | MAPQ | `30` | Qualidade do **mapeamento**: confiança de que a posição está certa (mesma lógica logarítmica do Phred) |
| 6 | CIGAR | `8M2I4M1D3M` | "Receita" de como a read se encaixou |
| 7 | RNEXT | `=` | Referência onde a parceira (mate) caiu; `=` significa "a mesma" |
| 8 | PNEXT | `37` | Posição da parceira |
| 9 | TLEN | `39` | Tamanho estimado do fragmento |
| 10 | SEQ | `TTAGATAAAGGATACTA` | A sequência da read |
| 11 | QUAL | `*` | Qualidade por base (`*` = não informada neste exemplo; normalmente é a mesma string do FASTQ) |

Depois da 11ª coluna podem vir **tags opcionais**, como `NM:i:1` (número de diferenças entre a read e a referência).

**Decifrando o CIGAR.** É uma sequência de pares "número + letra":

| Letra | Significado | Em palavras |
|---|---|---|
| **M** | Match/mismatch | A read está "alinhada" nessas posições (as bases podem ou não ser idênticas) |
| **I** | Insertion | A read tem bases a mais que a referência não tem |
| **D** | Deletion | A referência tem bases que a read não tem |
| **N** | Skipped region | Região "pulada" na referência (em RNA-seq, tipicamente um íntron; revisitaremos na Unidade 4) |
| **S** | Soft clip | Bases da ponta da read que **não** foram alinhadas, mas continuam no arquivo |
| **H** | Hard clip | Bases da ponta da read que não foram alinhadas e foram **removidas** do arquivo |

Exemplo: `8M2I4M1D3M` = 8 bases alinhadas, 2 bases inseridas, 4 alinhadas, 1 base deletada, 3 alinhadas. Conferência: 8 + 2 + 4 + 3 = 17, que é o comprimento de `TTAGATAAAGGATACTA`. ✅

**Decifrando o FLAG.** O número é a **soma** de "bandeiras" (cada uma vale uma potência de 2):

| Valor | Significado |
|---|---|
| 1 | A read faz parte de um par |
| 2 | O par foi alinhado "corretamente" (distância e orientação esperadas) |
| 4 | A read **não** alinhou |
| 8 | A parceira não alinhou |
| 16 | A read está na fita reversa |
| 32 | A parceira está na fita reversa |
| 64 | É a primeira do par (R1) |
| 128 | É a segunda do par (R2) |

Exemplos: `163 = 128 + 32 + 2 + 1` → é a segunda do par, a parceira está na fita reversa, o par está correto e a read é pareada. `83 = 64 + 16 + 2 + 1` → primeira do par, na fita reversa, par correto, pareada. `0` → read sem par, alinhada na fita direta.

> Você não precisa fazer essa conta à mão. É o que o **SAMtools** faz por você (ver a seguir). Entender o conceito ajuda a interpretar o que as estatísticas significam.

### SAMtools

O "canivete suíço" dos arquivos SAM/BAM: ordena, filtra, indexa, converte e **resume**. Na aula usaremos o **`samtools flagstat`**, que lê o FLAG de todas as reads e devolve um resumo estatístico (quantas alinharam, quantas formaram pares corretos etc.).

---

## 1.6 Mudando de escala: busca e comparação de sequências

Até aqui: milhões de reads contra **um** genoma. Agora, outra pergunta:

> Tenho **uma** sequência (um gene, uma proteína). **Que sequências parecidas já existem no mundo?**

### BLAST

**BLAST** (*Basic Local Alignment Search Tool*) é a ferramenta mais usada em bioinformática para comparar uma sequência contra bancos mundiais (GenBank, RefSeq, UniProt) e devolver os "parentes" mais parecidos.

| Variante | Compara | Quando usar |
|---|---|---|
| `blastn` | DNA contra DNA | Gene, fragmento de genoma |
| `blastp` | Proteína contra proteína | Sequência de aminoácidos |
| `blastx` | DNA traduzido nos 6 quadros de leitura contra proteínas | DNA de função desconhecida |

**Como ler um resultado do BLAST:**

| Coluna | O que significa | Como interpretar |
|---|---|---|
| **Query cover** | % da **sua** sequência que foi alinhada | Alto = o hit cobre quase toda a sua sequência |
| **Percent identity** | % de posições idênticas no trecho alinhado | Alto = muito parecido; mas veja sempre junto com a cobertura |
| **E-value** | Quantos acertos desse nível se esperaria por **puro acaso** num banco desse tamanho | **Quanto menor, mais significativo.** Valores muito próximos de 0 (ex.: `1e-50`) indicam similaridade real |
| **Bit score** | Pontuação do alinhamento, independente do tamanho do banco | Maior = melhor |

> **Cuidado clássico:** 100% de identidade em 20 aminoácidos de uma sequência de 500 **não** é um bom hit. Olhe sempre identidade **e** cobertura **e** E-value juntos.

### Clustal Omega e MUSCLE

O BLAST compara uma sequência contra um banco gigante. Estas ferramentas fazem **alinhamento múltiplo**: alinham **várias sequências ao mesmo tempo**, lado a lado, mostrando regiões **conservadas** (iguais) e **variáveis**. Útil, por exemplo, para comparar a mesma proteína em espécies diferentes.

Legenda de símbolos na linha abaixo do alinhamento do Clustal:

| Símbolo | Significado |
|---|---|
| `*` | Posição idêntica em todas as sequências |
| `:` | Substituições muito conservativas (aminoácidos de propriedades parecidas) |
| `.` | Substituições pouco conservativas |
| (vazio) | Sem conservação |

**Por que isso importa:** regiões muito conservadas ao longo da evolução costumam ser **funcionalmente importantes** (sítios ativos de enzimas, por exemplo). Esse raciocínio volta na Unidade 8 (filogenética) e na Unidade 7 (estrutural).

---

# PARTE 2 — PRÁTICA GUIADA NO GALAXY

## Visão geral

O **Galaxy** (https://usegalaxy.org) é uma plataforma web gratuita de bioinformática. Você escolhe a ferramenta em um menu, preenche um formulário e clica em **Run Tool/Execute**. Os resultados aparecem no painel **History** (à direita). Cada resultado é um "dataset" numerado.

**Cores do histórico:**

| Cor do dataset | Significa |
|---|---|
| Cinza | Na fila de espera |
| Amarelo | Em execução |
| Verde | Concluído ✅ |
| Vermelho | Erro ❌ (veja a Parte 4) |

> O Galaxy público é compartilhado. Em horários de pico, a fila pode demorar alguns minutos. **Isso é normal.** Aproveite a espera para ler o passo seguinte.

## Amostras utilizadas

| | Procarioto | Eucarioto |
|---|---|---|
| **Organismo** | *Escherichia coli* O157 | Levedura (*Saccharomyces cerevisiae*) |
| **Accession SRA** | `SRR957824` | `SRR453566` |
| **Genoma de referência** | *E. coli* O157:H7 str. Sakai, `NC_002695.1` (RefSeq) | Assembly R64-1-1 (Ensembl Fungi) |
| **Papel na aula** | Amostra **principal** da prática | Amostra opcional, para comparar com um eucarioto |

Pontos de atenção sobre as amostras:

- **Os dados do SRA são brutos.** Chegam exatamente como saíram do sequenciador, sem limpeza. É por isso que o exercício FastQC → Trimmomatic faz sentido.
- **`SRR957824` pertence a um projeto de investigação de surto de *E. coli* O157** (BioProject PRJNA215830). A referência Sakai **não é** o genoma exato dessa amostra: é um parente próximo, de alta qualidade, usado como gabarito padrão do sorotipo O157:H7.
- **Procedência exata da amostra.** Não afirme em trabalhos que ela veio "de um paciente" sem conferir. O campo *isolation source* do **BioSample** associado ao run informa se a origem foi clínica, alimentar ou ambiental. Verificar metadados faz parte do trabalho do bioinformata.
- **Para a levedura:** antes de usar `SRR453566`, confira no SRA o tipo de biblioteca (RNA-seq ou DNA) e se é *single-end* ou *paired-end*, pois isso muda as opções do Passo 5.
- **Tamanho dos dados.** O professor poderá fornecer uma versão reduzida (subconjunto de reads) para caber no tempo de aula.

---

## Passo 0 — Preparação (10 min)

1. Entre em **https://usegalaxy.org** e faça login.
2. No painel **History** (à direita), clique no nome "Unnamed history" e renomeie para `Aula2_SeuNome`.

**Você deve ver:** o título do histórico atualizado e a lista vazia.

---

## Passo 1 — Baixar os dados do SRA (15 min)

O **SRA** (*Sequence Read Archive*) é o repositório mundial de dados brutos de sequenciamento, mantido pelo NCBI. No Galaxy, usamos a ferramenta **fastq-dump** para baixar um *run* e convertê-lo em FASTQ.

1. Abra a ferramenta **"Faster Download and Extract Reads in FASTQ" (fastq-dump)** (link na Seção 6).
2. No campo de *accession*, digite `SRR957824`.
3. Clique em **Run Tool**.

**Você deve ver:** um item no histórico, que fica verde ao terminar, com nome parecido com **"Paired-end data (fastq-dump)"**. Ele é uma **coleção pareada** que contém R1 e R2 juntos. Clique nele e depois nos arquivos internos para espiar o conteúdo.

**Exercício rápido:** clique no "olho" 👁 de um dos arquivos. Identifique as 4 linhas de uma read e o que significam os caracteres da linha 4.

---

## Passo 2 — Diagnóstico com FastQC (15 min)

1. Abra a ferramenta **FastQC** (link na Seção 6).
2. Em *Short read data from your current history*, selecione os arquivos FASTQ brutos (R1 e R2).
3. Clique em **Run Tool**.
4. Ao terminar, abra o resultado **"Web page"** (👁).

**Perguntas-guia (anote as respostas, vão para o seu relatório):**

1. Quantas reads há na amostra (*Total Sequences*)? Qual o comprimento?
2. Em *Per base sequence quality*: a qualidade cai nas pontas? A partir de qual posição?
3. Em *Adapter Content*: há sinal de adaptadores? Em qual percentual?
4. Em *GC content*: a curva parece um "sino" único? (Para a *E. coli*, o %GC esperado é de cerca de 50%.)
5. Algum módulo ficou ⚠️ ou ❌? Você considera que isso é preocupante? Por quê?

📸 **Tire um print do gráfico *Per base sequence quality* (ANTES da limpeza).** Você vai precisar dele.

---

## Passo 3 — Limpeza com Trimmomatic (15 min)

1. Abra a ferramenta **Trimmomatic** (link na Seção 6).
2. Em *Single-end or paired-end reads*, escolha **Paired-end** e selecione R1 e R2.
3. Configure as operações (o professor confirmará os valores):
   - **ILLUMINACLIP** (remoção de adaptadores Illumina): ativada.
   - **SLIDINGWINDOW**: janela de 4 bases, qualidade média mínima 20.
   - **MINLEN**: 36.
4. Clique em **Run Tool**.

**Você deve ver:** quatro novos datasets:

- `Trimmomatic on ... (R1 paired)` ✅ usar
- `Trimmomatic on ... (R2 paired)` ✅ usar
- `... (R1 unpaired)` ❌ não usar agora
- `... (R2 unpaired)` ❌ não usar agora

> Atenção: confirme que o tipo dos arquivos é **fastqsanger** ou **fastqsanger.gz**. O sufixo `.gz` apenas indica compressão; o Galaxy descomprime sozinho. É a forma recomendada de guardar dados de sequenciamento.

---

## Passo 4 — Comparar antes × depois (10 min)

1. Rode o **FastQC novamente**, agora sobre os arquivos **R1 paired** e **R2 paired**.
2. Compare com o relatório do Passo 2.

**Perguntas-guia:**

1. A qualidade nas pontas melhorou?
2. O sinal de adaptadores diminuiu ou sumiu?
3. Quantas reads sobraram? Perdemos muitas? Isso é aceitável?

📸 **Print do gráfico *Per base sequence quality* (DEPOIS da limpeza)**, para comparar com o anterior.

---

## Passo 5 — Alinhamento com Bowtie2 (20 min)

### 5.1 Preparar os arquivos limpos como coleção pareada

O Bowtie2 espera os dois arquivos do par amarrados em uma **coleção pareada** (*paired collection*). Isso evita que R1 e R2 sejam trocados ou desalinhados por engano. Os arquivos limpos do Trimmomatic estão soltos, então é preciso montar a coleção.

1. No painel do histórico, clique no ícone de **seleção múltipla** (caixinhas, no topo do histórico).
2. Marque **somente** os dois datasets **R1 paired** e **R2 paired** do Trimmomatic.
3. No menu *"n of m selected"*, escolha **Build List of Pairs** (construir lista de pares).
4. Na tela **Auto Pairing**, preencha os dois filtros de nome:
   - campo da esquerda (forward): `forward`
   - campo da direita (reverse): `reverse`

   Você deve ver *"Auto-matched 1 pair(s)"*. Se aparecer *0*, tente `R1`/`R2` ou `_1`/`_2`, conforme o nome dos seus arquivos.
5. **Se o pareamento automático não funcionar**, avance para a etapa **Builder** e pareie manualmente: clique no ícone de **corrente** (🔗) de cada um dos dois datasets "UNPAIRED" para uni-los em um par. Confirme que R1 está como *forward* e R2 como *reverse*.
6. Dê um nome à coleção (ex.: `SRR957824_trimmed_paired`) e clique em **Create list**.

**Você deve ver:** uma nova coleção no histórico, contendo um par.

### 5.2 Obter o genoma de referência

Opção A (recomendada): verifique se o genoma já está disponível como **índice embutido** (*built-in index*) na lista da ferramenta Bowtie2. Se estiver, pule para 5.3.

Opção B: envie o FASTA da referência.

1. Menu **Upload** → aba **Paste/Fetch data**.
2. Cole esta URL (acesso direto ao FASTA do cromossomo `NC_002695.1` no NCBI):
   `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002695.1&rettype=fasta&retmode=text`
3. Em *Type*, escolha **fasta**. Clique em **Start**.

**Você deve ver:** um dataset FASTA com **uma** sequência (o cromossomo da *E. coli* Sakai, de ~5,5 milhões de bases).

> Os plasmídeos da cepa Sakai não estão incluídos neste arquivo. Para este exercício introdutório isso não é problema, mas explica por que uma pequena fração das reads pode não alinhar.

### 5.3 Rodar o Bowtie2

1. Abra a ferramenta **Bowtie2** (link na Seção 6).
2. *Is this single or paired library*: **Paired-end**.
3. *FASTQ Paired Dataset*: selecione a **coleção** criada em 5.1 (`SRR957824_trimmed_paired`).
4. *Will you select a reference genome from your history or use a built-in index?*
   - se usar índice embutido: escolha a opção *built-in* e selecione o genoma;
   - se enviou o FASTA: escolha *use a genome from history* e selecione-o.
5. Deixe os **demais parâmetros no padrão**.
6. Clique em **Run Tool**.

**Você deve ver:** um dataset no formato **BAM** (as reads já posicionadas no genoma) e, em geral, um resumo de alinhamento da própria ferramenta.

> **"Não aparece o meu arquivo limpo na lista."** O campo exige uma **coleção pareada**, não dois arquivos soltos. Volte ao passo 5.1. Não é necessário compactar nem zipar nada.

**Reflexão:** se um par R1/R2 fosse trocado ou desalinhado, o que aconteceria com o alinhamento? (Dica: pense na distância esperada entre as duas pontas do fragmento.)

---

## Passo 6 — Resumo estatístico com SAMtools flagstat (10 min)

1. Abra a ferramenta **Samtools flagstat** (link na Seção 6).
2. Em *BAM file*, selecione o BAM gerado pelo Bowtie2.
3. Clique em **Run Tool** e abra o resultado.

### Exemplo de saída e como interpretá-la

Resultado obtido numa execução de referência desta aula (os seus números podem diferir um pouco, conforme o subconjunto de dados e os parâmetros):

```
2895102 + 0 in total (QC-passed reads + QC-failed reads)
2895102 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2765970 + 0 mapped (95.54% : N/A)
2765970 + 0 primary mapped (95.54% : N/A)
2895102 + 0 paired in sequencing
1447551 + 0 read1
1447551 + 0 read2
2632294 + 0 properly paired (90.92% : N/A)
2755020 + 0 with itself and mate mapped
10950 + 0 singletons (0.38% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

| Linha | Significado | Leitura neste exemplo |
|---|---|---|
| `in total` | Total de reads processadas | ~2,9 milhões. O "+ 0" são reads reprovadas no controle de qualidade do sequenciador |
| `primary` | Alinhamentos "principais" (o melhor para cada read) | Igual ao total |
| `secondary` / `supplementary` | Alinhamentos alternativos / parciais | Zero: dado limpo, sem ambiguidade registrada |
| `duplicates` | Reads marcadas como duplicatas de PCR | Zero **porque não rodamos a etapa de marcação**, não porque não existam duplicatas |
| **`mapped`** | Reads que alinharam em algum lugar | **95,54%** ✅ |
| `paired in sequencing` | Reads de sequenciamento pareado | Todas |
| `read1` / `read2` | Quantas são R1 e quantas são R2 | Iguais (1.447.551 cada): pareamento íntegro ✅ |
| **`properly paired`** | Pares alinhados com distância e orientação esperadas | **90,92%** ✅ |
| `with itself and mate mapped` | Read e parceira, ambas alinhadas | Quase todas |
| `singletons` | Read alinhou, parceira não | 0,38% (muito baixo) ✅ |
| `with mate mapped to a different chr` | Parceira em outra sequência | Zero: sem sinal de rearranjos estranhos |

### Discussão

**Por que não deu 100%?** Alguns motivos plausíveis:

1. A amostra é um **parente próximo**, e não idêntica à referência: existem regiões próprias da cepa que a referência não tem.
2. Os **plasmídeos** da amostra não estão na referência usada.
3. Resquício de reads de baixa qualidade que passaram pela limpeza.

**O que 95% de alinhamento sugere?** Que amostra e referência são geneticamente muito próximas, como esperado para duas cepas do mesmo sorotipo. Uma taxa muito baixa levantaria a hipótese de **referência errada** ou **contaminação**. Isso é um diagnóstico valioso.

---

## Passo 7 — Busca com BLAST (15 min)

1. Acesse **https://blast.ncbi.nlm.nih.gov** e clique em **Protein BLAST (blastp)**.
2. No campo de sequência, cole a sequência da **sua** proteína de interesse em formato FASTA (comece com `>nome`).
   *Se você não tiver uma ainda, use como exemplo a lisozima de galinha (UniProt **P00698**): baixe o FASTA na página dela no UniProt e cole.*
3. Em *Database*, deixe **Standard databases (nr)** (ou *Non-redundant protein sequences*).
4. Clique em **BLAST** e aguarde (pode levar de segundos a minutos).

**Você deve ver:** uma tabela de *Descriptions* com as sequências mais parecidas, ordenadas por similaridade.

**Perguntas-guia:**

1. Qual foi o primeiro hit **que não é a sua própria sequência**? De qual organismo?
2. Qual a *Percent identity*, a *Query cover* e o *E-value* dele?
3. O hit é convincente? (Lembre: identidade, cobertura **e** E-value juntos.)
4. Em quantos organismos diferentes essa proteína aparece? O que isso sugere sobre a conservação evolutiva dela?

📸 **Print da tabela de resultados**, com os hits à vista.

---

## Passo 8 — Alinhamento múltiplo com Clustal Omega (10 min)

1. No resultado do BLAST, marque a caixa de **3 a 4 hits** de organismos diferentes.
2. Clique em **Download → FASTA (complete sequence)** para obter as sequências.
3. Acesse **https://www.ebi.ac.uk/jdispatcher/msa/clustalo**, cole as sequências (incluindo a sua) e clique em **Submit**.
4. Abra a aba **Alignments**.

**Perguntas-guia:**

1. Existem trechos longos com `*` (idênticos em todas)? Onde? O que isso poderia significar funcionalmente?
2. Existem trechos muito variáveis? Em quais sequências?
3. As espécies mais "próximas" evolutivamente parecem ter mais trechos em comum?

📸 **Print do alinhamento**, mostrando uma região conservada.

---

## Passo 9 (complementar/opcional) — Sequenciamento Sanger: arquivos `.ab1`

**Contexto.** Antes do NGS, o padrão era o **sequenciamento Sanger**, que produz **uma leitura longa por vez** (até ~900 a 1.000 bases) em um arquivo **`.ab1`**. Esse arquivo guarda o **cromatograma**: os picos coloridos de fluorescência, um por base. O conceito de Phred Score **nasceu aqui**, no programa Phred, que analisava esses cromatogramas, e depois foi reaproveitado no NGS.

1. **Obter os dados.** Acesse **https://zenodo.org/records/7104640** e baixe um par de arquivos `.ab1` (leitura *forward* e *reverse* do mesmo fragmento) listados no registro.
2. **Ver o cromatograma e o Phred por base.** Acesse **https://conductscience.com/tools/file-formats/ab1-chromatogram#tool-calculator**, envie um dos `.ab1` e observe:
   - os picos coloridos (um pico bem definido = base confiável);
   - o Phred por posição;
   - as extremidades da leitura (início e fim costumam ter qualidade pior). Use a função de *trimming* da ferramenta para ver onde ela cortaria e exporte a sequência em FASTA.
3. **Montar o contig.** Repita o passo 2 para a outra leitura. Depois, acesse **https://doua.prabi.fr/cgi-bin/run_cap3**, envie as duas sequências FASTA e execute o **CAP3**. Ele procura a região de sobreposição entre as leituras (testando também a orientação reversa) e gera a sequência **consenso** (*contig*).

**Perguntas-guia:**

1. Em que posições a qualidade do cromatograma começa a cair?
2. A sobreposição entre as duas leituras foi grande ou pequena? O consenso ficou mais longo do que cada leitura isolada?
3. Em posições onde as duas leituras discordam, qual delas tinha melhor qualidade? (Esse é o raciocínio que o CAP3 usa quando recebe também a qualidade.)

**Conexão com o que vimos:** montar um contig de 2 leituras Sanger é, em miniatura, a lógica de "montagem" que veremos em larga escala na Unidade 3 (montagem de genomas).

---

# PARTE 3 — ATIVIDADE ASSÍNCRONA (avaliativa)

## O que entregar

Um **resumo de 1 página** (PDF ou Word) com prints de tela e interpretação **em texto próprio**.

### Roteiro

1. **Dados e FastQC.** Escolha um dataset pequeno do SRA (pode ser o da aula, `SRR957824`, ou outro). Rode o FastQC pelo Galaxy.
   Escreva **3 a 4 linhas** interpretando o relatório: qualidade geral, onde a qualidade cai, presença de adaptadores.
2. **BLAST da sua pesquisa.** Escolha uma sequência (proteína ou gene) **ligada à sua linha de pesquisa** (química de proteínas/macromoléculas) e rode a busca BLAST.
3. **Registre:**
   - o parente mais próximo encontrado (outro que não a própria sequência);
   - identidade (%), cobertura (%) e **E-value**;
   - **uma frase** interpretando biologicamente o que isso significa.
4. Anexe os **prints** de cada etapa.

### Modelo de estrutura

| Seção | Conteúdo | Tamanho sugerido |
|---|---|---|
| Identificação | Nome, data, accessions usados | 2 linhas |
| 1. Qualidade dos dados | Print do FastQC + interpretação | ½ página |
| 2. Busca de similaridade | Print do BLAST + tabela com identidade, cobertura, E-value | ½ página |
| 3. Interpretação | O que os resultados significam para a sua pesquisa | 3 a 5 linhas |

### Critério de avaliação

O que será avaliado é a **clareza da interpretação**, não a perfeição técnica. O objetivo é demonstrar que você entendeu o **significado** de cada resultado, e não apenas que "rodou o programa".

### Checklist final

- [ ] Citei os *accessions* utilizados (SRA e/ou RefSeq/UniProt).
- [ ] Interpretei com minhas palavras, em vez de colar saídas brutas.
- [ ] Os prints estão legíveis e legendados.
- [ ] Relacionei o resultado à minha linha de pesquisa.
- [ ] Mencionei limitações (por exemplo: "a referência é um parente próximo, não a cepa exata").

---

# PARTE 4 — SOLUÇÃO DE PROBLEMAS

| Sintoma | Causa provável | O que fazer |
|---|---|---|
| Dataset **vermelho** | Entrada errada ou erro na ferramenta | Clique no dataset → ícone de **"bug"/detalhes** e leia a mensagem. Confira se o tipo de dado (*datatype*) está correto |
| Dataset **cinza** por muito tempo | Fila do Galaxy público | Aguarde. Não relance a mesma execução várias vezes, isso piora a fila |
| Não aparece o arquivo limpo no Bowtie2 | O campo exige **coleção pareada** | Monte a coleção (Passo 5.1) |
| *Auto-matched 0 pair(s)* | Os filtros não reconhecem o nome dos arquivos | Tente `forward`/`reverse`, `R1`/`R2`, `_1`/`_2`; se não funcionar, pareie manualmente na etapa *Builder* |
| Dúvida sobre `.gz` | Arquivo comprimido | Pode usar normalmente. O Galaxy descomprime sozinho; `fastqsanger.gz` é aceito |
| Taxa de alinhamento muito baixa | Referência errada, contaminação ou amostra muito divergente | Confira o organismo da amostra e a referência; rode o FastQC de novo e olhe o %GC |
| Resultados diferentes dos do colega | Subconjuntos de dados ou parâmetros diferentes | Compare os parâmetros usados; pequenas diferenças são normais |
| BLAST demorando | Servidor do NCBI ocupado | Aguarde e não feche a aba; o resultado tem um *RID* (identificador) para recuperar depois |
| Perdi meu progresso | Histórico não salvo | O Galaxy salva automaticamente ao estar logado. Verifique em *User → Histories* |

---

# PARTE 5 — GLOSSÁRIO

| Termo | Explicação simples |
|---|---|
| **Adaptador** | Sequência artificial adicionada no preparo da amostra; não pertence ao organismo |
| **Alinhamento** | Posicionar reads (ou sequências) lado a lado para ver o que coincide |
| **ASCII** | Tabela que associa cada caractere a um número |
| **BAM / SAM** | Formatos de arquivo de alinhamento (BAM = versão binária do SAM) |
| **Cobertura (*coverage*)** | Quantas reads cobrem cada posição do genoma |
| **Contig** | Sequência contínua obtida pela sobreposição de leituras |
| **CIGAR** | "Receita" que descreve como uma read se encaixou na referência |
| **E-value** | Número esperado de acertos por acaso; menor = mais significativo |
| **FASTA** | Formato de texto para sequências (`>nome` seguido da sequência) |
| **FASTQ** | Formato de texto para reads com qualidade (4 linhas por read) |
| **FLAG** | Número do SAM que codifica propriedades da read (pareada, reversa etc.) |
| **Genoma de referência** | Genoma já montado, usado como gabarito |
| **Homólogas** | Sequências com ancestral comum, portanto parecidas |
| **Mapeamento** | Sinônimo de alinhar reads a uma referência |
| **MAPQ** | Confiança do mapeamento de uma read |
| **NGS** | *Next Generation Sequencing*: sequenciamento em larga escala |
| **Paired-end** | Sequenciamento das duas pontas de cada fragmento |
| **Phred Score** | Nota de confiança por base, em escala logarítmica |
| **Read** | Fragmento de sequência lido pelo sequenciador |
| **Run** | Uma corrida de sequenciamento depositada no SRA (começa com SRR) |
| **Singleton** | Read alinhada cuja parceira não alinhou |
| **SRA** | Repositório mundial de dados brutos de sequenciamento |
| **Trimming** | Corte das partes ruins das reads |

---

# PARTE 6 — LINKS

| Recurso | Endereço |
|---|---|
| Galaxy (plataforma) | https://usegalaxy.org |
| SRA (NCBI) | https://www.ncbi.nlm.nih.gov/sra |
| fastq-dump (baixar SRA no Galaxy) | https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fsra_tools%2Ffastq_dump%2F3.1.1%2Bgalaxy1&version=latest |
| FastQC | https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fdevteam%2Ffastqc%2Ffastqc%2F0.74%2Bgalaxy1&version=latest |
| Trimmomatic | https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fpjbriggs%2Ftrimmomatic%2Ftrimmomatic%2F0.39%2Bgalaxy2&version=latest |
| Bowtie2 | https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fdevteam%2Fbowtie2%2Fbowtie2%2F2.5.5%2Bgalaxy0&version=latest |
| SAMtools flagstat | https://usegalaxy.org/?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fdevteam%2Fsamtools_flagstat%2Fsamtools_flagstat%2F2.0.8&version=latest |
| Genoma de referência *E. coli* (FASTA direto) | https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_002695.1&rettype=fasta&retmode=text |
| BLAST (NCBI) | https://blast.ncbi.nlm.nih.gov |
| Clustal Omega (EBI) | https://www.ebi.ac.uk/jdispatcher/msa/clustalo |
| Dados `.ab1` de exemplo (Zenodo) | https://zenodo.org/records/7104640 |
| Visualizador de `.ab1` com Phred (ConductScience) | https://conductscience.com/tools/file-formats/ab1-chromatogram#tool-calculator |
| CAP3 — montagem de contig (PRABI-Doua) | https://doua.prabi.fr/cgi-bin/run_cap3 |
| Tutorial Galaxy sobre `.ab1` (GTN) | https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/Manage_AB1_Sanger/tutorial.html |

---

## Bibliografia

- Pevsner, J. *Bioinformatics and Functional Genomics*. Wiley-Blackwell.
- Compeau, P.; Pevzner, P. *Bioinformatics Algorithms: An Active Learning Approach*.
- Lesk, A. M. *Introduction to Bioinformatics*. Oxford University Press.
- Mariano, D. *Python para Bioinformática: Fundamentos de Programação para Bioinformática e Biologia Computacional*. Novatec, 2025.
- Especificação do formato SAM/BAM (documentação oficial do SAMtools/hts-specs).
- Documentação oficial: NCBI BLAST, SRA, FastQC, Trimmomatic, Bowtie2, Galaxy Training Network.

---

*Dica final: bioinformática se aprende fazendo. Se algo não funcionar de primeira, leia a mensagem de erro, confira o tipo de arquivo e tente de novo. Isso faz parte do processo.*
