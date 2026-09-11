# Aula 3 — Genômica: Montagem e Anotação de Genomas

**Disciplina:** Introdução à Bioinformática — PPGED/IME
**Carga horária da aula:** 3 horas (teoria + prática guiada)
**Público-alvo:** Pós-graduandos em Engenharia de Defesa, sem experiência prévia em programação ou biologia molecular aprofundada
**Pré-requisito direto:** Aula 2 (qualidade de NGS, alinhamento, BLAST) — todos os conceitos de "reads" e ferramentas do Galaxy já foram apresentados

---

## Ideia-guia da aula

Continuamos a metáfora da Aula 2: se o genoma sequenciado é um **livro despedaçado em milhões de pedacinhos (reads)**, a Aula 2 nos ensinou a checar a qualidade desses pedaços e alinhá-los contra um livro de referência já pronto. Hoje damos um passo além:

> **E se não existir um "livro de referência" pronto?** Precisamos remontar o livro do zero, juntando os pedaços que se sobrepõem (**montagem**), e depois "traduzir" esse livro remontado, marcando onde estão os capítulos e o que cada um significa (**anotação**).

---

## Parte 1 — TEORIA (1h30)

### 1.1 Gancho inicial — por que genômica importa para Engenharia de Defesa (10 min)

- Biovigilância e biodefesa: identificar rapidamente um patógeno desconhecido a partir de amostra ambiental ou clínica.
- Identificação forense biológica e verificação de agentes biológicos (cepa natural x manipulada).
- Ponte com pesquisa em química de proteínas: antes de estudar a estrutura de uma proteína (Unidade 7), é preciso saber que gene a codifica e onde ele está no genoma.

### 1.2 Como genes e genomas se organizam: procariontes x eucariontes (25 min)

Esse bloco é o alicerce conceitual da aula — sem ele, "anotar um genoma" vira uma caixa-preta.

**Procariontes (bactérias e arqueias):**
- Genoma geralmente **circular**, sem núcleo delimitado (fica solto no citoplasma).
- Genes ficam **compactados**, quase sem espaço entre eles — pouquíssimo DNA "não-codificante".
- Frequentemente organizados em **óperons**: vários genes relacionados são transcritos juntos, como uma linha de produção controlada por um único interruptor.
- Consequência prática: genoma pequeno (poucos milhões de "letras") e **fácil de anotar automaticamente** — é por isso que nossa prática de hoje usa uma bactéria.

**Eucariontes (do fungo ao ser humano):**
- Genoma **linear**, organizado dentro do núcleo da célula, em cromossomos.
- Muito mais DNA não-codificante (em humanos, menos de 2% do genoma codifica proteínas diretamente).
- Genes **fragmentados** em **éxons** (partes que codificam a proteína) e **íntrons** (trechos removidos antes da tradução) — como um texto útil intercalado com parágrafos de rascunho que são cortados antes da versão final (processo chamado *splicing*).
- Consequência prática: genoma muito maior e **anotação muito mais complexa**, porque o programa precisa "adivinhar" onde começa e termina cada éxon.

> **Por que isso importa agora:** o Prokka (ferramenta que usaremos na prática) funciona rápido e bem justamente porque trabalha com procarionte — genoma pequeno, sem íntrons, genes colados uns nos outros. Isso não é acaso, é a razão da nossa escolha didática.

### 1.3 Anatomia de um gene: as "placas de sinalização" do genoma (20 min)

Antes de falar em anotação automática, é importante entender **o que exatamente** um programa como o Prokka está procurando quando "lê" um genoma. Um gene não é só a sequência que vira proteína — ele vem cercado de **regiões regulatórias**, como um texto que tem instruções de formatação antes e depois do conteúdo em si.

Analogia-guia: pense num gene como um **parágrafo de um manual técnico**, com marcações de onde começar a ler, onde parar, e (em eucariontes) trechos que precisam ser "cortados" antes da versão final.

- **Promotor**: região *antes* do gene onde a maquinaria de transcrição (RNA polimerase) se encaixa para começar a "copiar" o DNA em RNA — é o "aqui começa a leitura", o interruptor liga/desliga da expressão do gene.
- **Sítio de início da transcrição**: o ponto exato onde a cópia em RNA começa a ser produzida (geralmente logo após o promotor).
- **Códon de início (ATG)** e **códons de parada (stop codons)**: marcam, dentro do RNA, onde começa e termina a informação que efetivamente vira proteína — como o "abre aspas" e "fecha aspas" de uma citação dentro do texto.
- **Terminador**: região *depois* do gene que sinaliza à maquinaria de transcrição para parar de copiar — o "ponto final" do parágrafo.
- **Sítios de splicing (doador e aceptor)** — *só em eucariontes*: marcações exatas nas bordas entre éxons e íntrons, reconhecidas por um complexo molecular chamado *spliceossomo*, que "corta" os íntrons e "cola" os éxons entre si. É justamente prever essas bordas com precisão que torna a anotação eucariótica tão mais difícil (é o que ferramentas como AUGUSTUS e MAKER tentam fazer).
- **Regiões não traduzidas (UTRs, 5' e 3')**: trechos que são transcritos em RNA, mas não viram proteína — ficam nas pontas do gene, com papel regulatório (ex.: controlar a estabilidade ou a velocidade de tradução do RNA).

**No caso específico de procariontes** (relevante para a prática de hoje), vale destacar duas peças a mais:
- **Operador**: sequência próxima ao promotor onde proteínas regulatórias podem se ligar, funcionando como uma "trava" que libera ou bloqueia a transcrição.
- **Sequência Shine-Dalgarno**: um pequeno trecho no RNA, um pouco antes do códon de início, que sinaliza ao ribossomo onde encaixar para começar a "traduzir" a proteína — o equivalente bacteriano de "comece a ler aqui".

> **Conexão com a prática:** quando o Prokka anota o genoma bacteriano, ele está essencialmente procurando por esses sinais — sobretudo o par início/parada e padrões característicos de promotor/Shine-Dalgarno — para delimitar cada gene. Como não há íntrons em procariontes, ele não precisa lidar com sítios de splicing, o que simplifica (e acelera) muito o processo. Já em eucariontes, prever corretamente os sítios de splicing é o principal motivo de MAKER/AUGUSTUS serem tão mais complexos e demorados que o Prokka.

> **Nota complementar (material opcional, não é exercício obrigatório):** o Prokka identifica os genes, mas **não prevê promotor nem terminador** — isso fica a cargo de ferramentas separadas. Para quem quiser explorar por conta própria, ambas são gratuitas e 100% via navegador (sem instalar nada):
> - **BPROM** (previsão de promotor bacteriano): http://linux1.softberry.com/berry.phtml?topic=bprom&group=programs&subgroup=gfindb — cole a sequência da região até ~150 pb antes do gene de interesse; o resultado traz a posição prevista do promotor e as caixas -10/-35 reconhecidas pela RNA polimerase.
> - **ARNold** (previsão de terminador rho-independente): https://rna.igmors.u-psud.fr/toolbox/arnold/ — cole a sequência logo após o gene; o resultado traz a posição do terminador previsto e a energia livre (ΔG) da estrutura em grampo, quanto mais negativa, mais confiável.

### 1.4 Montagem de genomas (20 min)

- **Sequenciamento gera fragmentos, não o genoma inteiro** (retomada rápida da Aula 2).
- **Montagem (assembly)**: remontar o "livro rasgado" juntando pedaços que se sobrepõem (*overlap*), formando pedaços maiores chamados **contigs**, que por sua vez se agrupam em **scaffolds**.
- **De novo x referência**: montar do zero (sem mapa/genoma parecido) x usar um genoma parecido já conhecido como guia — hoje faremos *de novo*, o cenário mais realista quando se descobre um organismo novo ou desconhecido.
- **Qualidade da montagem — métrica N50**: é o tamanho de contig tal que metade do genoma está contida em pedaços iguais ou maiores que ele. Quanto maior o N50, melhor (menos fragmentado) está a montagem.

### 1.5 Anotação de genomas (15 min)

- **O que é**: depois de montado, é preciso "etiquetar" onde estão os genes e qual é a função provável de cada um — como anotar num mapa onde ficam as casas, ruas e praças.
- **Ferramenta em procariontes — Prokka**: recebe o genoma montado e devolve uma lista de genes encontrados, com nome e função provável, comparando com bancos de referência já conhecidos. Funciona como um "leitor de OCR" especializado em DNA bacteriano.
- **E em eucariontes?** O processo é bem mais complexo por causa dos íntrons/éxons. As ferramentas de referência são:
  - **MAKER / MAKER2**: equivalente conceitual ao Prokka, mas para eucariontes — combina várias fontes de evidência (similaridade com proteínas conhecidas, dados de RNA-seq, previsão estatística) para decidir onde estão os genes.
  - **AUGUSTUS**: prediz genes com base em modelos estatísticos treinados para reconhecer padrões de éxon/íntron específicos de cada espécie (tem versão web, sem precisar instalar nada).
  - **BRAKER**: combina o AUGUSTUS com dados reais de RNA-seq para melhorar a precisão — em vez de só "prever" onde estão os genes, usa evidência real de que aquele trecho é de fato expresso.
  - **Na prática do dia a dia**, poucas vezes alguém anota um genoma eucarioto do zero: para a maioria dos organismos já estudados (humano, camundongo, etc.), a anotação já existe pronta em bancos como **Ensembl** e **NCBI Genome** — usa-se a referência publicada.

---

## Intervalo (10 min)

---

## Parte 2 — PRÁTICA GUIADA (1h20)

Toda a prática síncrona continua na plataforma **Galaxy** (https://usegalaxy.org), sem sair do padrão de "buscar ferramenta no menu → configurar → rodar → ver resultado" já usado na Aula 2. SPAdes, Quast e Prokka estão nativamente disponíveis no Galaxy.

### Dataset sugerido

Um genoma bacteriano **fictício, mas realista**, de *Staphylococcus aureus* em miniatura, criado especificamente para ensino pela comunidade Galaxy Training Network (GTN): genoma pequeno, cobertura baixa (~19x), reads pareadas de 150 pb — roda em minutos, não em horas, e é o mesmo dataset usado nos tutoriais oficiais de montagem e anotação do GTN.

**Links para importar direto no Galaxy** (FASTQ, Illumina paired-end):
- `https://zenodo.org/record/582600/files/mutant_R1.fastq`
- `https://zenodo.org/record/582600/files/mutant_R2.fastq`

> **Nota de preparação do professor:** a montagem consome mais tempo de processamento do que FastQC/BWA. Recomenda-se rodar o SPAdes previamente com esse dataset e deixar o resultado pronto, evitando 15-20 min de espera ao vivo. Os alunos replicam o processo em paralelo ou acompanham a execução já em andamento.
>
> Esse é o mesmo dataset usado nos tutoriais oficiais do GTN (`training.galaxyproject.org`, buscar "Genome Assembly" e "Genome annotation with Prokka") — útil como material de apoio ou gabarito de comparação. Para uma versão com dados reais (não simulados), o GTN também tem um tutorial com sequenciamento real de MRSA em Illumina MiSeq, mas é mais pesado computacionalmente e não recomendado para uma prática síncrona de poucas horas.

### Passo 1 — Criar a history e importar as reads (10 min)
1. No Galaxy, clique no ícone "+" no painel de History e crie uma nova, chamada `Aula3_Turma_PPGED`.
2. Vá em **Upload Data → Paste/Fetch data** e cole as duas URLs (`mutant_R1.fastq` e `mutant_R2.fastq`).
3. Marque o tipo como `fastqsanger` se o Galaxy não detectar automaticamente.

### Passo 2 — Renomear os arquivos (5 min)
1. Por padrão o Galaxy usa a URL como nome do arquivo. Clique no lápis (editar atributos) de cada dataset e renomeie para algo claro, como `reads_R1` e `reads_R2`, para não se perder depois.

### Passo 3 — Montagem com SPAdes (25 min)
1. No menu de ferramentas do Galaxy, buscar por "SPAdes".
2. Em "Input reads", escolher o modo paired-end e selecionar `reads_R1` e `reads_R2`. Deixar os parâmetros de k-mer no padrão (o SPAdes escolhe automaticamente vários tamanhos).
3. Executar (ou acompanhar a execução pré-rodada pelo professor) — com esse dataset pequeno, leva poucos minutos.
4. Explicar, enquanto o job roda ou já com o resultado pronto: por trás do SPAdes há um método matemático (grafo de De Bruijn) que organiza a sobreposição dos fragmentos — não é necessário entender a matemática, só saber que ele é o "motor" que resolve o quebra-cabeça.
5. Abrir o arquivo de saída (`contigs.fasta`) e observar a lista de contigs gerados.

### Passo 4 — Avaliar a montagem com Quast (10 min)
1. Buscar a ferramenta "Quast" no Galaxy.
2. Selecionar o arquivo de contigs gerado pelo SPAdes como entrada e executar.
3. Interpretar juntos o relatório: número de contigs, tamanho total da montagem, e o **N50** — discutir se a montagem ficou "boa" ou "fragmentada".

### Passo 5 — Anotação com Prokka (15 min)
1. Buscar a ferramenta "Prokka" no Galaxy.
2. Em "contigs to annotate", selecionar o arquivo de contigs do SPAdes. Em "Genetic code", escolher **11** (Bacterial, Archaeal and Plant Plastid Code).
3. Executar. As saídas incluem arquivos `.gff`, `.gbk` e uma tabela de genes anotados.
4. Abrir a tabela de saída e, em turma, identificar 3 a 5 genes anotados (nome, posição, função prevista). Opcional: usar o JBrowse do Galaxy para visualizar as anotações sobre o contig mais longo.
5. Relacionar com a teoria: por que essa anotação foi rápida e direta — genoma procarionte, sem necessidade de prever splicing.

### Passo 6 — Comparação rápida com um genoma eucarioto já anotado (10 min)
1. Abrir, em outra aba do navegador, o **Ensembl** (ensembl.org) e buscar um gene humano ou de camundongo qualquer já anotado.
2. Mostrar visualmente a diferença de escala (tamanho do genoma, número de cromossomos) e a presença de éxons/íntrons na estrutura do gene.
3. Reforçar: essa anotação não foi feita "na hora" — é fruto de anos de trabalho colaborativo mundial (MAKER/AUGUSTUS/BRAKER + curadoria manual), diferente do que acabamos de fazer em minutos com uma bactéria.

### Fechamento (5 min)
- Recapitular o fluxo do dia: reads → montagem (SPAdes) → avaliação (Quast) → anotação (Prokka) → contraste com a complexidade eucariótica.
- Conectar com a próxima aula (Transcriptômica): "o que fazer quando, além do DNA, queremos saber quais genes estão sendo ativamente usados pela célula?"

---

## Atividade prática assíncrona (entre a Aula 3 e a Aula 4)

**Entrega:** resumo de 1 página (PDF ou Word), com prints de tela e interpretação em texto próprio.

**Roteiro para o aluno:**
1. Reaproveitar (ou gerar novamente) a montagem feita em aula no Galaxy com o SPAdes, a partir do dataset `mutant_R1.fastq` / `mutant_R2.fastq` (Zenodo, imaginary *S. aureus*).
2. Rodar o Quast sobre essa montagem e registrar: número de contigs, tamanho total, N50 — com uma frase de interpretação sobre a qualidade da montagem.
3. Rodar o Prokka sobre o genoma montado e listar 3 genes anotados, indicando nome e função prevista.
4. Pesquisar, no NCBI Genome ou Ensembl, o mesmo organismo (ou um próximo) já montado e anotado profissionalmente, e comparar brevemente (nº de genes, tamanho do genoma) com o resultado obtido em aula.

**Critério de avaliação:** clareza da interpretação (não a perfeição técnica), demonstrando que o aluno entendeu o *significado* de cada resultado — não apenas que "rodou o programa".

---

## Materiais e links de apoio

| Ferramenta / Recurso | Acesso | Necessita instalação? |
|---|---|---|
| Galaxy | https://usegalaxy.org | Não (web) |
| Dataset de exemplo (reads *S. aureus* miniatura) | https://zenodo.org/record/582600 | Não (web) |
| SPAdes | Via Galaxy | Não |
| Quast | Via Galaxy | Não |
| Prokka | Via Galaxy | Não |
| NCBI Genome | https://www.ncbi.nlm.nih.gov/genome | Não (web) |
| Ensembl | https://www.ensembl.org | Não (web) |
| AUGUSTUS (mencionado na teoria) | https://bioinf.uni-greifswald.de/augustus/ | Não (web) |
| BPROM (previsão de promotor — material opcional) | http://linux1.softberry.com/berry.phtml?topic=bprom&group=programs&subgroup=gfindb | Não (web) |
| ARNold (previsão de terminador — material opcional) | https://rna.igmors.u-psud.fr/toolbox/arnold/ | Não (web) |

## Bibliografia de apoio para esta aula

- Pevsner, J. *Bioinformatics and Functional Genomics*. Wiley-Blackwell.
- Compeau, P.; Pevzner, P. *Bioinformatics Algorithms: An Active Learning Approach*.
- Mariano, D. *Python para Bioinformática: Fundamentos de Programação para Bioinformática e Biologia Computacional*. Novatec Editora, 2025.
- Documentação oficial: SPAdes, Prokka, Quast, MAKER2, AUGUSTUS, NCBI Genome, Ensembl.
