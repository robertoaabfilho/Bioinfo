# Material de Apoio — Aula 3: Genômica, Montagem e Anotação de Genomas

Oct 8, 2026 · @Roberto

## Como usar este material

Este texto acompanha a **Aula 3 de Introdução à Bioinformática (PPGED/IME)**. Ele foi escrito para quem **nunca programou** e **nunca estudou biologia molecular a fundo**. Cada termo técnico é explicado quando aparece pela primeira vez, e cada ferramenta usada na prática é descrita em detalhe: o que faz, o que recebe, o que devolve e como interpretar o resultado.

### O que você vai aprender

Ao final da aula você será capaz de:

1. Explicar do que é feito o material genético (nucleotídeos, DNA, RNA) e como ele vira proteína.
2. Comparar como **genes** e **genomas** se organizam em **procariotos** (bactérias) e **eucariotos** (animais, plantas, fungos).
3. Reconhecer as partes de um gene: promotor, início, região codificadora, terminador e, nos eucariotos, íntrons e sítios de *splicing*.
4. Explicar o que é **montagem** de um genoma (reads, contigs, scaffolds) e o que é **anotação**.
5. Usar, no **Galaxy**, as ferramentas **SPAdes** (montagem), **Quast** (qualidade) e **Prokka** (anotação) em um genoma bacteriano real.
6. Dizer o que essas ferramentas **não** fazem, e quais ferramentas complementares existem (BPROM, ARNold, MAKER, AUGUSTUS).

### Cronograma da aula (3 horas)

| Bloco | Duração | O que acontece |
| --- | --- | --- |
| Teoria 1 | \~35 min | Biologia molecular: nucleotídeos, DNA, RNA |
| Teoria 2 | \~30 min | Organização de genomas e de genes: procariotos × eucariotos |
| Teoria 3 | \~25 min | Montagem, anotação e apresentação das ferramentas |
| **Intervalo** | 10 min |  |
| Prática | \~80 min | Galaxy: dados → SPAdes → Quast → Prokka → JBrowse |
| Fechamento | \~10 min | Discussão e atividade assíncrona |

### O que preparar antes da aula

- Um computador com **navegador** (Chrome ou Firefox) e internet. Não é preciso instalar nada.
- Uma **conta gratuita no Galaxy**: acesse *usegalaxy.org*, clique em *Login or Register* e crie a conta com seu e-mail. Confirme o e-mail **antes** da aula, pois sem isso o Galaxy limita o uso.
- Nenhum conhecimento prévio de programação. Em nenhum momento você escreverá código: tudo será feito clicando em menus e formulários.

### Como ler este material

Os termos importantes aparecem em **negrito** na primeira vez. Há um glossário ao final. As perguntas de revisão (com gabarito) servem para você testar se entendeu antes da próxima aula.

## Parte 1 — Biologia molecular: do nucleotídeo à proteína

Para entender um genoma é preciso saber de que ele é feito. A boa notícia: a ideia central cabe em poucas linhas. **Toda a informação genética de um ser vivo é um texto escrito com um alfabeto de apenas quatro letras.** O resto é como esse texto é guardado, copiado e lido.

### 1.1 O nucleotídeo: a "letra" do alfabeto genético

Um **nucleotídeo** é uma pequena molécula com três partes ligadas entre si:

- **Um grupo fosfato**: carrega carga negativa e faz a "ligação" entre um nucleotídeo e o próximo.
- **Um açúcar de cinco carbonos (pentose)**: é a **desoxirribose** no DNA e a **ribose** no RNA. A diferença é apenas um átomo de oxigênio a menos na desoxirribose (o prefixo "desoxi" significa exatamente isso).
- **Uma base nitrogenada**: é a parte que funciona como "letra". Existem cinco, em duas famílias:
  - **Purinas** (duas anéis): **A**denina e **G**uanina.
  - **Pirimidinas** (um anel): **C**itosina, **T**imina (só no DNA) e **U**racila (só no RNA).

Pense num trem: o fosfato e o açúcar formam o "vagão" e a base é o "passageiro" que viaja pendurado. Em engenharia, seria o equivalente a uma sequência de símbolos em um fita, onde o **suporte** (a fita) é sempre igual e a **informação** está no símbolo.

### 1.2 O DNA: a dupla hélice

O **DNA** (ácido desoxirribonucleico) é formado por **duas fitas** de nucleotídeos enroladas uma na outra, formando uma **dupla hélice**. Pontos essenciais:

- **Complementaridade**: A sempre pareia com T, e G sempre pareia com C. O pareamento é feito por **ligações de hidrogênio**: A–T tem 2 ligações e G–C tem 3. Por isso trechos ricos em G–C são mais estáveis (à temperatura mais alta para separar as fitas).
- **Antiparalelismo**: as duas fitas correm em sentidos opostos. Cada fita tem uma ponta chamada **5’** (cinco linha) e outra chamada **3’** (três linha), nomes que vêm dos carbonos do açúcar. Uma fita vai de 5’→3’ e a parceira de 3’→5’. Isso importa porque as enzimas só "leem" e "escrevem" em um sentido.
- **Consequência prática**: se você conhece uma fita, conhece a outra (a **fita complementar reversa**). Os programas de bioinformática fazem isso o tempo todo, porque um gene pode estar em qualquer uma das duas fitas.

Exemplo: fita 5’-ATGGCA-3’ tem como parceira 3’-TACCGT-5’. Lida no sentido 5’→3’, a parceira é 5’-TGCCAT-3’ (é o "reverso complementar").

O DNA pode assumir formas diferentes (as mais conhecidas são **B**, **A** e **Z**) e, nas células, costuma estar **super-enovelado** (*supercoiling*), como um fio de telefone torcido que forma voltas sobre si mesmo. Isso permite que moléculas de milhões de letras caibam dentro de uma célula microscópica.

### 1.3 O RNA: a cópia de trabalho

O **RNA** (ácido ribonucleico) difere do DNA em três aspectos: usa **ribose**, usa **uracila (U)** no lugar da timina, e normalmente tem **fita simples**. Por ter fita simples, uma fita de RNA pode dobrar sobre si mesma e parear com outro trecho da mesma fita, formando estruturas em **grampo (*hairpin*)**. Guarde essa palavra: os grampos aparecem mais adiante, nos terminadores de transcrição.

Há três tipos principais:

| Tipo | Função |
| --- | --- |
| **mRNA** (RNA mensageiro) | Cópia de um gene, leva a receita da proteína até o ribossomo |
| **tRNA** (RNA transportador) | "Adaptador" que traz o aminoácido certo para cada trecho de 3 letras |
| **rRNA** (RNA ribossomal) | Principal componente do ribossomo, a "máquina" que monta a proteína |

### 1.4 O dogma central: replicação, transcrição e tradução

O fluxo da informação na célula acontece em três processos:

1. **Replicação** (DNA → DNA): a célula copia todo o genoma antes de se dividir.
2. **Transcrição** (DNA → RNA): uma enzima chamada **RNA polimerase** lê um gene e produz uma cópia em RNA.
3. **Tradução** (RNA → proteína): o ribossomo lê o mRNA de três em três letras e junta aminoácidos formando a proteína.

Cada trio de letras é um **códon**. Existem 64 códons possíveis (4³) e 20 aminoácidos, portanto vários códons significam o mesmo aminoácido (o código é "redundante"). O códon **AUG** (ATG no DNA) significa "começar" (e codifica o aminoácido metionina). Três códons — **UAA, UAG e UGA** — significam "parar" (**códons de parada ou *stop***).

A região do gene que vai do códon de início ao de parada, sem interrupção, chama-se **CDS** (*coding sequence*, sequência codificadora). Encontrar CDSs é a tarefa central da anotação de genomas bacterianos.

Em uma analogia de engenharia: o DNA é o **arquivo mestre** guardado em local seguro; o mRNA é uma **cópia temporária** impressa para uso na linha de produção; a proteína é o **produto** montado.

### 1.5 Por que isso importa para a prática?

O que você receberá do sequenciador é só texto: milhões de pedaços de A, C, G e T. Montar esse texto na ordem certa e depois descobrir onde estão os genes (trechos que começam com códon de início e terminam com códon de parada) é exatamente o que faremos.

## Parte 2 — Organização dos genomas: procariotos × eucariotos

O **genoma** é o conjunto completo do material genético de um organismo. Os seres vivos celulares se dividem em dois grandes grupos, e a organização do genoma é muito diferente entre eles.

- **Procariotos** (bactérias e arqueias): células **sem núcleo**. Exemplo: *Escherichia coli*.
- **Eucariotos** (animais, plantas, fungos, protozoários): células **com núcleo** delimitado por membrana. Exemplo: o ser humano.

### 2.1 Genoma procarioto

- **Onde fica**: numa região chamada **nucleoide**, que não é separada do resto da célula por membrana.
- **Cromossomo**: em geral **um único cromossomo circular** (uma fita de DNA com as pontas unidas, como um anel). Em *E. coli*, são cerca de 4,6 milhões de pares de bases (4,6 Mpb).
- **Plasmídeos**: pequenos círculos extras de DNA, independentes do cromossomo. Costumam carregar genes de **resistência a antibióticos** ou de toxinas, e podem ser passados entre bactérias.
- **Genoma compacto**: cerca de **88%** do DNA bacteriano codifica proteínas. Há pouco "espaço vazio".
- **Empacotamento**: o DNA é super-enovelado e organizado em alças por proteínas associadas ao nucleoide, sem a hierarquia complexa dos eucariotos.
- **Transferência horizontal de genes**: bactérias trocam DNA entre si, mesmo de espécies diferentes, por três mecanismos: **conjugação** (contato direto, como uma "ponte"), **transformação** (capturar DNA solto do ambiente) e **transdução** (DNA levado por vírus que infectam bactérias). É por isso que a resistência a antibióticos se espalha tão rápido.

### 2.2 Genoma eucarioto

- **Onde fica**: dentro do **núcleo**, protegido por uma membrana dupla. Parte menor do DNA fica nas **mitocôndrias** (e, nas plantas, nos **cloroplastos**), que têm DNA circular próprio. A teoria da **endossimbiose** explica isso: essas organelas descendem de bactérias que passaram a viver dentro de células maiores.
- **Cromossomos lineares e múltiplos**: o ser humano tem 46 cromossomos (23 pares) no núcleo.
- **Genoma grande e "espaçado"**: o genoma humano tem cerca de 3 bilhões de pares de bases, mas apenas \~1% a 2% codifica proteínas.
- **Tamanho do genoma não é sinal de complexidade**: algumas plantas e anfíbios têm genomas muito maiores que o humano, e o número de genes é parecido entre organismos bem diferentes.

Composição aproximada do genoma humano:

| Componente | Parcela aproximada |
| --- | --- |
| Éxons (partes que codificam proteínas) | \~1% |
| Íntrons (trechos internos dos genes, removidos depois) | \~24% |
| Regiões entre genes (intergênicas) | \~22% |
| Elementos transponíveis ("genes saltadores") | \~45% |
| Grandes duplicações | \~5% |
| Repetições simples | \~3% |

Os valores variam um pouco conforme a fonte, mas a mensagem é clara: **a maior parte do genoma humano não codifica proteína**. Isso torna a anotação de eucariotos muito mais difícil do que a de bactérias.

### 2.3 Como o DNA cabe na célula? A cromatina

Uma célula humana tem cerca de 2 metros de DNA esticado dentro de um núcleo de poucos micrômetros. O empacotamento é feito em níveis:

1. **Nucleossomo**: o DNA dá cerca de duas voltas em torno de um conjunto de oito proteínas chamadas **histonas** (duas cópias de H2A, H2B, H3 e H4). É o "carretel" básico, uma estrutura em "contas de um colar". O DNA entre dois carretéis chama-se **DNA ligante** (15 a 55 pares de bases) e é estabilizado pela histona **H1**.
2. **Fibra de 30 nm (solenoide)**: o colar de nucleossomos se enrola formando um cilindro.
3. **Alças presas a um esqueleto (scaffold)**: o solenoide forma alças ancoradas em uma matriz proteica que inclui a enzima **topoisomerase II**.
4. **Cromossomo condensado**: na divisão celular, tudo se compacta ainda mais (cerca de 700 nm e depois \~1400 nm), formando os cromossomos visíveis ao microscópio.

**Código das histonas**: modificações químicas nas "caudas" das histonas regulam se o DNA está acessível ou não. De forma simplificada, a metilação da lisina 9 da histona H3 (H3K9me) está associada a **silêncio** (gene desligado), enquanto marcas como H3K4me e acetilações estão associadas a **expressão** (gene ligado). Isso mostra que, nos eucariotos, a regulação não depende só da sequência, mas também de como o DNA está empacotado.

### 2.4 Regiões especiais do cromossomo eucarioto

- **Centrômero**: região onde as fibras que puxam os cromossomos se prendem na divisão celular. Muitas vezes é rico em sequências repetitivas.
- **Telômeros**: as pontas dos cromossomos, formadas por repetições curtas (nos humanos, TTAGGG). Funcionam como "capas" que protegem as extremidades. A enzima **telomerase** as reconstrói.
- **Origens de replicação**: pontos onde a cópia do DNA começa. Cromossomos eucariotos têm muitas; o cromossomo bacteriano tem em geral só uma.

### 2.5 Resumo comparativo

| Característica | Procarioto | Eucarioto |
| --- | --- | --- |
| Núcleo | Ausente (nucleoide) | Presente |
| Cromossomos | Geralmente 1, circular | Vários, lineares |
| DNA extra | Plasmídeos | Mitocôndria / cloroplasto |
| Tamanho típico | 0,5–10 Mpb | De milhões a bilhões de pb |
| Densidade de genes | Muito alta (\~88% codificante) | Baixa (humano: \~1–2%) |
| Empacotamento | Super-enovelamento e alças | Nucleossomos, histonas, cromatina |
| Íntrons | Praticamente ausentes | Muito comuns |

## Parte 3 — Organização dos genes: anatomia de um gene

Um **gene** não é só a parte que codifica a proteína. Ao redor dela existem trechos de DNA que funcionam como **"sinais de trânsito"**, dizendo à maquinaria da célula onde começar a copiar, onde a tradução inicia e onde tudo termina. A figura abaixo compara um gene de bactéria com um gene de eucarioto.

&#91;embedded content: gene procarioto vs eucarioto · 4 e 6 elementos\]

Observe que, no procarioto, a região colorida é uma peça única. No eucarioto, ela é dividida em éxons, com íntrons removidos no meio.

### 3.1 O gene procarioto

Lendo da esquerda para a direita (sentido 5’→3’ da fita):

1. **Promotor**: sequência onde a RNA polimerase se liga para começar a transcrição. Em bactérias, a RNA polimerase usa uma subunidade auxiliar chamada **fator sigma**, que reconhece duas "caixas" de sequência no DNA: a **caixa −35** (consenso TTGACA) e a **caixa −10**, também chamada *Pribnow box* (consenso TATAAT). O sinal negativo indica que ficam *antes* (a montante) do ponto onde a transcrição começa (posição +1). Quanto mais parecido com o consenso, mais forte costuma ser o promotor.
2. **Operador**: um trecho curto, na região do promotor, onde proteínas **repressoras** podem se ligar e bloquear a transcrição. É um "interruptor".
3. **Sítio de ligação do ribossomo (Shine-Dalgarno)**: após o início da transcrição, há no mRNA uma sequência rica em purinas (consenso AGGAGG), cerca de 5 a 10 letras antes do códon de início. Ela se pareia com o rRNA do ribossomo e o posiciona no lugar correto.
4. **Códon de início (ATG)**, seguido da **CDS** (região codificadora) e do **códon de parada**. Não há íntrons: a CDS é contínua.
5. **Terminador**: sinal que manda a RNA polimerase parar. Em muitas bactérias é um **terminador independente de Rho**: uma sequência que forma um **grampo** no RNA (trecho rico em G–C que dobra e pareia consigo mesmo), seguido de uma sequência de **uracilas (U)**. O grampo faz a polimerase "travar" e a cauda de U, que pareia fracamente, deixa o RNA se soltar.

**Operon**: nas bactérias, genes com funções relacionadas costumam ficar lado a lado e ser transcritos juntos, sob **um único promotor**, gerando um mRNA **policistrônico** (que carrega vários genes). O exemplo clássico é o **operon lac** de *E. coli*, com três genes (*lacZ*, *lacY*, *lacA*) para digerir lactose. Quando não há lactose, um repressor se liga ao operador e desliga o operon; quando há, o repressor sai e os três genes são ligados juntos. Isso é muito eficiente em organismos que precisam reagir rápido ao ambiente.

### 3.2 O gene eucarioto

Os eucariotos têm o mesmo princípio, com mais camadas:

1. **Promotor**: em geral inclui a **TATA box** (consenso TATAAA, cerca de 25–30 pb antes do início), onde se liga um conjunto de proteínas chamadas **fatores gerais de transcrição**, que recrutam a RNA polimerase II. Nem todos os promotores têm TATA box.
2. **Regiões reguladoras distantes (*enhancers*, "intensificadores")**: sequências que podem estar a milhares de pares de bases do gene e aumentam a transcrição, por dobramento do DNA que aproxima as regiões.
3. **5’ UTR e 3’ UTR**: *untranslated regions* (regiões não traduzidas). São trechos do mRNA antes do códon de início e depois do códon de parada. Não viram proteína, mas controlam estabilidade e eficiência da tradução.
4. **Éxons e íntrons**: o gene é interrompido. **Éxons** são os trechos que permanecem no mRNA maduro; **íntrons** são trechos intercalados que são **removidos**. No gene humano médio, existem cerca de 9 éxons, e os íntrons somam muito mais DNA que os éxons.
5. **Sítios de *splicing***: o corte e a emenda (*splicing*) são feitos pelo **spliceosomo**, um complexo de RNA e proteínas. Ele reconhece sinais nas bordas de cada íntron: o **sítio doador** (5’ do íntron, quase sempre começa com **GT**), o **sítio aceptor** (3’ do íntron, quase sempre termina com **AG**) e um ponto de ramificação interno com uma adenina. Isso é conhecido como "regra GT–AG".
6. **Processamento do mRNA**: ao mRNA acrescenta-se um **cap** ("capuz") na ponta 5’ e uma **cauda poli-A** na 3’, que protegem a molécula. Tudo isso ocorre no núcleo, antes do mRNA ir ao citoplasma.
7. **Splicing alternativo**: o mesmo gene pode gerar várias proteínas, escolhendo quais éxons manter. É uma das razões pelas quais humanos têm apenas \~20 mil genes, mas centenas de milhares de proteínas diferentes.
8. **Monocistrônico**: normalmente cada mRNA eucarioto carrega um só gene.

### 3.3 Por que isso complica a anotação de eucariotos?

Em bactérias, um programa procura simplesmente trechos longos que começam com ATG e terminam em códon de parada sem interrupção (ORFs). Em eucariotos, o gene está "picotado" em éxons separados por íntrons enormes, e são necessários modelos estatísticos e evidências de RNA para localizar cada éxon e cada sítio de *splicing*. Por isso, para eucariotos usam-se ferramentas diferentes (MAKER, AUGUSTUS, BRAKER), vistas adiante.

### 3.4 Comparativo

| Elemento | Procarioto | Eucarioto |
| --- | --- | --- |
| Promotor | Caixas −35 e −10, reconhecidas pelo fator sigma | TATA box e outros elementos, fatores gerais de transcrição |
| Regulação | Operador e repressores | Enhancers, cromatina, muitos fatores |
| Início da tradução | Shine-Dalgarno + ATG | Cap + sequência ao redor do ATG |
| Estrutura do gene | CDS contínua | Éxons e íntrons |
| Splicing | Não | Sim (GT–AG, spliceosomo) |
| Terminador | Grampo + poli-U (indep. de Rho) ou Rho-dependente | Sinal de poliadenilação (AAUAAA) |
| mRNA | Frequentemente policistrônico (operons) | Monocistrônico |
| Transcrição e tradução | Simultâneas, no mesmo compartimento | Separadas: núcleo e citoplasma |

## Parte 4 — Montagem e anotação de genomas

Até aqui falamos de como o genoma é organizado *dentro* da célula. Agora vamos ao que o bioinformata faz com ele *no computador*.

### 4.1 O problema: sequenciadores não leem cromossomos inteiros

Nenhuma tecnologia lê um cromossomo de ponta a ponta. Os sequenciadores modernos (por exemplo, **Illumina**) quebram o DNA em milhões de pedaços pequenos e leem cada pedaço. Cada pedaço lido é um ***read*** ("leitura"), com algo entre 100 e 300 letras. Uma analogia: é como pegar 1.000 cópias do mesmo livro, passar todas no triturador e depois tentar reconstruir o texto original juntando os pedaços que se sobrepõem.

**Cobertura (*coverage*)**: quantas vezes, em média, cada letra do genoma foi lida. Se a soma de todos os reads for 100 Mpb e o genoma tiver 5 Mpb, a cobertura é de 20×. Cobertura alta ajuda a corrigir erros de leitura e a garantir que todas as regiões foram vistas.

**Sequenciamento *paired-end* (pares de extremidades)**: o equipamento lê o mesmo fragmento de DNA pelas **duas pontas**, gerando dois arquivos: **R1** (a leitura da ponta esquerda) e **R2** (a da ponta direita). Sabendo que os dois reads estão separados por uma distância aproximada, o programa de montagem consegue atravessar trechos repetitivos e ordenar melhor os pedaços. Na prática você usará exatamente esse par de arquivos.

### 4.2 Formatos de arquivo

- **FASTA** (.fasta, .fa, .fna): arquivo de texto simples com sequências. Cada sequência começa com uma linha de cabeçalho iniciada por `>`, seguida das letras. Exemplo:

```
>contig_1 comprimento=1250
ATGGCGTTAGCTAGCTAGGATCGATCG...
```

- **FASTQ** (.fastq, .fq): como o FASTA, mas **inclui a qualidade** de cada letra lida. Cada read ocupa 4 linhas: (1) identificador iniciado por `@`; (2) a sequência; (3) um sinal `+`; (4) uma linha de qualidade, em que cada símbolo codifica a confiança naquela letra (escala **Phred**: Q30 significa 1 erro em 1.000 letras). Os arquivos que você receberá para a aula estão neste formato.
- **GFF/GFF3** (.gff): tabela que descreve **onde estão os elementos anotados** (gene, CDS, tRNA...) em cada sequência: nome do contig, início, fim, fita (+ ou −), tipo e atributos.
- **GenBank** (.gbk): formato ricamente anotado, que traz a sequência e as anotações juntas, no padrão do banco NCBI.

### 4.3 Montagem (*assembly*)

**Montar** é juntar os reads em sequências cada vez mais longas. A saída depende do nível:

- **Contigs** (*contiguous sequences*): trechos contínuos obtidos por sobreposição de reads, sem lacunas. São o produto direto da montagem.
- **Scaffolds**: contigs ordenados e orientados com base na informação dos pares (R1–R2), com **lacunas de tamanho estimado** (representadas por Ns) entre eles.
- **Genoma completo (*closed*)**: um único cromossomo, sem lacunas. Exige tecnologias complementares (leituras longas) e trabalho adicional.

Existem duas estratégias principais:

| Estratégia | Como funciona | Quando usar |
| --- | --- | --- |
| **Montagem por referência** | Alinha os reads sobre um genoma já conhecido da mesma espécie | Quando já existe um genoma de referência bom (ex.: humano) |
| **Montagem *de novo*** | Reconstrói o genoma apenas a partir dos reads, sem referência | Organismos novos ou cepas muito diferentes (ex.: nossa prática) |

O truque dos montadores modernos é o **grafo de De Bruijn**. Em vez de comparar todos os reads com todos (inviável com milhões de reads), o programa corta cada read em "palavrinhas" de tamanho fixo chamadas **k-mers** (por exemplo, k = 55 letras) e constrói um **grafo**, ou seja, uma rede de nós e setas, em que cada k-mer se liga ao seguinte que se sobrepõe a ele. O genoma vira, então, um **caminho** nesse grafo. Onde o caminho se bifurca (por causa de repetições ou erros), o programa precisa decidir ou cortar o contig. É por isso que genomas com muitas repetições (como os eucariotos) são difíceis de montar: o mesmo trecho aparece em muitos lugares e o grafo fica emaranhado.

### 4.4 Como saber se a montagem é boa? As métricas

- **Número de contigs**: quanto menos, melhor (o ideal seria 1 por cromossomo).
- **Tamanho total**: deve ser próximo ao tamanho esperado do genoma da espécie.
- **Maior contig**: quanto maior, melhor.
- **N50**: a métrica mais citada. Ordene os contigs do maior para o menor e some seus tamanhos até atingir **50% do total** montado. O tamanho do contig em que isso acontece é o N50. Interpretação: *"metade do genoma está em contigs de pelo menos esse tamanho"*. **Quanto maior o N50, mais contínua é a montagem.**
- **L50**: quantos contigs foram necessários para chegar a esses 50%.
- **GC (%)**: porcentagem de G e C no genoma. É característica de cada espécie e serve como "impressão digital"; se a montagem tiver GC muito diferente do esperado, pode haver contaminação.

#### Exercício resolvido: calculando o N50

Suponha uma montagem com 5 contigs de 40, 30, 15, 10 e 5 kb (kb = mil letras).

1. Soma total = 40 + 30 + 15 + 10 + 5 = **100 kb**. Metade = 50 kb.
2. Ordenados: 40, 30, 15, 10, 5. Soma acumulada: 40 → 70 → ...
3. A soma ultrapassa 50 kb quando chegamos ao segundo contig (40 + 30 = 70). Logo **N50 = 30 kb** e **L50 = 2**.

**Tente você:** contigs de 60, 20, 10, 6 e 4 kb. Qual o N50? *(Resposta: soma total = 100 kb, metade = 50; o primeiro contig já tem 60 kb, então N50 = 60 kb e L50 = 1.)*

### 4.5 Anotação (*annotation*)

Ter a sequência montada é como ter um livro sem capítulos nem índice. **Anotar** é adicionar a informação biológica sobre onde estão os elementos e para que servem. Divide-se em dois passos:

1. **Anotação estrutural** ("onde está cada coisa"): localizar genes, regiões codificadoras (CDS), RNAs (tRNA, rRNA), etc. Em bactérias, o programa procura **ORFs** (*open reading frames*, fases abertas de leitura): trechos longos que começam em códon de início e terminam em códon de parada sem interrupção. Como o DNA tem duas fitas e cada fita pode ser lida em três "molduras" (frames), são seis possibilidades por posição.
2. **Anotação funcional** ("para que serve"): comparar cada proteína prevista com bancos de dados de proteínas já conhecidas (por similaridade de sequência). Se a sequência é parecida com uma proteína de função conhecida, atribui-se a ela uma função provável (ex.: "DNA polymerase III subunit beta"). Se não há correspondência, o gene recebe o nome **"hypothetical protein"** (proteína hipotética), pois sabemos que provavelmente é um gene, mas não sabemos o que faz.

**Para eucariotos**, a anotação estrutural é mais complexa (éxons, íntrons, splicing) e depende de ferramentas que combinam **modelos estatísticos** (como Modelos Ocultos de Markov) com **evidências experimentais** (RNA-seq, proteínas de espécies parentes). As mais usadas são **AUGUSTUS** (preditor *ab initio*, treinado para cada espécie), **BRAKER** (que treina o AUGUSTUS automaticamente a partir de evidências) e **MAKER** (pipeline que integra várias fontes de evidência). Elas não serão usadas na prática, mas é importante saber que existem.

### 4.6 Um ponto de atenção importante

As ferramentas da nossa prática (SPAdes, Quast e Prokka) **encontram genes que codificam proteínas, tRNAs e rRNAs**. Elas **não** localizam promotores, terminadores ou sítios de *splicing*. Para essas regiões reguladoras existem ferramentas específicas, que mostraremos ao final (BPROM e ARNold).

## Parte 5 — As ferramentas, uma a uma

Cada ferramenta abaixo é descrita com o mesmo roteiro: **o que é**, **o que recebe**, **o que devolve** e **como interpretar**.

### 5.1 Galaxy (a plataforma)

- **O que é**: uma plataforma web gratuita, de código aberto, que permite rodar programas de bioinformática **sem instalar nada e sem digitar comandos**. Cada programa vira um formulário com campos e botões. Você já a utilizou na Aula 2.
- **Como funciona**: o Galaxy roda os programas em servidores remotos. Você não precisa de computador potente. Três conceitos: (1) **Histórico** — a coluna à direita, onde ficam todos os arquivos (de entrada e resultados); (2) **Ferramentas** — a coluna à esquerda, com caixa de busca; (3) **Painel central** — onde você preenche o formulário da ferramenta e vê os resultados.
- **Cores no histórico**: cinza = na fila; amarelo = executando; **verde = terminou com sucesso**; vermelho = erro.
- **Por que usar**: reprodutibilidade (o Galaxy registra exatamente o que foi feito) e acessibilidade. Em pesquisa real, o mesmo fluxo pode ser rodado em linha de comando.

### 5.2 SPAdes (montagem)

- **O que é**: um **montador de genomas *de novo*** projetado para bactérias e outros genomas pequenos. O nome vem de *St. Petersburg genome assembler*. Usa grafos de De Bruijn com vários tamanhos de k-mer.
- **O que recebe**: os reads em FASTQ (no nosso caso, os dois arquivos *paired-end*, R1 e R2).
- **O que devolve**: o arquivo de **contigs** (FASTA) e, normalmente, de **scaffolds**; também o grafo de montagem e um log.
- **Parâmetros relevantes**: o *modo* (usaremos o padrão para isolados bacterianos) e os tamanhos de k-mer (deixaremos em automático). Não precisa mexer em mais nada.
- **Como interpretar**: a saída não diz se a montagem é boa. Para isso, usamos o Quast.

### 5.3 Quast (avaliação da montagem)

- **O que é**: *QUality ASsessment Tool*. Calcula métricas que resumem a qualidade de uma montagem.
- **O que recebe**: um ou mais arquivos de contigs (FASTA). Opcionalmente, um genoma de referência e anotações; sem eles, o Quast ainda calcula as métricas internas.
- **O que devolve**: um relatório (HTML e tabela) com: número de contigs, tamanho total, maior contig, **N50**, **L50**, %GC e gráficos de distribuição de tamanhos (como o gráfico de Nx e o *Cumulative length*).
- **Como interpretar**: queremos **poucos contigs**, **N50 alto**, tamanho total plausível e GC compatível com a espécie. Contigs muito pequenos (< 500 pb) costumam ser lixo e podem ser descartados.

### 5.4 Prokka (anotação de genomas bacterianos)

- **O que é**: um programa que **anota genomas de bactérias, arqueias e vírus** de forma automatizada. Internamente, coordena várias ferramentas: **Prodigal** (prediz CDSs), **Aragorn** (acha tRNAs), **Barrnap** (acha rRNAs) e **bancos de dados de proteínas** para dar nome a cada proteína prevista por similaridade.
- **O que recebe**: o arquivo FASTA com os contigs. Você informa o **código genético**: usaremos o **11** (o código de bactérias e arqueias, que difere ligeiramente do código "padrão" quanto aos códons de início).
- **O que devolve**: vários arquivos: **GFF** (as anotações, para ver em um navegador de genoma), **GBK** (GenBank), **FNA** (nucleotídeos dos contigs), **FAA** (sequências das proteínas), **FFN** (nucleotídeos dos genes), **TSV** (tabela de features) e **TXT** (resumo com contagens).
- **Como interpretar**: o resumo mostra quantos CDS, tRNA e rRNA foram achados. Para uma bactéria de \~4 Mpb esperam-se cerca de 4.000 CDSs. Muitos serão "hypothetical protein", e isso é normal.
- **O que o Prokka NÃO faz**: não prediz promotores, terminadores de transcrição nem sítios de *splicing*, e não serve para eucariotos.

### 5.5 JBrowse (visualização)

- **O que é**: um **navegador de genoma** (*genome browser*) interativo, disponível dentro do Galaxy. Mostra o genoma como uma "régua" horizontal, e as anotações como caixas coloridas sobre ela, com zoom.
- **O que recebe**: o FASTA dos contigs e o GFF do Prokka.
- **Como usar**: escolha um contig, aproxime o zoom até enxergar genes individuais (setas indicam a fita) e clique em um gene para ver o nome e a função.

### 5.6 NCBI Genome e Ensembl (bancos de referência)

- **NCBI Genome / GenBank**: acervo público do *National Center for Biotechnology Information* (EUA) com genomas de milhares de espécies. Permite baixar sequências (FASTA), anotações (GFF, GBK) e ver metadados.
- **Ensembl**: banco europeu (EMBL-EBI) especializado em **genomas de eucariotos**, em especial vertebrados. Tem navegador gráfico onde é possível ver um gene com seus éxons, íntrons, transcritos alternativos e regiões reguladoras. Vamos apenas mostrar a página de um gene humano, para você ver na prática a estrutura éxon–íntron de que falamos.

### 5.7 BPROM e ARNold (regiões reguladoras de bactérias)

São dois serviços web, **demonstrados ao final** (não fazem parte do fluxo principal):

- **BPROM** (Softberry): prediz **promotores bacterianos**, buscando as caixas −10 e −35 do **fator sigma 70** (o mais comum em *E. coli*). Recebe um trecho de DNA em formato FASTA e devolve as posições candidatas com uma pontuação (**LDF**, *linear discriminant function*); quanto maior, mais confiável. **Cuidados**: ele deve ser usado com trechos curtos, em geral os **até 150 pb antes** do início de um gene, e não com o genoma inteiro.
- **ARNold**: prediz **terminadores independentes de Rho**, procurando grampos de RNA seguidos de poli-U. Recebe uma sequência FASTA e devolve as posições dos grampos candidatos com a energia livre de dobramento (**ΔG**); quanto **mais negativo** o ΔG, mais estável o grampo e mais provável que seja um terminador real.

Isso fecha o círculo com a Parte 3: a anotação completa de um genoma combina ferramentas para encontrar os **genes** (Prokka) e ferramentas para encontrar os **sinais que os regulam** (BPROM e ARNold). Todas essas previsões são **candidatas**: são hipóteses a serem confirmadas experimentalmente.

## Parte 6 — Prática guiada no Galaxy (80 minutos)

**Objetivo**: montar e anotar o genoma de uma bactéria, partindo de reads brutos, e comparar o resultado com um gene eucarioto. A figura mostra o caminho completo; cada passo está detalhado abaixo.

&#91;embedded content: fluxo da prática · reads, montagem, avaliação e anotação\]

Os contigs são o ponto central: tanto o Quast quanto o Prokka partem do mesmo arquivo.

**Os dados**: o conjunto de dados é um isolado mutante de *Staphylococcus aureus* (bactéria comum na pele, que pode causar infecções), sequenciado em *paired-end*, com reads de 150 letras e cobertura de \~19×. Foi escolhido por ser pequeno e processar em poucos minutos. Ele está publicado no repositório **Zenodo** (registro 582600), com duas URLs:

```
https://zenodo.org/record/582600/files/mutant_R1.fastq
https://zenodo.org/record/582600/files/mutant_R2.fastq
```

**Pontos de atenção antes de começar**

- Os nomes de campos podem variar um pouco conforme a versão das ferramentas no servidor. Se algo não estiver exatamente como descrito, procure o campo equivalente ou chame o professor.
- Clique em *Run Tool* **uma única vez**. Cliques repetidos criam jobs duplicados e deixam o servidor mais lento.
- Seus números podem diferir levemente dos de colegas por causa de versões e de processos aleatórios nos algoritmos. Isso é normal.

### Passo 1 — Criar a *history* e importar os reads (10 min)

1. Acesse **usegalaxy.org** e faça login.
2. No painel *History* (à direita), clique no **+** para criar uma nova history. Clique em *Unnamed history* e renomeie para **Aula3\_Genomica**.
3. No canto superior esquerdo, clique em **Upload** (seta para cima) → **Paste/Fetch data**.
4. Cole as duas URLs acima, uma por linha.
5. No seletor *Type*, escolha **fastqsanger** (o formato FASTQ que o Galaxy espera).
6. Clique em **Start** e depois em **Close**.

**O que deve acontecer**: dois itens aparecem na history, passando de cinza (fila) para amarelo (baixando) e depois verde.

**Checkpoint**: a history tem exatamente dois datasets verdes; o formato mostrado é *fastqsanger*; ao clicar no olho, você vê linhas começando por `@`.

### Passo 2 — Renomear os arquivos (5 min)

Por padrão o Galaxy usa a URL inteira como nome, o que atrapalha depois.

1. Clique no **lápis** (editar atributos) do dataset terminado em *mutant\_R1.fastq*.
2. No campo *Name*, escreva **reads\_R1** e clique em *Save*.
3. Repita para o outro, chamando-o **reads\_R2**.

**Atenção à ordem**: R1 e R2 formam um par. Se trocar os nomes, o SPAdes receberá as pontas invertidas e a montagem piora ou falha.

*Opcional se sobrar tempo*: rode o **FastQC** (da Aula 2) para ver a qualidade dos reads.

### Passo 3 — Montar o genoma com o SPAdes (25 min)

1. Na caixa de busca de ferramentas, digite **SPAdes** e abra *SPAdes genome assembler for regular and single-cell projects*.
2. Escolha o modo **paired-end**, com datasets individuais.
3. Em *Forward reads*, selecione **reads\_R1**; em *Reverse reads*, **reads\_R2**.
4. Mantenha os demais parâmetros no padrão (inclusive k-mers automáticos).
5. Clique em **Run Tool**, uma vez só.

Com este conjunto pequeno, a montagem leva poucos minutos; o professor pode deixar uma execução pronta para você explorar enquanto a sua roda.

**O que o SPAdes entrega**

| Saída | O que é | Vamos usar? |
| --- | --- | --- |
| Contigs (FASTA) | Sequências montadas | Sim, no Quast e no Prokka |
| Contig stats | Tamanho e cobertura de cada contig | Sim, para observar |
| Scaffolds | Contigs ordenados com lacunas | Não hoje |
| Assembly graph | Grafo da montagem | Opcional |
| Log | Registro da execução | Só se der erro |

**Olhando os contigs**: clique no olho do arquivo de contigs. Cada cabeçalho segue o padrão `NODE_1_length_XXXXX_cov_YY.YYY`: o número depois de *length* é o tamanho em letras e o depois de *cov* é a cobertura média daquele contig.

**Perguntas para anotar**: (1) Quantos contigs foram gerados? Quanto mede o maior? (2) Algum contig tem cobertura muito maior que os demais? O que isso pode indicar? *(Dica: trechos repetidos ou plasmídeos, presentes em várias cópias.)*

**Checkpoint**: a saída de contigs está verde e não vazia; você reconheceu os cabeçalhos `>NODE_`.

### Passo 4 — Avaliar a montagem com o Quast (10 min)

1. Busque **Quast** e abra *Quast genome quality assessment tool*.
2. Em *Contigs/scaffolds file*, selecione o arquivo de **contigs** do SPAdes.
3. Tipo de montagem: **Genome**. Usar genoma de referência: **No**.
4. Clique em **Run Tool** e abra o relatório **HTML** (ícone do olho).

**Como ler o relatório**

| Métrica | Como interpretar |
| --- | --- |
| # contigs | Menos é melhor (menos fragmentada) |
| Largest contig | Maior é melhor |
| Total length | Deve ser próximo ao tamanho esperado do genoma |
| N50 | Maior é melhor |
| L50 | Menor é melhor |
| GC (%) | Característico da espécie; valor muito fora indica possível contaminação |

O Quast costuma ignorar, nas estatísticas principais, contigs com menos de 500 letras e mostra o total completo em linhas separadas. Lembre que a cobertura aqui é de só \~19×: é comum a montagem sair mais fragmentada do que sairia com mais dados.

**Anote seus resultados**

| Métrica | Meu resultado |
| --- | --- |
| # contigs |  |
| Largest contig |  |
| Total length |  |
| N50 |  |
| L50 |  |

**Checkpoint**: você abriu o relatório, localizou o N50 e consegue explicar o que ele diz sobre a montagem.

### Passo 5 — Anotar o genoma com o Prokka (15 min)

1. Busque **Prokka** e abra *Prokka prokaryotic genome annotation tool*.
2. Em *Contigs to annotate*, selecione o arquivo de contigs do SPAdes.
3. Em *Genetic code*, escolha **11** (bactérias, arqueias e plastídeos).
4. Se houver campos de gênero e espécie, preencha *Staphylococcus* e *aureus* (opcional).
5. Clique em **Run Tool**.

**Os arquivos de saída**: `tsv` (tabela de features, uma linha por gene — o que vamos ler), `txt` (resumo em números), `gff` (posições), `gbk` (GenBank), `fna` (contigs), `faa` (proteínas) e `ffn` (genes).

**Lendo a tabela `tsv`**: `locus_tag` (identificador do gene), `ftype` (CDS, rRNA, tRNA ou tmRNA), `length_bp` (tamanho), `gene` (nome, quando há), `EC_number` (número da enzima), `COG` (categoria funcional) e `product` (função provável). Não estranhe os **hypothetical protein**: muitos genes não têm função conhecida.

**Atividade em sala**: abra o `txt` e o `tsv` e responda: (1) Quantos CDS, rRNA e tRNA foram encontrados? (2) Quantas linhas são *hypothetical protein*, e é uma fração grande ou pequena? (3) Escolha 3 a 5 genes com nome e função e registre `locus_tag`, `ftype`, `gene`, `product` e `length_bp`.

**Conecte com a teoria**: no `gff`, compare início e fim de dois genes vizinhos. Há muito espaço entre eles? E por que o Prokka não precisou prever sítios de *splicing*?

**Opcional — JBrowse**: busque *JBrowse genome browser*, use o `fna` do Prokka como sequência de referência e o `gff` como trilha de genes. Escolha o contig mais longo e aproxime o zoom; clique com o botão direito em um gene e use *View Details*.

**Checkpoint**: o Prokka terminou (verde); você registrou o total de CDS, rRNA e tRNA e listou de 3 a 5 genes.

### Passo 6 — Contrastar com um gene eucarioto no Ensembl (10 min)

Aqui não rodamos nada: apenas navegamos por uma anotação pronta.

1. Em outra aba, abra **ensembl.org**.
2. Na busca, selecione a espécie **Human**, digite **TP53** e clique em *Go*.
3. Abra a página do gene TP53 (ID ENSG00000141510), e depois o **transcrito principal**.
4. Procure a visualização da estrutura do gene (éxons como caixas, íntrons como linhas) e a lista de **Exons**. Se não achar, procure no menu lateral por *Exons* ou *Splice variants*.

**Preencha**

| Característica | Sua bactéria (Prokka) | Gene humano (Ensembl) |
| --- | --- | --- |
| Tamanho do genoma |  | \~3 bilhões de pb |
| Número de genes codificadores |  | \~20 mil |
| Gene contínuo ou em éxons/íntrons? |  |  |
| Tempo para anotar | Minutos | Anos de trabalho colaborativo |

**Para discutir**: quantos éxons tem o TP53 no transcrito aberto? Que fração do gene é de fato codificadora? Por que o Prokka não anotaria corretamente um genoma humano? (A anotação do Ensembl combina programas de predição, dados de RNA-seq e revisão manual de especialistas.)

### Fechamento (5 min)

Discussão em grupo sobre os resultados do N50, a proporção de *hypothetical protein* e o contraste procarioto × eucarioto. Se houver tempo, o professor demonstra o **BPROM** e o **ARNold** sobre um trecho de DNA a montante de um gene.

| Passo | Tempo |
| --- | --- |
| 1. History e importação | 10 min |
| 2. Renomear | 5 min |
| 3. SPAdes | 25 min |
| 4. Quast | 10 min |
| 5. Prokka | 15 min |
| 6. Ensembl | 10 min |
| Fechamento | 5 min |
| **Total** | **80 min** |

## Problemas comuns e como resolver

A maioria dos problemas vem de três causas: **formato de arquivo errado**, **arquivo errado selecionado** e **servidor ocupado**.

| O que você vê | Causa provável | O que fazer |
| --- | --- | --- |
| Dataset do upload ficou vermelho | URL errada ou Zenodo lento/indisponível | Confira a URL, tente de novo; se persistir, peça os arquivos ao professor |
| Formato não é fastqsanger | Galaxy detectou o tipo errado | Lápis → aba *Datatypes* → escolha fastqsanger → salvar |
| Os reads não aparecem na lista da ferramenta | Formato incompatível | Corrija o formato como acima |
| SPAdes ou Quast cinza por muito tempo | Servidor público ocupado | Aguarde e não reenvie; use a execução pronta do professor |
| Vários jobs iguais na history | Clicou em executar mais de uma vez | Apague os duplicados com o X e espere o primeiro |
| SPAdes ficou vermelho | R1 e R2 trocados, ou o mesmo arquivo nos dois campos | Confira os campos e leia as últimas linhas do erro |
| Quast ou Prokka não aceitam seu arquivo | Selecionou o arquivo errado (log, grafo) | Escolha o arquivo de contigs em FASTA |
| Não acho meus arquivos | Está em outra history | Menu do usuário → lista de histories → Aula3\_Genomica |
| O Galaxy avisa que o espaço acabou | Limite de armazenamento | Apague datasets que não usa mais |

**Regra de ouro**: se um item ficar vermelho, não apague nem refaça tudo. Clique no nome do dataset, abra os detalhes e leia a mensagem. Depois chame o professor com a mensagem aberta na tela.

## Perguntas de revisão (com gabarito)

Tente responder antes de olhar as respostas.

1. **Cite as três partes de um nucleotídeo e diga qual delas funciona como "letra".**
2. **Se uma fita de DNA é 5’-ATGGCA-3’, qual é a fita complementar lida em 5’→3’?**
3. **O que é um códon e quantos códons de parada existem?**
4. **Cite três diferenças entre o genoma de uma bactéria e o de um humano.**
5. **O que é um operon? Por que não há operons típicos em humanos?**
6. **Para que servem o promotor e o terminador de um gene bacteriano?**
7. **O que são éxons, íntrons e sítios de splicing?**
8. **Diferencie read, contig e scaffold.**
9. **Uma montagem tem contigs de 50, 25, 15, 6 e 4 kb. Qual o N50?**
10. **O que o Prokka faz e o que ele não faz?**
11. **O que significa "hypothetical protein"? Isso é um erro?**
12. **Por que a anotação de um genoma eucarioto é mais difícil do que a de uma bactéria?**

### Gabarito

1. Fosfato, açúcar (pentose) e base nitrogenada. A **base** funciona como "letra".
2. 5’-TGCCAT-3’ (reverso complementar).
3. Códon é um trio de nucleotídeos que especifica um aminoácido (ou sinal de início/parada). Há **três** códons de parada: UAA, UAG e UGA.
4. Exemplos: bactéria tem cromossomo circular único, humano tem vários lineares; bactéria não tem núcleo; genoma bacteriano é compacto (\~88% codificante), o humano tem \~1–2% codificante; bactéria quase não tem íntrons; bactéria tem plasmídeos, e o eucarioto tem DNA em mitocôndrias; o DNA eucarioto está empacotado em histonas.
5. Operon é um conjunto de genes de funções relacionadas transcritos juntos sob um único promotor, em um mRNA policistrônico. Eucariotos normalmente têm mRNA monocistrônico, cada gene com sua regulação.
6. O **promotor** indica onde a RNA polimerase se liga para iniciar a transcrição; o **terminador** indica onde ela deve parar.
7. Éxons são trechos que permanecem no mRNA maduro; íntrons são removidos; os sítios de splicing (GT no início e AG no fim do íntron) são os sinais que o spliceosomo reconhece.
8. *Read*: fragmento lido pelo sequenciador. *Contig*: trecho contínuo montado por sobreposição de reads. *Scaffold*: contigs ordenados e orientados, com lacunas entre eles.
9. Total = 100 kb; metade = 50 kb. O primeiro contig (50 kb) já atinge 50%, então **N50 = 50 kb**.
10. **Faz**: prediz CDSs, tRNAs e rRNAs e atribui funções por similaridade com bancos de dados. **Não faz**: não prediz promotores, terminadores nem sítios de splicing, e não é para eucariotos.
11. É um gene provável cuja função não foi identificada (sem correspondência confiável nos bancos). Não é erro: é comum, mesmo em bactérias bem estudadas.
12. Por causa dos éxons e íntrons (o gene é "picotado"), das grandes regiões não codificantes, das repetições e do splicing alternativo; exige modelos estatísticos e dados de RNA.

## Atividade assíncrona, glossário e referências

### Atividade assíncrona

Entregue um **resumo de uma página** (PDF ou Word) interpretando os resultados da prática, até a data combinada com o professor, antes da Aula 4. O objetivo não é provar que você rodou os programas, e sim mostrar que **entendeu o significado de cada resultado**.

O resumo deve conter:

1. **Montagem (SPAdes + Quast)**: número de contigs, tamanho total e N50, com uma frase interpretando a qualidade (a montagem está fragmentada? por quê? qual o papel da cobertura de \~19×?).
2. **Anotação (Prokka)**: total de CDS, rRNA e tRNA; proporção de *hypothetical protein*; três genes anotados, indicando nome, função e tamanho.
3. **Comparação com referência**: localize no **NCBI Genome** um genoma de *Staphylococcus aureus* anotado profissionalmente e compare tamanho e número de genes com os seus resultados. Explique as diferenças.
4. **Procarioto × eucarioto**: um parágrafo curto explicando por que o Prokka funciona para bactérias, mas não para o genoma humano.

### Glossário

| Termo | Significado |
| --- | --- |
| **Anotação** | Identificar onde estão os elementos do genoma e para que servem |
| **Base nitrogenada** | A, C, G, T (DNA) ou A, C, G, U (RNA): a "letra" |
| **Cobertura** | Número médio de vezes que cada letra do genoma foi lida |
| **Contig** | Trecho contínuo montado por sobreposição de reads |
| **CDS** | Região codificadora de proteína, do códon de início ao de parada |
| **Códon** | Trio de nucleotídeos que especifica um aminoácido ou sinal |
| **Éxon / Íntron** | Trecho mantido / removido no mRNA maduro (eucariotos) |
| **FASTA / FASTQ** | Formatos de sequências sem / com qualidade |
| **GFF** | Tabela com as posições das anotações |
| **Grafo de De Bruijn** | Rede de k-mers usada pelos montadores |
| **Histona** | Proteína em torno da qual o DNA eucarioto se enrola |
| **k-mer** | Palavra de k letras obtida cortando um read |
| **N50 / L50** | Medidas de continuidade de uma montagem |
| **Nucleoide** | Região da célula procarionte onde fica o DNA |
| **Operon** | Genes transcritos juntos sob um promotor |
| **ORF** | Fase aberta de leitura: trecho sem códon de parada interno |
| **Paired-end** | Leitura das duas pontas de um fragmento (R1 e R2) |
| **Plasmídeo** | Pequeno DNA circular extra, comum em bactérias |
| **Promotor** | Sequência onde a RNA polimerase inicia a transcrição |
| **Read** | Fragmento de DNA lido pelo sequenciador |
| **Scaffold** | Contigs ordenados, com lacunas entre eles |
| **Splicing** | Remoção dos íntrons e união dos éxons no mRNA |
| **Terminador** | Sinal que encerra a transcrição |

### Para ir além (opcional)

- Rode o **FastQC** nos reads antes da montagem.
- Importe o genoma de referência do mesmo registro do Zenodo e rode o **Quast com referência** para comparar.
- Experimente o **BPROM** e o **ARNold** com os 150 pb a montante de um gene anotado pelo Prokka.
- Visualize o grafo da montagem no programa **Bandage**.

### Referências e links

- Galaxy: usegalaxy.org · Galaxy Training Network (training.galaxyproject.org), tutoriais de montagem e anotação.
- Dados da prática: Zenodo, registro 582600.
- SPAdes: Bankevich et al., *J Comput Biol*, 2012. Quast: Gurevich et al., *Bioinformatics*, 2013. Prokka: Seemann, *Bioinformatics*, 2014.
- NCBI Genome: ncbi.nlm.nih.gov/genome · Ensembl: ensembl.org.
- BPROM (Softberry) · ARNold (serviço web do I2BC, França).
- Material de aula do Prof. Roberto Filho sobre material genético em procariotos e eucariotos (2016).
