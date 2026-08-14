# Ementa – Introdução à Bioinformática

**Programa de Pós-Graduação em Engenharia de Defesa (PPGED) – Instituto Militar de Engenharia (IME)**

## 1. Identificação da Disciplina

| Campo | Descrição |
| --- | --- |
| **Disciplina** | Introdução à Bioinformática |
| **Carga horária** | 30 horas |
| **Modalidade** | Online (síncrona) |
| **Docentes** | Prof. Mariana Quezado (IME) e [seu nome] |
| **Público-alvo** | Discentes do Programa de Pós-Graduação em Engenharia de Defesa (PPGED/IME) |
| **Pré-requisitos** | Conhecimento básico das linguagens Python e R; noções fundamentais de biologia molecular (estrutura de ácidos nucleicos e proteínas, dogma central) e de biotecnologia |

## 2. Objetivo Geral

Capacitar o aluno de pós-graduação ao uso crítico e prático das principais ferramentas de bioinformática, com ênfase em análise de sequências, proteômica, epigenômica, bioinformática estrutural e integração de dados ômicos, aplicados à pesquisa em química de proteínas e macromoléculas.

## 3. Objetivos Específicos

Ao final da disciplina, o aluno deverá ser capaz de:

- Operar ambientes computacionais típicos de bioinformática (linha de comando, gerenciadores de pacotes, Python/R, plataformas em nuvem);
- Avaliar e processar dados brutos de sequenciamento de nova geração (NGS);
- Realizar buscas e comparações de sequências biológicas utilizando ferramentas de alinhamento;
- Interpretar resultados de análises de genomas e transcriptomas;
- Processar e interpretar dados de proteômica obtidos por espectrometria de massas, identificando e quantificando proteínas de interesse;
- Interpretar dados epigenômicos (metilação do DNA, modificações de histonas e acessibilidade da cromatina) e sua relação com a regulação da expressão gênica;
- Aplicar técnicas de bioinformática estrutural (visualização molecular, modelagem por homologia, docking e dinâmica molecular);
- Construir e interpretar árvores filogenéticas no contexto da evolução molecular;
- Utilizar ferramentas de mineração de texto e bases de vias metabólicas para integração de dados biológicos;
- Apresentar e defender oralmente, perante banca docente, os resultados de um trabalho de pesquisa que integre ferramentas de bioinformática;
- Relacionar as ferramentas estudadas a problemas de química de proteínas e macromoléculas de interesse em pesquisa aplicada.

## 4. Ementa (Conteúdo Programático)

| Unidade | Tópico | Carga horária |
| --- | --- | --- |
| 1 | Introdução ao ambiente computacional para bioinformática e mineração de textos com dados biológicos (ênfase em biotecnologia) | 3h |
| 2 | Qualidade e processamento de dados de NGS; identificação e comparação de sequências biológicas | 3h |
| 3 | Genômica: montagem e anotação de genomas | 3h |
| 4 | Transcriptômica: RNA-seq e análise de expressão diferencial | 3h |
| 5 | Proteômica: identificação e quantificação de proteínas por espectrometria de massas | 3h |
| 6 | Epigenômica: metilação do DNA, modificações de histonas e acessibilidade da cromatina | 3h |
| 7 | Bioinformática estrutural: visualização, modelagem por homologia, docking e dinâmica molecular | 3h |
| 8 | Análise filogenética e evolução molecular | 3h |
| 9 | Integração de vias metabólicas; aplicações em química de proteínas e macromoléculas | 3h |
| 10 | Defesa dos artigos e encerramento da disciplina | 3h |
| | **Total** | **30h** |

## 5. Metodologia

Aulas expositivo-dialogadas síncronas, com demonstrações práticas guiadas ("hands-on") e estudos de caso baseados em dados públicos. Ferramentas e bases de dados utilizadas em cada unidade:

| Unidade | Ferramentas / Bases de dados |
| --- | --- |
| **1. Ambiente computacional e mineração de texto** | Terminal Linux/Shell; Anaconda/Conda; Jupyter Notebook; Python (NumPy, Pandas, Biopython) e R (Bioconductor); Google Colab; repositórios NCBI, ENA e Ensembl; PubMed/PMC e ferramentas de mineração de textos com dados biológicos, com ênfase em biotecnologia (ex.: PubTator) |
| **2. Qualidade/processamento NGS e comparação de sequências** | FastQC, MultiQC (controle de qualidade); Trimmomatic/Cutadapt (trimming); BWA e Bowtie2 (alinhamento a referência); SAMtools; dados de exemplo do SRA (Sequence Read Archive); BLAST (blastn, blastp, blastx) e BLAST online (NCBI); Clustal Omega e MUSCLE (alinhamento múltiplo); bancos GenBank, UniProt, RefSeq |
| **3. Genômica** | SPAdes (montagem de genomas); Prokka (anotação de genomas); bases NCBI Genome e Ensembl |
| **4. Transcriptômica** | Trinity (montagem de transcriptoma); HISAT2/STAR (mapeamento de RNA-seq); DESeq2/edgeR em R (expressão diferencial); base Ensembl |
| **5. Proteômica** | Espectrometria de massas (fundamentos e fluxo de análise "bottom-up"); MaxQuant e mecanismo de busca Andromeda (identificação e quantificação de peptídeos/proteínas); Perseus (tratamento estatístico e visualização); PatternLab for Proteomics (ferramenta nacional, Fiocruz/UFRJ); Skyline (proteômica direcionada); bancos UniProt e ProteomeXchange/PRIDE (repositório de dados brutos de espectrometria de massas) |
| **6. Epigenômica** | Bismark e methylKit (análise de metilação do DNA por sequenciamento com bissulfito); MACS2 (chamada de picos em ChIP-seq para modificações de histonas e fatores de transcrição); conceitos de ATAC-seq (acessibilidade da cromatina); IGV e UCSC Genome Browser (visualização de tracks epigenômicos); bases ENCODE, Roadmap Epigenomics e GEO |
| **7. Bioinformática estrutural** | PyMOL e ChimeraX (visualização); SWISS-MODEL e AlphaFold/ColabFold (modelagem por homologia); banco de dados PDB (Protein Data Bank); AutoDock Vina e AutoDock Tools (docking molecular); GROMACS (dinâmica molecular); base PubChem (ligantes) |
| **8. Filogenética e evolução molecular** | MEGA e/ou IQ-TREE (construção de árvores); Jalview (visualização de alinhamentos); bases GenBank e TimeTree |
| **9. Vias metabólicas e integração de dados** | KEGG e Reactome (vias metabólicas); STRING (redes de interação proteína-proteína) |
| **10. Defesa dos artigos e encerramento** | Apresentação oral dos artigos científicos finais (ferramentas de apresentação como PowerPoint/Google Slides, a critério do aluno); sessão de perguntas e respostas; discussão integrada dos temas da disciplina |

Atividades práticas assíncronas complementam os encontros síncronos, com exercícios computacionais aplicados às ferramentas de cada unidade.

## 6. Avaliação

- Lista de exercícios práticos aplicando as ferramentas apresentadas (análise de sequência, proteômica, epigenômica, estrutural e/ou filogenética);
- **Trabalho final, em formato de artigo científico**, integrando ao menos duas técnicas estudadas na disciplina, aplicadas a um problema relacionado à linha de pesquisa do aluno (ex.: análise estrutural de uma proteína de interesse combinada à busca de sequências homólogas, a dados de proteômica ou à análise filogenética). O artigo deverá seguir estrutura padrão (introdução, metodologia, resultados, discussão e conclusão) e formato compatível com submissão a periódico ou evento da área;
- Defesa oral do artigo final na unidade 10, com apresentação sintética do trabalho seguida de arguição e discussão.

## 7. Bibliografia

- Pevsner, J. *Bioinformatics and Functional Genomics*. Wiley-Blackwell.
- Compeau, P.; Pevzner, P. *Bioinformatics Algorithms: An Active Learning Approach*.
- Lesk, A. M. *Introduction to Bioinformatics*. Oxford University Press.
- Mariano, D. *Python para Bioinformática: Fundamentos de Programação para Bioinformática e Biologia Computacional*. Novatec Editora, 2025.
- Aebersold, R.; Mann, M. *Mass-spectrometry-based proteomics*. Nature (artigo de referência sobre fundamentos de proteômica).
- Allis, C. D.; Jenuwein, T.; Reinberg, D. (Eds.). *Epigenetics*. Cold Spring Harbor Laboratory Press (referência sobre fundamentos de epigenética e epigenômica).
- Documentação oficial de ferramentas: NCBI BLAST, PDB, UniProt, KEGG, AlphaFold/SWISS-MODEL, GROMACS/AutoDock, MaxQuant/Perseus, Bismark/MACS2.
- Artigos científicos complementares indicados pelos docentes conforme a evolução da turma.
