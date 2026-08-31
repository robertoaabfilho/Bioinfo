#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
============================================================================
 PIPELINE DIDÁTICO: Alinhamento (MUSCLE) -> Consenso -> BLAST (blastn/blastx/tblastx)
 Disciplina: Introdução à Bioinformática - PPGED/IME
============================================================================

O QUE ESTE SCRIPT FAZ (visão geral, sem jargão):

  1) Lê um arquivo com várias sequências de DNA no formato FASTA
     (o formato de texto mais usado em bioinformática: cada sequência
     começa com uma linha ">nome" seguida da sequência propriamente dita).

  2) Alinha essas sequências com o programa MUSCLE. "Alinhar" significa
     colocar as sequências uma embaixo da outra, inserindo "buracos"
     (representados por "-") de forma que posições homólogas (que vieram
     de um ancestral comum) fiquem na mesma coluna. É como comparar
     várias frases parecidas e alinhar as palavras equivalentes.

  3) A partir do alinhamento, calcula uma "sequência consenso": em cada
     coluna do alinhamento, escolhe a base mais frequente. O resultado é
     uma sequência "representativa" de todo o conjunto - útil, por
     exemplo, quando você tem várias leituras/isolados de um mesmo gene
     e quer uma sequência única para representá-los.

  4) Usa essa sequência consenso para pesquisar bancos de dados públicos
     do NCBI com três variantes do BLAST (explicadas mais abaixo):
     blastn, blastx e tblastx.

PRÉ-REQUISITOS (fazer UMA VEZ, antes de rodar o script):

  a) Instalar o Python (>=3.8) e o pacote Biopython:
         pip install biopython

  b) Instalar o MUSCLE (programa de alinhamento) e deixá-lo acessível
     pelo terminal (na "PATH" do sistema). Duas opções:
       - Baixar o executável em https://drive5.com/muscle/ (versão 5,
         a mais atual); ou
       - No Linux/Conda:  conda install -c bioconda muscle
     Este script assume a sintaxe do MUSCLE versão 5. Se você tiver a
     versão 3.8 (mais antiga), veja a observação na função `rodar_muscle`.

  c) Ter conexão com a internet: as buscas BLAST deste script são feitas
     online, nos servidores do NCBI (não é necessário instalar nenhum
     banco de dados localmente).

  d) Ter um arquivo FASTA de entrada com as sequências de DNA que você
     quer alinhar (ex.: "minhas_sequencias.fasta").

COMO EXECUTAR:

    python pipeline_alinhamento_blast.py minhas_sequencias.fasta

============================================================================
"""

import sys
import os
import subprocess
import urllib.error

from collections import Counter

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML


# ----------------------------------------------------------------------
# ETAPA 1: RODAR O MUSCLE PARA ALINHAR AS SEQUÊNCIAS
# ----------------------------------------------------------------------
def rodar_muscle(arquivo_entrada, arquivo_alinhado):
    """
    Chama o programa MUSCLE por fora do Python (via linha de comando)
    para alinhar as sequências do arquivo_entrada e salvar o resultado
    em arquivo_alinhado (também em formato FASTA, mas agora "alinhado",
    ou seja, com os símbolos "-" indicando os buracos).

    O Biopython NÃO faz o alinhamento sozinho: ele apenas "conversa" com
    o programa MUSCLE, que é quem executa o algoritmo de alinhamento.
    Por isso é preciso ter o MUSCLE instalado separadamente.
    """
    print(f"[1/4] Rodando MUSCLE para alinhar as sequências de '{arquivo_entrada}'...")

    # Sintaxe do MUSCLE versão 5 (a mais recente):
    #   muscle -align entrada.fasta -output saida.fasta
    #
    # Se você usa a versão antiga (3.8), troque a linha abaixo por:
    #   comando = ["muscle", "-in", arquivo_entrada, "-out", arquivo_alinhado]
    comando = ["muscle", "-align", arquivo_entrada, "-output", arquivo_alinhado]

    try:
        subprocess.run(comando, check=True, capture_output=True, text=True)
    except FileNotFoundError:
        sys.exit(
            "ERRO: o programa 'muscle' não foi encontrado no sistema.\n"
            "Instale o MUSCLE e certifique-se de que ele está na PATH "
            "(veja as instruções no início deste script)."
        )
    except subprocess.CalledProcessError as erro:
        sys.exit(f"ERRO ao rodar o MUSCLE:\n{erro.stderr}")

    print(f"      Alinhamento salvo em '{arquivo_alinhado}'.")


# ----------------------------------------------------------------------
# ETAPA 2: CALCULAR A SEQUÊNCIA CONSENSO A PARTIR DO ALINHAMENTO
# ----------------------------------------------------------------------
def gerar_consenso(arquivo_alinhado, limiar=0.5):
    """
    Lê o alinhamento múltiplo (arquivo_alinhado) e devolve a sequência
    consenso: para cada coluna do alinhamento, escolhe a base (A, T, C,
    G) mais frequente. Se a base mais comum não atingir o 'limiar' de
    frequência definido (padrão: 50%), a posição recebe um "N" (base
    indefinida). Colunas onde a maioria é um "buraco" (gap, "-") são
    ignoradas, ou seja, não entram na sequência consenso final.

    Pense nisso como "tirar a média" das sequências alinhadas, coluna
    por coluna.

    OBSERVAÇÃO TÉCNICA: versões mais novas do Biopython descontinuaram
    o método pronto `AlignInfo.SummaryInfo.dumb_consensus`. Por isso,
    calculamos o consenso manualmente aqui, usando apenas o indexador
    `alinhamento[:, i]` (que devolve todos os caracteres da coluna i)
    e a classe `Counter` do próprio Python para contar as ocorrências.
    Isso deixa o script menos dependente de detalhes internos do
    Biopython que podem mudar entre versões.
    """
    print("[2/4] Calculando a sequência consenso...")

    alinhamento = AlignIO.read(arquivo_alinhado, "fasta")
    num_sequencias = len(alinhamento)
    comprimento = alinhamento.get_alignment_length()

    bases_consenso = []
    for i in range(comprimento):
        # alinhamento[:, i] devolve uma string com o caractere de CADA
        # sequência naquela coluna i do alinhamento (ex.: "AAT-A").
        coluna = alinhamento[:, i]
        contagem = Counter(coluna)
        caractere_mais_comum, vezes = contagem.most_common(1)[0]
        frequencia = vezes / num_sequencias

        if frequencia < limiar:
            bases_consenso.append("N")          # sem consenso claro
        elif caractere_mais_comum == "-":
            continue                             # coluna é majoritariamente gap
        else:
            bases_consenso.append(caractere_mais_comum)

    consenso = "".join(bases_consenso)
    print(f"      Consenso gerado com {len(consenso)} posições (bases).")
    return consenso


# ----------------------------------------------------------------------
# ETAPA 3: SALVAR A SEQUÊNCIA CONSENSO EM UM ARQUIVO FASTA
# ----------------------------------------------------------------------
def salvar_fasta(sequencia_texto, nome_id, arquivo_saida):
    """
    Empacota a sequência consenso (que por enquanto é só uma string de
    texto) em um objeto SeqRecord do Biopython e grava em um arquivo
    FASTA, para que possa ser reaproveitada por outros programas
    (incluindo o próprio BLAST).
    """
    registro = SeqRecord(Seq(sequencia_texto), id=nome_id,
                          description="sequencia consenso gerada por MUSCLE")
    SeqIO.write(registro, arquivo_saida, "fasta")
    print(f"      Consenso salvo em '{arquivo_saida}'.")


def salvar_consenso_txt(sequencia_texto, arquivo_saida, largura_linha=60):
    """
    Salva a sequência consenso em um arquivo .txt "simples", formatada
    em blocos de 'largura_linha' caracteres (o mesmo padrão de quebra
    de linha usado em arquivos FASTA), junto com um pequeno resumo
    (comprimento total e quantas posições ficaram como "N", ou seja,
    sem consenso claro entre as sequências).
    """
    quantidade_n = sequencia_texto.count("N")

    with open(arquivo_saida, "w") as arquivo:
        arquivo.write("SEQUÊNCIA CONSENSO\n")
        arquivo.write("==================\n")
        arquivo.write(f"Comprimento total: {len(sequencia_texto)} posições\n")
        arquivo.write(
            f"Posições sem consenso claro (N): {quantidade_n} "
            f"({100 * quantidade_n / len(sequencia_texto):.1f}%)\n\n"
        )
        for inicio in range(0, len(sequencia_texto), largura_linha):
            arquivo.write(sequencia_texto[inicio: inicio + largura_linha] + "\n")

    print(f"      Consenso (texto legível) salvo em '{arquivo_saida}'.")


# ----------------------------------------------------------------------
# ETAPA 4: RODAR OS TRÊS TIPOS DE BLAST NA SEQUÊNCIA CONSENSO
# ----------------------------------------------------------------------
def rodar_blast(sequencia_consenso, programa, banco_dados, arquivo_saida_xml):
    """
    Envia a sequência consenso para os servidores do NCBI e roda uma
    busca BLAST online.

    O QUE É O BLAST?
    -----------------
    BLAST (Basic Local Alignment Search Tool) é a ferramenta mais usada
    em bioinformática para responder à pergunta: "essa sequência que eu
    tenho é parecida com alguma sequência já conhecida (armazenada em
    bancos de dados públicos)?". Ele varre bancos como o "nt" (sequências
    de nucleotídeos) ou o "nr" (proteínas não-redundantes) procurando
    trechos semelhantes à sua sequência de consulta ("query").

    Cada variante do BLAST muda o TIPO da sequência de consulta e do
    banco de dados contra o qual ela é comparada:

      - blastn  : consulta = DNA/RNA   | banco = DNA/RNA
                  (compara nucleotídeo com nucleotídeo, "letra com letra").

      - blastx  : consulta = DNA/RNA   | banco = PROTEÍNA
                  Antes de comparar, o programa TRADUZ sua sequência de
                  DNA nos 6 possíveis "quadros de leitura" (as 3 formas
                  de ler o DNA na direção direta + as 3 formas na direção
                  reversa-complementar) e compara cada tradução com
                  proteínas conhecidas. Muito usado para descobrir qual
                  proteína um trecho de DNA (por exemplo, um gene ainda
                  não caracterizado) provavelmente codifica.

      - tblastx : consulta = DNA/RNA   | banco = DNA/RNA
                  (mas comparado como PROTEÍNA)
                  Aqui, tanto a sua sequência quanto TODAS as sequências
                  do banco de dados são traduzidas nos 6 quadros de
                  leitura, e a comparação é feita entre as proteínas
                  resultantes. É o mais lento dos três (traduz dos dois
                  lados), mas consegue detectar semelhanças "escondidas"
                  que só aparecem no nível da proteína, mesmo quando o
                  DNA já mudou bastante (porque o código genético é
                  degenerado: vários "trincas" de DNA podem virar o
                  mesmo aminoácido).

    Parâmetros importantes:
      - programa: "blastn", "blastx" ou "tblastx"
      - banco_dados: geralmente "nt" (nucleotídeos) para blastn e
        tblastx, e "nr" (proteínas) para blastx.

    IMPORTANTE: buscas online no NCBI podem levar de alguns segundos a
    vários minutos, dependendo do tamanho da fila de outros usuários e
    da estabilidade da conexão. Buscas com o banco "nr" (proteínas),
    usado pelo blastx, costumam ser as mais demoradas.

    RESILIÊNCIA A FALHAS DE REDE: uma busca BLAST online pode levar
    minutos para responder, e nesse meio tempo a conexão pode cair
    (erro de SSL, timeout, etc.) - isso não é um bug do script, é uma
    instabilidade da rede/servidor. Por isso, esta função tenta a busca
    até 'tentativas' vezes antes de desistir; se todas falharem, ela
    avisa e devolve False (em vez de travar o programa inteiro),
    permitindo que o pipeline continue para o próximo BLAST.
    """
    print(f"[3/4] Rodando {programa} (banco de dados: {banco_dados})... "
          f"isso pode demorar alguns minutos.")

    tentativas = 3
    for tentativa in range(1, tentativas + 1):
        try:
            resultado_handle = NCBIWWW.qblast(
                program=programa,
                database=banco_dados,
                sequence=sequencia_consenso,
            )

            # Salva o resultado bruto (formato XML) em disco, para não
            # perder a busca caso algo dê errado na etapa de leitura.
            with open(arquivo_saida_xml, "w") as arquivo:
                arquivo.write(resultado_handle.read())

            print(f"      Resultado do {programa} salvo em '{arquivo_saida_xml}'.")
            return True

        except (urllib.error.URLError, ConnectionError, TimeoutError) as erro:
            print(f"      [Aviso] Falha de conexão na tentativa {tentativa}/"
                  f"{tentativas} para {programa}: {erro}")
            if tentativa == tentativas:
                print(f"      [Aviso] Desistindo de {programa} após "
                      f"{tentativas} tentativas. Pulando para a próxima busca.")
                return False


def mostrar_melhores_hits(arquivo_xml, programa, arquivo_saida_txt,
                           n_hits=5, e_value_maximo=1e-5):
    """
    Lê o arquivo XML gerado pelo BLAST, imprime no terminal e também
    grava em um arquivo .txt um resumo com os 'n_hits' melhores
    resultados (os mais parecidos com a nossa sequência consenso).

    Dois conceitos-chave para interpretar os resultados:

      - E-value (valor esperado): estima quantos "achados" tão bons
        quanto aquele apareceriam por puro acaso, num banco de dados
        daquele tamanho. Quanto MENOR o E-value, mais confiável é a
        semelhança encontrada (E-value perto de 0 = semelhança muito
        provavelmente real e não coincidência).

      - Identidade (% identity): de todas as posições comparadas entre
        a sua sequência e a sequência do banco, qual fração é
        exatamente igual.
    """
    linhas = [f"===== Principais resultados do {programa} =====\n"]

    with open(arquivo_xml) as arquivo:
        registros = NCBIXML.parse(arquivo)
        for registro in registros:
            if not registro.alignments:
                linhas.append("Nenhum resultado significativo encontrado.\n")
                continue

            contador = 0
            for alinhamento in registro.alignments:
                for hsp in alinhamento.hsps:  # HSP = High-scoring Segment Pair
                    if hsp.expect <= e_value_maximo:
                        identidade_pct = 100 * hsp.identities / hsp.align_length
                        linhas.append(f"- {alinhamento.title[:80]}")
                        linhas.append(
                            f"    E-value: {hsp.expect:.2e} | "
                            f"Identidade: {identidade_pct:.1f}% | "
                            f"Comprimento do alinhamento: {hsp.align_length}"
                        )
                        contador += 1
                        break
                if contador >= n_hits:
                    break

    texto_resultado = "\n".join(linhas)
    print("\n" + texto_resultado)

    with open(arquivo_saida_txt, "w") as arquivo:
        arquivo.write(texto_resultado + "\n")
    print(f"      Resumo do {programa} salvo em '{arquivo_saida_txt}'.")


# ----------------------------------------------------------------------
# PROGRAMA PRINCIPAL: liga todas as etapas em sequência
# ----------------------------------------------------------------------
def main():
    if len(sys.argv) != 2:
        sys.exit("Uso: python pipeline_alinhamento_blast.py <arquivo_multifasta.fasta>")

    arquivo_entrada = sys.argv[1]
    if not os.path.exists(arquivo_entrada):
        sys.exit(f"ERRO: arquivo '{arquivo_entrada}' não encontrado.")

    # Nomes dos arquivos que serão gerados ao longo do pipeline
    arquivo_alinhado = "alinhamento_muscle.fasta"
    arquivo_consenso_fasta = "sequencia_consenso.fasta"
    arquivo_consenso_txt = "sequencia_consenso.txt"

    # --- Etapa 1: alinhamento com MUSCLE ---
    rodar_muscle(arquivo_entrada, arquivo_alinhado)

    # --- Etapa 2: geração da sequência consenso ---
    sequencia_consenso = gerar_consenso(arquivo_alinhado)

    # --- Etapa 3: salvar o consenso em FASTA e em TXT ---
    salvar_fasta(sequencia_consenso, "consenso", arquivo_consenso_fasta)
    salvar_consenso_txt(sequencia_consenso, arquivo_consenso_txt)

    # --- Etapa 4: rodar os três tipos de BLAST ---
    buscas = [
        # (programa, banco_de_dados, arquivo_xml, arquivo_txt_resumo)
        ("blastn", "nt", "resultado_blastn.xml", "resultado_blastn.txt"),
        ("blastx", "nr", "resultado_blastx.xml", "resultado_blastx.txt"),
        ("tblastx", "nt", "resultado_tblastx.xml", "resultado_tblastx.txt"),
    ]

    buscas_ok = []
    buscas_com_falha = []

    for programa, banco, arquivo_xml, arquivo_txt in buscas:
        sucesso = rodar_blast(sequencia_consenso, programa, banco, arquivo_xml)
        if sucesso:
            mostrar_melhores_hits(arquivo_xml, programa, arquivo_txt)
            buscas_ok.append(programa)
        else:
            buscas_com_falha.append(programa)

    print("\n[4/4] Pipeline concluído! Arquivos gerados:")
    print(f"  - Alinhamento MUSCLE  : {arquivo_alinhado}")
    print(f"  - Consenso (FASTA)    : {arquivo_consenso_fasta}")
    print(f"  - Consenso (TXT)      : {arquivo_consenso_txt}")
    for programa in buscas_ok:
        print(f"  - Resultado {programa:<8}: resultado_{programa}.xml "
              f"e resultado_{programa}.txt")
    if buscas_com_falha:
        print(f"\n  [Atenção] Não foi possível concluir: "
              f"{', '.join(buscas_com_falha)}. "
              f"Rode o script novamente mais tarde para tentar essas buscas "
              f"(o consenso e as buscas que já funcionaram estão salvos).")


if __name__ == "__main__":
    main()
