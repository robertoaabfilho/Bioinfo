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

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
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
    G) mais frequente. Se nenhuma base atingir o 'limiar' de frequência
    definido (padrão: 50%), a posição recebe um "N" (base indefinida).

    Pense nisso como "tirar a média" das sequências alinhadas, coluna
    por coluna.
    """
    print("[2/4] Calculando a sequência consenso...")

    alinhamento = AlignIO.read(arquivo_alinhado, "fasta")
    resumo = AlignInfo.SummaryInfo(alinhamento)

    # dumb_consensus: método simples de consenso "por maioria de votos"
    # disponível no Biopython. threshold = fração mínima de concordância
    # entre as sequências para aceitar uma base naquela coluna.
    consenso = resumo.dumb_consensus(threshold=limiar, ambiguous="N")

    print(f"      Consenso gerado com {len(consenso)} posições (bases).")
    return str(consenso)


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
    vários minutos, dependendo do tamanho da fila de outros usuários.
    """
    print(f"[3/4] Rodando {programa} (banco de dados: {banco_dados})... "
          f"isso pode demorar alguns minutos.")

    resultado_handle = NCBIWWW.qblast(
        program=programa,
        database=banco_dados,
        sequence=sequencia_consenso,
    )

    # Salva o resultado bruto (formato XML) em disco, para não perder
    # a busca caso algo dê errado na etapa de leitura/interpretação.
    with open(arquivo_saida_xml, "w") as arquivo:
        arquivo.write(resultado_handle.read())

    print(f"      Resultado do {programa} salvo em '{arquivo_saida_xml}'.")


def mostrar_melhores_hits(arquivo_xml, programa, n_hits=5, e_value_maximo=1e-5):
    """
    Lê o arquivo XML gerado pelo BLAST e imprime, de forma resumida, os
    'n_hits' melhores resultados (os mais parecidos com a nossa
    sequência consenso).

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
    print(f"\n===== Principais resultados do {programa} =====")

    with open(arquivo_xml) as arquivo:
        registros = NCBIXML.parse(arquivo)
        for registro in registros:
            if not registro.alignments:
                print("Nenhum resultado significativo encontrado.")
                continue

            contador = 0
            for alinhamento in registro.alignments:
                for hsp in alinhamento.hsps:  # HSP = High-scoring Segment Pair
                    if hsp.expect <= e_value_maximo:
                        identidade_pct = 100 * hsp.identities / hsp.align_length
                        print(f"- {alinhamento.title[:80]}")
                        print(f"    E-value: {hsp.expect:.2e} | "
                              f"Identidade: {identidade_pct:.1f}% | "
                              f"Comprimento do alinhamento: {hsp.align_length}")
                        contador += 1
                        break
                if contador >= n_hits:
                    break


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
    arquivo_consenso = "sequencia_consenso.fasta"

    # --- Etapa 1: alinhamento com MUSCLE ---
    rodar_muscle(arquivo_entrada, arquivo_alinhado)

    # --- Etapa 2: geração da sequência consenso ---
    sequencia_consenso = gerar_consenso(arquivo_alinhado)

    # --- Etapa 3: salvar o consenso em FASTA ---
    salvar_fasta(sequencia_consenso, "consenso", arquivo_consenso)

    # --- Etapa 4: rodar os três tipos de BLAST ---
    buscas = [
        # (programa, banco_de_dados, arquivo_de_saida)
        ("blastn", "nt", "resultado_blastn.xml"),
        ("blastx", "nr", "resultado_blastx.xml"),
        ("tblastx", "nt", "resultado_tblastx.xml"),
    ]

    for programa, banco, arquivo_saida in buscas:
        rodar_blast(sequencia_consenso, programa, banco, arquivo_saida)
        mostrar_melhores_hits(arquivo_saida, programa)

    print("\n[4/4] Pipeline concluído! Arquivos gerados:")
    print(f"  - Alinhamento MUSCLE : {arquivo_alinhado}")
    print(f"  - Sequência consenso : {arquivo_consenso}")
    print(f"  - Resultados BLAST   : resultado_blastn.xml, "
          f"resultado_blastx.xml, resultado_tblastx.xml")


if __name__ == "__main__":
    main()
