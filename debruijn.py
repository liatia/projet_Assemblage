#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Lara Herrmann"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Lara Herrmann"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lara Herrmann"
__email__ = "lara.herrma@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    with open(fastq_file) as filin:
        for line in enumerate(filin):
            yield next(filin)[:-1]
            next(filin)       
            next(filin)


def cut_kmer(seq, kmer_size):
    for i in range(len(seq) - kmer_size + 1):
        yield seq[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    kmer_dict = {}
    kmer_list = []
    seq_list = read_fastq(fastq_file)
    for seq in seq_list : 
        kmer_list = kmer_list + (list(cut_kmer(seq, kmer_size)))
    for kmer in set(kmer_list):
        kmer_dict[kmer] = kmer_list.count(kmer)
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph.add_node(prefix)
        graph.add_node(suffix)
        graph.add_edge(prefix, suffix, weight = kmer_dict[kmer])
    # nx.draw(graph, with_labels = True)
    # plt.show()
    return graph


def get_starting_nodes(graph):
    nodes_list_in = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []: 
            nodes_list_in.append(node)
    return nodes_list_in


def get_sink_nodes(graph):
    nodes_list_out = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            nodes_list_out.append(node)
    return nodes_list_out


def get_contigs(graph, nodes_list_in, nodes_list_out):
    list_contigs = []
    for node_in in nodes_list_in: 
        for node_out in nodes_list_out:
            for paths in nx.all_simple_paths(graph, source=node_in, target=node_out):
                contig = paths[0] 
                for path in paths[1:]:
                    contig = contig + path[-1]
                list_contigs.append((contig, len(contig)))
    print(list_contigs)
    return list_contigs


def fill(text, width=80):    
    """Split text with a line return to respect fasta format"""    
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(list_contigs, output_file):
    # écrit un fichier de sortie contenant les contigs selon le format fasta 
    # (retour chariot tous les 80 charactères) à l’aide de la fonction fill
    with open(output_file, "w") as filout:
        for i, contig in enumerate(list_contigs):
            filout.write(">" + str(i) + " len=" + str(contig[1]) + "\n")
            print(fill(contig[0]))
            filout.write(fill(contig[0])+"\n")


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # 1. Création du graphe de Bruijn
    # 1.a. Identification des k-mer uniques
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)

    # 1.b. Construction de l'arbre de Bruijn
    graph = build_graph(kmer_dict)

    # 2. Parcours du graphe de de Bruijn
    nodes_list_in = get_starting_nodes(graph)
    nodes_list_out = get_sink_nodes(graph)
    list_contigs = get_contigs(graph, nodes_list_in, nodes_list_out)
    save_contigs(list_contigs, args.output_file)

    # 3. Simplification du graphe de de Bruijn
    # 3.a. Résolution des bulles



if __name__ == '__main__':
    main()

