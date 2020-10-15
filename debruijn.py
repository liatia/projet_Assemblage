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
    G = nx.DiGraph()
    for kmer in kmer_dict:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        G.add_nodes_from(prefix)
        G.add_nodes_from(suffix)
        G.add_edge(prefix, suffix, weight = kmer_dict[kmer])
    nx.draw(G, with_labels = True)
    plt.show()
    return G


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


if __name__ == '__main__':
    main()

