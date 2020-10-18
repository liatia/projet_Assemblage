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
import statistics
from operator import itemgetter
import random
from random import randint
import networkx as nx
import matplotlib.pyplot as plt
random.seed(9001)

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
    """ Lecture du fichier fastq en entree
        Retourne : les sequences d'acides nucleiques
    """
    with open(fastq_file) as filin:
        for line in enumerate(filin):
            yield next(filin)[:-1]
            next(filin)
            next(filin)


def cut_kmer(seq, kmer_size):
    """ Coupe les sequences en kmers
        Retourne : tous les kmers trouves dans les sequences
    """
    for i in range(len(seq) - kmer_size + 1):
        yield seq[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """ Contruction d'un dict de kmers
        Retourne : un dictionnaire avec cle : kmer et valeur : le nombre
        d'occurence du kmer
    """
    kmer_dict = {}
    kmer_list = []
    seq_list = read_fastq(fastq_file)
    for seq in seq_list:
        kmer_list = kmer_list + (list(cut_kmer(seq, kmer_size)))
    for kmer in set(kmer_list):
        kmer_dict[kmer] = kmer_list.count(kmer)
    return kmer_dict


def build_graph(kmer_dict):
    """Creation d'un graph a partir des kmers
    Les noeuds correspondent aux prefix et au suffix des kmers
        Retourne : un graph
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph.add_node(prefix)
        graph.add_node(suffix)
        graph.add_edge(prefix, suffix, weight = kmer_dict[kmer])
    return graph


def get_starting_nodes(graph):
    """Obtention de noeuds d'entree sans predecesseurs
        Retourne : une liste de noeuds d'entree
    """
    nodes_list_in = []
    for node in graph.nodes:
        if list(graph.predecessors(node)) == []:
            nodes_list_in.append(node)
    return nodes_list_in


def get_sink_nodes(graph):
    """Obtention de noeuds de sortie sans successeurs
        Retourne : une liste de noeuds de sortie
    """
    nodes_list_out = []
    for node in graph.nodes:
        if list(graph.successors(node)) == []:
            nodes_list_out.append(node)
    return nodes_list_out


def get_contigs(graph, nodes_list_in, nodes_list_out):
    """Creation d'une liste de contigs en concatenant les kmers d'un chemin
        Retourne d'une lite de tuple(contig, longueur du contig)
    """
    list_contigs = []
    for node_in in nodes_list_in:
        for node_out in nodes_list_out:
            for paths in nx.all_simple_paths(graph, source=node_in, target=node_out):
                contig = paths[0]
                for path in paths[1:]:
                    contig = contig + path[-1]
                list_contigs.append((contig, len(contig)))
    return list_contigs


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def save_contigs(list_contigs, output_file):
    """ Ecrit un fichier de sortie contenant les contigs selon le format fasta
    """
    with open(output_file, "w") as filout:
        for i, contig in enumerate(list_contigs):
            filout.write(">contig_" + str(i) + " len=" + str(contig[1]) + "\n")
            filout.write(fill(contig[0])+"\n")


def std(val_list):
    """ Calcul de la variance
    """
    return statistics.stdev(val_list)


def path_average_weight(graph, path):
    """ Calcul d'un poids moyen dans un graph
    """
    weight_list = []
    for i in range(len(path)-1):
        weight_list.append(graph.get_edge_data(path[i],
            path[i+1])['weight'])
    return statistics.mean(weight_list)


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """ Detection des chemins indesirables et des noeuds a supprimer
    """
    for path in path_list:
        for i in range(len(path)):
            if (i == 0 and delete_entry_node == False) or \
            (i == len(path)-1 and delete_sink_node == False):
                continue
            graph.remove_node(path[i])
    return graph


def select_best_path(graph, path_list, path_size_list, average_weight_list,
    delete_entry_node = False, delete_sink_node = False):
    """ Selection du meilleur chemin
    """
    if average_weight_list[0] > average_weight_list[1]:
        remove_paths(graph, [path_list[1]], delete_entry_node, delete_sink_node)
    elif average_weight_list[0] < average_weight_list[1]:
        remove_paths(graph, [path_list[0]], delete_entry_node, delete_sink_node)
    else:
        if path_size_list[0] > path_size_list[1]:
            remove_paths(graph, [path_list[1]], delete_entry_node, delete_sink_node)
        elif path_size_list[0] < path_size_list[1]:
            remove_paths(graph, [path_list[0]], delete_entry_node, delete_sink_node)
        else:
            remove_paths(graph, [path_list[random.randint(0,1)]], delete_entry_node,
                delete_sink_node)
    return graph


def solve_bubble(graph, ancetre, descendant):
    paths = list(nx.all_simple_paths(graph, source=ancetre, target=descendant))
    while len(paths) >= 2:
        path_size_list = [len(paths[0]), len(paths[1])]
        average_weight_list = [path_average_weight(graph, paths[0]),
                               path_average_weight(graph, paths[1])]
        graph = select_best_path(graph, paths, path_size_list, average_weight_list)
        paths = list(nx.all_simple_paths(graph, ancetre, descendant))
    return graph


def simplify_bubbles(graph):
    flag_anc = 0
    flag_des = 0
    ancetre_list = []
    descendant_list = []
    nodes_list_in = get_starting_nodes(graph)
    nodes_list_out = get_sink_nodes(graph)
    for node_in in nodes_list_in:
        for node_out in nodes_list_out:
            paths = list(nx.all_simple_paths(graph, source=node_in, target=node_out))
            while len(paths) > 1:
                for i in range(len(paths[0])):
                    if paths[0][i] != paths[1][i] and flag_anc == 0:
                        ancetre_list.append(paths[0][i-1])
                        flag_anc = 1
                    elif paths[0][i] == paths[1][i] and flag_anc == 1:
                        flag_anc = 0
                    if paths[0][len(paths[0])-1-i] != paths[1][len(paths[1])-1-i] and flag_des == 0:
                        descendant_list.append(paths[0][-i])
                        flag_des = 1
                    elif paths[0][len(paths[0])-1-i] == paths[1][len(paths[1])-1-i] and flag_des == 1:
                        flag_des = 0
                if len(ancetre_list) == 1 and len(descendant_list) == 1:
                    graph = solve_bubble(graph, ancetre_list[0], descendant_list[0])
                else :
                    for j in range(len(ancetre_list)):
                        graph = solve_bubble(graph, ancetre_list[j], descendant_list[len(descendant_list)-1-j])
                paths = list(nx.all_simple_paths(graph, source=node_in, target=node_out))
    return graph


def solve_entry_tips(graph, nodes_list_in):
    entry_list = []
    nb_successors = -1
    actual_node = list(graph.successors(nodes_list_in[0]))[0]
    while nb_successors != 0:
        nb_predecesseurs = len(list(graph.predecessors(actual_node)))
        if nb_predecesseurs > 1:
            entry_list.append(actual_node)
        actual_node = list(graph.successors(actual_node))[0]
        nb_successors = len(list(graph.successors(actual_node)))
    last_entry_node = entry_list[-1]
    print(last_entry_node)
    while len(nodes_list_in) != 1:
        path0 = list(nx.all_simple_paths(graph, source=nodes_list_in[0], target=last_entry_node))[0]
        path1 = list(nx.all_simple_paths(graph, source=nodes_list_in[1], target=last_entry_node))[0]
        path_list = [path0, path1]
        path_size_list = [len(path0), len(path1)]
        average_weight_list = [path_average_weight(graph, path0), path_average_weight(graph, path1)]
        graph = select_best_path(graph, path_list, path_size_list, average_weight_list, True, False)
        nodes_list_in = get_starting_nodes(graph)
    return graph


def solve_out_tips(graph, nodes_list_out):
    # nx.draw(graph, with_labels = True)
    # plt.show()
    out_list = []
    nb_predecesseurs = -1
    actual_node = list(graph.predecessors(nodes_list_out[0]))[0]
    while nb_predecesseurs != 0:
        nb_successors = len(list(graph.successors(actual_node)))
        if nb_successors > 1:
            out_list.append(actual_node)
        actual_node = list(graph.predecessors(actual_node))[0]
        nb_predecesseurs = len(list(graph.predecessors(actual_node)))
    last_out_node = out_list[-1]
    while len(nodes_list_out) != 1:
        path0 = list(nx.all_simple_paths(graph, source=last_out_node, target=nodes_list_out[0]))[0]
        path1 = list(nx.all_simple_paths(graph, source=last_out_node, target=nodes_list_out[1]))[0]
        path_list = [path0, path1]
        path_size_list = [len(path0), len(path1)]
        average_weight_list = [path_average_weight(graph, path0), path_average_weight(graph, path1)]
        graph = select_best_path(graph, path_list, path_size_list, average_weight_list, False, True)
        nodes_list_out = get_starting_nodes(graph)
    return graph

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
    # 3.b. Detection des pointes


if __name__ == '__main__':
    main()