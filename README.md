# masterPATH
MasterPATH is an exploratory network analysis method that employs the shortest path approach and centrality measure to uncover members of active molecular pathways leading to the studied phenotype based on the results of functional genomics screening data.

## Databases

The method works with an integarted network that consists of interaction from the following databases:
Human Integrated Protein-Protein Interaction rEference database (HIPPIE):

http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/download.php

Human Protein Reference Database (HPRD):

http://hprd.org/download

SignaLink 2.0 database:

http://signalink.org/download

SIGnaling Network Open Resource database (SIGNOR):

http://signor.uniroma2.it/downloads.php

tFactS database

http://www.tfacts.org/TFactS-new/TFactS-v2/index1.html

TransmiR database:

http://www.cuilab.cn/transmir

miRTarBase database

http://mirtarbase.mbc.nctu.edu.tw/php/download.php


## Installation

The easiest way to use masterPATH is to download masterPATH.jar from JAR/ folder.

Otherwise the repo can be cloned into a new e.g. Netbeans IDE project and can be built.

The source files are available in the src/masterPATH folder.

Dependency libraries: -- commons-lang3-3.3.2

## Usage

First, Wrapper and Network classes should be imported:

    import masterPATH.Wrapper;
    import masterPATH.Network;
    
Next, create a Wrapper and a Network objects :
    
    Network nw ;
    Wrapper wr = new Wrapper();
    
Next, load network :

     nw = wr.load_network(String file_nodes, file_interactions);
     
where 
_String file_nodes_ is a full path to a file with network nodes,  
_String file_interactions_ is a full path to a file with network interactions.
The prebuilt network files are available in the Networks/ folder.


Finally, perform the computaion for mixed directed and undirected network:

    wr.find_shortest_paths_and_calculate_centrality(
            nw,
            file_hitlist,
            file_fimplementers,
            file_output,
            prefix,
            max_length_for_shortest_path,
            min_len_for_path,
            max_len_for_path,
            folder_for_random_paths,
            prefix_for_random_paths,
            number_of_permutations)


or for undirected network :
  
        wr.find_shortest_paths_and_calculate_centrality_ppi(
            nw,
            file_hitlist,
            file_fimplementers,
            file_output,
            prefix,
            max_length_for_shortest_path,
            min_len_for_path,
            max_len_for_path,
            folder_for_random_paths,
            prefix_for_random_paths,
            number_of_permutations)
            
where :

_Network nw_ is a network object, 

_String file_hitlist_ is a full path to a file with hit genes, 

_String file_fimplementers_ is a full path to a file with "final implemmenters",

_String file_output_ is a full path to an output file,

_String prefix_ is a prefix for the shortest paths ids,

_int max_length_for_shortest_path_ is maximum length for the breadth-first algorithm,

_int min_len_for_path_ is minimum length of the paths for which centrality will be calculated, 

_int max_len_for_path_ is maximim length of the paths for which centrality will be calculated,

_String folder_for_random_paths_ is a full path to folder where files for permutation analysis will be stored,

_String prefix_for_random_paths_ is a prefix for the shortest paths ids for permuted hit lists,

_int number_of_permutations_ is number of permutations.

## Output files

Output files from the method are :

__file_output + _paths_centrality__ file with paths information is a tab separated text file. File format : _path id_, _path as a list of interaction ids (each interaction is separated by a tab)_, _reserved field_, _centrality_, _hit gene-final implementer pairs that yield this path as a list of HGNC ids separated by semicolon_, _list of the shortest paths ids that yield this path separated by semicolon_, _reserved field_,  _reserved field_, _HGNC id of the node in the path with maximum centrality_, _official symbol of the node in the path with maximum centrality_, _hit gene-final implementer pairs that yield this path as a list of official symbols separated by semicolon_, _path as a list of intercator1-interactor2 pairs  (each interaction is separated by a tab)_, _centrality_, _p-value_.


__file_output + _nodes_centrality__ file with nodes information is a tab separated text file. File format : _node id_, _node ofiicial symbol_, _centrality_, _p-value_.
