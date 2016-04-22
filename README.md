# Sparse Bilingual Word Representations 

Code for Sparse Bilingual Embeddings as described in [Sparse Bilingual Word Representations for Cross-lingual Lexical Entailment](http://cs.umd.edu/~yogarshi/publications/2016/naacl2016.pdf).

## Prerequisites

* MATLAB

## Getting Embeddings

Run 'sh fasta\_biling.sh' with the following parameters (in order):

* En Vocab File : One word per line (|e| lines)
* Fr Vocab File : One word per line (|f| lines)
* Dense En embeddings : One vector per line, each vector a space seperated list of floats (|e| lines)
* Dense Fr embeddings : One vector per line, each vector a space seperated list of floats (|f| lines)
* Alignment matrix : .mat file containing the crosslingual statistics matrix S (of size |e| x |f|)

Example files are available [here](http://www.umiacs.umd.edu/~yogarshi/data/bisparse/sample_data.tar.gz).

There are other hyperparameters in the script which you should consider adjusting.

## Data

The data folder contains

* final\_dataset.tsv - The French-English crosslingual lexical entailment dataset
* bisparse\_{en,fr}.txt - The French-English bilingual sparse vectors used to obtained results in the paper

## Utils

This folder contains Some other useful code :

* top\_dims.py - Interpret the dimensions given a (sparse) vector file

If you use this code or the associated dataset, please cite the paper!

	@InProceedings{VyasCarpuat2015,
        	Title = {Sparse Bilingual Word Representations for Cross-lingual Lexical Entailment},
        	Booktitle = {Proceedings of NAACL},
        	Author = {Vyas, Yogarshi and Carpuat, Marine},
        	Year = {2016},
        	Location = {San Diego, United States of America}
   	}

