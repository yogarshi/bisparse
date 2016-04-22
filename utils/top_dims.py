__author__ = 'yogarshi'


import argparse, sys
import math
import numpy as np

def load_vectors(vector_file):

	f = open(vector_file)

	#ig = f.readline()
	#ig = f.readline()
	words = []
	embeddings = []
	for each_line in f:
		x = each_line.strip().split()
		words.append(x[0])
		embeddings.append([float(y) for y in x[1:]])
	embeddings_array = np.array(embeddings)
	return words, embeddings_array


def get_top_words(words, embeddings_array, k):

	top_words = []
	embeddings_arrayT = embeddings_array.transpose()
	for each_dimension in embeddings_arrayT:
		#print each_dimension
		top10 = [words[x] for x in each_dimension.argsort()[-k:][::-1] if '_' not in  words[x]]
		top_words.append(top10)
		#all_words = [words[x] for x in range(len(each_dimension)) if each_dimension[x] > 0]
		#top_words.append(all_words)


	return top_words


def main(argv):

	parser = argparse.ArgumentParser(description="Interpret sparse embeddings")

	parser.add_argument("vector_file", help='The file from which sparse vectors have to be loaded')

	args = parser.parse_args()

	#Load the vectors
	words, embeddings_array = load_vectors(args.vector_file)
	top_words = get_top_words(words, embeddings_array, 10)

	count =0
	for each in top_words:
		print count, "\t", each
		count += 1

if __name__ == "__main__":
	main(sys.argv[1:])
