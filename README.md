# Part-of-Speech-Tagging
Part of Speech Tagging on Gene Sequences using the Viterbi Algorithm

Files:

gene.TRAIN: Train file rewritten to include informative word classes

gene2.TRAIN: Train file used as reference to preserve original words after rewriting gene.TRAIN

gene.COUNTS: File used to store original counts without informative word classes

geneIWC.COUNTS: File used to store counts with informative word classes

gene_dev.p1.OUT: File used to output predictions on dev set from baseline tagger

gene_dev.p2.OUT: File used to output predictions on dev set from Trigram HMM

gene_train.p1.OUT: File used to output predictions on train set from baseline tagger

gene_train.p2.OUT: File used to output predictions on train set from Trigram HMM

geneunlabeled.TRAIN: Words from gene.TRAIN with labels removed. Used when predicting on test set

updatecounts.py:

	- Purpose:
		- Takes original train data (gene2.TRAIN)
		- 
		- Separates words into informative word classes
		- Returns new gene.TRAIN file with informative word classes 
	- Methods:
		- return_informative_class(word):
			- Takes word and decides if it belongs to informative word class
			- Returns corresponding word class
			- To use no informative word classes, comment out every line and uncomment last line
	- Usage: 
		1. Run count_freqs.py gene2.train > gene.counts to get original gene.COUNTS file
		2. Run updatecounts.py to get new gene.COUNTS and gene.TRAIN file with informative word classes
		3. Run count_freqs.py gene.train > gene.counts to get new gene.COUNTS file with informative word classes 

trigramhmm.py: 

	- Purpose:
		- Computes emission parameters
		- Contains both baseline model and TrigramHMM model with Viterbi Algorithm
		- Writes output predictions to gene_dev.p1.OUT, gene_dev.p2.OUT, and gene_train.p1.OUT
		- Contains Extensions (Interpolation)
	- Methods:
		- compute_trigram_prob(trigram):
			- Computes probability of seeing trigram
		- compute_bigram_probs(bigram):
			- Computes probability of seeing bigram
		- compute_unigram_probs(unigram):
			- Computes probability of seeing unigram
		- baseline_tagger(word):
			- Takes word and returns tag from emission parameter with highest prob
			- Incorporates return_informative_class() method
		- Viterbi(words, tags, Interpolation = False, gammas = [.7, .2, .1])
			- Takes in sentence in the form of a 1D array of strings
			- Outputs optimal sequence of tags
			- If Interpolation = True: Incorporates simple interpolation
			- gammas = Coefficients for interpolation
		- interpolation(trigram, gammas):
			- Calculates interpolation probability given trigram
			- gammas: Coefficients used in interpolation calculation
			- Helper method inside Viterbi() method
		- return_possible_states(i):
			- Returns possible states given position i in sentence
			- Helper method inside Viterbi() method

	- Usage:
		- Run trigramhmm.py to generate gene_dev.p1.OUT, gene_dev.p2.OUT, gene_train.p1.OUT, and gene_train.p2.OUT


