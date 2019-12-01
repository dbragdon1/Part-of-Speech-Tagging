
from collections import defaultdict
from math import log
from nltk import ngrams
from updatecounts import return_informative_class

unigramdict = defaultdict(int)
bigramdict = defaultdict(int)
trigramdict = defaultdict(int)
emission_dict = defaultdict(int)
word_dict = defaultdict(int)
tag_dict = defaultdict(int)

#Create dictionaries from Informative Word Class counts
for line in open('geneIWC.COUNTS'):
	if 'WORDTAG' in line:
		data = line.strip().split(maxsplit = 3)
		emission_dict[(data[2], data[3])] += int(data[0])
		word_dict[data[3]] += int(data[0])
		tag_dict[data[2]] += int(data[0])
	if "1-GRAM" in line:
		data = line.strip().split(maxsplit = 3)
		unigramdict[data[2]] += int(data[0])
		
	elif "2-GRAM" in line:
		data = line.strip().split(maxsplit = 4)
		bigram = (data[-2], data[-1])
		bigramdict[bigram] += int(data[0])

	elif "3-GRAM" in line:
		data = line.strip().split(maxsplit = 5)
		#print(data)
		trigram = (data[-3], data[-2], data[-1])
		#print(trigram)
		trigramdict[trigram] += int(data[0])

for bigram in bigramdict:
	if "*" in bigram:
		unigramdict['*'] += bigramdict[bigram]
	if "STOP" in bigram:
		unigramdict['STOP'] += bigramdict[bigram]

def compute_emission_parameters(tag, word):
	return emission_dict[(tag, word)] / tag_dict[tag]

emission_params = defaultdict(float)
for emission in emission_dict:
	emission_params[emission] = compute_emission_parameters(emission[0], emission[1])

def compute_trigram_prob(trigram):
	bigram = (trigram[0], trigram[1])
	return trigramdict[trigram] / bigramdict[bigram]

def compute_unigram_probs(unigram):
	return unigramdict[unigram] / sum(unigramdict.values())

def compute_bigram_probs(bigram):
	return bigramdict[bigram] / unigramdict[bigram[0]]

unigramprobs = defaultdict(float)
for unigram in unigramdict:
	unigramprobs[unigram] = compute_unigram_probs(unigram)

bigramprobs = defaultdict(float)
for bigram in bigramdict:
	bigramprobs[bigram] = compute_bigram_probs(bigram)

trigram_probs = defaultdict(float)
for trigram in trigramdict:
	trigram_probs[trigram] = compute_trigram_prob(trigram)

#Create original count dictionaries
#Used when predicting tags in baseline tagger
orig_word_dict = defaultdict(int)
orig_tag_dict = defaultdict(int)

for line in open('gene.COUNTS'):
	if 'WORDTAG' in line:
		data = line.strip().split(maxsplit = 3)
		orig_word_dict[data[3]] += int(data[0])
		orig_tag_dict[data[2]] += int(data[0])

def baseline_tagger(word):
	#stores emission for each tag
	possible_ems = []

	#stores possible tags
	possible_tags = []

	#Check if word is rare
	if orig_word_dict[word] < 5:
		#Iterate over possible tags
		for tag in ['O', 'I-GENE']:
			#Assign word to appropriate informative class
			
			new_word = return_informative_class(word)
			#Append probability of (tag, _CLASS_)
			possible_ems.append(emission_params[(tag, new_word)])
			#Append tag
			possible_tags.append(tag)
	#If word is not rare
	else: 
		for tag in ['O', 'I-GENE']:
			#Assign word to appropriate emission (tag, word)
			possible_ems.append(emission_params[(tag, word)])
			possible_tags.append(tag)
	#Find highest probability in emissions
	max_idx = possible_ems.index(max(possible_ems))
	#Find tag corresponding to highest probability
	argmax_tag = possible_tags[max_idx]
	#returns emission with original word, rather than informative class
	#Preserves word to be written to dev predictions file 
	emission = (argmax_tag, word)

	return emission

#Write Baseline Predictions to gene_dev.p1.out 
out1 = open("gene_dev.p1.out", 'w')
for line in open('gene.DEV'):
	word = line.strip()
	if word == "":
		out1.write("\n")
	else:
		pred_em = baseline_tagger(word)
		out1.write(pred_em[1] + " " + pred_em[0] + "\n")
out1.close()

trainout1 = open('gene_train.p1.out', 'w')
for line in open('geneunlabeled.TRAIN'):
	word = line.strip()
	if word == "":
		trainout1.write("\n")
	else:
		pred_em = baseline_tagger(word)
		#print(pred_em)
		trainout1.write(pred_em[1] + " " + pred_em[0] + "\n")
trainout1.close()


#words: sequence of words in the form of a 1D array
#gammas: gammas to use for interpolation
#		 [.7, .2, .1] are best combination I have found so far
#Tags: tags that can exist in position K if K is not -1, 0, or n + 1

def Viterbi(words, tags, Interpolation = False, gammas = [.7, .2, .1]):
	
	#initialize tags to empty array
	y = [None] * len(words)

	#probability dictionary
	pi = defaultdict(float)

	#backpointer dictionary
	backpointer = defaultdict(float)

	#####Helper functions#####

	#Interpolate trigram probabilities for transition parameters
	def interpolation(trigram, gammas):
		u,v,w = trigram[0], trigram[1], trigram[2]
		#trigram
		first_term = (gammas[0] * trigram_probs[(u,v,w)])
		#bigram
		second_term = (gammas[1] * bigramprobs[(v,w)]) 
		#unigram
		third_term = (gammas[2] * unigramprobs[w]) 
    	
		return first_term + second_term + third_term

	#find possible states at current time step
	def return_possible_states(i):
		#if time step is less than zero, tag can only be '*'
		if i <= 0:
			return ['*']
		else:
			#return ['O', 'I-GENE']
			return [tag for tag in tags]

	###### begin algorithm ######

	#Initialize prob dict
	pi[(0, "*", "*")] = 1

	n = len(words)

	#create indices where gene tags can exist (1...n)
	K = range(1, n + 1)

	for k in K:
		#iterate over S_i-1
		for u in return_possible_states(k - 1):
			#iterate over S_i
			for v in return_possible_states(k):
				#initialize best probability to -1 so that all successive probs will be bigger
				best_prob = -1
				best_w = None
				#iterate over S_i-2
				for w in return_possible_states(k - 2):
					word = words[k - 1]

					if (word_dict[word] < 5):
						new_word = return_informative_class(word)
	                	#assign _RARE_ tag to word
						e_p = emission_params[(v, new_word)]

					else:
	                	#assign normal tag to word
						e_p = emission_params[(v, word)]

	            	#interpolation step
					if Interpolation == True:
						interpol = interpolation((w, u, v), gammas)
						p_curr_w = pi[(k - 1, w, u)] * interpol * e_p
					else:
						p_curr_w = pi[(k - 1), w, u] * trigram_probs[(w, u, v)] * e_p

					#calculate pi(k, u, v) = pi(k - 1, u, v) + q(v | w, u) + e(x_k | v)
					#for current w
					#p_curr_w = pi[(k - 1, w, u)] * interpol * e_p
					#p_curr_w = pi[(k - 1), w, u] * trigram_probs[(w, u, v)] * e_p
					#if probability from current w is higher than previous w's
					if p_curr_w > best_prob:
						best_prob = p_curr_w
						best_w = w

	            #store best prob at kth state in pi
				pi[(k, u, v)] = best_prob
				#store best w at kth state in backpointer matrix
				backpointer[(k, u, v)] = best_w

	#Proposition 2 in Collins Notes
    #return maximum probability for ending trigram
    #initialize probability to low number
	
	endingmax = -1

    #iterate over position S_{n-1}
	for u in return_possible_states(n - 1):
		#iterate over position S_{n}
		for v in return_possible_states(n):
			#assign STOP to S_{n + 1}

			#get (u,v,STOP) trigram and interpolate the probability
			if Interpolation == True:
				interpol = interpolation((u, v, "STOP"), gammas)
				p_y = pi[(n, u, v)] * interpol

			else:
				p_y = pi[(n, u, v)] * trigram_probs[(u, v, "STOP")]

			if p_y > endingmax:

            	#update max for higher probabilities
				endingmax = p_y

				#u = S_{n - 2}
				#v = S_{n - 1}
				y[n - 2] = u
				y[n - 1] = v
    
    #traverse backpointer matrix for tags
	for k in range(n - 3, -1, -1):

		#y_k = bp(k + 2, y_{k+1}, y__{k+2})
		y[k] = backpointer[(k + 3, y[k + 1], y[k + 2])]
	return y

sent = []
out2 = open('gene_dev.p2.out', 'w')
for line in open('gene.DEV', 'r'):
	word = [line.strip()]
	if word == ['']:
		tags = Viterbi(sent, tags = ['O', 'I-GENE'], Interpolation = True, gammas = [.8, .1, .1])
		for word, tag, in zip(sent, tags):
			out2.write(word + " " + tag + "\n")
		out2.write("\n")
		sent = []
		
	else:
		sent += word
out2.close()

words = []
for line in open('geneunlabeled.TRAIN'):
	word = [line.strip()]
	words.append(word)
size = len(words)
idx_list = [idx + 1 for idx, val in enumerate(words) if val == ['']]
res = [words[i: j] for i, j in zip([0] + idx_list, idx_list + ([size] if idx_list[-1] != size else []))]

new_sents = []
for sent in res:
	new_sent = []
	for word in sent:
		if word == ['']:
			continue
		else: new_sent += word
	new_sents.append(new_sent)
	new_sent = []

trainout2 = open('gene_train.p2.out', 'w')
for sent in new_sents:
	tags = Viterbi(sent, tags = ['O', 'I-GENE'], Interpolation = True, gammas = [.7, .2, .1])
	for word, tag, in zip(sent, tags):
		trainout2.write(word + " " + tag + "\n")
	trainout2.write("\n")

trainout2.close()
