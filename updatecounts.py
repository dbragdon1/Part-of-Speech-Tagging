from collections import defaultdict
original_file = open('gene2.train', 'r')

word_dict = defaultdict(int)
tag_dict = defaultdict(int)

#Get word counts from original file without Informative Word Classes
for line in open('gene.COUNTS'):
	if 'WORDTAG' in line:
		data = line.strip().split(maxsplit = 3)
		word_dict[data[3]] += int(data[0])
		tag_dict[data[2]] += int(data[0])

#Given a word, assign it to an appropriate informative class
#Assumes word shows up less than 5 times in training data
def return_informative_class(word):
	if True in [char.isdigit() for char in word]:
		return '_DIGIT_'
	elif [char.isupper() for char in word] == ([True] * len(word)):
		return '_CAPS_'
	if word[-1].isupper():
		return '_LASTUPPER_'
	else: return '_RARE_'
	#return '_RARE_'

#Generate unlabeled train file for predicting on training data
unlabeled_train = open('geneunlabeled.TRAIN', 'w')
def create_train_unlabeled():
	for line in open('gene2.TRAIN'):
		if line.strip() == "":
			unlabeled_train.write('\n')
		else:
			word = line.strip().split()[0]
			unlabeled_train.write(word + "\n")
create_train_unlabeled()
unlabeled_train.close()


#Rewrite gene.TRAIN to contain informative word classes
new_file = open('gene.train', 'w')
i = 0
for line in original_file:
	data = line.strip().split()
	#print(data)
	if data == []:
		#new_file.write('\n')
		new_file.write("\n")

	else: 
		if word_dict[data[0]] < 5:
			new_word = return_informative_class(data[0])
			new_file.write(new_word + " " + data[1] + "\n")
		else: 
			new_file.write(data[0] + " " + data[1] + "\n")

new_file.close()
original_file.close()
