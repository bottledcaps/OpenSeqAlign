class dict_of_counting_dicts():
	def __init__(self, encoded_reads = None):
		
		self.dict = {}
		if encoded_reads != None:	
			for encoded_read in encoded_reads:
				for key, v in encoded_read.insertion_dict.items():
					self.add_term_at_index(v, key)
	def add_term_at_index(self, term, index):
		if index in self.dict:
			self.dict[index].add_term(term)
		else:
			self.dict[index] = counting_dict()
			self.dict[index].add_term(term)


class counting_dict():
	def __init__(self):
		self.dict = {}
		self.total = 0
	def add_term(self, term):
		if term == '_':
			pass
		elif term not in self.dict.keys():
			self.dict[term] = 1
			self.total += 1
		else:
			self.dict[term] += 1
			self.total += 1
	def get_highest_pair(self):
		max_value = max(self.dict.values())
		for key, v in self.dict.items():
			if v == max_value:
				return key, v
	def get_highest_key(self):
		max_value = max(self.dict.values())
		for key, v in self.dict.items():
			if v == max_value:
				return key
	
	def calc_ins(self, loc_total):
		#if the highest-occurring insertion string is at least half of the total number of reads which pass by this index
		#this process can be more complicated if needed?? Such as a hamming-distance relationship between strings of identical length supporting each other
		#since aca and act are ""similar"" in some level which would be nice to implement
		string, value = self.get_highest
		if value >= (loc_total // 2) - 15:
			#resolve the counting dictionary into that insertion string
			return string
		#if the insertions don't meet the minimum barrier to manifest a real insertion, delete the 
		else:
			return None
			
