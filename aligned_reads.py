
"""
class consensus_char(counting_dict):
	def __init__(self):
		self.chars_dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '=': 0}
	def get_highest(self):
		max_value = max(self.chars_dict.values())
		for key in self.chars_dict.keys():
			if self.chars_dict[key] == max_value:
				return self.chars_dict
	def add_char(self, char):
		if char == '_':
			pass
		if char not in self.chars_dict.keys():
			raise InputError('character not a valid consensus_char type (must be A,C,G,T, or =)')
		self.chars_dict[char] += 1"""

class aligned_reads():
	def __init__(self, ref_string):
		self.ref_seq = ref_string
		self.encoded_reads = []
	def add_read(self, encoded_read):
		self.encoded_reads.append(encoded_read)
	def get_alignment(self, index):
		return decode(self.encoded_reads[index], self.ref_seq)
	def assemble_consensus_read(self):
		if hasattr(self, cons_read):
			print("Already created a consensus read for this sample, but redoing and overwriting")


		def resolve(self, combined_ins_dicts):
			#manually handle insertion at front and 0th index together, since each should compare to the 0th index of the encoded read strings
			
			first_index_char_dict = counting_dict()

			new_string = ''
			new_ins_dict = {}

			for encoded_read in self.encoded_reads:
				first_index_char_dict.add_term(encoded_read.str[0])

			if -1 in combined_ins_dicts.dict:
				ins_str = combined_ins_dicts.dict[-1].calc_ins(first_index_char_dict.total)
				if ins_str != None:
					new_ins_dict[-1] = ins_str

			if 0 in combined_ins_dicts.dict:	
				ins_str = combined_ins_dicts.dict[0].calc_ins(first_index_char_dict.total)
				if ins_str != None:
					new_ins_dict[0] = ins_str

			new_string = new_string + first_index_char_dict.get_highest_key()


			for index in range(1, len(self.ref_seq)):
				cons_char = counting_dict()
				for encoded_read in self.encoded_reads:
					cons_char.add_term(encoded_read.str[index])
				new_string = new_string + cons_char.get_highest_key()
				ins_str = combined_ins_dicts.dict[index].calc_ins(cons_char.total)
				if ins_str != None:
					new_ins_dict[index] = ins_str
			return new_string, new_ins_dict



		self.cons_read = encoded_read()
		combined_ins_dicts = dict_of_counting_dicts(self.encoded_reads)
		cons_read.str, cons_read.consensus_read = self.resolve(combined_insertion_dicts)

		#at each index
		"""for index in range(len(self.ref_seq)):
			#create a consensus char object
			consensus_char_dict = counting_dict()
			for encoded_read in aligned_reads.encoded_reads:
				consensus_char_dict.add_term(encoded_read.str[index])
			#and then add the highest char to the new encoded read
			new_encoded_read.str = new_encoded_read.str + consensus_char_dict.get_highest()

		insertions_dict = dict_of_counting_dicts()
		#for every possible insertion index...
		for insertion_index in range(-1, len(self.ref_seq)):
			#iterate through the encoded reads...
			for encoded_read in aligned_reads.encoded_reads:
				#if the current examined index is an insertion index for this particular encoded read,
				if insertion_index in encoded_read.insertion_dict.keys():
					#add the term to the insertion dict of counting dicts
					insertions_dict.add_term_at_index(encoded_read[insertion_index], insertion_index)
"""
class encoded_read():
	def __init__(self, aligned_read = None, aligned_ref = None):
		#separate way to initialize, with blank data, for when building a new encoded read
		if aligned_read == None or aligned_ref == None:
			self.str = ''
			self.insertion_dict = {}
		else:
			def encode(aligned_read, aligned_ref):
				def transform(aligned_read, aligned_ref):
					"""Takes a side-by-side alignment and transforms it into just an transformed read, which can be reversed by the translate function. 
					Gaps in the aligned ref corresponding to insertions in the aligned read become lowercase characters in the transformed read
					"""
					length = len(aligned_ref)
					transformed_read = ''

					for index in range(length):
						if (aligned_ref[index] == '=' or aligned_ref[index] == '_') and aligned_read[index] != '=' and aligned_read[index] != '_':
							transformed_read = transformed_read + str.lower(aligned_read[index])
						else:
							transformed_read = transformed_read + (aligned_read[index])
					return transformed_read		
				def compress(transformed_read):
					"""Initializes an insertion dictionary where keys are location, pair is the insertion characters. 
					The location key maps to the index it is right after in the reference sequence. 
					For example {0: 'acc', -1: 'agtc'} maps to a starting gap of 'agtc' and a gap of 'acc' right after the 0th index of the reference
					>>> from aligned_reads import *
					>>> compress('GcATG=CU')
					('GATG=CU', {0: 'c'})
					"""	
					insertion_dict = {}
					compressed_read = ''
					index = 0
					ref_index = -1
					while index < len(transformed_read):
						#if the current character is a lowercase character...
						if str.islower(transformed_read[index]):
							insertion = ''
							#record the continuous string of lowercase letters until it ends
							while index < len(transformed_read) and str.islower(transformed_read[index]):
								insertion = insertion + transformed_read[index]
								index += 1
							#then put it into the dictionary with the start index as the key
							insertion_dict[ref_index] = insertion
						else:
							compressed_read = compressed_read + transformed_read[index]
							index += 1
							ref_index += 1
					return compressed_read, insertion_dict 
				return compress(transform(aligned_read, aligned_ref))
			self.str, self.insertion_dict = encode(aligned_read, aligned_ref)
	def get_alignment_to(self, ref_string):
		def decode(compressed_read, insertion_dict, ref_string):
			"""Reverses encoding. Decompresses then translates
			>>> from aligned_reads import *
			>>> decode('GATG=CU', {0:c}, 'GATTACA')
			('GCATG=CU', 'G=ATTACA')
			"""
			def translate(transformed_read, reference_string):
				"""Reverses the process of transformation. Takes a translated read and its original reference and puts them into side-by-side alignment form
				>>> from aligned_reads import *
				>>> translate('GcATG=CU', 'GATTACA')
				('GCATG=CU', 'G=ATACCA')
				>>> translate(transform('GCATG=CU', 'G=ATTACA'), 'GATTACA') #translation is the opposite of transforming
				('GCATG=CU', 'G=ATACCA')
				"""
				length = len(transformed_read)
				translated_ref = ''
				translated_read = ''
				ref_index = 0
				#for every character in the transformed read...
				for index in range(length):
					#if the current char is lowercase...
					if str.islower(transformed_read[index]):
						#add a gap to reference ...
						translated_ref = translated_ref + '=' 
						#and add an uppercase char to the read. Don't increase the ref_index.
						translated_read = translated_read + str.upper(transformed_read[index])
					else:
						#Otherwise add a character from each to the new strings, and increment the ref_index.
						translated_ref = translated_ref + reference_string[ref_index]
						translated_read = translated_read + transformed_read[index]
						ref_index += 1
				#additionally, post process the translated_ref for gaps at the beginning and end (these should be _'s)
				"""index = 0
				while translated_ref[index] == '=':
					translated_ref[index] = '_'
					index += 1
				index = len(translated_ref) - 1
				while translated_ref[index] == '=':
					translated_ref[index] = '_'
					index -= 1"""
				return translated_read, translated_ref
			def decompress(compressed_read, insertion_dict):
				"""Takes a commpressed read and turns it into a translated read by inserting terms from dictionary into string. Refer to figure for full diagram.
				>>> from aligned_reads import *
				>>> decompress('GATG=CU', {0, 'c'})
				'GcATG=CU'
				"""
				decompressed_read = ''
				keys_list = insertion_dict.keys()
				#first check if there is a beginning insertion
				if -1 in keys_list:
					decompressed_read = decompressed_read + insertion_dict[-1]
				#for each char in the compressed read...
				for index in range(len(compressed_read)):
					#make that char part of the new string
					decompressed_read = decompressed_read + compressed_read[index]
					#afterwards, check if this index has a insertion right after it..
					if index in keys_list:
						#and if it does, add that insertion into the decompressed read.
						decompressed_read = decompressed_read + insertion_dict[index]
				return decompressed_read
			return translate(decompress(compressed_read, insertion_dict), ref_string)
		return decode(self.str, self.insertion_dict, ref_string)