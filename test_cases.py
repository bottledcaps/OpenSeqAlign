from OpenSeqAlign import *
from aligned_reads import *
from counting_dict_objs import *

class TestFunctions():
	def test_nwalign(self):
		string1, string2 = OpenSeqAlign.nwalign('GCATGCU', 'GATTACA')
		assert string1 == '====GCATGCU'
		assert string2 == 'GATTACA===='
		string1, string2 = OpenSeqAlign.nwalign('GCATGCU', 'GATTACA', no_end_gap_penalty = False)
		assert string1 == 'GCATGCU'
		assert string2 == 'GATTACA'
class TestObjects():
	def test_encoded_reads_object(self):
		string1, string2 = ('GCATG=CU', 'G=ATTACA')
		new_encoded_read = encoded_read(string1, string2)
		assert new_encoded_read.str == 'GATG=CU'
		assert new_encoded_read.insertion_dict == {0: 'c'}
		string1, string2 = new_encoded_read.get_alignment_to('GATTACA')
		assert string1 == 'GCATG=CU'
		assert string2 == 'G=ATTACA'
		string1 = '____GCATGCU'
		string2 = 'GATTACA____'
		new_encoded_read = encoded_read(string1, string2)
		assert new_encoded_read.str == '____GCA' 
		assert new_encoded_read.insertion_dict == {6: 'tgcu'}
		new_encoded_read = encoded_read()
		assert new_encoded_read.str == ''
		assert new_encoded_read.insertion_dict == {}
	def test_counting_dict(self):
		counters = counting_dict()
		assert counters.total == 0
		assert counters.dict == {}
		counters.add_term('a')
		assert counters.total == 1
		assert counters.dict == {'a': 1}
		counters.add_term('c')
		assert counters.total == 2
		assert counters.dict == {'a': 1, 'c': 1}
		counters.add_term('_')
		assert counters.total == 2
		assert counters.dict == {'a': 1, 'c': 1}
		counters.add_term('a')
		assert counters.total == 3
		assert counters.dict == {'a': 2, 'c': 1}
		assert counters.get_highest() == 'a'
	def test_dict_of_counting_dicts(self):
		ins_dict = dict_of_counting_dicts()
		for index in range(0, 7):
			for number in range(4, 7):
				ins_dict.add_term_at_index(number, index)
		for index in range(0, 7):
			assert ins_dict.dict[index].dict == {4: 1, 5: 1, 6: 1}
		ins_dict.add_term_at_index(5, 0)
		assert ins_dict.dict[0].dict == {4: 1, 5: 2, 6: 1}
		assert ins_dict.dict[0].total == 4
	def test_resolve(self):
		temp_count_dict = counting_dict()
		for _ in range(10):
			temp_count_dict.add_term('acattgtta')
		assert temp_count_dict.dict == {'acattgtta': 10}
		assert temp_count_dict.resolve() == 'acattgtta'
	#def test_resolve_all(self):
	def test_dict_of_counting_dicts_init(self):

		string1, string2 = ('GCATG=CU', 'G=ATTACA')
		new_encoded_read = encoded_read(string1, string2)
		new_encoded_read2 = encoded_read(string1, string2)
		new_encoded_read3 = encoded_read(string1, string2)
		encoded_reads_list = [new_encoded_read, new_encoded_read2, new_encoded_read3]
		new_dict_of_counting_dicts = dict_of_counting_dicts(encoded_reads_list)
		assert new_dict_of_counting_dicts.dict[0].dict == {'c': 3}
		assert new_dict_of_counting_dicts.dict[0].total == 3
		assert new_dict_of_counting_dicts.dict.keys() == [0]

