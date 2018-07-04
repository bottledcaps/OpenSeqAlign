class OpenSeqAlign:
	def default_gap_penalty(open_penalty=10, extend_penalty = 1): 
		"""the default gap penalty calculation, always outputs positive values
		to use a different gap penalty, pass a different gap_penalty function into needleman_wunsch
		default_gap_penalty returns a penalty function of the form penalty = open_penalty + extend_penalty*(length-1) where length is an input
		"""

		def gap_penalty(length):
			if length <= 0:
				return 0
			return open_penalty + extend_penalty*(length-1)
		return gap_penalty
	
	def nwalign(string1, string2, no_end_gap_penalty = True, match_score = 1, mismatch_score = -1, gap_penalty = default_gap_penalty()):
		"""Aligns 2 strings, string1 and string2. Default values can be modified by passing in different values.
		This needleman_wunsch implementation does not penalize gaps at the beginning and ends, by default. 
		For example, the following strings align as shown since gaps at the beginning and end do not count against the alignment by default.
		>>>nwalign('GATTACA', 'GCATGCU')
		('GATTACA====', '====GCATGCU')
		"""
		def create_path_grid():
			def init_grid(value, y_length = len(string2) + 1, x_length = len(string2) + 1):
				a = []
				for _ in range(y_length):
					a.append([])
				for x in range(y_length):
					for y in range(x_length):
						a[x].append(value)
				return a
			def init_score_grid(y_length = len(string2) + 1, x_length = len(string2) + 1):
				a = init_grid(None)
				for x in range(x_length):
					if no_end_gap_penalty:
						a[0][x] = 0
					"""removes starting gap penalties
					to keep starting gap penalties set g[0][x] to be gap_penalty(x)
					"""
					else:
						a[0][x] = gap_penalty(x)
				for y in range(y_length):
					if no_end_gap_penalty:
						a[0][x] = 0
					else:
						a[y][0] = gap_penalty(y)
				return a
			def init_path_grid(y_length = len(string2) + 1, x_length = len(string2) + 1	):
				a = init_grid(None)
				#elements along top and left edge point to top left entry
				for x in range(1, x_length):
					a[0][x] = -1
				for x in range(1, y_length):
					a[x][0] = 1
				return a

			def calc_matrices():
				def subst_score(char1, char2):
					if char1 == 'N' or char2 == 'N':
						return match_score
					elif char1 == char2:
						return match_score
					else:
						return mismatch_score
				
				def score(y_loc, x_loc):
					def calc_gaps():
						#key is value, pair indicates origin point
						#if pair is 0, no gap and up left
						gap_dict = {}
						#go up
						k = 1
						while (y_loc-k) >= 0 and scoring_matrix[y_loc-k][x_loc] != None:
							gap_dict[scoring_matrix[y_loc-k][x_loc] - gap_penalty(k)] = k

							k = k + 1
						#go left
						j = 1
						while (x_loc-j) >= 0 and scoring_matrix[y_loc][x_loc - j] != None:
							gap_dict[scoring_matrix[y_loc][x_loc-j] - gap_penalty(j)] = -j

							j = j + 1
						return gap_dict
					

					score_path = calc_gaps()
					calc_sub = scoring_matrix[y_loc-1][x_loc-1] + subst_score(string1[x_loc-1], string2[y_loc-1])
					score_path[calc_sub] = 0
					max_key = max(score_path.keys())
					path_matrix[y_loc][x_loc] = score_path[max_key]
					scoring_matrix[y_loc][x_loc] = max_key
				
				def fix_end():
					#changes value in bottom right of path_grid to point to maximum value in same row or column, 
					#or if the value from up-left is actually best with this alt calculaution
					y_loc = len(string2)
					x_loc = len(string1)
					end_dict = {}
					end_dict[scoring_matrix[y_loc - 1][x_loc - 1] + subst_score(string1[x_loc-1], string2[y_loc-1])] = 0
					for k in range(1, y_loc):
						end_dict[scoring_matrix[y_loc-k][x_loc]] = k
					for j in range(1, x_loc):
						end_dict[scoring_matrix[y_loc][x_loc-j]] = -j
					max_key = max(end_dict.keys())
					path_matrix[y_loc][x_loc] = end_dict[max_key]

				for y in range(1, length2+1):
					for x in range(1, length1+1):
						score(y, x)
				if no_end_gap_penalty:
					fix_end()
			scoring_matrix = init_score_grid()
			path_matrix = init_path_grid()
			calc_matrices()

			return path_matrix

		def traceback():
			new_string1 = ''
			new_string2 = ''
			current_x = len(string1) 
			current_y = len(string2)
			while(path_grid[current_y][current_x] != None):
				#k is current_path_pointer
				k = path_grid[current_y][current_x]

				if k == 0:
					new_string1 = string1[current_x-1] + new_string1
					new_string2 = string2[current_y-1] + new_string2
					current_y -= 1
					current_x -= 1
				elif k > 0:

					new_string2 = string2[current_y-k:current_y] + new_string2
					new_string1 = k*'=' + new_string1
					current_y -= k
				elif k < 0:
					new_string1 = string1[current_x - (-k):current_x] + new_string1
					new_string2 = (-k)*'=' + new_string2
					current_x -= (-k)
			return new_string1, new_string2

		path_grid = create_path_grid()
		new_string1, new_string2 = traceback()
		return new_string1, new_string2

	def heur_nwalign(string1, string2, max_displacement, match_score = 1, mismatch_score = -1, gap_penalty = default_gap_penalty()):
		def get_indexes(row_number, x_loc): #converts from diagonal matrix to string indexes
			if (x_loc - max_displacement) < 0:
				first = row_number - 1
				second = row_number - 1 + (max_displacement - x_loc)
			else:
				first = row_number - 1 + (x_loc - max_displacement)
				second = row_number - 1
			return first, second
		def indexes_to_row_x(string1_index, string2_index):
			if string1_index > string2_index:
				row_number = string2_index + 1
				x_loc = (string1_index - string2_index) + max_displacement
				return row_number, x_loc
			else:
				row_number = string1_index + 1
				x_loc = max_displacement - (string2_index - string1_index) 
				return row_number, x_loc
		def create_heur_path_grid():
			def abs(input):
				if input < 0:
					return -input
				return input
			def init_grid(value, y_length = shorter_len + 1, x_length = row_width):
			#given two strings, create blank scoring matrix with dimensions 1 greater than each
				a = []
				for _ in range(y_length):
					a.append([])
				for x in range(y_length):
					for y in range(x_length):
						a[x].append(value)
				return a
			def init_scoring_grid():
				a = init_grid(None)
				for x in range(len(g[0])):
					a[0][x] = -gap_penalty(abs(x - max_displacement))
				return a
			def init_path_grid():
				g = init_grid(None)
				for x in range(len(g[0])):
					if (x - max_displacement) < 0:
						g[0][x] = -1
					if (x - max_displacement) > 0:
						g[0][x] = 1
				return g
			def calc_matrices():
				def score(row_number, x_loc):
					#Can do non-linear 
					#DICTIONARIES: score is key, location is value
					#FIX GAP PENALTY

					string1_index, string2_index = get_indexes(row_number, x_loc)
					def calc_gaps():
						gap_dict = {}
						#go_down, gaps in string2
						k = 1
						while (string2_index - k) >= 0 and (x_loc + k) < row_width:
							#print(row_number, x_loc)
							#print(string1_index, string2_index)
							a, b = indexes_to_row_x(string1_index, string2_index - k)
							#print(a,b)
							gap_dict[scoring_matrix[a][b] - gap_penalty(k)] = k
							k += 1
						j = 1
						while (string1_index - j) >= 0 and (x_loc - j) >= 0:
							a, b = indexes_to_row_x(string1_index - j, string2_index)
							
							gap_dict[scoring_matrix[a][b] - gap_penalty(j)] = -j
							j += 1
						return gap_dict

					if string1_index >= length1 or string2_index >= length2:
						return None
					score_path_dict = calc_gaps()

					calc_sub = scoring_matrix[row_number - 1][x_loc] + sub_score(string1[string1_index] == string2[string2_index])
					score_path_dict[calc_sub] = 0
					max_key = max(score_path_dict.keys())
					path_grid[row_number][x_loc] = score_path_dict[max_key]

					return max_key		
					for row in range(shorter_len):
				for a in range(max_displacement+1):
					scoring_matrix[row+1][max_displacement+a] = score(row+1, max_displacement + a)
					scoring_matrix[row+1][max_displacement-a] = score(row+1, max_displacement - a)

			scoring_matrix = init_scoring_grid()
			path_matrix = init_path_grid()
			calc_matrices
			return path_matrix
		def traceback():
			new_string1 = ''
			new_string2 = ''
			current_row = shorter_len
			current_x_loc = max_displacement

			while(path_grid[current_row][current_x_loc+length1_flag]!=None):
				current_x_loc += length1_flag

			current_string1_index, current_string2_index = get_indexes(current_row, current_x_loc)

			
			while(path_grid[current_row][current_x_loc]) != None:
				#k is current_path_pointer
				
				k = path_grid[current_row][current_x_loc]
				if k == 0:
					new_string1 = string1[current_string1_index] + new_string1
					new_string2 = string2[current_string2_index] + new_string2
					current_row -= 1
					current_string1_index -= 1
					current_string2_index -= 1
				elif k > 0:
					new_string2 = string2[(current_string2_index - k)+1:current_string2_index+1] + new_string2
					new_string1 = k*'_' + new_string1
					current_string2_index -= k
					current_row, current_x_loc = indexes_to_row_x(current_string1_index, current_string2_index)
				elif k < 0:
					new_string1 = string1[(current_string1_index - (-k))+1:current_string1_index+1] + new_string1
					new_string2 = (-k)*'_' + new_string2
					current_string1_index -= (-k)
					current_row, current_x_loc = indexes_to_row_x(current_string1_index, current_string2_index)

			return new_string1, new_string2
	#initializations of values needed throughout program
	length1 = len(string1)
	length2 = len(string2)
	length1_flag = -1
	if length1 > length2:
		longer_len = length1
		shorter_len = length2
		length1_flag = 1
	else: 
		longer_len = length2
		shorter_len = length1

	row_width = (max_displacement * 2) + 1

	#function body / "doing stuff" 
	path_grid = creat_heur_path_grid()
	new_string1, new_string2 = traceback()
	return new_string1, new_string2

	def smith_waterman_align(string1, string2, match_score = 3, mismatch_score = -3):
		#strings 1 and 2 are input strings, presumably genomes
		#have to determine rules?
		length1 = len(string1)
		length2 = len(string2)
		matrix_max = 0
		x_max = None
		y_max = None
		def initialize_grid(value):
		#given two strings, create blank scoring matrix with dimensions 1 greater than each
			A = []
			for _ in range(length2 + 1):
				A.append([])
			for x in range(length2 + 1):
				for y in range(length1+1):
					A[x].append(value)
			return A
		def print_grid(grid):
			a = len(grid)
			for x in range(a):
				#if x != 0:
				#	print(string2[x-1])
				print(grid[x])
		def sub_score(tf):
			if tf:
				return match_score
			else:
				return mismatch_score
		def gap_penalty(length):
			return length
		def score(y_loc, x_loc):
			#return x_loc + y_loc
			nonlocal matrix_max, x_max, y_max
			max_val = 0

			calc_sub = g[y_loc-1][x_loc-1] + sub_score((string1[x_loc-1] == string2[y_loc-1]))
			calc_string1_gap = g[y_loc-1][x_loc] - 2
			calc_string2_gap = g[y_loc][x_loc-1] - 2

			subst = False
			up = False
			left = False


			if calc_sub > max_val:
				max_val = calc_sub
				#path_x, path_y = x_loc - 1, y_loc - 1
				subst = True

			if calc_string1_gap > max_val:
				max_val = calc_string1_gap
				#path_x = x_loc
				#path_y = y_loc-1
				up = True

			if calc_string2_gap > max_val:
				max_val = calc_string2_gap
				#path_x = x_loc-1
				#path_y = y_loc
				left = True

			if max_val > matrix_max:
				matrix_max = max_val
				x_max = x_loc
				y_max = y_loc



			if left:
				path_grid[y_loc][x_loc] = -1
			elif up:
				path_grid[y_loc][x_loc] = 1
			elif subst:
				path_grid[y_loc][x_loc] = 0 

			return max_val
		def align_strings():
			new_string1 = ''
			new_string2 = ''
			current_x = x_max
			current_y = y_max
			#update messes up the if statements
			while(path_grid[current_y][current_x] != None):
				if path_grid[current_y][current_x] == 0:

					new_string1 = string1[current_x-1] + new_string1
					new_string2 = string2[current_y-1] + new_string2
					current_y -= 1
					current_x -= 1
				elif path_grid[current_y][current_x] == 1:
					new_string2 = string2[current_y-1] + new_string2
					new_string1 = '-' + new_string1
					current_y -= 1
				elif path_grid[current_y][current_x] == -1:
					current_x -= 1
			return new_string1, new_string2


		g = initialize_grid(0)
		path_grid = initialize_grid(None) #negative x means to the left by x, positive x means up by x, 0 means up, left


		#x indexes spot in string 1, y indexes spot in string 2, [y+1][x+1] is current spot
		for y in range(length2):
			for x in range(length1):
				g[y+1][x+1] = score(y+1, x+1)
				
		new_string1, new_string2 = align_strings()
		return new_string1, new_string2