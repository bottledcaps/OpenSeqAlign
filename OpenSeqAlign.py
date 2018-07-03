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
	
	def nwalign(string1, string2, match_score = 1, mismatch_score = -1, gap_penalty = default_gap_penalty()):
		"""Aligns 2 strings, string1 and string2. Default values can be modified by passing in different values.
		This needleman_wunsch implementation does not penalize gaps at the beginning. 
		For example, the following strings align as shown since gaps at the beginning and end do not count against the alignment.
		>>>nwalign('GATTACA', 'GCATGCU')
		('GATTACA====', '====GCATGCU')
		...
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
					"""removes starting gap penalties
					to keep starting gap penalties set g[0][x] to be gap_penalty(x)
					"""
					a[0][x] = gap_penalty(x)
				for y in range(y_length):
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
				

				for y in range(1, length2+1):
					for x in range(1, length1+1):
						score(y, x)
						
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

	def no_end_gap_penalty_nwalign(string1, string2, match_score = 1, mismatch_score = -1, gap_penalty = default_gap_penalty()):
		"""Aligns 2 strings, string1 and string2. Default values can be modified by passing in different values.
		This needleman_wunsch implementation does not penalize gaps at the beginning. 
		For example, the following strings align as shown since gaps at the beginning and end do not count against the alignment.
		>>>nwalign('GATTACA', 'GCATGCU')
		('GATTACA====', '====GCATGCU')
		...
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
					"""removes starting gap penalties
					to keep starting gap penalties set g[0][x] to be gap_penalty(x)
					"""
					a[0][x] = 0
				for y in range(y_length):
					a[y][0] = 0
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
	
