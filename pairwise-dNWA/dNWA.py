# Compute the optimal (best-scoring) global alignment of two sequences with NWA for all sequences in the provided input file.
###############################################################

# import modules 
import itertools
import argparse
import math
import logging
import re
import os

# global variable definitions
number_of_alignments = 0
number_of_lines_in_infile = 0
number_of_lines_in_interfile = 0
number_of_lines_in_outfile = 0

# GAP values
PENALTY_GAP_OPENING = -1
PENALTY_GAP_EXTENSION = -0.1

# Scoring Matrix
WEIGHT_MATCH_F_GROUP = 8
WEIGHT_MATCH_T_GROUP = 4
WEIGHT_MATCH_H_GROUP = 2
WEIGHT_MATCH_X_GROUP = 1
WEIGHT_MATCH_NO_FID = 0
WEIGHT_MISMATCH_FIDS = -5
WEIGHT_MISMATCH_NOFID_FID = -1

# parsing arguments from commandline 
parser = argparse.ArgumentParser(description='aligns all input sequences pairwise with a modified Needleman-Wunsch algorithm.',epilog='Generates two output files per default with provided filename at the end. Output format: >ClusID of seq1,accession of seq1,>ClusID of seq2,accession of seq2, Alignmentscore, length of alignment')
parser.add_argument('filename', help='please provide an input file.')
parser.add_argument('-v', '--verbose', help='writes a second output file that is more human-readable.',action='store_true') 
parser.add_argument('-l', '--log', help='set loglevel',action='store', const='INFO', nargs='?',choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'])
parser.add_argument('-t', '--temp', help='keep temporary file with reformatted input',action='store_true')
parser.add_argument('-s','--noselfalignment', help='Deactivate the alignment of every Sequence with itself',action='store_true')
parser.add_argument('-g','--gapextension', help='Differentiate between gap openings and gap extensions',action='store_true')
parser.add_argument('-o', '--outfolder',help="specifiy location for output files")
args = parser.parse_args()

# activate logging if -log is set. Set log level if provided, otherwise set INFO as default
if args.log:
	numeric_level = args.log
	substring = '.txt'
	filename_wo_extension = os.path.basename(args.filename).replace(substring,'')
	logging.basicConfig(filename=args.outfolder+'/'+'basic_log_'+filename_wo_extension+'.log', encoding='utf-8', filemode='w', format='%(levelname)s:%(message)s (%(asctime)s)', datefmt='%m/%d/%Y %I:%M:%S %p', level=numeric_level)

# define functions
def transformRawInput():
	"""Transform the raw input into a format, that NWA can handle, and store it in a txt-file."""

	# open input file and close it in the end
	with open(args.filename, "r") as infile:
	
		# Create an empty array to hold the parsed data
		data = []
	
		# Initialize variables
		clusID = None
		specie_info = None
		cluster_count = 0
		numberOfLinesRead = 0
		totalAccessions = 0
		countedAccessions = 0
		global number_of_lines_in_infile
		filename = args.outfolder+'/'+'formatted_input_'+os.path.basename(args.filename)
		
		# Read the file line by line
		line = infile.readline()
		while line:
			numberOfLinesRead = numberOfLinesRead +1
			if line == '':
				continue
				
			# if line contains cluster infos, store them
			if line.startswith(">"):
				if clusID is not None:
					if args.log:
						logging.info('The number of read lines, that should be accessions, in cluster ' + str(clusID) + ' is '+str(countedAccessions)+'.')
				
				# set counter variables
				cluster_count = cluster_count+1
				countedAccessions = 0
				
				# If the line starts with ">", store the clusID in the clusID variable
				temp = re.findall(r'\d+', line)
				numbersInLine = list(map(int, temp))
				clusID = numbersInLine[0] # clusclusID
				totalAccessions = numbersInLine[1]
				if args.log:
					logging.info('The number of accessions in cluster ' + str(clusID) + ' is indicated in the file as '+str(totalAccessions)+'.')
			
			# else split the line into parts using the clusIDentifier as separator and write a transformed line into an intermediate file
			else:
				
				#assertion, that cluster infos are provclusIDed 
				if clusID is None:
					if args.log:
						logging.error('Cluster clusID is empty.')
					print('A cluster clusID is empty, please check input file')

				# count variable for checking number of accessions
				countedAccessions = countedAccessions+1
				
				parts = line.split(" ")			
				first_entry = parts[0] # extracts unique clusID of sequence
				last_entry = parts[-1] # extracts length of sequence
				
				# store the length of the sequence, the x as found in (x,y,END)
				parts_length = int(last_entry.split(',')[0][1:])
				
				# Write the output to a file
				try:
					with open(filename, 'a') as outfile:
					
						# Write the first line with "clusID,name of the sequence" to store this info with the sequence
						outfile.write('>' + str(clusID) + ',' + first_entry + '\n')

						# declare array for the second line with zeros (as int)
						sequence_list = [0]* parts_length
						
						# Write the numbers for domains into the list "parts"
						for entry in parts[2:-1]:
							parts = entry.split(',')
							motif_number = parts[-1][:-1]
							first_number = int(parts[0][1:])					
							second_number = int(parts[1])
							
						# write the domains	into the sequence_list
							for index in range(first_number,second_number):
								sequence_list[index] = motif_number

						# write data into a file
						output_string = ','.join(map(str, sequence_list))
						outfile.write(output_string + '\n')
				
				#handle exception
				except IOError:
					if args.log:
						logging.critical("The file cannot be opened")
					print("The file cannot be opened")
			
			# Read the next line
			line = infile.readline()
	
	#store number global for checks
	number_of_lines_in_infile = numberOfLinesRead		
			
	# write infos to log if logging activated
	if args.log:
		logging.info('The number of read lines, that should be accessions, in cluster ' + str(clusID) + ' is '+str(countedAccessions)+'.')
		logging.info('Input transformed')
		logging.info('#cluster= '+str(cluster_count))
		logging.info('#lines read = '+str(numberOfLinesRead))
		if totalAccessions != countedAccessions:
			logging.warning('In cluster ' + str(clusID) + ' the number of accessions read does not match the number of accessions specified in the file')
		if numberOfLinesRead <2:
			logging.error('Less than 2 lines read')

	# if file contains less than 2 lines it is not a useful input for the algorithm
	assert numberOfLinesRead >1, "The provided file contains less than 2 lines."
	
def readInput():
	"""Read the file into a list."""
	filename = args.outfolder+'/'+'formatted_input_'+os.path.basename(args.filename)
	with open(filename, "r") as f:
		
		# read the file into a list (with newlines) so that every line of the read file is one list entry
		lines_from_file = f.readlines()
		
		# strip newline from entries and store entries in a list
		lines = [lines_from_file[x].strip('\n') for x in range(len(lines_from_file))]

		# read all transformed sequences into a 2D list, so that NWA can align all of them in all combinations
		list_with_all_sequences = [lines[i:i+2] for i in range(0, len(lines), 2)]

		return list_with_all_sequences
	
def scoringFun(i,j):
	"""Scoring function, assign the values for match calculation and return the value."""
	
	# variables
	score = 0
	
	# global weight values
	global WEIGHT_MATCH_F_GROUP
	global WEIGHT_MATCH_T_GROUP
	global WEIGHT_MATCH_H_GROUP
	global WEIGHT_MATCH_X_GROUP
	global WEIGHT_MATCH_NO_FID
	global WEIGHT_MISMATCH_FIDS
	global WEIGHT_MISMATCH_NOFID_FID
	
	# split label to compare the different groups
	list_with_splitted_groups_i = i.split('.')
	list_with_splitted_groups_j = j.split('.')
	if len(list_with_splitted_groups_i)<4:
		list_with_splitted_groups_i.append('None')
	if len(list_with_splitted_groups_j)<4:
		list_with_splitted_groups_j.append('None')

	# scoring logic
	
	# 0 indicates "no labeled proteindomain" 
	if i == '0' and j == '0':
		score = WEIGHT_MATCH_NO_FID
	# match
	elif i.strip() == j.strip():
		score = WEIGHT_MATCH_F_GROUP
	# one sequence with NOFID, but not both
	elif i == '0' or j == '0':
		score = WEIGHT_MISMATCH_NOFID_FID
	#compare X-group, then H-group, ...
	elif list_with_splitted_groups_i[0] == list_with_splitted_groups_j[0]:
		if list_with_splitted_groups_i[1] == list_with_splitted_groups_j[1]:
			if list_with_splitted_groups_i[2] == list_with_splitted_groups_j[2]:
				if list_with_splitted_groups_i[3] == list_with_splitted_groups_j[3] and list_with_splitted_groups_i[3] != None:
					score = WEIGHT_MATCH_F_GROUP
				else:
						score = WEIGHT_MATCH_T_GROUP
			else:
				score = WEIGHT_MATCH_H_GROUP
		else:
			score = WEIGHT_MATCH_X_GROUP
	else:
		score = WEIGHT_MISMATCH_FIDS

	return score

def build_matrix(sequence_length1,sequence_length2):
	"""Build a matrix of the correct length for the scorematrix function, initiate it with 0 and return the matrix as 2d list."""

	# compute matrix size
	row = sequence_length2+1
	column = sequence_length1+1
	
	# initialize matrix with 0
	try:
		matrix = [[0 for column in range(column)] for row in range(row)]
	except:
		print('An error occured.')
		if args.log:
			logging.critical('unknown error in def build_matrix')
		sys.exit()
	if args.log:
		logging.debug('Matrix built')
	return matrix

def scorematrix(sequence_length1,sequence_length2,value_list1,value_list2):
	"""Build and return a scorematrix. Stores computed scores in the scorematrix with dimensions m + 1, n + 1, if m and n are the lengths of the sequences."""
	
	# for gap opening and extension
	global PENALTY_GAP_OPENING
	global PENALTY_GAP_EXTENSION
	previous_is_a_insertion = False
	previous_is_a_deletion = False
	
	matrix = build_matrix(sequence_length1,sequence_length2)
	
	# initialize indizes. i and j have to be 1 at the beginning, because M[0][0]=0
	i = 1
	j = 1
	
	#initialize matrix (+1 because of the 0 at 0,0)
	MATRIX_COLUMN_N = sequence_length2+1
	MATRIX_ROW_N = sequence_length1+1
	for i in range(MATRIX_COLUMN_N):
		matrix[i][0] = PENALTY_GAP_OPENING * i
	for j in range(MATRIX_ROW_N):
		matrix[0][j] = PENALTY_GAP_OPENING * j

	# assign new indizes
	i = 1
	j = 1
	
	# fill scorematrix
	if not args.gapextension:
		for i in range(1,MATRIX_COLUMN_N):
			for j in range(1,MATRIX_ROW_N):
				match = matrix[i-1][j-1] + scoringFun(value_list2[i-1], value_list1[j-1])
				insert = matrix[i][j-1] + PENALTY_GAP_OPENING
				delete = matrix[i-1][j] + PENALTY_GAP_OPENING
				matrix[i][j] = max(match, insert, delete)
				
	# fill scorematrix gap opening and extension
	if args.gapextension:
		for i in range(1,MATRIX_COLUMN_N):
			for j in range(1,MATRIX_ROW_N):
				match = matrix[i-1][j-1] + scoringFun(value_list2[i-1], value_list1[j-1])
				
				# check if insertion is extending a gap or not
				if previous_is_a_insertion == False:
					insert = matrix[i][j-1] + PENALTY_GAP_OPENING
				else:
					insert = matrix[i][j-1] + PENALTY_GAP_EXTENSION
					
				# check if deletion is extending a gap or not
				if previous_is_a_deletion == False:
					delete = matrix[i-1][j] + PENALTY_GAP_OPENING
				else:
					delete = matrix[i-1][j] + PENALTY_GAP_EXTENSION
				
				# write entry to matrix
				operation_dict = {"insertion":insert,"deletion":delete, "match":match}
				new_matrix_entry = max(operation_dict.values())
				matrix[i][j] = new_matrix_entry
				
				# reset value if needed 
				max_key = max(operation_dict, key=lambda k: operation_dict[k])
				if max_key == "deletion":
					previous_is_a_deletion = True
					previous_is_a_insertion = False
				elif max_key == "insertion":
					previous_is_a_insertion = True
					previous_is_a_deletion = False
				else:
					previous_is_a_deletion = False
					previous_is_a_insertion = False
	
	if args.log:
		logging.debug('Scorematrix done')		
	
	return matrix
	
def traceback(name1,name2,value_list1,value_list2, matrix):
	"""Go back through the matrix with information about the matches and align the sequences."""
	
	# initialize variables
	AlignmentA = ''
	AlignmentB = ''
	GAP_CHARACTER = '-'
	list_length1 = len(value_list1)
	list_length2 = len(value_list2)
	i = list_length2
	j = list_length1
	listIndex1 = i-1
	listIndex2 = j-1

	# compute alignment as in Needleman-Wunsch algorithm
	while i>0 or j>0:
		if i>0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + scoringFun(value_list2[listIndex1], value_list1[listIndex2]):
			AlignmentA = str(value_list2[listIndex1]) + ',' + AlignmentA
			AlignmentB = str(value_list1[listIndex2]) + ',' + AlignmentB
			listIndex1 = listIndex1-1	
			i = i - 1
			listIndex2 = listIndex2-1
			j = j - 1
		elif i>0 and (matrix[i][j] == matrix[i-1][j] + PENALTY_GAP_OPENING or matrix[i][j] == matrix[i-1][j] + PENALTY_GAP_EXTENSION):
			AlignmentA = str(value_list2[listIndex1]) + ',' + AlignmentA
			AlignmentB = GAP_CHARACTER + ',' + AlignmentB
			listIndex1 = listIndex1-1
			i = i - 1
		elif j>0 and (matrix[i][j] == matrix[i][j-1] + PENALTY_GAP_OPENING or matrix[i][j] == matrix[i][j-1] + PENALTY_GAP_EXTENSION):
			AlignmentA = GAP_CHARACTER + ',' + AlignmentA
			AlignmentB = str(value_list1[listIndex2]) + ',' + AlignmentB
			j = j - 1
			listIndex2 = listIndex2-1
		else:
			print('error in traceback')
			if args.log:
				logging.critical('error in traceback')
			quit()
	
	# prepare Strings for writing
	AlignmentA = AlignmentA.strip(',')
	AlignmentB = AlignmentB.strip(',')
	score_alignment = str(matrix[list_length2][list_length1])
	
	# build lists for verbose output and length computation
	list_of_AlignmentA = AlignmentA.split(',')
	list_of_AlignmentB = AlignmentB.split(',')
	alignmentlength = len(list_of_AlignmentA)	

	# build a string with names, score and length of alignment, set filename
	first_line = name1 + ',' + name2 + ',' + score_alignment + ',' + str(alignmentlength)
	filename = args.outfolder+'/'+'alignments_'+os.path.basename(args.filename)
	
	# write 3 lines in an outputfile (for every alignment)
	with open(filename, 'a', encoding="utf-8") as file_alignments:
		write_data = file_alignments.write(first_line + '\n' + AlignmentB + '\n' + AlignmentA + '\n')
	
	# check if option -verbose is set and if so, create the second output file
	if args.verbose:
		writeHumanreadableOutput(name1,name2,score_alignment,alignmentlength,list_of_AlignmentA,list_of_AlignmentB)

def needleman_wunsch(subset):
	"""Function with all calls for NWA, except selfalignments."""
	
	# global variables
	global number_of_alignments
	
	#input variables
	name1 = subset[0][0]
	name2 = subset[1][0]
	value_string1 = subset[0][1]
	value_string2 = subset[1][1]
	value_list1 = value_string1.split(',')
	value_list2 = value_string2.split(',')
	
	sequence_length1 = len(value_list1)
	sequence_length2 = len(value_list2)

	# build scorematrix for this alignment
	matrix = scorematrix(sequence_length1,sequence_length2,value_list1,value_list2)
	
	# compute traceback for this alignment
	matrix2 = traceback(name1,name2,value_list1,value_list2, matrix)
	
	# global count for logging
	number_of_alignments = number_of_alignments+1
	
def needleman_wunsch_self(entry):
	"""Function with all selfalingment calls for NWA."""
	
	# global variables
	global number_of_alignments
	
	#input variables
	name1 = entry[0]
	name2 = name1
	value_string = entry[1]
	value_list1 = value_string.split(',')
	value_list2 = value_list1
	
	sequence_length1 = len(value_list1)
	sequence_length2 = sequence_length1

	# build scorematrix for this alignment
	matrix = scorematrix(sequence_length1,sequence_length2,value_list1,value_list2)
	
	# compute traceback for this alignment
	matrix2 = traceback(name1,name2,value_list1,value_list2, matrix)
	
	# global count for logging
	number_of_alignments = number_of_alignments+1
	
def alignAllPairs(list_with_all_sequences):
	"""Calls NWA for every possible pair of sequences."""
	
	for subset in itertools.combinations(list_with_all_sequences, 2):
		needleman_wunsch(subset)
	
	# if it is not explicitly stated that selfalignments should not be calculated, calculate selfalignments.
	if not args.noselfalignment:
		if args.log:
			logging.info('Compute selfalignments.')
		for entry in list_with_all_sequences:
			needleman_wunsch_self(entry) #aktuell hängt traceback das einfach hintendran
	
def writeHumanreadableOutput(name1,name2,score_alignment,alignmentlength,list_of_AlignmentA,list_of_AlignmentB): #TODO
	"""Writes a more human readable output if -v is active and inserts padding with underscore."""
	
	# variables
	entry_splitted=[]
	transformedListA=[0]*alignmentlength
	transformedListB=[0]*alignmentlength
	iterate_A = 0
	iterate_B = 0
	
	# add padding to both Strings
	for entry in list_of_AlignmentB:
		if entry == '0':
			transformedListB[iterate_B] = "0000.0000.0000.0000"
		elif entry == "-":
			transformedListB[iterate_B] = "----.----.----.----"
		else:
			entry_splitted = entry.split('.')
			splitList = [0]*len(entry_splitted)
			j = 0
			for label in entry_splitted:
				splitList[j] = label.ljust(4, '_')
				j = j+1
				
			a_string = '.'.join(splitList)
			if len(splitList)==3:
				a_string = a_string + '.____'
			elif len(splitList)==2:
				a_string = a_string + '.____.____'
			elif len(splitList)==1:
				a_string = a_string + '.____.____.____'
			transformedListB[iterate_B] = a_string
			entry_splitted=[]
		iterate_B = iterate_B +1
		
	for entry in list_of_AlignmentA:
		if entry == '0':
			transformedListA[iterate_A] = "0000.0000.0000.0000"
		elif entry == "-":
			transformedListA[iterate_A] = "----.----.----.----"
		else:
			entry_splitted = entry.split('.')
			splitList = [0]*len(entry_splitted)
			j = 0
			for label in entry_splitted:
				if len(label)>4:
					logging.warning("Found a label with len >4, code asserts len 4, please update code accordingly.")
				splitList[j] = label.ljust(4, '_')
				j = j+1
				
			a_string = '.'.join(splitList)
			if len(splitList)==3:
				a_string = a_string + '.____'
			elif len(splitList)==2:
				a_string = a_string + '.____.____'
			elif len(splitList)==1:
				a_string = a_string + '.____.____.____'
			transformedListA[iterate_A] = a_string
			entry_splitted=[]
		iterate_A = iterate_A +1				

	# write sequences as strings
	hAlignmentA = ','.join(map(str, transformedListA))
	hAlignmentB = ','.join(map(str, transformedListB))
	
	# build annotation string with names and score
	first_line = name1 + ',' + name2 + ',' + score_alignment
	
	# write data in an output file
	with open(args.outfolder+'/'+'file_alignments_verbose_'+os.path.basename(args.filename), 'a', encoding="utf-8") as file_alignments:
		write_data = file_alignments.write(first_line + '\n' + hAlignmentB + '\n' + hAlignmentA + '\n')
	
	if args.log:
		logging.debug('second file with more humanreadable output created')
		
def inputChecks():
	"""Check if given file exists, is a file and not empty. Check if output files already existing."""
	
	# check if given file exists
	fileExists = os.path.exists(args.filename)
	if not fileExists:
		logging.error('The input file does not exist.')
	assert fileExists, "The provided file does not exist."
	
	# check if given path points to a file
	pathIsFile = os.path.isfile(args.filename)
	if not pathIsFile:
		logging.error('The specified path does not point to a file.')
	assert pathIsFile, "The specified path does not point to a file."
	
	# check if given file is not empty
	if (os.stat(args.filename).st_size == 0):
		logging.error('The file is empty.')
	assert os.stat(args.filename).st_size > 0, "The file is empty."
	
	# check if output file aleady exists from previous run and if so, set warning
	filename_output = 'alignments_'+os.path.basename(args.filename)
	fileExists = os.path.exists(filename_output)
	if fileExists:
		logging.warning('The alignment file already existed from a previous run.')
	
	# check if output file aleady exists from previous run and if so, set warning
	filename_intermediate = 'formatted_input_'+os.path.basename(args.filename)
	fileExists = os.path.exists(filename_intermediate)
	if fileExists:
		logging.warning('The intermediate file already existed from a previous run.')
	
def deleteTempFiles():
	"""Delete intermediate file."""
	
	if os.path.exists(args.outfolder+'/'+'formatted_input_'+os.path.basename(args.filename)):
		os.remove(args.outfolder+'/'+'formatted_input_'+os.path.basename(args.filename))
	else:
		print('No file found to delete')
		if args.log:
			logging.error('The intermediate file could not be found.')

def main():
	"""Main function which calls all functions"""
	try:
		# check if input file is okay and check if output files already existing from previous run
		inputChecks()
		
		# transform the given input in proper input for needleman wunsch algorithm
		transformRawInput()
		
		# read input data and align all sequences
		alignAllPairs(readInput())

		# check if option -temp is set and if so, delete file with transformed input
		if not args.temp:
			deleteTempFiles()
		
		# log info with gapextension value
		if args.log:
			if args.gapextension:
				global PENALTY_GAP_EXTENSION
				global PENALTY_GAP_OPENING
				infostring = 'Different penalties for gap extensions and openings are used. Gap extension penalty is ' + str(PENALTY_GAP_EXTENSION) + 'and gap opening penalty is ' + str(PENALTY_GAP_OPENING) + '.'
				logging.info(infostring)
			
	except:
		print("unknown error occurred")
	
# call main function to start program
main()
