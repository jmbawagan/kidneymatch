'''
	Author:
	Juan Miguel Bawagan
	Institute of Computer Science, UPLB

	Creates a solution (or solutions) from a csv data file. 
	See create_test_data.py for the format.

	(1)
	For multi-way matching, we'll treat the problem 
	as an assignment problem.
	We'll use the Hungarian algorithm to solve it.

	(2)
	If a strictly pair-wise solution is needed, we'll use maximum weight
	matching using the Blossom algorithm for finding augmenting paths
	and the primal-dual method for finding a matching of maximum weight
	(both invented by Jack Edmonds).
'''
#!/usr/bin/env python

import PySimpleGUI as sg
import sys
import csv
import numpy as np
from scipy.optimize import linear_sum_assignment
import networkx as nx
import networkx.algorithms.matching as matching



def read_csv(filename):
	''' Reads the csv file and returns a patient list, donor list,
	the bipartite matrix, the inverted bipartite matrix, and the 
	maximum score.
	'''
	
	file_handle = open(filename, "r")

	# read first line, patient list information
	patient_list = file_handle.readline()[:-1].split(",")
	patient_list.pop(0)
	n = len(patient_list)

	# create the max_score, donor list, and matrices
	max_score = 0
	donor_list = []
	matrix = np.zeros((n,n), dtype=int)
	inverted = np.zeros((n,n), dtype=int)

	# next n lines, donor and matrix information
	for i in range(n):
		temp_list = file_handle.readline()[:-1].split(",")
		donor_list.append(temp_list[0])

		# put values into matrix
		for j in range(1,n+1):
			#check if it's a blank, if yes, make it zero
			if temp_list[j] == "":
				temp_list[j] = 0
			
			#typecast to int temporarily
			matrix[i][j-1] = int(temp_list[j])

			if matrix[i][j-1] > max_score:	#get max_score while at it
				max_score = matrix[i][j-1]

		# create inverted matrix for hungarian method
		for j in range(1, n+1):
			inverted[i][j-1] = max_score - int(temp_list[j])

	file_handle.close()

	return patient_list, donor_list, matrix, inverted, max_score

def pairwise_match(matrix, max_score, regular=True, zero_out=False):
	''' Creates an complete graph G from the given donor-patient matrix.
	
	Each original donor-patient pair is set as a node with edges going to
	other donor-patient pairs. Edge weights are derived from the average score
	of the pairwise exchange or our modified average.
	score = AVG(value1,value2) * (MAXPOSSIBLESCORE - ABS_DIFF(value1,value2))
	Max weight matching is then done to get the pairwise matching.

	Returns the row_ind and col_ind lists indicating the matching
	(Similar to linearsumassignment).
	'''
	n = matrix.shape[0]
	
	# create graph
	# C(6,2) [6 taken 2] no of edges
	# ignore self loops for now
	G = nx.Graph()
	for i in range(0,n-1): 
		for j in range(i+1,n):
			if zero_out and (matrix[i][j] == 0 or matrix[j][i] == 0):
				# if zero out is True, we automatically make the score
				# zero if there is one incompatible pair
				# score = 0
				# EDIT:
				# we now instead skip adding the edge
				# instead of adding an edge with zero weight
				continue
			else:
				average = (matrix[i][j] + matrix[j][i]) / 2
				if regular:
					# regular average
					score = average
				else:
					# our modified average
					score = average * (max_score - abs(matrix[i][j]-matrix[j][i]))

			G.add_edge(i,j,weight=score)

	M = matching.max_weight_matching(G,maxcardinality=True)
	
	# initialize row and col indices list
	row_ind = [x for x in range(n)]
	col_ind = [-1 for x in range(n)]	# -1 for unmatched
	for pair in M:
		col_ind[pair[0]] = pair[1]
		col_ind[pair[1]] = pair[0]
		

	return row_ind, col_ind

def print_mapping(matrix, donor_ind, patient_ind, donor_list, patient_list, filename):
	''' Writes a .csv files of the donor-patient mappings.
	Format is donor,patient,score.
	'''
	# This can be modified later.
	# file_handle = open(filename+".csv", "w")
	file_handle = open(filename, "w")

	# file_handle.write(",")
	for i in range(len(donor_ind)):
		# donor,patient,score
		row = donor_ind[i]
		col = patient_ind[i]
		file_handle.write(str(donor_list[row]) + ",")
		if col == -1:
			file_handle.write("NONE,")
			file_handle.write("0\n")
		else:
			file_handle.write(str(patient_list[col]) + ",")
			file_handle.write(str(matrix[row][col]) + "\n")
			

	file_handle.close()


sg.ChangeLookAndFeel('TealMono')	# a background bug occurs when 
									# put after the layout
sg.SetOptions(element_padding=(0, 0))


# initialize empty table
data = [['donor1'," "," "," "," "," "," "," "," "],
		['donor2'," "," "," "," "," "," "," "," "],
		['donor3'," "," "," "," "," "," "," "," "],
		['donor4'," "," "," "," "," "," "," "," "],
		['donor5'," "," "," "," "," "," "," "," "],
		['donor6'," "," "," "," "," "," "," "," "],
		['donor7'," "," "," "," "," "," "," "," "],
		['donor8'," "," "," "," "," "," "," "," "]]
#header_list = ['','patient1','patient2','patient3','patient4','patient5','patient6','patient7','patient8']
header_list = ['','1','2','3','4','5','6','7','8']

data2 = [[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]]

# layout
layout = [	#[sg.Button('Load CSV File')],
			[sg.Text('File:'),
     		sg.InputText('',disabled=True,key='_inputtext_'),
			sg.FileBrowse('Browse CSV File', 
				tooltip=None,
				key="_file_",
				file_types=(("CSV Files","*.csv"),)),
			sg.Open(key="_open_", disabled=False)],
			[sg.Text(' ' * 30,font=("Helvetica", 2))],	# spacing

			[sg.Text("Score matrix")],
			[sg.Table(values=data,
				headings=header_list,
				key="_table1_",
				def_col_width=5,
				max_col_width=30,
				auto_size_columns=False,
				justification='right',
				alternating_row_color='#d7ebfa',
				#vertical_scroll_only=False,
				hide_vertical_scroll=True,
				#display_row_numbers=True,
				num_rows=min(len(data), 20))],


			[sg.Text(' ' * 30,font=("Helvetica", 2))],	# spacing
			[sg.Frame('Choose a matching option',
				layout=[
					[sg.Radio('Pair-wise matching', 
						"OPTION", 
						default=True,
						key="_pairwise_",
						change_submits=True), 
					sg.Checkbox('Use average', 
						default=True, 
						key="_average_", 
						change_submits=True),
					sg.Checkbox("Don't add edges with zero weight", 
						default=False, 
						key="_zero_", 
						change_submits=True)],
					[sg.Radio('Multi-way matching', 
						"OPTION", 
						key="_multi-way_",
						change_submits=True)]
				],
				title_color='dark blue', 
				#relief=sg.RELIEF_SUNKEN, 
				#border_width=0,
				tooltip=None)],
			[sg.Text(' ' * 30,font=("Helvetica", 2))], 	# spacing
			[sg.Button('RUN MATCHING')],


			[sg.Text('_' * 70,font=("Helvetica", 10))], # line
			[sg.Text(' ' * 30,font=("Helvetica", 2))], 	# spacing

			# RESULTS
			[sg.Text("RESULTS", 						
				key='_results_text',
				visible=False)],		
			[sg.Table(values=data2,
						headings=["Donor", "Patient", "Score"],
						key='_results_',
						def_col_width=7,
						max_col_width=30,
						auto_size_columns=False,
						justification='right',
						alternating_row_color='#d7ebfa',
						#vertical_scroll_only=False,
						hide_vertical_scroll=True,
						visible=False,
						#display_row_numbers=True,
						num_rows=min(len(data2), 20))],


			[sg.Text("Total matching score:", key='_score_', size=(30, 1), visible=False)],
			[sg.Text(' ' * 30,font=("Helvetica", 2))], # spacing
			
			#[sg.Button('Save results to file',key='_save_button_',visible=False)],

			[sg.Frame('',
				layout=[[
					#sg.Text('Save to:'),
					sg.SaveAs('Save results to..',
						key='_saveas_button_',
						#enable_events=True,
						target='_savelocation_',
						file_types=(("CSV Files","*.csv"),)),
					sg.InputText('',disabled=True,key='_savelocation_',size=(42, 1)),
					sg.Button('Save',size=(8, 1))]
				],
				border_width=0,
				key='_saveframe_',
				visible=False,
				tooltip=None)],

 
			[sg.Text(' ' * 30,font=("Helvetica", 2))], # spacing
			[sg.Text("Note: Due to GUI limitations, only 8 columns can be shown. Don't worry, everything still works!",font=("Helvetica", 8))]
			
			
]


# window
window = sg.Window("Kidney matching program v0.01").Layout(layout)


# initializations
file_open_flag = False
average_flag = True
zero_out_flag = False
matching_option = "pairwise"
results_flag = False

# event loop
while True:
	event, values = window.Read()
	
	#print(event)

	if event is None:
		break
	
	# -- set average or the special average --
	if event == '_average_':
		average_flag = not average_flag

	# -- set normal pair matching or zero_out --
	if event == '_zero_':
		zero_out_flag = not zero_out_flag

	# -- set matching option --
	if event == "_pairwise_":
		matching_option = "pairwise"

	if event == "_multi-way_":
		matching_option = "multi-way"

	# -- open CSV file --
	# from demo_table_csv.py
	if event == '_open_':
		fname = values['_file_']
		if fname == "":
			sg.Popup("Browse a file first!")
		else:
			#print("open:",values['_file_'])		
			button = sg.PopupYesNo('Does this file have column names already?')
			if fname is not None:
				# put into the Score Matrix
				with open(fname, "r") as infile:
					reader = csv.reader(infile)
					if button == 'Yes':
						header_list = next(reader)
					try:
						data = list(reader)  # read everything else into a list of rows
						if button == 'No':
							header_list = ['column' + str(x) for x in range(len(data[0]))]

						# read the input csv file
						patient_list, donor_list, matrix, inverted, max_score = read_csv(fname)
						file_open_flag = True
					except:
						sg.PopupError("Error reading, check if it's in the proper format!")
						file_open_flag = False
						#sys.exit(69)
			
			if file_open_flag == True:
				window.FindElement('_table1_').Update(num_rows=min(len(data),20))
				window.FindElement('_table1_').Update(values=data)  


	# -- run matching --
	if event == 'RUN MATCHING':
		if file_open_flag == False:
			sg.Popup("Open a CSV file first!")
		else:
			if matching_option == "pairwise":
				# pair-wise matching
				'''
				if average_flag == True:
					# regular average
					row_ind, col_ind = pairwise_match(matrix,max_score)
				else:
					#modified average
					row_ind, col_ind = pairwise_match(matrix,max_score,regular=False)
				'''
				row_ind, col_ind = pairwise_match(matrix,max_score,regular=average_flag,zero_out=zero_out_flag)
			else:
				#multi-way matching
				row_ind, col_ind = linear_sum_assignment(inverted)
				
			#print(row_ind,col_ind)	
			
			# put into results table!
			results_matrix = []
			total_match_score = 0
			for i in range(len(row_ind)):
				pair_list = []
				# donor,patient,score
				row = row_ind[i]
				col = col_ind[i]
				pair_list.append(donor_list[row])
				if col == -1:
					pair_list.append("NONE")
					pair_list.append(0)
				else:
					pair_list.append(patient_list[col])
					pair_list.append(matrix[row][col])
					total_match_score += matrix[row][col]
				
				results_matrix.append(pair_list)

			window.FindElement('_results_').Update(num_rows=min(len(row_ind),20))
			window.FindElement('_results_').Update(values=results_matrix)
			#window.FindElement('_score_').Update(("Total matching score: "+str(total_match_score)))
			new_text = "Total matching score: " + str(total_match_score)
			window.FindElement('_score_').Update(new_text)

			# make results elements visible
			window.FindElement('_results_text').Update(visible=True)
			window.FindElement('_results_').Update(visible=True)
			window.FindElement('_score_').Update(visible=True)
			window.FindElement('_saveframe_').Update(visible=True)

			results_flag = True
	
	
	# -- save results to file -- 
	if event == 'Save':
		if results_flag == False:
			sg.Popup("Run a matching first!")
		else:
			# write to file
			fname = values['_savelocation_']
			
			if fname[-4:len(fname)] != ".csv":
			 	fname = fname + '.csv'

			window.FindElement('_savelocation_').Update(fname)
			print_mapping(matrix, row_ind, col_ind, donor_list, patient_list, fname)
			sg.Popup("File saved!")
			# print(fname)
