import pandas as pd
import operator
import random
import numpy as np
import sys
import getopt
"""


"""
CATEGORICAL_DATA = '675161f1c87ff2648c61ff1c57c780f2'


def generate_row_query(keys, deltas, tr):
	q = []
	for i,dt in enumerate(deltas):
		key = keys[i]
		is_categorical = dt == CATEGORICAL_DATA
		if is_categorical:
			part = '=='.join([key, tr[key].__repr__()])
		else:
			structure = ['abs(', key, '-',tr[key],')', '<', dt]
			part = ''.join([str(x) for x in structure])
		q.append(part)
	return '&'.join(q)

def get_options(targets, controls, keys, deltas):
	tt = targets[keys]
	c = controls[keys]
	matching = {}
	if len(c) > len(tt):
		for i in tt.index:
			tr = tt.loc[i]
			control_query = generate_row_query(keys, deltas, tr)
			matches = c.query(control_query).index
			matching[i] = matches.drop_duplicates().tolist()
	else:
		for i in c.index:
			tr = c.loc[i]
			target_query = generate_row_query(keys, deltas, tr)
			matches = tt.query(target_query).index
			matching[i] = matches.drop_duplicates().tolist()
	return matching

def generate_matches(matching, goal):
	# Sort the targets by the number of controls they match
	frequency = { k : len(v) for k,v in matching.items() }
	frequency = sorted(frequency.items(), key=operator.itemgetter(1))
	success = True

	# Keep track of the already used controls
	used = []

	# The final mapping of targets : [control list]
	final = {}

	for key,_ in frequency:
		final[key] = []
		viable = matching[key]
		random.shuffle(viable)
		for control in viable:
			if len(final[key]) == goal:
				break
			if control not in used:
				used.append(control)
				final[key].append(control)
		if len(final[key]) < goal:
			success = False
	return (final, used, success, goal)

def maximize_matches(matching):
	prev = generate_matches(matching, 1)
	while prev[2] == False:
		return prev

	# If 1-1 matching was successful, attempt to maximize starting from 2
	success = prev[2]
	goal = 2

	while success:
		curr = generate_matches(matching, goal)
		success = curr[2]
		if success:
			prev = curr
			goal += 1
	
	return prev

def output_matches(path, outputfile, data, all_used, success, matched):
	new_data = data[data.index.isin(all_used)]

	if not success:
		print("Could not match 1-1, using the maximum number of matches found by the approximation algorithm")
		if '%s' in outputfile:
			outputfile = outputfile % ('max')
	else:
		print("Matched data 1-%s" % (matched))
		if '%s' in outputfile:
			outputfile = outputfile % (matched)

	new_data.to_csv(path + outputfile)
	print("Data in %s" % (path + outputfile))

def control_match(path, inputfile, outputfile, keys, deltas, condition='genotype',goal=-1):
	# Reformat arguments into Python format
	keys = keys.split('+')
	deltas = deltas.split(',')
	deltas = [CATEGORICAL_DATA if x == '' else int(x) for x in deltas]

	# Read data from the provided input file
	data = pd.read_csv(path + inputfile)

	# Assert that all of the provided keys are present in the data
	for key in keys:
		assert key in data.columns, '%s not a column in the input file (%s)' % (key, inputfile)

	# Assign new value for outputfile
	if not outputfile:
		outputfile = '1-%s_' + inputfile

	# Separate patients and controls
	targets = data[data[condition] == 1]
	controls = data[data[condition] == 0]

	match_by_control = len(targets) > len(controls)

	matching = get_options(targets, controls, keys, deltas)
	if goal != -1:
		final, used, success, matched = generate_matches(matching, goal)
		if success:
			if match_by_control:
				all_used = used + controls.index.tolist()
			else:
				all_used = used + targets.index.tolist()
			output_matches(path, outputfile, data, all_used, success, matched)
			return
		else:
			print("Failed to perform 1-%s, attempting to maximize..." % (goal))
			while not success:
				goal = 1
				print(deltas)
				deltas = [element + 1 if element != CATEGORICAL_DATA else element for element in deltas]
				matching = get_options(targets, controls, keys, deltas)
				final, used, success, matched = generate_matches(matching, goal)
			print("Used %s as delta values across keys. Generated a 1-%s match." % (deltas, goal))
	final, used, success, matched = maximize_matches(matching)
	if match_by_control:
		all_used = used + controls.index.tolist()
	else:
		all_used = used + targets.index.tolist()
	output_matches(path, outputfile, data, all_used, success, matched)
	if goal==-1:
		final, used, success, matched = maximize_matches(matching)
		#all_used = used + targets.index.tolist()
		if match_by_control:
			all_used = used + controls.index.tolist()
		else:
			all_used = used + targets.index.tolist()
		output_matches(path, outputfile, data, all_used, success, matched)

