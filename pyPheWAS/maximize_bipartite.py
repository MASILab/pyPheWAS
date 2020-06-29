"""
**Case/Control Matching functions**

Contains all functions that implement the :ref:`maximizeControls` data preparation tool.
"""

import pandas as pd
from hopcroftkarp import HopcroftKarp
import numpy as np

CATEGORICAL_DATA = '675161f1c87ff2648c61ff1c57c780f2'


def generate_row_query(keys, deltas, tr):
	q = []
	for i, dt in enumerate(deltas):
		key = keys[i]
		is_categorical = dt == CATEGORICAL_DATA
		if is_categorical:
			part = '=='.join([key, tr[key].__repr__()])
		else:
			structure = ['abs(', key, '-', tr[key], ')', '<=', dt]
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
			# matching[i] = matches.drop_duplicates().tolist()
			matching[i] = set(matches)
	else:
		for i in c.index:
			tr = c.loc[i]
			target_query = generate_row_query(keys, deltas, tr)
			matches = tt.query(target_query).index
			# matching[i] = matches.drop_duplicates().tolist()
			matching[i] = set(matches)

	return matching


def output_matches(path, outputfile, data, all_used, success, goal, matched):
	new_data = data[data.index.isin(all_used)]
	print('---')

	if not success:
		print("Could not match data 1-%d, using the maximum number of matches found by the approximation algorithm" % goal)
		print("Matched data 1-%0.3f" % matched)
	else:
		print("Matched data 1-%s" % (matched))

	new_data.to_csv(path / outputfile, index=False)
	print("New group file in %s" % (path / outputfile))

	return


def control_match(path, input, output, keys, deltas, condition='genotype', goal=1):
	"""
	Estimate an optimal case-control mapping.

	Match cases/controls (defined by value of ``condition``\ ) with a ``goal``\ :1 ratio
	based on similarity of ``keys``\ . Specify tolerance interval for ``keys`` via ``deltas``\ .
	For an exact match, specify a delta of 0. The order of
	``deltas`` values must match the order of the ``keys``\ .

	Optimal matches are estimated via the HopcroftKarp algorithm. A new group file
	with all matched subjects is saved to ``output``\ . The explicit matching of
	cases to controls is saved to ``output``\ __matched_pairs.csv

	:param path: path to input and output
	:param input: name of file containing group data
	:param output: name of output matched group file
	:param keys: comma-separated list of matching criteria (columns in input file)
	:param deltas: comma-separated list of tolerance intervals for the matching criteria
	:param condition: Field denoting group assignments [default: genotype]
	:param goal: n, target matching ratio (control:case => n:1) [default: 1]

	:type path: pathlib Path
	:type input: str
	:type output: str
	:type keys: str
	:type deltas: str
	:type condition: str
	:type goal: int

	:returns: None

	.. note:: If the matching algorithm cannot achieve the specified matching ratio, it will issue a warning and report the achieved ratio.

	"""
	# Reformat arguments
	keys = keys.replace(" ", "").split(',')
	deltas = deltas.replace(" ", "").split(',')
	deltas = [CATEGORICAL_DATA if x == '' else float(x) for x in deltas]

	# save original goal value
	orig_goal = goal

	# Assign new value for outputfile
	if output is None:
		inname = input.strip('.csv')
		output = inname + '__matched.csv'
		match_file = path / (inname + '__matched_pairs.csv')
	else:
		outname = output.strip('.csv')
		match_file = path / (outname + '__matched_pairs.csv')

	# Read data from the provided input file
	data = pd.read_csv(path / input)

	# Assert that all of the provided matching keys are present in the data
	for key in keys:
		assert key in data.columns, '%s is not a column in the input file (%s)' % (key, input)
	# Assert that condition column is present in the data
	assert condition in data.columns, 'Specified condition (%s) is not a column in the input file (%s)' % (condition, input)
	# Assert that condition column contains only '1' and '0'
	condition_vals = np.unique(data[condition])
	assert len(condition_vals) == 2, 'There are %d values (should only be 2) in the specified condition column (%s) in the input file (%s)' % (len(condition_vals), condition, input)
	for val in [0, 1]:
		assert val in condition_vals, 'The value %d is missing from the condition column (%s) in the input file (%s)' % (val, condition, input)


	# Separate patients and controls
	match_by_group0 = len(data[data[condition] == 1]) > len(data[data[condition] == 0])
	if match_by_group0:
		print('There are more cases (%s=1) than controls (%s=0) -- matching by controls' %(condition, condition))
		targets = data[data[condition] == 0].copy()
		controls = data[data[condition] == 1].copy()
	else:
		print('There are more controls (%s=0) than cases (%s=1) -- matching by cases' %(condition, condition))
		targets = data[data[condition] == 1].copy()
		controls = data[data[condition] == 0].copy()
	# save original number of targets (used to check ifany are dropped)
	orig_num_target = targets.shape[0]

	# create dictionary to store matching pairs
	targets['matching_ix'] = [[] for _ in range(targets.shape[0])]
	pairing = targets[['matching_ix']].to_dict(orient='index')

	cid = set()
	tid = set()
	set_num = 0
	while goal > 0:
		set_num += 1
		print('Getting match set %d' % set_num)
		matching = get_options(targets, controls, keys, deltas) # get all possible matches for each target
		matched = HopcroftKarp(matching).maximum_matching() # find optimal pairings
		# store matches
		for i in targets.index:
			if i in matched:
				cid.add(matched[i])
				tid.add(i)
				pairing[i]['matching_ix'].append(matched[i])
		# remove matched IDs from control pool
		rem_ids = set(controls.index).difference(cid)
		controls = controls.loc[rem_ids,:]
		goal = goal - 1

	final_ratio = float(len(cid)) / float(len(tid))
	all_used = cid.union(tid)

	print('Formatting matches')
	# convert pairing to dataframe & get matching ids from ix
	pairing_df = pd.DataFrame.from_dict(pairing, orient='index')
	# copy matching ix and remove extra columns
	targets.loc[pairing_df.index, 'matching_ix'] = pairing_df['matching_ix']
	# get match IDs from index
	get_ids = lambda i: list(data.loc[i['matching_ix'], 'id'].values)
	targets['matching_ids'] = targets.apply(get_ids, axis=1)
	# separate list of subject IDs & get matching subject's info
	cols = keys[:]
	for i in range(0,orig_goal):
		match_col = 'match'+str(i+1)
		cols.append(match_col)
		expand = lambda x: pd.Series(x['matching_ids'])
		try:
			targets[match_col] = targets.apply(expand, axis=1)[i]
			for key in keys:
				match_info = pd.merge(targets[[match_col]], data[['id', key]], left_on=match_col, right_on='id')
				match_info.rename(columns={key: match_col + '_' + key}, inplace=True)
				targets = pd.merge(targets, match_info.drop(columns='id'), on=match_col,how='left')
				cols.append(match_col + '_' + key)
		except: # no matches found for set i
			targets[match_col] = np.nan

	cols.insert(0, 'id')
	cols.insert(1, 'genotype')
	# export matching pairs
	print('Saving case/control mapping to %s' %match_file)
	targets.to_csv(match_file,index=False, columns=cols)

	if len(tid) != orig_num_target:
		g = 'controls' if match_by_group0 else 'cases'
		print('WARNING: Some %s were dropped during matching (%d out of %d %s remain)' %(g,len(tid),orig_num_target,g))

	if final_ratio == orig_goal:
		matching_success = 1
	else:
		matching_success = 0

	output_matches(path, output, data, all_used, matching_success, orig_goal, final_ratio)

	return
