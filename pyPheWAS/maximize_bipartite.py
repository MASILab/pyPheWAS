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
		if '%s' in outputfile:
			outputfile = outputfile % ('max')
	else:
		print("Matched data 1-%s" % (matched))
		if '%s' in outputfile:
			outputfile = outputfile % (matched)

	new_data.to_csv(path + outputfile, index=False)
	print("New group file in %s" % (path + outputfile))


def control_match(path, inputfile, outputfile, keys, deltas, condition='genotype', goal=1):
	# Reformat arguments into Python format
	keys = keys.replace(" ", "").split('+')
	deltas = deltas.replace(" ", "").split(',')
	deltas = [CATEGORICAL_DATA if x == '' else float(x) for x in deltas]

	# save original goal value
	orig_goal = goal

	# Read data from the provided input file
	data = pd.read_csv(path + inputfile)

	# Assert that all of the provided keys are present in the data
	for key in keys:
		assert key in data.columns, '%s is not a column in the input file (%s)' % (key, inputfile)
	# Assert that condition column is present in the data
	assert condition in data.columns, 'Specified condition (%s) is not a column in the input file (%s)' % (condition, inputfile)
	# Assert that condition column contains only '1' and '0'
	condition_vals = np.unique(data[condition])
	assert len(condition_vals) == 2, 'There are %d values (should only be 2) in the specified condition column (%s) in the input file (%s)' % (len(condition_vals), condition, inputfile)
	for val in [0, 1]:
		assert val in condition_vals, 'The value %d is missing from the condition column (%s) in the input file (%s)' % (val, condition, inputfile)

	# Assign new value for outputfile
	if not outputfile:
		outputfile = '1-%s_' + inputfile

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
		targets[match_col] = targets.apply(expand, axis=1)[i]
		for key in keys:
			match_info = pd.merge(targets[[match_col]], data[['id', key]], left_on=match_col, right_on='id')
			match_info.rename(columns={key: match_col + '_' + key}, inplace=True)
			targets = pd.merge(targets, match_info.drop(columns='id'), on=match_col,how='left')
			cols.append(match_col + '_' + key)

	cols.insert(0,'id')
	cols.insert(1, 'genotype')
	# export matching pairs
	print('Saving case/control mapping to ' + path + 'matches__' + outputfile)
	targets.to_csv(path + 'matches__' + outputfile,index=False, columns=cols)

	if final_ratio == orig_goal:
		matching_success = 1
	else:
		matching_success = 0

	output_matches(path, outputfile, data, all_used, matching_success, orig_goal, final_ratio)


