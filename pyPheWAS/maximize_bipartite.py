import pandas as pd
import operator
import random
import numpy as np
import sys
import getopt
from hopcroftkarp import HopcroftKarp
from tqdm import tqdm

"""


"""
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

    if not success:
        print("Could not match 1-%d, using the maximum number of matches found by the approximation algorithm" %goal)
        print("Matched data 1-%0.3f" %matched)
        if '%s' in outputfile:
            outputfile = outputfile % ('max')
    else:
        print("Matched data 1-%s" % (matched))
        if '%s' in outputfile:
            outputfile = outputfile % (matched)

    new_data.to_csv(path + outputfile, index=False)
    print("Data in %s" % (path + outputfile))


def control_match(path, inputfile, outputfile, keys, deltas, condition='genotype',goal=1):
    # Reformat arguments into Python format
    keys = keys.replace(" ","").split('+')
    deltas = deltas.replace(" ","").split(',')
    deltas = [CATEGORICAL_DATA if x == '' else float(x) for x in deltas]

    # save original goal value
    orig_goal = goal

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

    if match_by_control: # more targets than controls
        print('There are more targets (1) than controls (0) -- matching by controls')
        cid = set()
        tid = set()
        set_num = 0
        while goal>0:
            set_num += 1
            print('Getting match set %d' %set_num)
            matching = get_options(targets, controls, keys, deltas)
            matched = HopcroftKarp(matching).maximum_matching()
            for i in tqdm(controls.index,total=len(controls),desc="Checking matches"):
                if i in matched:
                    tid.add(matched[i])
                    cid.add(i)
            rem_ids = set(targets.index).difference(tid)
            targets=targets.ix[rem_ids]
            goal=goal-1
        final_ratio = float(len(tid)) / float(len(cid))
    else: # more controls than targets
        print('There are more controls (0) than targets (1) -- matching by targets')
        cid = set()
        tid = set()
        set_num = 0
        while goal>0:
            set_num += 1
            print('Getting match set %d' % set_num)
            matching = get_options(targets, controls, keys, deltas)
            matched = HopcroftKarp(matching).maximum_matching()
            for i in tqdm(targets.index,total=len(targets),desc="Checking matches"):
                if i in matched:
                    cid.add(matched[i])
                    tid.add(i)
                else:
                    print(data.loc[i,'id'])
            rem_ids = set(controls.index).difference(cid)
            controls=controls.ix[rem_ids]
            goal=goal-1
        final_ratio = float(len(cid)) / float(len(tid))
    all_used=cid.union(tid)

    if final_ratio == orig_goal:
        matching_success = 1
    else:
        matching_success = 0

    output_matches(path, outputfile, data, all_used, matching_success, orig_goal,final_ratio)


