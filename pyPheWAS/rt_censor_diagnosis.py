import pandas as pd
import numpy as np

def censor_diagnosis(genotype_file, phenotype_file, final_pfile, final_gfile, efield, delta_field, start_time=np.nan, end_time=np.nan):
	# read files & check field names
	print('Reading input files')
	genotypes = pd.read_csv(genotype_file)
	phenotypes = pd.read_csv(phenotype_file)
	mg = pd.merge(phenotypes, genotypes, on='id')
	# assert specified fields exist
	assert efield in phenotypes, 'Specified efield (%s) does not exist in phenotype file' % efield
	
	# censor the data
	if delta_field is not None:
		# censor with respect to the interval between efield and delta_field
		assert delta_field in mg, 'Specified delta_field (%s) does not exist in phenotype or genotype file' % delta_field
		print('Censoring %s with respect to the interval between %s and %s' %(efield, efield, delta_field))
		mg['diff'] = mg[delta_field] - mg[efield]
		if np.isfinite(start_time) and np.isnan(end_time):
			print('Start time specified: keeping events with (%s - %s) >= %0.2f' % (delta_field,efield, start_time))
			# final = mg[(mg['diff']>=start_time)|(np.isnan(mg['diff']))] # old behavior keeps nans - for when controls don't have the delta_column?
			final = mg[(mg['diff'] >= start_time)]
		elif np.isnan(start_time) and np.isfinite(end_time):
			print('End time specified: keeping events with (%s - %s) <= %0.2f' % (delta_field, efield, end_time))
			# final = mg[(mg['diff']<=end_time)|(np.isnan(mg['diff']))] # old behavior keeps nans - for when controls don't have the delta_column?
			final = mg[(mg['diff'] <= end_time)]
		else:
			print('Start & End times specified: keeping events with %0.2 <= (%s - %s) <= %0.2f' % (start_time, delta_field, efield, end_time))
			# final = mg[(mg['diff']>=start_time)&(mg['diff']<=end_time)|(np.isnan(mg['diff']))]  # old behavior keeps nans - for when controls don't have the delta_column?
			final = mg[(mg['diff'] >= start_time) & (mg['diff'] <= end_time)]
	else:
		# censor efield ages based on how start/end times are specified
		if np.isfinite(start_time) and np.isnan(end_time):
			print('Start time specified: keeping events with %s >= %0.2f' %(efield,start_time))
			final = mg[mg[efield] >= start_time]
		elif np.isnan(start_time) and np.isfinite(end_time):
			print('End time specified: keeping events with %s <= %0.2f' %(efield,end_time))
			final = mg[mg[efield] <= end_time]
		else:
			print('Start & End times specified: keeping events with %0.2f <= %s <= %0.2f' %(start_time,efield,end_time))
			final = mg[(mg[efield] >= start_time) & (mg[efield] <= end_time)]
	
	# save results
	final_gp = final.drop_duplicates('id')
	print('%d (out of %d) subjects remain after censoring' %(final_gp.shape[0], genotypes.shape[0]))
	if 'genotype' in final_gp:
		num_case = final_gp[final_gp['genotype']==1].shape[0]
		num_ctrl = final_gp[final_gp['genotype']==0].shape[0]
		print('Cases: %d\nControls: %d' %(num_case,num_ctrl))
		
	print("Saving new genotype file to %s" % final_gfile)
	final_gp.to_csv(final_gfile, columns=genotypes.columns, index=False)
	print("Saving new phenotype file to %s" % final_pfile)
	final.to_csv(final_pfile, columns=phenotypes.columns, index=False)



