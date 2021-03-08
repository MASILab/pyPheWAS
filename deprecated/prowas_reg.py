from pyProWAS import *
import os
import numpy as np
import pandas as pd

reg_type = 0
str_reg_type = "log"
path = "/nfs/share5/clineci/DownSyndrome/experiments/prowas_test/"
filename = "cpts_age.csv"
groupfile = "group.csv"
phewas_cov = ''
outfile = 'feature_matrix.csv'
covariates = 'SEX'
str_thresh_type = "fdr"
thresh_type = 1

"""
# gen_ftype = reg_type
phenotypes = get_input(path, filename,reg_type)
genotypes = get_group_file(path, groupfile)
fm = generate_feature_matrix(genotypes, phenotypes, reg_type)

print("Saving feature matrices to %s" % (path + outfile))

np.savetxt(path + 'agg_measures_' + outfile, fm[0],delimiter=',')
print("...")
np.savetxt(path + 'icd_age_' + outfile, fm[1],delimiter=',')
print("...")
np.savetxt(path + 'phewas_cov_' + outfile, fm[2],delimiter=',')

regressions = run_phewas(fm, genotypes, covariates,reg_type)

print("Saving regression data to %s" % (path + 'regressions.csv'))
header = ','.join(['str_reg_type', str_reg_type, 'group', groupfile]) + '\n'
f = open(os.sep.join([path, 'regressions.csv']), 'w')
f.write(header)
regressions.to_csv(f,index=False)
f.close()

"""

regressions = pd.read_csv(path + 'regressions.csv',dtype={'PheWAS Code':str},skiprows=1)

print("creating plots")

# Check if an imbalance will be used

imbalances = get_imbalances(regressions)

y = regressions['"-log(p)"']
pvalues = regressions['p-val'].values

# Get the threshold type
if thresh_type == 0:
    thresh = get_bon_thresh(pvalues, 0.05)
elif thresh_type == 1:
    thresh = get_fdr_thresh(pvalues, 0.05)

thresh = 0.5
print('%s threshold: %0.5f'%(str_thresh_type,thresh))

try:
    regressions[['lowlim', 'uplim']] = regressions['Conf-interval beta'].str.split(',', expand=True)
    regressions['uplim'] = regressions.uplim.str.replace(']', '')
    regressions['lowlim'] = regressions.lowlim.str.replace('[', '')
    regressions = regressions.astype(dtype={'uplim':float,'lowlim':float})
    yb = regressions[['beta', 'lowlim', 'uplim']].values
    yb = yb.astype(float)
except Exception as e:
    print('Error reading regression file:')
    print(e)
    sys.exit()

save = path + 'plot.png'
file_name, file_format = os.path.splitext(save)
saveb = file_name + '_beta' + file_format
file_format = file_format[1:] # remove '.' from from first index
print("Saving plot to %s" % (save))

plot_manhattan(regressions, -math.log10(thresh), save=save, save_format=file_format)
plot_odds_ratio(regressions, -math.log10(thresh), save=saveb, save_format=file_format)