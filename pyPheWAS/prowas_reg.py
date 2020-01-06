from pyProWAS import *
import os
import numpy as np

reg_type = 0
str_reg_type = "log"
path = "/nfs/share5/clineci/DownSyndrome/experiments/prowas_test/"
filename = "cpts_age.csv"
groupfile = "group.csv"
phewas_cov = ''
outfile = 'feature_matrix.csv'
covariates = 'SEX'

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

results = run_phewas(fm, genotypes, covariates,reg_type)

print("Saving regression data to %s" % (path + 'regressions.csv'))
header = ','.join(['str_reg_type', str_reg_type, 'group', groupfile]) + '\n'
f = open(os.sep.join([path, 'regressions.csv']), 'w')
f.write(header)
results.to_csv(f,index=False)
f.close()


