from pyProWAS import *
gen_ftype = reg_type
phenotypes = get_input(path, filename)
genotypes = get_group_file(path, groupfile)
fm = generate_feature_matrix(genotypes, phenotypes, phewas_cov)
results = run_phewas(fm, genotypes, covariates, phewas_cov=phewas_cov)
