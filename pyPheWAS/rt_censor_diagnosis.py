def censor_diagnosis(genotype_file,phenotype_file,final_pfile, final_gfile,field,start_time=float('nan'),end_time=float('nan')):
        import pandas as pd
        import numpy as np
        genotypes = pd.read_csv(genotype_file)
        phenotypes = pd.read_csv(phenotype_file)
        mg=pd.merge(phenotypes,genotypes,on='id')
        assert field in mg, 'Specified field (%s) does not exist in phenotype or genotype field' %field

        if np.isfinite(start_time) and np.isnan(end_time):
                final = mg[mg[field]>=start_time]
        elif np.isnan(start_time) and np.isfinite(end_time):
                final = mg[mg[field]<=end_time]
        else:
                final = mg[(mg[field]>=start_time)&(mg[field]<=end_time)]

        # for censoring with respect to the interval between AgeAtICD and a specified field (i.e. AgeAtDX or age at a scan time)
        # else:
        #         mg['diff']=mg[field]-mg['AgeAtICD']
        #         if np.isfinite(start_time) and np.isnan(end_time):
        #                 final = mg[(mg['diff']>=start_time)|(np.isnan(mg['diff']))]
        #         elif np.isnan(start_time) and np.isfinite(end_time):
        #                 final = mg[(mg['diff']<=end_time)|(np.isnan(mg['diff']))]
        #         else:
        #                 final = mg[(mg['diff']>=start_time)&(mg['diff']<=end_time)|(np.isnan(mg['diff']))]

        final.to_csv(final_pfile,columns=phenotypes.columns,index=False)
        print("New phenotype file in %s" % final_pfile)
        final_gp = final.drop_duplicates('id')
        final_gp.to_csv(final_gfile,columns=genotypes.columns,index=False)
        print("New group file in %s" % final_gfile)
        


