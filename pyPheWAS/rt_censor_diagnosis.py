
def censor_diagnosis(genotype_file,phenotype_file,final_pfile, field ='na',start_time=float('nan'),end_time=float('nan')):
        import pandas as pd
        import numpy as np
        genotypes = pd.read_csv(genotype_file)
        phenotypes = pd.read_csv(phenotype_file)
        mg=pd.merge(phenotypes,genotypes,on='id')
        if np.isnan(start_time) and np.isnan(end_time):
                print("Choose appropriate time period")
        if field=='na':
                if np.isfinite(start_time) and np.isnan(end_time):
                        final = mg[mg['AgeAtICD']>start_time]
                elif np.isnan(start_time) and np.isfinite(end_time):
                        final = mg[mg['AgeAtICD']<end_time]
                else:
                        final = mg[(mg['AgeAtICD']>start_time)&(mg['AgeAtICD']<end_time)]

        else:
                mg['diff']=mg[field]-mg['AgeAtICD']
                if np.isfinite(start_time) and np.isnan(end_time):
                        final = mg[(mg['diff']>start_time)|(np.isnan(mg['diff']))]
                elif np.isnan(start_time) and np.isfinite(end_time):
                        final = mg[(mg['diff']<end_time)|(np.isnan(mg['diff']))]
                else:
                        final = mg[(mg['diff']>start_time)&(mg['diff']<end_time)|(np.isnan(mg['diff']))]
        final[['id','icd9','AgeAtICD']].to_csv(final_pfile)

        


