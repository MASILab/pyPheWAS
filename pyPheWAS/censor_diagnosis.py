
def censor_diagnosis(path,genotype_file,phenotype_file,final_pfile, final_gfile, field ='na',type='ICD9',start_time=float('nan'),end_time=float('nan')):
        import pandas as pd
        import numpy as np
        genotypes = pd.read_csv(path+genotype_file)
        phenotypes = pd.read_csv(path+phenotype_file)
        mg=pd.merge(phenotypes,genotypes,on='id')
        if np.isnan(start_time) and np.isnan(end_time):
                print("Choose appropriate time period")
        if field=='na':
                if np.isfinite(start_time) and np.isnan(end_time):
                        final = mg[mg['AgeAt'+type]>=start_time]
                elif np.isnan(start_time) and np.isfinite(end_time):
                        final = mg[mg['AgeAt'+type]<=end_time]
                else:
                        final = mg[(mg['AgeAt'+type]>=start_time)&(mg['AgeAtICD']<=end_time)]

        else:
                mg['diff']=mg[field]-mg['AgeAt'+type]
                if np.isfinite(start_time) and np.isnan(end_time):
                        final = mg[(mg['diff']>=start_time)|(np.isnan(mg['diff']))]
                elif np.isnan(start_time) and np.isfinite(end_time):
                        final = mg[(mg['diff']<=end_time)|(np.isnan(mg['diff']))]
                else:
                        final = mg[(mg['diff']>=start_time)&(mg['diff']<=end_time)|(np.isnan(mg['diff']))]
        final['MaxAgeBeforeDx'] = final.groupby('id')['AgeAt'+type].transform('max')
        final.dropna(subset=['MaxAgeBeforeDx'],inplace=True)
        final[['id',type.lower(),'AgeAt'+type]].to_csv(path+final_pfile)
        final[['id', 'MaxAgeAtVisit','MaxAgeBeforeDx', 'AgeAtDx', 'genotype','sex']].drop_duplicates().to_csv(path+final_gfile)

        


