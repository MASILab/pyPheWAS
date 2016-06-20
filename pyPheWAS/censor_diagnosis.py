def censor_diagnosis(genotype_file,phenotype_file,final_gfile,final_pfile, field ='na',start_time='na',end_time='na'):
        genotypes = pd.read_csv(genotype_file)
        phenotypes = pd.read_csv(phenotype_file)
        mg=pd.merge(phenotypes,genotypes,on='IDNO')
        if start_time=='na' and end_time=='na':
                print "Choose appropriate time period"
        if field=='na':
                if start_time=='na' and end_time!='na': 
                        final = mg[mg['AgeAtICD']>start_time]
                elif start_time!='na' and end_time=='na': 
                        final = mg[mg['AgeAtICD']<end_time]
                else:
                        final = (mg[mg['AgeAtICD']>start_time])&(mg[mg['AgeAtICD']<end_time])

        else:
                mg['diff']=mg[field]-mg['AgeAtICD']
                if start_time=='na' and end_time!='na': 
                        final = mg[mg['diff']>start_time]
                elif start_time!='na' and end_time=='na': 
                        final = mg[mg['diff']<end_time]
                else:
                        final = (mg[mg['diff']>start_time])&(mg[mg['diff']<end_time])
        final[['IDNO','ICD','AgeAtICD']].to_csv(final_gpile)
        final[['IDNO','DXOI',field]].to_csv(final_gfile)

        


