import pandas as pd
from Bio import Entrez
import pickle
import time
import csv
import argparse
import os.path as osp
import os
from tqdm import tqdm
import numpy as np
import math

umls_cols = ['CUI', 'LAT', 'TS', 'LUI', 'STT', 'SUI', 'ISPREF', 'AUI', 'SAUI', 'SCUI', 'SDUI', 'SAB', 'TTY', 'CODE', 'STR', 'SRL', 'SUPPRESS', 'CVF', 'other']
dev_email = 'cailey.i.kerley@vanderbilt.edu'

def parse_args():
    parser = argparse.ArgumentParser(description="Search PubMed for PheCode-Related Articles")

    parser.add_argument('--umls', required=True, type=str, help='Path to UMLS Metathesaurus (MRCONSO.RRF)')
    parser.add_argument('--outdir', required=True, type=str, help='Path to output directory')
    parser.add_argument('--icd10', required=True, type=str, help ='Path to ICD10-PheCode lookup table')
    parser.add_argument('--icd9', required=True, type=str, help='Path to ICD9-PheCode lookup table')

    args = parser.parse_args()
    return args

def search(query,rs,f):
    Entrez.email = dev_email
    handle = Entrez.esearch(db='pubmed',
                           sort='relevance',
                           retstart=rs,
                           retmax=10000,
                           retmode='xml',
                           field=f,
                           # usehistory='y',
                           term=query)

    results = Entrez.read(handle)
    handle.close()
    return results


def main():
    start = time.time()
    ### Get Args ###
    args = parse_args()

    outdir = args.outdir
    if not osp.exists(outdir):
        os.mkdir(outdir)

    print('Reading ICD10 file...')
    icd10 = pd.read_csv(args.icd10)
    icd10.dropna(subset=['PheCode'],inplace=True)
    print('Reading ICD9 file...')
    icd9 = pd.read_csv(args.icd9)
    icd9.dropna(subset=['PheCode'],inplace=True)
    print('Reading UMLS file...')
    # TODO: check reading with read_csv instead
    umls = pd.read_table(args.umls,sep='|',header=None,names=umls_cols,low_memory=False)
    umls.drop(columns='other',inplace=True)
    umls = umls[umls['LAT'] == 'ENG'] # limit search to english terms

    pubmed_list = {}
    missed = []
    missed_phe = []
    weird = []
    weird_phw = []

    print('Making subset tables... ')
    phecode_list = icd9.drop_duplicates(subset=['PheCode']).sort_values(by=['PheCode']).reset_index(drop=True).copy()
    umls_icd10 = umls[umls['SAB'] == 'ICD10'].copy()
    umls_icd9 = umls[umls['SAB'] == 'ICD9CM'].copy()
    
    with open(osp.join(outdir,'errors.csv'),'w+') as e_file:
        e_file.write('PheCode,SearchType,SearchString,ErrorMsg\n')

        for ix, data in tqdm(phecode_list.iterrows(),total=phecode_list.shape[0]):
            uids = set() # unique identifiers of pubmed articles
            phe = data['PheCode']
            # print(phe)
    
            # Get all icd codes the map to current PheCode
            phe_icd10 = icd10[icd10['PheCode'] == phe]
            phe_icd9 = icd9[icd9['PheCode'] == phe]
            # Get all CUI codes the map those ICD codes
            if phe_icd10.shape[0] > 0:
                cui_icd10 = pd.merge(umls_icd10, phe_icd10, left_on='CODE', right_on='ICD10')
                cui_icd9 = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
                cui = cui_icd10.append(cui_icd9,sort=False)
            else:
                # print('Only ICD9')
                cui = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
            all_cui_str = pd.merge(umls, cui[['CUI']], on='CUI').drop_duplicates(subset=['CUI','STR'])
    
            num_cui = all_cui_str.shape[0]
            for k, group in tqdm(all_cui_str.groupby(np.arange(num_cui)//10),total=len(np.unique(np.arange(num_cui)//10))):
                # build search string from CUI strings
                ss = '('
                first = True
                for ix2, data2 in group.iterrows():
                    if not first: ss = ss + 'OR('
                    else: first = False
                    for term in data2['STR'].split():
                        if len(term)>0:
                            ss=ss+term.strip()+ ' AND '
                    ss = ss[:-5]+ ')'
                    
                try:
                    # search Titles & Abstracts
                    ss_r = search(ss, 0, 'TIAB')
    
                    # parse results
                    gx1 = ss_r['Count']
                    # print(gx1)
                    if int(gx1)>1000000:
                        # dump super big results to weird.pickle because this is weird
                        e_file.write('%s,%s,%s,%s\n' % (str(phe), 'TIAB', ss, 'Weird'))
                        weird_phe.append(phe)
                        weird.append(ss)
                        pickle_wd = open(osp.join(outdir,"weird.pickle"), "wb")
                        pickle.dump(weird_phe, pickle_wd)
                        pickle.dump(weird, pickle_wd)
                        pickle_wd.close()
                        continue
                    else:
                        uids = uids.union(ss_r['IdList'])
                        # only 10,000 results returned at a time, so keep going to get all of the uids
                        for n in range(10000, int(gx1), 10000):
                           uids = uids.union(search(ss, n, 'TIAB')['IdList'])
                           
                except Exception as e:
                    e_file.write('%s,%s,%s,%s\n' %(str(phe),'TIAB',ss,e.args[0]))
                    missed_phe.append(phe)
                    missed.append(ss)
                    pickle_ms = open(osp.join(outdir, "missed.pickle"), "wb")
                    pickle.dump(missed_phe, pickle_ms)
                    pickle.dump(missed, pickle_ms)
                    pickle_ms.close()
                    time.sleep(30)
                    
                try:
                    # repeat search, looking in Keywords field
                    ss_r = search(ss, 0, 'MESH')
    
                    # parse results
                    gx1 = ss_r['Count']
                    # print(gx1)
                    if int(gx1)>1000000:
                        # dump super big results to weird.pickle because this is weird
                        e_file.write('%s,%s,%s,%s\n' % (str(phe), 'MESH', ss, 'Weird'))
                        weird_phe.append(phe)
                        weird.append(ss)
                        pickle_wd = open(osp.join(outdir, "weird.pickle"), "wb")
                        pickle.dump(weird_phe, pickle_wd)
                        pickle.dump(weird, pickle_wd)
                        pickle_wd.close()
                        continue
                    else:
                        uids = uids.union(ss_r['IdList'])
                        # only 10,000 results returned at a time, so keep going to get all of the uids
                        for n in range(10000, int(gx1), 10000):
                            uids = uids.union(search(ss, n, 'MESH')['IdList'])
                            
                except Exception as e:
                    e_file.write('%s,%s,%s,%s\n' %(str(phe),'MESH',ss,e.args[0]))
                    missed_phe.append(phe)
                    missed.append(phe)
                    pickle_ms = open(osp.join(outdir, "missed.pickle"), "wb")
                    pickle.dump(missed_phe, pickle_ms)
                    pickle.dump(missed, pickle_ms)
                    pickle_ms.close()
                    time.sleep(30)
                    
            try:
                # Repeat search process, but search the PheCode string
                ss = '('
                for term in data['Phenotype'].split():
                    # build search string
                    if len(term) > 0:
                        ss = ss + term.strip() + ' '
                ss = ss + ')'
    
                # search Titles & Abstracts
                ss_r = search(ss, 0, 'TIAB')
                # parse results
                gx1 = ss_r['Count']
                # print(gx1)
                if int(gx1) > 1000000:
                    # dump super big results to weird.pickle because this is weird
                    e_file.write('%s,%s,%s,%s\n' % (str(phe), 'TIAB', ss, 'Weird'))
                    weird_phe.append(phe)
                    weird.append(ss)
                    pickle_wd = open(osp.join(outdir, "weird.pickle"), "wb")
                    pickle.dump(weird_phe, pickle_wd)
                    pickle.dump(weird, pickle_wd)
                    pickle_wd.close()
                else:
                    uids = uids.union(ss_r['IdList'])
                    # only 10,000 results returned at a time, so keep going to get all of the uids
                    for n in range(10000, int(gx1), 10000):
                        uids = uids.union(search(ss, n, 'TIAB')['IdList'])
                        
            except Exception as e:
                e_file.write('%s,%s,%s,%s\n' %(str(phe),'TIAB',ss,e.args[0]))
                missed_phe.append(phe)
                missed.append(phe)
                pickle_ms = open(osp.join(outdir, "missed.pickle"), "wb")
                pickle.dump(missed_phe, pickle_ms)
                pickle.dump(missed, pickle_ms)
                pickle_ms.close()
                time.sleep(30)
                
            try:
                # repeat search, looking in Keywords field
                ss_r = search(ss, 0, 'MESH')
                # parse results
                gx1 = ss_r['Count']
                # print(gx1)
                if int(gx1) > 1000000:
                    # dump super big results to weird.pickle because this is weird
                    e_file.write('%s,%s,%s,%s\n' % (str(phe), 'MESH', ss, 'Weird'))
                    weird_phe.append(phe)
                    weird.append(ss)
                    pickle_wd = open(osp.join(outdir, "weird.pickle"), "wb")
                    pickle.dump(weird_phe, pickle_wd)
                    pickle.dump(weird, pickle_wd)
                    pickle_wd.close()
                else:
                    uids = uids.union(ss_r['IdList'])
                    # only 10,000 results returned at a time, so keep going to get all of the uids
                    for n in range(10000, int(gx1), 10000):
                        uids = uids.union(search(ss, n, 'MESH')['IdList'])
                        
            except Exception as e:
                e_file.write('%s,%s,%s,%s\n' %(str(phe),'MESH',ss,e.args[0]))
                missed_phe.append(phe)
                missed.append(phe)
                pickle_ms = open(osp.join(outdir, "missed.pickle"), "wb")
                pickle.dump(missed_phe, pickle_ms)
                pickle.dump(missed, pickle_ms)
                pickle_ms.close()
                time.sleep(30)
    
    
            pubmed_list[phe] = list(uids)
            # print(pubmed_list.items())
            # print(len(uids))
    
            if np.mod(ix, 10) == 0:
                uid_df = pd.DataFrame(pubmed_list.items(),columns=['PheCode','IdsList'])
                uid_df.to_csv(osp.join(outdir,'phecode_pubmed_articles_'+str(ix)+'.csv'),index=False)
                time.sleep(60)
                pubmed_list = {}
    
        uid_df = pd.DataFrame(pubmed_list.items(), columns=['PheCode', 'IdsList'])
        uid_df.to_csv(osp.join(outdir, 'phecode_pubmed_articles_' + str(ix) + '.csv'), index=False)

    interval = time.time() - start
    hour = math.floor(interval / 3600.0)
    minute = math.floor((interval - hour * 3600) / 60)
    second = math.floor(interval - hour * 3600 - minute * 60)

    if hour > 0:
        time_str = '%dh:%dm:%ds' % (hour, minute, second)
    elif minute > 0:
        time_str = '%dm:%ds' % (minute, second)
    else:
        time_str = '%ds' % second

    print('---------------')
    print('PubMed Search Complete\nRuntime: %s' % time_str)
    print('---------------')


if __name__ == '__main__':
    main()
