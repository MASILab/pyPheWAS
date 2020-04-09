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

umls_cols = ['CUI', 'LAT', 'TS', 'LUI', 'STT', 'SUI', 'ISPREF', 'AUI', 'SAUI', 'SCUI', 'SDUI', 'SAB', 'TTY', 'CODE', 'STR', 'SRL', 'SUPPRESS', 'CVF', 'other']
dev_email = 'cailey.i.kerley@vanderbilt.edu'

def parse_args():
    parser = argparse.ArgumentParser(description="Search PubMed for PheCode-Related Articles")

    parser.add_argument('--umls', required=True, type=str, help='Path to UMLS Metathesaurus (MRCONSO.RRF)')
    parser.add_argument('--outdir', required=True, type=str, help='Path to output directory')
    parser.add_argument('--icd10', required=True, type=str, help ='Path to ICD10-PheCode lookup table')
    parser.add_argument('--icd9', required=True, type=str, help='Path to ICD9-PheCode lookup table')
    parser.add_argument('--phe', required=True, type=float, help='PheCode to search')

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
                           #usehistory='y',
                           term=query)

    results = Entrez.read(handle)
    handle.close()
    return results


def main():
    ### Get Args ###
    args = parse_args()

    outdir = args.outdir
    if not osp.exists(outdir):
        os.mkdir(outdir)
    phe = args.phe

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

    print('Making subset tables... ')
    phecode_list = icd9.drop_duplicates(subset=['PheCode']).reset_index(drop=True).sort_values(by=['PheCode']).copy()
    umls_icd10 = umls[umls['SAB'] == 'ICD10'].copy()
    umls_icd9 = umls[umls['SAB'] == 'ICD9CM'].copy()

    print('*****')
    print(phe)
    print('*****')
    tiab_list = {}
    mesh_list = {}
    pubmed_list = {}
    uids_all = set()
    missed = []
    weird = []
    

    # Get all icd codes the map to current PheCode
    phe_icd10 = icd10[icd10['PheCode'] == phe]
    phe_icd9 = icd9[icd9['PheCode'] == phe]
    # Get all CUI codes the map those ICD codes
    if phe_icd10.shape[0] > 0:
        cui_icd10 = pd.merge(umls_icd10, phe_icd10, left_on='CODE', right_on='ICD10')
        cui_icd9 = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
        cui = cui_icd10.append(cui_icd9,sort=False)
    else:
        cui = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
    all_cui_str = pd.merge(umls, cui[['CUI']], on='CUI').drop_duplicates(subset=['CUI','STR'])
    cui.drop_duplicates(subset=['CUI'],inplace=True)

    num_str = all_cui_str.shape[0]
    print('%d CUIs' % cui.shape[0])
    print('%d CUI Strings' % num_str)

    for k, group in tqdm(all_cui_str.groupby(np.arange(num_str) // 10), total=len(np.unique(np.arange(num_str) // 10))):
        # build search string from CUI strings
        ss = '('
        first = True
        for ix2, data2 in group.iterrows():
            if not first:
                ss = ss + 'OR('
            else:
                first = False
            for term in data2['STR'].split():
                if len(term) > 0:
                    ss = ss + term.strip() + ' AND '
            ss = ss[:-5]+ ')'
        
        try:
            # search Titles & Abstracts
            ss_r = search(ss, 0, 'TIAB')
            uids = set() # unique identifiers of pubmed articles

            # parse results
            gx1 = ss_r['Count']
            if int(gx1)>1000000:
                # dump super big results to weird.pickle because this is weird
                weird.append(ss)
                pickle_wd = open(osp.join(outdir,"weird.pickle"), "wb")
                pickle.dump(weird, pickle_wd)
                pickle_wd.close()
                continue
            else:
                uids = uids.union(ss_r['IdList'])
                # only 10,000 results returned at a time, so keep going to get all of the uids
                for n in range(10000, int(gx1), 10000):
                    uids = uids.union(search(ss, n, 'TIAB')['IdList'])
            
            tiab_list[ss] = list(uids)
            uids_all = uids_all.union(uids)
                    
        except Exception as e:
            print('Error searching TIAB with: ' + str(ss))
            print(e.args[0])
            missed.append(phe)
            pickle_ms = open(osp.join(outdir,"missed_TIAB.pickle"), "wb")
            pickle.dump(missed, pickle_ms)
            pickle_ms.close()
            time.sleep(60)

        try:
            # repeat search, looking in Keywords field
            ss_r = search(ss, 0, 'MESH')
            uids = set() # unique identifiers of pubmed articles

            # parse results
            gx1 = ss_r['Count']
            # print(gx1)
            if int(gx1)>1000000:
                # dump super big results to weird.pickle because this is weird
                weird.append(ss)
                pickle_wd = open(osp.join(outdir,"weird.pickle"), "wb")
                pickle.dump(weird, pickle_wd)
                pickle_wd.close()
                continue
            else:
                uids = uids.union(ss_r['IdList'])
                # only 10,000 results returned at a time, so keep going to get all of the uids
                for n in range(10000, int(gx1), 10000):
                    uids = uids.union(search(ss, n, 'MESH')['IdList'])
                    
            mesh_list[ss] = list(uids)
            uids_all = uids_all.union(uids)
                    
        except Exception as e:
            print('Error searching MESH with: ' + str(ss))
            print(e.args[0])
            missed.append(phe)
            pickle_ms = open(osp.join(outdir,"missed_MESH.pickle"), "wb")
            pickle.dump(missed, pickle_ms)
            pickle_ms.close()
            time.sleep(60)
                    
    try:       
        # Repeat search process, but search the PheCode string
        ss = '('
        phenotype = icd9.loc[icd9['PheCode'] == phe, 'Phenotype'].values[0]
        print(phenotype)
        for term in phenotype.split():
            # build search string
            if len(term) > 0:
                ss = ss + term.strip() + ' '
        ss = ss + ')'

        # search Titles & Abstracts
        ss_r = search(ss, 0, 'TIAB')
        uids = set() # unique identifiers of pubmed articles
        
        # parse results
        gx1 = ss_r['Count']
        # print(gx1)
        if int(gx1) > 1000000:
            # dump super big results to weird.pickle because this is weird
            weird.append(ss)
            pickle_wd = open(osp.join(outdir,"weird.pickle"), "wb")
            pickle.dump(weird, pickle_wd)
            pickle_wd.close()
        else:
            uids = uids.union(ss_r['IdList'])
            # only 10,000 results returned at a time, so keep going to get all of the uids
            for n in range(10000, int(gx1), 10000):
                uids = uids.union(search(ss, n, 'TIAB')['IdList'])
        tiab_list[ss] = list(uids)
        uids_all = uids_all.union(uids)
                
    except Exception as e:
        print('Error searching TIAB with: ' + str(ss))
        print(e.args[0])
        missed.append(phe)
        pickle_ms = open(osp.join(outdir,"missed_TIAB.pickle"), "wb")
        pickle.dump(missed, pickle_ms)
        pickle_ms.close()
        # time.sleep(60)
        
    try:
        # repeat search, looking in Keywords field
        ss_r = search(ss, 0, 'MESH')
        uids = set()  # unique identifiers of pubmed articles
        # parse results
        gx1 = ss_r['Count']
        # print(gx1)
        if int(gx1) > 1000000:
            # dump super big results to weird.pickle because this is weird
            weird.append(ss)
            pickle_wd = open(osp.join(outdir,"weird.pickle"), "wb")
            pickle.dump(weird, pickle_wd)
            pickle_wd.close()
        else:
            uids = uids.union(ss_r['IdList'])
            # only 10,000 results returned at a time, so keep going to get all of the uids
            for n in range(10000, int(gx1), 10000):
                uids = uids.union(search(ss, n, 'MESH')['IdList'])
                
        mesh_list[ss] = list(uids)
        uids_all = uids_all.union(uids)
                
    except Exception as e:
        print('Error searching MESH with: ' + str(ss))
        print(e.args[0])
        missed.append(phe)
        pickle_ms = open(osp.join(outdir,"missed_MESH.pickle"), "wb")
        pickle.dump(missed, pickle_ms)
        pickle_ms.close()
        # time.sleep(60)

    pubmed_list[str(phe)] = list(uids_all)
    print(len(uids_all))
        
    # save results to Excel workbook
    with pd.ExcelWriter(osp.join(outdir,'pubmed_search_' + str(phe) + '.xlsx')) as ex_f:
        cui.to_excel(ex_f,sheet_name='CUIs',index=False)
        all_cui_str.to_excel(ex_f,sheet_name='CUI Strings',index=False)
        
    tiab_df = pd.DataFrame(tiab_list.items(),columns=['Search String', 'IdsList'])
    tiab_df.to_csv(osp.join(outdir,'pubmed_search_' + str(phe) + '_TIAB_AND.csv'),index=False)
    mesh_df = pd.DataFrame(mesh_list.items(),columns=['Search String', 'IdsList'])
    mesh_df.to_csv(osp.join(outdir,'pubmed_search_' + str(phe) + '_MESH_AND.csv'),index=False)
    all_df = pd.DataFrame(pubmed_list.items(), columns=['PheCode', 'IdsList'])
    all_df.to_csv(osp.join(outdir,'pubmed_search_' + str(phe) + '_AND.csv'), index=False)
                                                                                

if __name__ == '__main__':
    main()
