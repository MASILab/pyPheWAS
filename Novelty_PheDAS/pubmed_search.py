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

    args = parser.parse_args()
    return args


# def WriteDictToCSV(csv_file,csv_columns,dict_data):
#     try:
#         with open(csv_file, 'w') as csvfile:
#             writer = csv.DictWriter(csvfile, fieldnames=csv_columns)
#             writer.writeheader()
#             for data in dict_data:
#                 writer.writerow(data)
#     except IOError as (errno, strerror):
#             print("I/O error({0}): {1}".format(errno, strerror))
#     return

def search(query,rs,f):
    Entrez.email = dev_email
    handle = Entrez.esearch(db='pubmed',
                           sort='relevance',
                           retstart=rs,
                           retmax=10000,
                           retmode='xml',
                           field=f,
                           usehistory='y',
                           term=query)

    results = Entrez.read(handle)
    handle.close()
    return results

# def searchf(query,rs,f,w):
#     Entrez.email = dev_email
#     handle = Entrez.esearch(db='pubmed',
#                            sort='relevance',
#                            retstart=rs,
#                            retmax=10000,
#                            retmode='xml',
#                            field=f,
#                            WebEnv=w,
#                            usehistory='y',
#                            term=query)
#
#     results = Entrez.read(handle)
#     handle.close()
#     return results


def main():
    ### Get Args ###
    args = parse_args()

    print('Reading ICD10 file...')
    icd10 = pd.read_csv(args.icd10)
    icd10.dropna(subset=['PheCode'],inplace=True)
    print('Reading ICD9 file...')
    icd9 = pd.read_csv(args.icd9)
    icd9.dropna(subset=['PheCode'],inplace=True)
    print('Reading UMLS file...')
    umls = pd.read_table(args.umls,sep='|',header=None,names=umls_cols,low_memory=False)
    umls.drop(columns='other',inplace=True)
    outdir = args.outdir
    if ~ osp.exists(outdir):
        os.mkdir(outdir)

    pubmed_list = {}
    missed = []
    weird = []

    print('Making subset tables... ')
    phecode_list = icd9.drop_duplicates(subset=['PheCode']).reset_index(drop=True).sort_values(by=['PheCode']).copy()
    umls_icd10 = umls[umls['SAB'] == 'ICD10'].copy()
    umls_icd9 = umls[umls['SAB'] == 'ICD9CM'].copy()

    for ix, data in tqdm(phecode_list.iterrows()):
        try:
            uids = set() # unique identifiers of pubmed articles
            phe = data['PheCode']

            # Get all icd codes the map to current PheCode
            phe_icd10 = icd10[icd10['PheCode'] == phe]
            phe_icd9 = icd9[icd9['PheCode'] == phe]
            # Get all CUI codes the map those ICD codes
            if phe_icd10.shape[0] > 0:
                cui_icd10 = pd.merge(umls_icd10, phe_icd10, left_on='CODE', right_on='ICD10')
                cui_icd9 = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
                cui = cui_icd10.append(cui_icd9)
            else:
                cui = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
            all_cui_str = pd.merge(umls, cui[['CUI']], on='CUI')

            for ix, data in tqdm(all_cui_str.iterrows()):
                # build search string from CUI str
                ss='('
                for term in data['STR'].split():
                    if len(term)>0:
                        ss=ss+term.strip()+ '+'
                ss=ss[:-1]+ ')'

                # search Titles & Abstracts
                ss_r = search(ss, 0, 'TIAB')

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


                # repeat search, looking in Keywords field
                ss_r = search(ss, 0, 'MESH')

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
                        uids = uids.union(search(ss, n, 'MESH')['IdList'])


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

            # repeat search, looking in Keywords field
            ss_r = search(ss, 0, 'MESH')
            # parse results
            gx1 = ss_r['Count']
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


            pubmed_list[phe] = list(uids)

            if np.mod(ix, 10) == 0:
               pd.DataFrame(pubmed_list.items(),columns=['PheCode','IdsList']).to_csv('phecode_pubmed_articles_'+str(ix)+'.csv',index=False)
               time.sleep(60)
               pubmed_list = {}

        except:
           missed.append(phe)
           pickle_ms = open(osp.join(outdir,"missed.pickle"), "wb")
           pickle.dump(missed, pickle_ms)
           pickle_ms.close()
           time.sleep(60)

        break # testing

    pd.DataFrame(pubmed_list.items(), columns=['PheCode', 'IdsList']).to_csv('phecode_pubmed_articles_' + str(ix) + '.csv', index=False)
                                                                                

if __name__ == '__main__':
    main()