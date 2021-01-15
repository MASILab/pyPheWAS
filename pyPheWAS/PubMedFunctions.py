from pyPheWAS.pyPhewasCorev2 import icd9_codes, icd10_codes
import pandas as pd
from pathlib import Path
import time
import numpy as np
from Bio import Entrez
import pickle
from tqdm import tqdm

dev_email = 'cailey.i.kerley@vanderbilt.edu'

def load_umls(umls_path):
    umls_cols = ['CUI', 'LAT', 'TS', 'LUI', 'STT', 'SUI', 'ISPREF', 'AUI', 'SAUI', 'SCUI', 'SDUI', 'SAB', 'TTY', 'CODE',
                 'STR', 'SRL', 'SUPPRESS', 'CVF', 'other']
    umls = pd.read_table(umls_path, sep='|', header=None, names=umls_cols, low_memory=False)
    umls.drop(columns='other', inplace=True)
    umls = umls[umls['LAT'] == 'ENG']  # limit search to english terms
    return umls


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

def run_PheWAS_PubMed_Query(umls, outdir):
    # make sets to hold publication IDs
    pubmed_list = {}
    missed = []
    missed_phe = []
    weird = []
    weird_phe = []

    print('Making ICD subset tables... ')
    phecode_list = icd9_codes.drop_duplicates(subset=['PheCode']).sort_values(by=['PheCode']).reset_index(drop=True).copy()
    umls_icd10 = umls[umls['SAB'] == 'ICD10'].copy()
    umls_icd9 = umls[umls['SAB'] == 'ICD9CM'].copy()

    with open(outdir / 'errors.csv', 'w+') as e_file:
        e_file.write('PheCode,SearchType,SearchString,ErrorMsg\n')

        for ix, data in tqdm(phecode_list.iterrows(), total=phecode_list.shape[0]):
            uids = set()  # unique identifiers of pubmed articles
            phe = data['PheCode']

            # Get all icd codes that map to current PheCode
            phe_icd10 = icd10_codes[icd10_codes['PheCode'] == phe]
            phe_icd9 = icd9_codes[icd9_codes['PheCode'] == phe]
            # Get all CUI codes the map those ICD codes
            if phe_icd10.shape[0] > 0:
                cui_icd10 = pd.merge(umls_icd10, phe_icd10, left_on='CODE', right_on='ICD10')
                cui_icd9 = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
                cui = cui_icd10.append(cui_icd9, sort=False)
            else:
                cui = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD9')
            all_cui_str = pd.merge(umls, cui[['CUI']], on='CUI').drop_duplicates(subset=['CUI', 'STR'])

            num_cui = all_cui_str.shape[0]
            # search CUI strings 10 at a time (STR0 or STR1 or .. or STR9)
            for k, group in tqdm(all_cui_str.groupby(np.arange(num_cui) // 10),
                                 total=len(np.unique(np.arange(num_cui) // 10)),
                                 desc='Searching for PheCode %s' %phe):
                # build search string from CUI strings
                ss = '('
                first = True
                for ix2, data2 in group.iterrows():
                    if not first:
                        ss = ss + 'OR('
                    else:
                        first = False
                    terms = data2['STR'].split()
                    ss = ss + ' AND '.join(terms) + ')'
                    # for term in data2['STR'].split():
                    #     if len(term) > 0:
                    #         ss = ss + term.strip() + ' AND '
                    ss = ss[:-5] + ')'

                try:
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
                        continue
                    else:
                        uids = uids.union(ss_r['IdList'])
                        # only 10,000 results returned at a time, so keep going to get all of the uids
                        for n in range(10000, int(gx1), 10000):
                            uids = uids.union(search(ss, n, 'TIAB')['IdList'])

                except Exception as e:
                    e_file.write('%s,%s,%s,%s\n' % (str(phe), 'TIAB', ss, e.args[0]))
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
                    if int(gx1) > 1000000:
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
                    e_file.write('%s,%s,%s,%s\n' % (str(phe), 'MESH', ss, e.args[0]))
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
                e_file.write('%s,%s,%s,%s\n' % (str(phe), 'TIAB', ss, e.args[0]))
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
                e_file.write('%s,%s,%s,%s\n' % (str(phe), 'MESH', ss, e.args[0]))
                missed_phe.append(phe)
                missed.append(phe)
                pickle_ms = open(osp.join(outdir, "missed.pickle"), "wb")
                pickle.dump(missed_phe, pickle_ms)
                pickle.dump(missed, pickle_ms)
                pickle_ms.close()
                time.sleep(30)

            pubmed_list[phe] = list(uids)

            if np.mod(ix, 10) == 0:
                uid_df = pd.DataFrame(pubmed_list.items(), columns=['PheCode', 'IdsList'])
                uid_df.to_csv(osp.join(outdir, 'phecode_pubmed_articles_' + str(ix) + '.csv'), index=False)
                time.sleep(60)
                pubmed_list = {}

        uid_df = pd.DataFrame(pubmed_list.items(), columns=['PheCode', 'IdsList'])
        uid_df.to_csv(osp.join(outdir, 'phecode_pubmed_articles_' + str(ix) + '.csv'), index=False)