from pyPheWAS.pyPhewasCorev2 import icd9_codes, icd10_codes
import pandas as pd
import time
import numpy as np
from Bio import Entrez
from tqdm import tqdm

dev_email = 'cailey.i.kerley@vanderbilt.edu'
umls_cols = ['CUI', 'LAT', 'TS', 'LUI', 'STT', 'SUI', 'ISPREF', 'AUI', 'SAUI', 'SCUI', 'SDUI', 'SAB', 'TTY', 'CODE',
                 'STR', 'SRL', 'SUPPRESS', 'CVF', 'other']
usable_cols = ['CUI','LAT','SAB','CODE','STR']

def load_umls(umls_path):
    umls = pd.read_table(umls_path, sep='|', header=None, names=umls_cols, skiprows=1, low_memory=False, usecols=usable_cols)
    umls = umls[umls['LAT']=='ENG']
    return umls

def load_search_terms(term_file):
    df = pd.read_csv(term_file, header=None, names=['search_terms'])
    return df


def entrez_query(query,rs,f):
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


def pubmed_search(search_str, search_type, error_fpath, phecode=None):
    """ Run PubMed search and parse results """
    try:
        ss_results = entrez_query(search_str, 0, search_type)
        gx1 = ss_results['Count']  # total number of results found by search
        if int(gx1) > 1000000:  # large magnitude of results is weird - raise exception
            msg = 'Weird - %s results' % gx1
            raise Exception(msg)
        else:
            result_ids = set(ss_results['IdList'])
            # only 10,000 results are returned at a time, so keep going to get all of the uids
            for n in range(10000, int(gx1), 10000):
                result_ids = result_ids.union(entrez_query(search_str, n, search_type)['IdList'])
    except Exception as e:
        with open(error_fpath, 'a+') as error_log:
            if phecode is not None:
                error_log.write('%s,%s,%s,%s\n' % (str(phecode), search_str, search_type, str(e)))
            else:
                error_log.write('%s,%s,%s\n' % (search_str, search_type, str(e)))
        time.sleep(10)
        return None
    return result_ids


def run_Custom_PubMed_Query(search_terms, outdir):
    """
    Search PubMed database for literature that mentions custom search terms.

    :param search_terms: dataframe of search terms
    :param outdir: path to directory in which to save results (custom_pubmed_search.csv) & error log
    :return: None
    """

    # init error log file
    error_fpath = outdir / 'errors.csv'
    with open(error_fpath, 'w+') as error_log:
        error_log.write('SearchString,SearchType,ErrorMsg\n')

    # search terms 10 at a time (STR0 or STR1 or .. or STR9)
    uids = set() # set to hold unique identifiers of pubmed articles
    num_terms = search_terms.shape[0]
    for k, group in tqdm(search_terms.groupby(np.arange(num_terms) // 10),
                         total=len(np.unique(np.arange(num_terms) // 10)),
                         desc='Running custom PubMed Search'):
        # build search string
        ss = '('
        first = True
        for ix2, data2 in group.iterrows():
            if not first:
                ss = ss + 'OR('
            else:
                first = False
            terms = data2['search_terms'].split()
            ss = ss + ' + '.join(terms) + ')'

        # run 2 different search types
        res = pubmed_search(ss, 'TIAB', error_fpath)  # search Titles & Abstracts for the constructed search string (ss)
        if res is not None:
            uids = uids.union(res)
        res = pubmed_search(ss, 'MESH', error_fpath)  # search Keywords for the constructed search string (ss)
        if res is not None:
            uids = uids.union(res)

    print('Custom Search Complete - Found %d Articles' % len(uids))

    out_f = outdir / 'custom_pubmed_search.csv'
    print('Saving results to %s' % str(out_f))
    uid_df = pd.DataFrame(columns=['Search', 'IdsList'])
    uid_df.loc[0,'Search'] = 'custom_search'
    uid_df.loc[0,'IdsList'] = list(uids)
    uid_df.to_csv(out_f, index=False)

    return


def run_PheWAS_PubMed_Query(umls, outdir):
    """
    Search PubMed database for literature concerning all PheCodes.

    :param umls: UMLs Metathesaurus dataframe
    :param outdir: pathlib Path to results directory
    :return: None
    """

    # init error log file
    error_fpath = outdir / 'errors.csv'
    with open(error_fpath, 'w+') as error_log:
        error_log.write('PheCode,SearchString,SearchType,ErrorMsg\n')

    print('Making ICD subset tables... ')
    phecode_list = icd9_codes.drop_duplicates(subset=['PheCode']).sort_values(by=['PheCode']).reset_index(drop=True).copy()
    umls_icd10 = umls[umls['SAB'] == 'ICD10'].copy()
    umls_icd9 = umls[umls['SAB'] == 'ICD9CM'].copy()

    # run pubmed search for one PheCode at a time
    pubmed_list = {} # dictionary to hold pubmed results for each phecode
    pubmed_counts = [] # list to hold number of pubmed search results for each phecode
    for ix, data in tqdm(phecode_list.iterrows(), total=phecode_list.shape[0]):
        uids = set()  # set to hold unique identifiers of pubmed articles for each PheCode
        phe = data['PheCode']

        # Get all icd codes that map to current PheCode
        phe_icd10 = icd10_codes[icd10_codes['PheCode'] == phe]
        phe_icd9 = icd9_codes[icd9_codes['PheCode'] == phe]
        # Get all CUI codes the map those ICD codes
        if phe_icd10.shape[0] > 0:
            cui_icd10 = pd.merge(umls_icd10, phe_icd10, left_on='CODE', right_on='ICD_CODE')
            cui_icd9 = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD_CODE')
            cui = cui_icd10.append(cui_icd9, sort=False)
        else:
            cui = pd.merge(umls_icd9, phe_icd9, left_on='CODE', right_on='ICD_CODE')
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
                ss = ss + ' + '.join(terms) + ')'

            # run 2 different search types
            res = pubmed_search(ss, 'TIAB', error_fpath, phecode=phe) # search Titles & Abstracts for the constructed search string (ss)
            if res is not None:
                uids = uids.union(res)
            res = pubmed_search(ss, 'MESH', error_fpath, phecode=phe) # search Keywords for the constructed search string (ss)
            if res is not None:
                uids = uids.union(res)

        # repeat, but this time search just for the PheCode string itself
        ss_phe = data['Phenotype']

        res = pubmed_search(ss_phe, 'TIAB', error_fpath, phecode=phe)  # search Titles & Abstracts for the constructed search string (ss)
        if res is not None:
            uids = uids.union(res)
        res = pubmed_search(ss_phe, 'MESH', error_fpath, phecode=phe)  # search Keywords for the constructed search string (ss)
        if res is not None:
            uids = uids.union(res)

        # done searching - add uids to pubmed result list
        pubmed_list[phe] = list(uids)
        pubmed_counts.append(len(uids))

        # write results to a file every 10 PheCodes
        if np.mod(ix+1, 10) == 0:
            uid_df = pd.DataFrame(pubmed_list.items(), columns=['PheWAS Code', 'IdsList'])
            uid_df['Publication_Count'] = pubmed_counts
            uid_df.to_csv(outdir / ('phecode_pubmed_articles_' + str(ix) + '.csv'), index=False)
            pubmed_list = {}
            pubmed_counts = []
            time.sleep(60) # pause before resuming so we don't kill PubMed
    # save the last dataframe
    uid_df = pd.DataFrame(pubmed_list.items(), columns=['PheCode', 'IdsList'])
    uid_df.to_csv(outdir /('phecode_pubmed_articles_' + str(ix) + '.csv'), index=False)
    return