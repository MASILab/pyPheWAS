import pandas as pd
import argparse
import os
import numpy as np
import os.path as osp
from ast import literal_eval
from tqdm import tqdm


def parse_args():
    parser = argparse.ArgumentParser(description="Add PubMED results to pyPheWAS stat file")

    parser.add_argument('--statfile', required=True, type=str, help='Name of the stat file (e.g. regressions.csv)')
    parser.add_argument('--dx_pm', required=True, type=str, help='Name of the Dx PubMED file (e.g. dx_PubMED_results.csv)')
    parser.add_argument('--pm_dir', required=True, type=str, help ='Path to PheCode PubMED directory')
    parser.add_argument('--path', required=False, default='.', type=str, help='Path to all input files and destination of output files')
    parser.add_argument('--outfile', required=False, default=None, type=str, help='Name of the updated regression output file (default: same as statfile)')

    args = parser.parse_args()
    return args


def main():
    ### Get Args ###
    args = parse_args()

    if args.outfile is None:
        outfile = osp.join(args.path, args.statfile)
    else:
        outfile = osp.join(args.path, args.outfile)

    reg_f = open(osp.join(args.path, args.statfile))
    reg_hdr = reg_f.readline()
    # TODO: Add str back in - fix saving of reg file
    reg = pd.read_csv(reg_f) # , dtype={"PheWAS Code":str})
    reg_f.close()

    reg.rename(columns={'Conf-interval beta':'beta.ci','std_error':'beta.se'},inplace=True)

    dx_pubmed = pd.read_csv(osp.join(args.path, args.dx_pm))
    pubmed_dir = args.pm_dir

    ### Set-up Dx list of UIDs ###
    dx_set = set(literal_eval(dx_pubmed.loc[0, "IdsList"]))
    dx_count = len(dx_set)
    reg["AD.PM.count"] = dx_count

    ### Get PheCode PubMED counts & joint PubMED counts ###
    reg["phe.PM.count"] = np.nan
    reg["joint.PM.count"] = np.nan
    reg["P.PM.phe.given.AD"] = np.nan
    reg["P.PM.AD.given.phe"] = np.nan

    tmp = open("tmp.txt","w+")

    pubmed_file_list = os.listdir(pubmed_dir)
    for j in tqdm(range(0,1856)):
        fname = "phewas_counts_%d.csv" % j
        if fname in pubmed_file_list:
            try:
                tmp.write("%s\n" % fname)
                phecode_pubmed = pd.read_csv(osp.join(pubmed_dir, fname)) # , dtype={"phewas_code": str})
            except Exception as e:
                ef = open("novelty_pubmed_errors.csv", 'a+')
                ef.write("error opening%s,%s" % (fname, e.args[0]))
                ef.close()

            for ix, data in phecode_pubmed.iterrows():
                phe = data["phewas_code"]
                tmp.write("%s\n" % phe)
                if phe in reg["PheWAS Code"].values:
                    phe_ix = reg["PheWAS Code"] == phe
                    phe_set = set(literal_eval(data["IdsList"]))
                    reg.loc[phe_ix, "phe.PM.count"] = len(phe_set)
                    joint_set = dx_set.intersection(phe_set)
                    reg.loc[phe_ix, "joint.PM.count"] = len(joint_set)
                    if len(dx_set) != 0:
                        reg.loc[phe_ix, "P.PM.phe.given.AD"] = len(joint_set) / len(dx_set)
                    if len(phe_set) != 0:
                        reg.loc[phe_ix, "P.PM.AD.given.phe"] = len(joint_set) / len(phe_set)

    for j in tqdm(['','2','3']):
        fname = "phewas_counts_missed%s.csv" % j
        if fname in pubmed_file_list:
            try:
                tmp.write("%s\n" % fname)
                phecode_pubmed = pd.read_csv(osp.join(pubmed_dir, fname)) # , dtype={"phewas_code": str})
            except Exception as e:
                ef = open("novelty_pubmed_errors.csv", 'a+')
                ef.write("error opening%s,%s" % (fname, e.args[0]))
                ef.close()

            for ix, data in phecode_pubmed.iterrows():
                phe = data["phewas_code"]
                tmp.write("%s\n" % phe)
                if phe in reg["PheWAS Code"].values:
                    phe_ix = reg["PheWAS Code"] == phe
                    phe_set = set(literal_eval(data["IdsList"]))
                    reg.loc[phe_ix, "phe.PM.count"] = len(phe_set)
                    joint_set = dx_set.intersection(phe_set)
                    reg.loc[phe_ix, "joint.PM.count"] = len(joint_set)
                    if len(dx_set) != 0:
                        reg.loc[phe_ix, "P.PM.phe.given.AD"] = len(joint_set) / len(dx_set)
                    if len(phe_set) != 0:
                        reg.loc[phe_ix, "P.PM.AD.given.phe"] = len(joint_set) / len(phe_set)

    tmp.close()
    reg.to_csv(outfile,index=False)


if __name__ == '__main__':
    main()
