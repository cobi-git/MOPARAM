#!/usr/bin/env python

import sys
from os import listdir
from pathlib import Path
import pickle
import os
from library.TEST_UTILS import TEST_OMICS, OLD_OMICS_NAME, DATA_DIR, GE, CNV, odf_stats, TCGA, pda_error, normalize
import pandas as pd
from os.path import join, isdir
from multiprocessing import Process

from preprocessing.Gene_level_preprocessing import preproc_omics, convert_to_gene_level, PDA


def makePDA(project, cancer, clinical, omics_string, matched):

    filedir = '%s/%s/PDA/' % (DATA_DIR,cancer)
    Path(filedir).mkdir(parents=True, exist_ok=True)
    filename = 'PDA_%s_%s_%s_%s_%s' % (project, cancer, clinical, omics_string, matched)

    path = join(filedir, filename + '.pda')
    if os.path.exists(path):
        print('***** Already exist ***** - ' + path)
        return

    print('Processing: ' + filename )
    pda = PDA()
    pda.project = project
    pda.cancer_type = cancer  # cancer type
    pda.clinical_feature = clinical  # target clinical feature
    pda.omics_types = omics_string.split('-')
    pda.omics_string = omics_string
    pda.aggregation_type = matched  # matched sample base (for MONTI and MOPA)

    # load clinical data
    clinical_dat = pd.read_csv('%s/%s/clinical/%s.csv' % (DATA_DIR, pda.cancer_type, pda.cancer_type),header=0, index_col='sampleID')
    cf_set = clinical_dat[pda.clinical_feature].dropna()  # .sort_values()
    cf_set = cf_set[cf_set != '-']

    # -------- DEBUG START --------# ## inukj: need to fix the ID's of either clinical or BDA/PDA
    new_barcodes = []
    for i in cf_set.index:
        new_barcodes.append(i)

    cf_set.index = new_barcodes
    # -------- DEBUG END   --------#

    # load omics data
    # parsing CDA files
    data_omics_level = {omics: [] for omics in pda.omics_types}
    common_samples = set()
    for o in data_omics_level.keys():
        print('loading %s omics data...' % (o))

        data_omics_level[o] = pd.read_csv('%s/%s/CDA/CDA_TCGA_%s_%s_count.csv' %
                                          (DATA_DIR, pda.cancer_type, pda.cancer_type, OLD_OMICS_NAME[o]), header=0,index_col='sample')
        if not len(common_samples):
            common_samples = data_omics_level[o].columns
        else:
            common_samples = common_samples.intersection(set(data_omics_level[o].columns))

    ## Sample ID matching
    # select matching barcodes
    matched_samples = set(cf_set.index).intersection(common_samples)

    # matched clinical feature
    pda.clinical_feature_labels = cf_set.loc[matched_samples].sort_values()
    pda.clinical_feature_groups = pda.clinical_feature_labels.value_counts()

    # filtering: clinical feature groups with size < max*0.2
    g_size= pda.clinical_feature_groups.max() * 0.2

    for idx, i in pda.clinical_feature_groups.items():
        if i < g_size and i < 10:
            pda.clinical_feature_labels=pda.clinical_feature_labels[pda.clinical_feature_labels != idx]

    pda.samplelist = pda.clinical_feature_labels.index.tolist()

    if len(matched_samples) == 0:
        pda_error('***** No matched samples *****', pda)

    pda.clinical_feature_groups = pda.clinical_feature_labels.value_counts()
    pda.clinical_feature_groups_n = len(pda.clinical_feature_groups)

    if pda.clinical_feature_groups_n < 2:
        pda_error('***** too few group size *****', pda)

    for o in data_omics_level.keys():
        data_omics_level[o] = data_omics_level[o][pda.clinical_feature_labels.index]
        if data_omics_level[o].empty:
            pda_error('***** too few samples in the group *****', pda)

    #preprocessing
    for o in data_omics_level.keys():
        pda.data_omics_level[o]=preproc_omics(o, data_omics_level)

    if CNV in data_omics_level.keys():
        pda.genelist = list(set(pda.data_omics_level[GE].index).intersection(set(pda.data_omics_level[CNV].index)))
    else:
        pda.genelist = list(pda.data_omics_level[GE].index)

    # filter NaN omics values, cond: NaN<80%
    for o in pda.data_omics_level.keys():
        pda.data_omics_level[o] = pda.data_omics_level[o].loc[pda.data_omics_level[o].isnull().mean(axis=1) < .8, :]
    for o in pda.data_omics_level.keys():
        if len(pda.data_omics_level[o]) == 0:
            pda_error('***** too much NAN *****', pda)

    # Finalizing PDA stats
    pda = odf_stats(pda)
    pda = convert_to_gene_level(pda) # convert omics-level to gene-level data
    pda = normalize(pda) # normalize omics data

    # save pda object - a single object with packed omics, clinical and gene list data
    print("Saving PDA object...")
    filename = '%s/PDA_%s_%s_%s_%s_%s_%d.pda' % (filedir,
                                                pda.project,
                                                pda.cancer_type,
                                                pda.clinical_feature,
                                                pda.omics_string,
                                                pda.aggregation_type,
                                                len(pda.samplelist))

    #inputPDAtoDB(project, cancer, clinical, omics_string, MATCHED, filename, pda.samplelist)
    pda_out = open(filename, 'wb')
    pickle.dump(pda, pda_out)
    pda_out.close()

if __name__ == "__main__":
    matched = 'matched'
    cancer_types = [f for f in listdir(DATA_DIR) if isdir(join(DATA_DIR, f))]
    proc = []
    proc_excution = []
    proc_end = []

    for cancer_type in cancer_types:
        print('------' + cancer_type + '------')
        clinical_features = pd.read_csv('%s/%s/clinical/%s.csv' % (DATA_DIR, cancer_type, cancer_type), header=0,	index_col='sampleID').columns.tolist()

        for clinical in clinical_features:
            if not clinical in ['pathologic-stage', 'gender', 'Age', 'Subtype']:
                continue
            print(clinical)
            for omics_string in TEST_OMICS:
                #makePDA(TCGA, cancer_type, clinical, omics_string, matched)
                p = Process(target=makePDA, args=(TCGA, cancer_type, clinical, omics_string, matched,))
                proc.append(p)

    thread_n = 20
    while len(proc) > 0:
        for i in proc:
            if len(proc_excution) < thread_n:
                p = proc.pop(0)
                p.start()
                proc_excution.append(p)
            else:
                break
        for excution_proc in proc_excution:
            if not excution_proc.is_alive():
                proc_excution.remove(excution_proc)
                proc_end.append(excution_proc)
                excution_proc.join()
                excution_proc.close()
else:
    project = sys.argv[1]
    cancer_type = sys.argv[2]
    clinical_feature =sys.argv[3]
    omics_string = sys.argv[4]
    aggregation_type = sys.argv[5]
    makePDA(project, cancer_type, clinical_feature, omics_string, aggregation_type)
