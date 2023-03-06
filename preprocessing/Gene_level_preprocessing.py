# Util functions for MOBD management
import pickle

import os
import pandas as pd
import numpy as np
from library.TEST_UTILS import DATA_DIR, GE, ME, MI, CNV

METADIR = '%s/annotation' % (os.path.dirname(os.path.abspath(__file__)))

# PDA() object
class PDA:
    def __init__(self):
        self.project = ''  # project tpye, eg) TCGA
        self.cancer_type = ''  # cancer type, eg) COAD
        self.clinical_feature = ''  # target clinical feature, eg) Molecular_Subtype
        self.clinical_feature_labels = None  # clinical tables to be used for analysis
        self.clinical_feature_groups = []  # clincial feature groups eg)
        self.clinical_feature_groups_n = 0
        self.conditions = []  # TBD
        self.omics_string = None
        self.aggregation_type = 'matched'  # 'any' | 'matched' (default) - 'any': as many samples possible, 'matched': common samples in all omics data
        self.omics_types = []  # list of omics types
        self.data_omics_level = {}  # omics level data
        self.data_gene_level = {}  # gene level data
        self.data_gene_level_norm = {}  # normalized gene level data
        self.data_n = {}  # number of samples in each omics
        self.normalization_info = {}  # normalization method per omics type --> TBD
        self.samplelist = []  # list of barcodes used for analysis
        self.genelist = []  # list of gene IDs (ENSGID) used for analysis


def gsymToendgid(dat_o, o):
    gsym_dict = {}
    gsym_set = set()
    for line in open('%s/geneid/ensgid_info.tsv' % (METADIR)):
        tok = line.strip().split('\t')
        ensgid = tok[0]
        enstid = tok[2]
        gsym = tok[4]
        gsym_dict[gsym] = (ensgid, enstid)
        gsym_set.add(gsym)

    matched_gsym = list(set(dat_o[o].index).intersection(gsym_set))
    dat_o[o] = dat_o[o].loc[matched_gsym]
    ensgid_index = [gsym_dict[i][0] for i in dat_o[o].index]
    dat_o[o].index = ensgid_index
    dat_o[o] = dat_o[o].groupby(dat_o[o].index).mean() #중복되는 id 평균값 처리
    return dat_o[o]

def preproc_omics(o, dat_o):
    if o == GE:
        gsymToendgid(dat_o, o)  ## Gene ID mapping
        return dat_o[o]
    elif o == ME:
        return dat_o[ME]
    elif o == MI:
        dat_o[o] = dat_o[o].fillna(0)
        # -- convert MIMAT IDs to miRNA IDs
        mirid_f = '%s/mirna/mirna_id_table.tsv' % (METADIR)
        mirid_dict = {}
        for line in open(mirid_f):
            mirid, mimatid = line.strip().split('\t')
            mirid_dict[mimatid] = mirid

        mirid_idx = []
        filter_out_mimats = []
        filtercnt = 0
        for idx in dat_o[MI].index.tolist():
            if idx not in mirid_dict:
                filter_out_mimats.append(idx)
                filtercnt += 1
                continue
            mirid_idx.append(mirid_dict[idx])
        # drop missing MIMAT IDs
        dat_o[MI].drop(filter_out_mimats, inplace=True)
        dat_o[MI].index = mirid_idx
        return dat_o[MI]
    elif o == CNV:
        dat_o[o] = dat_o[o].fillna(0)
        indexlist = [a.split('|')[0] for a in dat_o[o].index]
        dat_o[o].index = indexlist
        gsymToendgid(dat_o, o)
        return dat_o[o]

def convert_to_gene_level(pda):
    # load metadata
    gid_tid_file = '%s/geneid/gid_tid.txt' % (METADIR)  # ENSTID to ENSGID mapping information
    methylation_promoter_probes = '%s/methylation/promoter_probes_illumina450.txt' % (METADIR)  # promoter methylation probe information
    mirna_gene_target_file = '%s/mirna/mirna_target_gene.csv' % (METADIR)  # miRNA-target gene information

    # preparing gene-level data
    print('[%s-%s] converting to gene-centric data' % (pda.cancer_type, pda.clinical_feature))

    if GE in pda.data_omics_level.keys():    # gene level mRNA data
        pda.data_gene_level[GE] = pda.data_omics_level[GE].loc[pda.genelist]

    if ME in pda.data_omics_level.keys():# converting methylation data to gene-level
        pda.data_gene_level[ME] = make_methylation_gcentric(pda.data_omics_level[ME],
                                                            pda.genelist, gid_tid_file,
                                                            methylation_promoter_probes)
    if MI in pda.data_omics_level.keys():    # converting miRNA data to gene-level
        pda.data_gene_level[MI] = make_mir_gcentric(pda.data_omics_level[MI], pda.genelist, mirna_gene_target_file)
    if CNV in pda.data_omics_level.keys():
        pda.data_gene_level[CNV] = pda.data_omics_level[CNV].loc[pda.genelist]

    return pda

def make_methylation_gcentric(methf, genelist, gtid, methpromprob):
    meth_barcodes = methf.columns.tolist()
    methBetaValues = dict()
    print('reading methylation data...')
    for idx, val in methf.iterrows():
        tok = val.tolist()
        pid = idx  # probe id
        methBetaValues[pid] = [-1.0 if np.isnan(x) else float(x) for x in tok]

    # tid to gid converstion dictionary
    gid_dict = {}
    for line in open(gtid):
        gid, tid = line.strip().split('\t')
        gid_dict[tid] = gid

    print('reading TSS promoter probe annotation...')
    rf = open(methpromprob)
    targetRevDict = dict()
    while (True):
        line = rf.readline()
        if not line:
            break
        tok = line.strip().split('\t')
        pid = tok[0].strip()
        tid = tok[4].strip()
        if tid not in gid_dict:
            continue
        geneName = gid_dict[tid]
        if geneName not in targetRevDict:
            targetRevDict[geneName] = [pid]
        else:
            targetRevDict[geneName].append(pid)
    rf.close()

    meth_data = []

    def writeList(lst):
        row_data = []
        for i in range(len(lst)):
            row_data.append(float(lst[i]))
        meth_data.append(row_data)

    print('converting methylation to gene level data...')
    gene_list = []
    for geneNameIdx, geneName in enumerate(genelist):
        gene_list.append(geneName)
        if geneName not in targetRevDict:
            writeList([0.0 for _ in range(len(meth_barcodes))])
            continue

        pid = targetRevDict[geneName]
        pid = list(filter(lambda x: x in methBetaValues, pid))
        lst = []
        for i in range(len(meth_barcodes)):
            methAvg = 0.0
            num = 0.0
            for j in range(len(pid)):
                if methBetaValues[pid[j]][i] == -1: continue  # skip NA values
                num += 1.0
                methAvg = methAvg * (1.0 - 1.0 / num) + methBetaValues[pid[j]][i] / num  # running average

            lst.append(float(methAvg))

        writeList(lst)

    return pd.DataFrame(meth_data, columns=meth_barcodes, index=gene_list)


def make_mir_gcentric(mirf, genelist, mirtarget):
    # sep=','
    print('reading mirna file...')
    # rf = open(mirf)
    # mi_barcodes = rf.readline().strip().split(sep)[1:]
    mi_barcodes = mirf.columns.tolist()
    miExprs = dict()
    for idx, val in mirf.iterrows():
        token = val.tolist()
        miName = idx
        miExprs[miName] = list(map(lambda x: float(x), token))

    miNames = miExprs.keys()
    for miName in miNames:
        newLst = []
        oldLst = miExprs[miName]
        for barcode in mi_barcodes:
            newLst.append(oldLst[mi_barcodes.index(barcode)])
        miExprs[miName] = newLst

    print('reading mirna gene target file...')
    rf = open(mirtarget)
    targetRevDict = dict()
    while (True):
        line = rf.readline()
        if not line:
            break
        token = line.strip().split('\t')
        miName = token[0].strip()
        geneName = token[1].strip()
        if geneName not in targetRevDict:
            targetRevDict[geneName] = [miName]
        else:
            targetRevDict[geneName].append(miName)
    rf.close()

    mir_data = []

    def writeList(lst):
        row_data = []
        for i in range(len(lst)):
            row_data.append(float(lst[i]))
        mir_data.append(row_data)

    print('converting miRNAs to gene level data...')
    gene_list = []
    for geneNameIdx, geneName in enumerate(genelist):
        gene_list.append(geneName)
        if geneName not in targetRevDict:
            writeList([0.0 for _ in range(len(mi_barcodes))])
            continue

        miNames = targetRevDict[geneName]
        miNames = list(filter(lambda x: x in miExprs, miNames))
        lst = []
        for i in range(len(mi_barcodes)):
            miAvg = 0.0
            num = 0.0
            for j in range(len(miNames)):
                num += 1.0
                miAvg = miAvg * (1.0 - 1.0 / num) + miExprs[miNames[j]][i] / num  # running average

            lst.append(miAvg)
        arr = np.asarray(lst)
        lst = list(arr)
        writeList(lst)

    return pd.DataFrame(mir_data, columns=mi_barcodes, index=gene_list)
