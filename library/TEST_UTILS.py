from itertools import combinations

import os
from os import listdir
from os.path import join, isfile, isdir

import pickle

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import norm
import seaborn as sns

HOME_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_BASE_DIR = join(HOME_DIR, 'Result')
DATA_DIR = join(HOME_DIR, 'TCGA')

CLUSTERING_SCORE_DIR = join(OUT_BASE_DIR, 'Score')
CLUSTERING_RESULT_DIR = join(OUT_BASE_DIR, 'Clustering_results')
FEATURE_INFO_DIR = join(OUT_BASE_DIR, 'Feature_info')

TEST_CANCERS = ['BLCA', 'BRCA', 'COAD', 'HNSC', 'KIRP', 'LIHC', 'LUAD', 'SKCM', 'STAD','THCA']  # 10 종
TEST_OMICS = ['GE-ME', 'GE-MI', 'GE-CNV', 'GE-ME-MI', 'GE-ME-CNV', 'GE-MI-CNV', 'GE-ME-MI-CNV']
TEST_CLINICAL_FEATURES = ['Subtype', 'pathologic-stage', 'gender', 'Age']
BASIC_OMICS = 'GE-ME-MI'
MOVICS_TOOL=["PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF"]
TOTAL_TOOL= ["PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", 'mofa', 'iclusterplus','spectrum', 'snf']
CLUSTERING_TOOL = ["PINSPlus", "NEMO", "COCA", "LRAcluster", "ConsensusClustering", "IntNMF", 'mofa', 'iclusterplus','spectrum', 'snf']
TEST_ITERATION = 10
FEATURE_SELECTION_RATIO = 0.1
core_n = 10
BASE_DIR = 'Clustering_results'

ME = 'ME'
GE = 'GE'
MI ='MI'
CNV = 'CNV'

ARI = 'ari'
NMI = 'nmi'
F_MEASURE = 'f_measure'
CLUSTERING_SCORING_METHODS = [ARI, NMI, F_MEASURE]
CHANGE_TOOL_NAME = {'ConsensusClustering': 'CC', 'iclusterplus': 'Iclst', 'mofa': 'MOFA', 'snf': 'SNF', 'spectrum': 'Spectrum', 'LRAcluster': 'LRAclst', 'PINSPlus': 'PP'}
CHANGE_PREP_NAME = {'raw': 'Raw', 'norm': 'Norm', 'gene-centric': 'GC'}
CHANGE_TEST_NAME = {'omics_dimension_test':'feature_selection_test', 'omics_combination_test':'omics_combination_test_subtype'}
OLD_OMICS_NAME = {GE:GE, ME:'METH450', MI:'MIRNA', CNV:'CNVT'}

SAMPLE_SIZE_TEST = 'sample_size_test'
FEATURE_SELECTION_TEST = 'feature_selection_test'
BALANCE_TEST = 'balance_test'
ROBUSTNESS_TEST = 'robustness_test'
SUBTYPE_COMBINATION_TEST = 'subtype_combination_test'
OMICS_COMBINATION_TEST_AGE = 'omics_combination_test_age'
OMICS_COMBINATION_TEST_GENDER = 'omics_combination_test_gender'
OMICS_COMBINATION_TEST_STAGE = 'omics_combination_test_stage'
OMICS_COMBINATION_TEST_SUBTYPE = 'omics_combination_test'

TEST_LIST = [SAMPLE_SIZE_TEST, FEATURE_SELECTION_TEST, BALANCE_TEST, ROBUSTNESS_TEST, SUBTYPE_COMBINATION_TEST, OMICS_COMBINATION_TEST_SUBTYPE, OMICS_COMBINATION_TEST_GENDER, OMICS_COMBINATION_TEST_STAGE, OMICS_COMBINATION_TEST_AGE]



def check_exist(testname, testcase, cancer, norm, gencentric):
    norm_str = 'raw'

    score_path_dict = dict()
    if norm:
        norm_str = 'norm'
    if gencentric:
        norm_str = 'gene-centric'

    for score in CLUSTERING_SCORING_METHODS:
        dir_path = join(CLUSTERING_SCORE_DIR, testname.replace('-', '_'), cancer, norm_str, score)
        score_path_dict[score] = join(dir_path, '.'.join([testcase, norm_str, score]))

    if ((os.path.exists(score_path_dict[NMI]) and os.path.exists(score_path_dict[ARI]) and os.path.exists(
            score_path_dict[F_MEASURE]))):
        return True
    else:
        return False


def getPDA(pdapath):
    f = open(pdapath, 'rb')
    pda = pickle.load(f)
    f.close()
    return pda

def odf_stats(pda):  # getting stats of omics data
    for o in pda.data_omics_level.keys():
        pda.data_n[o] = pda.data_omics_level[o].shape
    return pda

def split_from_filename(filename):
    cancer = filename.split('/')[5]
    testname = filename.split('/')[7]
    testcase = filename.split('_')[-2]
    norm = filename.split('_')[-1]
    return  cancer, testname, testcase, norm

def f_measure(label, pred):
    p = make_comb(label)
    q = make_comb(pred)
    a = p&q
    b = p-q
    c = q-p
    f = (2*len(a)) / (2*len(a) + len(b) + len(c))
    return f

def flatten(t):
    return [item for sublist in t for item in sublist]

def make_comb(df):
    l = df.reset_index(drop=True)
    label_c = [l[l == a].index.tolist() for a in df.value_counts().index.tolist()]
    total_combination = set()
    for clist in label_c:
        total_combination.update(list(combinations(clist, 2)))
    return total_combination

def add_norm_tag(testcase, norm_bool, gene_centric):
    if gene_centric == True:
        tag ='_gene-centric'
        testcase = testcase + tag
    elif norm_bool == True:
        tag ='_norm'
        testcase = testcase + tag
    else:
        tag = '_raw'
        testcase = testcase + tag
    return testcase, tag

def getfilelist(dir):
    return [f for f in listdir(dir) if isfile(join(dir, f))]

def getdirlist(dir):
    return [f for f in listdir(dir) if isdir(join(dir, f))]



def plot_all_distributions(df, omics_name, prep, cancer):
    kwargs = dict(hist_kws={'alpha': .7}, kde_kws={'linewidth': 2})

    # Create a figure with multiple subplots
    plt.figure(figsize=(10, 7), dpi=80)
    plt.subplots_adjust(left=0.07, right=0.95, bottom=0.07, top=0.95)

    # Loop through all columns in the dataframe
    for i, column_name in enumerate(df.columns):
        # Select the current column from the dataframe
        column = df[column_name]

        # Plot the distribution of the current column
        sns.distplot(column, label=column_name, **kwargs, bins=20 )

    # Add a title and labels to the plot
    plt.title(cancer + '-' + prep + ": " + omics_name + '\'s ramdom 10 features')
    plt.xlabel("Value")
    plt.ylabel("Density")
    #plt.xlim(0, 150)
    plt.legend()

    # Display the plot
    plt.show()