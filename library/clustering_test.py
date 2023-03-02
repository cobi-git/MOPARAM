import csv
import os.path
import subprocess
from multiprocessing import Pool, Process
from os.path import join
from pathlib import Path
import numpy as np
import pandas as pd

# snf
import qnorm
from numpy.random import uniform
from sklearn.feature_selection import SelectKBest, chi2, f_classif
from sklearn.preprocessing import MinMaxScaler
from library.TEST_UTILS import *
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import normalized_mutual_info_score

def clustering_preprocessing(cancer, pdapath, testname, testcase, sample_size=0, group=[],
                             feature_selection_ratio=FEATURE_SELECTION_RATIO, balance=[],
                             noise_ratio=0, gene_centric_bool=True, norm_bool=True):
    # if check_exist(testname, testcase, cancer, norm, gene_centric):
    #     print('completed!!!')
    #     exit()
    #print(testname)


    odir = join(CLUSTERING_RESULT_DIR, cancer, testname.replace("-", "_"))
    Path(odir).mkdir(parents=True, exist_ok=True)
    pda = getPDA(pdapath)
    testcase, norm_tag = add_norm_tag(testcase, norm_bool, gene_centric_bool)
    nametag = '%s/%s_%s_%s_%s_%s' % (
        odir, pda.project, pda.cancer_type, pda.clinical_feature, pda.omics_string, testcase)

    #print(cancer + ":" + nametag.split('/')[-1])

    label_df = pda.clinical_feature_labels
    gene_list_by_omics_dict = dict()
    using_data = None

    df = label_df.value_counts()
    pda.clinical_feature_groups = df[df >= 10]
    #print(pda.clinical_feature_groups)
    if sample_size == 0:
        sample_size = pda.clinical_feature_groups.min()

    if len(group) > 0:
        sample_size = len(label_df[label_df == pda.clinical_feature_groups.idxmin()])
    else:
        group = pda.clinical_feature_groups.index

    if len(balance) > 0:
        sample_size_list = balance
    else:
        sample_size_list = [sample_size for a in group]

    if gene_centric_bool == False:
        using_data = pda.data_omics_level
    elif gene_centric_bool == True:
        using_data = pda.data_gene_level

    for k, v in using_data.items():
        using_data[k] = v.dropna() #drop features which contain na
    #     if k == CNV:
    #         using_data[k] = pd.DataFrame(MinMaxScaler().fit_transform(v.T).T, columns=v.columns, index=v.index)

    if feature_selection_ratio < 1:
        for k, v in using_data.items():
            X = v.transpose()
            method = chi2
            if k == CNV:
                method = f_classif
            feature_n = len(v)

            c = 1
            if k == ME:
                if gene_centric_bool != True:
                    c = 0.1
            select_feature_n = int(c * feature_selection_ratio * feature_n)
            selector = SelectKBest(method, k=select_feature_n)
            selector.fit(X, label_df)
            features = selector.get_support(indices=True)
            gene_list_by_omics_dict[k] = X.iloc[:, features].columns
    elif feature_selection_ratio == 1:
        gene_list_by_omics_dict = {k: using_data[k].index for k, v in using_data.items()}

    #make feature list
    feature_dict = {k:v.tolist() for k,v in gene_list_by_omics_dict.items()}
    feature_out_dir = join(FEATURE_INFO_DIR, testname.replace("-", "_"), cancer, norm_tag[1:])
    Path(feature_out_dir).mkdir(parents=True, exist_ok=True)
    feature_path = '%s/%s_%s_%s_%s_%s' % (
        feature_out_dir, pda.project, pda.cancer_type, pda.clinical_feature, pda.omics_string,
        testcase) + ".features.csv"

    # Find the length of the longest list in the dictionary
    max_length = max(len(x) for x in feature_dict.values())

    # Pad each list with None values to make them the same length
    for key in feature_dict:
        feature_dict[key] += [None] * (max_length - len(feature_dict[key]))

    Path(FEATURE_INFO_DIR).mkdir(parents=True, exist_ok=True)
    with open(feature_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(feature_dict.keys())  # write the header row
        rows = zip(*feature_dict.values())
        writer.writerows(rows)
    exit()
    if noise_ratio > 0:
        for k, v in using_data.items():
            mask = np.random.choice([True, False], size=v.shape, p=[noise_ratio, 1 - noise_ratio])
            mean = v.values.mean()
            std = v.values.std()
            noise = np.random.normal(mean, std, size=v.shape)
            noise_df = pd.DataFrame(noise, columns=v.columns, index=v.index)
            using_data[k] = v.mask(mask).fillna(noise_df)

    omics_n = len(using_data.keys())
    labels_list = []
    csvnames_list = []
    labelpath_list = []

    for i in range(1, TEST_ITERATION + 1, 1):
        sampled_index = []
        for idx, label in enumerate(group):
            r_df = label_df[label_df == label].sample(sample_size_list[idx])  #
            sampled_index.extend(r_df.index)

        data = {k: v.transpose().loc[sampled_index, gene_list_by_omics_dict[k]] for k, v in using_data.items()}
        labels = pda.clinical_feature_labels[sampled_index]
        labels_list.append(labels)

        if norm_bool == True:#normalization
            data = normalization(data)

        labelpath = nametag + ".labels" + str(i) + ".csv"
        labelpath_list.append(labelpath)
        labels.to_csv(labelpath)
        csvnames = []
        for k, v in data.items():
            name = nametag + '.' + k + '-' + str(i)
            csvnames.append(name)
            v.to_csv(name)
        csvnames_list.append(csvnames)

    return omics_n, labels_list, csvnames_list, nametag, labelpath_list

def getClusteringScore(labels, clst_df):
    label_true = np.array(labels)
    label_pred = np.array(clst_df)
    nmi = normalized_mutual_info_score(label_true, label_pred)
    ari = adjusted_rand_score(label_true, label_pred)
    f = f_measure(labels, clst_df)
    return nmi, ari, f

def calc_mean(v):
    v.loc['mean'] = v.mean()
    v.loc['max'] = v.max()
    v.loc['min'] = v.min()

def clustfile2score(clust_file_list, label_df_list, tool, score_dict):
    ari_list = []
    nmi_list = []
    f_measure_list = []
    for clust_file, label_df in zip(clust_file_list, label_df_list):
        clust_df = pd.read_csv(clust_file, index_col=0)
        nmi, ari, f = getClusteringScore(label_df, clust_df['cluster'])
        ari_list.append(ari)
        nmi_list.append(nmi)
        f_measure_list.append(f)
    ari_df = pd.DataFrame(ari_list, columns=[tool])
    nmi_df = pd.DataFrame(nmi_list, columns=[tool])
    f_measure_df = pd.DataFrame(f_measure_list, columns=[tool])
    calc_mean(ari_df)
    calc_mean(nmi_df)
    calc_mean(f_measure_df)
    score_dict[ARI][tool] = ari_df
    score_dict[NMI][tool] = nmi_df
    score_dict[F_MEASURE][tool] = f_measure_df

def module_thread(csvnames, module_filename, labelpath, module):
    print(module + ' thread...' )
    modulepath = '/data3/projects/2021_MODB/python/library/R/%s.r' % (module)
    query = "/usr/bin/Rscript %s %s %s" % (modulepath, labelpath, module_filename)
    for csvname in csvnames:
        query = "%s %s" % (query, csvname)
    subprocess.call(query, shell=True) #run R scripts

def write_score(score_dict, module_filename, module, label_df_list):
    clust_file_list = [module_filename + str(i) + ".cluster.csv" for i in range(1, 11)]
    if module_filename.find('movics') > 0:
        clust_file_list = [module_filename + str(i) + "." + module + ".cluster.csv" for i in range(1, 11)]
    clustfile2score(clust_file_list,label_df_list,module, score_dict)

def run_module(csvnames_list, label_df_list, filename, labelpath_list, module, score_dict):
    print(module + '...')
    module_filename = filename + '.' + module
    with Pool(TEST_ITERATION) as p:
        args = []
        for i in range(1, 11):
            args.append((csvnames_list[i-1], module_filename + str(i), labelpath_list[i-1], module))
        p.starmap(module_thread, args)
    if module == 'movics':
        for tool in MOVICS_TOOL:
            write_score(score_dict, module_filename, tool, label_df_list)
    else:
        write_score(score_dict, module_filename, module, label_df_list)

def normalization(data_dict):
    for k, v in data_dict.items():
        # quantile normalize, If axis=1, standardize each sample (column), if axis=0, standardize each feature(row)
        df_norm = qnorm.quantile_normalize(v, axis=1, ncpus=10)
        df_norm = pd.DataFrame(MinMaxScaler().fit_transform(df_norm), columns=df_norm.columns, index=df_norm.index)
        data_dict[k] = df_norm
    return data_dict

def clustering_test(omics_n, labels_list, csv_names_list, filename, labelpath_list):
    cancer, testname, testcase, norm = split_from_filename(filename)
    score_dict = dict()
    score_path_dict = dict()

    for score in CLUSTERING_SCORING_METHODS:
        dir_path = join(CLUSTERING_SCORE_DIR, testname, cancer, norm, score)
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        score_path_dict[score] = join(dir_path, '.'.join([testcase, norm, score]))
        print(score_path_dict[score])
        score_dict[score] = pd.DataFrame()

    if not (os.path.exists(score_path_dict[NMI]) and os.path.exists(score_path_dict[ARI]) and os.path.exists(score_path_dict[F_MEASURE]) ):
        run_module(csv_names_list, labels_list, filename, labelpath_list, 'SNF', score_dict)
        if omics_n != 1: #movics can't integrate single-omics data
            run_module(csv_names_list, labels_list, filename, labelpath_list, 'movics', score_dict)
        run_module(csv_names_list, labels_list, filename, labelpath_list, 'MOFA', score_dict)
        run_module(csv_names_list, labels_list, filename, labelpath_list, 'Spectrum', score_dict)
        run_module(csv_names_list, labels_list, filename, labelpath_list, 'IclusterPlus', score_dict)
        for score in CLUSTERING_SCORING_METHODS:
            score_dict[score].to_csv(score_path_dict[score])
        print(score_dict[NMI])
    else:
        print(score_path_dict[NMI])
        print(score_path_dict[ARI])
        print(score_path_dict[F_MEASURE])
        print('already done!')
        # for Score in CLUSTERING_SCORING_METHODS:
        #     score_path_dict[Score] = join(clustering_score_dir, testname, cancer, norm, Score, '.'.join([testcase, norm, Score]))
        #     if os.path.exists(score_path_dict[Score]):
        #         score_dict[Score] = pd.read_csv(score_path_dict[Score], index_col=0)





