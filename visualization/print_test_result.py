
import os.path
from glob import glob

from collections import Counter
from matplotlib import pyplot as plt
from textwrap import wrap
from numpy import arange
from pandas import read_csv, DataFrame

from library.clustering_test import *
import seaborn as sns

from visualization.visualize_utils import change_column_name

pd.set_option("display.max_rows", None, "display.max_columns", None, "display.max_colwidth", None)
pd.set_option('display.expand_frame_repr', False)
color_list = ['brown', 'orange', 'darkkhaki', 'limegreen', 'aqua', 'grey', 'blue', 'indigo', 'olive', 'magenta']


#color_list.extend(['k', ''])
x = [round(a,1) for a in arange(0, 1.1, 0.1)]

def save_score(testpath,score_dir, prep, score_method, test):
    testcase_dict = dict()
    resultlist = [a for a in glob(testpath + '/*'+prep+'*Score*'+score_method+'*.csv')]
    for result in resultlist:
        #testcase = result.split('-')[-4]
        testcase = result[result.rfind(test.replace('_', '-')) + len(test) + 1:result.rfind(prep) - 1]
        module = result.split('.')[-4]
        if module == 'movics':
            continue
        df = pd.read_csv(result, index_col=0).rename(columns={'Score': module})
        if testcase in testcase_dict:
            testcase_dict[testcase] = pd.concat([testcase_dict[testcase], df], axis=1)
        else:
            testcase_dict[testcase] = df

    for testcase, v in testcase_dict.items():
        path = score_dir + '.'.join([testcase, prep, score_method])
        v.to_csv(path)

def save_f1_score(testpath, prep, clinical_df, score_dir, test):
    # f-measure
    testcase_dict = dict()
    clusterlist = [a for a in glob(testpath + '/*' + prep + '*.cluster.csv')]
    testcaselist = sorted(set([gettestcase(a, test, prep) for a in clusterlist]))
    print(testcaselist)
    for testcase in testcaselist:
        path = score_dir + '.'.join([testcase, prep, 'f_measure'])
        if os.path.exists(path):
            continue
        print('testcase: ' + testcase)
        testcase_df = pd.DataFrame()
        for module in TOTAL_TOOL:
            if module == 'monti':
                continue
            filelist = sorted(
                [a for a in clusterlist if gettestcase(a, test, prep) == testcase and a.rfind(module) > 0])
            flist = []
            if len(filelist) != 10:
                flist = [np.NAN for a in range(10)]
            else:
                for file in filelist:
                    # print(file)
                    target = file.split('/')[-1].split('_')[2]
                    df_pred = pd.read_csv(file, index_col=0)
                    df_label = clinical_df[target]
                    df_label = df_label[df_pred.index]
                    flist.append(f_measure(df_label, df_pred.iloc[:,0]))
            testcase_df[module] = flist
        testcase_df.loc['mean'] = testcase_df.mean()
        testcase_df.loc['max'] = testcase_df.max()
        testcase_df.loc['min'] = testcase_df.min()
        testcase_dict[testcase] = testcase_df
        testcase_df.to_csv(path)

def gettestcase(a, test, prep):
    return a[a.rfind(test.replace('_','-')) + len(test) + 1:a.rfind(prep)-1]

def calc_f1_score():
    for cancer in TEST_CANCERS:
        print("**********"+cancer+"**********")
        testdir = join(DATA_DIR, cancer, BASE_DIR)
        testlist = getdirlist(testdir)
        clinical = join(DATA_DIR, cancer, 'clinical', cancer + '.csv')
        clinical_df = pd.read_csv(clinical, index_col=0)
        print(clinical_df)
        exit()
        for test in testlist:
            print('-------' + test + '-------')
            testpath = join(testdir, test)
            preplist = ['raw', 'norm', 'gene-centric']
            score_method_list = ['ari', 'nmi']
            for prep in preplist:
                score_dir = OUT_BASE_DIR + 'Score/' + test + '/' + cancer + '/' + prep + '/' + 'f_measure/'
                Path(score_dir).mkdir(parents=True, exist_ok=True)
                #f-measure
                save_f1_score(testpath, prep, clinical_df, score_dir, test)
                #NMI, ARI
                for score_method in score_method_list:
                    score_dir = OUT_BASE_DIR + 'Score/' + test + '/' + cancer + '/' + prep + '/' + score_method + '/'
                    Path(score_dir).mkdir(parents=True, exist_ok=True)
                    save_score(testpath,score_dir, prep, score_method, test)

def dict2df(score_dict, test, cancer):
    score_df = DataFrame()
    for testcase in score_dict.keys():
        for tool in sorted(score_dict[testcase].columns):
            df = DataFrame()
            df['Score'] = score_dict[testcase][tool]
            df['tool'] = tool
            df['testcase'] = testcase
            score_df = score_df.append(df, ignore_index=True)

    if test.find('combination') > -1:
        score_df.sort_values(by=['testcase'], key=lambda col: col.str.len(), inplace=True)
    elif test == 'balance_test':
        score_df.sort_values(by=['testcase'],
                                                key=lambda col: col.str.split('-', expand=True)[0].astype(int),
                                                inplace=True)
    else:
        score_df.sort_values(by=['testcase'], key=lambda col: col.astype(int), inplace=True)
    score_df['cancer'] = cancer

    return score_df

def violin_plot():
    for prep in ['raw', 'norm', 'gene-centric']:
        print(prep)
        for test in getdirlist(CLUSTERING_SCORE_DIR):
            print(test)
            plotdir_png = join(OUT_BASE_DIR, 'Violin_plot', prep)
            Path(plotdir_png).mkdir(parents=True, exist_ok=True)
            print('*' * 10 + test + '*' * 10)
            test_dir = join(CLUSTERING_SCORE_DIR, test)
            cluster_score_for_test_df = pd.DataFrame()
            ari_for_test_df = pd.DataFrame()
            nmi_for_test_df = pd.DataFrame()
            f_measure_for_test_df = pd.DataFrame()
            for cancer in sorted(getdirlist(test_dir)):
                #print(cancer)
                score_dir = join(test_dir, cancer, prep)
                score_dict = dict()
                for score in getdirlist(score_dir):
                    file_dir = join(score_dir, score)
                    file_list = getfilelist(file_dir)
                    #print(Score)
                    for file in file_list:
                        testcase = file.split('.')[0]
                        file_path = join(file_dir, file)
                        score_df = pd.read_csv(file_path, index_col=0)
                        score_df = score_df[:-3]
                        if testcase not in score_dict:
                            score_dict[testcase] = dict()
                        score_dict[testcase][score] = score_df

                cluster_score_dict = dict()
                ari_dict = dict()
                nmi_dict = dict()
                f_measure_dict = dict()
                for testcase in score_dict.keys():
                    cluster_score_mean = (score_dict[testcase]['ari'] + score_dict[testcase]['nmi'] + score_dict[testcase]['f_measure'])/3
                    cluster_score_dict[testcase] = round(cluster_score_mean, 3)
                    ari_dict[testcase] = round(score_dict[testcase]['ari'], 3)
                    nmi_dict[testcase] = round(score_dict[testcase]['nmi'], 3)
                    f_measure_dict[testcase] = round(score_dict[testcase]['f_measure'], 3)

                cluster_score_for_cancer_df = dict2df(cluster_score_dict, test, cancer)
                ari_for_cancer_df = dict2df(ari_dict, test, cancer)

                nmi_for_cancer_df = dict2df(nmi_dict, test, cancer)
                f_measure_for_cancer_df = dict2df(f_measure_dict, test, cancer)

                cluster_score_for_test_df = cluster_score_for_test_df.append(cluster_score_for_cancer_df)
                ari_for_test_df = ari_for_test_df.append(ari_for_cancer_df)
                nmi_for_test_df = nmi_for_test_df.append(nmi_for_cancer_df)
                f_measure_for_test_df = f_measure_for_test_df.append(f_measure_for_cancer_df)

            sns.set(font_scale=1.4)
            fig = plt.figure(figsize=(40, 10))
            sns.set(rc={"axes.unicode_minus": False, 'figure.figsize': (40, 10)}, style='whitegrid')

            xlabel = 'Test case'
            if test == SAMPLE_SIZE_TEST:
                xlabel = 'Number of samples per class'
            elif test == BALANCE_TEST:
                xlabel = 'Number of samples per class separated by \'-\''
            elif test == ROBUSTNESS_TEST:
                xlabel = 'Ratio of noise'
            elif test == FEATURE_SELECTION_TEST:
                xlabel = 'Ratio of feature selection'
            elif test.find('omics_combination')!= -1:
                xlabel = 'Omics combination'
            elif test == SUBTYPE_COMBINATION_TEST:
                xlabel = 'Subtype combination'


            for i, cancer in enumerate(sorted(getdirlist(test_dir))):
                ax = plt.subplot(2, 5, i+1)
                g = sns.violinplot(x='testcase', y='Score', hue=None, orient='v', data=cluster_score_for_test_df[cluster_score_for_test_df['cancer'] == cancer], palette='flare', cut=0)
                #mean
                sns.pointplot(x='testcase', y='Score',ax=ax, data=cluster_score_for_test_df[cluster_score_for_test_df['cancer'] == cancer], ci=None, color='black')
                legend_n = 1
                loc='best'
                if test == 'balance_test':
                    #ari
                    sns.pointplot(x='testcase', y='Score',ax=ax,label='ARI', data=ari_for_test_df[ari_for_test_df['cancer'] == cancer], ci=None, color='blue')
                    #nmi
                    sns.pointplot(x='testcase', y='Score',ax=ax,label='NMI', data=nmi_for_test_df[nmi_for_test_df['cancer'] == cancer], ci=None, color='green')
                    #f-measure
                    sns.pointplot(x='testcase', y='Score',ax=ax,label='F-measure', data=f_measure_for_test_df[f_measure_for_test_df['cancer'] == cancer], ci=None, color='red')
                    legend_n = 4
                    loc='lower right'
                if i == 0:
                    for curve, label in zip(ax.collections[-legend_n:], ['Cluster-Score', 'ARI', 'NMI', 'F-measure']):
                        curve.set_label(label)
                    ax.legend(title='Average Score',fancybox=True, framealpha=0.5, loc=loc, fontsize=7, title_fontsize=9)




                g.set(xlabel=None, ylabel=None)
                g.set_title(cancer)
                #g.set_titles(template='{col_name}')
                xlabels = g.get_xticklabels()

                g.set(ylim=(0, 1))
                g.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1])
                g.set_yticklabels(['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'])

                if test in ['feature_selection_test', 'robustness_test']:
                    for t in xlabels:
                        t.set_text(t.get_text() + '%')
                g.set_xticklabels(xlabels, rotation=40, fontsize=12, fontweight='bold')
                plt.tight_layout()


            if test=='balance_test':
                plt.subplots_adjust(bottom=0.16, left=0.025)
                fig.supxlabel(xlabel, y=0.02, fontweight='bold', fontsize='x-large')
            else:
                plt.subplots_adjust(bottom=0.13, left=0.025)
                fig.supxlabel(xlabel, y=0.02, fontweight='bold', fontsize='x-large')

            fig.supylabel('Cluster Score', x=0.005, fontweight='bold', fontsize='x-large')
            #plt.savefig(join(plotdir_pdf, test) + '.violinplot.pdf')

            plt.savefig(join(plotdir_png, test) + '.violinplot.png')


def data_load():
    print('Data loading...')
    result_filename = CLUSTERING_SCORE_DIR + 'total.p'
    if os.path.exists(result_filename):
       with open(result_filename, 'rb') as f:
           data = pickle.load(f)
           return data

    testdirs = getdirlist(CLUSTERING_SCORE_DIR)
    test_dict = dict()
    for test in testdirs:
        print('*'*10 + test +'*'*10)
        cancerdir = join(CLUSTERING_SCORE_DIR, test)
        cancerlist = getdirlist(cancerdir)
        cancer_dict = dict()
        for cancer in cancerlist:
            print('-' * 10 + cancer + '-' * 10)
            prepdir = join(cancerdir, cancer)
            preplist = getdirlist(prepdir)
            prep_dict = dict()
            for prep in preplist:
                print('#' * 10 + prep + '#' * 10)
                scoredir = join(prepdir, prep)
                scorelist = getdirlist(scoredir)
                score_dict = dict()
                for score in scorelist:
                    testcase_dir = join(scoredir, score)
                    testcase_list = getfilelist(testcase_dir)
                    testcase_df = pd.DataFrame()
                    print(testcase_list)
                    for testcase in testcase_list:
                        testcase_path = join(testcase_dir, testcase)
                        df = pd.read_csv(testcase_path, index_col=0)
                        df = df.loc['mean']
                        tc = testcase.split('.')[0]
                        df.name = int(tc) if tc.isnumeric() else tc
                        testcase_df = testcase_df.append(df)
                    testcase_df.sort_index(inplace=True)
                    score_dict[score] = testcase_df
                prep_dict[prep] = score_dict
            cancer_dict[cancer] = prep_dict
        test_dict[test] = cancer_dict

    with open(result_filename, 'wb') as f:
        pickle.dump(test_dict, f, pickle.HIGHEST_PROTOCOL)
    return test_dict

def best_case(df):
    cols = ['sample_size_test', 'feature_selection_test', 'balance_test', 'robustness_test', 'subtype_combination_test', 'omics_combination_test', 'omics_combination_test_gender', 'omics_combination_test_stage', 'omics_combination_test_age']
    preplist = ['raw', 'norm', 'gene-centric']
    total_mean_acc_dict = dict()
    for prep in preplist:
        best_testcase_result = pd.DataFrame(columns=cols)
        worst_testcase_result = pd.DataFrame(columns=cols)
        testcase_mean_acc_result = pd.DataFrame(columns=cols)
        best_tool_result = pd.DataFrame(columns=cols)
        best_tool_freq_total_dict = dict()
        for testname, t_v in df.items():
            #print('**************' + testname + '**************')
            bestcase_se = pd.Series()
            bestcase_se.name = testname
            worstcase_se = pd.Series()
            worstcase_se.name = testname

            mean_acc_se = pd.Series()
            mean_acc_se.name = testname
            besttool_se = pd.Series()
            besttool_se.name = testname

            for cancer, c_v in t_v.items():
                #print('-----------' + cancer +'-------------')
                #print(c_v)
                p_v = c_v[prep]
                cluster_score_df = (p_v['ari'] + p_v['nmi'] + p_v['f_measure'])/3
                if prep == 'norm':
                    p_v['ari']['Mean'] = p_v['ari'].mean(axis=1)
                    p_v['ari'].to_csv(join(CLUSTERING_SCORE_DIR, testname, cancer) + '/ARI.' + prep + '.csv')
                    p_v['nmi']['Mean'] = p_v['nmi'].mean(axis=1)
                    p_v['nmi'].to_csv(join(CLUSTERING_SCORE_DIR, testname, cancer) + '/NMI.' + prep + '.csv')
                    p_v['f_measure']['Mean'] = p_v['f_measure'].mean(axis=1)
                    p_v['f_measure'].to_csv(join(CLUSTERING_SCORE_DIR, testname, cancer) + '/F_MEASURE.' + prep + '.csv')

                cluster_score_df.rename(columns=CHANGE_TOOL_NAME, inplace=True)
                cluster_score_df = cluster_score_df.dropna(axis = 1)
                cluster_score_df['mean'] = cluster_score_df.mean(axis=1)

                besttool_str = '-'
                bestcase_str = '-'
                worstcase_str = '-'
                acc_max = 0
                mean_df = cluster_score_df['mean']

                if len(cluster_score_df) != 0:
                    acc_max = mean_df.max()
                    acc_min = mean_df.min()

                    maxidx_sr = mean_df[mean_df == acc_max]
                    minidx_sr = mean_df[mean_df == acc_min]

                    bestcase_str = ','.join([str(a) for a in maxidx_sr.index]) + '(' + str(round(acc_max,2)) + ')'
                    worstcase_str = ','.join([str(a) for a in minidx_sr.index]) + '(' + str(round(acc_min,2)) + ')'

                    testcase_df = cluster_score_df.drop(columns=['mean'])
                    testcase_len = len(testcase_df.index)
                    #testcase_df['min'] = cluster_score_df.min(axis=1)
                    testcase_df['Best'] = testcase_df.eq(testcase_df.max(axis=1), axis=0).apply(lambda x: ','.join(testcase_df.columns[x]), axis=1)
                    testcase_df['Mean'] = testcase_df.mean(axis=1)
                    #testcase_df['Worst'] = testcase_df.eq(testcase_df.min(axis=1), axis=0).apply(lambda x: ','.join(testcase_df.columns[x]), axis=1)
                    testcase_df.to_csv(join(CLUSTERING_SCORE_DIR, testname, cancer) + '/cluster_score.' + prep + '.csv')
                    #print(testcase_df)
                    best_freq = pd.Series(Counter(flatten([a.split(',') for a in testcase_df['Best']])))
                    best_freq = best_freq[best_freq == best_freq.max()]
                    besttool_list = [k for k, v in best_freq.items()]
                    for t in besttool_list:
                        if t in best_tool_freq_total_dict:
                            best_tool_freq_total_dict[t] = best_tool_freq_total_dict[t] + 1
                        else:
                            best_tool_freq_total_dict[t] = 1

                    besttool_str = ','.join(besttool_list) + '(' + str(best_freq.max()) + '/' + str(testcase_len) + ')'
                bestcase_se[cancer] = bestcase_str
                worstcase_se[cancer] = worstcase_str
                mean_acc_se[cancer] = acc_max
                besttool_se[cancer] = besttool_str
            best_testcase_result[testname] = bestcase_se
            worst_testcase_result[testname] = worstcase_se
            testcase_mean_acc_result[testname] = mean_acc_se
            best_tool_result[testname] = besttool_se

        best_testcase_result.rename(columns=CHANGE_TEST_NAME, inplace=True)
        change_column_name(best_testcase_result)
        worst_testcase_result.rename(columns=CHANGE_TEST_NAME, inplace=True)
        change_column_name(worst_testcase_result)
        testcase_mean_acc_result.rename(columns=CHANGE_TEST_NAME, inplace=True)
        change_column_name(testcase_mean_acc_result)
        best_tool_result.rename(columns=CHANGE_TEST_NAME, inplace=True)
        change_column_name(best_tool_result)
        best_tool_freq_total_dict = {k: v for k, v in sorted(best_tool_freq_total_dict.items(), key=lambda item: item[1], reverse=True)}
        best_tool_freq_total_str = ', '.join([k + ':' + str(v) for k, v in best_tool_freq_total_dict.items()])
        if prep == 'norm':
            print('Best test case:')
            print(best_testcase_result)
            print('Worst test case:')
            print(worst_testcase_result)
            print('Best tool')
            print(best_tool_result)
            print('Best tool frequency')
            print(best_tool_freq_total_str)
        best_testcase_result.to_csv(CLUSTERING_SCORE_DIR + '/bestcase.' + prep + '.csv')
        worst_testcase_result.to_csv(CLUSTERING_SCORE_DIR + '/worstcase.' + prep + '.csv')
        best_tool_result.to_csv(CLUSTERING_SCORE_DIR + '/besttool.' + prep + '.csv')
        total_mean_acc_dict[prep] = testcase_mean_acc_result

    total_mean_acc_dict = {v:total_mean_acc_dict[k] for k,v in CHANGE_PREP_NAME.items()}
    prep_compare = pd.DataFrame(data=None, columns=total_mean_acc_dict['Norm'].columns, index=total_mean_acc_dict['Norm'].index)

    prep_freq_list = []
    for i in prep_compare.index:
        for c in prep_compare.columns:
            ps = pd.Series({k:v.loc[i,c] for k,v in total_mean_acc_dict.items()})
            best_prep_str = np.nan
            if ps.max() != 0:
                best_prep_list = ps[ps == ps.max()].index.tolist()
                prep_freq_list.append(best_prep_list)
                if len(best_prep_list) == 3:
                    best_prep_str = 'All'
                else:
                    max_min = round(ps.max() - ps.mean(), 2)
                    if max_min == 0:
                        max_min = round(ps.max() - ps.mean(), 3)
                    best_prep_str = ','.join(best_prep_list) + '(+' + str(max_min) + ')'

            prep_compare.loc[i,c] = best_prep_str
    change_column_name(prep_compare)
    prep_freq_result = ', '.join([k + ':' + str(v) for k,v in Counter(flatten(prep_freq_list)).items()])
    prep_compare.to_csv(CLUSTERING_SCORE_DIR + '/best_preprocessing.csv')
    print(prep_compare)
    print(prep_freq_result)

def inner_plot(ax, test, cancer, df):
    #rtest = test.replace('_', '-')
    #df.index = [a.replace(rtest+'-','').replace('-'+prep,'') for a in df.index]

    df.index = ['\n'.join(wrap(str(l), 7)) for l in df.index]
    if test in ['sample_size_test', 'feature_selection_test', 'robustness_test']:
        df.index = df.index.astype(int)
    if test not in ['balance_test']:
        df = df.sort_index()
    if test in ['group_size_test', 'balance_test'] or test.find('combination') != -1:
        df = df.reindex(sorted(df.index,key=lambda d: (len(d), d)))
    if test == 'subtype_combination_test':
        df.sort_index(key=lambda x: x.str.len(), inplace=True)

    # if cancer=='LUAD':
    #     df.index = [a.replace('\n', '') for a in df.index]

    plt.sca(ax)
    x = df.index
    xi = range(len(x))
    mean_df = df.mean(axis=1)

    for i, c in enumerate(df.columns):
        ax.plot(xi, df[c], color=color_list[i], label=c, marker='o', linewidth=2, markersize=5, mec='black', alpha=0.5)

    ax.plot(range(len(mean_df)), mean_df, color='red', label='Mean', marker='D', linestyle='dashed', linewidth=2, markersize=5, mec='black')

    if test in ['feature_selection_test']:
        x = [str(a) + '%' for a in x]
    plt.xticks(xi, x)
    ax.title.set_text(cancer)
    ax.margins(x=0.025)

    # if (not test.endswith('age')) or (not test.endswith('stage')) or cancer == 'SKCM':
    #      plt.yticks(np.arange(0, 1.1, 0.1))

    plt.grid(True, axis='y', color='gray', alpha=0.5, linestyle='--')
    ax.tick_params(labelrotation=45, labelsize=8)

def line_plot_for_test():
    plotdir = join(OUT_BASE_DIR, 'Line_plot')
    for prep in ['raw', 'norm', 'gene-centric']:
        result_dir = join(plotdir, prep)
        Path(result_dir).mkdir(parents=True, exist_ok=True)
        print('-' * 10 + prep + '-' * 10)
        for test in getdirlist(CLUSTERING_SCORE_DIR):
            print('*'*10 + test + '*'*10)
            test_dir = join(CLUSTERING_SCORE_DIR, test)
            #output_dir = join(plotdir, test)
            fig = plt.figure(figsize=(18, 16))
            #fig.suptitle(test)
            fig.subplots(nrows=5, ncols=2)
            for cancer, a in zip(sorted(getdirlist(test_dir)), fig.axes):
                cluster_score_norm_path = join(test_dir,cancer,'cluster_score.'+prep+'.csv')
                df = pd.read_csv(cluster_score_norm_path, index_col=0).drop(columns=['Best', 'Mean'])
                inner_plot(a, test, cancer, df)
            if test == 'omics_combination_test_gender':
                fig.delaxes(fig.axes[-1])
            lines, labels = fig.axes[-1].get_legend_handles_labels()
            fig.legend(lines, labels, loc="lower center", ncol=5, mode='expand')
            fig.tight_layout()
            plt.subplots_adjust(bottom=0.12, left=0.035)
            fig.supxlabel('Testcase', y=0.06, fontweight='bold', fontsize='x-large')
            fig.supylabel('Cluster Score', x=0.005, fontweight='bold', fontsize='x-large')

            #plt.savefig(join(result_dir, test + '.' + prep + '.lineplot.pdf'))
            plt.savefig(join(result_dir, test + '.' + prep + '.lineplot.png'))

def average_testcase():
    testlist = getdirlist(CLUSTERING_SCORE_DIR)
    for test in testlist:
        # if test.find('omics_combination_test') != -1:
        #     continue
        print(test)
        test_dir = join(CLUSTERING_SCORE_DIR, test)
        for cancer in getdirlist(test_dir):
            cancer_dir = join(test_dir, cancer)
            for prep in getdirlist(cancer_dir):
                prep_dir = join(cancer_dir, prep)
                score_dict = dict()
                for score in CLUSTERING_SCORING_METHODS:
                    score_dir = join(prep_dir, score)
                    for filename in getfilelist(score_dir):
                        file_path = join(score_dir, filename)
                        testcase = filename.split('.')[0]
                        testcase_df = read_csv(file_path, index_col=0)
                        if testcase in score_dict:
                            score_dict[testcase][score] = testcase_df[:-3]
                        else:
                            score_dict[testcase] = dict()
                            score_dict[testcase][score] = testcase_df[:-3]

                cluster_score_dict = {k:pd.concat(v.values()).groupby(level=0).mean() for k,v in score_dict.items() if len(v.values()) == 3}
                sr_list = []
                for testcase, t_df in cluster_score_dict.items():
                    output_filename = join(prep_dir, testcase + '.' + prep + '.cluster_score.csv')
                    add_minmax(t_df)
                    t_df.to_csv(output_filename)
                    cluster_score_mean_sr = t_df.loc['mean']
                    cluster_score_mean_sr.name = testcase
                    sr_list.append(cluster_score_mean_sr)

                cluster_score_df = pd.concat(sr_list, axis=1)
                cluster_score_df = cluster_score_df.T
                cluster_score_df.rename(columns=CHANGE_TOOL_NAME, inplace=True)
                cluster_score_df = cluster_score_df.reindex(sorted(cluster_score_df.index))

                cluster_score_df['Best'] = cluster_score_df.eq(cluster_score_df.max(axis=1), axis=0).apply(
                    lambda x: ','.join(cluster_score_df.columns[x]), axis=1)
                cluster_score_df['Mean'] = cluster_score_df.mean(axis=1)
                cluster_score_df = cluster_score_df.round(3)

                cluster_score_df.sort_index(key=lambda x: x.str.len(), inplace=True)
                print(cluster_score_df)
                cluster_score_df.to_csv(join(cancer_dir, 'cluster_score.'+ prep +'.csv'))

def add_minmax(df):
    df.loc['max'] = df.max()
    df.loc['min'] = df.min()
    df.loc['mean'] = df.mean()

if __name__ == '__main__':
    #average_testcase()
    #calc_f1_score()
    #total_dict = data_load()
    #best_case(total_dict)
    #line_plot_for_test()
    violin_plot()
