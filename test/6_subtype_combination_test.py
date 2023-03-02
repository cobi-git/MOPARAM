from library.TEST_UTILS import TEST_CANCERS, TEST_CLINICAL_FEATURES, core_n, BASIC_OMICS, check_exist
from library.clustering_test import *

testname = 'subtype-combination-test'
proc = []
proc_excution = []
proc_end = []

def thread(cancer, pdapath, arg, gene_centric, norm, samplesize):
        testcase = '%s' % ('-'.join([str(a) for a in arg]))
        data_raw, labels, csvnames, filename, labelpath = \
            clustering_preprocessing(cancer, pdapath, testname, testcase, group=arg, gene_centric_bool=gene_centric, norm_bool=norm, sample_size=samplesize)
        clustering_test(data_raw, labels, csvnames, filename, labelpath)

for cancer in TEST_CANCERS:
    pdadir = join(DATA_DIR, cancer, 'PDA')
    pdalist = [join(pdadir, a) for a in getfilelist(pdadir) if a.split('_')[3] in TEST_CLINICAL_FEATURES]
    for pdapath in pdalist:
        if not pdapath.split('_')[5] == BASIC_OMICS or not (pdapath.split('_')[4] == 'Subtype'):
            continue
        pda = getPDA(pdapath)
        groups = [a for a in pda.clinical_feature_groups.index.to_list()]
        label_df = pda.clinical_feature_labels
        sample_size = len(label_df[label_df == groups[-1]])
        for i in range(2, len(groups) + 1, 1):
            #group = groups[0:i]
            for group in combinations(groups, i):
                group = list(group)
                #thread(cancer, pdapath, group, gene_centric=False, norm=True, samplesize=sample_size)
                p = Process(target=thread, args=(cancer, pdapath, group, False, False, sample_size)) # raw-not-norm
                proc.append(p)
                p = Process(target=thread, args=(cancer, pdapath, group, False, True, sample_size))  # raw-not-norm
                proc.append(p)
                p = Process(target=thread, args=(cancer, pdapath, group, True, True, sample_size))  # gene-centric-norm
                proc.append(p)

while len(proc) > 0:
    for i in proc:
        if len(proc_excution) < core_n:
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