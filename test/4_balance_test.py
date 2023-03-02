from math import floor

from library.TEST_UTILS import TEST_CANCERS, TEST_CLINICAL_FEATURES, BASIC_OMICS, core_n, check_exist
from library.clustering_test import *

testname = 'balance-test'

proc = []
proc_excution = []
proc_end = []

def thread(cancer, pdapath, arg, gene_centric, norm):
        testcase = '%s' % ('-'.join([str(a) for a in arg]))

        data_raw, labels, csvnames, filename, labelpath = \
            clustering_preprocessing(cancer, pdapath, testname, testcase, balance=arg, gene_centric_bool=gene_centric, norm_bool=norm)
        clustering_test(data_raw, labels, csvnames, filename, labelpath)

for cancer in TEST_CANCERS:
    pdadir = join(DATA_DIR, cancer, 'PDA')
    pdalist = [join(pdadir, a) for a in getfilelist(pdadir) if a.split('_')[3] in TEST_CLINICAL_FEATURES]
    for pdapath in pdalist:
        if not pdapath.split('_')[5] == BASIC_OMICS or not (pdapath.split('_')[4] == 'Subtype'):#
            continue

        pda = getPDA(pdapath)
        group = pda.clinical_feature_groups
        min = group.min()

        test_sample_list = []
        for i in range(0, 11, 1):
            test_sample = []
            for g in group.values:
                w = (g - min) / 10
                s = floor(w * i + min)
                temp = g
                if temp > s:
                    temp = s
                test_sample.append(temp)
            # print(test_ratio)
            test_sample_list.append(test_sample)
            # ratio = [round(a/min, 2) for a in test_sample]
        for test_balance in test_sample_list:
            #thread(cancer, pdapath, test_balance)
            p = Process(target=thread, args=(cancer, pdapath, test_balance, False, False))
            proc.append(p)
            p = Process(target=thread, args=(cancer, pdapath, test_balance, False, True))
            proc.append(p)
            p = Process(target=thread, args=(cancer, pdapath, test_balance, True, True))
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


