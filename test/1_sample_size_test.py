from library.TEST_UTILS import TEST_CANCERS, BASIC_OMICS, TEST_CLINICAL_FEATURES, core_n, check_exist
from library.clustering_test import *
testname = 'sample-size-test'

proc = []
proc_excution = []
proc_end = []

def thread(testname, cancer, pdapath, arg, gene_centric, norm):
    testcase = '%d' % (arg)
    data_raw, labels, csvnames, filename, labelpath = \
        clustering_preprocessing(cancer, pdapath, testname, testcase, sample_size=arg, gene_centric_bool=gene_centric, norm_bool=norm)
    clustering_test(data_raw, labels, csvnames, filename, labelpath)

for cancer in TEST_CANCERS:
    pdadir = join(DATA_DIR, cancer, 'PDA')
    pdalist = [join(pdadir, a) for a in getfilelist(pdadir) if a.split('_')[3] in TEST_CLINICAL_FEATURES]
    for pdapath in pdalist:
        if not pdapath.split('_')[5] == BASIC_OMICS or not (pdapath.split('_')[4] == 'Subtype'):
            continue

        pda = getPDA(pdapath)
        min = pda.clinical_feature_groups.min()
        max_sample_size = 40

        if max_sample_size > min:
            max_sample_size = min

        start = 5
        if start > min:
            start = min
        for sample_size in range(start, max_sample_size, 3):
            #thread(testname, cancer, pdapath, sample_size, False, True)
            p = Process(target=thread, args=(testname, cancer, pdapath, sample_size, False, False))  # RAW
            proc.append(p)
            p = Process(target=thread, args=(testname, cancer, pdapath, sample_size, False, True))  # NORM
            proc.append(p)
            p = Process(target=thread, args=(testname, cancer, pdapath, sample_size, True, True))  # GC
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