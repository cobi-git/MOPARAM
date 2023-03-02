from library.TEST_UTILS import TEST_CANCERS, TEST_CLINICAL_FEATURES, core_n, check_exist
from library.clustering_test import *

testname = 'omics-combination-test'

proc = []
proc_excution = []
proc_end = []

def thread(cancer, pdapath, arg,gene_centric, norm):
    testcase = '%s'%(arg)
    data_raw, labels, csvnames, filename, labelpath = \
        clustering_preprocessing(cancer, pdapath, testname, testcase, gene_centric_bool=gene_centric, norm_bool=norm)
    clustering_test(data_raw, labels, csvnames, filename, labelpath)

for cancer in TEST_CANCERS:
    pdadir = join(DATA_DIR, cancer, 'PDA')
    pdalist = [join(pdadir, a) for a in getfilelist(pdadir) if a.split('_')[3] in TEST_CLINICAL_FEATURES]
    for pdapath in pdalist:
        if not pdapath.split('_')[4] == 'Subtype':
            continue
        if pdapath.endswith('error.pda'):
            continue
        omics_string = pdapath.split('_')[5]
        if omics_string == 'GE':
            continue
        #thread(cancer, pdapath, omics_string, False, True)
        p = Process(target=thread, args=(cancer, pdapath, omics_string, False, False))
        proc.append(p)
        p = Process(target=thread, args=(cancer, pdapath, omics_string, False, True))
        proc.append(p)
        p = Process(target=thread, args=(cancer, pdapath, omics_string, True, True))
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