from os.path import join

import pandas as pd
import pickle
import seaborn as sns, numpy as np
from matplotlib import pyplot as plt

from library.TEST_UTILS import TEST_CANCERS, DATA_DIR, TEST_CLINICAL_FEATURES
from library.modb_utils import getfilelist, getPDA

# fig = plt.figure(figsize=(18, 14))
# fig.subplots(nrows=5, ncols=2)

for cancer in TEST_CANCERS:
    #if cancer != 'COAD':
    #    continue
    print(cancer)
    pdadir = join(DATA_DIR, cancer, 'PDA')
    pdalist = [join(pdadir, a) for a in getfilelist(pdadir) if a.split('_')[3] in TEST_CLINICAL_FEATURES]
    for pdapath in pdalist:
        if not pdapath.split('_')[4] == 'Age':
            continue
        # if not pdapath.split('_')[-3] == 'GE-METH450-MIRNA':
        #     continue
        if pdapath.split('_')[-3] == 'GE':
            continue
        print(pdapath.split('_')[-3])
        #ax.set_title(cancer)

        pda = getPDA(pdapath)
        df = pda.clinical_feature_labels
        new_df = pd.Series(index=df.index)
        if 'Middle' in df.value_counts().index:
            continue

        xpos = [40, 65]
        new_df[df < xpos[0]] = 'Young'
        new_df[ (df >= xpos[0]) & (df < xpos[1])] = 'Middle'
        new_df[df >= xpos[1]] = 'Old'

        pda.clinical_feature_labels = new_df
        pda.clinical_feature_groups = new_df.value_counts()
        pda.clinical_feature_groups_n = len(new_df.value_counts())

        # sns.histplot(df, ax=ax, kde=True)
        pda_out = open(pdapath, 'wb')
        pickle.dump(pda, pda_out)
        pda_out.close()

        # for xc in xpos:
        #     ax.axvline(x=xc, color='r', linestyle='--')

# fig.tight_layout()
# plt.savefig('../Report2/Age.dist.png')