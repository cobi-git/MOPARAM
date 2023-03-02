from os.path import join
import seaborn as sns, numpy as np
from matplotlib import pyplot as plt

from library.TEST_UTILS import TEST_CANCERS, DATA_DIR, TEST_CLINICAL_FEATURES
from library.modb_utils import getfilelist, getPDA



fig = plt.figure(figsize=(18, 14))
fig.subplots(nrows=5, ncols=2)
for cancer, ax in zip(TEST_CANCERS, fig.axes):
    #if cancer != 'COAD':
    #    continue
    print(cancer)
    pdadir = join(DATA_DIR, cancer, 'PDA')
    pdalist = [join(pdadir, a) for a in getfilelist(pdadir) if a.split('_')[3] in TEST_CLINICAL_FEATURES]
    for pdapath in pdalist:
        if not pdapath.split('_')[4] == 'pathologic-stage':
            continue
        if not pdapath.split('_')[-3] == 'GE-ME-MI':
            continue
        ax.set_title(cancer)
        pda = getPDA(pdapath)
        age_df = pda.clinical_feature_labels
        sns.histplot(age_df, ax=ax,  kde=True)
fig.tight_layout()
plt.savefig('../Report2/Stage.dist.png')