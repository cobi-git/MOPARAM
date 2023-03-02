from library.TEST_UTILS import TEST_CLINICAL_FEATURES
from library.modb_utils import *
import pandas as pd

cancerlist = sorted(getdirlist(PROJECTDIR))
clinicalset = TEST_CLINICAL_FEATURES

def viewPDA():
    for cancer in cancerlist:
        print("**********"+cancer+"**********")
        dir = '/'.join([PROJECTDIR, cancer, 'PDA'])
        pdalist = getfilelist(dir)
        for pdaname in pdalist:
            if pdaname.endswith('error.pda'):
                continue
            clincal = pdaname.split('_')[3]
            if clincal == 'Age':
                continue
            print('-------' + clincal + '-------')
            pda = getPDA(join(dir, pdaname))
            print(pda.clinical_feature_groups)

columns = [a for a in TEST_CLINICAL_FEATURES]
columns.insert(0, 'Cancer')
print(columns)

omics_columns = [a for a in PERMITTED_OMICS]
omics_columns.insert(0, 'Cancer')
print(omics_columns)
total = []
omics = []
for cancer in cancerlist:
    row = []
    row.append(cancer)
    omics_info = []
    omics_info.append(cancer)
    path = '/'.join([PROJECTDIR, cancer, 'PDA'])
    filelist = getfilelist(path)
    print("**********" + cancer + "**********")
    # if cancer != 'KIRP':
    #     continue
    for c in TEST_CLINICAL_FEATURES:
        print(c)
        file = [filename for filename in filelist if c in filename.split('_')[3]
                and not 'error' in filename and '-'.join(PERMITTED_OMICS) == filename.split('_')[4]]
        if len(file) != 0:
            filename = file[0]
            f = open(join(path, filename), 'rb')
            pda = pickle.load(f)
            f.close()
            if c == 'Subtype':
                d = pda.data_omics_level
                for a in d:
                    omics_info.append(str(d[a].shape[0]))
                omics.append(omics_info)
            input = pda.clinical_feature_groups
            input = input[input >= 10]
            input_sum =input.sum()
            input = '(' + input.astype(str) + ')'
            input = input.to_string().replace(' ','')
            input = input + '\n'+str(input_sum)
            print(input)
            row.append(input)
        else:
            row.append(None)
    total.append(row)
df_clinical = pd.DataFrame(data=total, columns=columns)
df_clinical.set_index('Cancer', inplace=True)

df_omics = pd.DataFrame(data=omics, columns=omics_columns)
df_omics.set_index('Cancer', inplace=True)

df_total = pd.concat([df_omics, df_clinical], axis=1)
df_total.to_csv('../Report2/TOTAL_Clinical_OMICS_INFO.csv')

#df_clinical.to_csv('../Report2/Clinical_INFO.csv')
#df_omics.to_csv('../Report2/OMICS_INFO.csv')

