# MOPARAM
A general guideline on parameter selection for multi-omics intergrative analysis.


## Table of Contents
1. [Introduction](#introduction)
2. [Getting start](#getting-start)
    1. [Prerequisites](#prerequisites)
    2. [Recommended system spec](#recomended-spec)
    3. [Data generation](#data-generation)
    4. [Preprocessing](#preprocessing)
3. [Benchmarking test](#benchmarking-test)
    1. [Execute all test in once](#exec-all-test)
    2. [Sample size test](#sample-size-test)
    3. [Feature selection test](#feature-selection-test)
    4. [Balance test](#balance-test)
    5. [Robustness test](#robustness-test)
    6. [Subtype combination test](#subtype-combination-test)
    7. [Omics combination test](#omics-combination-test) 
4. [Contributing](#contributing)


## Introduction <a name="introduction"></a>
With the advent of high-throughput sequencing technologies, multi-omics datasets are becoming increasingly complex, making it difficult to extract meaningful relationships between different biological features. In this paper, we propose a general guideline for multi-omics data analysis based on nine crucial factors, including sample size, feature selection, preprocessing type, noise ratio, sample balance, number of groups, cancer subtype combination, omics combination, and clinical features. We evaluate the effectiveness of the proposed guidelines by conducting benchmark tests on various types of TCGA cancer datasets using ten clustering-based multi-omics methods. Our results demonstrate that our proposed guidelines can significantly enhance the accuracy of multi-omics analysis, especially when it comes to feature selection, sample balance, and noise reduction. Overall, this paper provides a valuable resource for researchers and practitioners in the field, helping them to better navigate the challenges of multi-omics data analysis and extract meaningful relationships between different biological features.


## Getting start <a name="getting-start"></a>

### Prerequisites <a name="prerequisites"></a>
The following link is the data used in this study.

TCGA DATA: http://cobi.knu.ac.kr/tools/moparam/TCGA_DATA.tar.gz

*All compressed files must be extracted before use.

### Recommended system spec <a name="recomended-spec"></a>

Basically, the functions provided by the MOParam project are mainly data generation, preprocessing, and benchmarking tests, which operate in a parallel processing manner to handle large-scale data. To complete a single test within a day, it is recommended to perform them on a server-class computer with around 100 cores and 300GB of memory.

### Data generation <a name="data-generation"></a>
To use data generated by oneself in addition to the 10 TCGA cancer datasets refined in this study, the following data should be needed:

1. Clinical data
2. Omics data
<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/data_structure.jpg" alt="Data structure" width="800">

If the necessary files are prepared, it is recommended to store them in the following directory structure.

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/file_tree.png" alt="File tree" width="250">

Once everything is ready, you now need to integrate each omics data per sample.
This process is implemented through '/preprocessing/make_pda.py'.


### Preprocessing <a name="preprocessing"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/preprocessing.jpg" alt="Preprocessing" width="800">

1. RAW - The RAW dataset consists of omics data that has not undergone any preprocessing but were $log_2$ transformed.
2. NORM - The NORM dataset has undergone quantile normalization for each RAW omics data, as well as 0-1 scaling.
3. GL(Gene level) - The GL dataset utilizes normalized data that has undergone multi-stage integration to transform the features of the entire omics dataset into those of GE(Gene expression).

*Please refer to the /preprocessing/Gene_level_preprocessing.py code for information on generating the GL dataset.

## Benchmarking test <a name="benchmarking-test"></a>

###  Execute all test in once <a name="exec-all-test"></a>

Here is the code for exec_all_test.py. If you execute this code, you can perform all tests sequentially.

```python
from os.path import join
from library.TEST_UTILS import getfilelist
for file in sorted(getfilelist("./test")):
    file_path = join('./test', file)
    with open(file_path) as f:
        code_to_run = f.read()
    exec(code_to_run)
    


###  Sample size test <a name="sample-size-test"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/sample_size_test.jpg" alt="" width="700">

###  Feature selection test <a name="feature-selection-test"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/feature_selection_test.jpg" alt="" width="700">

###  Balance test <a name="balance-test"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/balance_test.jpg" alt="" width="700">

###  Robustness test <a name="robustness-test"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/robustness_test.jpg" alt="" width="700">

###  Subtype combination test <a name="subtype-combination-test"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/subtype_combination_test.jpg" alt="" width="700">

###  Omics combination test <a name="omics-combination-test"></a>

<img src="https://github.com/cobi-git/MOPARAM/blob/main/Images/omics_combination_test.jpg" alt="" width="700">

## Contributing <a name="contributing"></a>

Hwijun Kwon and Inuk Jung conceived the study and designed the experiments. Hwijun Kwon developed the multi-omics integration methods and performed the data analysis. Hwijun Kwon and Inuk Jung contributed to the experimental design and interpretation of the results. Hwijun Kwon and Inuk Jung wrote the manuscript. IlKon Kim provided critical feedback on the manuscript and helped with the interpretation of the results. All authors reviewed and approved the final version of the manuscript.

