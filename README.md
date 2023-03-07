# MOPARAM
A general guideline on parameter selection for multi-omics intergrative analysis.


## Table of Contents
1. [Introduction](#introduction)
2. [Getting start](#getting-start)
    1. [Prerequisites](#prerequisites)
    2. [Installing Dependencies](#installing-dependencies)
    3. [Installing the App](#installing-the-app)
3. [Usage](#usage)
    1. [Starting the App](#starting-the-app)
    2. [Using the App](#using-the-app)
4. [Contributing](#contributing)



## Introduction <a name="introduction"></a>
With the rapid development of high-throughput sequencing technologies, omics features can now be measured with greater detail and accurate biological context. To extract meaningful multi-omics relationships specific to a particular biological phenomena, various multi-omics integration (MOI) analysis methods have been proposed. While the biological and technical characteristics of each single-omics dataset have been studied in depth, many aspects of integrating two or more omics types still need to be investigated for efficient and robust results. In particular, the need to map multi-omics analysis results to clinical outcomes is increasing, with the expectation of an accurate explanation for biological problems. Aggregating different omics types results in a single dataset that is highly heterogeneous where the units of measurement, number of samples, and features per omics differ. One difficulty of analyzing such data is the absence of a generalized guideline, that entails the burden for making a number of important but difficult decisions, such as selection of samples, omics features. However, there is currently no standard guideline on how to perform multi-omics data analysis, including preprocessing and integration. While a guideline would be specific to a certain objective, a general guideline can be proposed by testing a wide range of different cohorts and objectives with multi-omics data. To provide a general guideline for multi-omics analysis, we identified nine crucial factors, including sample size, feature selection, preprocessing type, noise ratio, sample balance, number of groups, cancer subtype combination, omics combination, and clinical features. These were the factors that we considered to have impact on the performance of MOI analysis from both data and biological perspectives. Based on these factors, we devised seven benchmark tests and applied ten clustering-based multi-omics methods using various types of TCGA cancer datasets. The performance of the benchmark tests was evaluated with a range of different settings for each factor. The benchmark test results showed that a robust performance was achieved when using a minimum of 26 samples per group, selecting less than 10% of features, restricting the sample balance under a ratio of 3:1, and keeping the noise level below 30%. Among the factors, feature selection had the greatest impact on performance, resulting up to 34% enhancement in clustering performance with correct feature selection. Additionally, determining the optimal combination of omics and comparing clearly defined sample groups were important for improving the accuracy. We find that such general guidelines can be commonly applied to various practices on multi-omics data analysis.


## Getting start <a name="getting-start"></a>

### Prerequisites <a name="prerequisites"></a>

...

### Installing Dependencies <a name="installing-dependencies"></a>

...

### Installing the App <a name="installing-the-app"></a>

...

## Usage <a name="usage"></a>

### Starting the App <a name="starting-the-app"></a>

...

### Using the App <a name="using-the-app"></a>

...

## Contributing <a name="contributing"></a>

...

