import numpy as np
import pandas as pd
import scipy.stats.stats

###################################################################################

# filtering = "all_genes"
filtering = "protein_coding"
window_size = int(500)

###################################################################################

# Creating a dictionary of ENSG - gene symbol
with open("hgnc_complete_set.txt")as hgnc_complete_set:
    ensg_gene_symbol_dict = {}
    gene_symbol_ensg_dict = {}
    next(hgnc_complete_set)  # skip the first line (headlines)
    for line in hgnc_complete_set:
        ensg_name = ''
        for word in line.split():
            if word.startswith('ENSG'):
                ensg_name = word
        gene_symbol = line.split("\t")[1]
        ensg_gene_symbol_dict.update({ensg_name: gene_symbol})
        gene_symbol_ensg_dict.update({gene_symbol: ensg_name})
# print(ensg_gene_symbol_dict)  # test

###################################################################################

# creates a dictionary tissue_index_dict = {"tissue name: [list of sample indexes]"}
# creates a dictionary sample_tissue_dict = {"sample_name: tissue_name"}
tissue_index_dict = {}
tissue_samples_dict = {}  # maybe redundant
index = -1
with open("dr_1_sample_tissues_440.txt") as f:
    for sample_name in f:
        index += 1
        line = sample_name.split()
        tissue = ''.join(line[1:])
        sample = line[0]
        # print sample
        if tissue in tissue_index_dict:
            tissue_index_dict[tissue].append(index)
        else:
            tissue_index_dict[tissue] = [index]
        if tissue in tissue_samples_dict:
            tissue_samples_dict[tissue].append(sample)
        else:
            tissue_samples_dict[tissue] = []
            tissue_samples_dict[tissue].append(sample)

# print tissue_samples_dict  # test

######################################################################################################

# CREATE A DICT OF ALL GENES AS EMPTY DICTS
with open("dr1_samples_440_normalized_counts_cpm_edgeR_above_0.txt")as all_tissue_counts:
    gene_noise_dict = {}
    gene_percentile_dict = {}
    gene_expression_dict = {}
    gene_mean_expression_dict = {}
    gene_std_dict = {}
    gene_cv_dict = {}

    next(all_tissue_counts)  # skip the first line (samples names)
    for line in all_tissue_counts:
        gene_name = line.split()[0]
        gene_name = gene_name.split(".")[0]  # removing "." and numbers that follow in gene name
        gene_noise_dict.update({gene_name: {}})
        gene_percentile_dict.update({gene_name: {}})
        gene_expression_dict.update({gene_name: {}})
        gene_mean_expression_dict.update({gene_name: {}})
        gene_std_dict.update({gene_name: {}})
        gene_cv_dict.update({gene_name: {}})

########################################################################################################

# CREATE A LIST OF ALL GENES
gene_list = []
with open("dr1_samples_440_normalized_counts_cpm_edgeR_above_0.txt")as all_tissue_counts:
    next(all_tissue_counts)  # skip the first line (samples names)
    for line in all_tissue_counts:
        gene_name = line.split()[0]
        gene_name = gene_name.split(".")[0]  # removing "." and numbers that follow in gene name
        gene_list.append(gene_name)

########################################################################################################

# CREATE A LIST OF ALL TISSUES
tissues = []

########################################################################################################

# filter not protein-coding genes

protein_coding = []
with open("biomart_protein_coding.txt")as protein_coding_genes:
    next(protein_coding_genes)
    for gene in protein_coding_genes:
        protein_coding.append(gene.split(",")[0])
        # print gene.split(",")[0]
intersect_list = []
intersect_list = [val for val in gene_list if val in protein_coding]


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

if filtering == "protein_coding":
    gene_list = intersect_list  # just protein-coding.

########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################

# Initiating correlation files
# write headlines + 'w' deleting data from previous run in output file

# cv-counts
target4 = open("cv-median_counts-corr_%s.txt" % filtering, 'w')
target4.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target4.close()

target14 = open("cv-mean_counts-corr_%s.txt" % filtering, 'w')
target14.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target14.close()


# percentiles-counts
target8 = open("percentile-median_counts-corr_%s.txt" % filtering, 'w')
target8.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target8.close()

target18 = open("percentile-mean_counts-corr_%s.txt" % filtering, 'w')
target18.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target18.close()


# std-counts
target12 = open("std-median_counts-corr_%s.txt" % filtering, 'w')
target12.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target12.close()

target22 = open("std-mean_counts-corr_%s.txt" % filtering, 'w')
target22.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target22.close()


# noise-percentiles
target10 = open("cv-percentile-corr_%s.txt" % filtering, 'w')
target10.write("Tissue" + "\t" + "Spearman_corr_(rho)" + "\t" + "p-value" + "\n")
target10.close()


#######################################################################################

# MAIN LOOP: FOR EVERY TISSUE - CALCULATES GENE NOISE FOR ALL GENES THAT PASSED THE CUTOFF

for tissue in tissue_index_dict:
    first = min(tissue_index_dict[tissue])
    last = max(tissue_index_dict[tissue])
    # print first  # test
    # print last  # test
    if ((last - first) < 9) or (tissue == 'Cells-Transformedfibroblasts'): # only for tissues with more than 10 samples
        continue
    tissues.append(tissue)
    with open("dr1_samples_440_normalized_counts_cpm_edgeR_above_0.txt")as all_tissue_counts:
        tissue_counts_dict = {}
        next(all_tissue_counts)  # skip the first line (samples names)
        for gene in all_tissue_counts:
            line = gene.split()  # split() turns a line of text into a list of words
            gene_name = line[0].split(".")[0]  # removing "." and numbers that follow in gene name

            #######################################################################################
            #######################################################################################
            #######################################################################################
            #######################################################################################

            # if there are more than 20% '0' count/non-protein-coding - the gene will be filtered out
            # if (sum(1.0 for v in line[first+1:last+2] if v != '0') / len(line[first+1:last+2])) >= 0.8:
            if (sum(1.0 for v in line[first+1:last+2] if float(v) > 7.0) / len(line[first+1:last+2])) >= 0.8:

                #######################################################################################
                #######################################################################################
                #######################################################################################
                #######################################################################################
                if gene_name in gene_list:
                    tissue_counts_dict[gene_name] = line[first+1:last+2]  # dict: {gene_name: [count1, count2, count3, ...]}
        # print tissue_counts_dict   #test

        # turning string counts to float counts in "tissue_counts_dict"
        for gene in tissue_counts_dict:
            new_array = []
            for count in tissue_counts_dict[gene]:
                new_array.append(float(count))
            tissue_counts_dict[gene] = new_array
        # print tissue_counts_dict   #test

        '''
        # removing outliers:
        for gene in tissue_counts_dict:
            gene = sorted(gene)
            n = len(gene)
            # print math.ceil(0.1*n)
            # print n-math.ceil(0.1*n)
            gene = gene[int(math.ceil(0.1*n)): int((n-math.ceil(0.1*n)))]
        '''

        # calculating std for counts of each gene
        gene_std = {}
        for gene in tissue_counts_dict:
            gene_std[gene] = np.nanstd(tissue_counts_dict[gene])
        # print gene_std   # test

        # calculating mean for counts of each gene
        gene_mean = {}
        for gene in tissue_counts_dict:
            gene_mean[gene] = np.nanmean(tissue_counts_dict[gene])
        # print gene_mean   # test

        # calculating median for counts of each gene
        gene_median = {}
        for gene in tissue_counts_dict:
            gene_median[gene] = np.nanmedian(tissue_counts_dict[gene])
        # print gene_median   # test

        # calculating 75 and 25 percentiles, and iqr
        gene_q75 = {}
        gene_q25 = {}
        gene_iqr = {}
        for gene in tissue_counts_dict:
            gene_q75[gene] = np.percentile(tissue_counts_dict[gene], 75)
            gene_q25[gene] = np.percentile(tissue_counts_dict[gene], 25)
            gene_iqr[gene] = gene_q75[gene] - gene_q25[gene]


        # calculating noise(cv=std/mean)for counts of each gene
        gene_cv = {}
        for gene in tissue_counts_dict:
            if gene_mean[gene] != 0:  # only if the mean != 0
                gene_cv[gene] = (gene_std[gene]) / (gene_mean[gene])
        # print gene_noise  # test

        # calculating noise2(std/median)for counts of each gene
        gene_mcv = {}
        for gene in tissue_counts_dict:
            if gene_median[gene] != 0:  # only if the mean != 0
                gene_mcv[gene] = (gene_std[gene]) / (gene_median[gene])
        # print gene_noise  # test

        # calculating noise3(Quartile coefficient of dispersion)for counts of each gene
        gene_cqv = {}
        for gene in tissue_counts_dict:
            if (gene_q75[gene] + gene_q25[gene]) != 0:
                gene_cqv[gene] = (gene_q75[gene] - gene_q25[gene]) / (gene_q75[gene] + gene_q25[gene])
        # print gene_noise  # test

        # calculating noise4(IQR/median)for counts of each gene
        gene_iqr_median = {}
        for gene in tissue_counts_dict:
            if gene_median[gene] != 0:  # only if the mean != 0
                gene_iqr[gene] = (gene_q75[gene] - gene_q25[gene]) / (gene_median[gene])
        # print gene_noise  # test

        ########################################################################################################
        ########################################################################################################
        ########################################################################################################
        ########################################################################################################

        gene_noise = gene_cv

        ########################################################################################################
        ########################################################################################################
        ########################################################################################################
        ########################################################################################################

        # updating the gene`s entry in gene_noise_dict with it`s noise in this tissue
        for gene in gene_list:
            if gene in gene_noise:
                gene_noise_dict[gene].update({tissue: gene_noise[gene]})
            else:
                gene_noise_dict[gene].update({tissue: np.nan})

        # updating the gene`s entry in gene_std_dict with it`s noise in this tissue
        for gene in gene_list:
            if gene in gene_noise:
                gene_std_dict[gene].update({tissue: gene_std[gene]})
            else:
                gene_std_dict[gene].update({tissue: np.nan})

        # updating the gene`s entry in gene_cv_dict with it`s noise in this tissue
        for gene in gene_list:
            if gene in gene_cv:
                gene_cv_dict[gene].update({tissue: gene_cv[gene]})
            else:
                gene_cv_dict[gene].update({tissue: np.nan})
########################################################################################################
        # per-tissue gene expression matrix with all tissues

        target_samples = open('df_expression_samples_%s.csv' % tissue, 'w')
        target_samples.close()
        target_samples = open('df_expression_samples_%s.csv' % tissue, 'a+')
        for gene in tissue_counts_dict:
            target_samples.write(gene + ",")
            for count in tissue_counts_dict[gene]:
                target_samples.write(str(count) + ',')
            target_samples.write(',' + "\n")
        target_samples.close()


########################################################################################################

    # initializing output files
    target9 = open("sorted_gene_percentile_%s_%s.txt" % (tissue, filtering), 'w')
    target9.close()
    target5 = open("gene_percentile_noise_%s_%s.txt" % (tissue, filtering), 'w')
    target5.close()
    target11 = open("sorted_gene_percentile_%s_%s_reverse.txt" % (tissue, filtering), 'w')
    target11.close()

    # Sorting genes by noise values

    # sorting gene_noise key by values
    Sorted_Dict_Values = sorted(gene_noise.values(), reverse=True)
    Sorted_Dict_Keys = sorted(gene_noise, key=gene_noise.get, reverse=True)

    # reverse sorting gene_noise key by values
    Sorted_Dict_Values_reverse = sorted(gene_noise.values(), reverse=False)
    Sorted_Dict_Keys_reverse = sorted(gene_noise, key=gene_noise.get, reverse=False)

    # sorting by median count values:
    sorted_by_median_values = sorted(gene_median.values(), reverse=True)
    sorted_by_median_keys = sorted(gene_median, key=gene_median.get, reverse=True)


#######################################################################################

    # Sliding window

    gene_percentile = {}
    for i in range(0, len(sorted_by_median_values)):
        gene_name = sorted_by_median_keys[i]
        window_keys_list = []
        if i < (window_size / 2):
            window_keys_list = sorted_by_median_keys[0:window_size]
        elif (window_size / 2) <= i < (len(sorted_by_median_values) - (window_size / 2)):
            window_keys_list = sorted_by_median_keys[(i - (window_size / 2)):(i + (window_size / 2))]
        else:
            window_keys_list = sorted_by_median_keys[(len(sorted_by_median_values) - window_size):]

        window_noise_dict = {}
        for key in window_keys_list:
            window_noise_dict[key] = gene_noise[key]

        # sort window keys by noise values
        # bottom up
        sorted_window_keys = sorted(window_noise_dict, key=window_noise_dict.get, reverse=False)
        gene_percentile[gene_name] = ((sorted_window_keys.index(gene_name)+1.0)/window_size) * 100.0

    # updating the gene`s entry in gene_percentile_dict with it`s noise percentile in this tissue
    for gene in gene_list:
        if gene in gene_percentile:
            gene_percentile_dict[gene].update({tissue: gene_percentile[gene]})
            gene_expression_dict[gene].update({tissue: gene_median[gene]})
            gene_mean_expression_dict[gene].update({tissue: gene_mean[gene]})
        else:
            gene_percentile_dict[gene].update({tissue: np.nan})
            gene_expression_dict[gene].update({tissue: np.nan})
            gene_mean_expression_dict[gene].update({tissue: np.nan})


#######################################################################################

    # sorting gene according to percentiles

    sorted_gene_percentile = sorted(gene_percentile.values(), reverse=True)
    sorted_gene_percentile_keys = sorted(gene_percentile, key=gene_percentile.get, reverse=True)
    sorted_gene_percentile_keys_reverse = sorted(gene_percentile, key=gene_percentile.get, reverse=False)


#######################################################################################

    # gene_name-noise_percentile output files

    for i in range(len(sorted_gene_percentile_keys)):
        target9 = open("sorted_gene_percentile_%s_%s.txt" % (tissue, filtering), 'a+')
        target9.write(sorted_gene_percentile_keys[i])
        # target9.write(ensg_gene_symbol_dict[sorted_gene_percentile_keys[i]])
        target9.write("\n")
        target9.close()
        # if sorted_gene_percentile_keys[i] not in ensg_gene_symbol_dict:
        #     continue
        target5 = open("gene_percentile_noise_%s_%s.txt" % (tissue, filtering), 'a+')
        target5.write(sorted_gene_percentile_keys[i])
        # target5.write(ensg_gene_symbol_dict[sorted_gene_percentile_keys[i]])
        target5.write("\t")
        target5.write(str(sorted_gene_percentile[i]))
        target5.write("\n")
        target5.close()
        # just gene names

    # reverse sorted lists
    for i in range(len(sorted_gene_percentile_keys_reverse)):
        # just gene names
        target11 = open("sorted_gene_percentile_%s_%s_reverse.txt" % (tissue, filtering), 'a+')
        target11.write(sorted_gene_percentile_keys_reverse[i])
        target11.write("\n")
        target11.close()

#######################################################################################

    # correlation test to see if there is correlation between median_counts and noise percentile
    percentile_corr = []
    count_corr = []
    for gene in gene_list:
        if gene in gene_percentile:
            percentile_corr.append(gene_percentile[gene])
            count_corr.append(gene_median[gene])
    rho, pval = scipy.stats.spearmanr(percentile_corr, count_corr)
    target8 = open("percentile-median_counts-corr_%s.txt" % filtering, 'a+')
    target8.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target8.close()


    # correlation test to see if there is correlation between mean_counts and noise percentile
    percentile_corr = []
    count_corr = []
    for gene in gene_list:
        if gene in gene_percentile:
            percentile_corr.append(gene_percentile[gene])
            count_corr.append(gene_mean[gene])
    rho, pval = scipy.stats.spearmanr(percentile_corr, count_corr)
    target18 = open("percentile-mean_counts-corr_%s.txt" % filtering, 'a+')
    target18.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target18.close()


    # correlation test to see if there is correlation between median_counts and cv
    cv_corr = []
    count_corr = []
    for gene in gene_list:
        if gene in gene_cv:
            cv_corr.append(gene_cv[gene])
            count_corr.append(gene_median[gene])
    rho, pval = scipy.stats.spearmanr(cv_corr, count_corr)
    target4 = open("cv-median_counts-corr_%s.txt" % filtering, 'a+')
    target4.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target4.close()


    # correlation test to see if there is correlation between mean_counts and cv
    cv_corr = []
    count_corr = []
    for gene in gene_list:
        if gene in gene_cv:
            cv_corr.append(gene_cv[gene])
            count_corr.append(gene_mean[gene])
    rho, pval = scipy.stats.spearmanr(cv_corr, count_corr)
    target14 = open("cv-mean_counts-corr_%s.txt" % filtering, 'a+')
    target14.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target14.close()


    # correlation test to see if there is correlation between median_counts and std
    std_corr = []
    count_corr = []
    for gene in gene_list:
        if gene in gene_std:
            std_corr.append(gene_std[gene])
            count_corr.append(gene_median[gene])
    rho, pval = scipy.stats.spearmanr(std_corr, count_corr)
    target12 = open("std-median_counts-corr_%s.txt" % filtering, 'a+')
    target12.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target12.close()


    # correlation test to see if there is correlation between median_counts and std
    std_corr = []
    count_corr = []
    for gene in gene_list:
        if gene in gene_std:
            std_corr.append(gene_std[gene])
            count_corr.append(gene_mean[gene])
    rho, pval = scipy.stats.spearmanr(std_corr, count_corr)
    target22 = open("std-mean_counts-corr_%s.txt" % filtering, 'a+')
    target22.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target22.close()


    # correlation test to see if there is correlation between std/mean(cv)-based and percentile-based lists
    percentile_corr = []
    cv_corr = []
    for gene in gene_list:
        if gene in gene_percentile:
            percentile_corr.append(gene_percentile[gene])
            cv_corr.append(gene_cv[gene])
    rho, pval = scipy.stats.spearmanr(percentile_corr, cv_corr)
    target10 = open("cv-percentile-corr_%s.txt" % filtering, 'a+')
    target10.write(tissue + "\t" + str(rho) + "\t" + str(pval) + "\n")
    target10.close()

####################################################################################

    # creating a graph of gene noise against expression level

    # sorting by count values:
    sorted_by_median_values = sorted(gene_median.values(), reverse=True)
    sorted_by_median_keys = sorted(gene_median, key=gene_median.get, reverse=True)

    count = sorted_by_median_values
    noise = []
    log_noise_dict = {}
    for key in sorted_by_median_keys:
        noise.append(gene_noise[key])
        log_noise_dict[key] = np.log(gene_noise[key])

########################################################################################

#  HERE WE CREATE OUR PANDAS DATA FRAME

gene_noise_list = []
for gene in gene_list:
    gene_noise_list.append(gene_percentile_dict[gene])

df = pd.DataFrame.from_records(gene_noise_list, columns=tissues, index=gene_list)


# remove genes with no noise data

df_all = df.dropna(how='all')   # drop only if ALL columns are NaN
df_all.to_csv('gene_noise_%s_%s.csv' % (window_size, filtering))

df_no_nan = df.dropna()     # drop all rows that have any NaN values
# df_no_nan = df.dropna(thresh=10)   # Drop row if it does not have at least # values that are **not** NaN
# print df_no_nan  # test
# del df_no_nan['Cells-Transformedfibroblasts']
df_no_nan.to_csv('gene_noise_no_nan_%s_%s.csv' % (window_size, filtering))


#####################################################################

gene_expression_list = []
for gene in gene_list:
    gene_expression_list.append(gene_expression_dict[gene])

df_expression = pd.DataFrame.from_records(gene_expression_list, columns=tissues, index=gene_list)

# remove genes with no noise data
df_expression_all = df_expression.dropna(how='all')   # drop only if ALL columns are NaN
df_expression_all.to_csv('gene_expression_%s.csv' % filtering)



gene_mean_expression_list = []
for gene in gene_list:
    gene_mean_expression_list.append(gene_mean_expression_dict[gene])

df_mean_expression = pd.DataFrame.from_records(gene_mean_expression_list, columns=tissues, index=gene_list)

# remove genes with no noise data
df_mean_expression_all = df_mean_expression.dropna(how='all')   # drop only if ALL columns are NaN
df_mean_expression_all.to_csv('gene_mean_expression_%s.csv' % filtering)




gene_std_list = []
for gene in gene_list:
    gene_std_list.append(gene_std_dict[gene])

df_std = pd.DataFrame.from_records(gene_std_list, columns=tissues, index=gene_list)

# remove genes with no noise data
df_std_all = df_std.dropna(how='all')   # drop only if ALL columns are NaN
df_std_all.to_csv('gene_std_%s_%s.csv' % (window_size, filtering))



gene_cv_list = []
for gene in gene_list:
    gene_cv_list.append(gene_cv_dict[gene])

df_cv = pd.DataFrame.from_records(gene_cv_list, columns=tissues, index=gene_list)

# remove genes with no noise data
df_cv_all = df_cv.dropna(how='all')   # drop only if ALL columns are NaN
df_cv_all.to_csv('gene_cv_%s_%s.csv' % (window_size, filtering))
