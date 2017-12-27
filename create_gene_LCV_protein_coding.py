import numpy as np
import pandas as pd
import scipy.stats.stats

###################################################################################

# filtering = "all_genes"
filtering = "protein_coding"
window_size = int(500)

###################################################################################

# creates a dictionary tissue_index_dict = {"tissue name: [list of sample indexes]"}
# creates a dictionary sample_tissue_dict = {"sample_name: tissue_name"}
tissue_index_dict = {}
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

if filtering == "protein_coding":
    gene_list = intersect_list  # just protein-coding.

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

            if (sum(1.0 for v in line[first+1:last+2] if float(v) > 7.0) / len(line[first+1:last+2])) >= 0.8:

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

        # calculating noise(cv=std/mean)for counts of each gene
        gene_cv = {}
        for gene in tissue_counts_dict:
            if gene_mean[gene] != 0:  # only if the mean != 0
                gene_cv[gene] = (gene_std[gene]) / (gene_mean[gene])
        # print gene_noise  # test

        gene_noise = gene_cv

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
        else:
            gene_percentile_dict[gene].update({tissue: np.nan})

#######################################################################################

#  HERE WE CREATE OUR PANDAS DATA FRAME

gene_noise_list = []
for gene in gene_list:
    gene_noise_list.append(gene_percentile_dict[gene])

df = pd.DataFrame.from_records(gene_noise_list, columns=tissues, index=gene_list)

# remove genes with no noise data

df_all = df.dropna(how='all')   # drop only if ALL columns are NaN
df_all.to_csv('gene_noise_github%s_%s.csv' % (window_size, filtering))