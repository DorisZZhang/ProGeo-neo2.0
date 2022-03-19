#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2022/1/5 10:18
# @Author  : Yu Zhang

import pandas as pd
import os
from neoantigen_utils import NeoantigenUtils


def find_pep(pep, sequence):
    for seq in sequence:
        if pep in seq:
            return True
    return False


def handle01_tpm_filter(mut_type, hla_class=1):
    file_name = '/MHC-I' if hla_class == 1 else '/MHC-II'
    input1 = pd.read_table("reference_file/Tran_id_gene.txt", names=['gene', 'EM'], sep=',')
    input1.dropna(axis=0, how='any', inplace=True)
    em_dic = dict(zip(input1['EM'], input1['gene']))
    input2 = pd.read_table(utils.out_file + "outputs/abundance.tsv", sep='\t')
    input2 = input2[input2['tpm'] > 0]  # 将abundance中tpm>0的留下
    gene_list = []  # 保留dict中存在的gene及其对应tpm
    for i in input2.index:
        target_id = input2.at[i, 'target_id']
        if target_id in em_dic:
            gene_list.append([em_dic[target_id], input2.at[i, 'tpm']])
    gene_df = pd.DataFrame(gene_list, columns=['gene', 'tpm'])
    gene_df.drop_duplicates(inplace=True)
    gene_df.to_csv(utils.temp_file + '/kallisto_gene_expression.txt', sep='\t', index=False)  # 保存输出文件1为txt
    tpm_dict = dict(gene_list)
    for mut in mut_type:
        input3 = pd.read_table('reference_file/neoantigen_' + mut + '_1.txt', sep='\t', dtype=str)
        if mut_type == 'fusion':  # 融合基因：若有一个tpm>0则保留
            input3[['gene1', 'gene2']] = input3['gene'].str.split('--', expand=True)
            input3['tpm1'] = input3['gene1'].apply(lambda g: tpm_dict[g] if g in tpm_dict else None)
            input3['tpm2'] = input3['gene2'].apply(lambda g: tpm_dict[g] if g in tpm_dict else None)
            input3 = input3[~input3[['tpm1', 'tpm2']].isnull().T.all()]  # 保留tpm不全为空的行
            input3.drop('gene', axis=1, inplace=True)
        else:
            input3['tpm'] = input3['gene'].apply(lambda g: tpm_dict[g] if g in tpm_dict else None)
            input3.dropna(axis=0)
        input3.to_csv(utils.out_file + '/' + file_name +'/tpm_filter1-neoantigen_' + mut + '_1.txt', sep='\t', index=False)  # 保存输出文件2为txt


def handle02_MS_filter():
    input1 = pd.read_table("zq1227_tpm_snv_1.txt", sep='\t')
    input2 = pd.read_table("filter-mass/Maxpep-fdr0.01.txt", sep='\t')
    find_result = []
    for pep in input1['Peptide']:
        find_result.append(find_pep(pep, input2['Sequence'].values.tolist()))
    input1 = input1[find_result].copy()
    input1.drop_duplicates(inplace=True)
    input1.to_csv("outfile2/candid1_max_filter2-neoantigen.txt", index=False, sep='\t')


def handle03_threshold_filter():
    mut_type = 'fs'
    input1 = pd.read_csv('outfile2/tpm_filter1-neoantigen_' + mut_type + '_1.txt', sep='\t')
    # Aff(nM)≤34且tpm≥33
    input1 = input1[input1.iloc[:, 5] <= 34]
    input1 = input1[input1.iloc[:, -1] >= 33]
    input1.to_csv('outfiles/HI-Q_filter-neoantigen_' + mut_type + '_1.txt', sep='\t', index=False)


if __name__ == '__main__':
    utils = NeoantigenUtils(None, "outfile-filter-neoantigens")
    mut_type = ['snv', 'del', 'ins', 'fs', 'fusion']
    if not os.path.exists(utils.out_file):
        os.mkdir(utils.out_file)
        dir_list = ['/MHC-I', '/MHC-II']
        for d in dir_list:
            for mut in mut_type:
                os.mkdir(utils.out_file + d + '/' + mut)
    handle01_tpm_filter(mut_type)
    handle02_MS_filter()
    handle03_threshold_filter()
