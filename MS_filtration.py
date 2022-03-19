#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2022/1/4 18:19
# @Author  : Yu Zhang

import os
import pandas as pd
import re
from neoantigen_utils import NeoantigenUtils


def run_command(command):
    for c in command:
        print(c)
        os.system(c)


def handle01_create_db():
    command = ["cat " + utils.temp_file + "/*fasta > " + utils.temp_file + "/file/ref+var+pep.fasta"]
    run_command(command)


def handle02_xml():
    command = [
        "python3 software/gen_mqpar.py software/templates/labelfree.xml input/mass -o " + utils.temp_file + "/mass.xml"]
    run_command(command)


def handle03_maxquant():
    command = [
        "mono software/MaxQuant/bin/MaxQuantCmd.exe " + utils.temp_file + "/mass.xml > " + utils.temp_file + "/mass.txt"]
    # 得到filter-mass/mass/combined/txt/peptides.txt
    run_command(command)


def handle04_maxpep():
    peptide = pd.read_table(utils.temp_file + "/mass/combined/txt/peptides.txt", sep='\t', header=0)
    peptide = peptide.dropna(subset=['Proteins', 'Start position', 'End position'])
    peptide = peptide.reset_index(drop=True)
    output_list = []
    for line in range(len(peptide)):
        pep = peptide.loc[line, 'Proteins'].split(';')
        # before = peptide.loc[line, 'Amino acid before']
        # after = peptide.loc[line, 'Amino acid after']
        seq = peptide.loc[line, 'Sequence']
        start = int(peptide.loc[line, 'Start position'])
        end = int(peptide.loc[line, 'End position'])
        for p in pep:
            s = p.split('|')
            if len(s) == 2:
                mut_list = s[1].split('.')
                mutation = mut_list[1]  # 'A337V'  'H296_E297insY'  'G2403delinsDR'  # delins删除并插入
                if 'delins' in mutation:
                    mut_bef = mutation[0]
                    mut_aft = mutation.split('delins')[1][0]
                    mut_pos = int(re.findall(re.compile(r'\d{1,8}'), mutation)[0])
                elif 'ins' in mutation:
                    mut_bef = mutation.split('_')[1][0]
                    mut_aft = mutation.split('ins')[1][0]  # 取insert的首位字符
                    mut_pos = int(mutation.split('_')[1].split('ins')[0][1:])
                elif 'fs' in mutation:
                    fs = mutation.split('fs')[0]
                    mut_bef = fs[0]
                    mut_aft = fs[-1]
                    mut_pos = int(mutation[1:-1])
                elif re.findall(re.compile(r'[A-Z]\d{1,8}[A-Z]'), mutation):
                    mut_bef = mutation[0]
                    mut_aft = mutation[-1]
                    mut_pos = int(mutation[1:-1])
                else:
                    mut_pos = -1
                if 0 <= mut_pos <= end and start <= mut_pos:
                    pos = mut_pos - start
                    if seq[pos] == mut_aft:
                        final = mut_bef + str(pos + 1) + mut_aft
                        output_list.append([seq, s[0], final, (pos + 1), s[1]])
    output = pd.DataFrame(output_list, columns=['Sequence', 'Gene', 'AA Change', 'Mut Position1', 'Mut'])
    output.to_csv(utils.temp_file + "/Maxpep.txt", index=False, sep='\t')


if __name__ == '__main__':
    utils = NeoantigenUtils("filter-mass")
    handle01_create_db()
    handle02_xml()
    handle03_maxquant()
    handle04_maxpep()
