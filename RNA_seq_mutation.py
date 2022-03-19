#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2021/12/24 16:01
# @Author  : Yu Zhang
"""
RNA data processing(tumor RNA-seq):
"""

import os
import pandas as pd
from neoantigen_utils import NeoantigenUtils


def run_command(command):
    for c in command:
        print(c)
        os.system(c)


def activate_env():
    os.system("source ~/.bashrc")


def handle01_bwa():
    if len(tumor_files) == 2:
        command = [r"bwa index reference_file/hg38/hg38.fa",
                   r"bwa mem -t 8 -M -q -5 -R '@RG\tID:foo\tSM:bar\tLB:library1' reference_file/hg38/hg38.fa input/rna/" +
                   tumor_files[0] + " input/rna/" + tumor_files[1] + " > " + utils.temp_file + "/sample.sam"]
        run_command(command)
    else:
        print("please input files into the right position!")


def handle02_sam():
    """
    Map reads to reference genome  BWA（sam）
    """
    command = ["samtools fixmate -O bam " + utils.temp_file + "/sample.sam " + utils.temp_file + "/sample.bam",
               "samtools sort -O bam -o " + utils.temp_file + "/sorted_sample.bam -T temp " + utils.temp_file + "/sample.bam"]
    run_command(command)


def handle03_gatk():
    command = [
        "gatk MarkDuplicates -I " + utils.temp_file + "/sorted_sample.bam -O " + utils.temp_file + "/sorted_markedup.bam -M " + utils.temp_file + "/sorted_markedup_metric.bam",
        "gatk AddOrReplaceReadGroups -I  " + utils.temp_file + "/sorted_markedup.bam -O  " + utils.temp_file + "/sorted_markedup_Add.bam -ID 4 -LB lib1 -PL illumina -PU unit1 -SM Cancer",
        "gatk BaseRecalibrator -I " + utils.temp_file + "/sorted_markedup_Add.bam  -R  reference_file/hg38/hg38.fa  --known-sites reference_file/dbsnp_146.hg38.vcf -O " + utils.temp_file + "/recal_data-sample.table",
        "gatk ApplyBQSR -R reference_file/hg38/hg38.fa -I  " + utils.temp_file + "/sorted_markedup_Add.bam --bqsr-recal-file  " + utils.temp_file + "/recal_data-sample.table -O  " + utils.temp_file + "/recal-sample.bam"
    ]
    run_command(command)


def handle04_sam():
    command = [
        "samtools index " + utils.temp_file + "/recal-sample.bam",
        "bcftools mpileup -Ou -f reference_file/hg38/hg38.fa " + utils.temp_file + "/recal-sample.bam | bcftools call -vmO z -o " + utils.temp_file + "/sample.vcf.gz",
        "bcftools tabix -p vcf " + utils.temp_file + "/sample.vcf.gz",
        "bcftools filter -O z -o " + utils.temp_file + "/sample.filtered.vcf.gz -s LOWQUAL -i '%QUAL>20' " + utils.temp_file + "/sample.vcf.gz"
    ]
    run_command(command)


def handle05_annovar():
    command = ["gunzip " + utils.temp_file + "/sample.filtered.vcf.gz",
               "perl /data/wgz/software/annovar/table_annovar.pl " + utils.temp_file + "/sample.filtered.vcf /data/wgz/software/annovar/humandb1/ -buildver hg38 -out " +
               utils.out_file + "/rna_annovar_out -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation g,r,f,f,f -nastring . -vcfinput"]
    run_command(command)


def handle06_var_seq():
    utils.handle1(sample, reuniprot, 'n')
    utils.handle2(sample, reuniprot, 'n')
    utils.handle3(sample, reuniprot, 'n')
    utils.handle4(sample, reuniprot, 'n')
    utils.save_var_pro()
    if not os.path.exists('filter-mass'):    # 创建质谱过滤文件夹
        os.mkdir('filter-mass')
    command = ["cp " + utils.out_file + "/Varsequence.fasta filter-mass/Varsequence.fasta"]
    run_command(command)


def handle07_fusion():
    # 融合新抗原的检测
    if not os.path.exists(utils.temp_file + '/fusion_outdir'):
        os.mkdir(utils.temp_file + '/fusion_outdir')
    if not os.path.exists(utils.out_file + '/fusion'):
        os.mkdir(utils.out_file + '/fusion')
    # command = ["STAR-Fusion --genome_lib_dir reference_file/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir --left_fq input/rna/" +
    #            tumor_files[0] + " --right_fq input/rna/" + tumor_files[1] + " --examine_coding_effect --FusionInspector inspect --extract_fusion_reads --output_dir " +
    #            utils.temp_file + "/fusion_outdir"]
    # run_command(command)
    fusion = pd.read_table(utils.temp_file + '/fusion_outdir/star-fusion.fusion_predictions.abridged.coding_effect.tsv', sep='\t')
    fusion = fusion[fusion['FUSION_TRANSL'].str.contains(r'^[A-Z]')]  # 取出含氨基酸序列的行
    fusion['title'] = fusion.apply(lambda row: '>' + row['#FusionName'], axis=1)
    fusion_fasta = utils.output_peps(fusion[['title', 'FUSION_TRANSL']].copy())
    fusion_fasta.to_csv(utils.out_file + '/fusion/fusion-pep.fasta', header=None, index=False)   # 用于新抗原预测


def handle08_HLA_I_pred():
    if not os.path.exists(utils.temp_file + '/hla'):
        os.mkdir(utils.temp_file + '/hla')
    # MHC Ⅰ:Optitype  预测HLA分型
    command1 = [
        "razers3 -i 95 -m 1 -dr 0 -tc 8 -o " + utils.temp_file + "/Tumor_1.bam ./software/OptiType-1.3.5/data/hla_reference_rna.fasta input/rna/"+ tumor_files[0],
        "razers3 -i 95 -m 1 -dr 0 -tc 8 -o " + utils.temp_file + "/Tumor_2.bam ./software/OptiType-1.3.5/data/hla_reference_rna.fasta input/rna/" + tumor_files[1],
        "samtools fastq " + utils.temp_file + "/Tumor_1.bam > " + utils.temp_file + "/Tumor_1_fished.fastq",
        "samtools fastq " + utils.temp_file + "/Tumor_2.bam > " + utils.temp_file + "/Tumor_2_fished.fastq",
        # 得到预测分型的文件夹存于outfile-rna/hla，用户需要将预测得到的分型手动输入到.py文件中
        "python ./software/OptiType-1.3.5/OptiTypePipeline.py  -i " + utils.temp_file + "/Tumor_1_fished.fastq " + utils.temp_file + "/Tumor_2_fished.fastq  -r -v -o " + utils.out_file + "/hla",
    ]
    run_command(command1)


def handle08_HLA_II_pred():
    # NHC Ⅱ:HLAminer
    # 1. bwa比对  2. 分型
    command2 = [
        "bwa mem -a /data/wgz/software/HLAminer-1.4/HLAminer-1.4/database_bwamem/HLA-I_II_CDS.fasta input/rna/" + tumor_files[0] + " input/rna/" + tumor_files[1] + " > " + utils.temp_file + "/HLA_II.sam",
        "perl /data/wgz/software/HLAminer-1.4/HLAminer_1.4/bin/HLAminer.pl -a " + utils.temp_file + "HLA-II.sam -h /data/wgz/software/HLAminer-1.4/HLAminer_1.4/database/HLA-I_II_CDS.fasta -p /data/wgz/software/HLAminer-1.4/HLAminer_1.4/database/hla_nom_p.txt"
    ]
    run_command(command2)


def handle09_tpm():
    # 预测TPM
    # 得到的文件存在:outfile-rna/tpm/abundance.tsv
    if not os.path.exists(utils.out_file + '/tpm'):
        os.mkdir(utils.out_file + '/tpm')
    command = ["kallisto quant -i reference_file/hg38/h38.cdna.idx -o " + utils.out_file + "/tpm -b 100 input/rna/" + tumor_files[0] + " input/rna/" + tumor_files[1] + " > tpm.out"]
    run_command(command)


# not been used
def txt2fasta():
    txt_file = pd.read_table(utils.temp_file + 'HIQ_1223_snv_1.txt', sep='\t')
    txt_file.drop_duplicates(subset=['Peptide'], inplace=True)
    txt_file.reset_index(inplace=True)
    txt_file['num'] = txt_file['MHC'].index + 1
    txt_file['title'] = txt_file.apply(lambda row: '>' + str(row['num']) + '|' + row['MHC'] + '|' + row['gene'] + '|' + row['AAchange'] + '|' + str(row['Pos']), axis=1)
    fasta_file = utils.output_peps(txt_file[['title', 'Peptide']].copy())
    fasta_file.to_csv('HIQ_1223_snv_1.fasta', index=False, header=None, sep='\t')


if __name__ == '__main__':
    utils = NeoantigenUtils("outfile2", "outfile-rna")
    # 使用dirname函数访问上级目录
    current_path = os.path.abspath(os.path.dirname(__file__))
    print(current_path)
    if not os.path.exists(utils.out_file):
        os.mkdir(utils.out_file)
        dir_list = ['/txt_files', '/fasta_files', '/csv_files']
        for d in dir_list:
            os.mkdir(utils.out_file + d)
    tumor_files = os.listdir('input/rna')
    activate_env()
    # handle01_bwa()
    # handle02_sam()
    # handle03_gatk()
    # handle04_sam()
    # handle05_annovar()
    reuniprot = utils.read_fasta()
    # 读取样本文件
    sample = pd.read_table(utils.out_file + '/rna_annovar_out.hg38_multianno.txt', sep="\t")
    sample = sample[sample["Otherinfo10"].str.contains("PASS")]
    # handle06_var_seq()
    handle07_fusion()
    # handle08_HLA_I_pred()
    # handle08_HLA_II_pred()
    # handle09_tpm()
