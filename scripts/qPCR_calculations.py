#! /usr/bin/env python3

import argparse
import os
import glob
from subprocess import call
import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import linregress
import seaborn as sns
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
sns.set_style('white')


def qPCR_calculator(original_qPCR_Ct_dir, output_dir):

    # Pre-Step 1a: rename
    original_qPCR_Ct_dir=os.path.abspath(original_qPCR_Ct_dir)

    if len(glob.glob("{}/*Quantitation Ct*".format(original_qPCR_Ct_dir))) > 1:
        print()
        print("!!!")
        print("Error message: there are more than one Quantitiation Ct Results file in the given directory. Please make sure there is one.")
        print("!!!")
        print()
        return

    else:
        original_name = glob.glob("{}/*Quantitation Ct*".format(original_qPCR_Ct_dir))[0]
        original_file_name=os.path.basename(original_name)
        new_name = "_".join(os.path.basename(original_name).split(" ")[4:7])
        new_name_full_path = os.path.dirname(original_name) + '/intermediats/' + new_name

        mk_new_dir_cmd='mkdir {}/intermediats'.format(original_qPCR_Ct_dir)

        if os.path.isdir('{}/intermediats'.format(original_qPCR_Ct_dir)) == True:
            print("intermediates directory already exists! Skipping...")
        else:
            call(mk_new_dir_cmd, shell=True)


        rename_cmd = 'cp "{}" {}'.format(original_name, new_name_full_path)
        call(rename_cmd, shell=True)

        mk_output_result_cmd = 'mkdir {}/results'.format(output_dir)

        if os.path.isdir('{}/results'.format(output_dir)) == True:
            print("the results directory already exists! Skipping...")
        else:
            call(mk_output_result_cmd, shell=True)

        try:
            CTdb = pd.read_csv(new_name_full_path)

        except Exception:
            # ODS to XLSX conversion

            import jpype
            import asposecells
            jpype.startJVM()
            from asposecells.api import Workbook

            workbook = Workbook(new_name_full_path)
            workbook.save(new_name_full_path.replace('xlsx', 'xls'))
            jpype.shutdownJVM()

            CTdb = pd.read_excel(new_name_full_path.replace('xlsx', 'xls'))

        CTdb = CTdb.drop(columns=list(CTdb.filter(regex='Unnamed')))
        CTdb['Log Starting Quantity'] = [float(x) for x in CTdb['Log Starting Quantity']]
        CTdb['original_file_name'] = original_file_name

        # Step 1: Calculate intercept & slope

        print("Step 1: Calculate intercept & slope")

        CTdb_stds = CTdb[(CTdb['Content'].isin(['Std', 'NTC'])) & \
                         (CTdb['Fluor'].isin(['SYBR']))].reset_index(drop=True)
        CTdb_stds['Sample'] = [x.upper() if y != 'NTC' else 'S7' for x, y in zip(CTdb_stds['Sample'], CTdb_stds['Content'])]

        Qdb = CTdb_stds[['Content', 'Sample', 'Starting Quantity (SQ)', 'Log Starting Quantity']].drop_duplicates()
        Qdb = Qdb[-Qdb['Content'].isin(['NTC'])].reset_index(drop=True)

        std2ct_mean = dict()
        for sample, sdb in CTdb_stds.groupby('Sample'):
            ave_ct = np.mean(sdb['C(t) Mean'])
            std2ct_mean[sample] = ave_ct

        Qdb['Average Cq'] = Qdb['Sample'].map(std2ct_mean)

        std2delta_ct = {}
        for i, k in enumerate(std2ct_mean):
            # print(i,k)
            if i < 5:
                delta_ct = std2ct_mean[list(std2ct_mean.keys())[i + 1]] - std2ct_mean[list(std2ct_mean.keys())[i]]
                std2delta_ct[list(std2ct_mean.keys())[i + 1]] = delta_ct
        Qdb['Delta Cq'] = Qdb['Sample'].map(std2delta_ct)

        res = linregress(Qdb['Log Starting Quantity'], \
                         Qdb['Average Cq'])

        table = defaultdict(list)
        for i, item in enumerate(['slope', 'intercept', 'rvalue', 'pvalue', 'stderr']):
            table[item].append(res[i])

        Rdb = pd.DataFrame(table)

        r_sq = round(res.rvalue ** 2, 6)
        efficiency = (10 ** (-1 / res.slope) - 1) * 100

        Rdb['R-squared'] = r_sq
        Rdb['efficiency (%)'] = efficiency

        Rdb.to_csv(os.path.abspath(output_dir)+'/results/qPCR_Quantification_metrics.csv',index=False)

        # plot

        from matplotlib.ticker import MultipleLocator
        from matplotlib.ticker import MaxNLocator

        fig = plt.figure(figsize=(6.5, 4.5))
        ax = plt.subplot(111)

        ax.xaxis.set_major_locator(MaxNLocator(6))
        ax.yaxis.set_major_locator(MaxNLocator(4))

        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_minor_locator(MultipleLocator(1))

        plt.plot(Qdb['Log Starting Quantity'],Qdb['Average Cq'], 'o')

        for idx, val in enumerate(Qdb['Average Cq']):
            ax.text(-idx + 1.1, val + 1.5, str(round(val, 2)))

        plt.plot(Qdb['Log Starting Quantity'],res.intercept + res.slope * ((Qdb['Log Starting Quantity'])),
                 'r', label='y={:.4f}x+{:.4f}; \nR2={:.6f}'.format(res.slope, res.intercept, r_sq))
        plt.axvline(x=0, color='gray', linestyle='--')
        plt.legend(fontsize=11, frameon=True, bbox_to_anchor=(1.01, 1.025))

        plt.xlabel('Log (Conc in pM)')
        plt.ylabel('Average Cq')
        plt.ylim(2, 35)
        plt.xlim(-4.5, 2.5)
        plt.tight_layout()
        plt.savefig(os.path.abspath(output_dir)+'/results/qPCR_Quantification_metrics.pdf',bbox_inches='tight',format='pdf')

        # Step 2: Calculate library concentration

        print("Step 2: Calculate library concentration")

        CTdb_sp = CTdb[(CTdb['Content'].isin(['Unkn'])) & \
                       (CTdb['Fluor'].isin(['SYBR']))].reset_index(drop=True)

        Cdb = CTdb_sp[['Sample', 'C(t) Mean','original_file_name']]
        Cdb['Sample'] = [float(".".join(x.split("-"))) for x in Cdb['Sample']]
        Cdb.insert(1, 'Dilution', 10 ** 4)
        Cdb.insert(2, 'Average fragment length (bp)', 650)

        sp_std2ct_mean = dict()
        for sample, sdb in Cdb.groupby('Sample'):
            ave_ct = np.mean(sdb['C(t) Mean'])
            sp_std2ct_mean[sample] = ave_ct

        Cdb['Average Cq'] = Cdb['Sample'].map(sp_std2ct_mean)
        Cdb = Cdb.drop_duplicates(subset=['Sample', 'Average Cq']).drop(columns=['C(t) Mean']). \
            sort_values('Sample').reset_index(drop=True)
        Cdb['log (concentration)'] = [(x - res.intercept) / res.slope for x in Cdb['Average Cq']]
        Cdb['Average concentration (pM)'] = [10 ** x for x in Cdb['log (concentration)']]
        Cdb['Size-adjusted concentration (pM)'] = [x * (452 / y) for x, y in zip(Cdb['Average concentration (pM)'], \
                                                                                 Cdb['Average fragment length (bp)'])]
        Cdb['Concentration  of undiluted library (nM)'] = [x * y / 1000 for x, y in
                                                           zip(Cdb['Size-adjusted concentration (pM)'], \
                                                               Cdb['Dilution'])]
        Cdb_final = Cdb[['Sample', 'Concentration  of undiluted library (nM)','original_file_name']]

        rm_intermediate_dir_cmd = 'rm -r {}/intermediats'.format(original_qPCR_Ct_dir)
        call(rm_intermediate_dir_cmd,shell=True)

        Cdb_final.to_csv(os.path.abspath(output_dir)+'/results/qPCR_LibConc.csv', index=False)
        return Cdb_final


# if __name__ == '__main__':
#
#     parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,add_help=False,
#                                      description='Calculate library concentrations from qPCR data!')
#
#     # General Arguments
#     GenArgs = parser.add_argument_group('GENERAL ARGUMENTS')
#     GenArgs.add_argument('-h', action="help",help="show this help message and exit")
#     GenArgs.add_argument('-i','--qPCR_dir', help='the directory path where the qPCR Ct file was stored')
#     GenArgs.add_argument('-o','--output_dir', default='./', help='Output file directory')
#     args = parser.parse_args()
#
#
#     qPCR_calculator(args.qPCR_dir, args.output_dir)
