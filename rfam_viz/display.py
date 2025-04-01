#########################################################
#                        std Lib                        #
#########################################################
import os
from typing import List, Dict

#########################################################
#                      Dependencies                     #
#########################################################
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bpdf
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.signal
#########################################################
#                      Own modules                      #
#########################################################

def gap_fraction(df : pd.DataFrame, pdf):
    gf, gf_std = df['gap_freq'].mean(), df["gap_freq"].std()
    bf, bf_std = 1-gf, np.std([1-i for i in df["gap_freq"]])
    fig, axes = plt.subplots(1, 1, figsize=(12, 8))
    axes.bar(x = ["gaps", "bases"], height=[gf, bf], color=["blue", "green"])
    axes.errorbar(x = ["gaps", "bases"], y = [gf, bf], yerr=[gf_std, bf_std], fmt="o", color='black')
    axes.set_ylabel("Frequency")

    axes.set_title('Rfam Seed alignment overall gap fraction')

    fig.savefig(pdf, format="pdf")
    plt.close(fig)

def max_var_lenseq(df : pd.DataFrame, pdf):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    sns.scatterplot(data=df, x='avg_size', y="max_var", ax=axes[0])
    sns.scatterplot(data=df, x='len_msa', y="max_var", ax=axes[1])

    axes[0].set_title('Maximum variation depending on the average sequence size')
    axes[1].set_title('Maximum variation depending on the MSA length')

    fig.suptitle("Maximum variation in gaps' number per sequence depending on MSA's and sequence average size")
    fig.subplots_adjust(wspace=0.2)
    fig.savefig(pdf, format='pdf')
    plt.close(fig)

def gap_distribution_differences(matdist : np.matrix, pdf):
    fig, axes = plt.subplots(1, 1, figsize=(12, 10))
    sns.heatmap(data=matdist, cmap='jet', cbar=True, ax=axes)

    axes.set_title('L2 norm of the gap distribution differences in Rfam')
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format="pdf")
    plt.close(fig)

def show_outliers_dist(outliers : Dict[str, np.ndarray], pdf):
    col = len(outliers)//2
    fig, axes = plt.subplots(2, col, figsize=(18, 6))
    
    axes = axes.flat
    for index, (key,val) in enumerate(outliers.items()):
        axes[index].plot(val, marker='.')
        axes[index].set_xlabel(f"{100//len(val)}% quantiles")
        axes[index].set_ylabel("Fraction of gaps in this index")
        axes[index].set_title(key)
        axes[index].set_ylim((0.0, 1.0))
    
    fig.legend(["gap fraction"])
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format="pdf")
    plt.close(fig)

def show_consensus_dist(homogenous : List[np.ndarray], pdf):
    fig, axes = plt.subplots(2, 2, figsize=(18, 6))
    axes = axes.flat
    
    for index, mat in enumerate(homogenous):
        data = mat.mean(axis=0)
        std = mat.std(axis=0)
        #sns.scatterplot(data=mat.T, ax=axes[index], color='r', alpha=0.4, legend=False)
        axes[index].plot(list(range(len(data))), data, marker='.', color='blue')
        axes[index].errorbar(list(range(len(data))), data, std)
        axes[index].set_title(f"Gaps average distribution quartile number {index}")
        axes[index].set_xlabel(f"{100//len(data)}% quantiles")
        axes[index].set_ylabel(f"Fraction of gaps in this index")
    
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format="pdf")
    plt.close(fig)

def scatter_corr(df : pd.DataFrame, pdf):
    fig, axes = plt.subplots(2, 2, figsize=(13, 5))
    axes= axes.flat

    sns.scatterplot(data=df, x='avg_size', y="gap_freq", ax=axes[0])
    sns.scatterplot(data=df, x='len_msa', y="gap_freq", ax=axes[1])
    sns.scatterplot(data=df, x='avg_size', y="gap_freq_std", ax=axes[2])
    sns.scatterplot(data=df, x='var_len', y="gap_freq_std", ax=axes[3])

    axes[0].set_title('Gap fraction depending on the average sequence size')
    axes[1].set_title('Gap fraction depending on the MSA length')
    axes[2].set_title('Variation of Gap fraction depending on average sequence size')
    axes[3].set_title('Variation of Gap fraction depending on sequences size standard deviation')

    """    
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    axes[2].set_xscale('log')
    axes[3].set_xscale('log')
    """
    fig.suptitle("Correlation of gap frequency with sequences and MSA size depending on MSA's and sequence average size")
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format='pdf')
    plt.close(fig)


def disruption_score_display(disruption_dist_outliers : Dict[str, List[float]], pdf):
    col = len(disruption_dist_outliers)//2
    fig, axes = plt.subplots(2, col, figsize=(18, 6))
    
    axes = axes.flat
    for index, (key,val) in enumerate(disruption_dist_outliers.items()):
        axes[index].plot(val, marker='.')
        axes[index].set_xlabel(f"")
        axes[index].set_ylabel("")
        axes[index].set_title(key)
    
    fig.legend(["gap fraction"])
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format="pdf")
    plt.close(fig)

def gaps_origin(df : pd.DataFrame, pdf):
    fig, axes = plt.subplots(2, 2, figsize=(13, 5))
    axes= axes.flat

    sns.scatterplot(data=df, x='gap_per_seq', y="avg_size", ax=axes[0])
    sns.scatterplot(data=df, x='gps_std', y="avg_size", ax=axes[1])
    sns.scatterplot(data=df, x='prominent_fraction', y="avg_size", ax=axes[2])
    sns.scatterplot(data=df, x='gaps_per_seq_pf', y="avg_size", ax=axes[3])

    axes[0].set_title('Correlation between the number of gaps created by a sequence in a MSA\nand the average sequence size')
    axes[1].set_title('Correlation of the standard deviation')
    axes[2].set_title('Correlation between the size of a prominent fraction\nand the average sequence size')
    axes[3].set_title('Correlation between the number of gaps created by the sequences\nof the prominent fraction and the average sequence size')
    axes[0].set_xscale('log')
    axes[1].set_xscale('log')
    axes[2].set_xscale('log')
    axes[3].set_xscale('log')
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[2].set_yscale('log')
    axes[3].set_yscale('log')
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format='pdf')
    plt.close(fig)

def seq_disruption(df : pd.DataFrame, pdf):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes= axes.flat

    sns.scatterplot(data=df, x='disrupt_score', y="avg_size", ax=axes[0])
    sns.scatterplot(data=df, x='prominent_disrupt', y="avg_size", ax=axes[1])


    axes[0].set_title('Correlation between a family disruption score\nand the average sequence size')
    axes[1].set_title('Correlation between a family disruption score\nand the gaps creation in the sequences prominent fraction')

    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format='pdf')
    plt.close(fig)

def make_overview(df : pd.DataFrame, 
                  matdist : np.matrix, 
                  outliers : Dict[str, np.ndarray], 
                  homogeneous : List[np.ndarray], 
                  figdir : str):
    pdf = bpdf.PdfPages(os.path.join(figdir, "EDA_1_global_barplot.pdf"))
    
    for col in df.columns:
        if col == 'family':
            continue
        if col == "concomitance_score":
            fig, axes = plt.subplots(1, 1, figsize=(13, 5))

            sns.barplot(df,x = df.index, y = col, ax=axes, errorbar=('ci', 95))

            axes.set_yscale("log")
            fig.suptitle(f'Numeric Feature : {col}', fontsize=16, fontweight='bold')
            fig.subplots_adjust(wspace=0.2)
            fig.savefig(pdf, format='pdf')
            plt.close(fig)
            continue

        fig, axes = plt.subplots(1, 1, figsize=(13, 5))

        sns.barplot(df,x = df.index, y = col, ax=axes, errorbar=('ci', 95))

        fig.suptitle(f'Numeric Feature : {col}', fontsize=16, fontweight='bold')
        fig.subplots_adjust(wspace=0.2)
        fig.savefig(pdf, format='pdf')
        plt.close(fig)

    pdf.close()

    pdf = bpdf.PdfPages(os.path.join(figdir, "EDA_2_correlation.pdf"))
    df = df.drop('family', axis=1)
    fig, axes = plt.subplots(1, 1, figsize=(12, 8))
    sns.heatmap(data=df.corr('pearson'), annot=True, cmap='jet', cbar=True, ax=axes)
    fig.suptitle("Pearson correlation (linear)")
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format="pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 1, figsize=(12, 8))
    sns.heatmap(data=df.corr('spearman'), annot=True, cmap='jet', cbar=True, ax=axes)
    fig.suptitle("Spearman correlation (non linear)")
    fig.subplots_adjust(wspace=0.5, hspace=0.5)
    fig.savefig(pdf, format="pdf")
    plt.close(fig)

    gap_fraction(df, pdf=pdf)
    max_var_lenseq(df, pdf=pdf)
    scatter_corr(df=df, pdf=pdf)
    gaps_origin(df=df, pdf=pdf)
    seq_disruption(df=df, pdf=pdf)
    pdf.close()

    pdf = bpdf.PdfPages(os.path.join(figdir, "EDA_3_distance.pdf"))
    gap_distribution_differences(matdist=matdist, pdf=pdf)
    show_outliers_dist(outliers=outliers, pdf=pdf)
    show_consensus_dist(homogenous=homogeneous, pdf=pdf)
    pdf.close()

