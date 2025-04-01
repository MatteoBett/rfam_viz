#########################################################
#                        std Lib                        #
#########################################################
import os, sys, argparse
from dataclasses import dataclass
from typing import List

#########################################################
#                      Dependencies                     #
#########################################################
import pandas as pd

#########################################################
#                      Own modules                      #
#########################################################
import rfam_viz.loader as seqload
import rfam_viz.stats as stats
import rfam_viz.utils as utils
import rfam_viz.display as display

"""
Environment of use: rfam_viz
"""

if __name__ == "__main__":
    data = "./data/Rfam_own.seed"
    outfig = "./output/figures"

    do_one = False

    init_fam = stats.FamilyStats(
        gap_freq=[],
        gap_freq_std=[],
        gap_per_seq=[],
        gps_std=[],
        prominent_fraction=[],
        gaps_per_seq_pf=[],
        num_seq=[],
        len_msa=[],
        family=[],
        var_len=[],
        max_var=[],
        avg_size=[],
        disrupt_score=[],
        prominent_disrupt=[],
        concomitance_score=[]
    )
    init_gapdist = stats.GapDist(
        all_dist=[],
        family=[],
        
    )
    for index, (family_record, num_fam) in enumerate(seqload.load_fam(fam_ali_file=data)):
        utils.progressbar(iteration=index, total=num_fam)
        fam_stats, gapdist = stats.do_stats(family_record=family_record, famstats=init_fam, gapdist=init_gapdist)
        if do_one:
            exit(0)

    matdist = stats.distance_gap(gapdist=gapdist)
    outliers, homogeneous= stats.get_outliers(matdist=matdist, gapdist=gapdist)

    df = pd.DataFrame(data=fam_stats.__dict__)
    df = df[df["num_seq"] > 50]

    df = df.sort_values(by='avg_size', ascending=False)
    display.make_overview(df=df, 
                          matdist=matdist, 
                          outliers=outliers, 
                          homogeneous=homogeneous, 
                          figdir=outfig)
