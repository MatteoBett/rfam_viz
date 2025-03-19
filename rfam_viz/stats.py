#########################################################
#                        std Lib                        #
#########################################################
import re
from collections import Counter
from typing import List

#########################################################
#                      Dependencies                     #
#########################################################
import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
import scipy.signal

#########################################################
#                      Own modules                      #
#########################################################
from rfam_viz.loader import Family

@dataclass
class FamilyStats:
    gap_freq: List[float]
    gap_freq_std : List[float]
    gap_per_seq : List[float]
    gps_std : List[float]
    prominent_fraction: List[float]
    gaps_per_seq_pf : List[float]
    num_seq : List[int]
    len_msa: List[int]
    family: List[str]
    var_len : List[str]
    max_var : List[int]
    avg_size : List[float]
    

@dataclass
class GapDist:
    all_dist : List[List[float]]
    family : List[str]
    disrupt_score : List[float]

def get_stats(msa):

    pass

def freq_gaps_msa(family_record : Family, famstats : FamilyStats, gapdist : GapDist) -> float:
    all_val = []
    all_freq = []
    seqlen = []
    size = len(family_record.consensus_seq)
    gapcount = np.array([0]*size)
    all_disrupt_score = np.zeros((len(family_record.msa), size), dtype=np.float32)
    for i, seq in enumerate(family_record.msa): 
        hamming = 0
        
        for index, elt in enumerate(seq):
            if elt == '-':
                gapcount[index] += 1
            if elt != family_record.consensus_seq[index] and family_record.consensus_seq[index] == '-':
                hamming += 1
            elif elt != '-':
                all_disrupt_score[i, index] = 1
        
        all_freq.append(seq.count('-')/size)
        all_val.append(hamming)
        seqlen.append(len(re.sub('-', "", str(seq))))

    if gapcount.sum() == 0 or len(family_record.msa) <= 50:
        return famstats, gapdist
    
    gapcount = gapcount/gapcount.sum() 
    all_disrupt_score = all_disrupt_score * gapcount
    all_disrupt_score = all_disrupt_score.sum(axis=1)

    gapcount.sort()
    gapdist.all_dist.append(gapcount)
    gapdist.family.append(family_record.family)
    gapdist.disrupt_score.append(all_disrupt_score)

    avg = np.average(seqlen)
    maxlen = max(seqlen)
    seqlen = [elt/maxlen for elt in seqlen]
    std = np.std(all_val)
    prominent = len([i for i in all_val if abs(i - np.mean(all_val)) > 2*std])

    famstats.family.append(family_record.family)
    famstats.num_seq.append(len(family_record.msa))
    famstats.len_msa.append(len(family_record.msa[0]))
    famstats.gap_freq.append(np.mean(all_freq))
    famstats.gap_freq_std.append(np.std(all_freq))
    famstats.gap_per_seq.append(np.mean(all_val))
    famstats.gps_std.append(std)
    famstats.prominent_fraction.append(prominent/len(family_record.msa))
    famstats.gaps_per_seq_pf.append(np.mean(prominent))
    famstats.var_len.append(np.std(seqlen))
    famstats.max_var.append(max(seqlen) - min(seqlen))
    famstats.avg_size.append(avg)
    #print(family_record.family, max(seqlen), min(seqlen))
    return famstats, gapdist

def split_list(lst : List[np.ndarray], num_lists):
    # Calculate the size of each chunk, and any leftover items
    avg_len = len(lst) // num_lists
    remainder = len(lst) % num_lists

    result = []
    start_index = 0
    
    for i in range(num_lists):
        # Calculate the chunk size, adding an extra item if there's a remainder
        chunk_size = avg_len + (1 if i < remainder else 0)
        
        result.append(np.sum(lst[start_index:start_index + chunk_size]))
        start_index += chunk_size

    return result

def distance_gap(gapdist : GapDist, num_chunks : int = 20):
    """ 
    Norme L2 entre les distributions de gaps et tailles de séquences 
    entre toutes les familles.

    Pas de 5% par défaut
    """
    chunks_list = []
    for dist in gapdist.all_dist:
        chunks_list.append(split_list(lst=dist, num_lists=num_chunks))
    
    matdist = np.full((len(gapdist.all_dist), len(gapdist.all_dist)), 0.0, dtype=np.float32)
    l2_norm = lambda x,y : np.sqrt((x - y)**2)

    for i in range(len(gapdist.all_dist)):
        for j in range(len(gapdist.all_dist)):
            if i != j:
                for piece in range(num_chunks):
                    matdist[i,j] += l2_norm(x=chunks_list[i][piece], y=chunks_list[j][piece])
    return matdist

def get_disruption_dist(disruption_score : List[List[float]]):
    pass

def get_outliers(matdist : np.matrix, gapdist : GapDist, min_prominence : float = 0.75):
    all_mean = matdist.mean(axis=1)
    highest_values_indices, _ = scipy.signal.find_peaks(all_mean, prominence=min_prominence, height=None)
    outliers = {gapdist.family[i]:split_list(gapdist.all_dist[i], 20) for i in highest_values_indices}

    homogeneous_indices = all_mean.argsort()[:-len(highest_values_indices)]
    chunks_list = []
    for dist in gapdist.all_dist:
        chunks_list.append(split_list(lst=dist, num_lists=20))    

    chunks_list = np.array(chunks_list)
    piece = len(homogeneous_indices)//4
    homogeneous = np.split(homogeneous_indices, [piece, piece*2, piece*3])

    homogeneous = [chunks_list[split] for split in homogeneous]
    
    disruption_dist_homogenous = get_disruption_dist(disruption_score=gapdist.disrupt_score)
    disruption_dist_outliers = {gapdist.family[i]:gapdist.disrupt_score[i] for i in highest_values_indices}
    return outliers, homogeneous, disruption_dist_homogenous, disruption_dist_outliers


def do_stats(family_record : Family, famstats : FamilyStats, gapdist : GapDist):
    famstats, gapdist = freq_gaps_msa(family_record=family_record, famstats=famstats, gapdist=gapdist)
    return famstats, gapdist

