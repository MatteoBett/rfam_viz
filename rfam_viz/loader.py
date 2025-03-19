#########################################################
#                        std Lib                        #
#########################################################
import os
from dataclasses import dataclass
from typing import List, Generator, Any

#########################################################
#                      Dependencies                     #
#########################################################
from Bio import SeqIO
from Bio.Align import Alignment
import regex
#########################################################
#                      Own modules                      #
#########################################################

@dataclass
class Family:
    """ Dataclass resulting from the parsing of Stockholm MSA """
    msa: List[str]
    family: str
    consensus_ss: str
    consensus_seq: str


def load_fam(fam_ali_file : str) -> Generator[str, Any, Any]:
    with open(fam_ali_file, 'r', encoding="utf-8") as stock_msa:
        full_file = stock_msa.read()

    split_file = full_file.split('# STOCKHOLM 1.0')
    for piece in split_file:
        if len(piece) < 10:
            continue
        
        fam = regex.findall(r'#=GF DE\s+(.*?)(?=\n)', piece)
        consensus_ss = regex.findall(r'#=GC SS_cons\s+(.*?)(?=\n)', piece)
        consensus_seq = regex.findall(r'#=GC RF\s+(.*?)(?=\n)', piece)
        msa = [elt for elt in regex.findall(r'\s+([AUCG-]+)?(?=\n)', piece) if len(elt) > 5]
        
        yield Family(
            family=fam[0],
            consensus_seq=regex.sub("\.", "-", consensus_seq[0]),
            consensus_ss=consensus_ss[0],
            msa=msa
        ), len(split_file)
    