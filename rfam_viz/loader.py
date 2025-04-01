#########################################################
#                        std Lib                        #
#########################################################
import os
from dataclasses import dataclass
from typing import List, Generator, Any
from pathlib import Path


#########################################################
#                      Dependencies                     #
#########################################################
from Bio import SeqIO
import Bio
from Bio.Align import Alignment
import regex
import torch
from torch.utils.data import Dataset
#########################################################
#                      Own modules                      #
#########################################################

@dataclass
class Family:
    """ Dataclass resulting from the parsing of Stockholm MSA """
    path: str
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
            path=fam_ali_file,
            family=fam[0],
            consensus_seq=regex.sub("\.", "-", consensus_seq[0]),
            consensus_ss=consensus_ss[0],
            msa=msa
        ), len(split_file)

class DatasetDCA(Dataset):
    """Dataset class for handling multi-sequence alignments data."""
    def __init__(
        self,
        alphabet: str = '-AUCG',
        device : str = "cpu"
    ):
        """
        Initialize the dataset.
        """
        self.device = device

        self.alphabet = alphabet
        self.dtype = torch.float32

        self.msa = []
        self.mat = []

    def get_format(self, seqlist : List[str]):
        for record in seqlist:
            self.msa.append(record)
            encoded_seq = encode_sequence(record)
            
            if encoded_seq:
                self.mat.append(encoded_seq)
            else:
                print(f"Invalid character in sequence. Sequence has been removed.")
        
        self.encoded = self.mat.copy()       
        self.mat = one_hot(mat=self.mat, device=self.device, num_classes=len(self.alphabet))
        return self.mat

        """
            nseq = number of sequences in the MSA
            nnuc = sequences lengths in the MSA (i.e., number of nucleotides)
            nval = number of different nucleotides elements        
        """

    
def encode_sequence(seq : str | torch.Tensor):
    dico = {'-':0,'A':1,'U':2,'C':3,'G':4}
    new = []
    for nuc in seq:
        try :
            new.append(dico[nuc])
        except KeyError:
            return False
    
    return new
    
@torch.jit.script
def _one_hot(x: torch.Tensor, num_classes: int = -1, dtype: torch.dtype = torch.float32):
   
    if x.dim() != 2:
        raise ValueError(f"Input tensor x must be 2D, currently {x.shape}")
    
    if num_classes < 0:
        num_classes = x.max() + 1
    res = torch.zeros(x.shape[0], x.shape[1], num_classes, device=x.device, dtype=dtype)
    tmp = torch.meshgrid(
        torch.arange(x.shape[0], device=x.device),
        torch.arange(x.shape[1], device=x.device),
        indexing="ij",
    )
    index = (tmp[0], tmp[1], x)
    values = torch.ones(x.shape[0], x.shape[1], device=x.device, dtype=dtype)
    res.index_put_(index, values)
    
    return res


def one_hot(mat: List[List[int]] | torch.Tensor, device : str,  num_classes: int = -1, dtype: torch.dtype = torch.float32):
    """
    A fast one-hot encoding function faster than the PyTorch one working with torch.int32 and returning a float Tensor.
    Works only for 2D tensors.
    """
    if isinstance(mat, List):
        mat = torch.tensor(mat, device=device, dtype=torch.int32)
    return _one_hot(mat, num_classes, dtype)