import os, re

import numpy as np

tokens = ["-", "A", "U", "C", "G"]

normal_seq = ["A"]*25 + ['G']*25 + ["C"]*25 + ["U"]*25
begin_seq = ["-"]*25 + ['G']*25 + ["C"]*25 + ["U"]*25
end_seq = ["A"]*25 + ['G']*25 + ["-"]*25 + ["U"]*25
both = ["-"]*25 + ['G']*25 + ["-"]*25 + ["U"]*25


def make_ends(dist : int, len_seq : int):
    term1_start = np.random.randint(1, len_seq//4)
    term1_end = term1_start + dist
    term2_start = np.random.randint(term1_end +len_seq//4, len_seq - len_seq//4)
    term2_end = term2_start + dist
    return term1_start, term1_end, term2_start, term2_end

def make_data(with_cc : bool = True, dist : int = 20, datadir : str = r'/home/mbettiati/LBE_MatteoBettiati/code/Rfam_viz/data', n_seq : int = 100, len_seq : int = 200):
    
    if with_cc:
        out = "with_cc"
    else:
        out = "non_cc"
    with open(os.path.join(datadir, f"test_{out}.fasta"), 'w') as writest:
        for i in range(n_seq):
            if with_cc:
                if i < n_seq//4:
                    seq = begin_seq
                elif i > n_seq*3/4:
                    seq = end_seq
                else:
                    seq = normal_seq
            else:
                seq = normal_seq

            writest.write(f">test_{i}\n{"".join(seq)}\n")

make_data(with_cc=True)
make_data(with_cc=False)

            

