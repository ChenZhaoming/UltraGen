import pandas as pd
from tqdm import tqdm
import os

# HTR-selex
primer_5 = 'GGGAUAUCCUCCACGGAGUCGGCAAGCAGAAGACGGCAUACGAU'
primer_3 = 'GAUCGGAAGAGCGUCG'

# primer_5 = 'GGGGAGUUCUACAGUCCGACGAUC'
# primer_3 = 'UGGAAUUCUCGGGUGUCAAGG'

# RBNS-78 proteins
# primer_5 = 'GGGGAGUUCUACAGUCCGACGAUC'
# primer_3 = 'UGGAAUUCUCGGGUGUCAAGG'

d_complement = {'A': 'U', 'T':'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def get_reverse_complement(seq):
    """Function to get reverse complementary of input sequence

    Args:
      seq: Input sequence string

    Returns:
      Reversed complementary string:
    """
    return "".join(d_complement.get(base, base) for base in reversed(seq))

def load_file(file_dir):
    return pd.read_csv(file_dir)

if __name__ == "__main__":

    # IGF2BP1, PUM1, TARDBP

    # barcode_5 = 'TGAAGA'
    # barcode_3 = 'CAG'
    # file_dir = '/share/project/UltraGen/data/raw/PRJEB25907/DAZ3.csv'

    # barcode_5 = 'TGAAGA'
    # barcode_3 = 'ATG'
    # file_dir = '/share/project/UltraGen/data/raw/PRJEB25907/RBM24.csv'

    # barcode_5 = 'TAGCCA'
    # barcode_3 = 'TAT'
    # file_dir = '/share/project/UltraGen/data/raw/PRJEB25907/HNRNPCL1_raw.csv'

    barcode_5 = 'TCAGTA'
    barcode_3 = 'ACC'
    file_dir = '/share/project/UltraGen/data/raw/PRJEB25907/TARDBP.csv'

    # barcode_5 = 'TCTTAA'
    # barcode_3 = 'GCA'
    # file_dir = '/share/project/UltraGen/data/raw/PRJEB25907/IGF2BP1_raw.csv'

    # barcode_5 = 'TCTCGT'
    # barcode_3 = 'GCC'
    # file_dir = '/share/project/UltraGen/data/raw/PRJEB25907/SNRNP70.csv'

    # file_dir = '/share/project/UltraGen/RNA-bind-n-seq/DAZ3.csv'
    # file_dir = '/share/project/UltraGen/RNA-bind-n-seq/RBM24.csv'
    # file_dir = '/share/project/UltraGen/RNA-bind-n-seq/HNRNPCL1.csv'

    # file_dir = '/share/project/UltraGen/RNA-bind-n-seq/FUS.csv'
    # file_dir = '/share/project/UltraGen/RNA-bind-n-seq/TIA1.csv'
    
    df_data = load_file(file_dir)
    for idx in tqdm(df_data.index.tolist()):
        seq = df_data.loc[idx, 'motif']
        # new_seq = primer_5 + get_reverse_complement(barcode_5 + seq + barcode_3) + primer_3
        new_seq = primer_5 + get_reverse_complement(barcode_3) + seq.replace('T', 'U') + get_reverse_complement(barcode_5) + primer_3
        ### reverse complementary
        # new_seq = primer_5 + get_reverse_complement(seq) + primer_3
        # new_seq = primer_5 + seq.replace('T', 'U') + primer_3
        df_data.loc[idx, 'motif'] = new_seq

    # df_data.to_csv('/share/project/UltraGen/data/wanghui/HNRNPCL1_HTR.csv', index=False)
    # df_data.to_csv('/share/project/UltraGen/data/wanghui/TIA1_RBNS.csv', index=False)
    df_data.to_csv('/share/project/UltraGen/data/processed/HTR/PRJEB25907/TARDBP_new.csv', index=False)
    # df_data.loc[:, df_data.columns.tolist()[:-1]].to_csv('/share/project/UltraGen/data/wanghui/HNRNPCL1_RBNS.csv', index=False)