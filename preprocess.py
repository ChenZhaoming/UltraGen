import os
import re
import math
import difflib
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import argparse
from multiprocessing import Pool

opj=os.path.join
ope=os.path.exists

data_dir = '/home/zmchen/project/pretrain/data/conventional_selex'


d_complement = {'A': 'T', 'T':'A', 'G': 'C', 'C': 'G', 'N': 'N'}


def get_reverse_complement(seq):
    """Function to get reverse complementary of input sequence

    Args:
      seq: Input sequence string

    Returns:
      Reversed complementary string:
    """
    #return ''.join([d_complement[base] for base in seq[::-1]])
    return "".join(d_complement.get(base, base) for base in reversed(seq))


def fastq_parse(file_path):
    """Function to parse fastq file

    Args:
      file_path: Path of fastq file

    Returns:
      pd.DataFrame containing:
        * 'seq': parsed sequence
        * 'len': length of parsed sequence
        * 'cnt': number of parsed sequence
    """
#     d_parse = {}
    d_parse = defaultdict(int)
    for seq_record in SeqIO.parse(file_path, 'fastq'):
        ## get sequence from SeqIO
        seq = str(seq_record.seq)
#         if seq not in d_parse.keys():
#             d_parse[seq] = {}
#             d_parse[seq]['len'] = len(seq)
#             d_parse[seq]['cnt'] = 0
        ## the same sequence may appear many times
        d_parse[seq] += 1
    # return file in dataframe format
    output_data_frame = pd.DataFrame.from_dict(d_parse, orient='index', columns=['cnt'])
    output_data_frame['len'] = output_data_frame.index.map(lambda x: len(x))
#     return pd.DataFrame.from_dict(d_parse, orient='index', columns=['len', 'cnt'])
    return output_data_frame


def batch_load(file_dir):
    """Function to load fastq files in batch

    Args:
      file_dir: Folder path that stores the fastq files

    Returns:
      pd.DataFrame containing:
        * 'seq': parsed sequence
        * 'len': length of parsed sequence
        * 'cnt': number of parsed sequence
    """
    d_parse = {}
    l_index = []
    l_lens = []

    ## load fastq file in dataframe format and store in a dict
    for file in os.listdir(file_dir):
        d_parse[file.split('.')[0]] = fastq_parse(opj(file_dir, file))

    ## get the full set of parsed sequence
    for k, v in d_parse.items():
        l_index += v.index.tolist() 

    ## all parsed unique sequences
    l_index = list(set(l_index))
    
    ## initialize the dataframe that stores the final result
    df_merge = pd.DataFrame(np.zeros((len(l_index), len(d_parse)), dtype=np.int32), columns=list(d_parse.keys()))
    df_merge['seq'] = l_index
    df_merge = df_merge.set_index(['seq'])

    ## merge the parsed dataframe into one final dataframe
    for k, v in tqdm(d_parse.items()):
        for seq in v.index.tolist():
            df_merge.loc[seq, k] = df_merge.loc[seq, k] + v.loc[seq, 'cnt']

    ## rank the column names according the round number
    df_merge = df_merge.loc[:, sorted(df_merge.columns.tolist())]

    ## get the sequence lenth
    l_lens = [len(seq) for seq in df_merge.index.tolist()]
    df_merge['len'] = l_lens

    return df_merge


def sort_match(l_match, pattern):
    """Function to rank motifs according to the similarity against pattern

    Args:
      l_match: list that stores the matched motifs
      pattern: the sequence motif to be matched

    Returns:
      list containing ranked motifs:
    """
    l_score = []
    for match in l_match:
        ## calculate the similarity between pattern and motif
        l_score.append(difflib.SequenceMatcher(None, match, pattern).quick_ratio())

    return [item[0] for item in sorted(dict(zip(l_match, l_score)).items(), key=lambda d: d[1], reverse=True)]


def get_Pattern(pattern, seq):
    """Function to get motifs that match pattern target pattern

    Args:
      pattern: the sequence motif to be matched
      seq: the full sequence

    Returns:
      pd.DataFrame containing:
        * 'seq': parsed sequence
        * 'len': length of parsed sequence
        * 'cnt': number of parsed sequence
    """
    l_match = []
    ## change sequence motif from string format to list format
    l_pattern = [pattern[i] for i in range(len(pattern))]
    ## replace each position of target motif with [ATGCN]
    ## when the position is not in the head or in the tail
    ## replace each position of target motif with 1 to 3 [ATGCN]
    for i in range(len(l_pattern)):
        l_tmp = l_pattern.copy()
        if i == 0 or i == len(l_pattern) - 1:
            l_tmp[i] = '[ATGCN]?'
            re_pattern = re.compile(''.join(l_tmp))
            l_match = l_match + list(set(re_pattern.findall(seq)).difference(l_match))
        else:
            for j in range(4):
                l_tmp[i] = '[ATGCN]?' * (j+1)
                re_pattern = re.compile(''.join(l_tmp))
                l_match = l_match + list(set(re_pattern.findall(seq)).difference(l_match))

    ## insert [ATGCN] among consecutive residues
    l_tmp = []
    for i in range(len(pattern)):
        l_tmp.append(pattern[i])
        if i != len(pattern) - 1:
            l_tmp.append('[ATGCN]?')
    re_pattern = re.compile(''.join(l_tmp))
    
    l_match = l_match + list(set(re_pattern.findall(seq)).difference(l_match))
    ## random replce two residues with [ATGCN]
    if len(l_match) == 0:
        for i in range(len(l_pattern)):
            for j in range(len(l_pattern)):
                l_tmp = l_pattern.copy()
                if i != j:
                    l_tmp[i] = '[ATGCN]?'
                    l_tmp[j] = '[ATGCN]?'
                    re_pattern = re.compile(''.join(l_tmp))
                    l_match = l_match + list(set(re_pattern.findall(seq)).difference(l_match))
    ## random replce two residues with [ATGCN][ATGCN]
    if len(l_match) == 0:
        for i in range(len(l_pattern)):
            for j in range(len(l_pattern)):
                l_tmp = l_pattern.copy()
                if i != j:
                    l_tmp[i] = '[ATGCN]?[ATGCN]?'
                    l_tmp[j] = '[ATGCN]?[ATGCN]?'
                    re_pattern = re.compile(''.join(l_tmp))
                    l_match = l_match + list(set(re_pattern.findall(seq)).difference(l_match))

    return l_match[0] if len(l_match) > 0 else seq


def motif_extract(df_parse, term_5, term_3, midmin=-1, midmax=-1, fulllen=None, motif_lens=10, len_cutoff=0.9, reverse='reverse'):
    """
    term_5: pattern sequence in the 5' terminal
    term_3: pattern sequence in the 3' terminal
    """
    ## patterns
    head = get_reverse_complement(term_3) if reverse == 'reverse' else term_5
    tail = get_reverse_complement(term_5) if reverse == 'reverse' else term_3
    head_min = head[-motif_lens:]
    tail_min = tail[:motif_lens]

    l_pre = []
    l_suf = []
    l_mid = []
    l_mid_lens = []
    l_flag = []

    df_result = df_parse.loc[df_parse['len'] >= fulllen * len_cutoff].loc[df_parse['len'] <= fulllen * (2-len_cutoff)].copy() if fulllen is not None else df_parse.copy()

    for seq in tqdm(df_result.index.tolist()):
        prefix = '-'
        suffix = '-'
        # remain = seq
        middle = seq
        flag=0
        ###### header match
        ## full match
        if re.search(head, seq):
            prefix = seq[:seq.index(head) + len(head)]
            middle = seq[seq.index(head) + len(head):]
        ## partial match
        else:
            pattern = get_Pattern(head_min, seq)
            prefix = seq[:seq.index(pattern) + len(pattern)]
            middle = seq[seq.index(pattern) + len(pattern):]
        
        ###### tailer match
        ## full match
        if re.search(tail, middle):
            suffix = middle[middle.rindex(tail):]
            middle = middle[:middle.rindex(tail)]
        ## partial match
        else:
            pattern = get_Pattern(tail_min, middle)
            suffix = middle[middle.rindex(pattern):]
            middle = middle[:middle.rindex(pattern)]

        if prefix == head and suffix == tail and len(middle) >= midmin and len(middle) <= midmax:
            flag = 3
        elif len(middle) >= midmin and len(middle) <= midmax:
            flag = 2 if motif_lens == 10 else 1
        else:
            flag = 0

        l_pre.append(prefix)
        l_suf.append(suffix)
        l_mid.append(middle)
        l_mid_lens.append(len(middle))
        l_flag.append(flag)

    df_result['prefix'] = l_pre
    df_result['middle'] = l_mid
    df_result['midlen'] = l_mid_lens
    df_result['suffix'] = l_suf
    df_result['quality'] = l_flag

    return df_result


def multi_extract(ite, input_data, args):
    df_rst = motif_extract(input_data, term_5=args.term_5, term_3=args.term_3, midmin=args.midmin, midmax=args.midmax, motif_lens=args.motif_lens, reverse=args.reverse)
    return {ite: df_rst}

def data_preparation(processes_count, df_input):
    l_split = np.array_split(df_input, processes_count)
    d_pool = {str(i): l_split[i] for i in range(len(l_split))}
    return d_pool

def main(args):
    assert args.file_parse is not None or args.input_dir is not None
    ## parse sequences from fastq file(s)
    if args.file_parse is not None:
        df_data = batch_load(args.file_parse) if os.path.isdir(args.file_parse) else fastq_parse(args.file_parse)
    ## load parsed sequences from csv file
    else:
        df_data = pd.read_csv(args.input_dir).set_index('seq')

    if args.term_5 is not None and args.term_3 is not None and args.midmin > -1 and args.midmax > -1:
        
        if args.cpu_num == 1:
            df_rst = motif_extract(df_data, term_5=args.term_5, term_3=args.term_3, midmin=args.midmin, midmax=args.midmax, motif_lens=args.motif_lens, reverse=args.reverse)
        else:
            processes_pool = Pool(args.cpu_num)
            
            d_pool = data_preparation(args.cpu_num, df_data)

            results = [processes_pool.apply_async(multi_extract, args=(ite, input_data, args)) for ite, input_data in d_pool.items()]
            
            results_new = [p.get() for p in results]

            l_result_list = []

            for result in results_new:
                for ite in result.keys():
                    l_result_list.append(result[ite])
            
            df_rst = pd.concat(l_result_list, axis=0)
        
        df_data = df_rst.copy()
    
    if args.output_dir is not None:
        df_data.to_csv(args.output_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file_parse", type=str, default=None,
        help='''File dir for batch load'''
    )
    parser.add_argument(
        "--input_dir", type=str, default=None,
        help='''File dir for batch load'''
    )
    parser.add_argument(
        "--output_dir", type=str,
        help='''Output dir'''
    )
    parser.add_argument(
        "--term_5", type=str, default=None,
        help='''motif in 5' terminal '''
    )
    parser.add_argument(
        "--term_3", type=str, default=None,
        help='''motif in 3' terminal'''
    )
    parser.add_argument(
        "--midmin", type=int, default=-1,
        help='''lower bound of the motif length'''
    )
    parser.add_argument(
        "--midmax", type=int, default=-1,
        help='''upper bound of the motif length'''
    )
    parser.add_argument(
        "--motif_lens", type=int, default=10,
        help='''minimal length of the recognized motif in 5'- or 3' terminal'''
    )
    parser.add_argument(
        "--reverse", type=str, 
        help='''if reverse'''
    )
    parser.add_argument(
        "--cpu_num", type=int, default=1,
        help='''number of cpu core'''
    )
    args = parser.parse_args()

    main(args)
