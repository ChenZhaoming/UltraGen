import os
import re
import math
import difflib
import datetime
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
        ## the same sequence may appear many times
        d_parse[seq] += 1
    # return file in dataframe format
    output_data_frame = pd.DataFrame.from_dict(d_parse, orient='index', columns=['cnt'])
    output_data_frame['len'] = output_data_frame.index.map(lambda x: len(x))
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


# def sort_match(l_match, pattern):
#     """Function to rank motifs according to the similarity against pattern

#     Args:
#       l_match: list that stores the matched motifs
#       pattern: the sequence motif to be matched

#     Returns:
#       list containing ranked motifs:
#     """
#     l_score = []
#     for match in l_match:
#         ## calculate the similarity between pattern and motif
#         l_score.append(difflib.SequenceMatcher(None, match, pattern).quick_ratio())

#     return [item[0] for item in sorted(dict(zip(l_match, l_score)).items(), key=lambda d: d[1], reverse=True)]


def get_optimal_match(pattern, seq, midmin, midmax, end):
    """Function to get motifs that match pattern best

    Args:
      pattern: the sequence motif to be matched
      seq: the full sequence
      midmin: the minimal length requirement of random_region
      midmax: the maximal length requirement of random_region
      end: head or tail

    Returns:
      best matched sequence
    """
    ## for primer_A, find longest match, while for primer_B, find match that meets the length requirement of random_region
    # match_pattern_end = ",})" if end == 'head' else ","+str(midmax)+"})([ATGCN]{0,})"
    match_pattern_end = ",})" if end == 'head' else ","+str(midmax)+"})"
    ## full match
    # re_pattern = re.compile(r"([ATGCN]{0,})("+pattern+")([ATGCN]{"+str(midmin)+","+str(midmax)+"})([ATGCN]{0,})")
    re_pattern = re.compile(r"([ATGCN]{0,})("+pattern+")([ATGCN]{"+str(midmin)+match_pattern_end)
    matchs = re.search(re_pattern, seq)
    if matchs:
        primer = matchs.group(1) + matchs.group(2)
        # randme = matchs.group(3) if end =='head' else matchs.group(3) + matchs.group(4)
        randme = matchs.group(3)
        return primer, randme
    ## two situations that causue unmatch:
    ## 1. the length of target sequence does not satisfy the minimal requirement of full match. (len(pattern) + midmin > len(seq))
    ##    this may due to truncation of pattern sequence
    ## 2. there are some mistake in pattern
    ## problem 1 can be solved by passing a subsequence of pattern into this function

    #### solve the problem 2
    ## replace each residue in pattern with [ATGCN]{0,2}, e.g. ATG -> A[ATGCN]{0,2}G (ATG -> ACG, ATG -> AG, ATG -> ATCG)
    ## matching the replacement, insertion, and deletion incidence of a single residue
    # pattern_mismatch = list(set([pattern[:i]+'[ATGCN]{0,2}'+pattern[i+1:] for i in range(len(pattern))]))
    pattern_mismatch = [pattern[:i]+'[ATGCN]{0,2}'+pattern[i+1:] for i in range(len(pattern))]
    for ptn in pattern_mismatch:
        # re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+",})")
        # re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+","+str(midmax)+"})([ATGCN]{0,})")
        re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+match_pattern_end)
        matchs = re.search(re_pattern, seq)
        if matchs:
            # return matchs.group(1) + matchs.group(2), matchs.group(3) if end =='head' else matchs.group(1) + matchs.group(2), matchs.group(3) + matchs.group(4)
            primer = matchs.group(1) + matchs.group(2)
            # randme = matchs.group(3) if end =='head' else matchs.group(3) + matchs.group(4)
            randme = matchs.group(3)
            return primer, randme
    # ## remove one residue in pattern, e.g. ATG -> AG
    # pattern_mismatch = list(set([pattern[:i]+pattern[i+1:] for i in range(len(pattern))]))
    # for ptn in pattern_mismatch:
    #     re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+",})")
    #     matchs = re.search(re_pattern, seq)
    #     if matchs:
    #         return matchs.group(1) + matchs.group(2), matchs.group(3)
    # ## insert a residue in pattern with ATGCN, e.g. ATG -> A[ATGCN]TG
    # pattern_mismatch = list(set([pattern[:i]+res+pattern[i:] for res in 'ATGCN' for i in range(len(pattern) + 1)]))
    # for ptn in pattern_mismatch:
    #     re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+",})")
    #     matchs = re.search(re_pattern, seq)
    #     if matchs:
    #         return matchs.group(1) + matchs.group(2), matchs.group(3)
    ## random replace two residues in pattern with [ATGCN]{0,2}, e.g. ATG -> [ATGCN]{0,2}C[ATGCN]{0,2}
    # pattern_mismatch = list(set([pattern[:j]+'[ATGCN]{0,2}'+pattern[j+1:i]+'[ATGCN]{0,2}'+pattern[i+1:] for i in range(1,len(pattern)-1) for j in range(1,i)]))
    pattern_mismatch = [pattern[:j]+'[ATGCN]{0,2}'+pattern[j+1:i]+'[ATGCN]{0,2}'+pattern[i+1:] for i in range(1,len(pattern)-1) for j in range(1,i)]
    for ptn in pattern_mismatch:
        # re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+",})")
        # re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+","+str(midmax)+"})([ATGCN]{0,})")
        re_pattern = re.compile(r"([ATGCN]{0,})("+ptn+")([ATGCN]{"+str(midmin)+match_pattern_end)
        matchs = re.search(re_pattern, seq)
        if matchs:
            # return matchs.group(1) + matchs.group(2), matchs.group(3) if end =='head' else matchs.group(1) + matchs.group(2), matchs.group(3) + matchs.group(4)
            primer = matchs.group(1) + matchs.group(2)
            # randme = matchs.group(3) if end =='head' else matchs.group(3) + matchs.group(4)
            randme = matchs.group(3)
            return primer, randme
    return False


def get_one_end_match(pattern, seq, midmin, midmax, cutoffmin, cutoffmax, end):
    """Function to get motifs that match pattern best

    Args:
      pattern: the sequence motif to be matched
      seq: the full sequence
      midmin: the minimal length requirement of random_region
      midmax: the maximal length requirement of random_region
      cutoffmin: the minimal length requirement of subsequence match
      cutoffmax: the maximal length requirement of subsequence match
      end: head or tail

    Returns:
      best matched sequence
    """
    matchs = False
    ## if the length of target sequence does not satisfy the minimal truncation length of pattern sequence
    if cutoffmin + midmin > len(seq):
        return matchs
    ## try to get optimal match for pattern sequence
    if len(pattern) + midmin <= len(seq):
        matchs = get_optimal_match(pattern, seq, midmin, midmax, end)
    ## seq length do not meet the requirement, or cannot get optimal match for pattern
    if not matchs or len(pattern) + midmin > len(seq):
        ## try to get match for the subsequence of pattern
        for cutoff in range(cutoffmax, cutoffmin, -1):
            sub_pattern = pattern[-cutoff:]
            matchs = get_optimal_match(sub_pattern, seq, midmin, midmax, end)
            if matchs:
                break
    ## matched motifs or False
    return matchs


def get_two_end_match(patternA, patternB, seq, midmin, midmax):
    """Function to get motifs that match pattern best

    Args:
      patternA: the sequence of primer_A to be matched
      patternB: the sequence of primer_B to be matched
      seq: the full sequence
      midmin: the minimal length requirement of random_region
      midmax: the maximal length requirement of random_region

    Returns:
      best matched sequence
    """
    ## full match in two ends, assume that a target sequence is consists of the following parts:
    ## system_seq，primerA，random_region，primerB，system_seq
    re_pattern = re.compile(r"([ATGCN]{0,})("+patternA+")([ATGCN]{"+str(midmin)+","+str(midmax)+"})("+patternB+")([ATGCN]{0,})")
    matchs = re.search(re_pattern, seq)
    if matchs:
        return matchs.group(1) + matchs.group(2), matchs.group(3), matchs.group(4) + matchs.group(5)
    
    ## random replace one residue in patternA or patternB with [ATGCN]{0,2}
    pattern_mismatch = ["("+patternA[:i]+"[ATGCN]{0,2}"+patternA[i+1:]+")([ATGCN]{"+str(midmin)+","+str(midmax)+"})("+patternB[:j]+"[ATGCN]{0,2}"+patternB[j+1:]+")" for i in range(len(patternA)) for j in range(len(patternB))]
    for ptn in pattern_mismatch:
        # re_pattern = re.compile(r"([ATGCN]{0,})"+ptn+"([ATGCN]{0,})")
        re_pattern = re.compile(r"([ATGCN]{0,})"+ptn)
        matchs = re.search(re_pattern, seq)
        if matchs:
            return matchs.group(1) + matchs.group(2), matchs.group(3), matchs.group(4)
            # return matchs.group(1) + matchs.group(2), matchs.group(3), matchs.group(4) + matchs.group(5)
    ## matched motifs or False
    return False


def motif_extract(df_parse, term_5, term_3, midmin=-1, midmax=-1, fulllen=None, motif_lens=10, len_cutoff=0.9, reverse='reverse'):
    """
    term_5: pattern sequence in the 5' terminal
    term_3: pattern sequence in the 3' terminal
    """
    ## patterns
    head = get_reverse_complement(term_3) if reverse == 'reverse' else term_5
    tail = get_reverse_complement(term_5) if reverse == 'reverse' else term_3

    l_pre = []
    l_suf = []
    l_mid = []
    l_mid_lens = []
    l_quality = []

    # df_parse = df_parse.loc[df_parse['len'] >= fulllen * len_cutoff].loc[df_parse['len'] <= fulllen * (2-len_cutoff)].copy() if fulllen is not None else df_parse.copy()

    for seq in tqdm(df_parse.index.tolist()):
        ## quality flag
        quality = 0

        ## try to find matches in two ends
        if len(seq) >= len(head) + midmin + len(tail):
            matchs = get_two_end_match(head, tail, seq, midmin, midmax)
            if matchs:
                ## system_seq + primerA
                prefix = matchs[0]
                ## random_region
                middle = matchs[1]
                ## primerB + system_seq
                suffix = matchs[2]
                ## best quality
                quality = 2 if head in prefix and tail in suffix else 1

        ## the length of target sequence does not meet the requirement, or can not find match in both ends
        if quality == 0:
            ## find match in primer_A, the minimal length of recognized motif is 8
            matchA = get_one_end_match(head, seq, midmin, midmax, cutoffmin=8, cutoffmax=motif_lens, end='head')
            prefix, middle = matchA if matchA else ('-', seq)
            ## find match in primer_B, the minimal length of recognized motif is 6
            matchB = get_one_end_match(tail[::-1], middle[::-1], midmin, midmax, cutoffmin=6, cutoffmax=motif_lens, end='tail')
            suffix, middle = (matchB[0][::-1], matchB[1][::-1]) if matchB else ('-', middle)
            ## quality is assigned 1 if both primer_A and primer_B were matched and random_motif meets the length requirement
            quality = 1 if matchA and matchB and len(middle) >= midmin and len(middle) <= midmax else 0

        l_pre.append(prefix)
        l_suf.append(suffix)
        l_mid.append(middle)
        l_mid_lens.append(len(middle))
        l_quality.append(quality)

    df_parse['prefix'] = l_pre
    df_parse['middle'] = l_mid
    df_parse['midlen'] = l_mid_lens
    df_parse['suffix'] = l_suf
    df_parse['quality'] = l_quality

    return df_parse


def multi_extract(ite, input_data, args):
    df_rst = motif_extract(input_data, term_5=args.term_5, term_3=args.term_3, midmin=args.midmin, midmax=args.midmax, motif_lens=args.motif_lens, reverse=args.reverse)
    return {ite: df_rst}

def data_preparation(processes_count, df_input):
    l_split = np.array_split(df_input, processes_count)
    d_pool = {str(i): l_split[i] for i in range(len(l_split))}
    return d_pool

def main(args):
    ## time record
    time_old = datetime.datetime.now()

    assert args.file_parse is not None or args.input_dir is not None
    ## parse sequences from fastq file(s)
    if args.file_parse is not None:
        df_data = batch_load(args.file_parse) if os.path.isdir(args.file_parse) else fastq_parse(args.file_parse)
    ## load parsed sequences from csv file
    else:
        df_data = pd.read_csv(args.input_dir).set_index('seq')
    ## run motif extraction
    if args.term_5 is not None and args.term_3 is not None and args.midmin > -1 and args.midmax > -1:
        ## run in single cpu
        if args.cpu_num == 1:
            df_rst = motif_extract(df_data, term_5=args.term_5, term_3=args.term_3, midmin=args.midmin, midmax=args.midmax, motif_lens=args.motif_lens, reverse=args.reverse)
        ## run in multiple cpus
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
    
    ## save file
    if args.output_dir is not None:
        df_data.to_csv(args.output_dir)

    ## time record
    time_new = datetime.datetime.now()
    time_diff = time_new - time_old

    print("{} of sequence to be processed. {} CPUs were used. {} minutes were passed.".format(df_data.shape[0], args.cpu_num, time_diff.total_seconds() / 60))

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