import random
import RNA
import math
import os
import re
import itertools
import Levenshtein
import multiprocessing as mp
import numpy as np
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO


NT_complement = {'A': 'U', 'U':'A', 'G': 'C', 'C': 'G'}

NTs = [k for k,v in NT_complement.items()]

paired_NTs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

desired_ss = {'5_term': '.....(((((((((((........)))))))))))', '3_term': '(((((((((((........))))))))))).....'}


def levenshteinDistance(qSeq, sSeq):
    """ Calculate levenshtein distance between query sequence and source sequence

    Args:
        qSeq: query sequence
        sSeq: source sequence

    Return:
        levenshteinDistance
    """
    return Levenshtein.distance(qSeq, sSeq)


def gen_combination(pools, repeats):
    """ Get combinations of items in source pools

    Args:
        pools: list of string that to be combined
        repeats: number of repeats

    Return: list of string
    """
    return [''.join(combination) for combination in list(itertools.product(pools, repeat=repeats))]


def gen_random(lens):
    return ''.join([random.sample(NTs, 1)[0] for i in range(lens)])


def gen_paired_seq(seq_lens):
    """ generate paired sequences with desired length

    Args:
        seq_lens: desired sequence length

    Return:
        generated sequence, reverse complementary of generated sequence
    """
    l_left_seq = []
    l_right_seq = []
    # randomly sample a pair of NTs at a time
    for i in range(seq_lens):
        pair = random.sample(paired_NTs, 1)[0]
        l_left_seq.append(pair[0])
        l_right_seq.append(pair[1])

    return ''.join(l_left_seq), ''.join(l_right_seq)[::-1]


def gen_loop(fixed, pad_lens):
    """ generate unstructured loop with [A/C/G/U]{0-3} fixed [A/C/G/U]{3-0}

    Args:
        fixed: 'CUUGA'
        pad_lens: number of random NTs to be padded in either left or right region of fixed sequece

    Return:
        list of seqs with target_seq_lens NTs
    """
    target_seq_lens = len(fixed) + pad_lens
    # generate non-redundant combinations for NTs
    random_seqs = gen_combination(NTs, pad_lens)
    
    candidate_seqs = []

    for seq in random_seqs:
        left_lens = random.randint(0, pad_lens)
        candidate_seqs.append(seq[:left_lens] + fixed + seq[left_lens:])
    # # generate flanking sequences
    # flanking_seqs = gen_combination(random_seqs, 2)

    # candidate_seqs = []
    # for flank_seq in flanking_seqs:
    #     # [A/C/G/U]{3} fixed [A/C/G/U]{3}
    #     combined_seq = flank_seq[:pad_lens] + fixed + flank_seq[-pad_lens:]
    #     # slice to get sequence with desired length
    #     sub_seqs = [combined_seq[i:i+target_seq_lens] for i in range(pad_lens + 1)]
    #     candidate_seqs.extend(sub_seqs)

    candidate_with_desired_ss = []
    # predict secondary structure for each candidate sequence
    for seq in list(set(candidate_seqs)):
        # sequence without any secondary structure is expected
        if RNA.fold(seq)[0] == '.' * target_seq_lens:
            candidate_with_desired_ss.append(seq)

    return candidate_with_desired_ss


def pad_seqs(seq):
    """ padding seq to satisfy length requirement
    """
    # total seq length
    lens = random.randint(50, 60)
    # left primer length
    left_lens = random.randint(0, lens-len(seq))
    # right primer length
    right_lens = lens-len(seq)-left_lens

    return gen_random(left_lens) + seq + gen_random(right_lens)


def rnafold(seq):
    """ get ss in dot-bracket format
    """
    return RNA.fold(seq)[0]


def filter_random_seqs(items):
    """ filter random generated seqs with secondary structure
    """
    
    seq = items[0]
    ss = items[1]

    re_pattern = re.compile("([\(]{4,})([\.]{5,})([\)]{4,})")

    iter = re.finditer(re_pattern, ss)

    for it in iter:
        i, j = it.span()
        # idxs of loop
        if 'CUUGA' in seq[i+len(it.group(1)):j-len(it.group(3))]:
            return '1', i, j
    return '0', 0, 0


def check_ss_valid(sec_struc):
    """ only one type of brackets expected. e.g. ..(... or ...)..

    Args:
        sec_struc: secondary structure in dot-bracket format

    Return:
        True if input ss is valid else False
    """
    items = list(set([char for char in sec_struc]))
    if ')' in items and '(' in items:
        return False
    if items == ['.']:
        return False
    return True


def mutation_event(seq):
    sourceNTs = [NT for NT in seq]
    # if exists cuuga motif, 70% probability to mutate this motif
    if 'CUUGA' in seq and random.random() > 0.20:
        idx = random.randint(seq.index('CUUGA'), seq.index('CUUGA')+4)
    else:
        idx = random.randint(0, len(sourceNTs)-1)
    
    event = random.randint(1, 3)
    # insert
    if event == 1:
        sourceNTs[idx] = sourceNTs[idx] + random.sample(NTs, 1)[0]
    # delete
    if event == 2:
        sourceNTs[idx] = ''
    # revise
    if event == 3:
        sourceNTs[idx] = random.sample(list(set(NTs).difference([sourceNTs[idx]])), 1)[0]
    return ''.join(sourceNTs)


def random_mutate_seq(seq):
    """ random mutate 1-5 NT(s) within the given sequence

    Args:
        seq: source sequence

    Return:
        mutated seq
    """
    numMuts = random.randint(1, 5)

    for i in range(numMuts):
        seq = mutation_event(seq)

    return seq


def seqPartition(seqs, sec_strucs):
    # 1,2,3 represent one positive set and two negative sets, in which group2 are seqs without CUUGA, and group3 are seqs with paired CUUGA
    group1, group2, group3 = [], [], []

    for idx in range(len(seqs)):
        seq = seqs[idx]
        ss = sec_strucs[idx]
        if 'CUUGA' not in seq:
            group2.append(seq)
        else:
            # may be more than one CUUGA motif
            iter = re.finditer(re.compile("CUUGA"), seq)
            for it in iter:
                i, j = it.span()
                if ss[i:j].count('.') < 4:
                    # group3
                    group3.append(seq)
                    break
            if len(group1) + len(group2) + len(group3) < idx+1:
                group1.append(seq)
    
    return group1, group2, group3


def mutate_seq(seq, dist):
    """ get mutated seqs with preset levenshteinDistance

    Args:
        seq: source sequence
        dist: desired levenshtein distance between source sequence and mutated sequence

    Return:
        list of all the possible sequences with desired levenshtein distance against source sequence
    """
    l_seqs = []
    # generate all the possible sequences with desired length
    candidates = gen_combination(NTs, len(seq))
    # measure levenshteinDistance
    for candidate in candidates:
        if levenshteinDistance(candidate, seq) == dist:
            l_seqs.append(candidate)
    return l_seqs


def get_single_mutant(seq, flanking_lens, loop_idx):
    """ loop_idx is 16 when 5-NT random primer is in 5' terminal and is 11 when 5-NT random primer is in 3' terminal

    Args:
        seq: source sequence
        flanking_lens: regions to be mutated
        loop_idx: 16 or 11 according to the position of 5-NT random region

    Return:
        list of single-NT-mutated sequences
    """
    loop_lens = 8
    # left flanking sequence
    left_seq = seq[loop_idx-flanking_lens:loop_idx]
    # right flanking sequence
    right_seq = seq[loop_idx+loop_lens:loop_idx+loop_lens+flanking_lens]

    l_seqs = []
    # mutate one NT in left flanking region
    left_candidates = mutate_seq(left_seq, 1)
    # mutate one NT in right flanking region
    right_candidates = mutate_seq(right_seq, 1)

    for candidate in left_candidates:
        l_seqs.append(seq[0:loop_idx-flanking_lens] + candidate + seq[loop_idx:])

    for candidate in right_candidates:
        l_seqs.append(seq[0:loop_idx+loop_lens] + candidate + seq[loop_idx+loop_lens+flanking_lens:])

    return l_seqs


def get_best_mutant(seq, mutants, loop_idx):
    """ keep the best mutant with the highest levenshteinDistance

    Args:
        seq: source sequence
        mutants: list of mutants, current one mutation only
        loop_idx: 16 or 11 according to the position of 5-NT random region

    Return:
        soruce sequence, mutated sequence
    """
    loop_lens = 8
    ls_dist = 0
    best_mutant = seq
    keyMotif = 'CUUG'
    gt_ss = RNA.fold(seq)[0]
    best_ss = gt_ss
    for mutant in mutants:
        ss, energy = RNA.fold(mutant)
        offset = seq[loop_idx:loop_idx+loop_lens].index(keyMotif)
        # both loop and motif region are unpaired in the origin
        loop_ss = ss[loop_idx:loop_idx+loop_lens]
        motif_ss = ss[loop_idx+offset:loop_idx+offset+len(keyMotif)]
        # calculate levenshteinDistance
        dist = levenshteinDistance(gt_ss, ss)
        # mutant with desired secondary structure and lower energy
        if check_ss_valid(motif_ss) and check_ss_valid(loop_ss) and dist > ls_dist:
            best_mutant = mutant
            best_ss = ss
            ls_dist = dist

    return seq, best_mutant, gt_ss, best_ss, ls_dist


def get_mutant(seq):
    """ function for multiprocessing
    """
    ## processing seqs with random region in 5' terminal
    loop_idx = 16
    ## processing seqs with random region in 3' terminal
    # loop_idx = 11
    l_single_mutant = get_single_mutant(seq, flanking_lens=5, loop_idx=loop_idx)
    return get_best_mutant(seq, l_single_mutant, loop_idx=loop_idx)


def seq2fasta(seqs, save_dir):
    """ save a batch of sequences into a fasta file
    """
    with open(save_dir, 'w') as fp:
        for i in range(len(seqs)):
            fp.write('>seq_' + str(i) + '\n')
            fp.write(seqs[i] + '\n')


def seq2df(results):
    """ change result into dataframe
    """
    # calculate Levenshtein distance for each row
    # dists = np.array([Levenshtein.distance(item[2], item[3]) for item in np.array(results)]).reshape(-1, 1)
    # change results into dataframe
    df = pd.DataFrame(results, columns=['Origin seq', 'Mutated Seq', 'Origin ss', 'Mutated ss', 'Levenshtein'])
    df = df.sort_values(by='Levenshtein', ascending=False, inplace=False)

    return df.loc[df['Levenshtein'] > 0, ['Origin seq', 'Mutated Seq', 'Mutated ss', 'Levenshtein']]


if __name__ == "__main__":
    # ########################################### first version of generated seqs ###########################################
    # # number of seqs to be generated
    # num_seqs = 10000000
    # # specify ss type,
    # ssType = '5_term'
    # # ssType = '3_term'
    # ## 5_term for .....(((((((((((........))))))))))) and 3_term for (((((((((((........))))))))))).....
    # ssTemplate = desired_ss[ssType]

    # print('prepare samping pools for each components')
    # ## 256 candidates
    # l_loop_candidate = gen_loop(fixed='CUUGA', pad_lens=3)
    # ## 1024 candidates
    # l_primer_candidate = gen_combination(pools=NTs, repeats=5)

    # l_pair_candidate = []
    # # the function gen_combination is not utilized as there will be too many seqs when desired sequence length is long
    # for i in tqdm(range(num_seqs)):
    #     l_pair_candidate.append(gen_paired_seq(seq_lens=11))

    # print('generate sequences with {} primer, {} pair, and {} loop candidates'.format(len(l_primer_candidate), len(l_pair_candidate), len(l_loop_candidate)))
    # # sample primers
    # primers = np.char.array(np.random.choice(l_primer_candidate, size=num_seqs, replace=True))
    # # sample pairs
    # pairs = random.sample(l_pair_candidate, num_seqs)
    # left_pairs = np.char.array(np.array(pairs)[:, 0])
    # right_pairs = np.char.array(np.array(pairs)[:, 1])
    # # sample loops
    # loops = np.char.array(np.random.choice(l_loop_candidate, size=num_seqs, replace=True))
    # # concatenate components to get full seqs
    # seqPools = list(primers + left_pairs + loops + right_pairs) if ssType == '5_term' else list(left_pairs + loops + right_pairs + primers)
    # # predict ss for generated seqs
    # print('run ss prediction')
    # with mp.Pool(100) as p:
    #     ssPools = list(tqdm(p.imap(rnafold, seqPools), total=len(seqPools)))

    # assert len(seqPools) == len(ssPools)
    # # screen seqs with desired ss
    # l_seq_candidate = []
    # for i in range(len(seqPools)):
    #     if ssPools[i] == ssTemplate:
    #         l_seq_candidate.append(seqPools[i])
    # ## save generated seqs into a fasta file
    # seq2fasta(l_seq_candidate, '/share/project/UltraGen/data/raw/random_generation/CUUGA/stem_loop_CUUGA_primer3_0926_v2.fasta')

    ########################################### second version of generated seqs ###########################################
    # number of seqs to be generated
    num_seqs = 1000000
    ## loop candidates
    print('prepare candidates')
    l_loop_candidate = []
    for pad in tqdm(range(4, 8)):
        l_loop_candidate.extend(gen_loop(fixed='CUUGA', pad_lens=pad))
    # num_seqs candidates
    l_pair_candidate = []
    # the function gen_combination is not utilized as there will be too many seqs when desired sequence length is long
    for lens in tqdm(range(6, 16)):
        for repeats in range(int(num_seqs / 10)):
            l_pair_candidate.append(gen_paired_seq(seq_lens=lens))
    # sample pairs
    pairs = random.sample(l_pair_candidate, num_seqs)
    left_pairs = np.char.array(np.array(pairs)[:, 0])
    right_pairs = np.char.array(np.array(pairs)[:, 1])
    # sample loops
    loops = np.char.array(np.random.choice(l_loop_candidate, size=num_seqs, replace=True))
    # seqs without primers
    seqPools = list(left_pairs + loops + right_pairs)
    # seqs with primers
    print('padding seqs')
    with mp.Pool(100) as p:
        seqPools = list(tqdm(p.imap(pad_seqs, seqPools), total=len(seqPools)))
    # run ss prediction
    print('run ss prediction')
    with mp.Pool(100) as p:
        ssPools = list(tqdm(p.imap(rnafold, seqPools), total=len(seqPools)))
    
    assert len(seqPools) == len(ssPools)
    # screen seqs with desired ss
    print('filter seqs')
    with mp.Pool(100) as p:
        filters = list(tqdm(p.imap(filter_random_seqs, np.vstack([seqPools, ssPools]).T), total=len(seqPools)))
    
    df_result = pd.DataFrame(np.hstack([np.vstack([seqPools, ssPools]).T, np.array(filters)]), columns=['seq', 'ss', 'filter', 'loop_start', 'loop_end'])
    ## save generated seqs into a fasta file
    # df_result.to_csv('/share/project/UltraGen/data/raw/random_generation/CUUGA/origin_v2.csv', index=False)
    df_result.groupby(['filter']).get_group('1').loc[:, ['seq', 'ss', 'loop_start', 'loop_end']].to_csv('/share/project/UltraGen/data/raw/random_generation/CUUGA/origin_v3.csv', index=False)


    # NOTE: sequences are pre-generated and can be load directly
    # target_seqs = [str(seq_record.seq) for seq_record in SeqIO.parse('/share/project/UltraGen/data/raw/random_generation/CUUGA/stem_loop_CUUGA_5term.fasta', 'fasta')]
    # # target_seqs = [str(seq_record.seq) for seq_record in SeqIO.parse('/share/project/UltraGen/data/raw/random_generation/CUUGA/stem_loop_CUUGA_3term.fasta', 'fasta')]

    # with mp.Pool(100) as p:
    #     results = list(tqdm(p.imap(get_mutant, target_seqs), total=len(target_seqs)))
    # # change results to dataframe
    # df_result = seq2df(results)
    # print('Saving generated {} seqs'.format(len(df_result)))
    # df_result.to_csv('/share/project/UltraGen/data/raw/random_generation/CUUGA/CUUG_5term.csv', index=False)

    df_data = pd.read_csv('/share/project/UltraGen/data/raw/random_generation/CUUGA/origin_v3.csv')
    mutated_seqs = []
    # random select original seqs to be mutated 500 times, which generate 50w seqs
    seqPool1 = random.sample(df_data['seq'].tolist(), 1000)
    for i in tqdm(range(500)):
        with mp.Pool(100) as p:
            mutated_seqs.extend(list(p.imap(random_mutate_seq, seqPool1)))
    # random select original seqs to be mutated 2 time, which generate 50w seqs
    seqPool2 = random.sample(list(set(df_data['seq']).difference(seqPool1)), 50000)
    for i in tqdm(range(10)):
        with mp.Pool(100) as p:
            mutated_seqs.extend(list(p.imap(random_mutate_seq, seqPool2)))

    print('run ss prediction')
    with mp.Pool(100) as p:
        mutated_ss = list(tqdm(p.imap(rnafold, mutated_seqs), total=len(mutated_seqs)))

    group1, group2, group3 = seqPartition(mutated_seqs, mutated_ss)

    print(len(group1), len(group2), len(group3))

    assert len(group1) + len(group2) + len(group3) == len(mutated_seqs)

    pos_seqs = group1 + seqPool1 + seqPool2
    neg_seqs = group2 + group3

    print(len(pos_seqs), len(neg_seqs))

    seq2fasta(pos_seqs, '/share/project/UltraGen/data/raw/random_generation/CUUGA/pos_v4.fasta')
    seq2fasta(neg_seqs, '/share/project/UltraGen/data/raw/random_generation/CUUGA/neg_v4.fasta')