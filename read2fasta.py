import os
import re
import pandas as pd
from tqdm import tqdm

opj=os.path.join
ope=os.path.exists

primer_5 = "GGAGCUCAGCCUUCACUGC".replace('U',"T")
primer_3 = "GGCACCACGGUCGGAUCCAC".replace("U","T")
fixed_middle = 'CTGCTTCGGCAG'

pattern = re.compile('^(TTCACTGC|TCACTGC|CACTGC|ACTGC|CTGC|TGC|GC|C)([ATGCN]{24,})')

# input_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left/split_reads'
# input_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left_primer/reads/right'
# input_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/RLB/RLB_right/reads'
# input_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left/reads_gt_2'
input_dir = '/share/project/UltraGen/UltraSelex'
# input_dir = '/home/zmchen/project/rfdiff/camsol'
# output_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left/reads_with_primer'
# output_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left/reads_in_fastq'
# output_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left_primer/fasta'
# output_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/RLB/RLB_right/fasta'
# output_dir = '/home/zmchen/project/pretrain/rbns_pipeline/results/NLB/NLB_left/fastq_gt_2'
output_dir = '/share/project/UltraGen/UltraSelex'
# output_dir = '/home/zmchen/project/rfdiff/camsol'

left_region = True
# add_primer = True
add_primer = False

if left_region:
    left_primer = primer_5
    right_primer = fixed_middle
else:
    left_primer = fixed_middle
    right_primer = primer_3

def read2fasta():
    for file in tqdm(os.listdir(input_dir)):
        l_write = []
        if not file.endswith('reads'):
            continue
        file_path = opj(input_dir, file)
        with open(file_path, 'r') as fi:
            lines = [line.rstrip() for line in fi.readlines()]
        
        if add_primer:
            if left_region:
                for idx in range(len(lines)):
                    l_write.append('>seq_' + str(idx) + '\n')
                    seq = lines[idx]
                    m = re.match(pattern, seq)
                    if m:
                        l_write.append(left_primer + m.group(2) + right_primer + '\n')
                    else:
                        l_write.append(left_primer + seq + right_primer + '\n')
            else:
                for idx in range(len(lines)):
                    l_write.append('>seq_' + str(idx) + '\n')
                    l_write.append(left_primer + lines[idx] + right_primer + '\n')
        else:
            for idx in range(len(lines)):
                l_write.append('>seq_' + str(idx) + '\n')
                l_write.append(lines[idx] + '\n')

        with open(opj(output_dir, file.split('.')[0] + '.fasta'), 'w') as fo:
            for line in l_write:
                fo.write(line)

def read2fastq():
    for file in tqdm(os.listdir(input_dir)):
        if ope(opj(output_dir, file.split('.')[0] + '.fastq')):
            continue

        l_write = []
        file_path = opj(input_dir, file)
        with open(file_path, 'r') as fi:
            lines = [line.rstrip() for line in fi.readlines()]
        
        if add_primer:
            if left_region:
                for idx in range(len(lines)):
                    # first row in fastq
                    l_write.append('@seq_' + str(idx) + '\n')
                    seq = lines[idx]
                    m = re.match(pattern, seq)
                    seq_FL = left_primer + m.group(2) + right_primer if m else left_primer + seq + right_primer
                    # second row in fastq
                    l_write.append(seq_FL + '\n')
                    # third row in fastq
                    l_write.append('+seq_' + str(idx) + '\n')
                    # forth row in fastq
                    l_write.append('@' * len(seq_FL) + '\n')
                    
            else:
                for idx in range(len(lines)):
                    # first row in fastq
                    l_write.append('@seq_' + str(idx) + '\n')
                    # second row in fastq
                    seq_FL = left_primer + lines[idx] + right_primer
                    l_write.append(seq_FL + '\n')
                    # third row in fastq
                    l_write.append('+seq_' + str(idx) + '\n')
                    # forth row in fastq
                    l_write.append('@' * len(seq_FL) + '\n')
        else:
            for idx in range(len(lines)):
                # first row in fastq
                l_write.append('@seq_' + str(idx) + '\n')
                # second row in fastq
                l_write.append(lines[idx] + '\n')
                # third row in fastq
                l_write.append('+seq_' + str(idx) + '\n')
                # forth row in fastq
                l_write.append('@' * len(lines[idx]) + '\n')

        with open(opj(output_dir, file.split('.')[0] + '.fastq'), 'w') as fo:
            for line in l_write:
                fo.write(line)

if __name__ == "__main__":
    read2fasta()