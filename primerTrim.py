import pandas as pd
import numpy as np
import Levenshtein
import argparse
import multiprocessing
import re
from tqdm import tqdm

# 编辑距离
def levenshteinDistance(qSeq, sSeq = 'CTGCTTCGGCAG'):
  return Levenshtein.distance(qSeq, sSeq)

def splitSeq(rawSeq, getSeq=True, startPos=25, endPos=42, maximumChange=2, subjectSeq ='CTGCTTCGGCAG'):
  PrimerA = 'GGAGCUCAGCCUUCACUGC'.replace('U','T') #'GTGGATCCGACCGTGGTGCC'
  PrimerB = 'GGCACCACGGUCGGAUCCAC'.replace('U','T') #'GCAGTGAAGGCTGAGCTCC'
  inputSeq = rawSeq[startPos:endPos]
  Score12min = 99
  for i in range(0, len(inputSeq)-12+1):
    if inputSeq[i:i+12] == subjectSeq:
      startIndex= startPos+i
      endIndex= startPos+i+12
      Score12min = 0
      if getSeq:
        return [rawSeq[:startIndex],rawSeq[startIndex:endIndex],rawSeq[endIndex:]]#,Score12min]
      else:
        return True
    else:
      temScore = levenshteinDistance(inputSeq[i:i+12]) 
      if temScore < Score12min:
        Score12min = temScore
        startIndex= startPos+i
        endIndex= startPos+i+12
      if not getSeq and temScore<=2:
        return True
  if not getSeq:
    return False
  if Score12min ==1:                                
    for j in range(0, len(inputSeq)-13+1):
      if levenshteinDistance(inputSeq[j:j+13]) == 1:
        startIndex = startPos+j
        endIndex = startPos+j+13
        Score12min = 10
        break
    return [rawSeq[:startIndex],rawSeq[startIndex:endIndex],rawSeq[endIndex:]]#,Score12min]
  elif Score12min ==2:
    for k in range(0, len(inputSeq)-14+1):
      if levenshteinDistance(inputSeq[k:k+14]) == 1:
        startIndex = startPos+k
        endIndex = startPos+k+14
        Score12min = 20
        break
    return [rawSeq[:startIndex],rawSeq[startIndex:endIndex],rawSeq[endIndex:]]#,Score12min]
  else:
    return False


def trimLeader(mergeSeqInfo):
  #leftSeq, middleSeq, rightSeq= re.sub('[()\']','',(mergeSeqInfo)).split(',')
  leftSeq, middleSeq, rightSeq = mergeSeqInfo #.split('%')
  pattern = re.compile('^(CTGC|TGC|GC|C)([ATGCN]{24,})')
  m = re.match(pattern, leftSeq)
  if m:
    numTrim=len(m.group(1))
    if numTrim==1 and len(leftSeq)>=26+numTrim:
      newleftSeq=m.group(2)
    elif numTrim==2 and len(leftSeq)>=26+numTrim:
      newleftSeq=m.group(2)
    elif numTrim==3 and len(leftSeq)>=26+numTrim:
      newleftSeq=m.group(2)
    elif numTrim==4 and len(leftSeq)>=26+numTrim:
      newleftSeq=m.group(2)
    else:
      newleftSeq=leftSeq
  else:
    newleftSeq=leftSeq
  # return newleftSeq+middleSeq+rightSeq
  return newleftSeq, middleSeq, rightSeq


def mergeSeq(mergeSeqInfo):
  #leftSeq, middleSeq, rightSeq= re.sub('[()\']','',(mergeSeqInfo)).split(',')
  leftSeq, middleSeq, rightSeq = mergeSeqInfo #.split('%')
  return leftSeq+rightSeq

def leftSeq(mergeSeqInfo):
  #leftSeq, middleSeq, rightSeq= re.sub('[()\']','',(mergeSeqInfo)).split(',')
  leftSeq, middleSeq, rightSeq = mergeSeqInfo #.split('%')
  return leftSeq

def rightSeq(mergeSeqInfo):
  #leftSeq, middleSeq, rightSeq= re.sub('[()\']','',(mergeSeqInfo)).split(',')
  leftSeq, middleSeq, rightSeq = mergeSeqInfo #.split('%')
  return rightSeq


def main(args):
  if args.input_dir.endswith('csv'):
    rawData = pd.read_csv(args.input_dir)
  if args.input_dir.endswith('csv.gz'):
    rawData = pd.read_csv(args.input_dir, compression='gzip')
  # if args.input_dir.endswith('reads'):
  #   with open(args.input_dir, 'r') as fi:
  #     lines = [line.split()[0] for line in fi.readlines()]
  #   rawData = pd.DataFrame(lines, columns=['Sequence'])

  ncpus = multiprocessing.cpu_count()-1
  p = multiprocessing.Pool(ncpus)
  rawData['split_seq'] = list(tqdm(p.imap(splitSeq, rawData['Sequence']), total=len(rawData)))
  p.close()

  judgeCondition2= rawData['split_seq'] == False
  print ("{} % sequence will be filtered due to the missing bridge region N12!".format(np.sum(judgeCondition2)/len(judgeCondition2)*100))
  rawData=rawData[~judgeCondition2]

  tqdm.pandas()
  rawData['split_seq'] = rawData['split_seq'].progress_apply(trimLeader)

  judgeCondition4= rawData['split_seq'].apply(lambda x: len(x[0])>=25 and len(x[1])==12 and len(x[2])>=25)
  print ("{} % sequence will be filtered due to the mismatch in N26+N12+N26!".format((1-np.sum(judgeCondition4)/len(judgeCondition4))*100))
  rawData=rawData[judgeCondition4]

  rawData['merged_seq'] = rawData['split_seq'].progress_apply(mergeSeq)
  rawData['left_seq'] = rawData['split_seq'].progress_apply(leftSeq)
  rawData['right_seq'] = rawData['split_seq'].progress_apply(rightSeq)

  if args.input_dir.endswith('csv'):
    rawData.to_csv(args.output_dir, index=False)
  if args.input_dir.endswith('csv.gz'):
    rawData.to_csv(args.output_dir, compression='gzip', index=False)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument(
      "--input_dir", type=str, default=None,
      help='''Input file dir'''
  )
  parser.add_argument(
      "--output_dir", type=str, default=None,
      help='''Output file dir'''
  )
  args = parser.parse_args()

  main(args)