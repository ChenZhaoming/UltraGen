import os
import requests
import pandas as pd
from tqdm import tqdm
opj=os.path.join
ope=os.path.exists

save_dir = '/home/zmchen/project/pretrain/data/SRP041098'

# source_url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc=SRR1284112'
source_url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc='

def request_on_url(URL, out_path):
	r = requests.get(URL)
	with open(out_path, "wb") as code:
		code.write(r.content)

def batch_download(source_url, file_nb, save_dir):
	for i in tqdm(range(file_nb)):
		url = source_url[:-3] + str(int(source_url[-3:]) + i)
		path = opj(save_dir, url.split('=')[-1] + '.fastq.gz')
		if ope(path):
			continue
		request_on_url(url, path)

def batch_download_excel(target_files, save_dir):
	for file in tqdm(target_files):
		url = source_url + file
		path = opj(save_dir, file + '.fastq.gz')
		if ope(path):
			continue
		request_on_url(url, path)

if __name__ == '__main__':
	# batch_download(source_url, 8, save_dir)
	df_data = pd.read_excel('/home/zmchen/project/pretrain/data/SRP041098/accession_ID.xlsx', sheet_name='Sheet1')
	targets = df_data['access_ID'].tolist()
	batch_download_excel(targets, save_dir)
