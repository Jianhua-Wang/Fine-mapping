#! /usr/bin/env python3

import pandas as pd
import numpy as np
from multiprocessing import Pool
from subprocess import call,check_output
import os,sys

def extract_overlap_and_cal_ld(pop,chr_id, start, stop):
	out_prefix = '{}/{}_{}_{}_{}'.format(out_dir,prefix,chr_id,start,stop)
	raw = pd.read_csv('{}.txt'.format(out_prefix),sep='\t')
	ref = pd.read_csv('{}/{}/{}_{}_{}_{}.txt.gz'.format(reference_dir,pop,pop,chr_id,start,stop),sep='\t')

	ref['snp'] = ref['CHROM'].astype(str)+':'+ref['POS'].astype(str)+':'+ref['REF']+':'+ref['ALT']
	raw['positive'] = raw['CHR'].astype(str)+':'+raw['BP'].astype(str)+':'+raw['NEA']+':'+raw['EA']
	raw['negative'] = raw['CHR'].astype(str)+':'+raw['BP'].astype(str)+':'+raw['EA']+':'+raw['NEA']

	negative = raw.merge(ref,left_on='negative',right_on='snp',how='inner')
	positive = raw.merge(ref,left_on='positive',right_on='snp',how='inner')

	negative['Zscore'] = -negative['Zscore']

	processed = pd.concat([positive,negative])

	if pd.isnull(processed.iloc[0,3]):
		processed['MAF'] = processed['{}_MAF'.format(pop)]
	processed['rsID'] = processed['ID']
	processed = processed.drop_duplicates('rsID')
	processed = processed.sort_values('BP')
	processed['coding'] = 1
	processed[['CHR', 'BP', 'rsID','MAF', 'EA', 'NEA', 'BETA', 'SE', 'P', 'Zscore','index']].to_csv('{}.processed'.format(out_prefix),sep=' ',index=False)
	paintor = processed[['CHR', 'BP', 'rsID','MAF', 'EA', 'NEA', 'BETA', 'SE', 'P', 'Zscore','index']]
	# paintor['PAINTOR'] = -1
	# paintor.to_csv('{}.processed.results'.format(out_prefix),sep=' ',index=False)
	finemap_z = processed[['rsID','CHR', 'BP','EA', 'NEA','MAF', 'BETA', 'SE']]
	finemap_z.columns = ['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']
	finemap_z.to_csv('{}.processed.z'.format(out_prefix),sep=' ',index=False)
	processed[['rsID','Zscore',]].to_csv('{}.processed.caviarbf'.format(out_prefix),sep=' ',index=False,header=False)
	processed['coding'].to_csv('{}.processed.annotations'.format(out_prefix),sep=' ',index=False,header=True)
	np.savetxt('{}.processed.ld'.format(out_prefix), np.corrcoef(processed.iloc[:,18:-2].values), fmt='%1.4e')

def run_block(chr_id, start, stop):
	out_prefix = '{}/{}_{}_{}_{}.processed'.format(out_dir,prefix,chr_id,start,stop)
	prefix_name = '{}_{}_{}_{}.processed'.format(prefix,chr_id,start,stop)

    # PAINTOR
	call('echo "{}" > {}.input'.format(prefix_name,out_prefix),shell=True)
	call('{} -input {}.input -out {} -Zhead Zscore -LDname ld -enumerate {} -in {}'.format(paintor,out_prefix,out_dir,max_causal,out_dir),shell=True)

    # CAVIARBF
	n_variants = check_output('wc -l {}.caviarbf'.format(out_prefix),shell=True).decode().split()[0]
	call('{} -z {}.caviarbf -r {}.ld -t 0 -a 0.1281429 -n {} -c {} -o {}.caviarbf.out'.format(caviarbf,out_prefix,out_prefix,sample_size,max_causal,out_prefix),shell=True)
	call('{} -i {}.caviarbf.out -m {} -p 0 -o {}.prior0'.format(model_search,out_prefix,n_variants,out_prefix),shell=True)

    # FINEMAP
	call('echo "z;ld;snp;config;cred;log;n_samples\n{}.z;{}.ld;{}.snp;{}.config;{}.cred;{}.log;{}" > {}.master'.format(out_prefix,out_prefix,out_prefix,out_prefix,out_prefix,out_prefix,sample_size,out_prefix),shell=True)
	call('{} --sss --in-files {}.master --n-causal-snps {}'.format(finemap,out_prefix,max_causal),shell=True)

def merge_results(chr_id, start, stop):
	out_prefix = '{}/{}_{}_{}_{}.processed'.format(out_dir,prefix,chr_id,start,stop)
	paintor = pd.read_csv('{}.results'.format(out_prefix),sep=' ',)
	paintor.columns = ['CHR', 'BP', 'rsID', 'MAF', 'EA', 'NEA', 'BETA', 'SE', 'P', 'Zscore','index', 'PAINTOR']
	caviarbf = pd.read_csv('{}.prior0.marginal'.format(out_prefix),sep=' ',names=['No','caviarbf']).sort_values('No')
	finemap = pd.read_csv('{}.config'.format(out_prefix),sep=' ',usecols=['config','prob'],index_col='config')

	merge = paintor.copy()
	merge['CAVIARBF'] = merge['PAINTOR']
	merge['FINEMAP'] = merge['PAINTOR']
	merge['CAVIARBF'] = caviarbf['caviarbf'].values
	merge['FINEMAP'] = [finemap.loc[rs,'prob'] for rs in merge['rsID']]
	merge['LD'] = check_output('head -n {} {}.ld | tail -n 1'.format(merge['P'].idxmin()+1,out_prefix),shell=True).decode().split()
	merge.index = merge['index']
	# merge.to_csv('{}.causal'.format(out_prefix),sep=' ',index=False)

	block = pd.read_csv('{}/{}_{}_{}_{}.txt'.format(out_dir,prefix,chr_id,start,stop),sep='\t')
	block['PAINTOR'] = -1
	block['CAVIARBF'] = -1
	block['FINEMAP'] = -1
	block['LD'] = -1
	block.index = block['index']
	block.update(merge.drop_duplicates('index'))
	block = block.astype({'CHR':int,'BP':int})
	block.to_csv('{}/{}_{}_{}_{}.causal.txt'.format(out_dir,prefix,chr_id,start,stop),sep='\t',index=False,header=False)

# read blocks file
paintor = '/f/jianhua/software/PAINTOR_V3.0/PAINTOR'
caviarbf = '/f/jianhua/software/caviarbf/caviarbf'
model_search = '/f/jianhua/software/caviarbf/model_search'
finemap = '/f/jianhua/software/finemap_v1.3.1_x86_64/finemap_v1.3.1_x86_64'
blocks_file = '/f/jianhua/causal_db/ref/blocks.txt'
reference_dir = '/f/jianhua/causal_db/ref/ld/txt'
blocks = pd.read_csv(blocks_file, sep='\t')
threads = 10
max_causal = sys.argv[3]
sample_size = sys.argv[4]

# read input file
input_file = sys.argv[1]
pop = sys.argv[2]
pwd = os.getcwd()
prefix = input_file.split('/')[-1].split('.')[0]
out_dir = '{}/{}'.format(pwd,prefix)
if not os.path.exists(out_dir):
	os.mkdir(out_dir)
print('loading input file')
df = pd.read_csv(input_file, sep='\t')
df = df.dropna(subset=['CHR','BP','EA','NEA','BETA','SE','P','Zscore'])
if pd.isnull(df.iloc[0,3]):
	pass
else:
	df = df[(df['MAF']>0) & (df['MAF']<=0.5)]
df['index'] = df.index
df = df[df['P']>0]
significant_df = df[df['P']<=5e-8]

if len(significant_df) == 0:
	significant_df = df[df['P'].idxmin():df['P'].idxmin()+1]

# split files
print('split files')
significant_blocks = []
for chr_id in range(1,23):
    df_chr = df[df['CHR']==chr_id]
    significant_chr = significant_df[significant_df['CHR']==chr_id]
    block_chr = blocks[blocks['chr']==chr_id]
    for significant_bp in significant_chr['BP']:
        for i in block_chr.index:
            if block_chr.loc[i,'start']<=significant_bp<=block_chr.loc[i,'stop']:
                if i not in significant_blocks:
                    significant_blocks.append(i)
                    start,stop = block_chr.loc[i,'start'],block_chr.loc[i,'stop']
                    df_chr[(df_chr['BP']>=start) & (df_chr['BP']<=stop)].to_csv('{}/{}_{}_{}_{}.txt'.format(out_dir,prefix,chr_id,start,stop), sep='\t',index=False)
significant_blocks = blocks.loc[significant_blocks]

# prepare input files of finemapping tools
print('prepare input files of finemapping tools')
p = Pool(threads)
for i in significant_blocks.index:
	chr_id, start, stop = significant_blocks.loc[i].values
	p.apply_async(extract_overlap_and_cal_ld,(pop,chr_id, start, stop))
p.close()
p.join()

# run finemapping tools
print('run finemapping tools')
p = Pool(threads)
for i in significant_blocks.index:
	chr_id, start, stop = significant_blocks.loc[i].values
	p.apply_async(run_block,(chr_id, start, stop))
p.close()
p.join()

# gather results from three tools
print('gather results')
p = Pool(threads)
for i in significant_blocks.index:
	chr_id, start, stop = significant_blocks.loc[i].values
	p.apply_async(merge_results,(chr_id, start, stop))
p.close()
p.join()

# remove intermediate files
for file in os.listdir(out_dir):
	if file.endswith('.causal.txt'):
		pass
	else:
		os.remove(out_dir +'/'+ file)
significant_blocks.to_csv('{}/{}/{}_significant_blocks.txt'.format(pwd,prefix,prefix),sep='\t',index=False)
# call('rm -rf {}'.format(out_dir),shell=True)

