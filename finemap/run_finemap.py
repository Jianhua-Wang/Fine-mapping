from wjhfinemap import *
import pandas as pd
import sys

sumstat = pd.read_csv(sys.argv[1],sep='\t')
sumstat = make_uniqueID(sumstat)
loci = get_distance_loci(sumstat)
cojo_loci = get_cojo_loci(sumstat, loci, '/h/jianhua/REF/1000G/20130502/pop_plink/EUR.chr{chrom}', 10000)
cojo_loci.to_csv(sys.argv[4], sep='\t', index=False)
sample_size = int(sys.argv[3])

bfile_pattern = '/h/jianhua/REF/1000G/20130502/pop_plink/EUR.chr{chrom}'
cojo_loci["sample_size"] = sample_size
cojo_loci["bfile"] = cojo_loci["chr"].apply(lambda x: bfile_pattern.format(chrom=x))
cojo_loci["sub_df"] = cojo_loci.index.map(
    lambda x: sumstat[
        (sumstat["CHR"] == cojo_loci.loc[x, "chr"])
        & (sumstat["BP"] >= cojo_loci.loc[x, "start"])
        & (sumstat["BP"] <= cojo_loci.loc[x, "end"])
    ]
)

cojo_loci["cond_snps"] = cojo_loci.index.map(
    lambda x: cojo_loci[
        (cojo_loci["chr"] == cojo_loci.loc[x, "chr"])
        & (cojo_loci["lead_bp"] >= cojo_loci.loc[x, "start"])
        & (cojo_loci["lead_bp"] <= cojo_loci.loc[x, "end"])
        & (cojo_loci["lead_bp"] != cojo_loci.loc[x, "lead_bp"])
    ]["lead_snp"].values
)

p = Pool(20)
cojo_loci["cond_sumstat"] = p.starmap(cojo_cond, cojo_loci[['bfile','sub_df','cond_snps','sample_size']].values)
p.close()

p = Pool(20)
cojo_loci["prob"] = p.starmap(run_finemap, cojo_loci[['cond_sumstat','bfile','sample_size']].values)
p.close()

all_prob = pd.concat(cojo_loci["prob"].values, ignore_index=True)
all_prob.to_csv(sys.argv[2], sep='\t', index=False)