import pandas as pd
import numpy as np
import uuid
from subprocess import call
import os
import tempfile
from multiprocessing import Pool
import time


def make_uniqueID(
    sumstat,
    cols={
        "chr_col": "CHR",
        "bp_col": "BP",
        "rsid_col": "rsID",
        "a1_col": "EA",
        "a2_col": "NEA",
    },
    inplace=False,
):
    df = sumstat.copy()
    allele_df = df[[cols["a1_col"], cols["a2_col"]]].copy()
    b = allele_df.values
    b.sort(axis=1)
    allele_df[[cols["a1_col"], cols["a2_col"]]] = b
    allele_df["SNP_ID"] = (
        df[cols["chr_col"]].astype(str)
        + "-"
        + df[cols["bp_col"]].astype(str)
        + "-"
        + allele_df[cols["a1_col"]]
        + "-"
        + allele_df[cols["a2_col"]]
    )
    if inplace:
        df[cols["rsid_col"]] = allele_df["SNP_ID"]
    else:
        df.insert(loc=0, column="SNP_ID", value=allele_df["SNP_ID"].values)
    return df


def merge_loci(sig_blocks):
    merge_blocks = sig_blocks.copy()
    while True:
        for i in merge_blocks.index:
            chrom, start, end = merge_blocks.loc[i, ["chr", "start", "end"]]
            sub_df = merge_blocks[
                (merge_blocks["chr"] == chrom)
                & (merge_blocks["end"] >= start)
                & (merge_blocks["start"] <= end)
            ]
            lead_idx = sub_df["lead_p"].idxmin()
            merge_blocks.loc[i] = (
                sub_df.loc[lead_idx, "lead_snp"],
                sub_df.loc[lead_idx, "lead_p"],
                chrom,
                sub_df["start"].min(),
                sub_df["end"].max(),
            )
        merge_blocks = merge_blocks.drop_duplicates()
        overlap = 0
        for i in merge_blocks.index:
            chrom, start, end = merge_blocks.loc[i, ["chr", "start", "end"]]
            sub_df = merge_blocks[
                (merge_blocks["chr"] == chrom)
                & (merge_blocks["end"] >= start)
                & (merge_blocks["start"] <= end)
            ]
            if len(sub_df) > 1:
                overlap = 1
                break
        if not overlap:
            break

    merge_blocks = merge_blocks.drop_duplicates()
    merge_blocks.reset_index(drop=True, inplace=True)
    return merge_blocks


def get_distance_loci(
    sumstat, sig_cutoff=5e-8, window_size=500000, merge_overlap=True
):
    sig_df = sumstat[sumstat["P"] <= sig_cutoff].copy().sort_values("P")
    # sig_df["rsID"] = sig_df["CHR"].astype(str) + "-" + sig_df["BP"].astype(str)
    sig_blocks = pd.DataFrame(
        columns=["lead_snp", "lead_p", "chr", "start", "end"]
    )
    ith = 0
    while len(sig_df) > 0:
        chrom, lead_bp, lead_snp, lead_p = sig_df.iloc[0][
            ["CHR", "BP", "SNP_ID", "P"]
        ].values
        sig_blocks.loc[ith] = (
            lead_snp,
            lead_p,
            chrom,
            lead_bp - window_size,
            lead_bp + window_size,
        )
        sig_df = sig_df[
            ~(
                (sig_df["CHR"] == chrom)
                & (sig_df["BP"] >= lead_bp - window_size)
                & (sig_df["BP"] <= lead_bp + window_size)
            )
        ]
        ith += 1
    sig_blocks = sig_blocks.sort_values(["chr", "start", "end"])
    sig_blocks.reset_index(drop=True, inplace=True)
    if merge_overlap:
        sig_blocks = merge_loci(sig_blocks)
    return sig_blocks


def ref_intersec(bim_prefix, sumstat, out_prefix, mac=10, annotate_eaf=True):
    chrom, start, end = (
        sumstat["CHR"].values[0],
        sumstat["BP"].min(),
        sumstat["BP"].max(),
    )
    with open(f"{out_prefix}.region", "w") as f:
        f.write(f"{chrom}\t{start}\t{end}\tregion")

    call(
        f"""plink --bfile {bim_prefix} \
             --extract range {out_prefix}.region \
             --keep-allele-order --mac {mac} \
             --make-bed --out {out_prefix}.orig \
             > {out_prefix}.orig.log""",
        shell=True,
    )
    bim = pd.read_csv(
        f"{out_prefix}.orig.bim", delim_whitespace=True, header=None
    )
    bim[1] = bim.index
    bim.to_csv(f"{out_prefix}.orig.bim", sep="\t", index=False, header=False)
    call(
        f"""plink --bfile {out_prefix}.orig \
             --keep-allele-order \
             --list-duplicate-vars ids-only suppress-first \
             --out {out_prefix}.orig \
             > {out_prefix}.orig.rmdup.log""",
        shell=True,
    )
    call(
        f"""plink --bfile {out_prefix}.orig \
             --keep-allele-order \
             --exclude {out_prefix}.orig.dupvar \
             --make-bed --out {out_prefix}.orig.rmdup \
             > {out_prefix}.orig.rmdup.log""",
        shell=True,
    )

    bim = pd.read_csv(
        f"{out_prefix}.orig.rmdup.bim", delim_whitespace=True, header=None
    )
    bim = make_uniqueID(
        bim,
        cols={
            "chr_col": 0,
            "bp_col": 3,
            "rsid_col": 1,
            "a1_col": 4,
            "a2_col": 5,
        },
        inplace=True,
    )
    bim.to_csv(
        f"{out_prefix}.orig.rmdup.bim", sep="\t", index=False, header=False
    )
    refalt_snps = bim[3].astype(str) + "-" + bim[5] + "-" + bim[4]
    refalt_snps = refalt_snps.values
    sub_df = sumstat.copy()
    sub_df["refalt"] = (
        sub_df["BP"].astype(str) + "-" + sub_df["NEA"] + "-" + sub_df["EA"]
    )
    sub_df["altref"] = (
        sub_df["BP"].astype(str) + "-" + sub_df["EA"] + "-" + sub_df["NEA"]
    )

    sub_df = sub_df[
        (sub_df["refalt"].isin(refalt_snps))
        | (sub_df["altref"].isin(refalt_snps))
    ]
    sub_df.reset_index(inplace=True, drop=True)
    b = sub_df[["EA", "NEA"]].copy()
    sub_df["EA"] = b["NEA"].where(sub_df["altref"].isin(refalt_snps), b["EA"])
    sub_df["NEA"] = b["EA"].where(sub_df["altref"].isin(refalt_snps), b["NEA"])
    sub_df["BETA"] = -sub_df["BETA"].where(
        sub_df["altref"].isin(refalt_snps), sub_df["BETA"]
    )
    sub_df["Zscore"] = -sub_df["Zscore"].where(
        sub_df["altref"].isin(refalt_snps), sub_df["Zscore"]
    )

    # sub_df["snp"] = sub_df["BP"].astype(str) + "-" + sub_df["NEA"] + "-" + sub_df["EA"]
    sub_df = sub_df.sort_values(["CHR", "BP"])
    sub_df.index = sub_df["SNP_ID"].values
    sub_df["SNP_ID"].to_csv(
        f"{out_prefix}.snplist", sep="\t", header=False, index=False
    )

    call(
        f"""plink --bfile {out_prefix}.orig.rmdup \
            --extract {out_prefix}.snplist \
            --make-bed \
            --keep-allele-order \
            --out {out_prefix} > {out_prefix}.log""",
        shell=True,
    )
    if annotate_eaf:
        call(
            f"plink --bfile {out_prefix} --freq --out {out_prefix} > {out_prefix}.frq.log",
            shell=True,
        )
        freq = pd.read_csv(
            f"{out_prefix}.frq", delim_whitespace=True, index_col="SNP"
        )
        freq["EA"] = sub_df["EA"]
        freq["EAF"] = freq["MAF"].where(
            freq["EA"] == freq["A1"], 1 - freq["MAF"]
        )
        sub_df["EAF"] = freq["EAF"]
        sub_df["MAF"] = freq["MAF"]
    sub_df = sub_df.drop(["refalt", "altref"], axis=1)
    sub_df.reset_index(inplace=True, drop=True)
    sub_df.to_csv(f"{out_prefix}.sumstat", sep="\t", index=False)
    # call(f"rm {out_prefix}.orig.*", shell=True)
    return sub_df


def cojo_slct(bfile, sumstat, sample_size, sig_cutoff=5e-8, maf=0.01):
    with tempfile.TemporaryDirectory(dir="./") as tmpdirname:
        out_prefix = f"{tmpdirname}/{str(uuid.uuid1())}"
        df = ref_intersec(bfile, sumstat, out_prefix, annotate_eaf=True)
        df["N"] = sample_size
        df = df[["SNP_ID", "EA", "NEA", "EAF", "BETA", "SE", "P", "N"]]
        df.columns = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
        df.to_csv(f"{out_prefix}.ma", sep=" ", index=False)
        cmd = f"""gcta64 --bfile {out_prefix} \
                --cojo-p {sig_cutoff} \
                --cojo-collinear 0.9 \
                --maf {maf} \
                --cojo-file {out_prefix}.ma \
                --cojo-slct \
                --out {out_prefix}.slct > {out_prefix}.slct.log"""
        call(cmd, shell=True)
        if os.path.exists(f"{out_prefix}.slct.jma.cojo"):
            slct_res = pd.read_csv(f"{out_prefix}.slct.jma.cojo", sep="\t")
            slct_res = slct_res[
                (slct_res["p"] <= sig_cutoff) & (slct_res["pJ"] <= sig_cutoff)
            ]
        else:
            slct_res = None
        # call(f'rm {out_prefix}*', shell=True)
    return slct_res


def get_cojo_loci(
    allsumstat,
    loci,
    bfile_pattern,
    sample_size,
    window_size=1000000,
    threads=10,
):
    run_loci = loci.copy()
    run_loci["sample_size"] = sample_size
    run_loci["bfile"] = run_loci["chr"].apply(
        lambda x: bfile_pattern.format(chrom=x)
    )
    run_loci["sub_df"] = run_loci.index.map(
        lambda x: allsumstat[
            (allsumstat["CHR"] == run_loci.loc[x, "chr"])
            & (allsumstat["BP"] >= run_loci.loc[x, "start"])
            & (allsumstat["BP"] <= run_loci.loc[x, "end"])
        ]
    )
    p = Pool(threads)
    slct_res = p.starmap(
        cojo_slct, run_loci[["bfile", "sub_df", "sample_size"]].values
    )
    p.close()
    slct_res = pd.concat(slct_res, ignore_index=True)
    sig_blocks = pd.DataFrame(
        columns=["lead_snp", "lead_p", "lead_bp", "chr", "start", "end"]
    )
    for i in slct_res.index:
        sig_blocks.loc[i] = (
            slct_res.loc[i, "SNP"],
            slct_res.loc[i, "p"],
            slct_res.loc[i, "bp"],
            slct_res.loc[i, "Chr"],
            slct_res.loc[i, "bp"] - window_size,
            slct_res.loc[i, "bp"] + window_size,
        )
    return sig_blocks


def cojo_cond(bfile, sumstat, cond_snps, sample_size):
    if len(cond_snps) == 0:
        return sumstat.copy()
    with tempfile.TemporaryDirectory(dir="./") as tmpdirname:
        out_prefix = f"{tmpdirname}/{str(uuid.uuid1())}"
        overlap_sumstat = ref_intersec(
            bfile, sumstat, out_prefix, annotate_eaf=True
        )
        cojo_input = overlap_sumstat.copy()
        cojo_input["N"] = sample_size
        cojo_input = cojo_input[
            ["SNP_ID", "EA", "NEA", "EAF", "BETA", "SE", "P", "N"]
        ]
        cojo_input.columns = ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]
        cojo_input.to_csv(f"{out_prefix}.ma", sep=" ", index=False)
        with open(f"{out_prefix}.snps", "w") as f:
            f.write("\n".join(cond_snps))
        cmd = f"""gcta64 --bfile {out_prefix} \
                --cojo-file {out_prefix}.ma \
                --cojo-cond {out_prefix}.snps \
                --out {out_prefix}.cond > {out_prefix}.cond.log"""
        call(cmd, shell=True)
        if os.path.exists(f"{out_prefix}.cond.cma.cojo"):
            cond_res = pd.read_csv(f"{out_prefix}.cond.cma.cojo", sep="\t")
            overlap_sumstat.index = overlap_sumstat["SNP_ID"].values
            cond_res.index = cond_res["SNP"].values
            overlap_sumstat["bC"] = cond_res["bC"]
            overlap_sumstat["bC_se"] = cond_res["bC_se"]
            overlap_sumstat["pC"] = cond_res["pC"]
            overlap_sumstat["Zscore"] = (
                overlap_sumstat["BETA"] / overlap_sumstat["SE"]
            )
            overlap_sumstat.reset_index(drop=True, inplace=True)
        else:
            print(out_prefix, "cojo_cond")
            time.sleep(1800)
            raise RuntimeError("conditional analysis failed")
    return overlap_sumstat


def run_finemap(
    sumstat,
    bfile,
    sample_size,
    max_causal=1,
    credible_set_cutoff=0.95,
    cojo=True,
):
    with tempfile.TemporaryDirectory(dir="./") as tmpdirname:
        out_prefix = f"{tmpdirname}/{str(uuid.uuid1())}"
        overlap_sumstat = ref_intersec(
            bfile, sumstat, out_prefix, annotate_eaf=True
        )
        finemap_z = overlap_sumstat.copy()
        # if 'bC' not in finemap_z.columns:
        #     print(out_prefix, 'bC not in ')
        #     time.sleep(1800)
        if cojo and "bC" in finemap_z.columns:
            finemap_z = finemap_z[
                ["SNP_ID", "CHR", "BP", "EA", "NEA", "MAF", "bC", "bC_se"]
            ].copy()
        else:
            finemap_z = finemap_z[
                ["SNP_ID", "CHR", "BP", "EA", "NEA", "MAF", "BETA", "SE"]
            ].copy()
        finemap_z["MAF"] = finemap_z["MAF"].replace(0, 0.5)
        finemap_z.columns = [
            "rsid",
            "chromosome",
            "position",
            "allele1",
            "allele2",
            "maf",
            "beta",
            "se",
        ]
        finemap_z = finemap_z.dropna()
        finemap_z.to_csv(f"{out_prefix}.finemap.z", sep=" ", index=False)
        finemap_z["rsid"].to_csv(
            f"{out_prefix}.ld.snps", index=False, header=False
        )

        call(
            f"plink --bfile {out_prefix} --extract {out_prefix}.ld.snps --r2 square spaces --out {out_prefix} > {out_prefix}.ld.log",
            shell=True,
        )

        with open(f"{out_prefix}.master", "w") as f:
            f.write("z;ld;snp;config;cred;log;n_samples\n")
            f.write(
                f"{out_prefix}.finemap.z;{out_prefix}.ld;{out_prefix}.finemap.snp;{out_prefix}.finemap.config;{out_prefix}.finemap.cred;{out_prefix}.finemap.log;{sample_size}"
            )
        call(
            f"finemap --sss --in-files {out_prefix}.master --n-causal-snps {max_causal} > {out_prefix}.finemap.log",
            shell=True,
        )
        pp_df = sumstat.copy()
        if os.path.exists(f"{out_prefix}.finemap.snp"):
            finemap_res = pd.read_csv(
                f"{out_prefix}.finemap.snp", index_col="rsid", sep=" "
            )
            pp_df["prob"] = pp_df["SNP_ID"].map(finemap_res["prob"])
            pp_df = pp_df.sort_values("prob", ascending=False)

            incld_idx = []
            for idx, cumsum in pp_df["prob"].cumsum().items():
                incld_idx.append(idx)
                if cumsum >= credible_set_cutoff:
                    break
            pp_df = pp_df[pp_df.index.isin(incld_idx)]
        else:
            raise RuntimeError("Run Finemap failed!")
    return pp_df


def get_clump_loci():
    pass


if __name__ == "__main__":
    sample_size = 10000
    in_file = "./normalized/GD09516.txt.gz"
    out_loci = "./GD09516_cojoloci.txt"
    out_prob = "./GD09516_credsets.txt"
    sumstat = pd.read_csv(in_file, sep="\t")
    sumstat = make_uniqueID(sumstat)
    loci = get_distance_loci(sumstat)
    cojo_loci = get_cojo_loci(
        sumstat,
        loci,
        "/h/jianhua/REF/1000G/20130502/pop_plink/EUR.chr{chrom}",
        10000,
    )
    cojo_loci.to_csv(out_loci, sep="\t", index=False)
    # sample_size = int(sys.argv[3])

    bfile_pattern = "/h/jianhua/REF/1000G/20130502/pop_plink/EUR.chr{chrom}"
    cojo_loci["sample_size"] = sample_size
    cojo_loci["bfile"] = cojo_loci["chr"].apply(
        lambda x: bfile_pattern.format(chrom=x)
    )
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
    cojo_loci["cond_sumstat"] = p.starmap(
        cojo_cond,
        cojo_loci[["bfile", "sub_df", "cond_snps", "sample_size"]].values,
    )
    p.close()

    p = Pool(20)
    cojo_loci["prob"] = p.starmap(
        run_finemap, cojo_loci[["cond_sumstat", "bfile", "sample_size"]].values
    )
    p.close()

    all_prob = pd.concat(cojo_loci["prob"].values, ignore_index=True)
    all_prob.to_csv(out_prob, sep="\t", index=False)

