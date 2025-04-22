import logging
import math
import re

import numpy as np
import pandas as pd
from scipy.stats import chi2

from gsMap.config import FormatSumstatsConfig

VALID_SNPS = {"AC", "AG", "CA", "CT", "GA", "GT", "TC", "TG"}
logger = logging.getLogger(__name__)

default_cnames = {
    "chromosome": "CHR",
    # RS NUMBER
    "variant_id": 'SNP',
    "SNP": "SNP",
    "RS": "SNP",
    "RSID": "SNP",
    "rs_id": "SNP",
    "RS_NUMBER": "SNP",
    "RS_NUMBERS": "SNP",
    # P-VALUE
    "P": "P",
    "PVALUE": "P",
    "p_value": "P",
    "P_VALUE": "P",
    "Pvalue": "P",
    "PVAL": "P",
    "P_VAL": "P",
    "GC_PVALUE": "P",
    "p": "P",
    # EFFECT_ALLELE (A1)
    "A1": "A1",
    "ALLELE1": "A1",
    "ALLELE_1": "A1",
    "EFFECT_ALLELE": "A1",
    "effect_allele": "A1",
    "REFERENCE_ALLELE": "A1",
    "INC_ALLELE": "A1",
    "EA": "A1",
    # NON_EFFECT_ALLELE (A2)
    "A2": "A2",
    "ALLELE2": "A2",
    "ALLELE_2": "A2",
    "OTHER_ALLELE": "A2",
    "other_allele": "A2",
    "NON_EFFECT_ALLELE": "A2",
    "DEC_ALLELE": "A2",
    "NEA": "A2",
    # N
    "N": "N",
    "NCASE": "N_CAS",
    "CASES_N": "N_CAS",
    "N_CASE": "N_CAS",
    "N_CASES": "N_CAS",
    "n_cases": "N_CAS",
    "N_CONTROLS": "N_CON",
    "n_controls": "N_CON",
    "N_CAS": "N_CAS",
    "N_CON": "N_CON",
    "NCONTROL": "N_CON",
    "CONTROLS_N": "N_CON",
    "N_CONTROL": "N_CON",
    "WEIGHT": "N",
    # SIGNED STATISTICS
    "ZSCORE": "Z",
    "Z-SCORE": "Z",
    "GC_ZSCORE": "Z",
    "Z": "Z",
    "OR": "OR",
    "B": "BETA",
    "BETA": "BETA",
    "LOG_ODDS": "LOG_ODDS",
    "EFFECTS": "BETA",
    "EFFECT": "BETA",
    "b": "BETA",
    "beta": "BETA",
    
    # SE
    "se": "SE",
    "StdErr": "SE",
    "standard_error": "SE",
    # INFO
    "INFO": "INFO",
    "Info": "INFO",
    # MAF
    "EAF": "FRQ",
    "FRQ": "FRQ",
    "MAF": "FRQ",
    "FRQ_U": "FRQ",
    "F_U": "FRQ",
    "frq_A1": "FRQ",
    "frq": "FRQ",
    "freq": "FRQ",
}


def get_compression(fh):
    """
    Read filename suffixes and figure out whether it is gzipped,bzip2'ed or not compressed
    """
    if fh.endswith("gz"):
        compression = "gzip"
    elif fh.endswith("bz2"):
        compression = "bz2"
    else:
        compression = None

    return compression


def gwas_checkname(gwas, config):
    """
    Iterpret column names of gwas
    """
    old_name = gwas.columns
    mapped_cnames = {}
    for col in gwas.columns:
        mapped_cnames[col] = default_cnames.get(col, col)
    gwas.columns = list(mapped_cnames.values())

    # When column names are provided by users
    name_updates = {
        "SNP": config.snp,
        "A1": config.a1,
        "A2": config.a2,
        "INFO": config.info,
        "BETA": config.beta,
        "SE": config.se,
        "P": config.p,
        "FRQ": config.frq,
        "N": config.n,
        "Z": config.z,
        "Chr": config.chr,
        "Pos": config.pos,
        "OR": config.OR,
        "SE_OR": config.se_OR,
    }

    for key, value in name_updates.items():
        if value is not None and value in gwas.columns:
            gwas.rename(columns={value: key}, inplace=True)
    new_name = gwas.columns
    # check the name duplication
    for head in new_name:
        numc = list(new_name).count(head)
        if numc > 1:
            raise ValueError(
                f"Found {numc} different {head} columns, please check your {head} column."
            )

    name_dict = {new_name[i]: old_name[i] for i in range(len(new_name))}

    # When at OR scale
    if "OR" in new_name and "SE_OR" in new_name:
        gwas["BETA"] = gwas.OR.apply(lambda x: math.log(x) if x > 0 else None)
        gwas["SE"] = gwas.SE_OR.apply(lambda x: math.log(x) if x > 0 else None)

    interpreting = {
        "SNP": "Variant ID (e.g., rs number).",
        "A1": "Allele 1, interpreted as the effect allele for signed sumstat.",
        "A2": "Allele 2, interpreted as the non-effect allele for signed sumstat.",
        "BETA": "[linear/logistic] regression coefficient (0 → no effect; above 0 → A1 is trait/risk increasing).",
        "SE": "Standard error of the regression coefficient.",
        "OR": "Odds ratio, will be transferred to linear scale.",
        "SE_OR": "Standard error of the odds ratio, will be transferred to linear scale.",
        "P": "P-Value.",
        "Z": "Z-Value.",
        "N": "Sample size.",
        "INFO": "INFO score (imputation quality; higher → better imputation).",
        "FRQ": "Allele frequency of A1.",
        "Chr": "Chromsome.",
        "Pos": "SNP positions.",
    }

    logger.info("\nIterpreting column names as follows:")
    for key, _value in interpreting.items():
        if key in new_name:
            logger.info(f"{name_dict[key]}: {interpreting[key]}")

    return gwas


def gwas_checkformat(gwas, config):
    """
    Check column names required for different format
    """
    if config.format == "gsMap":
        condition1 = np.any(np.isin(["P", "Z"], gwas.columns))
        condition2 = np.all(np.isin(["BETA", "SE"], gwas.columns))
        if not (condition1 or condition2):
            raise ValueError(
                "To munge GWAS data into gsMap format, either P or Z values, or both BETA and SE values, are required."
            )
        else:
            if "Z" in gwas.columns:
                pass
            elif "P" in gwas.columns:
                gwas["Z"] = np.sqrt(chi2.isf(gwas.P, 1)) * np.where(gwas["BETA"] < 0, -1, 1)
            else:
                gwas["Z"] = gwas.BETA / gwas.SE

    elif config.format == "COJO":
        condition = np.all(np.isin(["A1", "A2", "FRQ", "BETA", "SE", "P", "N"], gwas.columns))
        if not condition:
            raise ValueError(
                "To munge GWAS data into COJO format, either A1|A2|FRQ|BETA|SE|P|N, are required."
            )
        else:
            gwas["Z"] = np.sqrt(chi2.isf(gwas.P, 1)) * np.where(gwas["BETA"] < 0, -1, 1)

    return gwas


def filter_info(info, config):
    """Remove INFO < args.info_min (default 0.9) and complain about out-of-bounds INFO."""
    if type(info) is pd.Series:  # one INFO column
        jj = ((info > 2.0) | (info < 0)) & info.notnull()
        ii = info >= config.info_min
    elif type(info) is pd.DataFrame:  # several INFO columns
        jj = ((info > 2.0) & info.notnull()).any(axis=1) | ((info < 0) & info.notnull()).any(
            axis=1
        )
        ii = info.sum(axis=1) >= config.info_min * (len(info.columns))
    else:
        raise ValueError("Expected pd.DataFrame or pd.Series.")

    bad_info = jj.sum()
    if bad_info > 0:
        msg = "WARNING: {N} SNPs had INFO outside of [0,1.5]. The INFO column may be mislabeled."
        logger.warning(msg.format(N=bad_info))

    return ii


def filter_frq(frq, config):
    """
    Filter on MAF. Remove MAF < args.maf_min and out-of-bounds MAF.
    """
    jj = (frq < 0) | (frq > 1)
    bad_frq = jj.sum()
    if bad_frq > 0:
        msg = "WARNING: {N} SNPs had FRQ outside of [0,1]. The FRQ column may be mislabeled."
        logger.warning(msg.format(N=bad_frq))

    frq = np.minimum(frq, 1 - frq)
    ii = frq > config.maf_min
    return ii & ~jj


def filter_pvals(P, config):
    """Remove out-of-bounds P-values"""
    ii = (P > 0) & (P <= 1)
    bad_p = (~ii).sum()
    if bad_p > 0:
        msg = "WARNING: {N} SNPs had P outside of (0,1]. The P column may be mislabeled."
        logger.warning(msg.format(N=bad_p))

    return ii


def filter_alleles(a):
    """Remove alleles that do not describe strand-unambiguous SNPs"""
    return a.isin(VALID_SNPS)


def gwas_qc(gwas, config):
    """
    Filter out SNPs based on INFO, FRQ, MAF, N, and Genotypes.
    """
    old = len(gwas)
    logger.info("\nFiltering SNPs as follows:")
    # filter: SNPs with missing values
    drops = {"NA": 0, "P": 0, "INFO": 0, "FRQ": 0, "A": 0, "SNP": 0, "Dup": 0, "N": 0}

    gwas = gwas.dropna(
        axis=0, how="any", subset=filter(lambda x: x != "INFO", gwas.columns)
    ).reset_index(drop=True)

    drops["NA"] = old - len(gwas)
    logger.info(f"Removed {drops['NA']} SNPs with missing values.")

    # filter: SNPs with Info < 0.9
    if "INFO" in gwas.columns:
        old = len(gwas)
        gwas = gwas.loc[filter_info(gwas["INFO"], config)]
        drops["INFO"] = old - len(gwas)
        logger.info(f"Removed {drops['INFO']} SNPs with INFO <= 0.9.")

    # filter: SNPs with MAF <= 0.01
    if "FRQ" in gwas.columns:
        old = len(gwas)
        gwas = gwas.loc[filter_frq(gwas["FRQ"], config)]
        drops["FRQ"] += old - len(gwas)
        logger.info(f"Removed {drops['FRQ']} SNPs with MAF <= 0.01.")

    # filter: P-value that out-of-bounds [0,1]
    if "P" in gwas.columns:
        old = len(gwas)
        gwas = gwas.loc[filter_pvals(gwas["P"], config)]
        drops["P"] += old - len(gwas)
        logger.info(f"Removed {drops['P']} SNPs with out-of-bounds p-values.")

    # filter: Variants that are strand-ambiguous
    if "A1" in gwas.columns and "A2" in gwas.columns:
        gwas.A1 = gwas.A1.str.upper()
        gwas.A2 = gwas.A2.str.upper()
        gwas = gwas.loc[filter_alleles(gwas.A1 + gwas.A2)]
        drops["A"] += old - len(gwas)
        logger.info(f"Removed {drops['A']} variants that were not SNPs or were strand-ambiguous.")

    # filter: Duplicated rs numbers
    if "SNP" in gwas.columns:
        old = len(gwas)
        gwas = gwas.drop_duplicates(subset="SNP").reset_index(drop=True)
        drops["Dup"] += old - len(gwas)
        logger.info(f"Removed {drops['Dup']} SNPs with duplicated rs numbers.")

    # filter:Sample size
    n_min = gwas.N.quantile(0.9) / 1.5
    old = len(gwas)
    gwas = gwas[gwas.N >= n_min].reset_index(drop=True)
    drops["N"] += old - len(gwas)
    logger.info(f"Removed {drops['N']} SNPs with N < {n_min}.")

    return gwas


def variant_to_rsid(gwas, config):
    """
    Convert variant id (Chr, Pos) to rsid
    """
    logger.info("\nConverting the SNP position to rsid. This process may take some time.")
    unique_ids = set(gwas["id"])
    chr_format = gwas["Chr"].unique().astype(str)
    chr_format = [re.sub(r"\d+", "", value) for value in chr_format][1]

    dtype = {"chr": str, "pos": str, "ref": str, "alt": str, "dbsnp": str}
    chunk_iter = pd.read_csv(
        config.dbsnp,
        chunksize=config.chunksize,
        sep="\t",
        skiprows=1,
        dtype=dtype,
        names=["chr", "pos", "ref", "alt", "dbsnp"],
    )

    # Iterate over chunks
    matching_id = pd.DataFrame()
    for chunk in chunk_iter:
        chunk["id"] = chr_format + chunk["chr"] + "_" + chunk["pos"]
        matching_id = pd.concat(
            [matching_id, chunk[chunk["id"].isin(unique_ids)][["dbsnp", "id"]]]
        )

    matching_id = matching_id.drop_duplicates(subset="dbsnp").reset_index(drop=True)
    matching_id = matching_id.drop_duplicates(subset="id").reset_index(drop=True)
    matching_id.index = matching_id.id
    return matching_id


def clean_SNP_id(gwas, config):
    """
    Clean SNP id
    """
    old = len(gwas)
    condition1 = "SNP" in gwas.columns
    condition2 = np.all(np.isin(["Chr", "Pos"], gwas.columns))

    if not (condition1 or condition2):
        raise ValueError("Either SNP rsid, or both SNP chromosome and position, are required.")
    elif condition1:
        pass
    elif condition2:
        if config.dbsnp is None:
            raise ValueError("To Convert SNP positions to rsid, dbsnp reference is required.")
        else:
            gwas["id"] = gwas["Chr"].astype(str) + "_" + gwas["Pos"].astype(str)
            gwas = gwas.drop_duplicates(subset="id").reset_index(drop=True)
            gwas.index = gwas.id

            matching_id = variant_to_rsid(gwas, config)
            gwas = gwas.loc[matching_id.id]
            gwas["SNP"] = matching_id.dbsnp
            num_fail = old - len(gwas)
            logger.info(f"Removed {num_fail} SNPs that did not convert to rsid.")

    return gwas


def gwas_metadata(gwas, config):
    """
    Report key features of GWAS data
    """
    logger.info("\nSummary of GWAS data:")
    CHISQ = gwas.Z**2
    mean_chisq = CHISQ.mean()
    logger.info("Mean chi^2 = " + str(round(mean_chisq, 3)))
    if mean_chisq < 1.02:
        logger.warning("Mean chi^2 may be too small.")

    logger.info("Lambda GC = " + str(round(CHISQ.median() / 0.4549, 3)))
    logger.info("Max chi^2 = " + str(round(CHISQ.max(), 3)))
    logger.info(
        f"{(CHISQ > 29).sum()} Genome-wide significant SNPs (some may have been removed by filtering)."
    )


def gwas_format(config: FormatSumstatsConfig):
    """
    Format GWAS data
    """
    logger.info(f"------Formating gwas data for {config.sumstats}...")
    compression_type = get_compression(config.sumstats)
    gwas = pd.read_csv(
        config.sumstats,
        delim_whitespace=True,
        header=0,
        compression=compression_type,
        na_values=[".", "NA"],
    )

    if isinstance(config.n, int | float):
        logger.info(f"Set the sample size of gwas data as {config.n}.")
        gwas["N"] = config.n
        config.n = "N"

    logger.info(f"Read {len(gwas)} SNPs from {config.sumstats}.")

    # Check name and format
    gwas = gwas_checkname(gwas, config)
    gwas = gwas_checkformat(gwas, config)
    # Clean the snp id
    gwas = clean_SNP_id(gwas, config)
    # QC
    gwas = gwas_qc(gwas, config)
    # Meta
    gwas_metadata(gwas, config)

    # Saving the data
    if config.format == "COJO":
        keep = ["SNP", "A1", "A2", "FRQ", "BETA", "SE", "P", "N"]
        appendix = ".cojo"
    elif config.format == "gsMap":
        keep = ["SNP", "A1", "A2", "Z", "N"]
        appendix = ".sumstats"

    if "Chr" in gwas.columns and "Pos" in gwas.columns and config.keep_chr_pos is True:
        keep = keep + ["Chr", "Pos"]

    gwas = gwas[keep]
    out_name = config.out + appendix + ".gz"

    logger.info(f"\nWriting summary statistics for {len(gwas)} SNPs to {out_name}.")
    gwas.to_csv(out_name, sep="\t", index=False, float_format="%.3f", compression="gzip")
