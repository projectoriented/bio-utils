#!/usr/bin/env python3
"""
Usage: ./sex-estimator.py your.bam --threads 8
Check the sex of sample.
Author: Mei Wu, https://github.com/projectoriented
"""
import os
import sys
import logging
import argparse

import pandas as pd
import pysam

from datetime import datetime

from io import StringIO
import multiprocessing as mp

LOG = logging.getLogger()
logging.basicConfig(stream=sys.stderr, level="INFO", format='%(asctime)s - %(levelname)s - %(message)s')

# These are gathered from ASD cohort
TRUE_MALEY_RATIO = 0.009847771
TRUE_FEMALEY_RATIO = 0.002462753

datef = datetime.today().strftime('%Y%m%d')


def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    parser.add_argument("filepath", nargs=1, type=str)
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument("--maleY", type=int, default=TRUE_MALEY_RATIO,
                        help="The Y ratio estimated from a truely female individual")
    parser.add_argument("--femaleY", type=int, default=TRUE_FEMALEY_RATIO,
                        help="The Y ratio estimated from a truely male individual")
    parser.add_argument("--rl", type=int, default=10000,
                        help="The required minimum read length. This is technology-specific.")

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()

    filepath = args.filepath[0]
    threads = args.threads
    maleY = args.maleY
    femaleY = args.femaleY
    min_rl = args.rl

    autosomes_set = ['chr{}'.format(x) for x in list(range(1, 23))]
    genome_dict = {entry: "autosome" for entry in autosomes_set}
    genome_dict.update({"chrY": "chrY"})

    df = gather_dataframes(filepath=filepath, num_processes=threads, target_chroms=["chrY"] + autosomes_set, rl=min_rl)
    df["sample"] = filepath

    df["genome_type"] = df.apply(lambda row: genome_dict.get(row["#rname"], ""), axis=1)

    stats = df.groupby(["genome_type"], group_keys=False)["numreads"].sum().reset_index().set_index("genome_type").T
    stats["ratioY"] = stats["chrY"] / stats["autosome"]
    stats["Yestimate"] = (stats["ratioY"] - femaleY) / (maleY - femaleY)
    stats["assigned_sex"] = stats["Yestimate"].apply(lambda x: guesstimate_sex(x=x, male_ratio=maleY))
    stats["sample"] = filepath

    stats.reset_index(drop=True, inplace=True)
    df.columns.name = ""

    file_basename = os.path.basename(filepath)

    os.makedirs(datef, exist_ok=True)
    prefix = os.path.join(datef, file_basename)
    stats.to_csv(f"{prefix}-Ystats.tsv", header=True, index=False, sep="\t")
    df.to_csv(f"{prefix}-stats.tsv", header=True, index=False, sep="\t")


def guesstimate_sex(x, male_ratio):
    target = (male_ratio / 1.3)
    if x >= target:
        return "male"
    elif target > x > (male_ratio / 1.5):
        return "maybe-male"
    else:
        return "female"


def gather_dataframes(filepath, num_processes, target_chroms, rl):
    arglist = []
    for entry in target_chroms:
        arglist.append(make_adl_args(region=entry, min_read_len=rl))

    filepath_list = [filepath] * len(arglist)

    with mp.Pool(processes=num_processes) as pool:
        results = pool.starmap(get_coverage, zip(filepath_list, arglist))

    results = [result for result in results]
    final_result = pd.concat(results, ignore_index=True)

    return final_result


def make_adl_args(**kwargs) -> list:
    adl_args = []
    for k, v in kwargs.items():
        if "_" in k:
            k = k.replace("_", "-")
        adl_args.append(f"--{k}")
        adl_args.append(str(v))
    return adl_args


def get_coverage(filepath, adl_args: list) -> pd.DataFrame:
    try:
        out = pysam.coverage(filepath, *adl_args)
        df = pd.read_table(StringIO(out), header=0)
    except pysam.utils.SamtoolsError:
        LOG.warning(f"{adl_args} is invalid.")
        df = pd.DataFrame()
    return df


if __name__ == "__main__":
    sys.exit(main())
