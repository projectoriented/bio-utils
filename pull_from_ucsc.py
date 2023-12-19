#!/usr/bin/env python3

"""
Usage: pull_from_ucsc.py --track 14455_p1 --genome --chrom
This script gets tracks from UCSC.
"""

import requests
import pandas as pd
import sys
import argparse

# Check Python version
if not sys.version_info >= (3, 6):
    sys.exit("Please have python>=3.6 exec in your system")


def get_parser():
    """Get parser"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--track', '-t', type=str, required=True,
                        help='Desired track name')
    parser.add_argument('--genome', '-g', required=False, help='Desired genome build. (default: hg38)',
                        default='hg38')
    parser.add_argument('--chrom', '-c', type=str, required=False,
                        help='Desired chromosome')
    parser.add_argument('--output', '-o', type=str, required=False, default=sys.stdout,
                        help='Output')
    return parser


def main():
    """Get bed"""
    parser = get_parser()
    args = parser.parse_args()

    target_params = {
        'genome': args.genome,
        'track': args.track,
    }

    if args.chrom:
        target_params['chrom'] = args.chrom

    get_df(target_params).to_csv(args.output, sep='\t', index=False)


def get_df(target_params: dict):
    breakdown_params = ';'.join(map('='.join, target_params.items()))
    api_url = f'https://api.genome.ucsc.edu/getData/track?{breakdown_params}'

    response = requests.get(api_url)
    json_res = response.json()

    track = target_params['track']

    is_list = isinstance(json_res[track], list)
    if is_list:
        df = pd.json_normalize(json_res, track)
    else:
        df = pd.DataFrame()
        for k in json_res[track].keys():
            df = pd.concat([df, pd.json_normalize(json_res, [track, k])])

    return df


if __name__ == "__main__":
    sys.exit(main())
