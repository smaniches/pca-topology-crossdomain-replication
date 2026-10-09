#!/usr/bin/env python3
"""Render the archived CPTAC covariance-null distribution as an auditable SVG.

Uses only Python's standard library. This graph is a histogram of the real
499 null draw values; it is not a schematic of the original point clouds.

Usage from the repository root:
  python audit/cptac_covariance_null/plot_null_distribution.py
  python audit/cptac_covariance_null/plot_null_distribution.py --check
"""
import argparse
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / 'results' / 'cptac_tumor_covariance_null_20261009'
SOURCE = RESULTS / '499_draws.json'
TARGET = RESULTS / 'null_distribution.svg'


def render(payload):
    draws = payload['null_max_h1']
    obs = float(payload['observed_max_h1'])
    n = len(draws)
    if n != 499 or payload['draws'] != n or not math.isfinite(obs):
        raise ValueError('Expected the verified 499-draw CPTAC experiment')
    if not all(math.isfinite(float(x)) for x in draws):
        raise ValueError('Null distribution contains nonfinite values')
    exceed = sum(x >= obs for x in draws)
    p = (exceed + 1) / (n + 1)
    mean = statistics.mean(draws)
    if (exceed != payload['null_exceedances'] or
        abs(p - payload['p_one_sided_plus_one']) > 1e-12 or
        abs(mean - payload['null_mean']) > 1e-10):
        raise ValueError('Raw null draws disagree with archived summary')

    # Fixed-width bins have an exact boundary at the observed statistic.
    # This ensures amber bars represent exactly draws >= observed.
    bin_width = 0.25
    left = obs - 12 * bin_width
    total_bins = 20
    counts = [0] * total_bins
    for value in draws:
        j = math.floor((value - left) / bin_width)
        if j < 0 or j >= total_bins:
            raise ValueError('A draw falls outside the documented axis bounds')
        counts[j] += 1
    if sum(counts[12:]) != exceed or sum(counts) != n:
        raise ValueError('Histogram tail does not match the empirical p-value')

    width, height = 1000, 515
    x0, x1 = 82, 954
    y0, y1 = 132, 394
    ymax = max(100, 25 * math.ceil(max(counts) / 25))
    def xp(value):
        return x0 + (value - left) / (total_bins * bin_width) * (x1 - x0)
    def yp(value):
        return y1 - (value / ymax) * (y1 - y0)
    observed_x = xp(obs)
    elements = []
    add = elements.append
    add('<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="515" viewBox="0 0 1000 515" role="img" aria-labelledby="title description">')
    add('  <title id="title">CPTAC tumor-only covariance-null distribution</title>')
    add(f'  <desc id="description">Histogram of {n} real covariance-preserving null draws. Observed maximum H1 persistence is {obs:.3f}; {exceed} null draws meet or exceed it. The one-sided plus-one Monte Carlo p-value is {p:.3f}. This does not meet a threshold of 0.05.</desc>')
    add('  <rect width="1000" height="515" fill="#ffffff"/>')
    add('  <text x="82" y="44" font-family="Arial,Helvetica,sans-serif" font-size="24" font-weight="700" fill="#172b45">Does the observed loop exceed a covariance-preserving null?</text>')
    add('  <text x="82" y="77" font-family="Arial,Helvetica,sans-serif" font-size="17" fill="#4d6076">CPTAC-CCRCC · 110 tumor samples · PCA(50) · 499 null draws</text>')
    add(f'  <rect x="{observed_x:.2f}" y="{y0}" width="{x1-observed_x:.2f}" height="{y1-y0}" fill="#fff6e9"/>')
    for tick in range(0, ymax+1, 25):
        y = yp(tick)
        add(f'  <line x1="{x0}" y1="{y:.2f}" x2="{x1}" y2="{y:.2f}" stroke="#dbe4ec" stroke-width="1"/>')
        add(f'  <text x="68" y="{y+5:.2f}" text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="14" fill="#4d6076">{tick}</text>')
    for j, count in enumerate(counts):
        px = xp(left + j * bin_width) + 1.7
        w = xp(left + (j+1)*bin_width) - xp(left+j*bin_width) - 3.4
        y = yp(count)
        color = '#d58a32' if j >= 12 else '#4b81b5'
        add(f'  <rect x="{px:.2f}" y="{y:.2f}" width="{w:.2f}" height="{y1-y:.2f}" fill="{color}"/>')
    add(f'  <line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#738399" stroke-width="1.4"/>')
    for tick in (3, 4, 5, 6, 7):
        x = xp(tick)
        add(f'  <line x1="{x:.2f}" y1="{y1}" x2="{x:.2f}" y2="{y1+5}" stroke="#738399"/>')
        add(f'  <text x="{x:.2f}" y="{y1+24}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="15" fill="#34485d">{tick}</text>')
    mean_x = xp(mean)
    add(f'  <line x1="{mean_x:.2f}" y1="{y0}" x2="{mean_x:.2f}" y2="{y1}" stroke="#53657a" stroke-width="2.4" stroke-dasharray="7 6"/>')
    add(f'  <line x1="{observed_x:.2f}" y1="{y0}" x2="{observed_x:.2f}" y2="{y1}" stroke="#a65e18" stroke-width="3.2"/>')
    add(f'  <text x="{mean_x-9:.2f}" y="{y0+16}" text-anchor="end" font-family="Arial,Helvetica,sans-serif" font-size="15" fill="#34485d">Null mean {mean:.3f}</text>')
    add(f'  <text x="{observed_x+10:.2f}" y="{y0+16}" font-family="Arial,Helvetica,sans-serif" font-size="15" font-weight="700" fill="#925012">Observed {obs:.3f}</text>')
    add(f'  <text x="{observed_x+10:.2f}" y="{y0+39}" font-family="Arial,Helvetica,sans-serif" font-size="14" fill="#925012">{exceed}/{n} nulls as large or larger</text>')
    add(f'  <text x="{observed_x+10:.2f}" y="{y0+59}" font-family="Arial,Helvetica,sans-serif" font-size="14" fill="#925012">Monte Carlo p = {p:.3f}</text>')
    add(f'  <text x="{(x0+x1)/2:.0f}" y="452" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="17" fill="#34485d">Maximum finite H1 persistence (Euclidean PCA score space)</text>')
    add('  <text x="82" y="491" font-family="Arial,Helvetica,sans-serif" font-size="14" fill="#4d6076">Post-preregistration sensitivity. Did not meet the 0.05 threshold; biology remains unresolved.</text>')
    add('</svg>')
    return '\n'.join(elements) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='Verify the checked-in SVG is byte-identical to the raw-data rendering')
    args = parser.parse_args()
    svg = render(json.loads(SOURCE.read_text(encoding='utf-8')))
    if args.check:
        if not TARGET.exists() or TARGET.read_text(encoding='utf-8') != svg:
            parser.error(f'{TARGET} does not match the archived null draws')
        print('Null distribution SVG matches all 499 archived draws')
    else:
        TARGET.write_text(svg, encoding='utf-8')
        print(f'Wrote {TARGET}')


if __name__ == '__main__':
    main()
