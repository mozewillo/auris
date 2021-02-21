import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import os.path


def check_alphabet(sequence):
    """check alphabet correctness of a single sequence(domain)"""
    sequence = set(sequence)
    alphabet = {'A', 'C', 'G', 'T', 'N', '-'}
    if len(sequence.difference(alphabet)) != 0:
        raise IOError("Given sequence contains unsupported characters")


def check_genomes(genomes, reference_id):
    """check alphabet correctness of genomes"""
    genomes = genomes[genomes['genome'].str.contains('^[^BDEFHIJKLMOPRSUW]+$')]
    try:
        ref_genome = genomes.loc[reference_id]['genome']
        ref_len = len(ref_genome)
        genomes = genomes[genomes['genome'].str.len() == ref_len]
        return genomes, ref_genome
    except KeyError:
        print("Given reference sequence contains unsupported characters")


def group_genomes(genomes):
    """group entries by genomes (same sequence)"""
    genomes = genomes.reset_index()
    genomes_num = len(genomes)
    genomes['count'] = genomes.groupby('genome')['id'].transform('count')
    aggregation_functions = {'id': 'min', 'genome': 'first', 'count': 'first'}
    genomes = genomes.groupby('genome', as_index=False).aggregate(aggregation_functions)
    return genomes, genomes_num


def plot_variants(gvariants, out_vis, threshold, all=False):
    """plot variants frequency in population"""
    gvariants = gvariants[gvariants['count'] > threshold]
    gvariants = gvariants.sort_values(by='frequency', ascending=False)
    # plot all the variants above the given threshold
    if all and int(np.round(len(gvariants) / 50)) > 1:
        num_of_plots = int(np.round(len(gvariants) / 50))
        fig, axs = plt.subplots(num_of_plots, figsize=(14, 5 * num_of_plots), sharey=True)
        for i in range(num_of_plots):
            if i != num_of_plots - 1:
                vars = gvariants[i * 50:(i + 1) * 50]
            else:
                vars = gvariants[i * 50:]
            vars.plot(kind="bar", x="symbol", y="frequency", ax=axs[i])
    # plot the first 100 of the most frequent variants
    elif len(gvariants) > 60:
        mid = int(len(gvariants) / 2)
        vars1 = gvariants[:mid]
        vars2 = gvariants[mid:100]
        fig, (ax1, ax2) = plt.subplots(2, figsize=(14, 12), sharey=True)
        ax1.set_ylabel('frequency')
        ax2.set_ylabel('frequency')
        ax1.set_title('Most common variants frequencies')
        vars1.plot(kind="bar", x="symbol", y="frequency", ax=ax1)
        # ax1.set_ylim([0,1])
        vars2.plot(kind="bar", x="symbol", y="frequency", ax=ax2)
        plt.show()
    # data fitted in one readable plot
    else:
        gvariants.plot(kind="bar", x="symbol", y="frequency", figsize=(14, 10), ylim=(0, 1),
                       title="Variants frequencies")
    plt.savefig(out_vis)


def variant_call(gens, ref_gen, gens_num, out_vis, threshold, all):
    """call the nucleotide that differs from the reference genome"""
    gvariants = {'ref': [], 'pos': [], 'alt': [], 'symbol': [], 'count': [], 'genome_id': []}
    for ind, row in gens.iterrows():
        i = 1
        g = row['genome']
        id = row['id']
        c = row['count']
        for n, m in zip(g, ref_gen):
            if n != m:
                gvariants['ref'].append(m)
                gvariants['alt'].append(n)
                gvariants['pos'].append(i)
                s = m + str(i) + n
                gvariants['symbol'].append(s)
                gvariants['count'].append(c)
                gvariants['genome_id'].append(id)
            i += 1
    gvariants = pd.DataFrame.from_dict(gvariants)
    aggregation_functions = {'symbol': 'first', 'ref': 'first', 'pos': 'first', 'alt': 'first', 'count': 'sum',
                             "genome_id": 'first'}
    gvariants = gvariants.groupby(gvariants['symbol']).aggregate(aggregation_functions)
    gvariants['frequency'] = gvariants['count'] / gens_num
    print(f"{len(gvariants)} unique variants found")
    plot_variants(gvariants, out_vis, threshold, all)
    return gvariants


def domain_search(genomes, domain):
    """search for given domain and return 1 if it is present """
    regex = ''
    for n in domain:
        regex = regex + '[' + n + 'N](\-)*'
    is_present = []
    for ind, row in genomes.iterrows():
        genome = row['genome']
        domain_found = int(bool(re.search(regex, genome)))
        is_present.append(domain_found)
    genomes['isDomainPresent'] = is_present
    return genomes


def out_file(filename, ext):
    """check file extension"""
    if not filename.endswith(ext):
        filename = filename + ext
    return filename


