import pandas as pd
from analyze.auris import *

parser = argparse.ArgumentParser(description='Specify the parameters to run analyses.')
parser.add_argument('input_csv', type=str, help="Path to input CSV")
parser.add_argument('out_vis', type=str, help="Path to output pdf")
parser.add_argument('out_csv', type=str, help="Path to output CSV")
parser.add_argument('searched_domain', type=str, help="String to be searched")
parser.add_argument('reference_id', type=str, help="ID of the reference genome")
parser.add_argument('-count_threshold', '--ct', type=int, required=False, default=1,
                    help="Threshold for plotting variations")
parser.add_argument('-a', '--a', type=bool, default=False, help="plot all the variants above the threshold ,"
                                                                "default = the first 100 of the most frequent variants")
args = parser.parse_args()


if __name__ == '__main__':
    if os.path.isfile(args.input_csv) and args.input_csv.endswith('.csv'):  # check user parameters and validate data
        gens = pd.read_csv(args.input_csv, index_col='id')
    else:
        raise IOError("Given input file does not exist or does not have csv extension")
    domain = args.searched_domain
    check_alphabet(domain)
    csv_out = out_file(args.out_csv, '.csv')
    out_vis = out_file(args.out_vis, '.pdf')
    reference_id = args.reference_id
    counts_threshold = args.ct
    plot_all = args.a
    gens, ref_genome = check_genomes(gens, reference_id)
    print(f"Loaded {len(gens)} genomes with the correct alphabet [A, C, G, T, N, -]")
    gens, genomes_num = group_genomes(gens)
    print(f"{len(gens)} unique genomes found")                               # call variants and search for the domain
    variants = variant_call(gens, ref_genome, genomes_num, out_vis, counts_threshold, plot_all)
    gens = domain_search(gens, domain)
    gens.to_csv(csv_out)