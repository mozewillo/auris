
##AURIS
Package for identifying common mutations and analyzing short genomes for the existence of specific
substrings.

#### INSTALLATION

**Download**  
Download package from github: https://github.com/mozewillo/auris

**Install**  
pip install path_to/auris

**Usage**  
python -m auris.analyze <input.csv> <plot_output.pdf> <output.csv> <domain_sequence>
<reference_genome_id> [-count_threshold CT] [-a A]

**Help**  
python -m auris.analyze -h

```
Specify the parameters to run analyses.

positional arguments:
  input_csv             Path to input CSV
  out_vis               Path to output pdf
  out_csv               Path to output CSV
  searched_domain       String to be searched
  reference_id          ID of the reference genome

optional arguments:
  -h, --help            show this help message and exit
  -count_threshold CT, --ct CT
                        Threshold for plotting variations
  -a A, --a A           plot all the variants above the threshold ,default =
                        the first 100 of the most frequent variants
```

