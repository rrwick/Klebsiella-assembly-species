# _Klebsiella_ assembly species

This tool is designed to assign species to a _Klebsiella_ assembly. However, unlike other tools, like [Kleborate](https://github.com/katholt/Kleborate/) which assign a single species, this tool 'paints' each region with a species, allowing it to detect and characterise cross-species hybrids.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6462461.svg)](https://doi.org/10.5281/zenodo.6462461)



## Requirements

[BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK279690/) commands must be available on the command line (specifically the commands `makeblastdb`, `blastn` and `tblastn`). BLAST+ can usually be easily installed using a package manager such as [Homebrew](http://brew.sh/) (on Mac) or [apt-get](https://help.ubuntu.com/community/AptGet/Howto) (on Ubuntu and related Linux distributions).

You'll need reference genome assemblies for each _Klebsiella_ species: _pneumoniae_, _quasipneumoniae_subsp_quasipneumoniae_, _quasipneumoniae_subsp_similipneumoniae_, _variicola_, _quasivariicola_, _oxytoca_, _michiganensis_, _grimontii_ and _aerogenes_. Put these in a directory and make sure that the species name (exactly as shown here) is in the path or filename. It's better to have a good variety of genomes from each species, as this will better capture the pan genome. However, more genomes will cause the tool to run slower, so you may need to strike a balance. Also, ensure that no hybrids are in your reference set, e.g. a genome which is a blend of _pneumoniae_ and _variicola_.



## Installation

This tool is a single Python 3 script with no third-party dependencies. It will run without any installation:
```
git clone https://github.com/rrwick/Klebsiella-assembly-species
Klebsiella-assembly-species/klebsiella_assembly_species.py --help
```

If you want, you can copy the it to somewhere in your PATH for easy usage:
```
cp Klebsiella-assembly-species/klebsiella_assembly_species.py ~/.local/bin
klebsiella_assembly_species.py --help
```



## Example commands

__Test one assembly__:<br>
`klebsiella_assembly_species.py --reference_dir references klebs_assembly.fasta`

__Test many assemblies__:<br>
`klebsiella_assembly_species.py --reference_dir references klebs_assemblies/*.fasta`

__Save blocks for plotting in R__:<br>
`klebsiella_assembly_species.py --reference_dir references --save_blocks block_dir klebs_assemblies/*.fasta`



## Full usage

```
usage: klebsiella_assembly_species.py [-h] [--hit_length HIT_LENGTH]
                                      [--hit_id HIT_ID] [--min_diff MIN_DIFF]
                                      [--min_block_size MIN_BLOCK_SIZE]
                                      [--reference_dir REFERENCE_DIR]
                                      [--save_blocks SAVE_BLOCKS]
                                      [--save_match_vals SAVE_MATCH_VALS]
                                      queries [queries ...]

Klebsiella species assembly checker

positional arguments:
  queries               FASTA file(s) of assembly to check

optional arguments:
  -h, --help            show this help message and exit
  --hit_length HIT_LENGTH
                        BLAST hits shorter than this will be ignored (default:
                        100)
  --hit_id HIT_ID       BLAST hits with lower identity than this will be
                        ignored (default: 98.0)
  --min_diff MIN_DIFF   Bases will not be assigned to a species if they have
                        less than this much difference between the best and
                        second-best species (default: 0.5)
  --min_block_size MIN_BLOCK_SIZE
                        Blocks smaller than this size will be filtered out
                        (default: 10000)
  --reference_dir REFERENCE_DIR
                        Where to find the reference genomes (default: same
                        directory as script)
  --save_blocks SAVE_BLOCKS
                        Saves blocks files for each query to this directory
                        (default: do not save blocks)
  --save_match_vals SAVE_MATCH_VALS
                        Saves blocks files for each query to this directory
                        (default: do not save match values)
```



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
