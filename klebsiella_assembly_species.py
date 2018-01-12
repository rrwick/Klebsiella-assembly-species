#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/DASCRUBBER_wrapper

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not,
see <http://www.gnu.org/licenses/>.
"""


import os
import subprocess
import argparse
import collections
import statistics


def get_arguments():
    parser = argparse.ArgumentParser(description='Klebsiella species assembly checker',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--hit_length', type=int, default=100,
                        help='BLAST hits shorter than this will be ignored')
    parser.add_argument('--hit_id', type=float, default=98.0,
                        help='BLAST hits with lower identity than this will be ignored')
    parser.add_argument('--min_diff', type=float, default=0.5,
                        help='Bases will not be assigned to a species if they have less than this '
                             'much difference between the best and second-best species')
    parser.add_argument('--min_block_size', type=float, default=10000,
                        help='Blocks smaller than this size will be filtered out')
    parser.add_argument('--reference_dir', type=str, default='same directory as script',
                        help='Where to find the reference genomes')
    parser.add_argument('--save_blocks', type=str, default='do not save blocks',
                        help='Saves blocks files for each query to this directory')
    parser.add_argument('--save_match_vals', type=str, default='do not save match values',
                        help='Saves blocks files for each query to this directory')
    parser.add_argument('queries', type=str, nargs='+',
                        help='FASTA file(s) of assembly to check')

    args = parser.parse_args()
    if args.reference_dir == 'same directory as script':
        args.reference_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'references')

    if args.save_blocks == 'do not save blocks':
        args.save_blocks = None
    elif not os.path.exists(args.save_blocks):
        os.makedirs(args.save_blocks)

    if args.save_match_vals == 'do not save match values':
        args.save_match_vals = None
    elif not os.path.exists(args.save_match_vals):
        os.makedirs(args.save_match_vals)

    return args


def main():
    args = get_arguments()
    ref_fasta_files = get_reference_fasta_files(args.reference_dir)
    for ref_fasta in ref_fasta_files:
        build_blast_db(ref_fasta)
    ref_species = get_ref_species(ref_fasta_files)
    print_table_header(ref_species)

    # For memory efficiency, we'll store a contig's per-base consensus using a number index.
    index_1 = {0: 'no call',
               1: 'too close',
               2: 'pneumoniae',
               3: 'quasipneumoniae_subsp_quasipneumoniae',
               4: 'quasipneumoniae_subsp_similipneumoniae',
               5: 'variicola',
               6: 'quasivariicola',
               7: 'oxytoca',
               8: 'michiganensis',
               9: 'grimontii',
               10: 'aerogenes'}
    index_2 = {call: num for num, call in index_1.items()}

    for query in args.queries:
        names_lengths = get_contig_names_and_lengths(query)
        results = initialise_results(names_lengths, ref_species)

        for ref_fasta in ref_fasta_files:
            blast_hits = [h for h in blast_query_against_reference(ref_fasta, query)
                          if h.alignment_length >= args.hit_length
                          and h.percent_identity >= args.hit_id]
            add_results(results, blast_hits, ref_fasta)

        query_totals = {i: 0 for i in index_1.keys()}
        query_blocks = []

        contig_names_size_order = [x[0] for x in
                                   sorted(list(names_lengths.items()),
                                          key=lambda x: x[1], reverse=True)]

        # Loop through contigs in order of largest to smallest.
        for contig_name in contig_names_size_order:
            contig_length = names_lengths[contig_name]
            contig_results = results[contig_name]
            consensus = get_per_base_results(contig_results, contig_length, args, ref_species,
                                             index_2)
            blocks = group_results_into_blocks(consensus, contig_length)
            blocks = filter_blocks(blocks, args.min_block_size)
            for block in blocks:
                query_totals[block[2]] += block[1] - block[0]
                query_blocks.append([contig_name, block[0], block[1], index_1[block[2]]])

        if args.save_match_vals:
            save_match_vals_files(query, results, args.save_match_vals, contig_names_size_order,
                                  ref_species)

        if args.save_blocks is not None:
            save_blocks_file(query, query_blocks, args.save_blocks)

        # Get the unknown fraction ('no call' and 'too close').
        unknown_base_count = sum(query_totals[i] for i in [0, 1])
        total_base_count = sum(query_totals.values())
        unknown_percentage = 100.0 * unknown_base_count / total_base_count

        # Get percentages for the query, not counting 'no call' and 'too close' bases.
        ref_indices = [index_2[s] for s in ref_species]
        base_count = sum(query_totals[i] for i in ref_indices)
        percentages = {index_1[i]: 100.0 * query_totals[i] / base_count
                       for i in ref_indices}
        print_table_row(query, ref_species, percentages, unknown_percentage)


def initialise_results(names_lengths, ref_species):
    """
    The results dictionary will hold the species identity per base per contig.
    Key = contig name, species name
    Value = list of identities
    """
    results = {}
    for contig_name, contig_length in names_lengths.items():
        results[contig_name] = {}
        for species in ref_species:
            results[contig_name][species] = [0.0] * contig_length
    return results


def add_results(results, blast_hits, ref_fasta):
    species = get_ref_species_from_filename(ref_fasta)
    for hit in blast_hits:
        contig_species_results = results[hit.contig_name][species]
        for i in range(hit.contig_start, hit.contig_end):
            if hit.percent_identity > contig_species_results[i]:
                contig_species_results[i] = hit.percent_identity


def get_per_base_results(contig_results, contig_length, args, ref_species, index_2):
    consensus = [0] * contig_length
    for i in range(contig_length):
        options = sorted([(contig_results[species][i], species) for species in ref_species],
                         reverse=True)

        # Don't assign a species if the hit is too weak.
        if options[0][0] < args.hit_id:
            continue

        # Don't assign a species if the first and second-best hits are too close.
        elif options[0][0] - options[1][0] < args.min_diff:
            consensus[i] = 1
            continue

        # Assign this base the species of the best match.
        consensus[i] = index_2[options[0][1]]

    return consensus


def group_results_into_blocks(consensus, contig_length):
    blocks = []
    block_start = 0
    for i in range(1, contig_length):
        if consensus[i] != consensus[i - 1]:
            blocks.append((block_start, i, consensus[i - 1]))
            block_start = i
    blocks.append((block_start, contig_length, consensus[-1]))
    return blocks


def filter_blocks(blocks, min_block_size):
    block_count = len(blocks)
    filtered_blocks = []
    for i, block in enumerate(blocks):
        block_size = block[1] - block[0]
        block_group = block[2]

        # 'no call' and 'too close' blocks always pass the filter.
        if block_group < 2:
            filtered_blocks.append(block)

        # Big blocks always pass the filter.
        elif block_size >= min_block_size:
            filtered_blocks.append(block)

        # Small blocks for a species will pass the filter if their neighbouring blocks are for the
        # same species. Specifically, we add the size of previous/following blocks to this block's
        # size if the species matches, then assess the adjusted block size.
        else:
            # Check previous blocks.
            for j in range(i-1, -1, -1):
                previous_block = blocks[j]
                previous_block_size = previous_block[1] - previous_block[0]
                previous_block_group = previous_block[2]
                if previous_block_group < 2:
                    continue
                elif previous_block_group != block_group:
                    break
                else:  # previous_block_group == block_group:
                    block_size += previous_block_size

            # Check following blocks.
            for j in range(i+1, block_count):
                following_block = blocks[j]
                following_block_size = following_block[1] - following_block[0]
                following_block_group = following_block[2]
                if following_block_group < 2:
                    continue
                elif following_block_group != block_group:
                    break
                else:  # following_block_group == block_group:
                    block_size += following_block_size

            if block_size >= min_block_size:
                filtered_blocks.append(block)
            else:  # Still too small
                filtered_blocks.append((block[0], block[1], 0))
    return filtered_blocks


def save_match_vals_files(query, results, directory, contig_names_size_order, ref_species):
    for contig_name in contig_names_size_order:
        contig_results = results[contig_name]
        for species in ref_species:
            species_results = contig_results[species]
            with open(get_values_filename(query, contig_name, species, directory), 'wt') as values_file:
                values_file.write('\t'.join(['start', 'end', 'match']))
                values_file.write('\n')
                range_start, range_end = 0, 100
                while range_end <= len(species_results):
                    val = statistics.mean(species_results[range_start:range_end])
                    values_file.write('\t'.join([str(range_start), str(range_end), '%.3f' % val]))
                    values_file.write('\n')
                    range_start += 100
                    range_end += 100


def save_blocks_file(query, query_blocks, directory):
    with open(get_block_filename(query, directory), 'wt') as blocks_file:
        blocks_file.write('\t'.join(['contig', 'start', 'end', 'classification']))
        blocks_file.write('\n')
        for block in query_blocks:
            blocks_file.write('\t'.join([str(x) for x in block]))
            blocks_file.write('\n')


def get_reference_fasta_files(reference_dir):
    return sorted([os.path.join(dir_path, f)
                   for dir_path, _, filenames in os.walk(reference_dir)
                   for f in filenames if f.endswith('.fasta') and not f.startswith('._')])


def get_ref_species(ref_fasta_files):
    ref_species = set(get_ref_species_from_filename(f) for f in ref_fasta_files)
    return [s for s in ['pneumoniae', 'quasipneumoniae_subsp_quasipneumoniae',
                        'quasipneumoniae_subsp_similipneumoniae', 'variicola', 'quasivariicola',
                        'oxytoca', 'michiganensis', 'grimontii', 'aerogenes']
            if s in ref_species]


def build_blast_db(ref_fasta):
    if not (os.path.isfile(ref_fasta + '.nhr') and
            os.path.isfile(ref_fasta + '.nin') and
            os.path.isfile(ref_fasta + '.nsq')):
        makeblastdb_cmd = 'makeblastdb -dbtype nucl -in ' + ref_fasta
        subprocess.check_output(makeblastdb_cmd, shell=True).decode()


def blast_query_against_reference(ref_fasta, query):
    blast_cmd = 'blastn -db ' + ref_fasta + ' -query ' + query + ' -outfmt 6'
    blast_output = subprocess.check_output(blast_cmd, shell=True).decode()
    blast_hits = [BlastHit(line) for line in blast_output.split('\n') if line]
    return sorted(blast_hits, key=lambda x: x.alignment_length, reverse=True)


def get_ref_species_from_filename(ref_fasta):
    if 'quasipneumoniae_subsp_quasipneumoniae' in ref_fasta.lower():
        return 'quasipneumoniae_subsp_quasipneumoniae'
    if 'quasipneumoniae_subsp_similipneumoniae' in ref_fasta.lower():
        return 'quasipneumoniae_subsp_similipneumoniae'
    elif 'pneumoniae' in ref_fasta.lower():
        return 'pneumoniae'
    elif 'quasivariicola' in ref_fasta.lower():
        return 'quasivariicola'
    elif 'variicola' in ref_fasta.lower():
        return 'variicola'
    elif 'oxytoca' in ref_fasta.lower():
        return 'oxytoca'
    elif 'michiganensis' in ref_fasta.lower():
        return 'michiganensis'
    elif 'grimontii' in ref_fasta.lower():
        return 'grimontii'
    elif 'aerogenes' in ref_fasta.lower():
        return 'aerogenes'
    else:
        assert False


def print_table_header(ref_species):
    print('\t'.join(['query', 'unknown'] + ref_species))


def print_table_row(query, ref_species, percentages, unknown_percentage):
    row = [os.path.basename(query), '%.2f' % unknown_percentage] + ['%.2f' % percentages[s] for s in ref_species]
    print('\t'.join(row))


def get_values_filename(query, contig, species, directory):
    query_name = os.path.basename(query)
    if query_name.endswith('.fasta'):
        query_name = query_name[:-6]
    if query_name.endswith('.fa'):
        query_name = query_name[:-3]
    if not os.path.exists(os.path.join(directory, query_name)):
        os.makedirs(os.path.join(directory, query_name))
    return os.path.join(directory, query_name, contig + '_' + species + '_values')


def get_block_filename(query, directory):
    block_filename = os.path.basename(query)
    if block_filename.endswith('.fasta'):
        block_filename = block_filename[:-6]
    if block_filename.endswith('.fa'):
        block_filename = block_filename[:-3]
    return os.path.join(directory, block_filename + '_blocks')


class BlastHit(object):
    def __init__(self, blast_line):
        line_parts = blast_line.split('\t')
        self.contig_name = line_parts[0]
        self.percent_identity = float(line_parts[2])
        self.alignment_length = int(line_parts[3])
        self.error_count = round(self.alignment_length * (100.0 - self.percent_identity) / 100.0)

        # Hit coordinates are stored in a Python-esque manner (0 based, exclusive end).
        self.contig_start = int(line_parts[6]) - 1
        self.contig_end = int(line_parts[7])

    def __repr__(self):
        return self.contig_name + ': ' + str(self.contig_start) + '-' + str(self.contig_end) + \
               ', ID: ' + str(self.percent_identity) + '%'


def get_contig_names_and_lengths(filename):
    """
    Returns a dictionary of contig name -> contig length.
    """
    names_lengths = collections.OrderedDict()
    with open(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    names_lengths[name] = len(sequence)
                    sequence = ''
                name = line[1:].split()[0]
            else:
                sequence += line
        if name:
            names_lengths[name] = len(sequence)
    return names_lengths


if __name__ == '__main__':
    main()
