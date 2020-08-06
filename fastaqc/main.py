#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import logging
from pbr.version import VersionInfo
from Bio import SeqIO
import itertools
from prettytable import PrettyTable
import tabulate
from xopen import xopen
from fastaqc.alphabet import Alphabet
import pprint

__version__ = VersionInfo('fastaqc').semantic_version().release_string()

def main():
    parser = argparse.ArgumentParser(description='Version ' + __version__ + '\nCheck fasta file', formatter_class=RawTextHelpFormatter)
    parser.set_defaults(func=help)
    parser.add_argument('--verbose', '-v', action="store_true")

    subparsers = parser.add_subparsers()
    info_parser = subparsers.add_parser('info')
    info_parser.add_argument('fasta', nargs=argparse.REMAINDER)
    info_parser.set_defaults(func=info)

    args = parser.parse_args()
    config = {} # implement when needed
    if ('verbose' in vars(args) and args.verbose):
        logging.basicConfig(level=logging.INFO)
    args.parser = parser
    args.func(args, config)

def help(args, config):
    args.parser.print_help()

def info(args, cfg):
  if not args.fasta:
    help(args, cfg)
    return
  # Concept
  #
  # The base method iterates through all sequences in the fasta file and applies multiple 
  # actions on it. Each action is modelled as a function with two parameters, the sequence
  # record and the stats-dictionary. The stats-dictionary holds the results of each action.
  # Subsequent actions can access the results from previous actions.
  #
  # Coonvetions for stats dictionary
  #  * Intermediary results for only one sequence (those that will not be presented to the user)
  #    start with an underscore.
  checks = [
    count, 
    collect_lengths,
    compute_character_distribution,
    compute_character_positions,
    detect_sequence_type,
    detect_ambiguous_and_special_characters,
    set_sequence_category_name,
    count_sequence_types,
    count_sequences_with_special_characters,
    count_sequences_with_ambiguous_characters,
    count_sequences_with_unknown_characters,
    clear_temporary_fields
  ]

  for filename in args.fasta:
    with xopen(filename) as fh:
      stats = {
        'filename': filename
      }
      for record in SeqIO.parse(fh, "fasta"):
        for c in checks:
          c(record, stats)
      print_stats(stats)

def clear_temporary_fields(record, stats):
  for_removal = []
  for k in stats.keys():
    if k.startswith('_'):
      for_removal.append(k)
  for k in for_removal:
    del stats[k]

def print_stats(stats):
  header = [stats['filename'], 'count']
  table = []
  table.append(['sequences', stats['sequences']])
  logging.info(pprint.pformat(stats))
  for k,v in sorted(stats['type_counts'].items()):
    table.append(['   ' + k, v])
    if 'special_char_count' in stats:
      if k in stats['special_char_count']:
        for k2,v2 in sorted(stats['special_char_count'][k].items()):
          table.append(['      ' + k2,v2])
    if 'ambiguous_char_count' in stats:
      if k in stats['ambiguous_char_count']:
        for k2,v2 in sorted(stats['ambiguous_char_count'][k].items()):
          table.append(['      ' + k2,v2])
    if 'unknown_char_count' in stats:
      if k in stats['unknown_char_count']:
        for k2,v2 in sorted(stats['unknown_char_count'][k].items()):
          table.append(['      ' + k2,v2])

  tabulate.PRESERVE_WHITESPACE = True
  print(tabulate.tabulate(table, headers=header, tablefmt='pretty', colalign=('left', 'right')))
  if merge(stats['unknown_char_count']):
    print("WARNING: The file contains unknown characters for DNA, RNA and AA sequences. ")
    print("         It will probably fail in applications with strict alphabet checking.")

  if 'seq_lenghts' in stats:
    import plotille
    print('')
    print('Sequence length distribution')
    print(plotille.histogram(stats['seq_lenghts'], height=25, x_min=0))
  print('')

def merge(dict_of_dicts):
  merged = {}
  for name, subdict in dict_of_dicts.items():
    for k,v in subdict.items():
      if k not in merged:
        merged[k] = v
      else:
        merged[k] = merged[k] + v
  return merged

def count(record, stats):
  '''counts the number of processed sequences in the field "sequences".'''
  if 'sequences' not in stats:
    stats['sequences'] = 0
  stats['sequences'] = stats['sequences'] + 1

def compute_character_distribution(record, stats):
  distribution = {}
  for c in record.seq.upper():
    if c not in distribution:
      distribution[c] = 0
    distribution[c] = distribution[c] + 1

  stats['_character_distribution'] = distribution

def compute_character_positions(record, stats):
  '''computes a dictionary with all positions per character of the current 
  sequence and stores it in "_character_positions"'''
  positions = {}
  for i, c in enumerate(record.seq.upper()):
    if c not in positions:
      positions[c] = []
    positions[c].append(i)
  stats['_character_positions'] = positions

def count_sequences_with_special_characters(record, stats):
  '''counts the sequences with special characters (depends on the alphabet)
  and stores them in "special_char_count.<sequence_category>.<character>"'''
  alphabet = assert_sequence_type_available(stats)
  _count(stats, 'special_char_count', alphabet.special_chars)

def count_sequences_with_ambiguous_characters(record, stats):
  '''counts the sequences with ambiguous characters (depends on the alphabet)
  and stores them in "ambiguous_char_count.<sequence_category>.<character>"'''
  alphabet = assert_sequence_type_available(stats)
  _count(stats, 'ambiguous_char_count', alphabet.ambiguous_chars)

def _count(stats, fieldname, chars):
  category_name = assert_sequence_category_name_available(stats)
  c_dist = assert_character_distribution_available(stats)
  if fieldname not in stats:
    stats[fieldname] = {}
  if category_name not in stats[fieldname]:
    stats[fieldname][category_name] = {}
  counts = stats[fieldname][category_name]
  for c in chars:
    if c in c_dist:
      if c not in counts:
        counts[c] = 0
      counts[c] = counts[c] + 1

def count_sequences_with_unknown_characters(record, stats):
  '''counts the sequences with unknown characters (depends on the alphabet)
  and stores them in "ambiguous_char_count.<sequence_category>.<character>"'''
  category_name = assert_sequence_category_name_available(stats)
  alphabet = assert_sequence_type_available(stats)
  c_dist = assert_character_distribution_available(stats)
  if 'unknown_char_count' not in stats:
    stats['unknown_char_count'] = {}
  if category_name not in stats['unknown_char_count']:
    stats['unknown_char_count'][category_name] = {}
  counts = stats['unknown_char_count'][category_name]
  chars = set(alphabet.all_chars)
  for c in c_dist.keys():
    if c not in chars:
      if c not in counts:
        counts[c] = 0
      counts[c] = counts[c] + 1

def detect_sequence_type(record, stats):
  c_dist = assert_character_distribution_available(stats)
  if _contains_only(c_dist, Alphabet.DNA.all_chars):
    stats['_type'] = Alphabet.DNA
  elif _contains_only(c_dist, Alphabet.RNA.all_chars):
    stats['_type'] = Alphabet.RNA
  elif _contains_only(c_dist, Alphabet.AA.all_chars):
    stats['_type'] = Alphabet.AA
  else:
    stats['_type'] = Alphabet.OTHER

def detect_ambiguous_and_special_characters(record, stats):
  c_dist = assert_character_distribution_available(stats)
  alphabet = assert_sequence_type_available(stats)
  type_flags = set()
  if _contains_only(c_dist, alphabet.unambiguous_chars):
    type_flags.add('unambiguous')
  else:
    if _contains(c_dist, alphabet.ambiguous_chars):
      type_flags.add('ambiguous')
    if _contains(c_dist, alphabet.special_chars):
      type_flags.add('special')
  stats['_type_flags'] = type_flags

def set_sequence_category_name(record, stats):
  alphabet = assert_sequence_type_available(stats)
  flags = assert_sequence_type_flags_available(stats)
  name = "{} ({})".format(alphabet.name, ",".join(sorted(flags)))
  stats['_category_name'] = name

def count_sequence_types(record, stats):
  alphabet = assert_sequence_type_available(stats)
  flags = assert_sequence_type_flags_available(stats)
  type = assert_sequence_category_name_available(stats)
  if 'type_counts' not in stats:
    stats['type_counts'] = {}
  if type not in stats['type_counts']:
    stats['type_counts'][type] = 0
  stats['type_counts'][type] = stats['type_counts'][type]  + 1

def collect_lengths(record, stats):
  if 'seq_lenghts' not in stats:
    stats['seq_lenghts'] = []
  stats['seq_lenghts'].append(len(record.seq))

def _contains(dist, characters):
  for c in characters:
    if c in dist:
      return True
  return False

def _contains_only(dist, characters):
  remaining = len(dist)
  for c in characters:
    if c in dist:
      remaining = remaining - 1
  return remaining == 0

def assert_character_distribution_available(stats):
  assert '_character_distribution' in stats, 'Sequence character distribution not availabe. It must be computed before this check.'
  return stats['_character_distribution']

def assert_sequence_type_available(stats):
  assert '_type' in stats, 'Sequence type information not availabe. It must be computed before this check.'
  return stats['_type']

def assert_sequence_type_flags_available(stats):
  assert '_type_flags' in stats, 'Sequence type flag information not availabe. It must be computed before this check.'
  return stats['_type_flags']

def assert_sequence_category_name_available(stats):
  assert '_category_name' in stats, 'Sequence category_name information not availabe. It must be computed before this check.'
  return stats['_category_name']

if __name__ == "__main__":
    main()
