#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import logging
from pbr.version import VersionInfo
import itertools
from prettytable import PrettyTable
import tabulate
from xopen import xopen

__version__ = VersionInfo('fastaqc').semantic_version().release_string()

def main():
    parser = argparse.ArgumentParser(description='Version ' + __version__ + '\nValidate fasta file', formatter_class=RawTextHelpFormatter)
    parser.set_defaults(func=help)

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
  from Bio import SeqIO
  checks = [
    count, 
    compute_character_distribution,
    check_sequence_type,
    check_has_star,
    check_has_dot,
    check_has_X,
    check_has_N,
    check_has_dash,
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

def print_stats(stats):
  header = [stats['filename'], 'count']
  table = []
  table.append(['sequences', stats['sequences']])
  _add_row_if_non_null2(table, stats['sequence_types'], 'aminoacid_unambiguous', '   AA')
  _add_row_if_non_null2(table, stats['sequence_types'], 'aminoacid_ambiguous', '   AA (ambiguous)')
  _add_row_if_non_null2(table, stats['sequence_types'], 'dna_unambiguous', '   DNA')
  _add_row_if_non_null2(table, stats['sequence_types'], 'dna_ambiguous', '    DNA (ambiguous)')
  _add_row_if_non_null2(table, stats['sequence_types'], 'rna_ambiguous', '   RNA')
  _add_row_if_non_null2(table, stats['sequence_types'], 'rna_unambiguous', '   RNA (ambiguous)')
  _add_row_if_non_null2(table, stats['sequence_types'], 'other', '   other')
  table.append(['',''])
  table.append(['sequences with special characters',''])
  _add_row_if_non_null2(table, stats, 'Xs', '   X')
  _add_row_if_non_null2(table, stats, 'Ns', '   N')
  _add_row_if_non_null2(table, stats, 'stars', '   *')
  _add_row_if_non_null2(table, stats, 'dots', '   .')
  _add_row_if_non_null2(table, stats, 'dash', '   -')
  _add_row_if_non_null2(table, stats, 'unhandled_characters', 'unknown characters')
  tabulate.PRESERVE_WHITESPACE = True
  print(tabulate.tabulate(table, headers=header, tablefmt='fancy_grid', colalign=('left', 'right')))
  print('')

def _add_row_if_non_null2(table, stats, field, fieldname=None):
  if field in stats:
    if stats[field]:
      if fieldname:
        table.append([fieldname, stats[field]])
      else:
        table.append([field, stats[field]])

def print_stats_with_PTable(stats):
  table = PrettyTable(['field', 'count'])
  table.align['field'] = 'l'
  table.align['count'] = 'r'
  table.add_row(['total sequences', stats['sequences']])
  _add_row_if_non_null(table, stats['sequence_types'], 'aminoacid_unambiguous', 'AA')
  _add_row_if_non_null(table, stats['sequence_types'], 'aminoacid_ambiguous', 'AA with X')
  _add_row_if_non_null(table, stats['sequence_types'], 'dna_unambiguous', 'DNA')
  _add_row_if_non_null(table, stats['sequence_types'], 'dna_ambiguous', 'DNA with N, ...')
  _add_row_if_non_null(table, stats['sequence_types'], 'rna_ambiguous', 'RNA')
  _add_row_if_non_null(table, stats['sequence_types'], 'rna_unambiguous', 'RNA with N, ...')
  _add_row_if_non_null(table, stats['sequence_types'], 'other')
  _add_row_if_non_null(table, stats, 'Xs', '# sequences containing X')
  _add_row_if_non_null(table, stats, 'Ns', '# sequences containing N')
  _add_row_if_non_null(table, stats, 'stars', '# sequences containing *')
  _add_row_if_non_null(table, stats, 'dots', '# sequences containing .')
  _add_row_if_non_null(table, stats, 'dash', '# sequences containing -')
  _add_row_if_non_null(table, stats, 'unhandled_characters', 'unknown characters')

  print(table.get_string(title=stats['filename']))
  print('')


def _add_row_if_non_null(table, stats, field, fieldname=None):
  if field in stats:
    if stats[field]:
      if fieldname:
        table.add_row([fieldname, stats[field]])
      else:
        table.add_row([field, stats[field]])


def _print_if_non_null(template, stats, field):
  if field in stats:
    if stats[field] :
      print(template.format(stats[field]))

def count(record, stats):
  if 'sequences' not in stats:
    stats['sequences'] = 0
  stats['sequences'] = stats['sequences'] + 1

def compute_character_distribution(record, stats):
  distribution = {}
  for c in record.seq.upper():
    if c not in distribution:
      distribution[c] = 0
    else:
      distribution[c] = distribution[c] + 1
  stats['character_distribution'] = distribution

def check_is_protein(record, stats):
  for c in record.seq:
    protein_alphabet = {''}

def check_has_X(record, stats):
  '''Checks if the sequence contains 'X'. This character denotes ambiguity in
  amino acid sequences. In blastp these characters won't be aligned.'''
  _process_character(stats, 'X', 'Xs')

def check_has_star(record, stats):
  '''Checks if the sequence contains '*' and increases the counter. 
  This usually denotes stop codons in amino acid sequences. Some tools 
  cannot handle these sequnces'''
  _process_character(stats, '*', 'stars')

def check_has_dot(record, stats):
  '''Checks if the sequence contains '.' and increases the counter. 
  While unproblematic for DNA and RNA sequences this may cause problems for proteins.
  The diamond aligner cannot handle sequences with dots'''
  _process_character(stats, '.', 'dots')

def check_has_dash(record, stats):
  '''Checks if the sequence contains '-'. This character is used in alignments saved
  in fasta files.'''
  _process_character(stats, '-', 'dash')

def check_has_N(record, stats):
  '''Checks if the sequence contains 'N'. This character denotes ambiguity in
  DNA sequences.'''
  _process_character(stats, 'N', 'Ns')

dna_unambiguous_alphabet = ['A','T','G','C']
dna_ambiguous_alphabet = ['A','T','G','C','N','R','Y','S','W','K','M','B','D','H','V','.','-']
rna_unambiguous_alphabet = ['A','U','G','C']
rna_ambiguous_alphabet = ['A','U','G','C','N','R','Y','S','W','K','M','B','D','H','V','.','-']
aa_unambiguous_alphabet=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']
# source https://www.ddbj.nig.ac.jp/ddbj/code-e.html
aa_ambiguous_alphabet=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','*','_','B','Z','J']

allchars = set(dna_unambiguous_alphabet + dna_ambiguous_alphabet + rna_unambiguous_alphabet + rna_ambiguous_alphabet + aa_unambiguous_alphabet + aa_ambiguous_alphabet)

##############
# is this really needed?

alphabets = {
  'dna_unambiguous': dna_unambiguous_alphabet,
  'dna_ambiguous': dna_ambiguous_alphabet,
  'rna_unambiguous': rna_unambiguous_alphabet,
  'rna_ambiguous': rna_ambiguous_alphabet,
  'aa_unambiguous': aa_unambiguous_alphabet,
  'aa_ambiguous': aa_ambiguous_alphabet,
}
################

def check_sequence_type(record, stats):
  c_dist = stats['character_distribution']
  if not 'sequence_types' in stats:
    stats['sequence_types'] = {
      'dna_unambiguous': 0,
      'rna_unambiguous': 0,
      'dna_ambiguous': 0,
      'rna_ambiguous': 0,
      'aminoacid_unambiguous': 0,
      'aminoacid_ambiguous': 0,
      'other': 0,
    }
  sequence_types = stats['sequence_types']
    
  if not 'unhandled_characters' in stats:
    stats['unhandled_characters'] = {
    }
  unhandled_characters = stats['unhandled_characters']

  if contains_only(c_dist, dna_unambiguous_alphabet):
    sequence_types['dna_unambiguous'] = sequence_types['dna_unambiguous'] + 1
  elif contains_only(c_dist, rna_unambiguous_alphabet):
    sequence_types['rna_unambiguous'] = sequence_types['rna_unambiguous'] + 1
  elif contains_only(c_dist, dna_ambiguous_alphabet):
    sequence_types['dna_ambiguous'] = sequence_types['dna_ambiguous'] + 1
  elif contains_only(c_dist, rna_ambiguous_alphabet):
    sequence_types['rna_ambiguous'] = sequence_types['rna_ambiguous'] + 1
  elif contains_only(c_dist, aa_unambiguous_alphabet):
    sequence_types['aminoacid_unambiguous'] = sequence_types['aminoacid_unambiguous'] + 1
  elif contains_only(c_dist, aa_ambiguous_alphabet):
    sequence_types['aminoacid_ambiguous'] = sequence_types['aminoacid_ambiguous'] + 1
  else:
    sequence_types['other'] = sequence_types['other'] + 1
    identify_unknown_character(unhandled_characters, c_dist)

def identify_unknown_character(unhandled_characters, character_distribution):
  characters = character_distribution.keys()
  for c in characters:
    if c not in allchars:
      if c not in unhandled_characters:
        unhandled_characters[c] = 0
      unhandled_characters[c] = unhandled_characters[c] + 1

def contains_only(dist, characters):
  remaining = len(dist)
  for c in characters:
    if c in dist:
      remaining = remaining - 1
  return remaining == 0

def _process_character(stats, character, field):
  if character in stats['character_distribution']:
    if field not in stats:
      stats[field] = 0
    stats[field] = stats[field] + 1

if __name__ == "__main__":
    main()
