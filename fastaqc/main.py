#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import os
import logging
from dbxref import resolver, config
from pbr.version import VersionInfo
import json

__version__ = VersionInfo('fastaqc').semantic_version().release_string()

def main():
    parser = argparse.ArgumentParser(description='Version ' + __version__ + '\nValidate fasta file', formatter_class=RawTextHelpFormatter)
    parser.set_defaults(func=help)

    subparsers = parser.add_subparsers()
    info_parser = subparsers.add_parser('info')
    info_parser.add_argument('fasta')
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
    print ('fastaqc Version ' + __version__)
    print ('')
    print ('Supported dbxref databases:')
    providers = config.load_providers()
    for provider in providers:
      print ('   ' + provider['name'])
      print ('     Prefixes: ' + str.join(', ', [x for x in provider['prefixes']]))
      print ('     Formats : ' + str.join(', ', [x for x in provider['resources']]))

if __name__ == "__main__":
    main()
