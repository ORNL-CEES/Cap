#!/usr/bin/env python

'''
Usage:
  diff-clang-format.py [--file-extension=<arg>...] [options] <path>...

Option:
  -h --help                    Show this screen
  -q --quiet                   Do not print the diff
  --file-extension=<arg>       Filename extension with a dot [default: .hpp .cpp]
  --style=<style>              Coding style supported by clang-format [default: LLVM]
  --configuration-file=<file>  Style configuation .clang-format file
  --apply-patch                Apply diff patch to the original source
'''

from docopt import docopt
from subprocess import Popen, PIPE
from os import walk, getcwd, remove, rename
from shutil import copy

def diff_with_formatted_source(original_file, style, patch):
    '''Compute the diff between the C++ source code and its formatted version
    using clang-format.

    Parameters
    ----------
    original_file : str
        Absolute path to the file to be formatted with clang-format.
    style : str
        Coding style.

    Returns
    -------
    bytes
        Output of diff.
    '''
    formatted_file = getcwd()+'/'+'.tmp'
    cmd = ['clang-format', '-style='+style, original_file]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode or stderr:
        stdout, stderr = p.communicate()
        raise RuntimeError('clang-format failed')
    with open(formatted_file, 'wb') as fout:
        fout.write(stdout)
    cmd = ['diff', formatted_file, original_file]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    if stderr:
        raise RuntimeError('diff failed')
    if patch:
        rename(formatted_file, original_file)
    else:
        remove(formatted_file)
    return stdout

def run(paths, file_extensions, style, config, patch):
    '''
    Parameters
    ----------
    paths : list of str
    file_extensions : list of str
    style : str
    config : str or None

    Returns
    -------
    dict
        Gather output of all the diffs.
    '''
    if config:
        copy(config, getcwd()+'/'+'.clang-format')

    diffs = {}
    for path in paths:
        for root, dirs, files in walk(path):
            for file in files:
                file_path = root+'/'+file
                if (any([file.endswith(extension) for extension in file_extensions])):
                    diff = diff_with_formatted_source(file_path, style, patch)
                    if diff:
                        diffs[file_path] = diff
    return diffs

if __name__ == '__main__':
    args = docopt(__doc__)

    paths = args['<path>']
    extensions = args['--file-extension']
    style = args['--style']
    config = args['--configuration-file']
    patch = args['--apply-patch']

    diffs = run(paths, extensions, style, config, patch)
    if diffs:
        if not args['--quiet']:
            for file, diff in diffs.items():
                print('####', file, '####')
                print(diff.decode('utf-8'))
        print('Bad format')
    else:
        print('OK')
    reformatted_files = len(diffs)
    if patch:
        exit(0)
    else:
        exit(reformatted_files)
