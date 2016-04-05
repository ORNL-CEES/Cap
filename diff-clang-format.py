#!/usr/bin/env python

'''
Usage:
  diff-clang-format.py [--file-extension=<arg>...] [options] <path>...

Option:
  -h --help                    Show this screen
  -q --quiet                   Do not print the diff
  --file-extension=<arg>       Filename extension with a dot [default: .hpp .cpp]
  --style=<style>              Coding style supportted by clang-format [default: LLVM]
  --configuration-file=<file>  Style configuation .clang-format file
'''

from docopt import docopt
from subprocess import Popen, PIPE
from os import walk, getcwd, remove
from shutil import copy

def diff_with_formatted_source(original_file, style):
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
    ostream = b''
    if stdout:
        headers = ''
        headers += '<<< ' + original_file + '\n'
        headers += '--- diff ---\n'
        headers += '>>>' + formatted_file + '\n'
        ostream += headers.encode('utf-8')
        ostream += stdout
    remove(formatted_file)
    return ostream

def run(paths, file_extensions, style, config):
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
                if (any([file.endswith(extension) for extension in file_extensions])):
                    diff = diff_with_formatted_source(root+'/'+file, style)
                    if diff:
                        diffs[root+'/'+file] = diff
    return diffs

if __name__ == '__main__':
    args = docopt(__doc__)

    paths = args['<path>']
    extensions = args['--file-extension']
    style = args['--style']
    config = args['--configuration-file']

    diffs = run(paths, extensions, style, config)
    if diffs:
        if not args['--quiet']:
            for file, diff in diffs.items():
                print('####', file, '####')
                print(diff.decode('utf-8'))
        print('Bad format')
        exit(1)
    else:
        print('OK')
        exit(0)
