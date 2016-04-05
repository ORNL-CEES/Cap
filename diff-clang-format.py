#!/usr/bin/env python

'''
Usage: diff-clang-format.py [--file-extension=<arg>]... [PATHS ...]

Options:
  -h --help               Show this screen
  --file-extension=<arg>  File extensions [default: .hpp .cpp]
'''

from docopt import docopt
from subprocess import Popen, PIPE
from os import walk, getcwd, remove

def print_diff_with_formatted_source(original_file):
    '''Print the diff between the C++ source code and its formatted version
    using clang-format.
    Parameters
    ----------
    original_file : str
        Absolute path to the file to be formatted with clang-format.
    '''
    formatted_file = getcwd()+'/'+'.tmp'
    cmd = ['clang-format', '-style=file', original_file]
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

if __name__ == '__main__':
    args = docopt(__doc__)
    for path in args['PATHS']:
        for root, dirs, files in walk(path):
            for file in files:
                if (any([file.endswith(extension) for extension in args['--file-extension']])):
                    ostream = print_diff_with_formatted_source(root+'/'+file)
                    if ostream:
                        print(ostream.decode('utf-8'))
