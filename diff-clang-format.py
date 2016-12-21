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
  --binary=<path>              Location of binary to use for clang-format [default: clang-format]
'''

from docopt import docopt
from subprocess import Popen, PIPE
from os import walk, chdir, getcwd, remove, path
from sys import exc_info
from shutil import copy, rmtree, move
from tempfile import mkdtemp

def diff_with_formatted_source(original_file, command, patch):
    '''Compute the diff between the C++ source code and its formatted version
    using clang-format.

    Parameters
    ----------
    original_file : str
        Absolute path to the file to be formatted with clang-format.
    command : list of str
        Clang-format executable with options.
    patch : bool
        Apply the diff patch to the original source file.

    Returns
    -------
    bytes
        Output of diff.
    '''
    formatted_file = getcwd() + '/' + path.basename(original_file)
    cmd = command[:]
    cmd.append(original_file)
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode or stderr:
        print(cmd)
        print(stderr.decode('utf-8'))
        raise RuntimeError('clang-format failed')
    with open(formatted_file, 'wb') as fout:
        fout.write(stdout)
    cmd = ['diff', '-u',
        original_file, '--label', 'a' + original_file,
        formatted_file, '--label', 'b' + original_file]
    p = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    if p.returncode > 1 or stderr:
        print(cmd)
        print(stderr.decode('utf-8'))
        raise RuntimeError('diff failed')
    if patch and stdout:
        move(formatted_file, original_file)
    else:
        remove(formatted_file)
    return stdout

def run(paths, file_extensions, command, patch):
    '''Search recursively for files ending with extension and check the format
    of the source code.

    Parameters
    ----------
    paths : list of str
        Search paths.
    file_extensions : list of str
        File extensions with a dot.
    command : list of str
        Clang-format executable with options.
    patch : bool
        Apply the diff patch to the original source files.

    Returns
    -------
    dict
        Gather output of all the diffs.
    '''
    diffs = {}
    for path in paths:
        for root, dirs, files in walk(path):
            for file in files:
                file_path = root+'/'+file
                if (any(file.endswith(extension) for extension in file_extensions)):
                    diff = diff_with_formatted_source(file_path, command, patch)
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
    command = [args['--binary']]
    if style:
      command.append('-style='+style)
    current_directory = getcwd()
    temporary_directory = mkdtemp()
    if config:
        copy(config, temporary_directory + '/' + '.clang-format')
    chdir(temporary_directory)
    try:
        diffs = run(paths, extensions, command, patch)
    except:
        e = exc_info()[0]
        print("<p>Error: %s</p>" % e)
        raise
    chdir(current_directory)
    rmtree(temporary_directory);
    if diffs:
        if not args['--quiet']:
            for file, diff in diffs.items():
                print(diff.decode('utf-8'))
        print('{0} file(s) not formatted properly:'.format(len(diffs)))
        for key in diffs.keys():
            print('    {0}'.format(key))
    else:
        print('OK')
    if patch:
        exit(0)
    else:
        reformatted_files = len(diffs)
        exit(reformatted_files)
