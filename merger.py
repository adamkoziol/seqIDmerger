#!/usr/bin/env python
import os

__author__ = 'adamkoziol'


class Merger(object):

    def idmerge(self):
        pass

    def __init__(self, args, start):
        """
        :param args: list of arguments passed to the script
        Initialises the variables required for this class
        """
        import sys
        from glob import glob
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args['path'], "")
        self.start = start
        if args['f']:
            self.idfile = args['f']
            if "/" not in self.idfile:
                self.idfile = u'{}{}'.format(self.path, self.idfile)
        else:
            self.idfile = map(lambda x: glob('{}*{}'.format(self.path, x)), ['.txt', '.csv', '.tsv'])
            filecount = 0
            for extension in self.idfile:
                if extension:
                    if len(extension) == 1:
                        self.idfile = extension[0]
                        filecount += 1
                if filecount == 0:
                    print u'Could not find a seq ID file with a .txt, .tsv, or .csv extension in the path'
                    sys.exit()
                elif filecount > 1:
                    print(u'Too many entries found for the ID file. Please check that there is only one'
                          u'.txt, .tsv, or .csv file in the path')
                    sys.exit()
        assert os.path.isfile(self.idfile), u'seqID file cannot be found {0!r:s}'.format(self.idfile)
        self.delimiter = args['d'].lower()
        if self.delimiter == 'space':
            self.delimiter = ' '
        elif self.delimiter == 'tab':
            self.delimiter = '\t'
        elif self.delimiter == 'comma' or self.delimiter == ',':
            self.delimiter = ','
        print u'{} {} {}'.format(self.path, self.idfile, self.delimiter)


# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    import subprocess
    from time import time
    # Get the current commit of the pipeline from git
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(__file__)[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
    from argparse import ArgumentParser
    # Parser for arguments
    parser = ArgumentParser(description='Merges seqIDs')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s commit {}'.format(commit))
    parser.add_argument('path',  help='Specify path')
    parser.add_argument('-f', metavar='idFile', help='The name and path of the file containing seqIDs to '
                        'merge and reassemble. If this file is in the path, then including the path is not necessary'
                        ' for this argument. Alternatively, as long as the file has a .txt, .csv, or. tsv file '
                        'extension, you can omit this argument altogether. Note: you don\'t supply the argument, and'
                        'there are multiple files with any of these extensions, the program will fail')
    parser.add_argument('-d', metavar='delimiter', default='space', help='The delimiter used to separate seqIDs. '
                        'Popular options are "space", "tab", and "comma". Default is space. Note: you can use custom'
                        'delimiters. Just be aware that a delimiter, such as "-" will break the program if there are'
                        'hyphens in your sample names')

    # Get the arguments into a list
    arguments = vars(parser.parse_args())

    starttime = time()
    # Run the pipeline
    output = Merger(arguments, starttime)
