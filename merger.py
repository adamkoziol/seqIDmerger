#!/usr/bin/env python
import os
import sys
import errno
from glob import glob
__author__ = 'adamkoziol'

# Initialise the count
count = 0


class GenObject(object):
    """Object to store static variables"""
    def __init__(self, x=None):
        start = (lambda y: y if y else {})(x)
        super(GenObject, self).__setattr__('datastore', start)

    def __getattr__(self, key):
        return self.datastore[key]

    def __setattr__(self, key, value):
        if value:
            self.datastore[key] = value
        else:
            self.datastore[key] = "NA"


class MetadataObject(object):
    """Object to store static variables"""
    def __init__(self):
        """Create datastore attr with empty dict"""
        super(MetadataObject, self).__setattr__('datastore', {})

    def __getattr__(self, key):
        """:key is retrieved from datastore if exists, for nested attr recursively :self.__setattr__"""
        if key not in self.datastore:
            self.__setattr__(key)
        return self.datastore[key]

    def __setattr__(self, key, value=GenObject(), **args):
        """Add :value to :key in datastore or create GenObject for nested attr"""
        if args:
            self.datastore[key].value = args
        else:
            self.datastore[key] = value

    def dump(self):
        """Print only the nested dictionary values; removes __methods__ and __members__ attributes"""
        metadata = {}
        for attr in self.datastore:
            metadata[attr] = {}
            if not attr.startswith('__'):
                if isinstance(self.datastore[attr], str):
                    metadata[attr] = self.datastore[attr]
                else:
                    metadata[attr] = self.datastore[attr].datastore
        return metadata


def make_path(inpath):
    """
    from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL
    :param inpath: string of the supplied path
    """
    import os
    try:
        # os.makedirs makes parental folders as required
        os.makedirs(inpath)
    # Except os errors
    except OSError as exception:
        # If the os error is anything but directory exists, then raise
        if exception.errno != errno.EEXIST:
            raise


class Merger(object):

    def idseek(self):
        # Create a list of all the lines in the file: open(self.idfile).readlines()
        # Create a lambda function
        # Map the list to the lambda function and split the list based on the delimiter: x.split(self.delimiter)
        # List comprehension of individual seq IDs without whitespace: [x.rstrip() for x in ...]
        self.seqids = map(lambda x: [x.rstrip() for x in x.split(self.delimiter)], open(self.idfile).readlines())

    def idfind(self):
        """Find the fastq files associated with the seq IDs pulled from the seq ID file. Populate a MetadataObject
        with the name of the merged files as well as the fastq file names and paths"""
        for merge in self.seqids:
            # Initialise the lists
            seqidfile = []
            outputname = []
            for ids in merge:
                # Ensure that the id exists. Dues to the way the ids were pulled from the file, newline characters
                # will be entered into the list. Skip them
                if ids:
                    # Glob for files in the path with the seq ID and 'fastq'
                    idfile = glob('{}{}*fastq*'.format(self.path, ids))
                    # Assertion to ensure that all the files specified in :self.idfile are present in the path
                    assert idfile, 'Cannot find files for seq ID: {}. Please check that the seqIDs ' \
                                   'provided in the seq ID file match the files present in the path'.format(ids)
                    # Append the fastq file and path and the seq ID to the appropriate list
                    seqidfile.append(idfile)
                    outputname.append(ids)
            # Create the MetadataObject to store variables associated with each file set to merge
            data = MetadataObject()
            # Set the name of the file set
            data.name = '-'.join(outputname)
            # Create the general category for the MetadataObject
            data.general = GenObject()
            # Populate the fastq files in the file set
            data.general.fastqfiles = seqidfile
            # Append the MetadataObjects to a list
            self.data.append(data)

    def idmerge(self):
        """Merge the files together"""
        from threading import Thread
        # As :self.data.general.fastqfiles is a list of lists, find the length of each child list with len(), add
        # these lengths together with sum(), and use range() to set the number of times to iterate
        for i in range([sum(len(sample.general.fastqfiles) for sample in self.data)][0]):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.merge, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.data:
            # Initialise strings to hold the forward and reverse fastq files
            forwardfiles = []
            reversefiles = []
            # Create the output directory
            sample.general.outputdir = '{}{}'.format(self.path, sample.name)
            make_path(sample.general.outputdir)
            # Iterate through the samples
            for files in sample.general.fastqfiles:
                # Find the forward and reverse files (forward files must have have either _R1_ or _1_
                forwardfiles.append((map(lambda x: x if '_R1_' in x or '_1_' in x else '', sorted(files)))[0])
                reversefiles.append((map(lambda x: x if '_R2_' in x or '_2_' in x else '', sorted(files)))[1])
            # Add the files to the processing queue
            sample.general.outputforward = '{}/{}_S1_L001_R1_001.fastq.gz'.format(sample.general.outputdir, sample.name)
            sample.general.outputreverse = '{}/{}_S1_L001_R2_001.fastq.gz'.format(sample.general.outputdir, sample.name)
            # Add the command object to self.data
            sample.commands = GenObject()
            sample.commands.forwardmerge = 'cat {} > {}'.format(' '.join(forwardfiles), sample.general.outputforward)
            sample.commands.reversemerge = 'cat {} > {}'.format(' '.join(reversefiles), sample.general.outputreverse)
            self.mergequeue.put((sample.commands.forwardmerge, sample.general.outputforward))
            self.mergequeue.put((sample.commands.reversemerge, sample.general.outputreverse))
        # Join the threads
        self.mergequeue.join()

    def merge(self):
        while True:  # while daemon
            # Unpack the merge command and the output file from the queue
            (mergecommand, outputfile) = self.mergequeue.get()
            # Don't run the command if the output file exists
            if not os.path.isfile(outputfile):
                try:
                    self.execute(mergecommand)
                except KeyboardInterrupt:
                    self.printtime(u'Keyboard interrupt! The system call will not stop until it is finished.')
                    self.mergequeue.empty()
                    try:
                        os.remove(outputfile)
                    except IOError:
                        pass
                    sys.exit()
            # Signal to mergequeue that job is done
            self.mergequeue.task_done()

    def filelink(self):
        # If the creation of a sample sheet is necessary
        if self.samplesheet:
            # Extract the path of the current script from the full path + file name
            samplesheet = open('{}/SampleSheet.csv'.format(os.path.split(os.path.abspath(__file__))[0])).readlines()
            # Iterate through each merged file
            for sample in self.data:
                # Append enough information to the list to allow the pipeline to work
                samplesheet.append('{},{},,,NA,NA,NA,NA,NA,NA\n'.format(sample.name, sample.name))
            # Initialise the name and path of the output sample sheet
            outsamplesheet = '{}/SampleSheet.csv'.format(self.assemblypath)
            # Don't overwrite a sample sheet already present in the directory
            if not os.path.isfile(outsamplesheet):
                # Open the file to write and write to it
                with open(outsamplesheet, 'wb') as writesheet:
                    writesheet.write(''.join(samplesheet))
        # Link the files to the assembly path
        for sample in self.data:
            try:
                # Link the files
                os.symlink(sample.general.outputforward, '{}/{}'.format(self.assemblypath,
                                                                        os.path.basename(sample.general.outputforward)))
                os.symlink(sample.general.outputreverse, '{}/{}'.format(self.assemblypath,
                                                                        os.path.basename(sample.general.outputreverse)))
            # Except os errors
            except OSError as exception:
                # If the os error is anything but directory exists, then raise
                if exception.errno != errno.EEXIST:
                    raise

    def execute(self, command, outfile=""):
        """
        Allows for dots to be printed to the terminal while waiting for a long system call to run
        :param command: the command to be executed
        :param outfile: optional string of an output file
        from https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
        """
        import time
        from subprocess import Popen, PIPE, STDOUT
        # Initialise the starting time
        start = int(time.time())
        maxtime = 0
        # Run the commands - direct stdout to PIPE and stderr to stdout
        process = Popen(command, shell=True, stdout=PIPE, stderr=STDOUT)
        # Create the output file - if not provided, then nothing should happen
        writeout = open(outfile, "ab+") if outfile else ""
        # Poll process for new output until finished
        while True:
            # If an output file name is provided
            if outfile:
                # Get stdout into a variable
                nextline = process.stdout.readline()
                # Print stdout to the file
                writeout.write(nextline)
            # Break from the loop if the command is finished
            if process.poll() is not None:
                break
            # Adding sleep commands slowed down this method when there was lots of output. Difference between the start
            # time of the analysis and the current time. Action on each second passed
            currenttime = int(time.time())
            # As each thread will print a dot at the same time, often the dots printed to the terminal do not look
            # even. Instead of 80 per line, there are sometimes around 78-82, or just one. Having this random number
            # seems to fix this
            from random import randint
            # Set the number to be a random integer between 0 and 999
            number = randint(0, 999)
            if currenttime - start > maxtime + number:
                # Set the max time for each iteration
                maxtime = currenttime - start
                # Print up to 80 dots on a line, with a one second delay between each dot
                if self.count <= 80:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    self.count += 1
                # Once there are 80 dots on a line, start a new line with the the time
                else:
                    sys.stdout.write('\n.')
                    sys.stdout.flush()
                    self.count = 1
        # Close the output file
        writeout.close() if outfile else ""

    def __init__(self, args, start):
        """
        :param args: list of arguments passed to the script
        :param start: the start time
        Initialises the variables required for this class
        """
        from Queue import Queue
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args['path'], "")
        self.start = start
        # Determine which seq ID file to use
        # If an argument for the file is provided, use it
        if args['f']:
            self.idfile = args['f']
            # If there is no path information present in the argument, use :path + the file name
            if '/' not in self.idfile:
                self.idfile = '{}{}'.format(self.path, self.idfile)
                assert os.path.isfile(self.idfile), 'Could not find the seq ID file. Please double check the supplied' \
                                                    'arguments'
        # If no argument is provided, try to find the file
        else:
            # Look for .txt, .tsv, or .csv files
            self.idfile = map(lambda x: glob('{}*{}'.format(self.path, x)), ['.txt', '.csv', '.tsv'])
            # Initialise the file count
            filecount = 0
            # Iterate through each extension
            for extension in self.idfile:
                if extension:
                    # If a single file with that extension was found, set :self.idfile to that file
                    extensiontype = extension[0].split('.')[1]
                    assert len(extension) == 1, u'Too many .{} entries found for the ID file.'.format(extensiontype)
                    self.idfile = extension[0]
                    # Increment the count
                    filecount += 1
            # Assertions to exit if there are too many or too few potential seq ID files
            assert filecount <= 1, u'Too many potential ID files found. Please check that there is only one .txt,' \
                                   u' .tsv, or .csv file in the path'
            assert filecount >= 1, u'Could not find a seq ID file with a .txt, .tsv, or .csv extension in the path'
        # Assertion to ensure that the seq ID file exists
        assert os.path.isfile(self.idfile), u'seqID file cannot be found {0!r:s}'.format(self.idfile)
        self.printtime(u'Using {} as the file containing seq IDs to be merged'.format(self.idfile))

        # Set the delimiter
        self.delimiter = args['d'].lower()
        if self.delimiter == 'space':
            self.delimiter = ' '
        elif self.delimiter == 'tab':
            self.delimiter = '\t'
        elif self.delimiter == 'comma' or self.delimiter == ',':
            self.delimiter = ','
        # Initialise class variables
        self.seqids = ""
        self.seqfiles = []
        self.data = []
        self.mergequeue = Queue()
        self.count = 0
        # Find which IDs need to be merged together from the text file
        self.idseek()
        # Find the files corresponding to the IDs
        self.idfind()
        # Merge the files together
        self.idmerge()
        # Exit
        self.printtime(u'Files have been successfully merged.')
        if args['linkFiles']:
            # Create the assembly folder and path from the supplied arguments
            self.assemblyfolder = args['o'] if args['o'] else self.path.split('/')[-2]
            self.assemblypath = os.path.join(args['a'], "") + self.assemblyfolder
            make_path(self.assemblypath)
            assert os.path.isdir(self.assemblypath), 'Could not create the destination folder. Please double-check ' \
                                                     'your supplied arguments'
            self.samplesheet = args['samplesheet']
            # Run the linking
            self.filelink()
            self.printtime(u'Files have been successfully linked to the assembly folder.')
        sys.exit()

    def printtime(self, string):
        """Prints a string in bold with the elapsed time
        :param string: a string to be printed in bold
        """
        import time
        print(u'\n\033[1m' + u'[Elapsed Time: {:.2f} seconds] {}'.format(time.time() - self.start, string) + u'\033[0m')

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    from time import time
    from argparse import ArgumentParser
    import subprocess
    # Extract the path of the current script from the full path + file name
    homepath = os.path.split(os.path.abspath(__file__))[0]
    # Find the commit of the script by running a command to change to the directory containing the script and run
    # a git command to return the short version of the commit hash
    commit = subprocess.Popen('cd {} && git rev-parse --short HEAD'.format(homepath),
                              shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()
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
    parser.add_argument('-l', '--linkFiles', action='store_true', help='Optionally link the files to the \'WGS_Spades\''
                        ' directory. Note that this is specific to the local setup here and is not recommended unless'
                        ' your set-up is similar')
    parser.add_argument('-a', metavar='assemblyLocation', default='/nas/akoziol/WGS_Spades/', help='Path to a folder'
                        ' where files are automatically assembled using a cluster. Only relevant if linking the files')
    parser.add_argument('-o', metavar='outputdirectory', help='A directory name to use when linking the merged .fastq'
                        'files to the WGS_Spades folder e.g. 2016-01-19_ListeriaMerged. If this is not provided, then'
                        'the program will use the name of lowest folder in the path e.g. \'files\' will be used if the'
                        'path is \'/path/to/files\'')
    parser.add_argument('-s', '--samplesheet', action='store_true', help='Depending on the version of the assembly '
                        'pipeline, a sample sheet is required. Including this option will populate a basic sample sheet'
                        'with enough information in order to allow the pipeline to proceed')

    # Get the arguments into a list
    arguments = vars(parser.parse_args())

    starttime = time()
    # Run the pipeline
    output = Merger(arguments, starttime)
