Seq ID Merger
==============
# Introduction

Merges .fastq.gz files based on an input list of seq IDs.

# Installation

From the command line:
git clone https://github.com/adamkoziol/seqIDmerger.git

Ideally add the scripts to your path. If you don't know how, that's fine. Just run the merger.py script from within 
the seqIDmerger folder

# Requirements
- Unix-like environment
- List of seq IDs: Ideally, list should be a .txt file with each row containing the seq IDs of the files to be merged, separated by whitespace
e.g.
2013-SEQ-0073 2014-SEQ-0625
2014-SEQ-0029 2014-SEQ-0033 2014-SEQ-0074
- The .fastq files corresponding to the seq IDs in the list


# Arguments

path
The path of the folder containing the .fastq(.gz) files and the list of seq IDs
-f 
The name and path of the file containing seqIDs to merge and reassemble. If this file is in the path, then including the path is not necessary
for this argument. Alternatively, as long as the file has a .txt, .csv, or. tsv file extension, you can omit this argument altogether. Note: if you don't supply the argument, and
there are multiple files with any of these extensions, the program will fail
-d
The delimiter used to separate seqIDs. Popular options are "space", "tab", and "comma". Default is space. Note: you can use custom delimiters. 
Just be aware that a delimiter, such as "-" will break the program if there are hyphens in your sample names
-l
Optionally link the files to the WGS_Spades directory. Note that this is specific to the local setup here and is not recommended unless your set-up is similar
-a
Path to a folder where files are automatically assembled using a cluster. Only relevant if linking the files. Default is /nas/akoziol/WGS_Spades/
-o
A directory name to use when linking the merged .fastq files to the WGS_Spades folder e.g. 2016-01-19_ListeriaMerged. If this is not provided, then
the program will use the name of lowest folder in the path e.g. 'files' will be used if the path is '/path/to/files'
-s
Depending on the version of the assembly pipeline, a sample sheet is required. Including this option will populate a basic sample sheet
with enough information in order to allow the pipeline to proceed
-c
Copies rather than symbolically linking the files to the destination folder

# Running 

Example command to run the script

- merger.py /nas0/Bioinformatic_Requests/6346_Ashley_20160118 -f IDs.txt -l -s
    - merger.py is in $PATH
    - path argument: /nas0/Bioinformatic_Requests/6346_Ashley_20160118
    - file of IDs: IDs.txt
    - link the files: -l
    - create a sample sheet: -s
    