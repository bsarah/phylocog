Program description
-------

This is a commandline tool written in Python, that accepts a specific domainmapper output file with sequence information as input. It transforms the input format and aligns all given sequences pairwise using a modified Needleman-Wunsch algorithm.
It generates an output files with the specified filename at the end, containing the alignments.  

How to use
-------

To use the program you have to provide an input file in the same directory. Example of use:

```sh
python programname.py inputfile.txt
```

Output format
-------

### intermediate file
This document is used as input for the modified Needleman-Wunsch algorithm. It contains the written out sequence information from the provided file.

### output file 
The output file contains all alignments and the information about the sequences in the following representation: 
> \>ClusID of seq1,accession of seq1,>ClusID of seq2,accession of seq2, alignmentscore, length of alignment
> Alignment

Options
-------

| Option | Usage |
| ------ | ------ |
| filenamne | you need to provide an input file |
| -h, -\-help | displays help |
| -v, -\-verbose | writes a second output file that is more human-readable |
| -s, -\-noselfalignment| deactivates the computation of alignments like (a,a) |
| -l, -\-log | sets the log level to [DEBUG, INFO, WARNING, ERROR, CRITICAL] |
| -t, -\-temp | keeps the intermediate file |
| -g, -\-gapextension | different penalties for gap openings and gap extensions |

Set the logging level option such as -\-log=INFO.

License
-------

The code is released under the MIT license. See [LICENSE](https://github.com/blaueste/bachelor_thesis/blob/main/LICENSE).
