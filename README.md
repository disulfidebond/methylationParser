# methylationParser
The software package MethylationParser and Formatter (MAAF), is used to parse and format genomic methylation data from Bismark.  It uses python 2.7, and is capable of running on any computer or server, however, using a computer with at least 16 GB of RAM is strongly recommended for human genome samples.  MAAF takes as input an output coverage file from Bismark, and requires an annotated gene file in GTF format to be in the current working directory named “genes.gtf”, then filters methylated regions that meet at least at coverage of 10, and greater than 50% methylation.  Next, it attempts to match the methylated regions to annotated genes, and finally outputs a comma separated list of gene name with each methylated site.  Alternatively, a specific pathway or list of genes can be provided, and in lieu of the entire transcriptome, MAAF will narrow down the search to the specified genes.  MAAF usually takes between 4-18 hours to process a sample, depending on the sample size, amount of genes to query, and amount of available RAM on the computer.

Notes for README and Usage
import bismarck coverage file
important to initially import and store in list if using argparse, optionally parse and store
if -c argument is sortcov: 
run parseMethylatedData()
strip then split line on ‘\t’
if coverage is > 10 and methylation > 50% adds to list, then returns list
parses genes.gtf file into tuple
scans methylation list and genes.gtf tuple to match methylation regions to genes, then returns tuple
formats data with format [geneName: name, geneCDS: [], Methylation: [{methylationEntry1}, {methylationEntry2}, {methylationEntryN}]}]
prints data
if -c argument is sortmethyl:
run parseFormattedDataAsList() to parse formatted data from format step above
prints data
if -c argument is sortmethylP, with -f geneList argument:
runs steps from sortcov with added step of parsing a gene list for a pathway and matching genes in the pathway to the genes.gtf file to use to match to methylation regions
output currently is parsed with a manual step by searching for a “####” line to find output of genes, then removing newlines with text edit, then running parseMethylationFile.py to get a count of the genes with methylation sites for each sample
