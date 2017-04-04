#!/usr/bin/python

# Author: Chris Baker
# Email: ccmbaker@princeton.edu

# The original version of this script was written while the author was in the
# Pierce Lab, Department of Organismic and Evolutionary Biology, Harvard University

# This update, which allows the script to work with NCBI accession.version numbers instead
# of GI numbers, was written while the author was in the Pringle and Tarnita Labs,
# Department of Ecology and Evolutionary Biology, Princeton University

# Last modified: 4 April 2017

#################################################

import os, argparse

from textwrap import dedent
from time import strftime
from time import localtime
from sys import exit

from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles
from cogent.parse.ncbi_taxonomy import NcbiTaxonNode
from cogent.core.tree import TreeNode

#################################################

# parse command line options

parser = argparse.ArgumentParser(description='Generates QIIME-compatible files from a suitable FASTA file or list of accession numbers, plus the NCBI taxonomy database.')

inputfilepath_args = parser.add_mutually_exclusive_group(required=True)

# Path including filename for input FASTA file. Expects deflines that look like ">XX######.# ..." where "XX######.# " is the NCBI accession.version number that uniquely identifies the sequence, and "..." is arbitrary text. Does NOT work with other defline formats, including the old deflines with the format ">gi|######|...". For files with these old-style deflines see v.1.0 of this script.
inputfilepath_args.add_argument('-i','--inputfasta', 
    metavar='path/file.ext', 
    dest='infile_fasta_path', 
    action='store', 
    help='The path, including filename, of an input FASTA file. Expects deflines that look like ">XX######.# ..." where "XX######.# " is the NCBI accession.version number that uniquely identifies the sequence, and "..." is arbitrary text. Does NOT work with other defline formats, including the old deflines with the format ">gi|######|...". [REQUIRED, unless an accession number list is supplied with -L or --List. If both are supplied, the script fails.]')

# The path, including filename, of an input accession number list as an alternative to an input FASTA file.
inputfilepath_args.add_argument('-L','--List', 
    metavar='path/file.ext', 
    dest='infile_list_path', 
    action='store', 
    help='The path, including filename, of an input accession number list (one accession number per line) as an alternative to supplying an input FASTA file. [REQUIRED, unless a FASTA file is supplied with -i or --inputfasta. If both are supplied, the script fails.]')

# Path to output file.
parser.add_argument('-o','--outputfile', 
    metavar='path/file.ext', 
    dest='outfile_path', 
    action='store', 
    help='The path, including filename, for the accession-number-to-taxonomy output file. Will be created if it does not exist. Defaults to same directory as input FASTA or list file, with file name derived from the input file\'s name, if not supplied here.')

# Path for output log file.
parser.add_argument('-g','--outputlog', 
    metavar='path/file.ext', 
    dest='logfile_path', 
    action='store', 
    help='The path, including filename, for the output log file. Will be created if it does not exist. Defaults to same directory as input FASTA or list file, with file name derived from the input file\'s name, if not supplied here.')

# Force overwrite if output file already exists?
parser.add_argument('-f','--force', 
    dest='force_overwrite', 
    action='store_true', 
    default=False, 
    help='If -f or --force option are passed, output file will overwrite any existing file with the same name. Without this option, output filename will be modified by the addition of a date/time string in the event that a file with the same name already exists.  [Default: False]')

# Filename and path for NCBI node-parent list (from taxdump.tar.gz)
parser.add_argument('-n','--nodes', 
    metavar='path/', 
    dest='ncbi_taxonomy_dir', 
    action='store', 
    default='./', 
    help='The directory where the files nodes.dmp, names.dmp, merged.dmp and delnodes.dmp are located. The files are typically downloaded from the NCBI in the compressed archive taxdump.tar.gz. [Default: ./]')

# Filename and path for NCBI accession number - taxid list (e.g. from nucl_gb.accession2taxid.gz)
parser.add_argument('-a','--acc2taxid', 
    metavar='path/file.ext', 
    dest='infile_acc2taxid_path', 
    action='store', 
    default='./nucl_gb.accession2taxid', 
    help='The path, including filename, of the (uncompressed) input accession number - taxid file, typically downloaded from the NCBI as nucl_gb.accession2taxid.gz and uncompressed to nucl_gb.accession2taxid. [Default: ./nucl_gb.accession2taxid]')

# lineage info to retain
parser.add_argument('-r','--ranks', 
    metavar='rank1,rank2,...', 
    dest='output_ranks', 
    action='store', 
    default='phylum,class,order,family,genus,species', 
    help=dedent('''\
        A comma-separated list (no spaces) of the taxonomic ranks you want to keep in the taxon names output. Ranks can be provided in any order, but names output will be generated in that order, so should typically be in descending order. Any taxonomic names assigned to one of the ranks in the list will be retained; other names will be discarded. If a taxon has no name for one of the ranks in the list, then NA will be inserted in the output. No sequences or taxa will be discarded -- this option only affects the names that are output for each taxon or sequence. The following ranks are available in the NCBI taxonomy database as at 2 June 2012 (in alphabetical order): 
            class, 
            family, 
            forma, 
            genus, 
            infraclass, 
            infraorder, 
            kingdom, 
            no_rank, 
            order, 
            parvorder, 
            phylum, 
            species, 
            species_group, 
            species_subgroup, 
            subclass, 
            subfamily, 
            subgenus, 
            subkingdom, 
            suborder, 
            subphylum, 
            subspecies, 
            subtribe, 
            superclass, 
            superfamily, 
            superkingdom, 
            superorder, 
            superphylum, 
            tribe, 
            varietas. 
        Note the use of underscores in three of the names. You should include these on the command line but they will be replaced by spaces when the script runs so that the names match the NCBI ranks. [Default: phylum,class,order,family,genus,species]'''))

args = parser.parse_args()

#################################################

def main():
    
    args.infile_nodesdmp_path, args.infile_namesdmp_path, args.infile_mergeddmp_path, \
        args.infile_delnodesdmp_path = input_files()
    
    args.outfile_path, args.logfile_path = output_files()
    
    args.output_ranks = args.output_ranks.replace('_',',').split(',')
    
    accessions_in_input_file = get_accessions_from_input_file()
    
    ncbi_full_taxonomy = NcbiTaxonomyFromFiles(open(args.infile_nodesdmp_path), open(args.infile_namesdmp_path))
    
    merged_taxids = get_merged_nodes()
    
    deleted_taxids = get_deleted_nodes()
        
    included_nodes, taxid, missing_accessions = obtain_nodes_for_each_accession(accessions_in_input_file, \
        args.infile_acc2taxid_path, ncbi_full_taxonomy, merged_taxids, deleted_taxids)
    
    taxid_taxonomy, missing_taxonomy = generate_taxonid_taxonomy(included_nodes, ncbi_full_taxonomy, args.output_ranks)
    
    generate_output_files(args.outfile_path, taxid_taxonomy, taxid, missing_accessions, missing_taxonomy)
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('Finished running entrez_qiime.py at ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()

#################################################

def input_files():
    
    if not os.path.isdir(args.ncbi_taxonomy_dir):
        exit('Error: path supplied via -n does not appear to be a valid directory (' + args.ncbi_taxonomy_dir + '))')
    
    infile_nodesdmp_path = os.path.join(args.ncbi_taxonomy_dir, 'nodes.dmp')
    
    if not os.path.isfile(infile_nodesdmp_path):
        exit('Error: cannot find file nodes.dmp at ' + infile_nodesdmp_path)
    
    infile_namesdmp_path = os.path.join(args.ncbi_taxonomy_dir, 'names.dmp')
    
    if not os.path.isfile(infile_namesdmp_path):
        exit('Error: cannot find file names.dmp at ' + infile_namesdmp_path)
    
    infile_mergeddmp_path = os.path.join(args.ncbi_taxonomy_dir, 'merged.dmp')
    
    if not os.path.isfile(infile_mergeddmp_path):
        exit('Error: cannot find file merged.dmp at ' + infile_mergeddmp_path)
    
    infile_delnodesdmp_path = os.path.join(args.ncbi_taxonomy_dir, 'delnodes.dmp')
    
    if not os.path.isfile(infile_delnodesdmp_path):
        exit('Error: cannot find file delnodes.dmp at ' + infile_delnodesdmp_path)
    
    if args.infile_fasta_path != None:
        if not os.path.isfile(args.infile_fasta_path):
            exit('Error: cannot find input FASTA file at ' + args.infile_fasta_path)
    elif args.infile_list_path != None:
        if not os.path.isfile(args.infile_list_path):
            exit('Error: cannot find input accession number list at ' + args.infile_list_path)
    else:
        print "Need to supply either FASTA file or accession number list as input."    # should never get here because either -i or -L is a required argument.
    
    if not os.path.isfile(args.infile_acc2taxid_path):
        exit('Error: cannot find input accession number - TaxonID file at ' + args.infile_acc2taxid_path)
    
    return(infile_nodesdmp_path, 
           infile_namesdmp_path, 
           infile_mergeddmp_path, 
           infile_delnodesdmp_path)

#################################################

def output_files():
    
    timestamp = strftime("%Y%m%d%H%M%S",localtime())
    
    # Construct output filenames if not supplied as arguments.
    
    if args.infile_fasta_path != None:
        inputfile_dir = os.path.dirname(args.infile_fasta_path)
        inputfile_filename = os.path.splitext(os.path.basename(args.infile_fasta_path))[0].rstrip(".")      # rstrip in case there are multiple periods at the end of the filename
    else:
        inputfile_dir = os.path.dirname(args.infile_list_path)
        inputfile_filename = os.path.splitext(os.path.basename(args.infile_list_path))[0].rstrip(".")       # rstrip in case there are multiple periods at the end of the filename

    if args.outfile_path == None:
        # if no output filepath supplied, save file to input directory and construct filename based on input filename
        outputfile_dir = inputfile_dir
        outputfile_filename = inputfile_filename + '_accession_taxonomy.txt'
        outfile_path = os.path.join(outputfile_dir, outputfile_filename)
    elif os.path.isdir(args.outfile_path):
        # if supplied output filepath is actually a directory, construct filename based on input filename
        outputfile_dir = args.outfile_path
        outputfile_filename = inputfile_filename + '_accession_taxonomy.txt'
        outfile_path = os.path.join(outputfile_dir, outputfile_filename)
    else:
        # if supplied output filepath is not an existing directory, then use it for output file
        outfile_path = args.outfile_path
    
    if args.logfile_path == None:
        # if no logfile filepath supplied, save logfile to input directory and construct filename based on input filename
        logfile_dir = inputfile_dir
        logfile_filename = inputfile_filename + '.log'
        logfile_path = os.path.join(logfile_dir, logfile_filename)
    elif os.path.isdir(args.logfile_path):
        # if supplied logfile filepath is actually a directory, construct filename based on input filename
        logfile_dir = args.logfile_path
        logfile_filename = inputfile_filename + '.log'
        logfile_path = os.path.join(logfile_dir, logfile_filename)
    else:
        # if supplied logfile filepath is not an existing directory, then use it for output logfile
        logfile_path = args.logfile_path
    
    # Try to deal with relative paths
   
    if not os.path.dirname(outfile_path):
        outfile_path = './' + outfile_path
    
    if not os.path.dirname(logfile_path):
        logfile_path = './' + logfile_path

    # Check that directories for output files exist, otherwise create them.
    
    if not os.path.exists(os.path.dirname(outfile_path)):
        os.makedirs(os.path.dirname(outfile_path))
    
    if not os.path.exists(os.path.dirname(logfile_path)):
        os.makedirs(os.path.dirname(logfile_path))

    # Checks that the output file and log file don't already exist. If they do, either delete them or modify name to be used.
    
    if args.force_overwrite:
        
        if os.path.isfile(outfile_path):
            os.unlink(outfile_path)
        
        if os.path.isfile(logfile_path):
            os.unlink(logfile_path)
    
    else:
        
        if os.path.isfile(outfile_path):
            outfile_path = os.path.splitext(outfile_path)[0] + '_' + timestamp + os.path.splitext(outfile_path)[1]
        
        if os.path.isfile(logfile_path):
            logfile_path = os.path.splitext(logfile_path)[0] + '_' + timestamp + os.path.splitext(logfile_path)[1]
            
    # Checks that the output file and log file aren't the same file. If they are, append different extensions.
    
    if outfile_path == outfile_path:
        print "yes"
        outfile_path = outfile_path + '.txt'
        logfile_path = logfile_path + '.log'

    # Initiate log file.
    
    log_file = open(logfile_path, 'w')
    log_file.write('Log file for entrez_qiime.py run ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.write('Using entrez_qiime v.2.0' + '\n')
    log_file.close()
    
    return(outfile_path, logfile_path)

#################################################

def get_accessions_from_input_file():
    
    # Create list of accession numbers in the input file (either FASTA or accession number list).
    
    included_accessions = {}
    
    if args.infile_fasta_path != None:
        
        input_file = open(args.infile_fasta_path, 'r')
        
        for curr_line in input_file:
            if curr_line[0] == '>':
                curr_line_accession = curr_line.split(' ',1)[0][1:]
                curr_line_accession = curr_line_accession.rstrip()
                included_accessions[curr_line_accession] = True
        
    else:
        
        input_file = open(args.infile_list_path, 'r')
        
        for curr_line in input_file:
            curr_line_accession = curr_line.rstrip()
            if curr_line_accession:
                included_accessions[curr_line_accession] = True
    
    input_file.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('get_accessions_from_input_file() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return included_accessions

#################################################

def get_merged_nodes():
    
    merged = {}
    
    mergeddmp = open(args.infile_mergeddmp_path,'r')
    
    for curr_line in mergeddmp:
        curr_line_old_taxid = curr_line.split('|')[0].strip()
        curr_line_new_taxid = curr_line.split('|')[1].strip()
        merged[curr_line_old_taxid] = curr_line_new_taxid
    
    mergeddmp.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('get_merged_nodes() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return(merged)

#################################################

def get_deleted_nodes():
    
    deleted = {}
    
    delnodesdmp = open(args.infile_delnodesdmp_path,'r')
    
    for curr_line in delnodesdmp:
        curr_line_old_taxid = curr_line.split('|')[0].strip()
        deleted[curr_line_old_taxid] = True
    
    delnodesdmp.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('get_deleted_nodes() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return(deleted)

#################################################

def obtain_nodes_for_each_accession(accessions_in_input_file, infile_acc2taxid_path, \
    ncbi_full_taxonomy, merged_taxids, deleted_taxids):
    
    included_nodes = []
    node_numbers = {}
    taxid = {}
    failed_taxids = {}
    missing_accessions = {}
    
    with open(infile_acc2taxid_path, 'r') as acc2taxid:
    
        discard_header_line = next(acc2taxid)
        
        for curr_line in acc2taxid:
        
            curr_accession, curr_taxid = curr_line.rstrip().split('\t')[1:3]
            
            if curr_accession in accessions_in_input_file:
                
                if not curr_taxid in node_numbers:
                    try:
                        included_nodes.append(ncbi_full_taxonomy[curr_taxid])
                    except KeyError:
                        if curr_taxid in merged_taxids:
                            old_taxid = curr_taxid
                            curr_taxid = merged_taxids[curr_taxid]
                            log_file = open(args.logfile_path, 'a')
                            log_file.write('The following TaxonID was changed:\n    acc2taxid gives accession=' \
                                           + str(curr_accession) + ' with taxid=' + str(old_taxid) + \
                                           ' --> changed to taxid=' + str(curr_taxid) + '.\n')
                            log_file.close()
                        elif curr_taxid in deleted_taxids:
                            old_taxid = curr_taxid
                            curr_taxid = 1 # assigns to root
                            log_file = open(args.logfile_path, 'a')
                            log_file.write('The following accession number was assigned to a deleted TaxonID:\n    acc2taxid gives accession=' \
                                           + str(curr_accession) + ' with taxid=' + str(old_taxid) + \
                                           ' --> changed to taxid=' + str(curr_taxid) + '.\n')
                            log_file.close()
                        else:
                            old_taxid = curr_taxid
                            curr_taxid = 1 # assigns to root
                            log_file = open(args.logfile_path, 'a')
                            log_file.write('The following accession number was assigned to a TaxonID that could not be found:\n    acc2taxid gives accession=' \
                                           + str(curr_accession) + '  with taxid=' + str(old_taxid) + \
                                           ' --> changed to taxid=' + str(curr_taxid) + '.\n')
                            log_file.close()
                    finally:
                        node_numbers[curr_taxid] = True
                        try:
                            included_nodes.append(ncbi_full_taxonomy[curr_taxid])
                        except KeyError:
                            failed_taxids[curr_taxid] = True
                            log_file = open(args.logfile_path, 'a')
                            log_file.write('The following taxid could not be added to included_nodes: taxid=' + str(curr_taxid) + '.\n')
                            log_file.close()
                taxid[curr_accession] = curr_taxid
    
    # determine any accession numbers from input FASTA or list file that do not have taxid information
    for curr_accession in accessions_in_input_file:
        if curr_accession not in taxid:
            missing_accessions[curr_accession] = True
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('obtain_nodes_for_each_accession() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return included_nodes, taxid, missing_accessions

#################################################

def generate_taxonid_taxonomy(included_nodes, ncbi_full_taxonomy, output_ranks):
    
    # Generate dictionary that matches taxonIDs to lineage information for all nodes represented by a sequence in the FASTA file.
    
    taxid_taxonomy = {}
    
    ranks_lookup = dict([(r,idx) for idx,r in enumerate(output_ranks)])
    
    for node in included_nodes:
        lineage = ['NA'] * len(ranks_lookup)
        curr = node
        lineage_complete = False
        while lineage_complete is False:
            if curr.Rank in ranks_lookup:
                lineage[ranks_lookup[curr.Rank]] = curr.Name
            curr = curr.Parent
            if curr is None:
                lineage_complete = True
            elif curr.TaxonId in taxid_taxonomy:
                for level in range(0,len(lineage)):
                    if (taxid_taxonomy[curr.TaxonId][level] is not 'NA') and (lineage[level] is 'NA'):
                        lineage[level] = taxid_taxonomy[curr.TaxonId][level]
                lineage_complete = True
        taxid_taxonomy[node.TaxonId] = lineage
    
    missing_taxonomy = ['NA'] * len(ranks_lookup)   # taxonomy string to be used in the output file for any accessions in input FASTA or list file that do not have corresponding TaxonID information
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('generate_taxonid_taxonomy() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return taxid_taxonomy, missing_taxonomy

#################################################

def generate_output_files(outfile_path, taxid_taxonomy, taxid, missing_accessions, missing_taxonomy):
    
    # Generate accession-to-taxonomy mapping file (e.g. for QIIME).
    
    accession_taxonomy = open(outfile_path, 'w')
    
    for curr_accession in taxid:
        accession_taxonomy.write(curr_accession + '\t' + ';'.join(taxid_taxonomy[int(taxid[curr_accession])]) + '\n')
    
    # fill in taxonomy file with NAs for accession numbers that are missing taxIDs
    log_file = open(args.logfile_path, 'a')
    for curr_accession in missing_accessions:
        accession_taxonomy.write(curr_accession + '\t' + ';'.join(missing_taxonomy) + '\n')
        log_file.write('The following accession number has no taxid associated with it: ' + curr_accession + '\n')
    log_file.close()
    
    accession_taxonomy.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('generate_output_files() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()

#################################################

if __name__ == "__main__":
    main()

#################################################
