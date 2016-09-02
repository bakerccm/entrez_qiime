#!/usr/local/bin/python

# Author: Chris Baker
# Affiliation: Pierce Lab, Department of Organismic and Evolutionary Biology, Harvard University
# Email: ccmbaker@fas.harvard.edu
# Last modified: 9 June 2012

# module load bio/qiime-1.5.0

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

parser = argparse.ArgumentParser(description='Generates QIIME-compatible files from an Entrez FASTA file and the NCBI taxonomy database.')

# Path including filename for input FASTA file.
# Expects deflines that look like ">gi|######|..." where "######" is
# the GI number we want, and "..." is arbitrary text.
parser.add_argument('-i','--inputfasta', 
    metavar='path/file.ext', 
    dest='infile_fasta_path', 
    action='store', 
    required=True, 
    help='The path, including filename, of the input FASTA file. Filename should have an extension, as it will be used as the basis for output file names e.g. file_stripped.ext. [REQUIRED]')

# Path to output directory. Output fasta file, output qiime taxonomy mapping file,
# output gi-taxid file and output list of gis to be included in database will be saved here,
# plus file with taxonomy data for each TaxonId if requested with -t option.
parser.add_argument('-o','--outputdir', 
    metavar='path/', 
    dest='outdir_path', 
    action='store', 
    default='./', 
    help='The directory where output files should be saved. Will be created if it does not exist. [Default: ./]')

# Force overwrite if any output files already exist?
parser.add_argument('-f','--force', 
    dest='force_overwrite', 
    action='store_true', 
    default=False, 
    help='If -f or --force option are passed, output files will overwrite any existing files with the same names. Without this option, output filenames will be modified by the addition of a date/time string in the event that a file with the same name already exists.  [Default: False]')

# Filename and path for NCBI node-parent list (from taxdmp.zip)
parser.add_argument('-n','--nodes', 
    metavar='path/', 
    dest='ncbi_taxonomy_dir', 
    action='store', 
    default='./', 
    help='The directory where the files nodes.dmp, names.dmp, merged.dmp and delnodes.dmp are located. The files are typically downloaded from the NCBI in the compressed archive taxdump.tar.gz. [Default: ./]')

# Filename and path for NCBI gi-taxid list (from nt database dump)
parser.add_argument('-g','--gitaxid', 
    metavar='path/file.ext', 
    dest='infile_gi_taxid_dmp_path', 
    action='store', 
    default='./gi_taxid_nucl.dmp', 
    help='The path, including filename, of the (uncompressed) input gi-taxid file e.g. from the nt database dump, typically downloaded from the NCBI as gi_taxid_nucl.dmp.gz and uncompressed to gi_taxid_nucl.dmp.  [Default: ./gi_taxid_nucl.dmp]')

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

# Should script generate file with taxonomy data for each TaxonId?
parser.add_argument('-t','--taxid', 
    dest='write_taxonid_taxonomy_tofile', 
    action='store_true', 
    default=False, 
    help='If -t or --taxid are passed, an additional output file will be generated that lists taxonomy for each TaxonID.  [Default: False]')

args = parser.parse_args()

#################################################

def main():
    
    args.infile_nodesdmp_path, args.infile_namesdmp_path, args.infile_mergeddmp_path, \
        args.infile_delnodesdmp_path = input_files()
    
    args.outfile_fasta_path, args.outfile_taxonid_taxonomy_path, args.outfile_gi_taxonomy_path, \
        args.outfile_gi_taxid_path, args.outfile_gilist_path, args.logfile_path = output_files()
    
    args.output_ranks = args.output_ranks.replace('_',',').split(',')
    
    gi_numbers_in_input_fasta = get_gi_numbers_from_fasta()
    
    ncbi_full_taxonomy = NcbiTaxonomyFromFiles(open(args.infile_nodesdmp_path), open(args.infile_namesdmp_path))
    
    merged_taxids = get_merged_nodes()
    
    deleted_taxids = get_deleted_nodes()
    
    included_nodes, taxid = obtain_nodes_for_GIs(gi_numbers_in_input_fasta, \
        args.infile_gi_taxid_dmp_path, args.outfile_gi_taxid_path, ncbi_full_taxonomy, \
        merged_taxids, deleted_taxids)
    
    generate_stripped_fasta(taxid)
    
    taxid_taxonomy = generate_taxonid_taxonomy(included_nodes, ncbi_full_taxonomy, args.output_ranks)
    
    if args.write_taxonid_taxonomy_tofile:
        write_taxonid_taxonomy(args.outfile_taxonid_taxonomy_path, taxid_taxonomy)
    
    generate_gi_files(args.outfile_gi_taxonomy_path, args.outfile_gi_taxid_path, \
        args.outfile_gilist_path, taxid_taxonomy, taxid)
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('Finished running entrez_qiime.py at ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()

#################################################

def input_files():
    
    if not args.ncbi_taxonomy_dir[-1:] == '/':
        args.ncbi_taxonomy_dir = args.ncbi_taxonomy_dir + '/'
    
    infile_nodesdmp_path = args.ncbi_taxonomy_dir + 'nodes.dmp'
    
    if not os.path.isfile(infile_nodesdmp_path):
        exit('Error: cannot find file nodes.dmp at ' + infile_nodesdmp_path)
    
    infile_namesdmp_path = args.ncbi_taxonomy_dir + 'names.dmp'
    
    if not os.path.isfile(infile_namesdmp_path):
        exit('Error: cannot find file names.dmp at ' + infile_namesdmp_path)
    
    infile_mergeddmp_path = args.ncbi_taxonomy_dir + 'merged.dmp'
    
    if not os.path.isfile(infile_mergeddmp_path):
        exit('Error: cannot find file merged.dmp at ' + infile_mergeddmp_path)
    
    infile_delnodesdmp_path = args.ncbi_taxonomy_dir + 'delnodes.dmp'
    
    if not os.path.isfile(infile_delnodesdmp_path):
        exit('Error: cannot find file delnodes.dmp at ' + infile_delnodesdmp_path)
    
    if not os.path.isfile(args.infile_fasta_path):
        exit('Error: cannot find input FASTA file at ' + args.infile_fasta_path)
    
    if not os.path.isfile(args.infile_gi_taxid_dmp_path):
        exit('Error: cannot find input GI - TaxonID file at ' + args.infile_gi_taxid_dmp_path)
    
    return(infile_nodesdmp_path, 
           infile_namesdmp_path, 
           infile_mergeddmp_path, 
           infile_delnodesdmp_path)

#################################################

def output_files():
    
    timestamp = strftime("%Y%m%d%H%M%S",localtime())
    
    ''' Checks that outdir_path exists, otherwise creates it. '''
    
    if not args.outdir_path[-1:] == '/':
        args.outdir_path = args.outdir_path + '/'
    
    if os.path.isfile(args.outdir_path):
        os.makedirs(args.outdir_path + timestamp)
    
    if not os.path.exists(args.outdir_path):
        os.makedirs(args.outdir_path)
    
    ''' Constructs output file names from command line options. '''
    
    inputfilenameroot = os.path.splitext(os.path.basename(args.infile_fasta_path))[0]
    
    inputfilenameextn = os.path.splitext(os.path.basename(args.infile_fasta_path))[1]
    
    outfile_fasta_path = args.outdir_path + inputfilenameroot + '_stripped' + inputfilenameextn # keeps the same extension as the input fasta file rather than forcing it to be .fasta or something like that
    
    outfile_taxonid_taxonomy_path = args.outdir_path + inputfilenameroot + '_taxid_taxonomy.txt'
    
    outfile_gi_taxonomy_path = args.outdir_path + inputfilenameroot + '_gi_taxonomy.txt'
    
    outfile_gi_taxid_path = args.outdir_path + inputfilenameroot + '_gi_taxid.txt'
    
    outfile_gilist_path = args.outdir_path + inputfilenameroot + '_gi.txt'
    
    logfile_path = args.outdir_path + inputfilenameroot + '.log'
    
    ''' Checks that the five output filenames constructed above don't already refer to files. '''
    ''' If they do, either deletes them or modifies name to be used.                          '''
    
    if args.force_overwrite:
        
        if os.path.isfile(outfile_fasta_path):
            os.unlink(outfile_fasta_path)
        
        if os.path.isfile(outfile_taxonid_taxonomy_path) and args.write_taxonid_taxonomy_tofile:
            os.unlink(outfile_taxonid_taxonomy_path)
        
        if os.path.isfile(outfile_gi_taxonomy_path):
            os.unlink(outfile_gi_taxonomy_path)
        
        if os.path.isfile(outfile_gi_taxid_path):
            os.unlink(outfile_gi_taxid_path)
        
        if os.path.isfile(outfile_gilist_path):
            os.unlink(outfile_gilist_path)
            
        if os.path.isfile(logfile_path):
            os.unlink(logfile_path)
    
    else:
        
        if os.path.isfile(outfile_fasta_path):
            outfile_fasta_path = os.path.splitext(outfile_fasta_path)[0] + '_' + timestamp + os.path.splitext(outfile_fasta_path)[1]
        
        if os.path.isfile(outfile_taxonid_taxonomy_path)  and args.write_taxonid_taxonomy_tofile:
            outfile_taxonid_taxonomy_path = os.path.splitext(outfile_taxonid_taxonomy_path)[0] + '_' + timestamp + os.path.splitext(outfile_taxonid_taxonomy_path)[1]
        
        if os.path.isfile(outfile_gi_taxonomy_path):
            outfile_gi_taxonomy_path = os.path.splitext(outfile_gi_taxonomy_path)[0] + '_' + timestamp + os.path.splitext(outfile_gi_taxonomy_path)[1]
        
        if os.path.isfile(outfile_gi_taxid_path):
            outfile_gi_taxid_path = os.path.splitext(outfile_gi_taxid_path)[0] + '_' + timestamp + os.path.splitext(outfile_gi_taxid_path)[1]
        
        if os.path.isfile(outfile_gilist_path):
            outfile_gilist_path = os.path.splitext(outfile_gilist_path)[0] + '_' + timestamp + os.path.splitext(outfile_gilist_path)[1]
        
        if os.path.isfile(logfile_path):
            logfile_path = os.path.splitext(logfile_path)[0] + '_' + timestamp + os.path.splitext(logfile_path)[1]
            
    ''' Initiate log file. '''
    
    log_file = open(logfile_path, 'w')
    log_file.write('Log file for entrez_qiime.py run ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return(outfile_fasta_path, 
           outfile_taxonid_taxonomy_path, 
           outfile_gi_taxonomy_path, 
           outfile_gi_taxid_path, 
           outfile_gilist_path, 
           logfile_path)

#################################################

def get_gi_numbers_from_fasta():
    
    """creates list of gi numbers in the input fasta file"""
    
    in_fasta_file = open(args.infile_fasta_path, 'r')
    
    included_GIs = {}
    
    for curr_line in in_fasta_file:
        if curr_line[0] == '>':
            curr_line_gi_number = int(curr_line.split('|',2)[1])
            included_GIs[curr_line_gi_number] = True
    
    in_fasta_file.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('get_gi_numbers_from_fasta() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return included_GIs

#################################################

def get_merged_nodes():
    
    merged = {}
    
    mergeddmp = open(args.infile_mergeddmp_path,'r')
    
    for curr_line in mergeddmp:
        curr_line_old_taxid = int(curr_line.split('|')[0].strip())
        curr_line_new_taxid = int(curr_line.split('|')[1].strip())
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
        curr_line_old_taxid = int(curr_line.split('|')[0].strip())
        deleted[curr_line_old_taxid] = True
    
    delnodesdmp.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('get_deleted_nodes() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return(deleted)

#################################################

def obtain_nodes_for_GIs(gi_numbers_in_input_fasta, infile_gi_taxid_dmp_path, outfile_gi_taxid_path, \
    ncbi_full_taxonomy, merged_taxids, deleted_taxids):
    
    included_nodes = []
    node_numbers = {}
    taxid = {}
    failed_taxids = {}
    
    gi_taxid_dmp = open(infile_gi_taxid_dmp_path, 'r')
    
    for curr_gi_taxid in gi_taxid_dmp:
        curr_gi, curr_taxid = curr_gi_taxid.rstrip().partition('\t')[::2]
        curr_gi = int(curr_gi)
        curr_taxid = int(curr_taxid)
        if curr_gi in gi_numbers_in_input_fasta:
            if not curr_taxid in node_numbers:
                try:
                    included_nodes.append(ncbi_full_taxonomy[curr_taxid])
                except KeyError:
                    if curr_taxid in merged_taxids:
                        old_taxid = curr_taxid
                        curr_taxid = merged_taxids[curr_taxid]
                        log_file = open(args.logfile_path, 'a')
                        log_file.write('The following TaxonID was changed:\n    gi_taxid_dmp gives gi=' \
                                       + str(curr_gi) + ' with taxid=' + str(old_taxid) + \
                                       ' --> changed to taxid=' + str(curr_taxid) + '.\n')
                        log_file.close()
                    elif curr_taxid in deleted_taxids:
                        old_taxid = curr_taxid
                        curr_taxid = 1 # assigns to root
                        log_file = open(args.logfile_path, 'a')
                        log_file.write('The following GI was assigned to a deleted TaxonID:\n    gi_taxid_dmp gives gi=' \
                                       + str(curr_gi) + ' with taxid=' + str(old_taxid) + \
                                       ' --> changed to taxid=' + str(curr_taxid) + '.\n')
                        log_file.close()
                    else:
                        old_taxid = curr_taxid
                        curr_taxid = 1 # assigns to root
                        log_file = open(args.logfile_path, 'a')
                        log_file.write('The following GI was assigned to a TaxonID that could not be found:\n    gi_taxid_dmp gives gi=' \
                                       + str(curr_gi) + '  with taxid=' + str(old_taxid) + \
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
            taxid[curr_gi] = curr_taxid
    
    gi_taxid_dmp.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('obtain_nodes_for_GIs() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return included_nodes, taxid

#################################################

def generate_stripped_fasta(taxid):

    """ creates stripped fasta file including only GIs that have taxon information """
    
    in_fasta_file = open(args.infile_fasta_path, 'r')
    out_fasta_file = open(args.outfile_fasta_path, 'w')
    
    for curr_line in in_fasta_file:
        if curr_line[0] == '>':
            curr_sequence_gi_number = int(curr_line.split('|',2)[1])
        if curr_sequence_gi_number in taxid:
            if curr_line[0] == '>':
                out_fasta_file.write('>' + str(curr_sequence_gi_number) + '\n')
            else:
                out_fasta_file.write(curr_line)
        elif curr_line[0] == '>':
            log_file = open(args.logfile_path, 'a')
            log_file.write('The sequence with gi=' + str(curr_sequence_gi_number) + \
                ' is in the original fasta file but is not in gi_taxid.dmp\n' + \
                '    and has been excluded from the output fasta file.\n')
            log_file.close()
    
    in_fasta_file.close()
    out_fasta_file.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('generate_stripped_fasta() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()

#################################################

def generate_taxonid_taxonomy(included_nodes, ncbi_full_taxonomy, output_ranks):
    
    """generates dictionary that matches taxonIDs to lineage information """
    """for all nodes represented by a GI in the masked fasta file        """
    
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
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('generate_taxonid_taxonomy() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()
    
    return taxid_taxonomy

#################################################

def write_taxonid_taxonomy(outfile_taxonid_taxonomy_path, taxid_taxonomy):
    
    """ writes taxonid and taxonomy info to file (i.e. one line per node)         """
    
    fout = open(outfile_taxonid_taxonomy_path, 'w')
    
    for key, value in taxid_taxonomy.items():
        fout.write(str(key) + '\t' + ';'.join(value) + '\n')
    
    fout.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('write_taxonid_taxonomy() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()

#################################################

def generate_gi_files(outfile_gi_taxonomy_path, outfile_gi_taxid_path, \
    outfile_gilist_path, taxid_taxonomy, taxid):
    
    """ generates GI-taxonomy mapping file (e.g. for QIIME), GI-taxid map file    """
    """ (e.g. for processing files for megan), and GI list (e.g. as an input to   """
    """ the blastdbaliastool)                                                     """
    
    gi_taxonomy = open(outfile_gi_taxonomy_path, 'w')
    gi_taxid = open(outfile_gi_taxid_path, 'w')
    gi_list = open(outfile_gilist_path, 'w')
    
    for curr_gi_number in taxid:
        gi_list.write(str(curr_gi_number) + '\n')
        gi_taxid.write(str(curr_gi_number) + '\t' + str(taxid[int(curr_gi_number)]) + '\n')
        gi_taxonomy.write(str(curr_gi_number) + '\t' + ';'.join(taxid_taxonomy[taxid[curr_gi_number]]) + '\n')
    
    gi_taxonomy.close()
    gi_taxid.close()
    gi_list.close()
    
    log_file = open(args.logfile_path, 'a')
    log_file.write('generate_gi_files() finished ' + strftime("%H:%M:%S on %d-%m-%Y",localtime()) + '\n')
    log_file.close()

#################################################

if __name__ == "__main__":
    main()

#################################################
