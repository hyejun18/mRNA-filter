import argparse

##################################################
# Args processing
##################################################

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--refflat', '-r', 
        type=str, 
        required=True, 
        dest='refFlat', 
        help='/path/to/RefFlat.txt'
        )
    parser.add_argument(
        '--fasta-directory', '-a', 
        type=str, 
        required=True, 
        dest='fa_dir', 
        help='/path/to/fasta-containing-directory'
        )
    parser.add_argument(
        '--outfile', '-o', 
        type=str, 
        required=True, 
        dest='out_file', 
        help='/path/to/results.txt'
        )
    parser.add_argument(
        '--fold-change', '-f', 
        type=str, 
        dest='fold_change_table', 
        help='/path/to/mRNA-fold-change.txt'
        )
    args = parser.parse_args()
    return args
# end: def getArgs