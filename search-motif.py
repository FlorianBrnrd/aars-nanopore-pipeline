import parasail
import argparse
from Bio import SeqIO


DEFAULT_ALIGN_PARAMS = {'match': 1,
                        'mismatch': -1,
                        'gap_open': 2,    # penalty to create a gap
                        'gap_extend': 1}  # penalty to extend a gap 


def semi_global_alignment(reference, query, params=DEFAULT_ALIGN_PARAMS):
    """ 
    Perform semi-global alignment between two sequences
    """

    subs_mat = parasail.matrix_create("ACGT", params['match'], params['mismatch'])
    alignment = parasail.sg_trace_striped_32(reference, query, params['gap_open'], params['gap_extend'], subs_mat)

    return alignment


def parse_args():
    """
    Parse user's arguments
    """
    
    parser = argparse.ArgumentParser(description='Find motif in sequences',
                                     add_help=True)

    main_group = parser.add_argument_group('Main options')

    main_group.add_argument('-i', '--input', required=True,
                            help='Fasta file')

    main_group.add_argument('-o', '--output', required=True,
                            help= 'Output file')

    main_group.add_argument('-m', '--motif', default='AACTAACGCTATATAAGTATCAGTTTCTGTACTTTATTG',
                            help='Sequence to be searched (default: SL sequence)')

    main_group.add_argument('-s', '--sensitivity', default=0.7,
                            help='Sensitivity value - must be between 0 and 1 (default: 0.7)')

    args = parser.parse_args()

    return args


def evaluate_motif(record, motif, sensitivity):
    
    """
    Perform semi-global alignment between reads and the query motif.
    Match is accepted given a sensitivity parameter reflecting the level of similitude between the match and the expected motif.
    This parameter ranges between 0 and 1 (default: 0.7).
    If no match is found, the query sequence is slightly shortened and searched again, until we reach a minimal sequence.
    """

    aln = semi_global_alignment(str(record.seq), motif)
    score = aln.score

    # No match with full sequence -> enter degenerated sequence search
    if score < len(motif) * sensitivity:
        
        #  Each loop, removes the two first bases of the previous query sequence
        for n in reversed(range(10, len(motif), 2)):

            aln = semi_global_alignment(str(record.seq), motif[-n:])
            score = aln.score

            # Match -> motif found. Loop terminated.
            if score >= sensitivity * n:
                match = 1
                break

            # No match + minimal query sequence -> No motif found. Loop terminated.
            elif n == 10 and score < sensitivity * n:
                match = 0

            # No match -> try shorter query sequence. Loop continues.
            else:
                continue

    # Match with full sequence.
    else:
        match = 1

    return match


def search_motif():
    """
    Process input .fasta file to extract reads containing a given motif
    """

    args = parse_args()

    input = args.input
    output = args.output
    motif = args.motif
    sensitivity = args.sensitivity

    with open(output, "w+") as outfile:
        
        # Loop over records in fasta file
        for record in SeqIO.parse(input, "fasta"):
            
            # Search query in each sequence
            if evaluate_motif(record, motif=motif, sensitivity=sensitivity) == 1:
                
                # Output record if a match is found
                outfile.write(f'>{record.id}\n{record.seq}\n')


if __name__ == '__main__':

    search_motif()
