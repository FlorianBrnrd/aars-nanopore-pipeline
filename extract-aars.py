import pysam
import argparse


# Reference name of aaRS transcripts
aars = ['LtaP27.1400.mRNA', 'LtaP35.1480.mRNA', 'LtaP22.1510.mRNA', 'LtaP36.5770.mRNA', 'LtaP30.3140.mRNA',
        'LtaP13.1000.mRNA', 'LtaP11.0110.mRNA', 'LtaP12.0270.mRNA', 'LtaP14.1430.mRNA', 'LtaP15.1380.mRNA',
        'LtaP21.0830.mRNA', 'LtaP30.0700.mRNA', 'LtaP34.2190.mRNA', 'LtaP30.3280.mRNA', 'LtaP36.3950.mRNA',
        'LtaP32.0940.mRNA', 'LtaP19.0950.mRNA', 'LtaP23.0360.mRNA', 'LtaP29.0060.mRNA', 'LtaP15.0260.mRNA',
        'LtaP30.0190.mRNA', 'LtaP30.0530.mRNA', 'LtaP21.1070.mRNA', 'LtaP29.2400.mRNA', 'LtaP18.1190.mRNA',
        'LtaP18.1200.mRNA']


def parse_args():
    """
    Input arguments for processing transcriptome alignment file and extracting nanopore reads mapped to aaRS transcripts into a .fasta file
    """

    parser = argparse.ArgumentParser(description='extract nanopore reads mapped to aaRS transcripts',
                                     add_help=True)

    main_group = parser.add_argument_group('Main options')

    main_group.add_argument('-i', '--input', required=True,
                            help='Transcriptome alignments file')

    main_group.add_argument('-o', '--output', required=True,
                            help= 'Output file name')

    args = parser.parse_args()

    return args


def extract_aars():
    """
    extract sequence of reads mapped to aaRS into a .fasta file
    """

    args = parse_args()

    input = args.input
    output = args.output

    # load alignment file
    transcriptome = pysam.AlignmentFile(input, 'rb')

    with open(output, 'w+') as fasta:

        for read in transcriptome:
                
            # Select only primary alignments
            if not read.is_unmapped and not read.is_supplementary and not read.is_secondary and read.seq is not None:
                
                # Get reference name
                ref = transcriptome.get_reference_name(read.reference_id)
                
                # Extract read sequence if the reference name corresponds to an aaRS transcript
                for gene in aars:
                        
                    if ref == gene:
                        # Write sequence as fasta format
                        fasta.write('>{}\n{}\n'.format(read.query_name, read.seq))


if __name__ == '__main__':

    extract_aars()
