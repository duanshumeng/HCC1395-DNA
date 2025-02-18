
import argparse, re, os
import genomic_file_handlers as genome

def run():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Variant Call Type, i.e., snp or indel
    parser.add_argument('-infile',  '--input-vcf',  type=str, help='Input VCF file', required=True)
    parser.add_argument('-outfile', '--output-vcf', type=str, help='Output VCF file', required=True)

    # Parse the arguments:
    args = parser.parse_args()
    infile = args.input_vcf
    outfile = args.output_vcf

    tumor_name = os.path.basename(infile).split('.')[0]+'.T'
    normal_name = os.path.basename(infile).split('.')[0]+'.N'

    return infile, outfile,tumor_name,normal_name


def convert(infile, outfile,tumor_name,normal_name):

    idx_chrom,idx_pos,idx_id,idx_ref,idx_alt,idx_qual,idx_filter,idx_info,idx_format,idx_SM1,idx_SM2 = 0,1,2,3,4,5,6,7,8,9,10

    with genome.open_textfile(infile) as vcf, open(outfile, 'w') as vcfout:

        line_i = vcf.readline().rstrip()

        # VCF header
        while line_i.startswith('##'):
            vcfout.write( line_i + '\n' )
            line_i = vcf.readline().rstrip()

        # This is the #CHROM line:
        header = line_i.split('\t')
        header[header.index('TUMOR')] = tumor_name
        header[header.index('NORMAL')] = normal_name
        print(header)
        vcfout.write('\t'.join(header) + '\n' )
        line_i = vcf.readline().rstrip()


        while line_i:

            # Print "SomaticSniper" into the INFO field if it is called so, otherwise never mind.
            item = line_i.split('\t')

            # In the REF field, non-GCTA characters should be changed to N to fit the VCF standard:
            item[idx_ref] = re.sub( r'[^GCTA]', 'N', item[idx_ref], flags=re.I )
            line_i = '\t'.join(item)

            vcfout.write( line_i + '\n' )

            line_i = vcf.readline().rstrip()


if __name__ == '__main__':
    infile, outfile,tumor_name,normal_name = run()
    convert(infile, outfile,tumor_name,normal_name)