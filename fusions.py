################################################################################################
## Filtering of set fusions in BAM file                                                        #
## Developed and maintained by: Jonathan Tang                                                  #
## Organization: Washington University in St. Louis                                            #
## Arguments: 1. BAM, CRAM, or SAM input, 2 List of fusions of interest,                       #
##      3. Reference BED file, 4. Output file path                                             #
## Example: fusions.py regions.bam ROIS.txt gene_ref.bed fusion_py.txt                         #
################################################################################################

# Libraries
import sys
import numpy as np
import pandas as pd
import pysam

# Functions
#   ROI Lookup
#       Looks up the chromosome, start point, and end point of a gene of interest.
#           Checks if any row of reference gene column contains a gene in the roi list,
#           then pulls output columns. Throws an error if an empty dataframe is returned.
#       Inputs:
#           roi: list of genes
#           ref: reference dataframe from reference file
#       Outputs:
#           goi: dataframe with columns:
#               gene: name of gene
#               chromosome: chromosome of gene
#               startcoord: starting basepair of gene
#               endcoord: ending basepair of gene
def roi_lookup(roi, ref):
    goi = ref[ref.gene.isin(roi)][['gene', 'chromosome', 'startcoord', 'endcoord', 'strand']]
    if goi.empty:
        raise ValueError("Genes not identified in reference file")
    else:
        return goi

#   Load BAM
#       Assigns pysam AlignmentFile, creates columns, and iterates through pysam AlignmentFile
#           while writing SA and MQ tags to a column specifically. Includes all tags in tags column
#           in case they are needed in the future.
#       Inputs:
#           bam_path: path to BAM file
#       Outputs:
#           dataframe with columns:
#               qname: query name
#               flag: read flags
#               rname: chromosome of read
#               pos: location of read
#               mapq: mapping quality of read
#               cigar: cigar of read
#               rnext: chromosome of mate
#               pnext: location of mate
#               tlen: length of transcript
#               seq: sequence of read
#               qual: quality score
#               sa: supplementary alignments
#               mateq: mate MAPQ
#               tags: all tags
def load_bam(bam_path):
    bam = pysam.AlignmentFile(bam_path)
    qname = []
    flag = []
    rname = []
    pos = []
    mapq = []
    cigar = []
    rnext = []
    pnext = []
    tlen = []
    seq = []
    qual = []
    sa = []
    mateq = []
    all_tags = []
    for read in bam:
        qname.append(read.query_name)
        flag.append(read.flag)
        rname.append(read.reference_id)
        pos.append(read.pos)
        mapq.append(read.mapping_quality)
        cigar.append(read.cigarstring)
        rnext.append(read.next_reference_id)
        pnext.append(read.next_reference_start)
        tlen.append(read.template_length)
        seq.append(read.query_sequence)
        qual.append(read.query_name)
        if (read.has_tag('SA')):
            sa.append(read.get_tag('SA'))
        else:
            sa.append('NA')
        if (read.has_tag('MQ')):
            mateq.append(read.get_tag('MQ'))
        else:
            mateq.append('NA')
        all_tags.append(read.get_tags())
    return(pd.DataFrame({
        'qname': qname,
        'flag': flag,
        'rname': rname,
        'pos': pos,
        'mapq': mapq,
        'cigar': cigar,
        'rnext': rnext,
        'pnext': pnext,
        'tlen': tlen,
        'seq': seq,
        'qual': qual,
        'sa': sa,
        'mateq': mateq,
        'tags': all_tags
    }))

# Collects file paths as passed by arguments
if len(sys.argv) == 5:
    bam_path = sys.argv[1]
    roi_path = sys.argv[2]
    ref_path = sys.argv[3]
    out_path = sys.argv[4]
else:
    raise Exception("Incorrect number of arguments: need 4")

# Load gene reference dataframe and name columns
ref_cols = [
    'chromosome',
    'startcoord',
    'endcoord',
    'ident',
    'gene',
    'strand'
]

ref = pd.read_table(ref_path, delimiter = '\t', lineterminator = '\n', header = None)
ref.set_axis(ref_cols, axis = 1, inplace = True)

# Load ROIs
#   Collect genes of interest and their locations per reference file by looking up
#   locations of gene, merging them into the ROIs frame, and removing original name columns
#   Rename expected strands to 0/1 from +/- respectively.
rois = pd.read_table(roi_path, delimiter = '_', lineterminator = '\n')
roi_locs_1 = roi_lookup(rois.gene1, ref)
roi_locs_2 = roi_lookup(rois.gene2, ref)
rois = rois.merge(roi_locs_1, left_on = 'gene1', right_on = 'gene')
rois = rois.merge(roi_locs_2, left_on = 'gene2', right_on = 'gene')
rois = rois.loc[:, ~rois.columns.isin(['gene1', 'gene2'])]
rois.strand_x = rois.strand_x.str.replace('+', '0')
rois.strand_x = rois.strand_x.str.replace('-', '1')
rois.strand_y = rois.strand_y.str.replace('+', '0')
rois.strand_y = rois.strand_y.str.replace('-', '1')

# Load BAM
bam = load_bam(bam_path)

# Filter readthroughs
#   Removes readthroughs within 1 Mb by filtering out results with an rnext of '=' (same
#   chromosome) and a distance between 'pos' and 'pnext' of <=1000000
bam_filter = (bam.rnext == '=') & (abs(bam.pos - bam.pnext <= 1000000))
bam.drop(bam[bam_filter].index, inplace = True)

# Rename '=' values in 'rnext' to corresponding chromosome in 'rname'
bam_rename = bam[bam.rnext == '=']
bam.loc[bam_rename.index, 'rnext'] = bam.rname[bam_rename.index]

# Strandedness
#   Formats flag to binary and pulls strand values (str[7] is mate reverse strand, str[8]
#   is read reverse strand). 0 is + and 1 is -.
binflags = pd.DataFrame([format(i, '013b') for i in bam.flag])
binflags.set_index(bam.index, inplace = True)
bam['strand1'] = binflags[0].str[7]
bam['strand2'] = binflags[0].str[8]

# Supplemental alignments
#   Splits sa column to list, returns to the bam, and explodes the list so that each
#   supplementary alignment has its own row and each main alignment has a row with ''
supp_align_split = bam.sa.str.split(';')
bam.loc[supp_align_split.index, 'sa'] = supp_align_split
bam = bam.explode('sa').fillna('')
bam.reset_index(drop = True, inplace = True)

#   Pulls SA rows, splits them for information, and creates a dataframe with them. Replaces
#   +/- strands with 0/1 respectively. Sets flag to 4096 and sets the correct index.
supp_align = bam.loc[bam.sa.str.startswith('chr'), 'sa'].str.split(',')
supp_align_df = pd.DataFrame(supp_align.tolist())
supp_align_df.set_axis(['rnext', 'pnext', 'strand2', 'cigar', 'mapq', 'NM'], axis = 1, inplace = True)
supp_align_df.strand2 = supp_align_df.strand2.str.replace('+', '0')
supp_align_df.strand2 = supp_align_df.strand2.str.replace('-', '1')
supp_align_df['flag'] = 4096
supp_align_df.set_index(supp_align.index, inplace = True)

#   Updates bam with SA dataframe and fixes some types.
bam.update(supp_align_df, overwrite = True)
bam.pos = pd.to_numeric(bam.pos)
bam.pnext = pd.to_numeric(bam.pnext)
bam.flag = pd.to_numeric(bam.flag)
bam.mapq = pd.to_numeric(bam.mapq)

# Pull relevant results
#   Iterates through rows of rois and pulls rows with appropriate chromosome locations,
#   adds columns with expected strand orientation, then concats them to dataframe out.
#   Note: expected strand orientation is for the fusion gene pairing, and is not accurate
#   to the gene itself, e.g. if GENE1-GENE2 is +/- but is read as GENE2-GENE1, it will still
#   be labeled as +/-. This shouldn't affect pair orientation logic though. Also note that MAPQ
#   filtering is performed to exclude supplementary alignments with MAPQs below threshold.
out = pd.DataFrame()
for fusion in rois.iterrows():
    bam_filter = \
    ((bam.rname == fusion[1].chromosome_x) & \
        (bam.rnext == fusion[1].chromosome_y) & \
        (bam.pos > fusion[1].startcoord_x) & \
        (bam.pos < fusion[1].endcoord_x) & \
        (bam.pnext > fusion[1].startcoord_y) & \
        (bam.pnext < fusion[1].endcoord_y) & \
        (bam.mapq >= 30) & \
        (bam.mateq >= 30)) | \
    ((bam.rname == fusion[1].chromosome_y) & \
        (bam.rnext == fusion[1].chromosome_x) & \
        (bam.pos > fusion[1].startcoord_y) & \
        (bam.pos < fusion[1].endcoord_y) & \
        (bam.pnext > fusion[1].startcoord_x) & \
        (bam.pnext < fusion[1].endcoord_x) & \
        (bam.mapq >= 30) & \
        (bam.mateq >= 30))
    
    bam_fusion = bam[bam_filter]
    bam_fusion.reset_index(drop = True, inplace = True)
    bam_fusion = bam_fusion.assign(strand1_exp = fusion[1].strand_x)
    bam_fusion = bam_fusion.assign(strand2_exp = fusion[1].strand_y)
    bam_fusion = bam_fusion.assign(gene1 = fusion[1].gene_x)
    bam_fusion = bam_fusion.assign(gene2 = fusion[1].gene_y)
    out = pd.concat([out, bam_fusion])

out = out.fillna(value = 'NA')
out['strand1'] = pd.to_numeric(out.strand1)
out['strand2'] = pd.to_numeric(out.strand2)
out['strand1_exp'] = pd.to_numeric(out.strand1_exp)
out['strand2_exp'] = pd.to_numeric(out.strand2_exp)

# Pair orientation filter
#   Checks orientations against expected strand directions. Supplementary alignments (flag >= 2048)
#   are checked for same or complementary directions and pairs (flag < 2048) are checked for
#   discordant directions. None or '' values are replaced with 'NA' for clean awk later.
pairs = out[(out.flag < 2048)]
pairs_actual = (pairs.strand1 + pairs.strand2) % 2
pairs_exp = (pairs.strand1_exp + pairs.strand2_exp) % 2
pairs_keep = (pairs_actual != pairs_exp)

supps = out[(out.flag >= 2048)]
supps_actual = (supps.strand1 + supps.strand2) % 2
supps_exp = (supps.strand1_exp + supps.strand2_exp) % 2
supps_keep = (supps_actual == supps_exp)

out_filtered = pd.concat([pairs[pairs_keep], supps[supps_keep]])
out_filtered = out_filtered.fillna(value = 'NA').replace('', 'NA').sort_values('qname').reset_index(drop = True)

# Write results to output
out_filtered.to_csv(out_path, sep = '\t', line_terminator = '\n', index = False)
