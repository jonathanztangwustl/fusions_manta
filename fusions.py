################################################################################################
## Filtering of set fusions in BAM file                                                        #
## Developed and maintained by: Jonathan Tang                                                  #
## Email: jonathanztang@wustl.edu                                                              #
## Phone: +1-6106390698                                                                        #
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

# For testing
#   Generally, bam_ will be raw data, roi_ will be regions of interest/targets,
#   ref_ will be the reference file with genes and locations, and out_ will be output.
bam_path = 'fusions.bam'
roi_path = '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/ROIs.txt'
ref_path = '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/gene_reference.bed'
out_path = '/storage1/fs1/timley/Active/aml_ppg/tmp/jonathanztang/breakpoint_reader/terra/output.txt'

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

# Open BAM and read data
#   Uses pysam to call samtools view on bam, then parses \n and \t to dataframe.
#   Renames columns, removes 'None' lines, fixes data types.
bam_cols = [
    'qname',
    'flag',
    'rname',
    'pos',
    'mapq',
    'cigar',
    'rnext',
    'pnext',
    'tlen',
    'seq',
    'qual',
    'tag1',
    'tag2',
    'tag3',
    'tag4',
    'tag5',
    'tag6',
    'tag7',
    'tag8',
    'tag9',
    'tag10'
]

bam = pysam.view(bam_path)
bam_list = bam.split('\n')
bam = pd.DataFrame([line.split('\t') for line in bam_list])
bam.set_axis(1, bam_cols[0:len(bam.columns)], inplace = True)
bam.dropna(subset = ['pos', 'pnext'], inplace = True)
bam.pos = pd.to_numeric(bam.pos)
bam.pnext = pd.to_numeric(bam.pnext)
bam.flag = pd.to_numeric(bam.flag)

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
#   Pulls values in 'tag' columns that start with 'SA' and then realigns them to a column.
#   Searches all tag columns, because tags aren't constant. Switches to series for update().
tag_col_names = [col for col in bam if col.startswith('tag')]
tag_cols = bam[tag_col_names].fillna('empty')
supp_align_filter = tag_cols[tag_cols.apply(lambda col: col.str.startswith('SA'), axis = 1)]
supp_align_raw = supp_align_filter.tag1
for tag in tag_col_names[1:]:
    supp_align_raw.update((supp_align_filter[[tag]]).squeeze())

# Supplemental alignments - explode
#   Removes SA tag, splits to list, returns to bam, explodes list to
#   rows so each SA has its own row and each main alignment has a row with '' in tag3.
supp_align_split = supp_align_raw.str.replace('SA:Z:', '').str.split(';')
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
#   be labeled as +/-. This shouldn't affect pair orientation logic though.
out = pd.DataFrame()
for fusion in rois.iterrows():
    bam_filter = \
    ((bam.rname == fusion[1].chromosome_x) & \
        (bam.rnext == fusion[1].chromosome_y) & \
        (bam.pos > fusion[1].startcoord_x) & \
        (bam.pos < fusion[1].endcoord_x) & \
        (bam.pnext > fusion[1].startcoord_y) & \
        (bam.pnext < fusion[1].endcoord_y) & \
        (bam.mapq >= 30)) | \
    ((bam.rname == fusion[1].chromosome_y) & \
        (bam.rnext == fusion[1].chromosome_x) & \
        (bam.pos > fusion[1].startcoord_y) & \
        (bam.pos < fusion[1].endcoord_y) & \
        (bam.pnext > fusion[1].startcoord_x) & \
        (bam.pnext < fusion[1].endcoord_x) & \
        (bam.mapq >= 30))
    
    bam_fusion = bam[bam_filter]
    bam_fusion.reset_index(drop = True, inplace = True)
    bam_fusion = bam_fusion.assign(strand1_exp = fusion[1].strand_x)
    bam_fusion = bam_fusion.assign(strand2_exp = fusion[1].strand_y)
    bam_fusion = bam_fusion.assign(gene1 = fusion[1].gene_x)
    bam_fusion = bam_fusion.assign(gene2 = fusion[1].gene_y)
    out = pd.concat([out, bam_fusion])

out = out.fillna(value = 'empty')
out['strand1'] = pd.to_numeric(out.strand1)
out['strand2'] = pd.to_numeric(out.strand2)
out['strand1_exp'] = pd.to_numeric(out.strand1_exp)
out['strand2_exp'] = pd.to_numeric(out.strand2_exp)

# Pair orientation filter
#   Checks orientations against expected strand directions. Supplementary alignments (flag >= 2048)
#   are checked for same or complementary directions and pairs (flag < 2048) are checked for
#   discordant directions. None or '' values are replaced with 'empty' for clean awk later.
pairs = out[(out.flag < 2048)]
pairs_actual = (pairs.strand1 + pairs.strand2) % 2
pairs_exp = (pairs.strand1_exp + pairs.strand2_exp) % 2
pairs_keep = (pairs_actual != pairs_exp)

supps = out[(out.flag >= 2048)]
supps_actual = (supps.strand1 + supps.strand2) % 2
supps_exp = (supps.strand1_exp + supps.strand2_exp) % 2
supps_keep = (supps_actual == supps_exp)

out_filtered = pd.concat([pairs[pairs_keep], supps[supps_keep]])
out_filtered = out_filtered.fillna(value = 'empty').replace('', 'empty').sort_values('qname').reset_index(drop = True)

# Filters results based on mate MAPQ
#   Works in the same way as the SA filter, where all columns are searched for 'MQ' and then updated into
#   a single series that is then split and filtered by mate MAPQ.
mate_mapq = out_filtered[tag_col_names]
mate_mapq_filter = mate_mapq[mate_mapq.apply(lambda col: col.str.startswith('MQ'), axis = 1)]
mate_mapq = mate_mapq_filter.tag1
for tag in tag_col_names[1:]:
    mate_mapq.update((mate_mapq_filter[[tag]].squeeze()))

mate_mapq = pd.DataFrame(mate_mapq.str.split(':').to_list(), columns = ['tag', 'i', 'mq'])
mate_mapq.mq = pd.to_numeric(mate_mapq.mq)
mate_filter = (mate_mapq.mq >= 30).reset_index(drop = True)

out_filtered = out_filtered[mate_filter].reset_index(drop = True)

# Filter complementary duplicates
#   E.g. AB and BA reads. Relies on previous sorting. May not filter triplicate reads
#   successfully if reads are identical but another read is in the middle, e.g. AB BC BA.
#   Also filters rows that have identical rname/rnext/pos/pnext.
# TODO: identify and test triplicate. Might just need multi-level sorting, e.g. sort
#   by qname, then by rname/rnext or something. Could just use chr in tag3? e.g. SA tags?
# TODO: Make it work for more than pairs - doesn't work for BCR-ABL reads, for example.
#   Maybe hash them? e.g. make them into a string of chr1_100_chr2_200 and chr2_200_chr1_100
#   and then compare each row to that instead?
# TODO: I'm sure there's a better way to do this.
if not(out_filtered.dropna().empty):
    discard = [True]
    prev_row = out_filtered.iloc[0]
    out_filt_iter = out_filtered.itertuples()
    next(out_filt_iter)
    for row in out_filt_iter:
        if (((row.rname == prev_row.rnext) &
            (row.rnext == prev_row.rname) &
            (row.pos == prev_row.pnext) &
            (row.pnext == prev_row.pos)) |
            ((row.rname == prev_row.rname) &
            (row.rnext == prev_row.rnext) &
            (row.pos == prev_row.pos) &
            (row.pnext == prev_row.pnext))):
                discard.append(False)
        else:
            discard.append(True)
        prev_row = row
    out_filtered = out_filtered[discard].reset_index(drop = True)

# Write results to output
out_filtered.to_csv(out_path, sep = '\t', line_terminator = '\n', index = False)
