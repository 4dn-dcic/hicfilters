#!/usr/bin/env python

"""
Functions for processing BAM files containing aligned Hi-C reads.

This file contains definitions for two classes, Fragment and HiCPair,
as well as a number of helper functions.
"""

from __future__ import print_function

import pysam
import re
import argparse
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Restriction import *

#############################################
class Fragment(object):
    """
    Minimal class for restriction fragments.
    """
    def __init__(self, record, delim='\t'):
        data = record.split(delim)
        self.chrom = data[0]
        self.site1 = int(data[1])
        self.site2 = int(data[2])
        self.index = int(data[3])

    ##############################################
    def __str__(self):
        return '%s - [%d, %d] (%d)' % (self.chrom, self.site1, self.site2, self.index)

    ##############################################
    def __repr__(self):
        return str(self)

    ##############################################
    def __eq__(self, other):
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return NotImplemented

    ##############################################
    def __ne__(self, other):
        if type(other) is type(self):
            return not self == other
        return NotImplemented

    ##############################################
    def __hash__(self):
        return hash(tuple(sorted(self.__dict__.items())))

    ##############################################
    def __sub__(self, other):
        if type(other) is type(self):
            if self.chrom != other.chrom:
                return None
            else:
                return abs(self.index - other.index)
        return NotImplemented

##############################################
def reverseComplement(seq):
    """
    Reverse a string and return as a string.
    """
    return str(Seq(seq).reverse_complement())

##############################################
def getRestrictionEnzyme(res):
    """
    Get RestrictionType object for given enzyme name.
    """
    b = RestrictionBatch()
    b.add(res)
    return b.get(res)

##############################################
def getRestrictionSequence(res):
    """
    Get the recognition site for the given restriction enzyme.
    """
    enzyme = getRestrictionEnzyme(res)
    return enzyme.elucidate().replace('N','').replace('^','').replace('_','')

##############################################
def getLigationSequence(res):
    """
    Get the ligation junction sequence expected from a
    restriction enzyme.
    """
    enzyme = getRestrictionEnzyme(res)
    cutseq = enzyme.elucidate()
    seq1 = re.sub('[\^]', '', cutseq.split('_')[0])
    seq2 = re.sub('[_]', '', cutseq.split('^')[1]) 
    return (seq1 + seq2).strip('N')

##############################################
def getLigationSites(seq, ligseq, mid=True, offset=0):
    """
    Given a DNA sequence and a ligation junction, find set of positions of
    ligation sites.
    """
    half = len(ligseq) / 2   # should be an integer!
    if mid:                  # return 1-based positions of middle of ligation sites
        return [ x.start() + half + offset + 1 for x in re.finditer(ligseq, seq) ]
    else:
        return [ x.start() + offset + 1 for x in re.finditer(ligseq, seq) ]

##############################################
def reverseCigarString(cigar):
    """
    Given a CIGAR string, reverse it.
    """
    words = re.findall('\d{1,3}[A-Z]', cigar)
    return ''.join(words[::-1])

##############################################
def sortReads(reads, mapq_threshold=30):
    """
    Given a list of pysam.AlignedSegment objects, check that the list
    contains precisely one representative, and if so, return the list
    with the representative as the first element.
    """
    # Iterate over alignments, removing all unmapped or low-MAPQ supplementaries
    reprs = []
    newreads = []
    for i, r in enumerate(reads):
        if not r.is_supplementary:
            reprs.append(i)
        elif r.is_unmapped or r.mapping_quality <= mapq_threshold:
            continue
        newreads.append(r)
    # Return None if there are no representatives
    if len(reprs) != 1:
        return None
    else:
        return [ newreads[reprs[0]] ] + newreads[:reprs[0]] + newreads[reprs[0]+1:]

#############################################
def switchChimericAlignments(repr, supp):
    """
    Given a representative and a supplementary alignment, switch so that the
    supplementary is now representative. Namely:
     - Change the supplementary bit in the SAM flags of both alignments
     - Change the CIGAR strings of both alignments (soft becomes hard, vice versa)
     - Change the sequences of both alignments
    """
    # Turn bit 2048 on for representative and off for supplementary
    repr.flag = repr.flag | 2048
    supp.flag = supp.flag & ~ 2048

    # Change soft clips in representative CIGAR to hard clips
    re.sub('\d{1,3}S', lambda s: s.replace('S', 'H'), repr.cigarstring)

    # Change hard clips in supplementary CIGAR to soft clips
    re.sub('\d{1,3}H', lambda s: s.replace('H', 'S'), supp.cigarstring)

    # Replace sequence in supplementary alignment with full sequence
    full_seq = repr.query_sequence
    full_quals = repr.query_qualities
    readlength = len(full_seq)
    if repr.is_reverse:    # Reverse to form presented in FASTQ file
        full_seq = reverseComplement(full_seq)
        full_quals = full_quals[::-1] 
    if supp.is_reverse:    # Reverse again if necessary
        supp.query_sequence = reverseComplement(full_seq)
        supp.query_qualities = full_quals[::-1]
    else:
        supp.query_sequence = full_seq
        supp.query_qualities = full_quals

    # Mask sequence in representative alignment according to hard clips
    if repr.is_reverse:    # Reverse again, to match representative alignment
        full_seq = reverseComplement(full_seq)
        full_quals = full_quals[::-1]
    clipleft = int(re.findall('^\d{1,3}H', repr.cigarstring)[0][:-1])
    clipright = int(re.findall('\d{1,3}H$', repr.cigarstring)[0][:-1])
    repr.query_sequence = full_seq[clipleft:readlength-clipright]
    repr.query_qualities = full_quals[clipleft:readlength-clipright]

    # Switch places when returning
    return supp, repr

#############################################
def getChimericBreakpoints(repr, supp, site=None):
    """
    Given a representative alignment and a supplementary alignment, determine
    the 1-based position at which the split(s) occurred.
     - If site is None, then simply return the position at which the
       representative alignment is clipped.
     - If a ligation site position is given, iterate through the range
       of possible breakpoints, [repr_clip_3end, ..., supp_clip_5end],
       and determine, if any, coincides with a ligation site. If none does,
       then return the split position that is closest to the given site position.

    NOTE: this method implicitly assumes that there is only one supplementary alignment.
          This is true for the great majority of chimeric reads.

    Parameters
    ----------
    repr : pysam.AlignedSegment
        representative alignment of chimeric read
    supp : pysam.AlignedSegment
        supplementary alignment of chimeric read
    site : int or NoneType, optional (default: None)
        location of known ligation site, used to correct breakpoint

    Returns
    -------
    breakpoint : int
        position of breakpoint (either the clipped end of the representative,
        or, if site is not None, the position of the known ligation site)
    """
    seq = repr.query_sequence          # query sequence
    cigar_repr = repr.cigarstring      # CIGAR string for representative alignment
    cigar_supp = supp.cigarstring      # CIGAR string for alignment
    readlength = sum([ int(x) for x in re.split('[A-Z]',cigar_repr) if x.isdigit() ])
                                       # infer read length from representative CIGAR
    breakpoint = 0

    # if one of the alignments are to the reverse strand, get
    # reverse complement of sequence and reverse CIGAR string.
    # this allows consideration of all alignment positions relative
    # to a single strand (the forward strand)
    if repr.is_reverse:
        seq = reverseComplement(seq)
        cigar_repr = reverseCigarString(cigar_repr)
    if supp.is_reverse:
        cigar_supp = reverseCigarString(cigar_supp)

    # determine positions of soft/hard clipping which determine the
    # chimeric breakpoint.
    # - In "Hi-C" mode, bwa mem returns the portion of the read
    #   that maps the 5' end of the read as representative (soft-clipping
    #   of 5' end is still possible)
    clip_repr = re.findall('\d{1,3}S$', cigar_repr)[0]
    break_repr = readlength - int(clip_repr.split('S')[0]) + 1
    clip_supp = re.findall('^\d{1,3}H', cigar_supp)[0]
    break_supp = int(clip_supp.split('H')[0]) + 1

    # if site is not given, return the breakpoint corresponding to the
    # representative alignment
    if site is None:
        breakpoint = break_repr
        difference = 0
    # otherwise, iterate over range between break_repr and break_supp, trying
    # to find ligation site position
    else:
        breaks = sorted([ break_repr, break_supp ])
        if site in range(breaks[0], breaks[1]+1):
            breakpoint = site
            difference = min([ abs(x - breakpoint) for x in breaks ])
        elif breaks[0] > site:
            breakpoint = breaks[0]
            difference = 0
        else:
            breakpoint = breaks[1]
            difference = 0

    return breakpoint, difference

#############################################
def combineChimericAlignment(repr, supp, readlength, breakpoint):
    """
    Combine representative and supplementary alignments to give a list
    of tuples that enumerate, from 5' to 3', the sequence of alignment
    positions along the genome.
    """
    combined_aln = [ (i, None) for i in range(readlength) ]

    # get aligned pairs from representative alignment, reversing read letter positions
    # if the read mapped to reverse strand
    repr_aln = [ x for x in repr.get_aligned_pairs() if x[0] is not None ]
    if repr.is_reverse:
        last = repr_aln[-1][0]
        repr_aln = [ (last - x[0], x[1]) for x in repr_aln ][::-1]
    # get aligned pairs from supplementary alignment, reversing read letter positions
    # if the read mapped to reverse strand
    supp_aln = [ x for x in supp.get_aligned_pairs() if x[0] is not None ]
    if supp.is_reverse:
        last = supp_aln[-1][0]
        supp_aln = [ (last - x[0], x[1]) for x in supp_aln ][::-1]

    # add None-entries for clipping events on either end of supplementary
    supp_cigar = supp.cigarstring
    if supp.is_reverse:
        supp_cigar = reverseCigarString(supp_cigar)
    supp_clip5 = [ int(x.split('H')[0]) for x in re.findall('^\d{1,3}H', supp_cigar) ]
    supp_clip3 = [ int(x.split('H')[0]) for x in re.findall('\d{1,3}H$', supp_cigar) ]
    if len(supp_clip5) > 0:
        supp_aln = list(enumerate([None] * supp_clip5[0] + [ x[1] for x in supp_aln ]))
    if len(supp_clip3) > 0:
        supp_aln = list(enumerate([ x[1] for x in supp_aln ] + [None] * supp_clip3[0]))

    # "extend" alignments to incorporate adjacent genomic positions
    idx, pos = [ x for x in repr_aln if x[1] is not None ][0]
    if repr.is_reverse:
        for i in range(idx, readlength):
            repr_aln[i] = (repr_aln[i][0], pos)
            pos -= 1
    else:
        for i in range(idx, readlength):
            repr_aln[i] = (repr_aln[i][0], pos)
            pos += 1
    idx, pos = [ x for x in supp_aln if x[1] is not None ][-1]
    if supp.is_reverse:
        for i in range(idx, -1, -1):
            supp_aln[i] = (supp_aln[i][0], pos)
            pos += 1
    else:
        for i in range(idx, -1, -1):
            supp_aln[i] = (supp_aln[i][0], pos)
            pos -= 1

    # combine both alignments
    for i in range(readlength):
        if i < breakpoint - 1:
            combined_aln[i] = repr_aln[i]
        else:
            combined_aln[i] = supp_aln[i]

    return [ x for x in combined_aln if x[1] is not None ]

#############################################
def getFragments(chrom, posleft, posright, tabix):
    """
    Given a chromosome name, the genomic positions of the ends of an
    alignmentm and a pysam.TabixFile of restriction fragments,
    return a list of the fragments that intersect with the alignment.
    """
    return [ Fragment(x) for x in tabix.fetch(chrom, posleft, posright) ]

#############################################
def getPositionsContiguous(read):
    """
    Get the 5' and 3' end positions of a contiguously aligned read.
    """
    if read.is_reverse:
        return [read.reference_end, read.reference_start + 1]
    else:
        return [read.reference_start + 1, read.reference_end]

#############################################
def getNearestRestrictionSite(read, fragments):
    """
    Given a pysam.AlignedSegment and a list of Fragment objects, return
    the nearest downstream restriction site to the 5' end of the read.
    """
    if read.is_reverse:
        sites = sorted([ f.site1 for f in fragments ], reverse=True)
        for s in sites:
            if read.reference_end >= s:
                return s
    else:
        sites = sorted([ f.site2 for f in fragments ])
        for s in sites:
            if read.reference_start + 1 <= s:
                return s

#############################################
def getFragmentsContiguous(read, tabix, cutlength):
    """
    Given a pysam.AlignedSegment and a pysam.TabixFile of restriction fragments,
    return a list of the fragments that intersect with the alignment. 
    """
    # Get chromosome name, left position, and right position
    chrom = read.reference_name
    posleft = read.reference_start + 1
    posright = read.reference_end
    # Correct by restriction enzyme cutter-length, depending
    # on alignment orientation
    if read.is_reverse:
        posleft += cutlength
    else:
        posright -= cutlength
    # Get restriction fragments
    return getFragments(chrom, posleft, posright, tabix)

#############################################
def getPositionsChimeric(repr, supp, ligseq, readlength):
    """
    Given two pysam.AlignedSegment objects, corresponding to the two pieces of
    a chimeric alignment, and a pysam.TabixFile of restriction fragments,
    get the (corrected) 5' and 3' end genomic positions of the two alignments,
    as follows:
      1) Locate any ligation sites along the read.
      2) If there is a ligation site, correct the breakpoint between
         the two alignments. If there is 
      3) Combine the two alignments at the corrected breakpoint
      4) Get the start/end positions of the two alignments from the
         combined alignment
    """
    # Locate all ligation sites along the representative sequence
    seq = repr.query_sequence
    if repr.is_reverse:
        seq = reverseComplement(seq)
    ligsites = getLigationSites(seq, ligseq, mid=True, offset=0)

    # If there are no ligation sites, return original breakpoint uncorrected
    if len(ligsites) == 0:
        breakpoint, _ = getChimericBreakpoints(repr, supp, site=None)
    # If there is a single ligation site, correct breakpoint
    elif len(ligsites) == 1:
        breakpoint, _ = getChimericBreakpoints(repr, supp, site=ligsites[0])
    # If there are >1 ligation sites, test all of them, and choose the
    # inferred breakpoint that is closest to the reported breakpoints
    # (i.e., has minimum difference)
    else:
        breakpoints = [ getChimericBreakpoints(repr, supp, site=site) for site in ligsites ]
        breakpoints.sort(key=lambda x: x[1])
        breakpoint = breakpoints[0][0]

    # Combine the two alignments
    combined_aln = combineChimericAlignment(repr, supp, readlength, breakpoint)

    # Obtain aligned positions of the two ends of each alignments
    repr_pos5 = combined_aln[0][1] + 1
    repr_pos3 = [ x[1] for x in combined_aln if x[0] == breakpoint-2 ][0] + 1
    supp_pos5 = [ x[1] for x in combined_aln if x[0] == breakpoint-1 ][0] + 1
    supp_pos3 = combined_aln[-1][1] + 1

    return repr_pos5, repr_pos3, supp_pos5, supp_pos3

#############################################
def getFragmentsChimeric(repr, supp, tabix, ligseq, readlength, cutlength,
                         return_positions=False):
    """
    Given a pair of pysam.AlignedSegment objects forming a chimeric alignment,
    along with a pysam.TabixFile of restriction fragments, determine the
    set of fragments that intersect with the two alignments. 
    """
    # Obtain 5' and 3' positions of both alignments
    repr_pos5, repr_pos3, supp_pos5, supp_pos3 =\
        getPositionsChimeric(repr, supp, ligseq, readlength)

    # Re-sort the aligned positions and get chromosome names
    repr_posleft, repr_posright = sorted([ repr_pos5, repr_pos3 ])
    supp_posleft, supp_posright = sorted([ supp_pos5, supp_pos3 ])
    repr_chrom = repr.reference_name
    supp_chrom = supp.reference_name

    # Correct positions by restriction enzyme cutter-length, depending on
    # alignment orientation
    if repr.is_reverse:
        repr_posleft += cutlength
    else:
        repr_posright -= cutlength
    if supp.is_reverse:
        supp_posleft += cutlength
    else:
        supp_posright -= cutlength

    # Get restriction fragments
    frags_repr = getFragments(repr_chrom, repr_posleft, repr_posright, tabix)
    frags_supp = getFragments(supp_chrom, supp_posleft, supp_posright, tabix)

    if return_positions:
        return frags_repr, frags_supp, repr_pos5, repr_pos3, supp_pos5, supp_pos3

    return frags_repr, frags_supp

#############################################
class HiCPair(object):
    """
    Class definitions for a HiCPair object.
    """
    def __init__(self, reads=None, chrs='all', ligseq=None, cutlength=None,
                 mapq_threshold=30, sep_threshold=1000, res_site_dist_threshold=750,
                 min_insert_size=0, max_insert_size=1000, fragments_file=None):
        """
        Constructor for the HiCPair object.

        PARAMETERS
        ----------
        reads : list of pysam.AlignedSegment objects
            Alignments of either read in pair.
        chrs : list of strings, or 'all'
            List of chromosomes.
        ligseq : str
            Name of restriction enzyme used for digestion. Use None
            for DNase Hi-C and Micro-C.
        cutlength : int
            Length of the restriction enzyme site motif. Use None
            for DNase Hi-C and Micro-C. 
        mapq_threshold : int
            Non-negative integer, ranging between 0 and 255. An
            alignment with MAPQ less than this number is flagged.
        sep_threshold : int
            Non-negative integer. A contact with separation less than
            this number is flagged.
        res_site_dist_threshold : int
            Non-negative integer. A pair in which either read has a
            5' end located further than this number from its nearest
            restriction site is flagged.
        min_insert_size : int
            Non-negative integer. A pair whose inferred insert size is
            less than this number is flagged.
        max_insert_size : int
            Non-negative integer. A pair whose inferred insert size is
            greater than this number is flagged.
        fragments_file : pysam.TabixFile
            TabixFile object for the .bed.gz file containing restriction
            fragment information.
        """
        self.reads1 = []
        self.reads2 = []
        self.pos1 = None
        self.pos2 = None
        self.readlength = None
        self.nearest_site1 = None
        self.nearest_site2 = None
        self.insert_size = None
        self.ligseq = ligseq
        self.cutlength = cutlength
        self.separations = None
        self.primary_sep = None

        # Initialize all flags to False, with exception of is_empty
        self.is_empty = True
        self.is_unmapped = False
        self.low_MAPQ = False
        self.other_contig = False
        self.multiple_loci = False
        self.close_pair = False
        self.multiple_fragments = False
        self.undigested_site = False
        self.same_fragment = False
        self.far_from_restriction_site = False
        self.wrong_insert_size = False

        # Open tabix file of restriction fragments, if specified
        self.fragments = None
        if self.ligseq is not None and fragments_file is None:
            raise Exception('Fragments file should be specified for regular Hi-C data')
        elif self.ligseq is None and fragments_file is not None:
            raise Exception('Fragments file specified without restriction enzyme')
        elif fragments_file is not None:
            self.fragments = fragments_file

        # Initialize list of standard chromosomes/contigs
        if chrs == 'all':
            chrs = [ 'chr%d' % i for i in range(1,23) ] + ['chrX', 'chrY']

        # Add alignments if specified
        if reads is not None:
            self.addReads(reads, mapq_threshold)

        # If nonempty, test for mappability
        if not self.is_empty:
            self.flagUnmappedPair()

            # If nonempty and unmapped, test for low-MAPQ and/or mapping to other contigs
            if not self.is_unmapped:
                self.flagLowMAPQPair(mapq_threshold=mapq_threshold)
                self.flagOtherContig(chrs)

                # If reads mapped with high MAPQ to a standard contig
                if not self.low_MAPQ and not self.other_contig:

                    # Compute genomic separations
                    self.computePositions()

                    # Flag pair if there are > 3 alignments, or there are 3 alignments
                    # and the minimum distance b/t each pair exceeds sep_threshold
                    self.flagMultipleLoci(sep_threshold=sep_threshold)

                    # Flag pair if separation is less than empirical threshold
                    # (<10kb for 6-cutter, <1kb for 4-cutter, respectively)
                    self.flagClosePair()

                    # Flag pair if reads map to multiple restriction fragments
                    if self.ligseq is not None:
                        self.assignFragments()

        # Update BAM tags
        try:
            self.updateBAMTags()
        except:
            pass

    #############################################
    def addReads(self, reads, mapq_threshold):
        """
        Add alignment data (lists of pysam.AlignedSegment objects).
        """
        reads1 = sortReads([ r for r in reads if r.is_read1 ],
                           mapq_threshold=mapq_threshold)
        reads2 = sortReads([ r for r in reads if r.is_read2 ],
                           mapq_threshold=mapq_threshold)
        if reads1 is not None and reads2 is not None:
            self.reads1 = reads1
            self.reads2 = reads2
            self.readlength = len(self.reads1[0].query_sequence)
            self.is_empty = False

    #############################################
    def clearReads(self):
        """
        Delete alignment data from HiCPair object and reset all flags and metadata.
        """
        self.reads1 = None
        self.reads2 = None
        self.pos1 = None
        self.pos2 = None
        self.readlength = None
        self.nearest_site1 = None
        self.nearest_site2 = None
        self.insert_size = None
        self.separations = None
        self.primary_sep = None
        self.is_empty = True
        self.is_unmapped = False
        self.other_contig = False
        self.multiple_loci = False
        self.close_pair = False
        self.multiple_fragments = False
        self.undigested_site = False
        self.same_fragment = False
        self.far_from_restriction_site = False
        self.wrong_insert_size = False

    #############################################
    def updateBAMTags(self):
        """
        Update tags for each alignment record in self.reads1 and self.reads2.

        COMPLETE LIST OF TAGS
        ---------------------
        BOOLEAN TAGS : YM - whether either read was unmapped or poorly mapped
                       YL - whether the reads mapped to multiple distal loci
                       YC - whether the reads represent a very close-range contact
                       YF - whether the reads mapped to multiple restriction fragments
                       YU - whether the reads contain an undigested restriction site
                       YS - whether the reads mapped to the same restriction fragment
                       YI - whether the reads come from an insert with anomalous length
                       YR - whether either read mapped very far from the nearest restriction site
        INTEGER TAGS : Z1 - position of read 1
                       Z2 - position of read 2
                       ZM - position of nearest restriction site downstream of read 1
                       ZN - position of nearest restriction site downstream of read 2
        """
        if self.is_empty:
            raise Exception('Cannot write empty HiCPair to file')
        for r in self.reads1 + self.reads2:
            r.set_tag('YM', int(self.is_unmapped or self.low_MAPQ or self.other_contig), value_type='i')
            r.set_tag('YL', int(self.multiple_loci), value_type='i')
            r.set_tag('YC', int(self.close_pair), value_type='i')
            r.set_tag('YF', int(self.multiple_fragments), value_type='i')
            r.set_tag('YU', int(self.undigested_site), value_type='i')
            r.set_tag('YS', int(self.same_fragment), value_type='i')
            r.set_tag('YI', int(self.wrong_insert_size), value_type='i')
            r.set_tag('YR', int(self.far_from_restriction_site), value_type='i')
            r.set_tag('Z1', self.pos1, value_type='i')
            r.set_tag('Z2', self.pos2, value_type='i')
            r.set_tag('ZM', self.nearest_site1, value_type='i')
            r.set_tag('ZN', self.nearest_site2, value_type='i')

    #############################################
    def writeToFile(self, outbam):
        """
        Write alignment records to a BAM file (pysam.AlignmentFile).
        """
        for r in self.reads1 + self.reads2:
            outbam.write(r)

    #############################################
    def computePositions(self):
        """
        Determine 5' positions of the alignments.
        """
        # If either read is empty, unmapped, low-MAPQ, or mapped to non-standard contig
        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.pos1 = None
            self.pos2 = None
            self.separations = None
        # If both reads mapped contiguously
        elif len(self.reads1) == 1 and len(self.reads2) == 1:
            self.pos1, _ = getPositionsContiguous(self.reads1[0])
            self.pos2, _ = getPositionsContiguous(self.reads2[0])
            self.separations = [ abs(self.pos1 - self.pos2) ]
        # If first read mapped chimerically
        elif len(self.reads1) == 2 and len(self.reads2) == 1:
            repr_pos5, _, supp_pos5, _ = getPositionsChimeric(self.reads1[0], self.reads1[1],
                                                              self.ligseq, self.readlength)
            self.pos1 = repr_pos5
            self.pos2, _ = getPositionsContiguous(self.reads2[0])
            self.separations = [ abs(self.pos1 - supp_pos5),
                                 abs(self.pos1 - self.pos2),
                                 abs(supp_pos5 - self.pos2) ]
        # If second read mapped chimerically
        elif len(self.reads1) == 1 and len(self.reads2) == 2:
            self.pos1, _ = getPositionsContiguous(self.reads1[0])
            repr_pos5, _, supp_pos5, _ = getPositionsChimeric(self.reads2[0], self.reads2[1],
                                                              self.ligseq, self.readlength)
            self.pos2 = repr_pos5
            self.separations = [ abs(self.pos1 - self.pos2),
                                 abs(self.pos1 - supp_pos5),
                                 abs(self.pos2 - supp_pos5) ]
        # If the read pair was otherwise (anomalously) aligned
        else:
            self.pos1 = None
            self.pos2 = None
            self.separations = None

        # Distance between representatives
        if self.pos1 is not None and self.pos2 is not None:
            self.primary_sep = abs(self.pos1 - self.pos2)
        else:
            self.primary_sep = None

    #############################################
    def flagUnmappedPair(self):
        """
        Sets self.is_unmapped to True if either read is unmapped. Only checks
        representatives.
        """
        if self.is_empty:
            self.is_unmapped = False
        else:
            unmapped1 = self.reads1[0].is_unmapped
            unmapped2 = self.reads2[0].is_unmapped
            # If the representative is unmapped while there is a supplementary,
            # replace representative with next supplementary
            if unmapped1 and len(self.reads1) > 1:
                repr, supp = switchChimericAlignments(self.reads1[0], self.reads1[1])
                self.reads1 = [repr, supp] + self.reads1[2:]
                unmapped1 = False
            # Do the same for the second read
            if unmapped2 and len(self.reads2) > 1:
                repr, supp = switchChimericAlignments(self.reads2[0], self.reads2[1])
                self.reads2 = [repr, supp] + self.reads2[2:]
                unmapped2 = False
            self.is_unmapped = (unmapped1 or unmapped2)

    #############################################
    def flagLowMAPQPair(self, mapq_threshold=30):
        """
        Sets self.low_MAPQ to True if both reads have been mapped, but either
        read has MAPQ at or less than <mapq_threshold>. Only checks representatives.
        """
        # Default to False if either read did not map
        if self.is_empty or self.is_unmapped:
            self.low_MAPQ = False
        # If both reads were mapped successfully
        else:
            low_mapq1 = (self.reads1[0].mapping_quality <= mapq_threshold)
            low_mapq2 = (self.reads2[0].mapping_quality <= mapq_threshold)
            # If either representative has low MAPQ, replace with
            # supplementary
            if low_mapq1 and len(self.reads1) > 1:
                repr, supp = switchChimericAlignments(self.reads1[0], self.reads1[1])
                self.reads1 = [repr, supp] + self.reads1[2:]
                low_mapq1 = False
            # Do the same for second read
            if low_mapq2 and len(self.reads2) > 1:
                repr, supp = switchChimericAlignments(self.reads2[0], self.reads2[1])
                self.reads2 = [repr, supp] + self.reads2[2:]
                low_mapq2 = False
            self.low_MAPQ = (low_mapq1 or low_mapq2)

    #############################################
    def flagOtherContig(self, chrs):
        """
        Sets self.other_contig to True if any one of the alignments mapped
        to a contig not in <chrs>.
        """
        if self.is_empty or self.is_unmapped:
            self.other_contig = False
        else:
            self.other_contig = (sum([ r.reference_name not in chrs for r in self.reads1 ]) > 0 or
                                 sum([ r.reference_name not in chrs for r in self.reads2 ]) > 0)

    #############################################
    def flagMultipleLoci(self, sep_threshold=1000):
        """
        Sets self.multiple_loci to True if there are > 3 alignments for the
        pair, or there are 3 alignments and the minimum distance b/t each
        pair of alignments exceeds <sep_threshold>.
        """
        # Default to False if read pair is empty or unmapped
        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.multiple_loci = False
        # If there are two alignments, set to False
        elif len(self.reads1) == 1 and len(self.reads2) == 1:
            self.multiple_loci = False
        # If there are three alignments, look for minimum separation
        # between each pair of alignments 
        elif (len(self.reads1) == 2 and len(self.reads2) == 1) or\
             (len(self.reads1) == 1 and len(self.reads2) == 2):
            if self.separations is None:
                self.computePositions()
            min_dist = min(self.separations)
            self.multiple_loci = (min_dist > sep_threshold)
        # If there are > 3 alignments, set to True
        else:
            self.multiple_loci = True

    #############################################
    def flagClosePair(self):
        """
        Sets self.close_pair to True if abs(self.pos1 - self.pos2) is less than
        or equal to a threshold based on the cutter length.
        """
        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.close_pair = False
        
        thresholds = { 4 : 1000, 6 : 10000 }
        if self.cutlength is None:
            self.close_pair = False
        elif self.pos1 is None or self.pos2 is None:
            self.close_pair = False
        else:
            self.close_pair = (self.primary_sep <= thresholds[self.cutlength])

    #############################################
    def assignFragments(self, res_site_dist_threshold=750, min_insert_size=0,
                        max_insert_size=1000):
        """
        Determines all restriction fragments that intersect with each
        alignment. Sets all fragment-related flags as necessary.
        """
        ###################
        # BOTH CONTIGUOUS #
        ###################
        def assignFragmentsBothContiguous():
            """
            Assign fragments in the case where both reads mapped contiguously.
            """
            self.multiple_fragments = False
            frags1 = getFragmentsContiguous(self.reads1[0], self.fragments, self.cutlength)
            frags2 = getFragmentsContiguous(self.reads2[0], self.fragments, self.cutlength)
            # If either read intersects with > 1 fragments, flag for undigested site
            if len(frags1) > 1 or len(frags2) > 1:
                self.undigested_site = True
                self.same_fragment = False
            # If both reads intersect with the same fragment, flag for same-fragment pair
            elif frags1[0] == frags2[0]:
                self.undigested_site = False
                self.same_fragment = True
            # If both reads intersect with unique distinct fragments
            else:
                self.undigested_site = False
                self.same_fragment = False
            return frags1, frags2

        ################
        # ONE CHIMERIC #
        ################
        def assignFragmentsOneChimeric(repr, supp, other):
            """
            Assign fragments in the case where one read mapped chimerically.
            """
            frags_repr, frags_supp = getFragmentsChimeric(repr, supp, self.fragments, self.ligseq,
                                                          self.readlength, self.cutlength)
            # Assign fragments for alignment of second read
            frags_other = getFragmentsContiguous(other, self.fragments, self.cutlength)
            # If any of the alignments map to more than one fragment
            if len(frags_repr) > 1 or len(frags_supp) > 1 or len(frags_other) > 1:
                self.undigested_site = True
                self.multiple_fragments = True
                self.same_fragment = False
            # If all of the alignments map to single fragments, and the supplementary
            # and contiguous alignments map to the same fragment
            elif frags_supp[0] == frags_other[0]:
                self.undigested_site = False
                self.multiple_fragments = False
                # If, in addition, the representative also maps to the same fragment
                if frags_repr[0] == frags_other[0]:
                    self.same_fragment = True
                else:
                    self.same_fragment = False
            # Otherwise, the alignments all map uniquely to distinct fragments
            else:
                self.undigested_site = False
                self.multiple_fragments = True
                self.same_fragment = False
            return frags_repr, frags_other

        #####################################
        # RESTRICTION SITES AND INSERT SIZE #
        #####################################
        def computeNearestSitesAndInsertSize(read1, read2, frags1, frags2):
            """
            Given alignments for the first and second reads, and their assigned
            fragments, locate the nearest restriction site to each read, and compute
            the expected insert size.
            """
            self.nearest_site1 = getNearestRestrictionSite(read1, frags1)
            self.nearest_site2 = getNearestRestrictionSite(read2, frags2)
            neardist1 = abs(self.nearest_site1 - self.pos1)
            neardist2 = abs(self.nearest_site2 - self.pos2)
            if neardist1 > res_site_dist_threshold or neardist2 > res_site_dist_threshold:
                self.far_from_restriction_site = True
            else:
                self.far_from_restriction_site = False
            # Compute insert size
            self.insert_size = neardist1 + neardist2
            if self.insert_size < min_insert_size or self.insert_size > max_insert_size:
                self.wrong_insert_size = True
            else:
                self.wrong_insert_size = False

        #####################################
        # If pair is empty, unmapped, low-MAPQ, or maps to non-standard contig,
        # don't bother with fragment assignment
        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.multiple_fragments = False
            self.undigested_site = False
            self.same_fragment = False
            self.far_from_restriction_site = False
            self.wrong_insert_size = False

        # If both reads mapped contiguously, then reads cannot intersect with
        # more than 2 distal fragments, save for the possibility of undigested
        # restriction sites
        elif len(self.reads1) == 1 and len(self.reads2) == 1:
            frags1, frags2 = assignFragmentsBothContiguous()
            computeNearestSitesAndInsertSize(self.reads1[0], self.reads2[0], frags1, frags2)

        # If the first read was chimerically aligned
        elif len(self.reads1) == 2 and len(self.reads2) == 1:
            frags1, frags2 = assignFragmentsOneChimeric(self.reads1[0], self.reads1[1], self.reads2[0])
            computeNearestSitesAndInsertSize(self.reads1[0], self.reads2[0], frags1, frags2)

        # If the second read was chimerically aligned
        elif len(self.reads1) == 1 and len(self.reads2) == 2:
            frags2, frags1 = assignFragmentsOneChimeric(self.reads2[0], self.reads2[1], self.reads1[0])
            computeNearestSitesAndInsertSize(self.reads1[0], self.reads2[0], frags1, frags2)

        # Otherwise, one or both reads map to multiple loci
        else:
            self.multiple_fragments = True
            self.undigested_site = False
            self.same_fragment = False
                
#############################################
####### MAIN SCRIPT FOR DATA ANALYSIS #######
#############################################
def parseBAM(bam, outbam, tabix, enzyme='MboI', mapq_threshold=30,
             sep_threshold=1000, verbose=True):
    """
    Parse input BAM file, process each pair and flag with necessary filters,
    and re-write each pair to an output BAM file.
    """
    curr_id = None
    curr_reads = []
    finished = False
    count = 0
    metadata = defaultdict(int)
    fragments_file = pysam.TabixFile(tabix)
    ligseq = getLigationSequence(enzyme)
    cutlength = len(getRestrictionEnzyme(enzyme))

    with pysam.AlignmentFile(bam, 'rb') as bf:
        outbf = pysam.AlignmentFile(outbam, 'wb', template=bf)
        while not finished:
            # Try to fetch the next read in the BAM file
            try:
                read = bf.next()
            except StopIteration:
                finished = True
            # Initialize curr_id if necessary
            if curr_id is None:
                curr_id = read.query_name
            # If we have reached the end of the file or we have reached
            # a new read ID, then process gathered alignment data
            if finished or curr_id != read.query_name:
                pair = HiCPair(reads=curr_reads, ligseq=ligseq,
                               cutlength=cutlength,
                               mapq_threshold=mapq_threshold,
                               sep_threshold=sep_threshold,
                               fragments_file=fragments_file)
                if pair.is_unmapped:               metadata['unmapped'] += 1
                if pair.low_MAPQ:                  metadata['low_MAPQ'] += 1
                if pair.other_contig:              metadata['other_contig'] += 1
                if pair.multiple_loci:             metadata['multiple_loci'] += 1
                if pair.close_pair:                metadata['close_pair'] += 1
                if pair.undigested_site:           metadata['undigested_site'] += 1
                if pair.multiple_fragments:        metadata['multiple_fragments'] += 1
                if pair.same_fragment:             metadata['same_fragment'] += 1
                if pair.far_from_restriction_site: metadata['far_from_restriction_site'] += 1
                if pair.wrong_insert_size:         metadata['wrong_insert_size'] += 1
                pair.writeToFile(outbf)
                # Update curr_id and re-initialize curr_reads
                curr_id = read.query_name
                curr_reads = []
                count += 1
                if count % 10000 == 0 and verbose:
                    print('Processed %d read pairs. Flag statistics:' % count)
                    print(dict(metadata))

            # Finally, gather alignment data from current read
            curr_reads.append(read)

        outbf.close()

    print('Finished processing %d read pairs. Flag statistics:' % count)
    print(dict(metadata))

    fragments_file.close()

#############################################
def parseInputs():
    """
    Parse input arguments (input BAM file and output BAM file).
    """
    parser = argparse.ArgumentParser(description='Filter Hi-C reads in BAM files.')
    parser.add_argument('-i', '--in', nargs=1, required=True,
                        help='Input BAM file')
    parser.add_argument('-o', '--out', nargs=1, required=True,
                        help='Output BAM file')
    parser.add_argument('-f', '--frag', nargs=1, default=None, required=False,
                        help='Tabix file for restriction fragments')
    parser.add_argument('-v', '--verbose', action='store_const',
                        const=True, default=False, required=False,
                        help='Verbose mode')
    args = vars(parser.parse_args())
    bam = args['in'][0]
    outbam = args['out'][0]
    tabix = args['frag'][0]
    verbose = args['verbose']
    return bam, outbam, tabix, verbose

#############################################
def main():
    bam, outbam, tabix, verbose = parseInputs()
    parseBAM(bam, outbam, tabix, verbose=verbose)

#############################################
if __name__ == '__main__':
    main()
