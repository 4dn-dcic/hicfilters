#!/usr/bin/env python

"""
Functions for processing BAM files containing aligned Hi-C reads.

This file contains definitions for three classes, Fragment, HiCPair,
and HiCDataset, as well as a number of helper functions.

Author: Kee-Myoung (Chris) Nam
        Department of Systems Biology
        Harvard Medical School
Last updated: 4/23/2017
"""

from __future__ import print_function

import os
import shutil
import pysam
import re
import math
import argparse
from collections import defaultdict
from Bio.Seq import Seq
from Bio.Restriction import *

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pyupset as pyu
from matplotlib_venn import venn2, venn3

chrom_lengths = { 'chr1'  : 249250621,
                  'chr2'  : 243199373,
                  'chr3'  : 198022430,
                  'chr4'  : 191154276,
                  'chr5'  : 180915260,
                  'chr6'  : 171115067,
                  'chr7'  : 159138663,
                  'chr8'  : 146364022,
                  'chr9'  : 141213431,
                  'chr10' : 135534747,
                  'chr11' : 135006516,
                  'chr12' : 133851895,
                  'chr13' : 115169878,
                  'chr14' : 107349540,
                  'chr15' : 102531392,
                  'chr16' : 90354753,
                  'chr17' : 81195210,
                  'chr18' : 78077248,
                  'chr19' : 59128983,
                  'chr20' : 63025520,
                  'chr21' : 48129895,
                  'chr22' : 51304566,
                  'chrX'  : 155270560,
                  'chrY'  : 155270560 }
genome_length = sum(chrom_lengths.values())

##############################################
class InvalidInputException(Exception):
    pass

class EmptyPairException(Exception):
    pass

class BadBreakpointException(Exception):
    pass

class NoNearestSiteException(Exception):
    pass

##############################################
class Fragment(object):
    """
    Minimal class for restriction fragments, defining the following
    binary operations:
      - Equality   : two fragments are "equal" if they are the same (i.e., same
                     index, same start position, same end position)
      - Difference : the "difference" of two fragments is the difference between
                     their indices (undefined if fragments lie on different
                     chromosomes)
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
########## VARIOUS HELPER FUNCTIONS ##########
##########  FOR MANIPULATING READS  ##########
##############################################
def reverse_complement(seq):
    """
    Get reverse complement of a given sequence.
    """
    return str(Seq(seq).reverse_complement())

##############################################
def get_restriction_enzyme(res):
    """
    Get RestrictionType object for given enzyme name.
    """
    b = RestrictionBatch()
    b.add(res)
    return b.get(res)

##############################################
def get_ligation_sequence(res):
    """
    Get the ligation junction sequence expected from a restriction enzyme.
    """
    enzyme = get_restriction_enzyme(res)
    cutseq = enzyme.elucidate()
    seq1 = re.sub('[\^]', '', cutseq.split('_')[0])
    seq2 = re.sub('[_]', '', cutseq.split('^')[1]) 
    return (seq1 + seq2).strip('N')

##############################################
def get_ligation_sites(seq, ligseq, mid=True, offset=0):
    """
    Given a DNA sequence and a ligation junction, find set of 1-based
    positions of ligation sites.
    """
    half = len(ligseq) / 2   # should be an integer!
    if mid:                  # return 1-based positions of middle of ligation sites
        return [ x.start() + half + offset + 1 for x in re.finditer(ligseq, seq) ]
    else:
        return [ x.start() + offset + 1 for x in re.finditer(ligseq, seq) ]

##############################################
def reverse_cigar_string(cigar):
    """
    Return reversed CIGAR string.
    """
    words = re.findall('\d{1,3}[A-Z]', cigar)
    return ''.join(words[::-1])

##############################################
def sort_reads(reads, mapq_threshold=30):
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
def switch_chimeric_alignments(repr, supp):
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
        full_seq = reverse_complement(full_seq)
        full_quals = full_quals[::-1] 
    if supp.is_reverse:    # Reverse again if necessary
        supp.query_sequence = reverse_complement(full_seq)
        supp.query_qualities = full_quals[::-1]
    else:
        supp.query_sequence = full_seq
        supp.query_qualities = full_quals

    # Mask sequence in representative alignment according to hard clips
    if repr.is_reverse:    # Reverse again, to match representative alignment
        full_seq = reverse_complement(full_seq)
        full_quals = full_quals[::-1]
    clipleft = int(re.findall('^\d{1,3}H', repr.cigarstring)[0][:-1])
    clipright = int(re.findall('\d{1,3}H$', repr.cigarstring)[0][:-1])
    repr.query_sequence = full_seq[clipleft:readlength-clipright]
    repr.query_qualities = full_quals[clipleft:readlength-clipright]

    # Switch places when returning
    return supp, repr

#############################################
def get_chimeric_breakpoints(repr, supp, site=None):
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
    difference : int
        difference between located breakpoint and the given ligation site
    """
    seq = repr.query_sequence          # query sequence
    cigar_repr = repr.cigarstring      # CIGAR string for representative alignment
    cigar_supp = supp.cigarstring      # CIGAR string for alignment
    readlength = sum([ int(x) for x in re.split('[A-Z]',cigar_repr) if x.isdigit() ])
                                       # infer read length from representative CIGAR
    breakpoint = 0

    # If one of the alignments are to the reverse strand, get
    # reverse complement of sequence and reverse CIGAR string.
    # This allows consideration of all alignment positions relative
    # to a single strand (the forward strand)
    if repr.is_reverse:
        seq = reverse_complement(seq)
        cigar_repr = reverse_cigar_string(cigar_repr)
    if supp.is_reverse:
        cigar_supp = reverse_cigar_string(cigar_supp)

    # Determine positions of soft/hard clipping which determine the
    # chimeric breakpoint.
    # - In "Hi-C" mode, bwa mem returns the portion of the read
    #   that maps the 5' end of the read as representative (soft-clipping
    #   of 5' end is still possible)
    clip_repr = re.findall('\d{1,3}S$', cigar_repr)[0]
    break_repr = readlength - int(clip_repr.split('S')[0]) + 1
    clip_supp = re.findall('^\d{1,3}H', cigar_supp)[0]
    break_supp = int(clip_supp.split('H')[0]) + 1

    # If site is not given, return the breakpoint corresponding to the
    # representative alignment
    if site is None:
        breakpoint = break_repr
        difference = 0
    # Otherwise, iterate over range between break_repr and break_supp, trying
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
def combine_chimeric_alignment(repr, supp, readlength, breakpoint):
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
        supp_cigar = reverse_cigar_string(supp_cigar)
    supp_clip5 = [ int(x.split('H')[0]) for x in re.findall('^\d{1,3}H', supp_cigar) ]
    supp_clip3 = [ int(x.split('H')[0]) for x in re.findall('\d{1,3}H$', supp_cigar) ]
    if len(supp_clip5) > 0:
        supp_aln = list(enumerate([None] * supp_clip5[0] + [ x[1] for x in supp_aln ]))
    if len(supp_clip3) > 0:
        supp_aln = list(enumerate([ x[1] for x in supp_aln ] + [None] * supp_clip3[0]))

    # "extend" alignments to incorporate adjacent genomic positions
    idx, pos = next( x for x in repr_aln if x[1] is not None )
    if repr.is_reverse:
        for i in range(idx, readlength):
            repr_aln[i] = (repr_aln[i][0], pos)
            pos -= 1
    else:
        for i in range(idx, readlength):
            repr_aln[i] = (repr_aln[i][0], pos)
            pos += 1
    idx, pos = next( x for x in supp_aln[::-1] if x[1] is not None )
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
def get_fragments(chrom, posleft, posright, tabix):
    """
    Given a chromosome name, the genomic positions of the ends of an
    alignmentm and a pysam.TabixFile of restriction fragments,
    return a list of the fragments that intersect with the alignment.
    """
    return set([ Fragment(x) for x in tabix.fetch(chrom, posleft, posright) ])

#############################################
def get_positions_contiguous(read):
    """
    Get the 5' and 3' end positions of a contiguously aligned read.
    """
    if read.is_reverse:
        return [read.reference_end, read.reference_start + 1]
    else:
        return [read.reference_start + 1, read.reference_end]

#############################################
def get_nearest_restriction_site(read, fragments):
    """
    Given a pysam.AlignedSegment and a list of Fragment objects, return
    the nearest downstream restriction site to the 5' end of the read.
    """
    # Determine position of 5' end of alignment
    if not read.is_reverse:
        pos = read.reference_start + 1
    else:
        pos = read.reference_end

    # If the read mapped to reverse strand, identify closest restriction site
    # to the left of the 5' alignment position
    if read.is_reverse:
        sites = [ f.site1 for f in fragments if f.site1 <= pos ]
        if sites is None or sites == []:
            raise NoNearestSiteException('No downstream restriction site: pos %d, fragments %s, reverse strand' %\
                                         (pos, fragments))
        return max(sites)
    # If the read mapped to forward strand, identify closest restriction site
    # to the right of the 5' alignment position
    else:
        sites = [ f.site2 for f in fragments if f.site2 >= pos ]
        if sites is None or sites == []:
            raise NoNearestSiteException('No downstream restriction site: pos %d, fragments %s, forward strand' %\
                                         (pos, fragments))
        return min(sites)

#############################################
def get_fragments_contiguous(read, tabix, cutlength, return_positions=True):
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
    posleft += cutlength
    posright -= cutlength

    # Get restriction fragments, returning positions if desired
    pos5 = posright if read.is_reverse else posleft
    pos3 = posleft if read.is_reverse else posright
    frags = get_fragments(chrom, posleft, posright, tabix)
    frag5 = get_fragments(chrom, pos5, pos5+1, tabix)
    if len(frag5) == 0:
        print(pos5, pos3, frags)
    (frag5,) = frag5

    if return_positions:
        return frags, frag5, pos5, pos3
    return frags, frag5

#############################################
def get_positions_chimeric(repr, supp, ligseq, readlength):
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
        seq = reverse_complement(seq)
    ligsites = get_ligation_sites(seq, ligseq, mid=True, offset=0)

    # If there are no ligation sites, raise Exception
    if len(ligsites) == 0:
        raise BadBreakpointException('Bad breakpoint found in chimeric read')
    # If there is a single ligation site, correct breakpoint
    elif len(ligsites) == 1:
        breakpoint, _ = get_chimeric_breakpoints(repr, supp, site=ligsites[0])
    # If there are >1 ligation sites, test all of them, and choose the
    # inferred breakpoint that is closest to the reported breakpoints
    # (i.e., has minimum difference)
    else:
        breakpoints = [ get_chimeric_breakpoints(repr, supp, site=site) for site in ligsites ]
        breakpoints.sort(key=lambda x: x[1])
        breakpoint = breakpoints[0][0]

    # Combine the two alignments
    combined_aln = combine_chimeric_alignment(repr, supp, readlength, breakpoint)

    # Obtain aligned positions of the two ends of each alignments
    repr_pos5 = combined_aln[0][1] + 1
    repr_pos3 = next(x[1] for x in combined_aln if x[0] == breakpoint-2) + 1
    supp_pos5 = next(x[1] for x in combined_aln if x[0] == breakpoint-1) + 1
    supp_pos3 = combined_aln[-1][1] + 1

    return repr_pos5, repr_pos3, supp_pos5, supp_pos3

#############################################
def get_fragments_chimeric(repr, supp, tabix, ligseq, readlength, cutlength,
                         return_positions=False):
    """
    Given a pair of pysam.AlignedSegment objects forming a chimeric alignment,
    along with a pysam.TabixFile of restriction fragments, determine the
    set of fragments that intersect with the two alignments. 
    """
    # Obtain 5' and 3' positions of both alignments
    try:
        repr_pos5, repr_pos3, supp_pos5, supp_pos3 =\
            get_positions_chimeric(repr, supp, ligseq, readlength)
    except BadBreakpointException:
        raise

    # Re-sort the aligned positions and get chromosome names
    repr_posleft, repr_posright = sorted([ repr_pos5, repr_pos3 ])
    supp_posleft, supp_posright = sorted([ supp_pos5, supp_pos3 ])
    repr_chrom = repr.reference_name
    supp_chrom = supp.reference_name

    # Correct positions by restriction enzyme cutter-length, depending on
    # alignment orientation
    repr_posleft += cutlength
    repr_posright -= cutlength
    supp_posleft += cutlength
    supp_posright -= cutlength

    # Get restriction fragments
    repr_pos5_1 = repr_pos5 - cutlength if repr.is_reverse else repr_pos5 + cutlength
    supp_pos5_1 = supp_pos5 - cutlength if supp.is_reverse else supp_pos5 + cutlength
    frags_repr = get_fragments(repr_chrom, repr_posleft, repr_posright, tabix)
    frag5_repr = get_fragments(repr_chrom, repr_pos5_1, repr_pos5_1+1, tabix)
    frags_supp = get_fragments(supp_chrom, supp_posleft, supp_posright, tabix)
    frag5_supp = get_fragments(supp_chrom, supp_pos5_1, supp_pos5_1+1, tabix)
    #print(frags_repr, frag5_repr, frags_supp, frag5_supp)

    try:
        (frag5_repr,) = frag5_repr
        (frag5_supp,) = frag5_supp
    except:
        print(frags_repr, frag5_repr, frags_supp, frag5_supp)
        raise

    if return_positions:
        return frags_repr, frag5_repr, frags_supp, frag5_supp,\
               repr_pos5, repr_pos3, supp_pos5, supp_pos3

    return frags_repr, frag5_repr, frags_supp, frag5_supp

#############################################
class HiCPair(object):
    """
    Class definitions for a HiCPair object.
    """
    def __init__(self, reads=None, chrs='all', ligseq=None, cutlength=None,
                 mapq_threshold=30, sep_threshold=1000,
                 res_site_mindist=None, res_site_maxdist=750,
                 min_insert_size=0, max_insert_size=1000,
                 dangling_ends_threshold=5, fragments_file=None):
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
        self.pair_id = None
        self.reads1 = []
        self.reads2 = []
        self.pos1 = None
        self.pos2 = None
        self.chrom1 = None
        self.chrom2 = None
        self.strand1 = None
        self.strand2 = None
        self.chimeric = None
        self.readlength = None
        self.nearest_site1 = None
        self.nearest_site2 = None
        self.nearest_site_abs1 = None
        self.nearest_site_abs2 = None
        self.insert_size = None
        self.ligseq = ligseq
        self.cutlength = cutlength
        self.separations = None
        self.primary_sep = None
        self.supp_other_frag_sep = None

        # Initialize all flags to False, with exception of is_empty
        self.is_empty = True
        self.is_unmapped = False
        self.low_MAPQ = False
        self.other_contig = False
        self.three_distal_loci = False
        self.close_pair = False
        self.many_alignments = False
        self.chimeric_three_fragments = False
        self.chimeric_supp_other_diff_chroms = False
        self.chimeric_three_chromosomes = False
        self.contiguous_undigested_site = False
        self.chimeric_undigested_site = False
        self.chimeric_bad_breakpoint = False
        self.chimeric_wrong_orientation = False
        self.chimeric_disordered = False
        self.same_fragment = False
        self.dangling_ends = False
        self.far_nearest_restriction_site = False
        self.wrong_insert_size = False

        # Open tabix file of restriction fragments, if specified
        self.fragments = None
        if self.ligseq is not None and fragments_file is None:
            raise InvalidInputException('Fragments file should be specified for regular Hi-C data')
        elif self.ligseq is None and fragments_file is not None:
            raise InvalidInputException('Fragments file specified without restriction enzyme')
        elif fragments_file is not None:
            self.fragments = fragments_file

        # Initialize list of standard chromosomes/contigs
        if chrs == 'all':
            chrs = [ 'chr%d' % i for i in range(1,23) ] + ['chrX', 'chrY']

        # Add alignments if specified
        if reads is not None:
            self.add_reads(reads, mapq_threshold)

        # If nonempty, test for mappability
        if not self.is_empty:
            self.flag_unmapped_pair()

            # If nonempty and unmapped, test for low-MAPQ and/or mapping to other contigs
            if not self.is_unmapped:
                self.flag_low_MAPQ_pair(mapq_threshold=mapq_threshold)
                self.flag_other_contig(chrs)

                # If reads mapped with high MAPQ to a standard contig
                if not self.low_MAPQ and not self.other_contig:
                    # Determine if either of the reads aligned chimerically
                    if len(self.reads1) == 1 and len(self.reads2) == 1:
                        self.chimeric = 0
                    elif len(self.reads1) == 2 and len(self.reads2) == 1:
                        self.chimeric = 1
                    elif len(self.reads1) == 1 and len(self.reads2) == 2:
                        self.chimeric = 2
                    else:
                        self.many_alignments = True
                        try:
                            self.update_BAM_tags()
                        except EmptyPairException:
                            raise
                        return

                    # Assign restriction fragments to alignments
                    if self.ligseq is not None:
                        try:
                            self.assign_fragments(res_site_maxdist=res_site_maxdist,
                                                  min_insert_size=min_insert_size,
                                                  max_insert_size=max_insert_size,
                                                  dangling_ends_threshold=dangling_ends_threshold)
                        except BadBreakpointException:
                            self.chimeric_bad_breakpoint = True
                            try:
                                self.update_BAM_tags()
                            except EmptyPairException:
                                raise
                            return
                    # If enzyme was not specified, compute positions and separation
                    elif self.chimeric == 0:
                        self.pos1, _ = get_positions_contiguous(self.reads1[0])
                        self.pos2, _ = get_positions_contiguous(self.reads2[0])
                        self.separations = [ abs(self.pos1 - self.pos2) ]
                    elif self.chimeric == 1:
                        self.pos1, _, supp_pos5, _ = get_positions_chimeric(self.reads1[0], self.reads1[1])
                        self.pos2, _ = get_positions_contiguous(self.reads2[0])
                        self.separations = [ abs(self.pos1 - self.pos2),
                                             abs(self.pos1 - supp_pos5),
                                             abs(self.pos2 - supp_pos5) ]
                    else:
                        self.pos2, _, supp_pos5, _ = get_positions_chimeric(self.reads2[0], self.reads2[1])
                        self.pos1, _ = get_positions_contiguous(self.reads1[0])
                        self.separations = [ abs(self.pos1 - self.pos2),
                                             abs(self.pos1 - supp_pos5),
                                             abs(self.pos2 - supp_pos5) ]

                    # Obtain alignment orientations, depending on read positions
                    if self.pos1 < self.pos2:
                        self.strand1 = self.reads1[0].is_reverse
                        self.strand2 = self.reads2[0].is_reverse
                    else:
                        self.strand1 = self.reads2[0].is_reverse
                        self.strand2 = self.reads1[0].is_reverse

                    self.primary_sep = abs(self.pos1 - self.pos2)

                    # Flag pair if there are 3 alignments and the minimum distance b/t
                    # each pair exceeds sep_threshold
                    if self.chimeric > 0:
                        self.flag_three_distal_loci(sep_threshold=sep_threshold)

                    # Flag pair if separation is less than empirical threshold
                    # (<10kb for 6-cutter, <1kb for 4-cutter, respectively)
                    self.flag_close_pair()
                   
        # Update BAM tags
        try:
            self.update_BAM_tags()
        except EmptyPairException:
            pass

    #############################################
    def add_reads(self, reads, mapq_threshold):
        """
        Add alignment data (lists of pysam.AlignedSegment objects).
        """
        reads1 = sort_reads([ r for r in reads if r.is_read1 ],
                           mapq_threshold=mapq_threshold)
        reads2 = sort_reads([ r for r in reads if r.is_read2 ],
                           mapq_threshold=mapq_threshold)
        if reads1 is not None and reads2 is not None:
            self.reads1 = reads1
            self.reads2 = reads2
            self.readlength = len(self.reads1[0].query_sequence)
            self.is_empty = False

            # Set self.pair_id if not yet set
            if self.pair_id is None:
                self.pair_id = self.reads1[0].query_name

    #############################################
    def clear_reads(self):
        """
        Delete alignment data from HiCPair object and reset all flags and metadata.
        """
        self.pair_id = None
        self.reads1 = None
        self.reads2 = None
        self.pos1 = None
        self.pos2 = None
        self.chrom1 = None
        self.chrom2 = None
        self.strand1 = None
        self.strand2 = None
        self.chimeric = None
        self.readlength = None
        self.nearest_site1 = None
        self.nearest_site2 = None
        self.nearest_site_abs1 = None
        self.nearest_site_abs2 = None
        self.insert_size = None
        self.separations = None
        self.primary_sep = None
        self.supp_other_frag_sep = None
        self.is_empty = True
        self.is_unmapped = False
        self.other_contig = False
        self.three_distal_loci = False
        self.close_pair = False
        self.many_alignments = False
        self.chimeric_three_fragments = False
        self.chimeric_supp_other_diff_chroms = False
        self.chimeric_three_chromosomes = False
        self.contiguous_undigested_site = False
        self.chimeric_undigested_site = False
        self.chimeric_bad_breakpoint = False
        self.chimeric_wrong_orientation = False
        self.chimeric_disordered = False
        self.same_fragment = False
        self.dangling_ends = False
        self.far_nearest_restriction_site = False
        self.wrong_insert_size = False

    #############################################
    def update_BAM_tags(self):
        """
        Update tags for each alignment record in self.reads1 and self.reads2.

        COMPLETE LIST OF TAGS
        ---------------------
        BOOLEAN TAGS : YM - whether either read was unmapped, mapped with MAPQ less than a
                            given threshold, or mapped to a non-standard contig
                       YL - whether the reads mapped to three or more distal loci
                       YC - whether the reads represent a close-range contact (< 1kb)
                       YF - whether the reads mapped to three restriction fragments
                       YG - whether the reads mapped to four or more restriction fragments
                       YU - whether the reads contain an undigested restriction site
                       YS - whether the reads mapped to the same restriction fragment
                       YD - whether the reads contain a restriction site close to
                            (within 111 bp of) the 5' end of any read
                       YI - whether the reads come from an insert with anomalous length
                       YR - whether either read mapped far (> 750 bp) from the nearest restriction site
        INTEGER TAGS : Z1 - position of read 1
                       Z2 - position of read 2
                       ZM - position of nearest restriction site downstream of read 1
                       ZN - position of nearest restriction site downstream of read 2
        """
        if self.is_empty:
            raise EmptyPairException('Cannot write empty HiCPair to file')
        for r in self.reads1 + self.reads2:
            r.set_tag('YM', int(self.is_unmapped or self.low_MAPQ or self.other_contig), value_type='i')
            r.set_tag('YL', int(self.three_distal_loci), value_type='i')
            r.set_tag('YC', int(self.close_pair), value_type='i')
            r.set_tag('YA', int(self.many_alignments), value_type='i')
            r.set_tag('YZ', int(self.chimeric_disordered), value_type='i')
            r.set_tag('YF', int(self.chimeric_three_fragments), value_type='i')
            r.set_tag('YH', int(self.chimeric_supp_other_diff_chroms), value_type='i')
            r.set_tag('YT', int(self.chimeric_three_chromosomes), value_type='i')
            r.set_tag('YB', int(self.chimeric_bad_breakpoint), value_type='i')
            r.set_tag('YU', int(self.contiguous_undigested_site or
                                self.chimeric_undigested_site), value_type='i')
            r.set_tag('YO', int(self.chimeric_wrong_orientation), value_type='i')
            r.set_tag('YS', int(self.same_fragment), value_type='i')
            r.set_tag('YD', int(self.dangling_ends), value_type='i')
            r.set_tag('YI', int(self.wrong_insert_size), value_type='i')
            r.set_tag('YR', int(self.far_nearest_restriction_site), value_type='i')
            r.set_tag('Z1', self.pos1, value_type='i')
            r.set_tag('Z2', self.pos2, value_type='i')
            r.set_tag('ZM', self.nearest_site1, value_type='i')
            r.set_tag('ZN', self.nearest_site2, value_type='i')

    #############################################
    def write_to_file(self, outbam):
        """
        Write alignment records to a BAM file (pysam.AlignmentFile).
        """
        for r in self.reads1 + self.reads2:
            outbam.write(r)

    #############################################
    def flag_unmapped_pair(self):
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
                repr, supp = switch_chimeric_alignments(self.reads1[0], self.reads1[1])
                self.reads1 = [repr, supp] + self.reads1[2:]
                unmapped1 = False
            # Do the same for the second read
            if unmapped2 and len(self.reads2) > 1:
                repr, supp = switch_chimeric_alignments(self.reads2[0], self.reads2[1])
                self.reads2 = [repr, supp] + self.reads2[2:]
                unmapped2 = False
            self.is_unmapped = (unmapped1 or unmapped2)

    #############################################
    def flag_low_MAPQ_pair(self, mapq_threshold=30):
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
                repr, supp = switch_chimeric_alignments(self.reads1[0], self.reads1[1])
                self.reads1 = [repr, supp] + self.reads1[2:]
                low_mapq1 = False
            # Do the same for second read
            if low_mapq2 and len(self.reads2) > 1:
                repr, supp = switch_chimeric_alignments(self.reads2[0], self.reads2[1])
                self.reads2 = [repr, supp] + self.reads2[2:]
                low_mapq2 = False
            self.low_MAPQ = (low_mapq1 or low_mapq2)

    #############################################
    def flag_other_contig(self, chrs):
        """
        Sets self.other_contig to True if any one of the alignments mapped
        to a contig not in <chrs>.
        """
        if self.is_empty or self.is_unmapped:
            self.other_contig = False
        else:
            self.chrom1 = self.reads1[0].reference_name
            self.chrom2 = self.reads2[0].reference_name
            self.other_contig = (sum([ r.reference_name not in chrs for r in self.reads1 ]) > 0 or
                                 sum([ r.reference_name not in chrs for r in self.reads2 ]) > 0)

    #############################################
    def flag_three_distal_loci(self, sep_threshold=1000):
        """
        Sets self.three_distal_loci to True if the pair has three alignments and
        the minimum distance b/t each pair of alignments exceeds <sep_threshold>.
        """
        # Default to False if read pair is empty or unmapped
        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.three_distal_loci = False
        # If there are 2 or >3 alignments, set to False
        elif self.chimeric is None or self.chimeric == 0:
            self.three_distal_loci = False
        # If there are three alignments, look for minimum separation
        # between each pair of alignments 
        else:
            min_dist = min(self.separations)
            self.three_distal_loci = (min_dist > sep_threshold)

    #############################################
    def flag_close_pair(self):
        """
        Sets self.close_pair to True if abs(self.pos1 - self.pos2) is less than
        or equal to a threshold based on the cutter length.
        """
        thresholds = { 4 : 1000, 6 : 10000 }

        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.close_pair = False
        elif self.cutlength is None:
            self.close_pair = False
        elif self.pos1 is None or self.pos2 is None:
            self.close_pair = False
        else:
            self.close_pair = (self.primary_sep <= thresholds[self.cutlength])

    #############################################
    def assign_fragments(self, res_site_maxdist=750,
                         min_insert_size=0, max_insert_size=1000,
                         dangling_ends_threshold=5):
        """
        Determines all restriction fragments that intersect with each
        alignment. Sets all fragment-related flags as necessary.
        """
        ###################
        # BOTH CONTIGUOUS #
        ###################
        def assign_fragments_both_contiguous():
            """
            Assign fragments in the case where both reads mapped contiguously.
            """
            # Initialize flags
            self.contiguous_undigested_site = False
            self.same_fragment = False
            self.chimeric_undigested_site = False
            self.chimeric_three_fragments = False
            self.chimeric_three_chromosomes = False
            self.chimeric_wrong_orientation = False
            self.chimeric_disordered = False

            # Assign fragments to both alignments
            frags1, _, self.pos1, _ = get_fragments_contiguous(self.reads1[0], self.fragments,
                                                                    self.cutlength, return_positions=True)
            frags2, _, self.pos2, _ = get_fragments_contiguous(self.reads2[0], self.fragments,
                                                                    self.cutlength, return_positions=True)
            # If either read intersects with > 1 fragments, flag for undigested site
            if len(frags1) > 1 or len(frags2) > 1:
                self.contiguous_undigested_site = True
            # If the two reads map to the same fragment (at least partially), flag for same-fragment 
            if len(frags1.intersection(frags2)) > 0:
                self.same_fragment = True
            return frags1, frags2

        ################
        # ONE CHIMERIC #
        ################
        def assign_fragments_one_chimeric(repr, supp, other):
            """
            Assign fragments in the case where one read mapped chimerically.
            """
            # Initialize flags. NOTE: The three alignments cannot map to the same fragment
            # because (at this stage) the representative and supplementary must map to 
            # distinct fragments
            self.contiguous_undigested_site = False
            self.same_fragment = False
            self.chimeric_undigested_site = False
            self.chimeric_three_fragments = False
            self.chimeric_supp_other_diff_chroms = False
            self.chimeric_three_chromosomes = False
            self.chimeric_wrong_orientation = False
            self.chimeric_disordered = False

            # Assign fragments to all three alignments
            try:
                frags_repr, _, frags_supp, frag5_supp, pos_repr, _, pos_supp, _ =\
                        get_fragments_chimeric(repr, supp, self.fragments, self.ligseq,
                                               self.readlength, self.cutlength, return_positions=True)
            except BadBreakpointException:
                raise

            frags_other, frag5_other, pos_other, _ =\
                    get_fragments_contiguous(other, self.fragments,
                                             self.cutlength, return_positions=True)
            self.supp_other_frag_sep = frag5_supp - frag5_other
            # If supplementary and contiguous alignments map to same fragment
            if len(frags_supp.intersection(frags_other)) > 0:
                # Require that supplementary and contiguous are inner-oriented
                # Otherwise, flag as having wrong orientation
                if not ((pos_supp < pos_other and not supp.is_reverse and other.is_reverse) or\
                        (pos_other < pos_supp and not other.is_reverse and supp.is_reverse)):
                    self.chimeric_wrong_orientation = True
            else:
                # Now test if representative and contiguous alignments map to same fragment
                if len(frags_repr.intersection(frags_other)) > 0:
                    self.chimeric_disordered = True
                # Otherwise, all alignments fall into distinct fragments. NOTE: The representative
                # and supplementary cannot fall into the same fragment, because they must be
                # (at this stage) separated by a ligation site
                else:
                    self.chimeric_three_fragments = True
                    # Check whether the supplementary and contiguous are on separate chromosomes
                    if supp.reference_name != other.reference_name:
                        self.chimeric_supp_other_diff_chroms = True
                        if (repr.reference_name != other.reference_name and
                                repr.reference_name != supp.reference_name):
                            self.chimeric_three_chromosomes = True
                    # Require that supplementary and contiguous are inner-oriented
                    # Otherwise, flag as having wrong orientation
                    if not ((pos_supp < pos_other and not supp.is_reverse and other.is_reverse) or
                            (pos_other < pos_supp and not other.is_reverse and supp.is_reverse)):
                        self.chimeric_wrong_orientation = True
            # Finally, if any of the alignments map to more than one fragment, mark
            # as having an undigested site
            if len(frags_repr) > 1 or len(frags_supp) > 1 or len(frags_other) > 1:
                self.chimeric_undigested_site = True

            return frags_repr, frags_supp, frags_other, pos_repr, pos_supp, pos_other

        #####################################
        # RESTRICTION SITES AND INSERT SIZE #
        #####################################
        def compute_nearest_sites_insert_size(read1, read2, frags1, frags2):
            """
            Given alignments for the first and second reads, and their assigned
            fragments, locate the nearest restriction site to each read, and compute
            the expected insert size.
            """
            # Nearest restriction site to a given read alignment is the furthest
            # restriction site downstream of the 5' end that serves as an endpoint
            # to a fragment intersecting the alignment
            try:
                self.nearest_site1 = get_nearest_restriction_site(read1, frags1)
            except NoNearestSiteException:
                raise
            try:
                self.nearest_site2 = get_nearest_restriction_site(read2, frags2)
            except NoNearestSiteException:
                raise
            neardist1 = abs(self.nearest_site1 - self.pos1)
            neardist2 = abs(self.nearest_site2 - self.pos2)
            if neardist1 > res_site_maxdist or neardist2 > res_site_maxdist:
                self.far_nearest_restriction_site = True
            else:
                self.far_nearest_restriction_site = False

            # Compute insert size
            if read1.reference_name != read2.reference_name:
                self.insert_size = neardist1 + neardist2
            else:
                d = abs(self.pos1 - self.pos2)
                self.insert_size = min([neardist1 + neardist2, d])
            if self.insert_size < min_insert_size or self.insert_size > max_insert_size:
                self.wrong_insert_size = True
            else:
                self.wrong_insert_size = False

        #################
        # DANGLING ENDS #
        #################
        def is_dangling_ends(frags1, frags2):
            """
            Determines if the 5' end of a read is very close (< 5 bp by default)
            to a restriction site. Should only be done for contiguously aligned reads
            or representatives in a chimeric alignment.
            """
            # Initialize distance to nearest site to be arbitrarily large
            self.nearest_site_abs1 = self.nearest_site_abs2 = 1e9
            self.dangling_ends = False
            for f in frags1:
                s = min(abs(f.site1 - self.pos1), abs(f.site2 - self.pos1))
                if s < self.nearest_site_abs1:
                    self.nearest_site_abs1 = s
            for f in frags2:
                s = min(abs(f.site1 - self.pos2), abs(f.site2 - self.pos2))
                if s < self.nearest_site_abs2:
                    self.nearest_site_abs2 = s
            if (self.nearest_site_abs1 < dangling_ends_threshold
                    or self.nearest_site_abs2 < dangling_ends_threshold):
                self.dangling_ends = True

        #####################################
        # If pair is empty, unmapped, low-MAPQ, or maps to non-standard contig,
        # don't bother with fragment assignment
        if self.is_empty or self.is_unmapped or self.low_MAPQ or self.other_contig:
            self.chimeric_three_fragments = False
            self.chimeric_supp_other_diff_chroms = False
            self.chimeric_three_chromosomes = False
            self.chimeric_undigested_site = False
            self.chimeric_wrong_orientation = False
            self.chimeric_disordered = False
            self.contiguous_undigested_site = False
            self.same_fragment = False
            self.dangling_ends = False
            self.far_nearest_restriction_site = False
            self.wrong_insert_size = False
            self.supp_other_frag_sep = None

        # If both reads mapped contiguously, then reads cannot intersect with
        # more than 2 distal fragments, save for the possibility of undigested
        # restriction sites
        elif self.chimeric == 0:
            frags1, frags2 = assign_fragments_both_contiguous()
            self.supp_other_frag_sep = None
            self.separations = [ abs(self.pos1 - self.pos2) ]
            compute_nearest_sites_insert_size(self.reads1[0], self.reads2[0], frags1, frags2)
            is_dangling_ends(frags1, frags2)

        # If the first read was chimerically aligned
        elif self.chimeric == 1:
            try:
                frags_repr, frags_supp, frags_other, self.pos1, supp_pos5, self.pos2 =\
                    assign_fragments_one_chimeric(self.reads1[0], self.reads1[1], self.reads2[0])
            except BadBreakpointException:
                raise
            self.separations = [ abs(self.pos1 - self.pos2),
                                 abs(self.pos1 - supp_pos5),
                                 abs(self.pos2 - supp_pos5) ]
            # Here, one must be careful with computing nearest restriction sites: the
            # contiguous alignment's 5' end should be compared with the fragments intersecting
            # the supplementary alignment
            if not (self.chimeric_supp_other_diff_chroms or self.chimeric_disordered or
                    self.chimeric_wrong_orientation):
                compute_nearest_sites_insert_size(self.reads1[0], self.reads2[0], frags_repr, frags_supp)
            is_dangling_ends(frags_repr, frags_other)

        # If the second read was chimerically aligned
        elif self.chimeric == 2:
            try:
                frags_repr, frags_supp, frags_other, self.pos2, supp_pos5, self.pos1 =\
                    assign_fragments_one_chimeric(self.reads2[0], self.reads2[1], self.reads1[0])
            except BadBreakpointException:
                raise
            self.separations = [ abs(self.pos1 - self.pos2),
                                 abs(self.pos1 - supp_pos5),
                                 abs(self.pos2 - supp_pos5) ]
            if not (self.chimeric_supp_other_diff_chroms or self.chimeric_disordered or
                    self.chimeric_wrong_orientation):
                compute_nearest_sites_insert_size(self.reads1[0], self.reads2[0], frags_supp, frags_repr)
            is_dangling_ends(frags_other, frags_repr)

        # Otherwise, one or both reads map to multiple loci
        else:
            self.many_alignments = True
            self.chimeric_three_fragments = False
            self.chimeric_supp_other_diff_chroms = False
            self.chimeric_three_chromosomes = False
            self.chimeric_wrong_orientation = False
            self.contiguous_undigested_site = False
            self.chimeric_undigested_site = False
            self.chimeric_disordered = False
            self.same_fragment = False
            self.dangling_ends = False
            self.far_nearest_restriction_site = False
            self.wrong_insert_size = False
    
#############################################
def bsearch(data, value, return_value=False):
    """
    Given a sorted list, return the index of the greatest element
    that is less than the given value, via binary search.
    """
    low = 0
    high = len(data) - 1
    if data[low] > value or data[high] <= value:
        return None

    while low <= high:
        mid = int((low + high) / 2)
        if mid == len(data) - 1:
            return None
        elif data[mid] <= value < data[mid+1]:
            return (mid,data[mid]) if return_value else mid
        elif data[mid] > value:
            high = mid - 1
        else:
            low = mid + 1

    # Check that value being returned is less than the given value
    return (mid,data[mid]) if return_value else mid

#############################################
class HiCDataset(object):
    """
    Class definitions for a HiCDataset class, which stores summary information
    regarding a set of Hi-C read pairs.
    """
    def __init__(self, outdir, resolution=1000, mindist=2e4, maxdist=1e9,
                 close_mindist=10, close_maxdist=1e5, binratio=1.01,
                 minsitedist=0, maxsitedist=2000):
        """
        A HiCDataset object stores a number of dictionaries that summarize the
        properties of the dataset.
        """
        self.outdir = outdir           # Directory for output PDFs
        self.resolution = resolution   # Resolution at which to compute scalings

        # Bins for computing the power-law scaling of contact probability
        # as a function of distance
        self.distbins = list(np.arange(mindist, maxdist, resolution))
        self.distances = [ 0 for _ in self.distbins ]

        # Log-scale bins for measuring the ratio of alignment orientations
        # (FF,FR,RF,RR) over close-range contacts (< close_maxdist).
        nclosebins = math.ceil((np.log10(close_maxdist) - np.log10(close_mindist))
                               / np.log10(binratio))
        self.closebins = [ int(close_mindist * (binratio ** i))
                           for i in range(int(nclosebins)) ]
        self.orientations  = [ { 'FF':0, 'FR':0, 'RF':0, 'RR':0 } for _ in self.closebins ]

        # Strand ratios for different sets of pairs, each with a different minimum distance
        self.distance_ratios = { d : { 'FF':0, 'FR':0, 'RF':0, 'RR':0 }
                                 for d in [200, 500, 1000, 2000, 5000, 10000, 20000] }

        # Log-scale bins for measuring distance dependence of dangling-ends contacts;
        # bins range from close_mindist to maxdist
        nlogdistbins = math.ceil((np.log10(maxdist) - np.log10(close_mindist))
                             / np.log10(binratio))
        self.logdistbins = [ int(close_mindist * (binratio ** i))
                             for i in range(int(nlogdistbins)) ]
        self.logdistances = [ 0 for _ in self.logdistbins ]

        # Base-resolution bins for distances to nearest restriction sites
        # and insert sizes
        self.sitebins = range(minsitedist, maxsitedist+1)
        self.nearest_sites = [ 0 for _ in self.sitebins ]
        self.contiguous_nearest_sites = [ 0 for _ in self.sitebins ]
        self.chimeric_nearest_sites   = [ 0 for _ in self.sitebins ]
        self.insert_sizes  = [ 0 for _ in self.sitebins ]

        # Dictionary of tallies for various filters
        self.filters = { 'unmapped' : 0,
                         'low_MAPQ' : 0,
                         'other_contig' : 0,
                         'three_distal_loci' : 0,
                         'many_alignments' : 0,
                         'close_pair' : 0,
                         'contiguous_undigested_site' : 0,
                         'chimeric_undigested_site' : 0,
                         'chimeric_three_fragments' : 0,
                         'chimeric_supp_other_diff_chroms' : 0,
                         'chimeric_three_chromosomes' : 0,
                         'chimeric_bad_breakpoint' : 0,
                         'chimeric_wrong_orientation' : 0,
                         'chimeric_disordered' : 0,
                         'same_fragment' : 0,
                         'dangling_ends' : 0,
                         'far_nearest_restriction_site' : 0,
                         'wrong_insert_size' : 0 }
        self.close_filters = { 'close_pair' : [],
                               'far_nearest_restriction_site' : [],
                               'same_fragment' : [],
                               'wrong_insert_size' : [],
                               'dangling_ends' : [] }
        self.triple_filters = { 'three_distal_loci' : [],
                                'chimeric_three_fragments' : [] }

        # Histograms of contact distance for same-fragment or dangling-ends contacts
        self.same_frag_distance = [ 0 for _ in self.closebins ]
        self.same_frag_control  = [ 0 for _ in self.closebins ]
        self.dang_ends_distance = [ 0 for _ in self.logdistbins ]
        self.dang_ends_control  = [ 0 for _ in self.logdistbins ]

        # Histograms of contact distance vs. n.r.s and insert size
        self.nrs_vs_distance = { strands : { 'nearest_site' : [], 'distance' : [] }
                                 for strands in ['FF','FR','RF','RR'] }
        self.nrsabs_vs_distance = { strands : { 'nearest_site_abs' : [], 'distance' : [] }
                                    for strands in ['FF','FR','RF','RR'] }
        self.insert_vs_distance = { 'insert_size' : [], 'distance' : [] }

        # Histograms of minimal contact distances over all three-fragment contacts
        self.triple_min_distance = [ 0 for _ in self.closebins ]
        self.triple_frag_distance = [ 0 for _ in range(200) ]

    #############################################
    def add_pair(self, pair):
        """
        Adds a read pair to the Hi-C dataset.
        """
        # Run through filters that apply to the pair
        if pair.is_empty:
            raise EmptyPairException('Empty pair should not be added to dataset')

        if pair.is_unmapped:                     self.filters['unmapped'] += 1
        if pair.low_MAPQ:                        self.filters['low_MAPQ'] += 1
        if pair.other_contig:                    self.filters['other_contig'] += 1
        if pair.many_alignments:                 self.filters['many_alignments'] += 1
        if pair.contiguous_undigested_site:      self.filters['contiguous_undigested_site'] += 1
        if pair.chimeric_undigested_site:        self.filters['chimeric_undigested_site'] += 1
        if pair.chimeric_disordered:             self.filters['chimeric_disordered'] += 1
        if pair.chimeric_supp_other_diff_chroms: self.filters['chimeric_supp_other_diff_chroms'] += 1
        if pair.chimeric_three_chromosomes:      self.filters['chimeric_three_chromosomes'] += 1
        if pair.chimeric_wrong_orientation:      self.filters['chimeric_wrong_orientation'] += 1
        if pair.chimeric_bad_breakpoint:         self.filters['chimeric_bad_breakpoint'] += 1

        if pair.close_pair:
            self.filters['close_pair'] += 1
            self.close_filters['close_pair'].append(pair.pair_id)
        if pair.same_fragment:
            self.filters['same_fragment'] += 1
            self.close_filters['same_fragment'].append(pair.pair_id)
        if pair.dangling_ends:
            self.filters['dangling_ends'] += 1
            self.close_filters['dangling_ends'].append(pair.pair_id)
        if pair.far_nearest_restriction_site:
            self.filters['far_nearest_restriction_site'] += 1
            self.close_filters['far_nearest_restriction_site'].append(pair.pair_id)
        if pair.wrong_insert_size:
            self.filters['wrong_insert_size'] += 1
            self.close_filters['wrong_insert_size'].append(pair.pair_id)
        if pair.three_distal_loci:
            self.filters['three_distal_loci'] += 1
            self.triple_filters['three_distal_loci'].append(pair.pair_id)
        if pair.chimeric_three_fragments:
            self.filters['chimeric_three_fragments'] += 1
            self.triple_filters['chimeric_three_fragments'].append(pair.pair_id)

        # Add distance, orientation, restriction site, and insert size information
        # for pairs that pass essential filters (non-empty, uniquely mapped to
        # standard contigs, either contiguously or in a 3-way chimeric alignment)
        if not pair.is_empty and not pair.is_unmapped and not pair.low_MAPQ\
                             and not pair.other_contig and not pair.many_alignments\
                             and not pair.chimeric_disordered\
                             and not pair.chimeric_supp_other_diff_chroms\
                             and not pair.chimeric_three_chromosomes\
                             and not pair.chimeric_wrong_orientation\
                             and not pair.chimeric_bad_breakpoint:
            # Exclude all cases for which positions are undefined
            if pair.pos1 is None or pair.pos2 is None:
                return
            # Compute distances and corresponding histogram bins for
            # all intra-chromosomal contacts
            if pair.chrom1 == pair.chrom2:
                dist = abs(pair.pos1 - pair.pos2)
                distbin = bsearch(self.distbins, dist)
                closebin = bsearch(self.closebins, dist)
                logdistbin = bsearch(self.logdistbins, dist)
            else:   # For all inter-chromosomal contacts, set distance to be very large
                dist = np.inf
                distbin = closebin = logdistbin = None
            # Tabulate distances for intra-chromosomal contacts passing filters
            if distbin is not None and not pair.close_pair and not pair.same_fragment and\
                    not pair.far_nearest_restriction_site and not pair.wrong_insert_size and\
                    not pair.three_distal_loci:
                self.distances[distbin] += (1. / (chrom_lengths[pair.chrom1] - dist - 1))
            # Alignment orientation as a function of distance, for close pairs (< 100kb)
            if closebin is None or pair.strand1 is None or pair.strand2 is None:
                pass
            elif pair.strand1 and pair.strand2:
                self.orientations[closebin]['RR'] += 1.
            elif pair.strand1:
                self.orientations[closebin]['RF'] += 1.
            elif pair.strand2:
                self.orientations[closebin]['FR'] += 1.
            else:
                self.orientations[closebin]['FF'] += 1.
            # Alignment orientation for different sets of reads with min. distance cutoffs
            if 0 < dist < np.inf:
                for d in [200, 500, 1000, 2000, 5000, 10000, 20000]:
                    if dist > d:
                        if pair.strand1 and pair.strand2:
                            self.distance_ratios[d]['RR'] += 1.
                        elif pair.strand1:
                            self.distance_ratios[d]['RF'] += 1.
                        elif pair.strand2:
                            self.distance_ratios[d]['FR'] += 1.
                        else:
                            self.distance_ratios[d]['FF'] += 1.
            # Tabulate contact distance for same-fragment pairs
            if closebin is not None:
                self.same_frag_control[closebin] += 1
                if pair.same_fragment:
                    self.same_frag_distance[closebin] += 1
            # Tabulate contact distance for dangling-ends pairs 
            if pair.dangling_ends:
                if logdistbin is not None:
                    self.dang_ends_distance[logdistbin] += 1
            if logdistbin is not None:
                self.dang_ends_control[logdistbin] += 1
            # Distances to absolute nearest restriction sites
            if pair.nearest_site_abs1 is not None and 0 < dist < np.inf:
                if pair.strand1 and pair.strand2:
                    self.nrsabs_vs_distance['RR']['nearest_site_abs'].append(pair.nearest_site_abs1)
                    self.nrsabs_vs_distance['RR']['distance'].append(math.log10(dist))
                elif pair.strand1:
                    self.nrsabs_vs_distance['RF']['nearest_site_abs'].append(pair.nearest_site_abs1)
                    self.nrsabs_vs_distance['RF']['distance'].append(math.log10(dist))
                elif pair.strand2:
                    self.nrsabs_vs_distance['FR']['nearest_site_abs'].append(pair.nearest_site_abs1)
                    self.nrsabs_vs_distance['FR']['distance'].append(math.log10(dist))
                else:
                    self.nrsabs_vs_distance['FF']['nearest_site_abs'].append(pair.nearest_site_abs1)
                    self.nrsabs_vs_distance['FF']['distance'].append(math.log10(dist))
            if pair.nearest_site_abs2 is not None and 0 < dist < np.inf:
                if pair.strand1 and pair.strand2:
                    self.nrsabs_vs_distance['RR']['nearest_site_abs'].append(pair.nearest_site_abs2)
                    self.nrsabs_vs_distance['RR']['distance'].append(math.log10(dist))
                elif pair.strand1:
                    self.nrsabs_vs_distance['RF']['nearest_site_abs'].append(pair.nearest_site_abs2)
                    self.nrsabs_vs_distance['RF']['distance'].append(math.log10(dist))
                elif pair.strand2:
                    self.nrsabs_vs_distance['FR']['nearest_site_abs'].append(pair.nearest_site_abs2)
                    self.nrsabs_vs_distance['FR']['distance'].append(math.log10(dist))
                else:
                    self.nrsabs_vs_distance['FF']['nearest_site_abs'].append(pair.nearest_site_abs2)
                    self.nrsabs_vs_distance['FF']['distance'].append(math.log10(dist))
            # Distances for three-fragment pairs
            if pair.chimeric_three_fragments:
                minbin = bsearch(self.closebins, min(pair.separations))
                if minbin is not None:
                    self.triple_min_distance[minbin] += 1
            # Distances to nearest restriction sites and insert sizes
            if pair.nearest_site1 is not None and pair.nearest_site2 is not None:
                neardist1 = abs(pair.pos1 - pair.nearest_site1)
                neardist2 = abs(pair.pos2 - pair.nearest_site2)
                nearbin1 = neardist1 - self.sitebins[0] if neardist1 <= self.sitebins[-1] else None
                nearbin2 = neardist2 - self.sitebins[0] if neardist2 <= self.sitebins[-1] else None
                if nearbin1 is not None:
                    self.nearest_sites[nearbin1] += 1
                    if pair.chimeric == 0 or pair.chimeric == 2:
                        self.contiguous_nearest_sites[nearbin1] += 1
                    else:
                        self.chimeric_nearest_sites[nearbin1] += 1
                if nearbin2 is not None:
                    self.nearest_sites[nearbin2] += 1
                    if pair.chimeric == 0 or pair.chimeric == 1:
                        self.contiguous_nearest_sites[nearbin2] += 1
                    else:
                        self.chimeric_nearest_sites[nearbin2] += 1
                # Append data for distance to nearest restriction site vs. contact distance
                if 0 < dist < np.inf:
                    if pair.strand1 and pair.strand2:
                        self.nrs_vs_distance['RR']['nearest_site'].append(neardist1)
                        self.nrs_vs_distance['RR']['nearest_site'].append(neardist2)
                        self.nrs_vs_distance['RR']['distance'].append(math.log10(dist))
                        self.nrs_vs_distance['RR']['distance'].append(math.log10(dist))
                    elif pair.strand1:
                        self.nrs_vs_distance['RF']['nearest_site'].append(neardist1)
                        self.nrs_vs_distance['RF']['nearest_site'].append(neardist2)
                        self.nrs_vs_distance['RF']['distance'].append(math.log10(dist))
                        self.nrs_vs_distance['RF']['distance'].append(math.log10(dist))
                    elif pair.strand2:
                        self.nrs_vs_distance['FR']['nearest_site'].append(neardist1)
                        self.nrs_vs_distance['FR']['nearest_site'].append(neardist2)
                        self.nrs_vs_distance['FR']['distance'].append(math.log10(dist))
                        self.nrs_vs_distance['FR']['distance'].append(math.log10(dist))
                    else:
                        self.nrs_vs_distance['FF']['nearest_site'].append(neardist1)
                        self.nrs_vs_distance['FF']['nearest_site'].append(neardist2)
                        self.nrs_vs_distance['FF']['distance'].append(math.log10(dist))
                        self.nrs_vs_distance['FF']['distance'].append(math.log10(dist))
                # Do the same for insert size
                insertbin = pair.insert_size - self.sitebins[0]\
                            if pair.insert_size <= self.sitebins[-1] else None
                if insertbin is not None:
                    self.insert_sizes[insertbin] += 1
                if 0 < dist < np.inf:
                    self.insert_vs_distance['insert_size'].append(pair.insert_size)
                    self.insert_vs_distance['distance'].append(math.log10(dist))

    #############################################
    def plot_summary_stats(self):
        """
        Plot the summary data and save to a PDF file with (relative) path given by
        self.plotfile.
        """
        # Initialize axes
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,9))

        # Plot various statistics
        self.plot_contact_distance(axes[0,0])
        self.plot_orientations_vs_distance(axes[0,1])
        #self.plot_distance_strand_ratios(axes[1,0])
        self.plot_nearest_restriction_sites(axes[1,0])
        self.plot_insert_sizes(axes[1,1])

        # Title each plot with a letter
        axes[0,0].set_title('A', loc='left', fontsize=18, fontweight='bold')
        axes[0,1].set_title('B', loc='left', fontsize=18, fontweight='bold')
        axes[1,0].set_title('C', loc='left', fontsize=18, fontweight='bold')
        axes[1,1].set_title('D', loc='left', fontsize=18, fontweight='bold')

        # Save figure to file
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'summary.pdf'))

    #############################################
    def plot_close_filter_stats(self):
        """
        Plot statistics for "close" filters:
         - Distance distribution for same-fragment pairs
         - Distribution of distances to nearest restriction sites vs. contact distance
           for close-range pairs
         - Distribution of insert sizes vs. contact distance for close-range pairs
        """
        # Initialize axes
        fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(10,14))

        # Plot various statistics
        self.plot_close_filter_venn(axes[0,0], axes[1,0], axes[2,0], axes[3,0])
        self.plot_same_fragment_distances(axes[0,1])
        self.plot_dangling_ends_distances(axes[1,1])
        self.plot_nrs_vs_distance(axes[2,1])
        self.plot_insert_size_vs_distance(axes[3,1])

        # Title each plot with a letter
        axes[0,0].set_title('A', loc='left', fontsize=18, fontweight='bold')
        axes[0,1].set_title('B', loc='left', fontsize=18, fontweight='bold')
        axes[1,0].set_title('C', loc='left', fontsize=18, fontweight='bold')
        axes[1,1].set_title('D', loc='left', fontsize=18, fontweight='bold')
        axes[2,0].set_title('E', loc='left', fontsize=18, fontweight='bold')
        axes[2,1].set_title('F', loc='left', fontsize=18, fontweight='bold')
        axes[3,0].set_title('G', loc='left', fontsize=18, fontweight='bold')
        axes[3,1].set_title('H', loc='left', fontsize=18, fontweight='bold')

        # Save figure to file
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'close_filters.pdf'))

    #############################################
    def plot_nrs_vs_distance_stats(self):
        """
        Plot a separate figure with the distance-to-NRS vs. contact distance plot
        and the strand orientation ratios of the respective categories.
        """
        # Initialize axes
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,4))

        # Plot various statistics
        self.plot_nrs_vs_distance_strand_ratios(axes[0])
        self.plot_nrs_correlations(axes[1])

        # Title each plot with a letter
        axes[0].set_title('A', loc='left', fontsize=18, fontweight='bold')
        axes[1].set_title('B', loc='left', fontsize=18, fontweight='bold')

        # Save figure to file
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'nrs_strand_ratios.pdf'))

    #############################################
    def plot_nrsabs_vs_distance_stats(self):
        """
        Plot a separate figure with the distance-to-absolute-NRS vs. contact distance plot
        and the strand orientation ratios of the respective categories.
        """
        # Initialize axes
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,4))

        # Plot various statistics
        self.plot_nrsabs_vs_distance_strand_ratios(axes[1])
        self.plot_nrsabs_vs_distance(axes[0])

        # Title each plot with a letter
        axes[0].set_title('A', loc='left', fontsize=18, fontweight='bold')
        axes[1].set_title('B', loc='left', fontsize=18, fontweight='bold')

        # Save figure to file
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'nrsabs_strand_ratios.pdf'))

    #############################################
    def plot_close_filter_intersect(self):
        """
        Plot intersections of filters using pyUpSet.
        """
        filters_dict = { ''.join([ x[0] for x in k.split('_') ]).upper() :
                         pd.DataFrame(self.close_filters[k],[ 0 for _ in self.close_filters[k] ])
                         for k in self.close_filters }
        fig = plt.figure()
        pyu.plot(filters_dict, inters_size_bounds=(1,np.inf))
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'close_filters_intersect.pdf'))

    #############################################
    def plot_triple_filter_stats(self):
        """
        Plot statistics for "close" filters:
         - Distance distribution for same-fragment pairs
         - Distribution of distances to nearest restriction sites vs. contact distance
           for close-range pairs
         - Distribution of insert sizes vs. contact distance for close-range pairs
        """
        # Initialize axes
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,4))

        # Plot various statistics
        self.plot_triple_filter_venn(axes[0])
        self.plot_triple_min_distance(axes[1])

        # Title each plot with a letter
        axes[0].set_title('A', loc='left', fontsize=18, fontweight='bold')
        axes[1].set_title('B', loc='left', fontsize=18, fontweight='bold')

        # Save figure to file
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'triple_filters.pdf'))

    #############################################
    def plot_triple_filter_intersect(self):
        """
        Plot intersections of filters using pyUpSet.
        """
        filters_dict = { ''.join([ x[0] for x in k.split('_') ]).upper() :
                         pd.DataFrame(self.triple_filters[k], [ 0 for _ in self.triple_filters[k] ])
                         for k in self.triple_filters }
        fig = plt.figure()
        pyu.plot(filters_dict, inters_size_bounds=(1,np.inf))
        plt.tight_layout()
        plt.savefig(os.path.join(self.outdir, 'triple_filters_intersect.pdf'))

    #############################################
    def plot_contact_distance(self, ax):
        """
        Plot distribution of contact distances.
        """
        total_dist = sum(self.distances)
        ax.loglog(self.distbins, [ x / total_dist for x in self.distances ])
        ax.set_xlabel('Contact distance')
        ax.set_ylabel('Probability')

    #############################################
    def plot_distance_strand_ratios(self, ax):
        """
        Plot a stacked 100% horizontal bar plot showing the relative ratios
        of the four strand orientations for various distance cutoffs.
        """
        # Normalize the values in each dict
        print(self.distance_ratios)
        for d in self.distance_ratios:
            total = sum(self.distance_ratios[d].values())
            for strand in self.distance_ratios[d]:
                self.distance_ratios[d][strand] /= total

        # Plot bar chart
        distances = [200, 500, 1000, 2000, 5000, 10000, 20000]
        y_cat = np.arange(len(distances))
        totals = { d : 1. for d in distances }
        for strand, strand_label, color in zip(['RR', 'RF', 'FR', 'FF'],
                                               ['Left', 'Outer', 'Inner', 'Right'],
                                               ["#8172B2", "#C44E52", "#55A868", "#4C72B0"]):
            ax.barh(y_cat, [ totals[d] for d in distances ],
                    align='center', color=color, label=strand_label,
                    edgecolor=color, height=0.4)
            for d in distances:
                totals[d] -= self.distance_ratios[d][strand]

        ax.set_xlabel('Fraction of read pairs')
        ax.set_yticks(y_cat)
        ax.set_yticklabels([ '> %d' % d for d in distances ])
        ax.legend(bbox_to_anchor=(0., 0.92, 1., .102), loc=3,
                  ncol=4, mode="expand", borderaxespad=0.)

    #############################################
    def plot_orientations_vs_distance(self, ax):
        """
        Plot distribution of alignment orientations (FF,FR,RF,RR) as
        a function of contact distance.
        """ 
        orientations_FF = [ self.orientations[i]['FF'] for i in range(len(self.closebins)) ]
        orientations_FR = [ self.orientations[i]['FR'] for i in range(len(self.closebins)) ]
        orientations_RF = [ self.orientations[i]['RF'] for i in range(len(self.closebins)) ]
        orientations_RR = [ self.orientations[i]['RR'] for i in range(len(self.closebins)) ]

        ax.semilogx(self.closebins, np.log10(orientations_FF), label='Right')
        ax.semilogx(self.closebins, np.log10(orientations_FR), label='Inner')
        ax.semilogx(self.closebins, np.log10(orientations_RF), label='Outer')
        ax.semilogx(self.closebins, np.log10(orientations_RR), label='Left')
        ax.legend(loc='upper right')
        ax.set_xlabel('Contact distance')
        ax.set_ylabel('Log frequency')

    #############################################
    def plot_nearest_restriction_sites(self, ax):
        """
        Plot distribution of distances to nearest restriction sites.
        """
        ax.plot(self.sitebins, self.nearest_sites, label='All')
        ax.plot(self.sitebins, self.contiguous_nearest_sites, alpha=0.7, label='Contig.')
        ax.plot(self.sitebins, self.chimeric_nearest_sites, alpha=0.7, label='Chim.')
        ax.set_xlim([0, 1000])
        ax.legend(loc='upper right')
        ax.set_xlabel('Distance to nearest restriction site')
        ax.set_ylabel('Frequency')

    #############################################
    def plot_insert_sizes(self, ax):
        """
        Plot distribution of insert sizes.
        """
        ax.plot(self.sitebins, self.insert_sizes)
        ax.set_xlabel('Insert size')
        ax.set_ylabel('Frequency')

    #############################################
    def plot_close_filter_venn(self, ax1, ax2, ax3, ax4):
        """
        Plot distribution of distances for same-fragment pairs.
        """
        same_fragment = set(self.close_filters['same_fragment'])
        dangling_ends = set(self.close_filters['dangling_ends'])
        far_from_nrs  = set(self.close_filters['far_nearest_restriction_site'])
        insert_size   = set(self.close_filters['wrong_insert_size'])
        close_pair    = set(self.close_filters['close_pair'])
        venn2([close_pair, same_fragment], ax=ax1, set_labels=('CP', 'SF'))
        venn2([close_pair, dangling_ends], ax=ax2, set_labels=('CP', 'DE'))
        venn2([close_pair, far_from_nrs], ax=ax3, set_labels=('CP', 'FNRS'))
        venn2([far_from_nrs, insert_size], ax=ax4, set_labels=('FNRS', 'WIS'))

    #############################################
    def plot_same_fragment_distances(self, ax):
        """
        Plot distribution of distances for same-fragment pairs.
        """
        ax.semilogx(self.closebins, self.same_frag_control, label='Ctl')
        ax.semilogx(self.closebins, self.same_frag_distance, label='SF')
        ax.legend(loc='upper right')
        ax.set_xlabel('Contact distance')
        ax.set_ylabel('Frequency')

    #############################################
    def plot_dangling_ends_distances(self, ax):
        """
        Plot distribution of distances for dangling-ends pairs.
        """
        ax.semilogx(self.logdistbins, self.dang_ends_control, label='Ctl')
        ax.semilogx(self.logdistbins, self.dang_ends_distance, label='DE')
        ax.legend(loc='upper right')
        ax.set_xlabel('Contact distance')
        ax.set_ylabel('Frequency')
    
    #############################################
    def plot_nrs_vs_distance(self, ax):
        """
        Plot a 2D histogram of contact distance vs. distance to nearest restriction site.
        """
        # Downsample if there are more than 40000 points, ensuring that NRSs
        # corresponding to a read pair are sampled together
        nrs_vs_distance_sample = { strands : { 'nearest_site' : [], 'distance' : [] }
                                        for strands in ['FF','FR','RF','RR'] }
        npoints = sum(len(self.nrs_vs_distance[strands]['nearest_site'])
                      for strands in ['FF','FR','RF','RR'])
        if npoints > 40000:
            # Assuming the four strand combinations are split evenly (each 25%)
            # downsample each set of points to 10000
            for strand in ['FF', 'FR', 'RF', 'RR']:
                idx = np.random.choice(range(0,len(self.nrs_vs_distance[strand]['nearest_site']),2),
                                       5000)
                for k in ['nearest_site', 'distance']:
                    for i in idx:
                        nrs_vs_distance_sample[strand][k].append(self.nrs_vs_distance[strand][k][i])
                        nrs_vs_distance_sample[strand][k].append(self.nrs_vs_distance[strand][k][i+1])
        else:
            nrs_vs_distance_sample = self.nrs_vs_distance

        # Plot each set of points
        for strand, strand_label, color in zip(['FF', 'FR', 'RF', 'RR'],
                                               ['Right', 'Inner', 'Outer', 'Left'],
                                               ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]):
            ax.scatter(nrs_vs_distance_sample[strand]['nearest_site'],
                       nrs_vs_distance_sample[strand]['distance'],
                       marker='.', color=color, linewidth=0, alpha=0.5, label=strand_label)
        ax.plot([750, 750], [0, 9], color='red')
        ax.plot([0, 3000],  [3, 3], color='red')
        ax.legend(loc='upper right')
        ax.set_xlim([0, 3000])
        ax.set_ylim([0, 9])
        ax.set_xlabel('Distance to nearest restriction site')
        ax.set_ylabel('Log contact distance')

    #############################################
    def plot_nrs_vs_distance_strand_ratios(self, ax):
        """
        Plot a stacked 100% horizontal bar plot showing the relative ratios
        of the four strand orientations in the four regions of the NRS vs. distance
        plot.
        """
        neither = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        CP      = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        FNRS    = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        both    = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        for strands in ['FF','FR','RF','RR']:
            for n, d in zip(self.nrs_vs_distance[strands]['nearest_site'],
                            self.nrs_vs_distance[strands]['distance']):
                if n > 750 and d < 3:         # If distance to NRS is large and contact distance is small
                    both[strands] += 1.
                elif n > 750:                 # If only the distance to NRS is large
                    FNRS[strands] += 1.
                elif d < 3:                   # If only the contact distance is small
                    CP[strands] += 1.
                else:                         # Otherwise
                    neither[strands] += 1.

        # Normalize the values in each dict
        for dct in [neither, CP, FNRS, both]:
            total = sum(dct.values())
            for k in dct:
                dct[k] /= total

        # Plot bar chart
        y_cat = np.arange(4)
        both_total = CP_total = FNRS_total = neither_total = 1.
        for strand, strand_label, color in zip(['RR', 'RF', 'FR', 'FF'],
                                               ['Left', 'Outer', 'Inner', 'Right'],
                                               ["#8172B2", "#C44E52", "#55A868", "#4C72B0"]):
            ax.barh(y_cat, [both_total, CP_total, FNRS_total, neither_total],
                    align='center', color=color, label=strand_label,
                    edgecolor=color, height=0.3)
            both_total -= both[strand]
            CP_total   -= CP[strand]
            FNRS_total -= FNRS[strand]
            neither_total -= neither[strand]

        ax.set_xlabel('Fraction of read pairs')
        ax.set_yticks(y_cat)
        ax.set_yticklabels(['Both', 'CP', 'FNRS', 'Neither'])
        ax.invert_yaxis()
        ax.legend(bbox_to_anchor=(0., 0.92, 1., .102), loc=3,
                  ncol=4, mode="expand", borderaxespad=0.)

    #############################################
    def plot_nrs_correlations(self, ax):
        """
        Plot the distances to NRS for all read pairs in a scatter plot.
        """
        # Use separate dictionaries to store close pairs vs. far pairs
        nrs1_vs_nrs2_close = { strands : { 'nearest_site1' : [], 'nearest_site2' : [] }
                               for strands in ['FF','FR','RF','RR'] }
        nrs1_vs_nrs2_far   = { strands : { 'nearest_site1' : [], 'nearest_site2' : [] }
                               for strands in ['FF','FR','RF','RR'] }
        # Downsample to 40000 points total if necessary
        npoints = sum(len(self.nrs_vs_distance[strands]['nearest_site'])
                      for strands in ['FF','FR','RF','RR'])
        for strand in ['FF', 'FR', 'RF', 'RR']:
            if npoints > 40000:
                idx = np.random.choice(range(0,len(self.nrs_vs_distance[strand]['nearest_site']),2),
                                       5000)
            else:
                idx = range(0,len(self.nrs_vs_distance[strand]['nearest_site']),2)
            for i in idx:
                if self.nrs_vs_distance[strand]['distance'][i] >= 3:
                    nrs1_vs_nrs2_far[strand]['nearest_site1'].append(
                            self.nrs_vs_distance[strand]['nearest_site'][i])
                    nrs1_vs_nrs2_far[strand]['nearest_site2'].append(
                            self.nrs_vs_distance[strand]['nearest_site'][i+1])
                else:
                    nrs1_vs_nrs2_close[strand]['nearest_site1'].append(
                            self.nrs_vs_distance[strand]['nearest_site'][i])
                    nrs1_vs_nrs2_close[strand]['nearest_site2'].append(
                            self.nrs_vs_distance[strand]['nearest_site'][i+1])

        # Plot each set of points
        for strand, strand_label, color in zip(['FF', 'FR', 'RF', 'RR'],
                                               ['Right', 'Inner', 'Outer', 'Left'],
                                               ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]):
            ax.scatter(nrs1_vs_nrs2_close[strand]['nearest_site1'],
                       nrs1_vs_nrs2_close[strand]['nearest_site2'],
                       marker='.', color=color, linewidth=0, alpha=0.5)
            ax.scatter(nrs1_vs_nrs2_far[strand]['nearest_site1'],
                       nrs1_vs_nrs2_far[strand]['nearest_site2'],
                       marker='x', color=color, linewidth=2, alpha=1., label=strand_label)
        ax.legend(loc='upper right')
        ax.set_xlim([0, 3000])
        ax.set_ylim([0, 3000])
        ax.set_xlabel('Distance to NRS, read 1')
        ax.set_ylabel('Distance to NRS, read 2')

    #############################################
    def plot_nrsabs_vs_distance(self, ax):
        """
        Plot a 2D histogram of contact distance vs. distance to absolute
        nearest restriction site.
        """
        # Downsample if there are more than 40000 points
        nrsabs_vs_distance_sample = { strands : { 'nearest_site_abs' : [], 'distance' : [] }
                                      for strands in ['FF','FR','RF','RR'] }

        npoints = sum(len(self.nrsabs_vs_distance[strands]['nearest_site_abs'])
                      for strands in ['FF','FR','RF','RR'])
        if npoints > 40000:
            # Assuming the four strand combinations are split evenly (each 25%)
            # downsample each set of points to 10000
            for strand in ['FF', 'FR', 'RF', 'RR']:
                idx = np.random.choice(range(0,len(self.nrsabs_vs_distance[strand]['nearest_site_abs']),2),
                                       5000)
                for k in ['nearest_site_abs', 'distance']:
                    for i in idx:
                        nrsabs_vs_distance_sample[strand][k].append(self.nrsabs_vs_distance[strand][k][i])
                        nrsabs_vs_distance_sample[strand][k].append(self.nrsabs_vs_distance[strand][k][i+1])
        else:
            nrsabs_vs_distance_sample = self.nrsabs_vs_distance

        # Plot each set of points
        for strand, strand_label, color in zip(['FF', 'FR', 'RF', 'RR'],
                                               ['Right', 'Inner', 'Outer', 'Left'],
                                               ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]):
            ax.scatter(nrsabs_vs_distance_sample[strand]['nearest_site_abs'],
                       nrsabs_vs_distance_sample[strand]['distance'],
                       marker='.', color=color, linewidth=0, alpha=0.5, label=strand_label)
        ax.legend(loc='upper right')
        ax.set_xlim([0, 3000])
        ax.set_ylim([0, 9])
        ax.set_xlabel('Distance to absolute NRS')
        ax.set_ylabel('Log contact distance')

    #############################################
    def plot_nrsabs_vs_distance_strand_ratios(self, ax):
        """
        Plot a stacked 100% horizontal bar plot showing the relative ratios
        of the four strand orientations in the four regions of the NRS vs. distance
        plot.
        """
        neither = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        close_pair = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        dang_ends = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        both = { 'FF' : 0, 'FR' : 0, 'RF' : 0, 'RR' : 0 }
        for strands in ['FF','FR','RF','RR']:
            for n, d in zip(self.nrsabs_vs_distance[strands]['nearest_site_abs'],
                            self.nrsabs_vs_distance[strands]['distance']):
                if n > 5 and d < 3:           # If distance to NRS is large and contact distance is small
                    both[strands] += 1.
                elif n > 5:                   # If only the distance to NRS is large
                    dang_ends[strands] += 1.
                elif d < 3:                   # If only the contact distance is small
                    close_pair[strands] += 1.
                else:                         # Otherwise
                    neither[strands] += 1.

        # Normalize the values in each dict
        for dct in [neither, close_pair, dang_ends, both]:
            total = sum(dct.values())
            for k in dct:
                dct[k] /= total

        # Plot bar chart
        y_cat = np.arange(4)
        both_total = close_pair_total = dang_ends_total = neither_total = 1.
        for strand, strand_label, color in zip(['RR', 'RF', 'FR', 'FF'],
                                               ['Left', 'Outer', 'Inner', 'Right'],
                                               ["#8172B2", "#C44E52", "#55A868", "#4C72B0"]):
            ax.barh(y_cat, [both_total,
                            close_pair_total,
                            dang_ends_total,
                            neither_total],
                    align='center', color=color, label=strand_label,
                    edgecolor=color, height=0.3)
            both_total -= both[strand]
            close_pair_total -= close_pair[strand]
            dang_ends_total -= dang_ends[strand]
            neither_total -= neither[strand]

        ax.set_xlabel('Fraction of read pairs')
        ax.set_yticks(y_cat)
        ax.set_yticklabels(['Both', 'CP', 'DE', 'Neither'])
        ax.invert_yaxis()
        ax.legend(bbox_to_anchor=(0., 0.92, 1., .102), loc=3,
                  ncol=4, mode="expand", borderaxespad=0.)

    #############################################
    def plot_insert_size_vs_distance(self, ax):
        """
        Plot a 2D histogram of contact distance vs. insert size.
        """
        # Downsample if there are more than 40000 points
        insert_vs_distance_sample = { 'insert_size' : [], 'distance' : [] }

        if len(self.insert_vs_distance['insert_size']) > 40000:
            idx = np.random.choice(range(len(self.insert_vs_distance['insert_size'])), 40000)
            for k in self.insert_vs_distance:
                insert_vs_distance_sample[k] = [ self.insert_vs_distance[k][i] for i in idx ]
        else:
            insert_vs_distance_sample = self.insert_vs_distance

        ax.scatter(insert_vs_distance_sample['insert_size'],
                   insert_vs_distance_sample['distance'],
                   marker='.', linewidth=0, alpha=0.5)
        ax.plot([1e3, 1e3], [0, 9], color='red')
        ax.plot([0, 6000],  [3, 3], color='red') 

        ax.set_xlim([0, 6000])
        ax.set_ylim([0, 9])
        ax.set_xlabel('Insert size')
        ax.set_ylabel('Log contact distance')

    #############################################
    def plot_triple_filter_venn(self, ax):
        """
        Plot distribution of distances for same-fragment pairs.
        """
        three_distal_loci = set(self.triple_filters['three_distal_loci'])
        three_fragments   = set(self.triple_filters['chimeric_three_fragments'])
        venn2([three_distal_loci, three_fragments], ax=ax, set_labels=('3DL', '3F'))

    #############################################
    def plot_triple_min_distance(self, ax):
        """
        Plot distribution of distances for three-fragment pairs.
        """
        ax.semilogx(self.closebins, self.triple_min_distance)
        ax.set_xlabel('Minimum contact distance')
        ax.set_ylabel('Frequency')

#############################################
####### MAIN SCRIPT FOR DATA ANALYSIS #######
#############################################
def parse_BAM(bam, outbam, tabix, enzyme='MboI', mapq_threshold=30,
              sep_threshold=1000, verbose=True):
    """
    Parse input BAM file, process each pair and flag with necessary filters,
    and re-write each pair to an output BAM file.

    Parameters
    ----------
    bam : string
        Path to input BAM file
    outbam : string
        Path to output BAM file, to which filter flags will be reported
    tabix : string
        Path to Tabix file of restriction fragments
    enzyme : string
        Name of restriction enzyme used to create dataset
    mapq_threshold : int
        MAPQ score threshold for distinguishing unique alignments
    sep_threshold : int
        Base-pair threshold for distinguishing alignments that are "close" together
    verbose : bool
        Whether or not to print verbose output to stdout
    """
    curr_id = None
    curr_reads = []
    finished = False
    npairs = 0
    
    # Initialize directory for output BAM file and figures
    outdir = outbam.split('.bam')[0] + '_figures/'
    outtxt = os.path.join(outdir, 'summary.txt')
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)

    dataset = HiCDataset(outdir, mindist=2e4, maxdist=1e9,
                         close_mindist=10, close_maxdist=1e5, binratio=1.01,
                         minsitedist=0, maxsitedist=2000, resolution=20000)
    fragments_file = pysam.TabixFile(tabix)
    ligseq = get_ligation_sequence(enzyme)           # Ligation sequence for enzyme
    cutlength = len(get_restriction_enzyme(enzyme))  # Cutter length of enzyme

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
                npairs += 1
                dataset.add_pair(pair)
                pair.write_to_file(outbf)
                # Update curr_id and re-initialize curr_reads
                curr_id = read.query_name
                curr_reads = []
                if npairs % 10000 == 0 and verbose:
                    print('Processed %d read pairs. Flag statistics:' % npairs)
                    for f in dataset.filters:
                        print('%s : %d' % (f, dataset.filters[f]))

            # Finally, gather alignment data from current read
            curr_reads.append(read)

        outbf.close()

    # Plot and print summary statistics
    print('Summarizing and plotting filter statistics ...')
    with open(outtxt, 'w') as f:
        f.write('Number of pairs\t%d\n' % npairs)
        for ft in dataset.filters:
            f.write('%s\t%d\n' % (ft, dataset.filters[ft]))
    dataset.plot_summary_stats()
    dataset.plot_close_filter_stats()
    dataset.plot_nrs_vs_distance_stats()
    dataset.plot_nrsabs_vs_distance_stats()
    dataset.plot_triple_filter_stats()

    # Plot intersections
    print('Plotting intersections of filters ...')
    dataset.plot_close_filter_intersect()
    dataset.plot_triple_filter_intersect()

    fragments_file.close()

#############################################
def parse_inputs():
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
    bam, outbam, tabix, verbose = parse_inputs()
    parse_BAM(bam, outbam, tabix, verbose=verbose)

#############################################
if __name__ == '__main__':
    main()
