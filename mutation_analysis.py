"""Script that analyzes occurrence of specific mutations in NA lineages.

For a sequence set, finds all occurrences of a specific mutation.

Then builds a phylogenetic tree (using codonphyml) of all sequences with the
mutation and a subset of those without.

Jesse Bloom, 2013."""


import os
import re
import random
import mutation_analysis_funcs


def RemoveOutliers(seqs, outlierfile):
    """Removes outlier sequences.

    *seqs* is a list of sequences as *(header, seq)* 2-tuples.

    *outlierfile* is a sequence that lists strain names to remove on separate lines.

    Returns a list that excludes any sequence in  *seqs* that has a header that
    contains a substring in *outlierfile*.
    """
    outliers = [head.strip() for head in open(outlierfile).readlines() if not head.isspace()]
    cleanseqs = []
    for (head, seq) in seqs:
        for outlier in outliers:
            if outlier in head:
                break
        else:
            cleanseqs.append((head, seq))
    return cleanseqs

def ClassifyByMutation(cds_alignment, prot_alignment, mutation):
    """Finds all sequences with a given protein mutation.

    Returns *(has_mut, no_mut, site_counts)*.

    *has_mut* lists cds sequences with mutation.

    *no_mut* lists cds sequences without mutation.

    *site_counts* keyed by amino acid identities and values are 
    counts of that identity at the site.
    """
    (wt, index, mut) = (mutation[0], int(mutation[1 : -1]) - 1, mutation[-1])
    has_mut = []
    no_mut = []
    site_counts = {}
    for (cds, prot) in zip(cds_alignment, prot_alignment):
        aa = prot[1][index]
        if aa in site_counts:
            site_counts[aa] += 1
        else:
            site_counts[aa] = 1
        if aa == mut:
            has_mut.append(cds)
        else:
            no_mut.append(cds)
    assert not has_mut or (site_counts[mut] == len(has_mut))
    return (has_mut, no_mut, site_counts)


def RemoveIdentical(seqs, identicalseqsfile):
    """Removes identical sequences / substrings.

    Goes through *seqs* and finds any sequences that are substrings or
    identical to another sequence. Returns a new list where all of these
    identical and substring sequences have been removed.

    Also writes the headers of identical sequence sets to *identicalseqsfile*.
    """
    seq_d = {}
    decorated_list = [(len(seq[1]), seq, seq[1].replace('-', '')) for seq in seqs]
    seqs.sort()
    seqs.reverse()
    cleanseqs = []
    for (n, (head, seq), nogapseq) in decorated_list:
        if nogapseq in seq_d:
            seq_d[nogapseq].append(head)
        else:
            substring = False
            for iseq in seq_d.iterkeys():
                if nogapseq in iseq:
                    substring = True
                    seq_d[iseq].append(head)
                    break
            if not substring:
                seq_d[nogapseq] = [head]
                cleanseqs.append((head, seq))
    f = open(identicalseqsfile, 'w')
    for heads in seq_d.itervalues():
        if len(heads) > 1:
            f.write("The following sequences are identical or substrings:\n%s\n\n" % ('\n'.join(heads)))
    f.close()
    return cleanseqs
        
    

def GetClosest(has_mut, no_mut, minclosest, printoutput=True):
    """For every sequence in *has_mut*, gets the closest sequence(s) in *no_mut*.

    *has_mut* and *no_mut* are lists of sequences as *(header, sequence)* tuples.

    For each sequence in *has_mut*, finds the single sequence (or set of sequences)
    with the highest fractional identity (excluding gaps) to the sequence in
    *has_mut*. Put all of these unique sequences into *closest*. The list
    *not_closest* is all sequences in *no_mut* that aren't in *closest*.
    *minclosest* is an integer specifying that we make sure that *closest*
    has at least this many entries.
    The returned variable is the tuple
    *(closest, not_closest)*

    If *printoutput* is *True*, prints to standard output some information.
    """
    closest = {}
    for (head, seq) in has_mut:
        decorated_list = []
        for (head2, seq2) in no_mut:
            ident = mutation_analysis_funcs.FracIdent(seq, seq2)
            decorated_list.append((ident, (head2, seq2)))
        decorated_list.sort()
        decorated_list.reverse()
        maxident = decorated_list[0][0]
        iclosest = [(head2, seq2) for (ident, (head2, seq2)) in decorated_list if ident >= maxident]
        ntoadd = minclosest - len(iclosest)
        if ntoadd > 0:
            iclosest += [(head2, seq2) for (ident, (head2, seq2)) in decorated_list[len(iclosest) : minclosest]]
        if printoutput:
            print "For the sequence with the mutation %s, the closest sequences without the mutation are:\n%s" % (head, '\n'.join(['\t%s' % head2 for (head2, seq2) in iclosest]))
        for tup in iclosest:
            closest[tup] = True
    not_closest = [tup for tup in no_mut if tup not in closest]
    closest = closest.keys()
    assert len(not_closest) + len(closest) == len(no_mut), "%d, %d, %d" % (len(not_closest), len(closest), len(no_mut))
    return (closest, not_closest)


def GetDate(head):
    """Returns date as float out of header. Returns *None* if date unknown."""
    date = head.split()[-2]
    if 'unknown' in date or 'NON' in date:
        return None
    else:
        return mutation_analysis_funcs.DateToOrdinal(date, 1968)


def GetNPerYear(seqs, nperyear):
    """Gets *nperyear* sequences per year for all seqs in *seqs*.
    """
    seqs_by_year = {}
    for (head, seq) in seqs:
        year = int(GetDate(head))
        if year not in seqs_by_year:
            seqs_by_year[year] = [(head, seq)]
        elif len(seqs_by_year[year]) < nperyear:
            seqs_by_year[year].append((head, seq))
    retained = []
    for yearseqs in seqs_by_year.itervalues():
        retained += yearseqs
    return retained


def main():
    """Main body of script."""
    # input / output variables
    # seqsets lists the sequence sets to analyze as tuples with these elements:
    #   * seqset : sequence set, should be name of sequence file with suffix .fasta
    #   * refseqhead : FASTA header of reference sequence for alignments in seqset file
    #   * mutation : mutation of interest in 1, 2, ... numbering
    #   * firstyear : only analyze seqs with year >= this
    #   * buildtree : build a tree, True or False
    #   * nperyear : maximum number per year of sequences without mutation of interest and not "closest" for sequence with mutation
    #   * minclosest: keep at least this many "closest" sequences for each with mutation
    seqsets = [
            ('pdmH1N1_NAs', 'cds:ACP41107 A/California/04/2009 H1N1 2009/04/01 NA', 'G147R', 2009, True, 10, 6),
            ('seasonalH1N1_NAs', 'cds:ADE28752 A/Brisbane/59/2007 H1N1 2007/07/01 NA', 'G147R', 2007, True, 6, 2),
            ('chickenH5N1_NAs', 'cds:BAM85823 A/chicken/Viet Nam/NCVD5/2003 H5N1 2003// NA', 'G147R', 2004, True, 10, 4),
            ]
    seqdir = './sequence_files/' # directory with sequence files
    maxgaps = 100 # maximum gaps to retain sequence if doesn't have mutation
    minseqlength = 1000 # only consider sequences if at least this long
    codonphyml = 'codonphyml' # codonPhyML program -- assumes in current path
    needlepath = '/Users/jbloom/EMBOSS-6.5.7/emboss/needle' # needle path
    outlierfile = '%s/outliers.txt' % seqdir # lists outlier strains
    strainmatch = re.compile('cds\:\w+ (?P<strain>[\?\:\'\/\-\(\)\. \w]+) ((H\d{1,2}){0,1}N\d{1,2}) ')
    protalignmentfile = 'protein_alignment.fasta'
    cdsalignmentfile = 'cds_alignment.fasta'
    sitecountsfile = 'site_counts.txt'
    identicalseqsfile = 'identical_sequences.txt'
    hasmutfile = 'seqs_with_mutation.txt'
    strainfile = 'strain_names.txt'
    codonphymlseqfile = 'sequences_for_tree.fasta'
    treeprefix = 'codonphyml'
    finaltreefile = '%s_tree.newick' % treeprefix

    # begin execution
    for (seqset, refseqhead, mutation, firstyear, buildtree, maxperyear, minclosest) in seqsets:
        print "\nPerforming analysis for mutation %s in %s" % (mutation, seqset)
        seqfile = "%s/%s.fasta" % (seqdir, seqset)
        seqs = mutation_analysis_funcs.ReadFASTA(seqfile)
        random.seed(1) # seed random number generator for reproducible output
        random.shuffle(seqs) # shuffle so there are not biases if only some are retained
        print "Read %d sequences from %s" % (len(seqs), seqfile)
        seqs = RemoveOutliers(seqs, outlierfile)
        print "After removing outliers specified in %s, there are %d sequences remaining." % (outlierfile, len(seqs))
        seqs = [(head, seq.upper()) for (head, seq) in seqs if len(seq) >= minseqlength]
        print "Retained %d sequences with length >= %d." % (len(seqs), minseqlength)
        dir = './%s_%s/' % (seqset, mutation)
        print "This analysis will be performed in directory %s" % dir
        if os.path.isdir(dir):
            print "Directory %s already exists, so the analysis will NOT be performed for this sequence set and mutation.\nIf you want to perform a new analysis, delete the existing directory." % dir
            continue
        print "Creating this directory."
        os.mkdir(dir)
        currdir = os.getcwd()
        os.chdir(dir)
        refseq = [(head, seq) for (head, seq) in seqs if head == refseqhead]
        if len(refseq) == 1:
            refseq = refseq[0]
            print "Pairwise alignments will be made to the reference sequence %s" % refseqhead
        else:
            raise ValueError("Failed to find single unique refseq")
        print "Now aligning the %d sequences using %s" % (len(seqs), needlepath)
        (prot_alignment, cds_alignment) = mutation_analysis_funcs.NeedleCDSandProtAlignments(refseq[1], seqs, needlepath)
        print "Aligned a total of %d sequences. Writing to %s and %s." % (len(seqs), protalignmentfile, cdsalignmentfile)
        mutation_analysis_funcs.WriteFASTA(prot_alignment, protalignmentfile)
        mutation_analysis_funcs.WriteFASTA(cds_alignment, cdsalignmentfile)
        (has_mut, no_mut, site_counts) = ClassifyByMutation(cds_alignment, prot_alignment, mutation)
        print "Overall, %d sequences have %s and %d do not.\nWriting the site counts to %s.\nWriting a list of the sequences with the mutation to %s" % (len(has_mut), mutation, len(no_mut), sitecountsfile, hasmutfile)
        open(sitecountsfile, 'w').write('Counts of identities at residue %s:\n%s' % (mutation[1 : -1], '\n'.join(['%d have %s' % (y, x) for (x, y) in site_counts.iteritems()])))
        open(hasmutfile, 'w').write('Sequences with %s:\n%s' % (mutation, '\n'.join([head for (head, seq) in has_mut])))
        has_mut = [(head, seq) for (head, seq) in has_mut if GetDate(head) >= firstyear]
        no_mut = [(head, seq) for (head, seq) in no_mut if GetDate(head) >= firstyear]
        print "After removing sequences from years before %d, have %d with mutation %s and %d without that mutation." % (firstyear, len(has_mut), mutation, len(no_mut))
        no_mut = RemoveIdentical(no_mut, identicalseqsfile)
        print "After removing identical sequences and substrings from those without %s, have %d unique sequences remaining without the mutation. Wrote information about the identical sequences to %s." % (mutation, len(no_mut), identicalseqsfile)
        print "Now finding the closest sequence with %s to every sequence with %s." % (mutation, mutation)
        (closest, not_closest) = GetClosest(has_mut, no_mut, minclosest, printoutput=True)
        print "Overall, categorized %d of the sequences without %s as being closest neighbors to a sequence with the mutation; the other %d are not closest." % (len(closest), mutation, len(not_closest))
        not_closest = GetNPerYear(not_closest, maxperyear)
        print "For the sequences without the mutation and not among the closest to a sequence with the mutation, retaining a maximum of %d per year. This gives a total of %d such sequences for all years." % (maxperyear, len(not_closest))
        combinedseqs = []
        strains = []
        has_mut_heads = dict([(head, seq) for (head, seq) in has_mut])
        for (head, seq) in has_mut + closest + not_closest:
            m = strainmatch.search(head)
            if not m:
                raise IOError("Failed to match strain in %s" % head)
            date = GetDate(head)
            if head in has_mut_heads:
                head += '_HAS_%s' % mutation
            head += '_%.2f' % date
            strain = m.group('strain')
            for char in [":", ",", ")", "(", ";", "]", "[", "'"]:
                head = head.replace(char, '')
            head = head.replace(' ', '_')
            combinedseqs.append((head, seq))
            strains.append((head, strain))
        print "Writing the strain names versus the full names to %s" % strainfile
        open(strainfile, 'w').write('\n'.join(['%s %s' % tup for tup in strains]))
        if not buildtree:
            print "Not building a tree, so done with this sequence set."
            continue
        print "Retained a total of %d sequences from all categories. Writing to %s" % (len(combinedseqs), codonphymlseqfile)
        mutation_analysis_funcs.WriteFASTA(combinedseqs, codonphymlseqfile)
        print "Now running %s to build a tree..." % codonphyml
        mutation_analysis_funcs.RunCodonphyml(codonphymlseqfile, codonphyml, 1, treeprefix, 'GY94_CF3x4')
        print "%s complete." % codonphyml
        if not os.path.isfile(finaltreefile):
            raise IOError("Failed to find the final tree of %s." % finaltreefile)
        print "The final tree has been written to the file %s.\n" % finaltreefile
        os.chdir(currdir)
    print "\nScript complete."


main() # run the script
