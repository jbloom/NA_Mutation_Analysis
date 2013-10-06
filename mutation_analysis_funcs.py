"""Module containing functions for mutation analysis.

Jesse Bloom, 2013."""


import os
import re
import sys
import shutil
import random
import cPickle
import cStringIO
import datetime
import tempfile
import subprocess



def WritePhylipSequenceFile(sequences, f, add_mm_line=False):
    """Writes a sequence file in sequential relaxed Phylip format.

    *sequences* is a list of 2-tuples *(head, seq)* where *head* is the
        sequence header and *seq* is the sequence (both strings). The sequences
        should be aligned, and so should all be the same length.
        
    *f* is a writeable file-like object to which the Phylip formatted
        sequences are written.

    *add_mm_line* specifies that we add the final line necessary for using
    ``codonPhyML`` with multi-models. If this option has its default value
    of *None*, then no such line is added. Otherwise, *add_mm_line* should
    be a string of length equal to the length of the sequences in *sequences*.
    The characters in this line are typically taken to specify the model 
    assigned to each site, and are written with the prefix ``#=GR mods``, 
    as in::

        #=GR mods AAAAABBAAAAAAAAAA

    Phylip file format is defined here: http://www.phylo.org/tools/phylip.html

    This function modifies that Phylip format as follows:

        * the sequence name (*head*) can be up to 100 characters long

        * there is always a space between the sequence name and the sequence

        * the following symbols are not allowed, blanks and any of :(),"

    If any of the headers violate these restrictions, or are not unique,
    an exception is raised.

    >>> sequences = [('head1', 'ATGACT'), ('head2', 'ATGCTA')]
    >>> f = cStringIO.StringIO()
    >>> WritePhylipSequenceFile(sequences, f)
    >>> f.seek(0)
    >>> print f.read()
    2 6
    head1 ATGACT
    head2 ATGCTA
    <BLANKLINE>

    >>> sequences = [('head1', 'ATGACT'), ('head1', 'ATGCTA')]
    >>> f = cStringIO.StringIO()
    >>> WritePhylipSequenceFile(sequences, f)
    Traceback (most recent call last):
        ...
    ValueError: Duplicate name of head1

    """
    maxnamelength = 100
    notallowed = ' :(),"'
    names_used = {}
    length = len(sequences[0][1])
    f.write("%d %d\n" % (len(sequences), length))
    for (name, seq) in sequences:
        if length != len(seq):
            raise ValueError("Sequences not all of the same length")
        if len(name) >= maxnamelength:
            raise ValueError("Name is too long:\n%s" % name)
        for x in notallowed:
            if x in name:
                raise ValueError("Name contains dis-allowed character %s" % x)
        if name in names_used:
            raise ValueError("Duplicate name of %s" % name)
        names_used[name] = True
        f.write("%s %s\n" % (name, seq))
    if add_mm_line:
        if not (isinstance(add_mm_line, str) and len(add_mm_line) == length):
            raise ValueError("add_mm_line is not a string of the same length as the sequences: %d and %d" % (len(add_mm_line), length))
        f.write("#=GR mods %s\n" % add_mm_line)


def RunCodonphyml(seqfile, path, seed, outprefix, model, add_mm_line=None):
    """Runs ``codonphyml``, returns log likelihood of tree.

    CALLING VARIABLES:

    * *seqfile* is a FASTA file with the sequences that we want to
      analyze. Each must have a unique header of less than 100 characters
      that does not include a blank or the characters :(),"

    * *path* is the path to the ``codonphyml`` executable. 

    * *seed* is the integer random number seed.

    * *outprefix* is the prefix for the files created by this function.
      These files have that prefix with the following suffixes. 
      Any existing files with these names are overwritten:

        - ``_config.drw`` is the darwin format input file used to run
          ``codonphyml``.

        - ``_output.txt`` is the output of ``codonphyml``.

        - ``_tree.newick`` is the inferred tree in Newick format.

        - ``_stats.txt`` contains the statistics from the inference,
          including all parameter values.

        - ``_loglikelihood.txt`` contains a single number representing
          the log likelihood.

    * *model* specifies the codon substitution model used. Valid values are:
    
        - *GY94_CF3x4* : The Goldman-Yang codon model with 
          the equilibrium codon frequencies estimated empirically using 
          the CF3x4 method, with kappa (the transition / transversion ratio) 
          estimated by maximum likelihood, and with omega (the dN/dS ratio) 
          estimated as a single parameter by maximum likelihood

        - *GY94_CF3x4_omega-gamma4* : The Goldman-Yang codon model with 
          the equilibrium codon frequencies estimated empirically using 
          the CF3x4 method, with kappa (the transition / transversion ratio) 
          estimated by maximum likelihood, and with omega (the dN/dS ratio) 
          drawn from four gamma-distributed categories with the distribution 
          shape parameter (alpha) estimated by maximum likelihood.

        - *KOSI07_F_omega-gamma4 : The Kosiol et al, 2007 empirical codon model 
          with the equilibrium codon frequencies estimated empirically using 
          the *F* method, with *kappa(tv)* (the relative decrease in transversions 
          versus transitions) estimated by maximum likelihood, and with *omega*
          (the elevation in nonsynonymous over synonymous) drawn from four 
          gamma-distributed categories with the distribution shape parameter 
          and mean estimated by maximum likelihood. This is based on 
          the *ECM+F+omega+1kappa(tv)* model described by Kosiol et al, 2007.

    * *add_mm_line* is *None* by default. If it is set to another
      value, it should be set to a string specifying a multi-model
      line in the format described in the documentation for the 
      parameter of the same name in *WritePhylipSequenceFile*.

    RETURN VARIABLE:

    This function returns a number giving the calculated log likelihood.

    It also creates the files specified by *outprefix*.

    """
    # define file names
    seqs = ReadFASTA(seqfile)
    phylipfile = tempfile.mkstemp(suffix='phylip')[1]
    WritePhylipSequenceFile(seqs, open(phylipfile, 'w'), add_mm_line=add_mm_line)
    statsfile = "%s_codonphyml_stats.txt" % phylipfile
    treefile = "%s_codonphyml_tree.txt" % phylipfile
    outfile = "%s_output.txt" % outprefix
    outtreefile = "%s_tree.newick" % outprefix
    outstatsfile = "%s_stats.txt" % outprefix
    llfile = "%s_loglikelihood.txt" % outprefix
    configfile = "%s_config.drw" % outprefix
    # create darwin format config file
    configlines = [
            "cpconfig := table();", # options in table called cpconfig
            "cpconfig['inputfile'] := '%s';" % phylipfile, # sequences
            "cpconfig['oformat'] := 'txt';", # output if txt format
            "cpconfig['sequential'] := true;", # sequences sequential phylip format
            "cpconfig['multiple'] := 1;", # number of datasets to analyze
            "cpconfig['bootstrap'] := -4;", # SH-aLRT branch supports
            "cpconfig['search'] := 'NNI';", # tree topology search
            "cpconfig['optimize'] := 'tlr';", # optimize topology (t), branch length (l), and rates (r)
            "cpconfig['quiet'] := true;", # no interactive questions asked
            "cpconfig['modrates'] := true;", # not sure what this option does
            "cpconfig['r_seed'] := %d;" % seed, # random number seed
            "cpconfig['datatype'] := 'codon';", # codon data
            ]
    if model == 'GY94_CF3x4_omega-gamma4':
        configlines += [
            "cpconfig['model'] := 'GY';", # Goldman-Yang model
            "cpconfig['frequencies'] := 'empirical';", # empirical codon frequencies
            "cpconfig['fmodel'] := 'CF3X4';", # CF3X4 codon frequency model
            "cpconfig['kappa'] := 'e';", # ML estimation of transition/transversion ratio
            "cpconfig['omega'] := 'DGAMMA';", # discrete gamma model for dN/dS
            "cpconfig['wclasses'] := 4;", # number of omega rate categories
            "cpconfig['alpha'] := 'e';", # ML estimation of gamma shape parameter
            ]
    elif model == 'GY94_CF3x4':
        configlines += [
            "cpconfig['model'] := 'GY';", # Goldman-Yang model
            "cpconfig['frequencies'] := 'empirical';", # empirical codon frequencies
            "cpconfig['fmodel'] := 'CF3X4';", # CF3X4 codon frequency model
            "cpconfig['kappa'] := 'e';", # ML estimation of transition/transversion ratio
            "cpconfig['omega'] := 'DM0';", # discrete gamma model for dN/dS
            "cpconfig['alpha'] := 'e';", # ML estimation of gamma shape parameter
            ]
    elif model == 'KOSI07_F_omega-gamma4':
        configlines += [
            "cpconfig['model'] := 'GY';", # Goldman-Yang model for codons
            "cpconfig['qrates'] := 'KOSI07';", # Kosiol 2007 empirical codon model
            "cpconfig['frequencies'] := 'empirical';", # empirical codon frequencies
            "cpconfig['fmodel'] := 'F1XCODONS';", # F codon frequency model
            "cpconfig['kappa'] := 'KAP3';", # ML estimation of kappa(tv)
            "cpconfig['omega'] := 'DGAMMA';", # discrete gamma model for dN/dS
            "cpconfig['wclasses'] := 4;", # number of omega rate categories
            "cpconfig['alpha'] := 'e';", # ML estimation of gamma shape parameter
            ]
    else:
        raise ValueError("Unrecognized model specification of %s" % model)
    open(configfile, 'w').write('\n'.join(configlines))
    # now run codonphyml
#    cmds = [path, '-i', phylipfile, '-q', '-d', 'codon', '--r_seed', str(seed), '-b', '0', '--quiet']
    cmds = [path, '--darwinconfig', configfile]
    try:
        try:
            p = subprocess.Popen(cmds, stdout=open(outfile, 'w'), stderr=subprocess.PIPE, stdin=subprocess.PIPE)
            (stdoutdata, stderrdata) = p.communicate('Y\n')
        except OSError:
            raise ValueError("codonphyml failed to run, probably the executable is not callable using the following path: %s\n%s" % (path, str(stderrdata)))
        # capture output
        if stderrdata:
            raise ValueError("codonphyml encountered errors:\n%s" % stderrdata)
        if not os.path.isfile(statsfile):
            raise IOError("Failed to find codonphyml stats file %s" % (statsfile))
        if not os.path.isfile(treefile):
            raise IOError("Failed to find codonphyml tree file %s" % (treefile))
        open(outtreefile, 'w').write(open(treefile).read())
        stats = open(statsfile).read()
        open(outstatsfile, 'w').write(stats)
        llmatch = re.compile('\. Log-likelihood: (?P<ll>\-\d+(\.\d+){0,1})')
        m = llmatch.search(stats)
        if not m:
            raise ValueError("Failed to find Log-likelihood in:\n%s" % stats)
        ll = m.group('ll')
        open(llfile, 'w').write(ll)
        ll = float(ll)
    finally:
        for f in [phylipfile, statsfile, treefile]:
            if os.path.isfile(f):
                os.remove(f)
    return ll


def AmbiguousNTCodes(nt):
    """Returns all possible nucleotides corresponding to an ambiguous code.

    This method takes as input a single nucleotide character, which is
    assumed to represent a nucleotide as one of the accepted
    codes for an ambiguous character.  Returns a list giving
    all possible codes for which a nucleotide might stand.  Raises
    an exception if `nt` is not a valid nucleotide code.

    EXAMPLES:

    >>> AmbiguousNTCodes('N')
    ['A', 'T', 'G', 'C']

    >>> AmbiguousNTCodes('R')
    ['A', 'G']

    >>> AmbiguousNTCodes('A')
    ['A']

    >>> AmbiguousNTCodes('-')
    ['-']

    >>> AmbiguousNTCodes('F')
    Traceback (most recent call last):
       ...
    ValueError: Invalid nt code of "F"
    """
    if nt in ['A', 'T', 'G', 'C', '-']:
        return [nt]
    elif nt == 'R':
        return ['A', 'G']
    elif nt == 'Y':
        return ['T', 'C']
    elif nt == 'K':
        return ['G', 'T']
    elif nt == 'M':
        return ['A', 'C']
    elif nt == 'S':
        return ['G', 'C']
    elif nt == 'W':
        return ['A', 'T']
    elif nt == 'B':
        return ['C', 'G', 'T']
    elif nt == 'D':
        return ['A', 'G', 'T']
    elif nt == 'H':
        return ['A', 'C', 'T']
    elif nt == 'V':
        return ['A', 'C', 'G']
    elif nt == 'N':
        return ['A', 'T', 'G', 'C']
    else: 
        raise ValueError('Invalid nt code of "%s"' % nt)


def WriteFASTA(headers_seqs, filename, writable_file=False):
    """Writes sequences to a FASTA file.

    CALLING VARIABLES:

    `headers_seqs` : list of 2-tuples specifying sequences and their
    corresponding headers.  Each entry is the 2-tuple `(header, seq)`
    where `header` is a string giving the header (without the leading ">"),
    and `seq` is the corresponding sequence.

    `filename` : string that specifies the name of the file to which the
     headers and sequences should be written.  If this file already exists,
     it is overwritten. 

    `writable_file` : Boolean switch specifying that rather than `filename`
    giving a string specifying the name of a file to which the sequences
    should be written, it instead specifies a writable file object to which
    the sequences should be written.

    RESULT OF THIS FUNCTION:

    The sequences are written to the file in the same order that they are specified
    in `headers_seqs`.
    """
    assert isinstance(writable_file, bool)
    if writable_file:
        f = filename
    else:
        f = open(filename, 'w')
    for (header, seq) in headers_seqs:
        f.write(">%s\n%s\n" % (header, seq))
    if not writable_file:
        f.close()


def ReadFASTA(fastafile):
    """Reads sequences from a FASTA file.

    CALLING VARIABLE:

    `fastafile` : specify the name of a FASTA file.

    RETURN VARIABLE: 
    
    This function reads all sequences from the FASTA file.  It returns the
    list `headers_seqs`.  This list is composed of a 2-tuple `(header, seq)`
    for every sequence entry in `fastafile`. `header` is the header for
    a sequence, with the leading ">" and any trailing spaces removed. `seq`
    is the corresponding sequence.
    """
    lines = open(fastafile).readlines()
    headers_seqs = []
    header = None
    seq = []
    for line in lines:
        if line[0] == '>':
            if (not header) and (not seq):
                pass # first sequence in file
            elif header and not seq:
                raise ValueError, "Empty sequence for %s" % header
            elif seq and not header:
                raise ValueError, "File does not begin with header."
            else:
                seq = ''.join(seq)
                seq = seq.replace(' ', '')
                headers_seqs.append((header, seq))
            header = line.strip()[1 : ]
            seq = []
        else:
            seq.append(line.strip())
    if (not header) and (not seq):
        pass # first sequence in file
    elif header and not seq:
        raise ValueError, "Empty sequence for %s" % header
    elif seq and not header:
        raise ValueError, "File does not begin with header."
    else:
        seq = ''.join(seq)
        seq = seq.replace(' ', '')
        headers_seqs.append((header, seq))
    return headers_seqs

def Needle(s1, s2, needlecmd, seqtype, outfile, gapopen=10.0,
           gapextend=0.5, outformat='fasta', overwriteoutfile=False, 
           endweight=False):
    """Needleman-Wunsch pairwise alignment using EMBOSS ``needle``.

    This function uses the Needleman-Wunsch algorithm as implemented in
    EMBOSS by the ``needle`` program to perform a pairwise alignment of one
    or more sequences. One target sequence (`s1`) is pairwise aligned with
    one or more other sequences (`s2`). The results are written to a
    specified output file.

    The return value is `None`, but `outfile` is created.

    This function is currently tested with ``needle`` from EMBOSS 6.5.7.0,
    although it should run with many versions.

    `s1` : string giving the first input sequence, to which all sequences
    in `s2` are aligned.

    `s2` : specifies the second input sequence(s) to which we are
    aligning `s1`. It can be:

        * a string giving the second input sequence to align to `s1` (if we
          are aligning just one pair)

        * a list of strings giving all input sequences to align to `s1`
          (if we are making multiple pairwise alignments).
        
        * a 2-tuple specifying that we read `s2` from an existing file. In
          this case, the 2-tuple should be `(filename, fileformat)`
          where `filename` is a string giving the filename, and
          `fileformat` is a string giving the file format (such as
          'fasta' or 'fastq').

    `needlecmd` : string giving executable path for ``needle``. For example,
    this can be a full pathname (such as ``/Users/jbloom/EMBOSS-6.5.7/emboss/needle``) 
    or if ``needle`` is installed in a directory in the current search path it might 
    just be the name of the program (i.e. 'needle').  
    
    `seqtype` : string specifying the sequence type. Can be either 
    'nucleotide' or 'protein'.
    
    `outfile` : string giving the name of the file to which we write the
    output alignment. The format is specified by `outformat`, and
    whether we overwrite any existing outfile is specified by
    `overwriteoutfile`.

    `gapopen` : the gap opening penalty, 10.0 by default.

    `gapextend` : the gap extension penalty, 0.5 by default.

    `outformat` : string specifying the format of the output alignment.
    Is 'fasta' by default.

    `overwriteoutfile` : Boolean switch specifying whether we overwrite
    any existing output files of this name. If `True`, we overwrite
    existing files. If `False`, we raise an `IOError` if there is
    already an existing file with the name of `outfile`.

    `endweight` : the weight applied against gaps at the end of the
    alignment. Is `False` by default, meaning that no end gap penalty
    is applied. If it is set to some other value, it should be the
    2-tuple `(endopen, endextend)` where `endopen` is the penalty for
    opening an end gap, and endextend is the penalty for extending
    an end gap.
    """
    try:
        arguments = ['-gapopen %f -gapextend %f' % (gapopen, gapextend)]
        if seqtype == 'nucleotide':
            arguments.append('-snucleotide1 -snucleotide2')
        elif seqtype == 'protein':
            arguments.append('-sprotein1 -sprotein2')
        else:
            raise ValueError("Invalid value of seqtype: %r" % seqtype)
        (s1fd, s1file) = tempfile.mkstemp() # temporary file
        WriteFASTA([('s1', s1)], s1file)
        arguments.append('-sformat1 fasta')
        if isinstance(s2, tuple) and len(s2) == 2:
            s2file = s2[0]
            if not os.path.isfile(s2file):
                raise IOError("Cannot find s2 file of %s" % s2file)
            arguments.append('-sformat2 %s' % s2[1])
        elif isinstance(s2, str):
            (s2fd, s2file) = tempfile.mkstemp() # temporary file
            WriteFASTA([('s2', s2)], s2file)
            arguments.append('-sformat2 fasta')
        elif isinstance(s2, list):
            (s2fd, s2file) = tempfile.mkstemp() # temporary file
            WriteFASTA(zip(['s2'] * len(s2), s2),
                    s2file)
            arguments.append('-sformat2 fasta')
        else:
            raise ValueError("Invalid argument for s2: %r" % s2)
        if endweight:
            if isinstance(endweight, tuple) and len(endweight) == 2:
                arguments.append('-endweight Y -endopen %f -endextend\
                    %f' % endweight)
            else:
                raise ValueError("Invalid value for endweight: %r" %
                    endweight)
        else:
            arguments.append('-endweight N')
        if os.path.isfile(outfile) and not overwriteoutfile:
            raise IOError("Needle outfile of %s already exists." %
            outfile)
        arguments.append('-aformat3 %s' % outformat)
        (stderrfd, stderrfile) = tempfile.mkstemp()
        (stdoutfd, stdoutfile) = tempfile.mkstemp()
        stdout = open(stdoutfile, 'w')
        stderr = open(stderrfile, 'w')
        cmdstr = '%s %s -asequence %s -bsequence %s -outfile %s' % \
                (needlecmd, ' '.join(arguments), s1file, s2file, outfile)
        returncode = subprocess.call(cmdstr.split(), stdout=stdout,
                stderr=stderr)
        stdout.close()
        stderr.close()
        if returncode != 0:
            raise IOError("Execution of needle failed. \n\nOutput: \n%s\
                \n\n Errors: \n%s" % ((open(stdoutfile).read()),
                open(stderrfile).read()))
    finally:
        # remove any temporary files created by function
        try:
            if os.path.isfile(s1file):
                os.remove(s1file) # remove the temporary file
                os.close(s1fd)
        except:
            pass
        try:
            if not isinstance(s2, tuple) and os.path.isfile(s2file):
                os.remove(s2file) # remove the temporary file
                os.close(s2fd)
        except:
            pass
        try:
            if os.path.isfile(stderrfile):
                os.remove(stderrfile) # remove the temporary file
                os.close(stderrfd)
        except:
            pass
        try:
            if os.path.isfile(stdoutfile):
                os.remove(stdoutfile) # remove the temporary file
                os.close(stdoutfd)
        except:
            pass

def StripGapsToFirstSequence(aligned_headers_seqs):
    """Strips gaps from a reference sequence, and all corresponding alignments.

    On input, *aligned_headers_seqs* is a list of two or more aligned sequences
    as *(header, sequence)* 2-tuples.

    The first sequence in this alignment is the reference sequence.
    The returned variable is a list similar to *aligned_headers_seqs*, but with
    all positions corresponding to gaps in this reference sequence stripped away.
    All gaps ('-') characters are removed from this reference sequence.  In addition,
    in all other aligned sequences in *aligned_headers_seqs*, every character at
    the same position as a gap in the reference sequence is removed.  Therefore,
    at the end of this procedure, all of the alignments have the same length
    as the reference sequence with its gaps stripped away.  The headers are 
    unchanged.  The order of sequences in this stripped alignment is also
    unchanged.

    >>> StripGapsToFirstSequence([('s1', '-AT-A-GC'), ('s2', 'AAT-TAGC'), ('s3', '--T-A-GC')])
    [('s1', 'ATAGC'), ('s2', 'ATTGC'), ('s3', '-TAGC')]
    """
    if not (isinstance(aligned_headers_seqs, list) and len(aligned_headers_seqs) >= 2):
        raise ValueError("aligned_headers_seqs does not specify at least two aligned sequences.")
    (ref_head, ref_seq) = aligned_headers_seqs[0]
    non_strip_positions = [] # positions not to strip away
    stripped_ref_seq = []
    for i in range(len(ref_seq)):
        r = ref_seq[i]
        if r != '-':
            non_strip_positions.append(i)
            stripped_ref_seq.append(r)
    stripped_headers_seqs = [(ref_head, ''.join(stripped_ref_seq))]
    for (ihead, iseq) in aligned_headers_seqs[1 : ]:
        istrippedseq = ''.join([iseq[i] for i in non_strip_positions])
        stripped_headers_seqs.append((ihead, istrippedseq))
    return stripped_headers_seqs



def DateToOrdinal(datestring, refyear):
    """Converts a date string to an ordinal date.

    *datestring* is a date given by a string such as '2007/2/13' (for
    Feb-13-2007), or '2007/2//' if no day is specified, or
    '2007//' if no day or month is specified. The '/' characters can
    also be '-'.

    *refdate* is an integer year from the approximate timeframe we are examining
    which is used to anchor the datestring date on the assumption
    that each year has 365.25 days.

    The returned value is a number (decimal) giving the date. If no
    day is specified, the 15th (halfway through the month) is chosen.
    If no month or day is specified, July 1 (halfway through the
    year) is chosen.

    >>> print "%.2f" % DateToOrdinal('2007/4/27', 1968)
    2007.32

    >>> print "%.2f" % DateToOrdinal('2007/4/', 1968)
    2007.29

    >>> print "%.2f" % DateToOrdinal('2007//', 1968)
    2007.50

    >>> print "%.2f" % DateToOrdinal('2007-4-27', 1968)
    2007.32
    """
    if not isinstance(refyear, int):
        raise ValueError('refyear is not an integer')
    refdate = datetime.date(refyear, 1, 1).toordinal()
    try:
        if '/' in datestring:
            (year, month, day) = datestring.split('/')
        else:
            (year, month, day) = datestring.split('-')
    except ValueError:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    if year and month and day:
        (year, month, day) = (int(year), int(month), int(day))
        date = datetime.date(year, month, day)
    elif year and month:
        (year, month) = (int(year), int(month))
        date = datetime.date(year, month, 15)
    elif year:
        year = int(year)
        date = datetime.date(year, 7, 1)
    else:
        raise ValueError("Invalid datestring of: %s" % str(datestring))
    return (date.toordinal() - refdate) / 365.25 + refyear



def FracIdent(s1, s2):
    """Returns fraction identities of two sequences, ignoring sites where either is '-'"""
    assert len(s1) == len(s2)
    nsame = ntotal = 0
    for (r1, r2) in zip(s1, s2):
        if r1 != '-' and r2 != '-':
            ntotal += 1
            if r1 == r2:
                nsame += 1
    return nsame / float(ntotal)


def Translate(headers_sequences, readthrough_n=False, readthrough_stop=False, truncate_incomplete=False, translate_gaps=False):
    """Translates a set of nucleotide sequences to amino acid sequences.

    CALLING VARIABLES:

    `headers_sequences` : list of tuples `(header, seq)` as would be returned
    by `Read`.  The sequences should all specify valid coding nucleotide
    sequences.  
    
    The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.

    `readthrough_n` : specifies that if any nucleotides
    in the sequence are equal to to an ambiguous nt code and cannot therefore
    be unambiguously translated into an amino acid, we simply translate through these 
    nucleotides by making the corresponding amino acid equal to "X".  By
    default, this option is `False`.  Note that even when this option is `False`,
    certain ambiguous nucleotides may still be translatable if they all lead to 
    the same amino acid.

    `readthrough_stop` : specifies that if we encounter any stop
    we simply translation them to 'X'.  By default,
    this option is `False`, meaning that we instead raise an error
    of an incomplete stop codon.

    `truncate_incomplete` : specifies that if the sequence
    length is not a multiple of three, we simply truncate off the one or two
    final nucleotides to make the length a multiple of three prior to translation.
    By default, this option is `False`, meaning that no such truncation is done.

    `translate_gaps` : specifies that a codon of '---' is translated
    to '-'. Codons with one '-' are also translated to gaps.

    RETURN VARIABLE:

    The returned variable is a new list in which
    all of the nucleotide sequences have been translated to their
    corresponding protein sequences, given by one letter codes.
    If any of the nucleotide sequences do not translate to
    valid protein sequences, an exception is raised.

    EXAMPLES:

    >>> Translate([('seq1', 'ATGTAA'), ('seq2', 'gggtgc')])
    [('seq1', 'M'), ('seq2', 'GC')]

    >>> Translate([('seq2', 'GGNTGC')])
    [('seq2', 'GC')]
    
    >>> Translate([('seq2', 'NGGTGC')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate codon NGG
    
    >>> Translate([('seq2', 'NGGTGC')], readthrough_n=True)
    [('seq2', 'XC')]

    >>> Translate([('seq2', 'TAATGC')])
    Traceback (most recent call last):
       ...
    ValueError: Premature stop codon

    >>> Translate([('seq2', 'TAATGC')], readthrough_stop=True)
    [('seq2', 'XC')]

    >>> Translate([('seq2', 'TGCA')])
    Traceback (most recent call last):
       ...
    ValueError: Sequence length is not a multiple of three

    >>> Translate([('seq2', 'TGCA')], truncate_incomplete=True)
    [('seq2', 'C')]

    >>> Translate([('seq2', 'TGC---')])
    Traceback (most recent call last):
       ...
    ValueError: Cannot translate gap.

    >>> Translate([('seq2', 'TGC---')], translate_gaps=True)
    [('seq2', 'C-')]
    """
    genetic_code = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
            'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'TAT':'Y', 'TAC':'Y', 'TAA':'STOP', 'TAG':'STOP',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N', 'AAC':'N',
            'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'TGT':'C', 'TGC':'C', 'TGA':'STOP', 'TGG':'W', 'CGT':'R',
            'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R',
            'AGG':'R', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}
    assert isinstance(headers_sequences, list)
    translated_headers_sequences = []
    for (head, seq) in headers_sequences:
        seq = seq.upper()
        if len(seq) % 3:
            if truncate_incomplete:
                seq = seq[ : -(len(seq) % 3)]
            else:
                raise ValueError, "Sequence length is not a multiple of three"
        prot_length = len(seq) // 3
        prot = []
        for i in range(prot_length):
            codon = seq[3 * i : 3 * (i + 1)]
            try:
                aa = genetic_code[codon]
            except KeyError:
                if '-' in codon:
                    if translate_gaps:
                        aa = '-'
                    else:
                        raise ValueError("Cannot translate gap.")
                else:
                    # see if we have an ambiguous nucleotide codon that doesn't matter in the translation
                    possible_nt1 = AmbiguousNTCodes(codon[0])
                    possible_nt2 = AmbiguousNTCodes(codon[1])
                    possible_nt3 = AmbiguousNTCodes(codon[2])
                    possible_codons = []
                    for nt1 in possible_nt1:
                        for nt2 in possible_nt2:
                            for nt3 in possible_nt3:
                                possible_codons.append("%s%s%s" % (nt1, nt2, nt3))
                    try:
                        aa = genetic_code[possible_codons[0]]
                    except KeyError:
                        raise KeyError("Cannot translate codon %s in %s" % (codon, head))
                    for possible_codon in possible_codons:
                        if genetic_code[possible_codon] != aa:
                            if readthrough_n:
                                aa = 'X'
                            else:
                                raise ValueError("Cannot translate codon %s" % codon)
            if aa == 'STOP' and i == prot_length - 1:
                aa = ''
            elif aa == 'STOP':
                if readthrough_stop:
                    aa = 'X'
                else:
                    raise ValueError("Premature stop codon")
            prot.append(aa)
        translated_headers_sequences.append((head, ''.join(prot)))
    return translated_headers_sequences


def FracIdentities(a):
    """Returns fraction identities among pairwise aligned sequences.

    *a* is a pairwise alignment of the form *[(head1, seq1), (head2, seq2)]*.
    In this tuple, *seq1* and *seq2* are of the same length and are aligned.

    Among all aligned sites (not gaps), returns the fraction of identities.
    """
    (s1, s2) = (a[0][1].upper(), a[1][1].upper())
    if len(s1) != len(s2):
        raise ValueError("The following two sequences differ in length:\n>%s\n%s\n>%s\n%s\nlengths are %d and %d" % (a[0][0], a[0][1], a[1][0], a[1][1], len(s1), len(s2)))
    nident = n = 0
    for (x, y) in zip(s1, s2):
        if x == '-' or y == '-':
            pass
        else:
            n += 1
            if x == y:
                nident += 1
    return nident / float(n)


def NeedleCDSandProtAlignments(refseq, seqs, needlecmd):
    """Uses EMBOSS needle to align sequences in *seqs* to *refseq*.

    The sequences in *seqs* and *refseq* should be coding DNA sequences.
    They must be translatable, but ambiguous nucleotides and truncation
    of incomplete sequences are allowed.

    *refseq* is a string giving the reference sequence.

    *seqs* is a list of *(header, sequence)* 2-tuples.

    *needlecmd* is the path to the EMBOSS needle program.

    Returns the 2-tuple *(prot_alignment, cds_alignment)*.

    *prot_alignment* contains each of the proteins encoded in 
    *seqs* as a *(header, protsequence)* 2-tuple, with the
    protein sequence aligned to that in *refseq* and with all
    gaps relative to *refseq* stripped away.

    *cds_alignment* is an alignment of the coding DNA sequences
    in *prot_alignment*, with the nucleotide alignments done
    according to the protein alignments.
    """
    prots = []
    heads = []
    refprot = Translate([('head', refseq)], truncate_incomplete=True, readthrough_n=True)[0][1]
    for (head, seq) in seqs:
        try:
            (head, prot) = Translate([(head, seq)], readthrough_n=True, truncate_incomplete=True)[0]
            heads.append(head)
            prots.append(prot)
        except:
            sys.stderr.write("PROBLEM translating sequence %s" % head)
            raise
    try:
        temp = tempfile.mkstemp()[1]
        os.remove(temp) # remove, needle with create again
        Needle(refprot, prots, needlecmd, 'protein', temp)
        alignments = ReadFASTA(temp)
    finally:
        if os.path.isfile(temp):
            os.remove(temp)
    assert len(alignments) == 2 * len(prots) == 2 * len(heads) == 2 * len(seqs)
    prot_alignment = []
    cds_alignment = []
    for i in range(len(prots)):
        prot = prots[i]
        head = heads[i]
        seq = seqs[i][1]
        assert seqs[i][0] == head
        (refa, prota) = (alignments[2 * i][1], alignments[2 * i + 1][1])
        assert len(refa) == len(prota)
        iref = iprot = 0
        alignedprot = []
        alignedcds = []
        for (aa_ref, aa_prot) in zip(refa, prota):
            assert (aa_ref == '-' or aa_ref == refprot[iref])
            assert (aa_prot == '-' or aa_prot == prot[iprot])
            if aa_ref == '-' and aa_prot != '-':
                iprot += 1
            elif aa_prot == '-' and aa_ref != '-':
                alignedprot.append(aa_prot)
                alignedcds.append('---')
                iref += 1
            elif aa_ref != '-' and aa_prot != '-':
                alignedprot.append(aa_prot)
                alignedcds.append(seq[3 * iprot : 3 * iprot + 3])
                iref += 1
                iprot += 1
            else:
                raise ValueError("Both prots in alignment have gap")
        alignedprot = ''.join(alignedprot)
        alignedcds = ''.join(alignedcds)
        assert alignedprot == Translate([(head, alignedcds)], readthrough_n=True, truncate_incomplete=True, translate_gaps=True)[0][1]
        assert len(alignedprot) == len(refprot)
        prot_alignment.append((head, alignedprot))
        cds_alignment.append((head, alignedcds))
    assert len(prot_alignment) == len(cds_alignment)
    return (prot_alignment, cds_alignment)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
