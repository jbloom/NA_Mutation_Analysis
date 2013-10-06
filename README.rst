=========================
NA_Mutation_Analysis
=========================

.. contents:: Contents
   :depth: 2

Summary
-----------

This is an analysis of the occurrences of the G147R mutation in influenza lineages with N1 neuraminidases (NAs). It examines the prevalence of this mutation in NA genes, and constructs phylogenetic trees. 

The G147R mutation is of interest because it confers receptor-binding activity on NA as described in `Hooper and Bloom, 2013`_. This analysis examines G147R in the three lineages that `Hooper and Bloom, 2013`_ reported to contain G147R:

* `pdmH1N1 NAs with G147R`_: human 2009 swine-origin pandemic H1N1 lineage.

* `seasonal H1N1 NAs with G147R`_: human seasonal H1N1 (up through 2009, then this lineage went extinct).

* `chicken H5N1 NAs with G147R`_

Note that the phylogenetic analysis here should be superior to that reported in `Hooper and Bloom, 2013`_. In that paper, the trees were built using protein sequences and only contained a subset of non-G147R sequences. Here the trees are built using codon sequences, and are constructed to ensure that the closest non-G147R neighbors of every G147R sequence is included.

The source code and data for this analysis are available `on GitHub`_. You can download the analysis by clipping the the *Download ZIP* button on the right-hand side `on GitHub`_.

This analysis was performed by `Jesse Bloom`_.


Software used
---------------

This analysis consists of a `Python`_ script. It has been tested with `Python`_ versions 2.6 and 2.7, and probably works with many other 2.* versions as well. The `Python`_ script is *mutation_analysis.py*. It uses functions defined in the module *mutation_analysis_funcs.py*.

The analysis uses the EMBOSS `needle`_ program for pairwise alignments, EMBOSS version 6.5.7.0.

The phylogenetic trees are built with `codonPhyML`_ version 1.00 201306.18

The trees are visualized using `FigTree`_ version 1.4.0.

The trees are re-rooted with `Path-O-Gen`_ version 1.4.


Input sequence files
-----------------------

Here is a summary of the input sequence data files used for the analysis. These are found in the subdirectory ``./sequence_files/`` in the main repository `on GitHub`_:

* *pdmH1N1_NAs.fasta* is the set of all human pandemic H1N1 influenza A NA protein-coding sequences as downloaded from the `Influenza Virus Resource`_ on October-4-2013.

* *seasonalH1N1_NAs.fasta* is the set of all human seasonal H1N1 influenza A NA protein-coding sequences as downloaded from the `Influenza Virus Resource`_ on October-4-2013.

* *chickenH5N1_NAs.fasta* is the set of all chicken H5N1 influenza A NA protein-coding sequences as downloaded from the `Influenza Virus Resource`_ on October-5-2013.

* *outliers.txt* is a list of sequences that are strong outliers on the phylogenetic tree due to branch placement or molecular clock analysis of trees with `Path-O-Gen`_. These sequences probably come from other lineages or are mis-annotated (for example, pdmH1N1 lineages that are labeled as seasonal H1N1). Any sequences in this file are excluded from analysis.


Analysis algorithm
--------------------------------------------

Here is a description of the algorithm used for the analysis. All of the steps of this analysis except for the final annotation of the phylogenetic trees for visual display are completely automated. So as long as you have installed the appropriate programs listed under `Software used`_, you can repeat this analysis by running the main `Python`_ script with the command::

    python mutation_analysis.py


For each of the sequence sets, the script performs the following steps. The first seven steps are algorithmic. The final step involves opening the tree in `FigTree`_ and manually annotating its appearance for better visual display.

The analysis for each sequence set is performed in a separate subdirectory.

The steps are:

    1) After removing any outliers specified in *outliers.txt*, all sequences of at least 1000 nucleotides are pairwise aligned to the reference sequence using `needle`_, and any gaps relative to the reference sequence are removed. The nucleotide sequences are aligned from protein alignments. The alignments are in the files *cds_alignment.fasta* and *protein_alignments.fasta*.

    2) The names of all strains with the mutation of interest are written to the file *seqs_with_mutation.txt*.

    3) The counts of different amino-acid identities at the site of interest are written to the file *site_counts.txt*.
    
    4) For all sequence sets (with and without the mutation of interest), all sequences with dates before the specified cutoff year are removed. 

    5) For all sequences **without** the mutation of interest, all redundant nucleotide sequences, sequences that are substrings of other sequences, or sequences with more than 100 gaps in the nucleotide alignment are removed. Lists of the names of sequences that are identical to other sequences are written to *identical_sequences.txt*.

    6) The sequences to retain for the phylogenetic analysis are then chosen. Retaining all of the sequences gives too many too build the trees, so just a subset is retained and written to *sequences_for_tree.txt*. The sequences that are retained are as follows:
    
        a) All sequences **with** the mutation of interest.
        
        b) The "closest" sequences **without** the mutation to each sequence **with** the mutation of interest. The purpose of this subset is to make sure that the closest sequences without the mutation are also in the trees. For each sequence with the mutation of interest, the closest sequences without the mutation based on fractional nucleotide identity are identified. If multiple sequences are tied for closest, we keep them all. We then also keep a set number (specified below for each analysis) of additional next-closest sequences.

        c) A random sampling of a set number of additional sequences per year of sequences **without** the mutation. The number per year for each analysis is specified below. The purpose here is to give an additional sampling of sequences without the mutation.

    7) `codonPhyML`_ is used to build a maximum-likelihood tree of the retained sequences. The tree is built using the `Goldman and Yang 1994`_ codon model with a single transition-transversion ratio estimated by maximum likelihood, the `CF3x4`_ method for estimating the equilibrium codon frequencies empirically, and a single omega (the dN/dS ratio) estimated by maximum likelihood. The final tree is saved with the name *codonphyml_tree.newick*. The tree includes branch supports calculated using the `SH-aLRT`_ method.

    8) The last step is manual annotation of the phylogenetic tree for visual appearance. This does not change the result, but does make it easier to look at. The formatted tree is saved as *codonphmyl_tree_formatted.trees* in a format that can be opened by `FigTree`_. The branches and tips are colored based on the presence / absence of the mutation. The tree is re-rooted based on where an analysis with `Path-O-Gen`_ suggests that the root should be placed. Some clades are collapsed, and key branch supports are shown (these are the `SH-aLRT`_ values calculated by `codonPhyML`_). Finally, a PDF of each tree is saved with the name *codonphyml_tree_formatted.pdf*.


pdmH1N1 NAs with G147R
--------------------------

The subdirectory ``./pdmH1N1_NAs_G147R/`` contains an analysis of G147R in viruses in the human 2009 pandemic H1N1 lineage. 

For the phylogenetic tree, in addition to the sequences with G147R, we retain the closest six sequences to each sequence with the mutation, and 10 randomly chosen sequences without the mutation per year. We keep sequences from 2009 and all subsequent years.

**Overall summary:** G147R is present in a few strains, but no evidence of evolutionary clusters. Some strains also have the intriguing G147E mutation.

Here are the counts of the amino-acid identities at site 147:

 * 8400 have G

 * 3 have R

 * 5 have E

 * 1 has  X

Here are all the sequences with G147R::

    cds:ADN26074 A/Finland/614/2009 H1N1 2009/07/24 NA
    cds:AFE11259 A/Tianjinhedong/SWL44/2011 H1N1 2011/01/08 NA
    cds:AFN20030 A/Singapore/SGH02/2011 H1N1 2011/02/09 NA


Here is the tree (the file ``./pdmH1N1_NAs_G147R/codonphyml_tree_formatted.pdf``). If your browser does not display the PDF embedded in HTML, click on the link to see the PDF alone:

.. figure:: ./pdmH1N1_NAs_G147R/codonphyml_tree_formatted.pdf
   :align: center
   :scale: 50%
   :alt: ./pdmH1N1_NAs_G147R/codonphyml_tree_formatted.pdf
   :target: ./pdmH1N1_NAs_G147R/codonphyml_tree_formatted.pdf

   Maximum-likelihood phylogenetic tree of G147R in pdmH1N1. Sequences with G147R are colored red. Support values for key branches are shown (`SH-aLRT`_ support). Some nodes are collapsed for visual display.


seasonal H1N1 NAs with G147R
-------------------------------

The subdirectory ``./seasonalH1N1_NAs_G147R/`` contains an analysis of G147R in viruses in the human seasonal H1N1 lineage (which circulated from 1918 to 1957, and then 1957 to 2009 in humans).

For the phylogenetic tree, in addition to the sequences with G147R, we retain the closest two sequences to each sequence with the mutation, and six randomly chosen sequences without the mutation per year. The tree includes only sequences from 2007 to 2009, because (with one exception) G147R only appeared starting in 2007.

**Overall summary:** G147R is present in 20 strains, with all but one of the occurrences in 2007 and later. There are some small phylogenetic clusters of G147R sequences.

Here are the counts of the amino-acid identities at site 147:

 * 4654 have G

 * 20 have R

 * 5 have X

Here are all the sequences with G147R::

    cds:ABD78030 A/South Canterbury/59/2000 H1N1 2000/09/06 NA

    cds:ABX58495 A/Tennessee/UR06-0238/2007 H1N1 2007/02/12 NA
    cds:ACY01424 A/Hamedan/117/2007 H1N1 2007/11/25 NA
    cds:ACA33659 A/Texas/74/2007 H1N1 2007/11/26 NA

    cds:ADA69512 A/Austria/404821/2008 H1N1 2008/01/21 NA
    cds:ACM17331 A/Austria/404811/2008 H1N1 2008/01/21 NA
    cds:ADA69518 A/Austria/405179/2008 H1N1 2008/01/23 NA
    cds:ACI94940 A/Austria/406109/2008 H1N1 2008/01/27 NA
    cds:BAH22142 A/Yokohama/30/2008 H1N1 2008/01/28 NA
    cds:ACM90850 A/Johannesburg/279/2008 H1N1 2008/07/22 NA
    cds:ADP89155 A/Thailand/Siriraj-05/2008 H1N1 2008// NA
    cds:ADZ53071 A/Hong Kong/01045/2008 H1N1 2008// NA
    cds:ADP89151 A/Thailand/Siriraj-01/2008 H1N1 2008// NA
    cds:ADP89152 A/Thailand/Siriraj-02/2008 H1N1 2008// NA

    cds:ADC45782 A/Niigata/08F188/2009 H1N1 2009/01/26 NA
    cds:AET84319 A/Iraq/WRAIR1683P/2009 H1N1 2009/03/ NA
    cds:ADA71159 A/Novosibirsk/3/2009 H1N1 2009/04/ NA
    cds:ACU44235 A/Kentucky/08/2009 H1N1 2009/05/12 NA
    cds:ACU44027 A/Kentucky/08/2009 H1N1 2009/05/12 NA
    cds:ADZ53099 A/Hong Kong/17566/2009 H1N1 2009// NA

Here is the tree (the file ``./seasonalH1N1_NAs_G147R/codonphyml_tree_formatted.pdf``). If your browser does not display the PDF embedded in HTML, click on the link to see the PDF alone:

.. figure:: ./seasonalH1N1_NAs_G147R/codonphyml_tree_formatted.pdf
   :align: center
   :scale: 90%
   :alt: ./seasonalH1N1_NAs_G147R/codonphyml_tree_formatted.pdf
   :target: ./seasonalH1N1_NAs_G147R/codonphyml_tree_formatted.pdf

   Maximum-likelihood phylogenetic tree of G147R in seasonal H1N1. Sequences with G147R are colored red. Support values for key branches are shown (`SH-aLRT`_ support). Some nodes are collapsed for visual display.


chicken H5N1 NAs with G147R
-------------------------------

The subdirectory ``./chickenH5N1_NAs_G147R/`` contains an analysis of G147R in viruses in chicken H5N1 viruses.

For the phylogenetic tree, in addition to the sequences with G147R, we retain the closest four sequences to each sequence with the mutation, and 10 randomly chosen sequences without the mutation per year. The tree includes only sequences from 2004 and later because the first G147R virus in chicken H5N1 is an a sequence from 2004.

**Overall summary:** G147R is present in 8 strains. There are some small phylogenetic clusters of G147R sequences.

Here are the counts of the amino-acid identities at site 147:

 * 1242 have G

 * 8 have R

 * 2 have E

Here are the sequences with G147R::

    cds:ADG59211 A/chicken/Gansu/44/2004 H5N1 2004// NA
    cds:ADG59204 A/chicken/Anhui/39/2004 H5N1 2004// NA
    cds:ADB26210 A/chicken/Nigeria/08RS848-93/2007 H5N1 2007/07/21 NA
    cds:AFH53768 A/chicken/Egypt/Kalyobia-18-CLEVB/2011 H5N1 2011/02/10 NA
    cds:AGG52920 A/chicken/Bangladesh/12VIR-7140-1/2011 H5N1 2011/12/19 NA
    cds:AGG52921 A/chicken/Bangladesh/12VIR-7140-2/2012 H5N1 2012/01/02 NA
    cds:AGG52922 A/chicken/Bangladesh/12VIR-7140-3/2012 H5N1 2012/01/08 NA
    cds:AGG52925 A/chicken/Bangladesh/12VIR-7140-6/2012 H5N1 2012/02/14 NA

Here is the tree (the file ``./chickenH5N1_NAs_G147R/codonphyml_tree_formatted.pdf``). If your browser does not display the PDF embedded in HTML, click on the link to see the PDF alone:

.. figure:: ./chickenH5N1_NAs_G147R/codonphyml_tree_formatted.pdf
   :align: center
   :scale: 75%
   :alt: ./chickenH5N1_NAs_G147R/codonphyml_tree_formatted.pdf
   :target: ./chickenH5N1_NAs_G147R/codonphyml_tree_formatted.pdf

   Maximum-likelihood phylogenetic tree of G147R in chicken H5N1. Sequences with G147R are colored red. Support values for key branches are shown (`SH-aLRT`_ support). Some nodes are collapsed for visual display.




.. _`on GitHub`: https://github.com/jbloom/NA_Mutation_Analysis
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Influenza Virus Resource`: http://www.ncbi.nlm.nih.gov/genomes/FLU/FLU.html
.. _`4JTV`: http://www.rcsb.org/pdb/explore.do?structureId=4JTV
.. _`HA_numbering`: https://github.com/jbloom/HA_numbering
.. _`4HMG`: http://www.rcsb.org/pdb/explore.do?structureId=4HMG
.. _`3TI6`: http://www.rcsb.org/pdb/explore/explore.do?structureId=3TI6
.. _`Zhang et al, 2013`: http://jvi.asm.org/content/87/10/5949
.. _`Vavricka et al, 2011`: http://www.plospathogens.org/article/info%3Adoi%2F10.1371%2Fjournal.ppat.1002249
.. _`PyMol`: http://www.pymol.org/
.. _`needle`: http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
.. _`Jensen-Shannon divergence`: https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
.. _`Python`: http://www.python.org/
.. _`codonPhyML`: http://sourceforge.net/projects/codonphyml/
.. _`RAxML`: https://github.com/stamatak/standard-RAxML
.. _`FigTree`: http://tree.bio.ed.ac.uk/software/figtree/
.. _`Path-O-Gen`: http://tree.bio.ed.ac.uk/software/pathogen/
.. _`CF3x4`: http://www.plosone.org/article/info%3Adoi/10.1371/journal.pone.0011230
.. _`Goldman and Yang 1994`: http://mbe.oxfordjournals.org/content/11/5/725.full.pdf
.. _`SH-aLRT`: http://sysbio.oxfordjournals.org/content/60/5/685.long
.. _`Hooper and Bloom, 2013`: http://www.ncbi.nlm.nih.gov/pubmed/24027333
