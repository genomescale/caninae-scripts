# Required software

Python3 and the biopython and python-docx packages. R and the ggplot2,
ggforce, ggthemes, ggstance, ggtree, deeptime, cowplot and Cairo packages.

# Downloading molecular loci

## Step 1

Enter the `molecular` folder.

## Step 2

Run `get-sequences.py` to download all sequences (except Lycaon pictus
sequences from Shang 2017). This script requires BioPython. You must change
the email address in the script to your own.

## Step 3

Download supplementary Lycaon pictus sequences `726_2017_2422_MOESM3_ESM.txt`
from https://doi.org/10.1007/s00726-017-2422-5

## Step 4

Run `append-supplementary.py` to combine the GenBank sequences with
supplementary sequences.

# Running PASTA on Tas2r sequences

## Step 1

Exerpt Tas2r sequences from combined.fasta into a separate file called
`tas2r.fasta` in the folder `tas2r`.

## Step 2

Enter the `tas2r` folder, and infer a gene tree (and alignment) of Tas2r
sequences using PASTA by running `run_pasta.py -i tas2r.fasta`

Clean up the labels by running `relabel.py`

Generate Supplementary Figure 12 by running `tas2r.R`

This figure was used to figure out what sequences to omit or to change
when writing the collate script in the next step.

# Preparing molecular loci

## Step 1

Re-enter the `molecular` folder and run `collate.py` to collate, align and
trim each locus. This script calls the programs `prank` and `trimal`, so both
must be installed and in your path. Errors in annotation are corrected by this
script.

## Step 2

Run `filter.py` to filter and phase loci. This script will, for loci in
`whitelist.csv`, choose one sequence per species per locus as directed in
`label-mapping.csv`. IUPAC ambiguity codes for two possible bases (e.g. R or
Y) are randomly resolved (e.g. R becomes A or G). The newly generated files
ending in `.filtered.fasta` should be used for further phylogenetic analysis.
This script also collects statistics on the alignments before and after
trimming, storing them in `msa-stats.csv`.

## Step 3

Run `plot-sources.R` to generate Supplementary Figure 13, which communicates
both the presence or absence of a sequence for a given species in the combined
data set, the source of the sequence when it is present, and the lengths of
multiple sequence alignments. The filename of the figure is
`filtered_sources.pdf`.

## Step 4

Run `make-alignment-table.py` to generate LaTeX code for the table of
alignment statistics.

# Preparing morphological character matrix

## Step 1

Enter the `morphology` folder.

## Step 2

Download `zsc12293-sup-0021-appendixs1.doc` and
`zsc12293-sup-0022-appendixs2.docx` supplementary files from
Zrzavy et al. (2018).

## Step 3

Run `make-nexus.py` to convert character state values and metadata in the
supplementary files into the character matrix `morphology.nex`. The species
names are corrected using the mapping in `label-mapping.csv`.

## Step 4

Add the missing state 'absent' before the 'present' state for the description
of state 197 'Caudal ventral iliac spine'. This state is omitted by the above
script due to an inconsistency in the original supplementary files.

# Preparing fossil taxa dates

## Step 1

Enter the `dates` folder. This contains the file
`raw-stratigraphic-ranges.csv`. That file lists the beginning and end X
coordinates of the stratigraphic ranges in Figure 4 of Zrzavy et al. (2018).

## Step 2

Run `process-stratigraphic-ranges.py`. This will generate a table of midpoint
dates for fossil taxa based on the previously measured ranges, called
`beast-dates.tsv`.

# Preparing XML files for FBD analyses

## Step 1

Use BEAUti to create a StarBEAST2 XML file. The data for this analysis should
be the filtered molecular loci, morphological matrix and midpoint dates. The
model should be as depicted in Supplementary Figure 5. Use a chain length of
1,073,741,824 and a log frequency of 131,072. The species tree file name
should be `sb2-total-evidence.x.trees`, the trace log file name should be
`sb2-total-evidence.x.log`, and the XML file name should be
`sb2-total-evidence.x.xml`. Make four copies of this file in the format in the
folder `sb2-total-evidence` replacing `x` with a, b, c or d to distinguish the
four xml, trace log and species tree files.

## Step 2

Comment the gene tree logging elements from each XML file. If you want to keep
the gene trees, make sure that every chain is writing gene tree samples to a
unique file path.

## Step 3

Create copies of each StarBEAST2 XML file in the folder
`concat-total-evidence`, replacing `sb2` with `concat` in the XML file names.
Edit the XML files to change `sb2` to `concat` for the species trees and trace
log file names, and to remove all aspects relating to gene trees including
StarBEAST2 operators.

# Preparing extant-only analysis

## Step 1

Run the `remove-extant-morphology.py` script to strip long-extinct species
from the morphology nexus file. This will create a new nexus file
`morphology-extant.nex` in the `morpology` folder containing only data from
extant or recently extinct taxa.

## Step 2

Again use BEAUti to create a StarBEAST2 XML file. This time tip dating is not
required, so only the filtered molecular loci and extant molecular data
should be used. Configure the model and priors to match the depiction in
Supplementary Figure 5, except using a Birth-Death tree prior conditioned on
the root instead of the origin, and without the FBD-specific parameters.
Again make four copies of the XML file, this time in the folder
`sb2-extant-only`using the filename `sb2-extant-only.x.xml`. Use the same
filename (with different extensions) for the tracelog and species trees
files. Replace `x` with a, b, c or d to distinguish the four chains.

## Step 3

BEAUti does not add morphology-only taxa to the `taxonsuperset` element. This
is not a problem, except when we want to include those taxa in a node
calibration. Since want to include *Canis rufus* and *Dusicyon australis* in
the calibration of Canini, add the following TaxonSet and Taxon elements for
those species to the `taxonsuperset` of each XML file:

```XML
<taxon id="Canis_rufus" spec="TaxonSet">
  <taxon id="Canis_rufus_x" spec="Taxon"/>
</taxon>
<taxon id="Dusicyon_australis" spec="TaxonSet">
  <taxon id="Dusicyon_australis_x" spec="Taxon"/>
</taxon>
```

Now add a node calibration element to the section of the XML containing
prior distributions:

```XML
<distribution id="StemCanis.prior" spec="beast.math.distributions.MRCAPrior" monophyletic="true" tree="@Tree.t:Species">
  <taxonset id="Canis" spec="TaxonSet">
    <taxon idref="Atelocynus_microtis"/>
    <taxon idref="Canis_aureus"/>
    <taxon idref="Canis_latrans"/>
    <taxon idref="Canis_lupus"/>
    <taxon idref="Canis_rufus"/>
    <taxon idref="Canis_simensis"/>
    <taxon idref="Cerdocyon_thous"/>
    <taxon idref="Chrysocyon_brachyurus"/>
    <taxon idref="Cuon_alpinus"/>
    <taxon idref="Dusicyon_australis"/>
    <taxon idref="Lupulella_adusta"/>
    <taxon idref="Lupulella_mesomelas"/>
    <taxon idref="Lycalopex_culpaeus"/>
    <taxon idref="Lycalopex_fulvipes"/>
    <taxon idref="Lycalopex_griseus"/>
    <taxon idref="Lycalopex_gymnocercus"/>
    <taxon idref="Lycalopex_sechurae"/>
    <taxon idref="Lycalopex_vetulus"/>
    <taxon idref="Lycaon_pictus"/>
    <taxon idref="Speothos_venaticus"/>
  </taxonset>
  <Exponential id="Exponential.1" name="distr" offset="5.0">
    <parameter id="RealParameter.135" spec="parameter.RealParameter" estimate="false" name="mean">1.66666666666</parameter>
  </Exponential>
</distribution>
```

And an appropriate log entry to the `tracelog` element, i.e.:

```XML
<log idref="StemCanis.prior"/>
```

# Running BEAST 2 XML files and checking ESS rates

## Step 1

Use the command `/usr/bin/time -o x beast sb2-total-evidence.x.xml`
(or `concat-total-evidence.x.xml` or `sb2-extant-only.x.xml`) when running
each XML file, in order to record the CPU time used. Again `x` is replaced
with a, b, c or d.

## Step 2

Run `sum-time.py sb2-total-evidence`, `sum-time.py concat-total-evidence` and
`sum-time.py sb2-extant-only` to get the total CPU time in hours for both
analyses.

## Step 3

Find the parameter with the lowest ESS for all analyses by running:

`loganalyser -b 0 sb2-total-evidence.combined.log > sb2-total-evidence.ess.txt`

`loganalyser -b 0 concat-total-evidence.combined.log > concat-total-evidence.ess.txt`

`loganalyser -b 0 sb2-extant-only.combined.log > sb2-extant-only.ess.txt`

The ESS per hour rate rate can then be calculated from this divided by the
times calculated in the previous step.

# Generating figures and tables from posterior distributions

## Step 1

Run `combine-logs.py sb2-total-evidence`,
`combine-logs.py concat-total-evidence` and `combine-logs.py sb2-extant-only`
to combine the four independent chains for each analysis.

## Step 2

Run `make-macroevo-table.py` to generate LaTeX code for the table of
macroevolutionary rates.

## Step 3

Run the following commands in order to generate Supplementary Figure 14:

`treeannotator -burnin 0 -heights mean sb2-total-evidence.combined.trees sb2-total-evidence.mcc.nexus`

`unfake.py sb2-total-evidence.mcc.nexus sb2-total-evidence.mcc.tree`

`Rscript sb2-sampled-ancestors.R`

The unfake script converts the leaf nodes with zero branch lengths used by
BEAST to encode sampled ancestors to degree-two internal nodes. It also
annotates those internal nodes as sampled ancestors for easier
visualization.

The R script includes manual customizations for whether sampled ancestor
labels are placed above or below their respective branches, and to flip the
orientation of individual nodes, both in order to stop labels from
overlapping in the figure. Due to the stochastic nature of MCMC, if you rerun
the MCMC chains, you will have to remove the existing customizations and add
your own as appropriate.

Run the following commands in order to generate Supplementary Figure 15:

`treeannotator -burnin 0 -heights mean concat-total-evidence.combined.trees concat-total-evidence.mcc.nexus`

`unfake.py concat-total-evidence.mcc.nexus concat-total-evidence.mcc.tree`

`Rscript concat-sampled-ancestors.R`

This script includes similar customizations to prevent labels from
overlapping, which will need to be removed in order to get the script to run
on your own MCMC chains. Run the following commands in order to generate
Supplementary Figure 16:

`treeannotator -burnin 0 -heights mean sb2-extant-only.combined.trees sb2-extant-only.mcc.nexus`

`unfake.py sb2-extant-only.mcc.nexus sb2-extant-only.mcc.tree`

`Rscript sb2-ultrametric.R`

There are no customizations needed for this figure as it is ultrametric
without sampled ancestors so labels never overlap.

## Step 4

Run `prune-trees.py sb2-total-evidence`, `prune-trees.py concat-total-evidence`
and `prune-trees.py sb2-extant-only` to generate posterior distributions of species
trees pruned to taxa with molecular sequence data.

## Step 5

Run `read-trees.py sb2-total-evidence-pruned.combined.trees sb2-total-evidence.combined.log`
to calculate the StarBEAST2 MCC tree along with node and branch expected values.

Run `read-trees.py concat-total-evidence-pruned.combined.trees concat-total-evidence.combined.log`
to do the same thing for the concatenated analysis.

Run `read-trees.py sb2-extant-only-pruned.combined.trees sb2-extant-only.combined.log`
to do the same thing for the concatenated analysis.

## Step 6

Add information on concatenation branch lengths to the StarBEAST2 MCC tree and generate Figure 3 by running the following commands:

`annotate-mcc-tree.py sb2-total-evidence-pruned.mcc.tree sb2-total-evidence-pruned.nodes.csv concat-total-evidence-pruned.branches.csv`

`fbd-length-differences.R`

## Step 7

Generate Figure 4 and Supplementary Figure 6 by running `clocks.R`

## Step 8

To generate Figure 5, run the following commands:

`ltt.py sb2-total-evidence.combined.trees`

`ltt.py concat-total-evidence.combined.trees`

`fbd-speciation-times.R`

## Step 9

To generate Supplementary Figure 7, run `msc-speciation-times.R`

## Step 10

To generate Figure 6 and Supplementary Figure 8, run
`fbd-crown-height-ellipses.R` and `msc-crown-height-ellipses.R`
respectively.

## Step 11

To generate Supplementary Figure 10, run `substitutions.R`
