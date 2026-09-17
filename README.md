#### Plotting gap compresses sequence identy across chromosomes using weighted averages to calculate sequence identity per window.  
Basic usage.  
```
bash window_identity/weighted_seq_ident_calc2.sh -threads 10 -mutation_rate 3e-8 -x asm20 -ref Zmays.fa -query Zsini.fa Bdact.fa Etef.fa Tgree.fa OkokoW.fa Eindi.fa
```

Run without options to pull up help page.
```
bash window_identity/weighted_seq_ident_calc2.sh
```

#### Synteny front-end
Anchors come from a protein liftover, matching `synLTR/module1.py`:

```
fasta_renamer_diploid -> liftover (miniprot + cd-hit + TEsorter) -> jcvi --prot
  -> anchor_builder -> gene_coords_extractor -> anchor_coord_subtracter (x2)
  -> anchor_coord_consolidator -> synmap_split (minimap2)
```

`liftover.py` maps the reference proteome onto each genome with miniprot under score
(`-outs`) and coverage (`-outc`) filters, de-duplicates the proteome with cd-hit, and strips
TE-derived peptides with a two-pass TEsorter + blastp screen. TE proteins seed false anchors
genome-wide, so this is what makes the anchors trustworthy. jcvi then anchors on **protein**
(diamond_blastp) rather than on nucleotide pseudo-CDS.

Relevant options (script defaults):

```
-peptide FILE       # reference proteome to lift over -- the biggest lever on synteny continuity
-outn N             # miniprot --outn, max alignments per protein (default 10; use 1 for haploid
                    #   or subgenome-split assemblies)
-outs FLOAT         # miniprot --outs, min score vs best hit    (default 0.95; do not lower)
-outc FLOAT         # miniprot --outc, min protein coverage     (default 0.1; little effect)
-tesorter yes|no    # two-pass TEsorter + blastp TE screen      (default yes; cost scales with proteome size)
-cscore FLOAT       # jcvi --cscore                             (default 0.99; lower for polyploids)
```

#### Getting continuous syntenic blocks
A 27-run sweep on *Poa annua* subgenomes, whose most size-discordant chromosome pair is 320 vs
98 Mb, measured how each option affects gaps between syntenic ribbons. Full results, mechanism
and a tuning recipe are in [TUNING.md](TUNING.md). In short:

- **The proteome is the dominant factor, and smaller is better.** A curated, evidence-supported
  proteome (models all three ab initio predictors agree on, plus stringtie support: 21k proteins)
  left 1.2% of chr1/chr2 in gaps between ribbons. The kitchen-sink merged proteome (901k) left
  46.9%.
  - Why: jcvi chains anchors within 20 *genes*, so every extra lifted model between two
    orthologs breaks chains.
  - cd-hit deduplication does not rescue a bloated proteome.
- **`-outn 1` and `-outs 0.99` help** on large proteomes (fewer paralog copies). `-outs 0.8` badly
  fragments blocks. `-tesorter yes` helps modestly.
- **`-outc` and `-min_block_size` do not change continuity.** `-min_block_size` only sets
  segment size under `genepair`, and does nothing under `block`.
- **`-stitch_gaps yes` closes small gaps between adjacent same-strand blocks.** On the best
  proteome it cut the chr1/chr2 gap fraction from 1.2% to 0.4%. The sweep first reported no
  effect because a bug had disabled stitching; see below.
- **Use `-partition block` for ribbon figures.**
- **Polyploids:** give each subgenome its own FASTA with `chr1`..`chrN` headers. Non-`chr`
  headers are treated as scaffolds and may be dropped.

Requires, in addition to the previous dependencies: `cd-hit`, `diamond`, and (unless
`-tesorter no`) `TEsorter`, `seqkit`, `blastp`, `makeblastdb`. The script checks for these
up front and fails fast rather than dying part-way through a long run.

#### Dotplot
Alongside jcvi's own dotplot (`{ID1}.{ID2}.pdf`, anchors placed by gene rank), Step 23 writes
`{ID1}.{ID2}.dotplot.pdf/.png` with `bin/dotplot.py`. It draws the same anchors at gene
midpoints in bp.

- **Shared bp scale.** Both axes use the same bp-per-inch scale, so one scale bar reads
  chromosome sizes on either genome, and local expansion shows as slope.
- **Labels.** Chromosomes are labelled `chr1`…`chrN`, genome names on the axes.
- **Same flips as the ribbon plot.** Query chromosomes are reversed exactly when riparian
  reverses them: same blocks, same decision code (`-flip no` for native orientation).
- **Ks colouring (`-kaks yes` only).** Dots are coloured by Ks:
  - viridis, on a square-root scale from 0 to that pair's 95th percentile;
  - higher values take the top colour (arrow on the colorbar), so one outlier cannot flatten
    the gradient;
  - anchors without a usable Ks are grey: never scored by ParaAT (e.g. an internal stop codon),
    or at or above `-ks_max`;
  - dots are drawn in random (seeded) order, because at genome scale they overlap, and a
    Ks-sorted order would let one end of the scale paint over the other.

```
-dotplot yes|no     # default yes
-flip yes|no        # default yes; applies to the riparian plot and the dotplot
```

#### Ribbon plot orientation
Assemblies often put homologous chromosomes in opposite orientations; for example, PiA chr2 is
reverse-complemented relative to PaA chr2. riparian draws each lower-track chromosome on the same
strand as its main partner above, so homologs do not draw as a twist.

- **Decision rule.** A chromosome is reversed when the bp-weighted strand of its blocks to that
  partner is mostly `-`, allowing for whether the partner is itself displayed reversed.
- **Display only.** `*.anchors.coords` is unchanged, and chromosome labels are not marked.
- **Orientation table.** `riparian.orientation.tsv` has one row per chromosome pair:
  - `upper`, `lower`: the lower chromosome and its main partner above;
  - `orientation`: `+` same strand, `-` opposite strand in the assemblies;
  - `opposite_frac`: share of their syntenic bp in the opposite orientation, so a clean `-`
    (≈1.00) can be told apart from a chromosome with internal inversions;
  - `flipped_in_plot`: whether the lower chromosome was drawn reverse-complemented.
- **Real rearrangements stay visible.** Only whole chromosomes are re-oriented, so an inversion
  *within* a chromosome still draws as a twist.
- **Native orientation.** `-flip no` in the pipeline, or `riparian.py --no-flip` standalone.
- **Older plots.** Earlier versions judged orientation from block positions, which needs at
  least two blocks. A chromosome carried by a single inverted block kept its native orientation
  and twisted, while a split chromosome was reversed. Reversed chromosomes were also labelled
  with a prime (`chr2′`). To redraw an existing run, delete `riparian.*` in its output
  directory and rerun; resume skips everything else.

#### Which ribbons get drawn

riparian discards two kinds of block before drawing: anything under 10 kb (`--min-block-len`), and
anything whose two sides differ in length by more than 5× (`--max-len-ratio`), which is the
signature of a chained-through-junk block that cannot be drawn as a ribbon. A block that fails the
ratio test is kept anyway ("rescued") if it is the only thing covering more than half of its span
on both genomes, so a genuinely expanded region is not silently dropped.

These are *block*-scale rules. Under `-partition genepair` a record is one gene-to-gene interval,
where extreme length skew is the normal signature of an expanded region rather than a chaining
artefact, so the pipeline turns both filters off (`--min-block-len 0 --max-len-ratio 0`) in that
mode. Left on, they dropped 13% of PaA on the *Poa* data — 144.7 Mb over 632 segments — and punched
holes the block plot did not have. Both modes now draw 100% of `*.anchors.coords`.

#### Synonymous divergence (Ks)
On by default. The syntenic anchors are already ortholog pairs, so the same liftover that
feeds the synteny front-end also yields in-frame CDS (`cds_walker.py`), and the two are
handed to ParaAT (mafft + pal2nal) and KaKs_Calculator:

```
liftover --outputs inframe -> {id}.cds.inframe
  + {ID1}.{ID2}.clean.anchors (already ParaAT's homolog format)
  -> ParaAT (mafft -> Epal2nal -> axt) -> KaKs_Calculator -> ks_summary -> matrix/tree/density
```

Per genome pair, the distance is the **median** Ks over gene pairs, keeping `0 <= Ks < -ks_max`.
The median (not the mean) because Ks has a long right tail from paralogous and saturated
anchors that drags a mean upward.

**Ks `NA` means Ks = 0, and it is kept.** KaKs_Calculator 3.0 prints any Ks estimate below 1e-6
as `NA` (`base.cpp`, `parseOutput`). That is a display convention, not a failure flag:

- no method ever assigns Ks the tool's internal NA sentinel;
- with zero synonymous differences (identical CDS, or nonsynonymous-only differences), YN falls
  back to Jukes–Cantor, which gives exactly 0;
- saturation yields large finite values, not `NA`.

Earlier versions of `ks_summary.py` dropped every `NA`. That discarded the most similar gene
pairs and biased recent divergences upward. On the *Poa* data it discarded 51% of PaA–PiA
anchors: the median moved from 0.016 to 0.000, and PaB–PsB from 0.0125 to 0.0061. Distant pairs
are barely affected (PaA–PaB 0.0907 → 0.0905).

A per-gene Ks of 0 still only means "below detection", roughly < 1/S for S synonymous sites, and
Ka/Ks is undefined there. That is a reason for caution about individual genes, not for dropping
them from a genome-wide distribution. When more than half of all pairs have Ks = 0, the median
is 0, and the Ks tree cannot resolve that divergence.

```
-kaks yes|no        # run the Ks branch                        (default yes)
-kaks_method NAME   # KaKs_Calculator method                   (default YN; 'ALL' unsupported)
-ks_rate RATE       # SYNONYMOUS rate for the Ks time tree     (default 1.5e-8)
-ks_max FLOAT       # saturation cutoff before the median      (default 2.0)
```

`-ks_rate` is deliberately separate from `-mutation_rate`: the latter is a genome-wide
nucleotide rate and calibrates the K2P tree, while Ks is a synonymous rate. Reusing one for
both would bias the Ks divergence times.

Outputs, mirroring the K2P set: `ks_density.pdf` (per-pair Ks density, dashed line at each
median), `ks_matrix.tsv`, `ks_matrix.nwk`(`.ape.pdf`), `ks_matrix.time.nwk`(`.ape.pdf`), plus
per-pair `{ID1}.{ID2}.kaks.tsv` and the codon alignments `{ID1}.{ID2}.axt`.

ParaAT and KaKs_Calculator are fetched and compiled on first use into `window_identity/tools/`
by `bin/setup_kaks_tools.sh` (needs `git`, `make`, `g++` that once; `perl` and `mafft` every
run). mafft is the aligner because ParaAT's `muscle` command line is muscle-v3 syntax and
breaks silently against muscle v5.

The Ks plot and tree need the R packages `optparse`, `ape` and `ggplot2`. A failing `Rscript`
does not stop the pipeline, so check that `ks_density.pdf` exists. Rerunning in the same output
directory regenerates only the plots.

#### Note....  
##### One #####
Currently, it uses UPGMA on distance matrix.   
We could also use minimum evoltuion approach for flexibility for varying mutation rates across phyla.   
Convert distance matrix to phy format.   
```
cat k2p_matrix.tsv 
	Zmays	Zsini	OkokoW	Bdact	Etef	Tgree	Eindi
Zmays	0.000000	0.142164	0.138620	0.141349	0.137039	0.144305	0.145037
Zsini	0.142164	0.000000	0.120776	0.127492	0.127879	0.125821	0.126729
OkokoW	0.138620	0.120776	0.000000	0.120185	0.121664	0.110176	0.117475
Bdact	0.141349	0.127492	0.120185	0.000000	0.126952	0.122811	0.128522
Etef	0.137039	0.127879	0.121664	0.126952	0.000000	0.121106	0.131794
Tgree	0.144305	0.125821	0.110176	0.122811	0.121106	0.000000	0.119423
Eindi	0.145037	0.126729	0.117475	0.128522	0.131794	0.119423	0.000000
```

```
cat k2p_matrix.phy
7
Zmays     0.000000 0.142164 0.138620 0.141349 0.137039 0.144305 0.145037
Zsini     0.142164 0.000000 0.120776 0.127492 0.127879 0.125821 0.126729
OkokoW    0.138620 0.120776 0.000000 0.120185 0.121664 0.110176 0.117475
Bdact     0.141349 0.127492 0.120185 0.000000 0.126952 0.122811 0.128522
Etef      0.137039 0.127879 0.121664 0.126952 0.000000 0.121106 0.131794
Tgree     0.144305 0.125821 0.110176 0.122811 0.121106 0.000000 0.119423
Eindi     0.145037 0.126729 0.117475 0.128522 0.131794 0.119423 0.000000
```

Run fastme for minimum evolution version of newick. 
```
fastme -i k2p_matrix.phy -o k2p_matrix.fastme.nwk
```

Might need to reroot to outgroup.
```
nw_reroot k2p_matrix.fastme.nwk Zmays > k2p_matrix.fastme.reroot.nwk
```

##### Two #####
~~Consider enhancing `all.anchors.coords.polished` to include inter-anchor interval.~~
Done. The synteny chain (Steps 6-9) now matches `synLTR/module1.py`:

```
anchor_builder -> gene_coords_extractor -> anchor_coord_subtracter (x2) -> anchor_coord_consolidator
```

`anchor_coord_consolidator.py --stitch-gaps` fills the interval between consecutive syntenic
blocks with a synthetic block, so that sequence reaches minimap2 in Step 10 instead of being
dropped. Stitching is suppressed where an opposite-strand block occupies the gap (an inversion)
or where two blocks are not adjacent in both genomes' orderings (a rearrangement).

Be precise about what this does and does not add. Both extractors already span the intergenic
sequence between anchors, including large gene deserts:

- `gene_coords_extractor_all4.py` (`-partition block`) cuts each jcvi cluster once, end to end.
- `gene_coords_extractor_all4_pairs.py` (`-partition genepair`) spans each consecutive anchor
  pair.

Column 4 is the block id. Merging (`-min_block_size`) is confined to one block, so a merge can
never run down a whole chromosome arm. Gaps, however, only exist *between* blocks, so stitching
works across blocks: it considers consecutive records of the same chromosome pair and strand.
A gap is stitched only if all of these hold:

- no opposite-strand block sits in it (an inversion);
- the two neighbours are consecutive in both genomes (not a rearrangement);
- its two sides differ in length by at most `-max_stitch_ratio`;
- on each genome it is no longer than `-stitch_flank_factor` × the smaller neighbouring block.
  Without this last guard, two tiny spurious blocks on non-homologous chromosomes were bridged
  across up to 75 Mb, and those fake blocks displaced real ones in the riparian plot.

Stitching judges **block extents**, not individual records: a block's extent is the span of its
records. Under `-partition block` a block is one record, so this is the same thing. Under
`genepair` a block is many gene-pair segments — only the block's outer edges bound a real gap,
and the segments bordering one are far too small and too locally scrambled to judge it by
themselves. Working on extents therefore fills the same gaps in both modes. (Where a genepair
extent stops short of the block's outer edge the flank guard can still refuse, so the genepair
run occasionally makes one stitch fewer; where both fire they agree exactly.)

(An earlier version grouped stitching by block id as well, so it could never find a gap and
silently did nothing. Re-scored with stitching restored and guarded, the 27-run sweep kept
every ranking, and no run lost coverage.)

Block orientation (column 3) is the majority vote of every anchor's relative gene strand; ties go
to anchor order. A genepair segment whose two anchors disagree takes its block's majority.
(Earlier versions used a block's first and last anchor only. That mislabelled 7–22% of block bp
on the *Poa* data, and riparian drew those blocks as inverted.)

The knobs:

```
-min_block_size N     # default 15000. Blocks with BOTH sides >= N are kept as-is; smaller
                      # blocks get merged into overlapping neighbours. Larger values merge
                      # more aggressively -> fewer, bigger blocks -> faster minimap2 in Step 10.
                      # Acts only under -partition genepair, and changes segment size, not
                      # which sequence is syntenic.
-stitch_gaps yes|no         # default yes
-max_stitch_ratio N         # default 5 (matches -max_len_ratio); 0 disables
-stitch_flank_factor N      # default 1: gap <= N x smaller neighbouring block; 0 disables
```

Note `-min_block_size` was previously hard-coded to 1000000, which merged almost every block.
The new default (15000, matching module1) is far more conservative and will produce more, smaller
blocks. If Step 10 becomes I/O-bound, raise it.

Consider flexibility to support wavefront alignment (WFA) in addition to minimap2 in `synmap_split.py`. 


