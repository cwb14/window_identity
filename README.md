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

Relevant options (defaults match module1):

```
-peptide FILE       # reference proteome to lift over
-outn N             # miniprot --outn, max alignments per protein (default 1; raise for polyploids)
-outs FLOAT         # miniprot --outs, min score vs best hit    (default 0.95)
-outc FLOAT         # miniprot --outc, min protein coverage     (default 0.9)
-tesorter yes|no    # two-pass TEsorter + blastp TE screen      (default yes; slow)
-cscore FLOAT       # jcvi --cscore                             (default 0.99; lower for polyploids)
```

Requires, in addition to the previous dependencies: `cd-hit`, `diamond`, and (unless
`-tesorter no`) `TEsorter`, `seqkit`, `blastp`, `makeblastdb`. The script checks for these
up front and fails fast rather than dying part-way through a long run.

#### Synonymous divergence (Ks)
On by default. The syntenic anchors are already ortholog pairs, so the same liftover that
feeds the synteny front-end also yields in-frame CDS (`cds_walker.py`), and the two are
handed to ParaAT (mafft + pal2nal) and KaKs_Calculator:

```
liftover --outputs inframe -> {id}.cds.inframe
  + {ID1}.{ID2}.clean.anchors (already ParaAT's homolog format)
  -> ParaAT (mafft -> Epal2nal -> axt) -> KaKs_Calculator -> ks_summary -> matrix/tree/density
```

Per genome pair, the distance is the **median** Ks over gene pairs, keeping `0 < Ks < -ks_max`.
The median (not the mean) because Ks has a long right tail from paralogous and saturated
anchors that drags a mean upward.

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

Be precise about what this does and does not add. `gene_coords_extractor_all4.py` emits one span
per *consecutive anchor pair within a jcvi cluster*, so a span already runs from gene i to gene
i+1 and already covers the intergenic sequence between anchors -- including large gene deserts.
A positive gap therefore only exists **between separate jcvi clusters (`###`)**, and that is the
only thing `--stitch-gaps` fills. On closely collinear genomes it may add nothing at all. Its
value shows up when synteny is fragmented into multiple collinear clusters on the same chromosome
pair and strand.

Two knobs control it:

```
-min_block_size N     # default 15000. Blocks with BOTH sides >= N are kept as-is; smaller
                      # blocks get merged into overlapping neighbours. Larger values merge
                      # more aggressively -> fewer, bigger blocks -> faster minimap2 in Step 10.
-stitch_gaps yes|no   # default yes.
```

Note `-min_block_size` was previously hard-coded to 1000000, which merged almost every block.
The new default (15000, matching module1) is far more conservative and will produce more, smaller
blocks. If Step 10 becomes I/O-bound, raise it.

Consider flexibility to support wavefront alignment (WFA) in addition to minimap2 in `synmap_split.py`. 


