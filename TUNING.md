# Tuning the synteny front-end for continuous syntenic blocks

What each synteny-related parameter of `weighted_seq_ident_calc2.sh` actually does to the
result, measured in a 27-run sweep. Read this before choosing `-peptide`; it matters far more
than any other option.

## Recommendations

| option | recommendation | effect in testing |
|---|---|---|
| `-peptide` | A **small, high-confidence, evidence-supported** proteome (roughly one model per locus). Bigger is worse. | **Dominant.** Gap fraction ranged 0.012–0.469 |
| `-outn` | `1` for haploid assemblies and for polyploids split into one FASTA per subgenome. Raise it only if homoeologs share one FASTA. | Strong with large proteomes; none with a curated one |
| `-outs` | Keep at `0.95` (default) or `0.99`. **Never lower it.** | 0.99 slightly better; 0.8 much worse |
| `-tesorter` | `yes` (default). | Modest gain; cheap with a curated proteome |
| `-outc` | Leave at the default. | None (0.1, 0.5 and 0.9 were indistinguishable) |
| `-min_block_size` | Sets segment size under `-partition genepair`; no role under `block`. | **No effect on continuity** in either mode |
| `-stitch_gaps` | `yes` (default). | Closes small gaps between adjacent same-strand blocks: best set 0.012 → 0.004 |
| `-partition` | `block` for ribbon (riparian) figures. | `genepair` segments are small, and riparian drops many of them |

## Why the proteome dominates

The anchors come from `jcvi.compara.catalog ortholog`:

1. `diamond blastp` between the two lifted-over proteomes.
2. **C-score filter.** `--cscore 0.99` keeps only near-reciprocal-best hits.
3. **Chaining.** An anchor joins a cluster if it lies within `--dist 20` **genes** of it in both
   genomes. Clusters need at least 4 anchors.

Distance in step 3 is counted in *gene ranks*, not base pairs. Every lifted gene model in
`{id}.bed` takes a rank, whether or not it has an ortholog. So each extra model between two true
orthologs (a paralog copy, a redundant isoform-level model, a fragmentary or TE-derived
prediction) pushes them apart. Extra models also make reciprocal-best hits rarer. Once 20 ranks
pass without an anchor, the chain breaks.

With `-partition block`, each jcvi cluster is drawn as one ribbon, so every break is a visible
gap. Gene-poor, repeat-rich regions such as expanded pericentromeres break first.

**Deduplication does not fix this.** cd-hit (90% identity, run inside `liftover.py`) cut a
456k-protein merged proteome to only 122k. The winning set lifted 14k. The harmful extra models
are distinct low-confidence predictions, not near-duplicates, so they have to be removed by
choosing the proteome on evidence.

`-outn` and `-outs` work through the same mechanism. Allowing more secondary alignments per
protein adds paralog copies to the gene order.

## The sweep

**Data.** *Poa annua* is an allotetraploid. It was split into two subgenome FASTAs: PaA
(1.12 Gb) and PaB (0.66 Gb). The target was the most size-discordant homoeologous pairs:
chr1 (320 vs 98 Mb) and chr2 (266 vs 146 Mb). Their pericentromeres are greatly expanded in A,
and that is where synteny usually fragments.

**Proteome.** The proteome was a merged annotation whose headers record which predictors support
each model: ab initio (annevo, tiberius, helixer), short-read stringtie, and long-read
isoquant/isoform collapsing.

**Fixed settings:** `-threads 6 -processes 5 -align no -kaks no -partition block -cscore 0.99`.

**Metric.** Blocks were filtered exactly as `riparian.py` filters them: at least 10 kb, and at
most 5-fold length skew unless the block fills otherwise uncovered sequence. The union of each
homoeologous pair's blocks was then taken on both chromosomes. **Gap fraction** is the uncovered
internal sequence, averaged over PaA/PaB chr1 and chr2 (lower is better). Genome-wide coverage
averages all 7 homoeologous pairs.

### Proteome

All rows use script defaults (`-outn 10 -outs 0.95 -outc 0.1 -tesorter yes`) unless noted.
"Lifted" is the protein count after cd-hit and the TE screen, i.e. what miniprot aligned.
Runtime is wall clock with `-threads 6 -processes 5`.

| proteome | input | lifted | gap fraction | chr1 cov A / B | chr2 cov A / B | genome cov | runtime |
|---|---|---|---|---|---|---|---|
| **all 3 ab initio tools ∩ stringtie**, `-outn 1` | 21k | 14k | **0.012** | 98.6 / 96.9 | 99.7 / 99.5 | **97.7** | 4 min |
| all 3 ab initio tools ∩ stringtie | 21k | 14k | 0.012 | 98.5 / 96.9 | 99.7 / 99.5 | 97.7 | 4 min |
| ≥2 ab initio tools ∩ stringtie | 31k | 19k | 0.041 | 98.8 / 94.0 | 98.8 / 91.6 | 97.1 | 6 min |
| annevo ∩ stringtie | 35k | 20k | 0.043 | 99.3 / 95.2 | 98.4 / 89.3 | 97.1 | 6 min |
| ≥1 ab initio tool ∩ stringtie | 37k | 21k | 0.048 | 99.3 / 91.3 | 98.0 / 91.5 | 96.7 | 8 min |
| ≥1 ab initio tool ∩ any transcript evidence | 46k | 25k | 0.074 | 98.6 / 92.5 | 95.5 / 83.5 | 95.6 | 8 min |
| ≥2 ab initio tools (no transcript requirement) | 47k | 28k | 0.091 | 89.8 / 87.0 | 96.2 / 90.1 | 95.4 | 8 min |
| annevo only | 63k | 35k | 0.161 | 82.7 / 85.8 | 82.4 / 84.0 | 91.4 | 11 min |
| stringtie only | 74k | 36k | 0.194 | 74.0 / 87.7 | 77.3 / 83.1 | 89.1 | 11 min |
| ≥1 ab initio tool | 128k | 59k | 0.286 | 63.0 / 82.6 | 64.5 / 75.2 | 81.4 | 23 min |
| progenitor species, ≥1 ab initio tool | 134k | 63k | 0.316 | 57.0 / 80.7 | 61.1 / 74.4 | 79.3 | 27 min |
| everything except long-read-only models | 165k | 70k | 0.331 | 58.0 / 80.3 | 57.0 / 71.6 | 79.2 | 27 min |
| progenitor species, all models | 445k | 116k | 0.390 | 48.2 / 75.1 | 52.4 / 67.7 | 73.2 | 45 min |
| all models | 456k | 117k | 0.391 | 49.0 / 75.8 | 51.2 / 67.2 | 74.5 | 43 min |
| all models, all 3 species | 901k | 169k | 0.469 | 38.2 / 67.9 | 44.1 / 61.9 | 66.1 | 81 min |

What this shows:

- **More proteins meant more gaps, almost monotonically.** The kitchen sink was the worst option
  and the slowest.
- **Requiring agreement helps more than any single source.** Intersections beat either source
  alone (annevo ∩ stringtie 0.043, annevo 0.161, stringtie 0.194), and stricter ab initio
  consensus kept improving the result.
- **Transcript support matters at every consensus level.** Dropping it roughly doubled gaps
  (0.041 with it vs 0.091 without, at the ≥2-tools level).
- **A related species' proteome was no better than the genome's own** (progenitors 0.316 vs
  0.286).

### Liftover and consolidation parameters

These were tested one at a time on the 128k "≥1 ab initio tool" proteome, which scores 0.286 at
defaults.

| change | gap fraction | reading |
|---|---|---|
| `-outn 1` | **0.207** | Fewer paralog copies; strongest non-proteome lever |
| `-outs 0.99` | 0.238 | Same mechanism, weaker |
| `-outc 0.9` | 0.280 | No effect |
| `-outc 0.5` | 0.284 | No effect |
| `-tesorter no` | 0.323 | TE-derived models hurt |
| `-outs 0.8` | 0.413 | Many more paralog copies; badly fragmented |

- **`-outn`/`-outs` on a curated proteome.** They barely matter there: `-outn 1` scored 0.0118
  vs 0.0120, and `-outs 0.99` on the 37k set scored 0.053 vs 0.048. A curated proteome has few
  paralog copies left to remove.
- **`-min_block_size` does nothing with `-partition block`.** Rerunning Step 8 with 0, 15 kb,
  1 Mb and 50 Mb gave byte-identical `*.anchors.coords`. Merging is confined to one block, and in
  block mode each block is a single record, so there is nothing to merge.
- **Under `-partition genepair`, `-min_block_size` only sets granularity.** On PaA–PaB, values
  of 0, 15 kb, 100 kb and 1 Mb gave 12,575, 7,879, 3,217 and 768 segments (median 22 kb, 51 kb,
  172 kb, 1.1 Mb) over the same 1,101.830 Mb. Raise it for fewer, larger Step 10 alignment jobs.
  It does not change which sequence is syntenic.
- **`-stitch_gaps` was broken during the sweep, and has since been fixed.** Stitching was
  grouped by block id, so it never saw the between-block gaps it exists to fill. Every gap
  fraction in this document was measured with stitching off.
  - **Restored:** stitching now spans blocks, with a new guard that a gap must be no longer than
    the smaller neighbouring block on each genome (`-stitch_flank_factor 1`).
  - **Re-scored sweep:** every ranking kept, and no run lost coverage. The best set dropped from
    0.012 to 0.004, i.e. the two small chr1 gaps (2.2 and 2.0 Mb) next to the pericentromere
    closed.
  - **Why the new guard matters:** without it, sparse spurious blocks on non-homologous
    chromosomes were stitched across up to 75 Mb (313 of 1,666 stitches). Those fake blocks
    displaced real blocks in 8 of 27 runs. With it, 800 stitches remain: 2 non-homologous (both
    under 0.3 Mb) and none over 10 Mb.
  - **Genepair:** block edges are locally scrambled, so the guards usually refuse; the final
    genepair runs stitched nothing.

### Is a sparse proteome hiding real structure?

A smaller proteome gives fewer anchors (2,068 on chr1 vs 12,430 with all models), so it was
checked for over-merging.

- **Same path.** Anchor dotplots of the best, 37k and annevo-only runs trace the same PaA–PaB
  path.
- **Real inversions persist.** A 17 Mb inversion on chr1 and the inverted pericentromeric
  segments appear in all three.
- **Best collinearity.** Within-block collinearity (anchor-weighted mean |Spearman ρ|,
  chr1 + chr2) was highest for the best set: 0.998, with 1 block below |ρ| 0.9. The 37k set
  scored 0.983 with 7 such blocks, and annevo-only 0.985 with 24. The small "inversions" the
  looser sets report inside pericentromeres are mostly scrambled noise (ρ ≈ 0).
- **Enough anchors for Ks.** It still yielded about 12.7k genome-wide anchors, and 98.6% of them
  were scored.

**Trade-off.** A curated proteome gives a clean syntenic backbone, but it omits genes absent from
the high-confidence set. Use a broader proteome when the question is about those genes, not about
block continuity.

## How to tune on your own data

The case study is one dataset. The mechanism (gene-rank chaining) is general, but the best
proteome depends on your annotation.

1. **Pick your hardest pair.** Choose your most size- or repeat-discordant homologous chromosome
   pair.
2. **Build 3–5 proteome variants of decreasing size.** Tighten by evidence agreement, e.g.
   ab initio ∩ transcript-supported, then stricter consensus. Use one model per locus where you
   can.
3. **Run each variant with `-align no -kaks no -partition block`.** With a curated proteome this
   takes minutes, and the variants run in parallel.
4. **Compare `riparian.pdf` for that pair.** Count the ribbon gaps, and check whether known
   inversions remain.
5. **Stop at diminishing returns.** Only then add `-outn 1` if the assembly is haploid or split
   by subgenome, and run the full job with `-kaks yes`.

## Practical notes

- **Polyploid inputs.** Give each subgenome its own FASTA. The FASTA filename becomes the genome
  ID. The renamer (`fasta_renamer_diploid.py`) only treats `Chr`/`chr`-style headers as
  chromosomes. Anything else (e.g. `Pa1A`) is handled as a scaffold, renamed `{id}_scaN`, and
  dropped when that name exceeds 13 characters. Rename headers to `chr1`…`chrN` first.
- **Changing only `-partition`, `-min_block_size` or `-stitch_gaps`.** Reuse a finished run instead of
  starting over. Copy its directory and delete:
  - `*.anchors.raw.coords`
  - `*.anchors.coords.polished`, `*.anchors.coords.polished2`
  - `*.anchors.coords`
  - `all.anchors.coords.polished`
  - `riparian.{pdf,png,html}`

  Then rerun with the new option. The resume guards reuse the liftover, anchors and Ks, which
  are all partition-independent. This takes under a minute.
- **Ks plots need R packages.** `ks_density_plotter.R` and `upgma.R` need `optparse`, `ape` and
  `ggplot2`. A failing `Rscript` does **not** stop the pipeline, so check that `ks_density.pdf`
  exists. Rerunning in the same directory regenerates only the plots.
- **Block orientation** (column 3 of `*.anchors.coords`) is the majority vote of every anchor's
  relative gene strand. Earlier versions used only a block's first and last anchor. On the Poa
  runs that mislabelled 7–22% of block bp, including perfectly collinear whole-chromosome
  blocks, and riparian drew them as inverted. For outputs from an older version, rebuild the
  coords as in the partition note above.
