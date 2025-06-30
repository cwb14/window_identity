#### Plotting gap compresses sequence identy across chromosomes using weighted averages to calculate sequence identity per window.  
Basic usage.  
```
bash window_identity/weighted_seq_ident_calc2.sh -threads 10 -mutation_rate 3e-8 -x asm20 -ref Zmays.fa -query Zsini.fa Bdact.fa Etef.fa Tgree.fa OkokoW.fa Eindi.fa
```

Run without options to pull up help page.
```
bash window_identity/weighted_seq_ident_calc2.sh
```


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
Consider enhancing `all.anchors.coords.polished` to include inter-anchor interval.  
Consider flexibility to support wavefront alignment (WFA) in addition to minimap2 in `synmap_split.py`. 


