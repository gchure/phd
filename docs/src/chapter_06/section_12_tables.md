## *E. coli* Primer and Strain List 
Here, we provide additional details about the genotypes of the strains
used, as well as the primer sequences used to generate them. *E. coli*
strains were derived from K12 MG1655. For those containing $R=22$, we
used strain HG104 which additionally has the *lacYZA* operon deleted
(positions 360,483 to 365,579), but still contains the native *lacI*
locus. All other strains used strain HG105, where both the *lacYZA* and
*lacI* operons have both been deleted (positions 360,483 to 366,637).

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;All 25x+11-yfp expression constructs were integrated at the *galK* locus
(between positions 1,504,078 and 1,505,112) while the 3\*1x-lacI
constructs were integrated at the *ybcN* locus (between positions
1,287,628 and 1,288,047). Integration was performed with $\lambda$ Red
recombineering [@sharan2009] as described in @garcia2011 using the
primers listed in Table 6.3. We follow the notation of Lutz and Bujard
[@lutz1997] for the nomenclature of the different constructs used.
Specifically, the first number refers to the antibiotic resistance
cassette that is present for selection (2 = kanamycin, 3 =
chloramphenicol, and 4 = spectinomycin) and the second number refers to
the promoter used to drive expression of either YFP or LacI (1 =
$P_{LtetO-1}$, and 5 = *lacUV5*). Note that in 25x+11-yfp, x refers to
the LacI operator used, which is centered at +11 (or alternatively,
begins at the transcription start site). For the different LacI
constructs, 3\*1x-lacI, x refers to the different ribosomal binding site
modifications that provide different repressor copy numbers and follows
from @garcia2011. The asterisk refers to the presence of FLP recombinase
sites flanking the chloramphenicol resistance gene that can be used to
lose this resistance. However, we maintained the resistance gene in our
constructs. A summary of the final genotypes of each strain is listed in
Table 6.4. In addition, each strain also contained the
plasmid pZS4\*1-mCherry and provided constitutive expression of the
mCherry fluorescent protein. This pZS plasmid is a low copy (SC101
origin of replication) where, like with 3\*1x-lacI, mCherry is driven by
a $P_{LtetO-1}$ promoter.

| **Primer** | **Sequence (5' $\rightarrow$ 3')** | **Notes**|
|:--|:--|:--|
| pZSForwSeq2 | TTCCCAACC TTACCAGAGG GC | Forward sequencing primer for 3\*1x-lacI |
|251 F | CCTTTCGTCT TCACCTCGA| Forward sequencing primer for 25x+11-YFP|
| YFP1 | ACTAGCAACAC CAGAACAGCCC | Reverse sequencing primer for 3\*1x-lacI and 25x+11-YFP | 
| HG 6.1 (*galK*) | gtttgcgcgc agtcagcgat atccattttc gcgaatccg gagtgtaag aaACTAGCAAC ACCAGAACA GCC | Reverse primer for 25x+11-YFP integration in to the *galK* locus (lowercase). |
| HG 6.3 (*galK*) | ttcatattgt tcagcgacag cttgctgtac ggcaggcac cagctcttc cgGGCTAATGC ACCCAGTAA GG | Forward integration primer for 25x+11-YFP with homology to the *galK* locus (lowercase.|
| HG11.1 (*ybcN*) | acctctgcgg aggggaagcg tgaacctctc acaagacgg catcaaatt acACTAGCAAC ACCAGAACA GCC| Reverse integration primer for 3\*1x-lacI with homology to the *ybcN* locus (lowercase). |
| HG11.3 (*ybcN*)| ctgtagatgtg tccgttcatg acacgaataa gcggtgtag ccattacgc cGGCTAATGCA CCCAGTAAG G| Forward integration primer for 3\*1x-lacI with homology to the *ybcN* locus (lowercase).|
| ybcN-control-upstream-1 | AGCGTTTGA CCTCTGCGGA | Sequencing primer to verify integration | 
| ybcN-control-downstream-1| GCTCAGGTT TACGCTTAC GACG | Sequencing primer to verify integration | 
: Primers used in this work.


|**Strain** | **Genotype** |
|:--|:--|
| O1, $R=0$ | HG105::*galK*<>25O1+11-YFP |
| O1, $R=22$ | HG104::*galK*<>25O1+11-YFP |
| O1, $R=60$ | HG105::*galK*<>25O1+11-YFP, *ybcN*<>3\*1RBS1147-lacI |
| O1, $R=124$ | HG105::*galK*<>25O1+11-YFP, *ybcN*<>3\*1RBS446-lacI |
| O1, $R=260$ | HG105::*galK*<>25O1+11-YFP, *ybcN*<>3\*1RBS1027-lacI |
| O1, $R=1220$ | HG105::*galK*<>25O1+11-YFP, *ybcN*<>3\*1RBS1-lacI |
| O1, $R=1740$ | HG105::*galK*<>25O1+11-YFP, *ybcN*<>3\*1RBS1L-lacI |
| O2, $R=0$ | HG105::*galK*<>25O2+11-YFP |
| O2, $R=22$ | HG104::*galK*<>25O2+11-YFP |
| O2, $R=60$ | HG105::*galK*<>25O2+11-YFP, *ybcN*<>3\*1RBS1147-lacI |
| O2, $R=124$ | HG105::*galK*<>25O2+11-YFP, *ybcN*<>3\*1RBS446-lacI |
| O2, $R=260$ | HG105::*galK*<>25O2+11-YFP, *ybcN*<>3\*1RBS1027-lacI |
| O2, $R=1220$ | HG105::*galK*<>25O2+11-YFP, *ybcN*<>3\*1RBS1-lacI |
| O2, $R=1740$ | HG105::*galK*<>25O2+11-YFP, *ybcN*<>3\*1RBS1L-lacI |
| O3, $R=0$ | HG105::*galK*<>25O3+11-YFP |
| O3, $R=22$ | HG104::*galK*<>25O3+11-YFP |
| O3, $R=60$ | HG105::*galK*<>25O3+11-YFP, *ybcN*<>3\*1RBS1147-lacI |
| O3, $R=124$ | HG105::*galK*<>25O3+11-YFP, *ybcN*<>3\*1RBS446-lacI |
| O3, $R=260$ | HG105::*galK*<>25O3+11-YFP, *ybcN*<>3\*1RBS1027-lacI |
| O3, $R=1220$ | HG105::*galK*<>25O3+11-YFP, *ybcN*<>3\*1RBS1-lacI |
| O3, $R=1740$ | HG105::*galK*<>25O3+11-YFP, *ybcN*<>3\*1RBS1L-lacI |
| Oid, $R=22$ | HG104::*galK*<>25Oid+11-YFP |
| Oid, $R=60$ | HG105::*galK*<>25Oid+11-YFP, *ybcN*<>3\*1RBS1147-lacI |
| Oid, $R=124$ | HG105::*galK*<>25Oid+11-YFP, *ybcN*<>3\*1RBS446-lacI |
: *E. coli* strains used in this work.


