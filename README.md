# _Gardnerella vaginalis_
_G. vaginalis_ is the primary pathogen involved in Bacterial Vaginosis (BV), the most common infection for reproductive aged women. G. vaginalis produces a sticky biofilm constructed of uncharacterized exopolysacchrides. This biofilm provides a home for other anaerobes to colonize the vagina, and further increases the risk of BV recurrence due to decreased antibiotic susceptibility.

By constructing a _G. vaginalis_ metabolic model we aim to understand essential biofilm metabolism, and then disrupt those newly characterized metabolic pathways to recoupe antibiotic efficacy. The following documentation outlines model construction and curration from start to finish, along with experimental data integration.

## Base Model Construction
Build metabolic model using [_Gardnerella vaginalis_ proteome](https://www.uniprot.org/proteomes/UP000035707) and [CarveMe](https://carveme.readthedocs.io/en/latest/index.html)

Did not require [gapfilling](https://github.com/jotech/gapseq) for a functional model.

### [_in Silico_ Media](https://github.com/lrd3uu/gardnerella_vaginalis/blob/main/insilico_media.py)
Defined media conditions will always be easier to simulate in silico, as compared to enriched media. I have included both enriched, KSFM & NYCIII, and defined media ([SVM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC265277/)) conditions within the provided code. [NYCIII media](https://www.atcc.org/~/media/FA8074C3B4B9450899EE2542D6AD7116.ashx) is the suggested growing conditions for _G. vaginalis_ and KSFM is the enriched media used for [vaginal epithelial cell](https://www.atcc.org/products/crl-2616) growth - a cell line later used for co-culture system to simulate host-pathogen interactions _in vitro_ 

Steps:
1. Identify media metabolites in the model; add missing extracellular metabolites if needed
2. Determine which metabolites are missing an exchange reaction; add missing exchange reactions
3. Locate transport reactions
4. Add missing transport reactions for metabolites that have both intracellular and extracellular metabolites
5. Simulate growth - Gapfill if necessary  
