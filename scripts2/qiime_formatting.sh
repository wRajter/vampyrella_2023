#!/bin/bash

# formatting reference sequeinceds and taxonomy data for qiime taxonomy assignement


# conda activate qiime2-2022.11

############ creating a FeatureData[Sequence] artifact ###############

# each sequence sould have the accession ID as the FASTA ID, e.g.:
# >KY676659.1 Eurycea nerea isolate MO42 cytochrome oxidase subunit 1 (CO1)
# which is perfect as the GenBank ID is the FASTA ID (everything after a space is the description)

# in our dataset we need to substitute the ; character for whitespace:
# >KC174719;tax=k:Eukaryota ----> >KC174719 tax=k:Eukaryota

sed 's/;/ /' ref_dummy.fasta > ref_dummy_formated.fasta


qiime tools import \
  --input-path ref_seq_formated.fasta \
  --output-path ref_sequences.qza \
  --type 'FeatureData[Sequence]'



############# creating a FeatureData[Taxonomy] artifact ###############

# we need to make a TSV mapping file with the genbank/accession ID to the taxonomic string you want to see
# two possible source-formats:
#   TSVTaxonomyFormat which has this header: Feature ID<tab>Taxon on the first line
#   HeaderlessTSVTaxonomyFormat which has no header.
# They both are just a TSV file with IDs, tab, then taxonomy string (usually delimited by ;)

# correct format example:
# Feature ID	Taxon
# 229854  k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
# 367523	k__Bacteria; p__Bacteroidetes; c__Flavobacteriia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium; s__
# 239330	k__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfuromonadales; f__Geobacteraceae; g__Geobacter; s__
# 203525	k__Bacteria; p__OP11; c__OP11-1; o__; f__; g__; s__



# my format:
# >MF039355.1.1540_U tax=k:Eukaryota,d:Apusozoa,p:Apusomonadidae,c:Apusomonadidae_Group-4,o:Apusomonadidae_Group-4_X,f:Apusomonadidae_Group-4_XX,g:Collodictyon,s:Collodictyon_triciliatum
# >KC174719 tax=k:Eukaryota,d:Rhizaria,p:Cercozoa,c:Endomyxa,o:Vampyrellida,f:Leptophryidae,g:Leptophrys,s:Leptophrys_solarihaesitans_strain_LAK2013
# >HE609038_LV03 tax=k:Eukaryota,d:Rhizaria,p:Cercozoa,c:Endomyxa,o:Vampyrellida,f:Leptophryidae,g:*,s:Leptophrys_vorax

# we need to:
# delete the > character
# delete *
# delete tax=
# replace whitespace with tab
# replace the : character with __
# replace , to ;
# in the second column add whitespace after ;

grep '>' ref_dummy_formated.fasta > taxonomy.tsv
# delete * and >
sed 's/[>*]//g' taxonomy.tsv
# delete tax=
sed 's/tax=//g' taxonomy.tsv
# replace whitespace with tab
sed 's/\s/\t/g' taxonomy.tsv
# replace the : character with __
sed 's/:/__/g' taxonomy.tsv
# replace , to ;
sed 's/,/;/g' taxonomy.tsv
# add whitespace after ;
sed 's/;/& /g' taxonomy2.tsv
TODO: # add Feature ID<tab>Taxon as first line
grep '>' ref_seq_formated.fasta | sed 's/[>*]//g' | \
sed 's/tax=//g' | \
sed 's/\s/\t/g' | \
sed 's/:/__/g' | \
sed 's/,/;/g' | \
sed 's/ //g' | \
sed 's/;/& /g' > taxonomy_formatted.tsv



qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path taxonomy_formatted.tsv \
  --output-path ref_taxonomy.qza
