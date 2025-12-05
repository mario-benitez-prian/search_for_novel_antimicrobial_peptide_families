#!/usr/bin/env python3
# cleaning_newick.py
#
# Lee la variable NEWICK (string), elimina etiquetas de nodos internos
# (todo lo que sigue a un ')' hasta la siguiente coma, ')' o ';'),
# y escribe el resultado en pantalla y en archivo cleaned_newick.nwk.
#
# Variables definidas al principio (sin argumentos).

NEWICK = """((((((((Anopheles gambiae)Anopheles)Culicidae,((Drosophila melanogaster)Drosophila)Drosophilidae)Diptera,(((Bombyx mori)Bombyx)Bombycidae)Lepidoptera,(((Tribolium castaneum)Tribolium)Tenebrionidae)Coleoptera,(((Apis mellifera)Apis)Apidae)Hymenoptera)Insecta,((((Daphnia pulex)Daphnia)Daphniidae)Diplostraca)Branchiopoda)Arthropoda,(((((Strongylocentrotus purpuratus)Strongylocentrotus)Strongylocentrotidae)Camarodonta)Echinoidea)Echinodermata,(((((Crassostrea gigas)Crassostrea)Ostreidae)Ostreida)Bivalvia,((((Lottia gigantea)Lottia)Lottiidae)LottiidaF+,(((Aplysia californica)Aplysia)Aplysiidae)Aplysiida)Gastropoda,((((Octopus bimaculoides)Octopus)Octopodidae)Octopoda)Cephalopoda)Mollusca,(((((Homo sapiens)Homo)Hominidae)Primates,(((Rattus norvegicus)Rattus,(Mus musculus)Mus)Muridae)Rodentia)Mammalia,((((Gallus gallus)Gallus)Phasianidae)Galliformes)Aves,((((Takifugu rubripes)Takifugu)Tetraodontidae)Tetraodontiformes,(((Danio rerio)Danio)Danionidae)Cypriniformes,(((Gasterosteus aculeatus)Gasterosteus)Gasterosteidae)Perciformes)Actinopteri,((((Anolis carolinensis)Anolis)Dactyloidae)Squamata)Lepidosauria,((((Latimeria chalumnae)Latimeria)Coelacanthidae)Coelacanthiformes)CoelacanthiformeO+,((((Callorhinchus milii)Callorhinchus)Callorhinchidae)Chimaeriformes)Chondrichthyes,((((Xenopus tropicalis)Xenopus)Pipidae)Anura)Amphibia,((((Petromyzon marinus)Petromyzon)Petromyzontidae)Petromyzontiformes)Hyperoartia,((((Ciona intestinalis)Ciona)Cionidae)Phlebobranchia)Ascidiacea)Chordata,(((((Nematostella vectensis)Nematostella)Edwardsiidae)Actiniaria)Anthozoa)Cnidaria,(((((Caenorhabditis elegans)Caenorhabditis)Rhabditidae)Rhabditida)Chromadorea)Nematoda)Metazoa,((((((Chlamydomonas reinhardtii)Chlamydomonas)Chlamydomonadaceae)Chlamydomonadales)Chlorophyceae)Chlorophyta,(((((Selaginella moellendorffii)Selaginella)Selaginellaceae)Selaginellales)Lycopodiopsida,((((Zea mays)Zea,(Sorghum bicolor)Sorghum,(Oryza sativa)Oryza,(Brachypodium distachyon)Brachypodium)Poaceae)Poales,(((Medicago truncatula)Medicago,(Glycine max)Glycine)Fabaceae)Fabales,(((Malus domestica)Malus)Rosaceae)Rosales,(((Arabidopsis thaliana)Arabidopsis)Brassicaceae)Brassicales,(((Populus trichocarpa)Populus)Salicaceae)Malpighiales,(((Cucumis sativus)Cucumis)Cucurbitaceae)Cucurbitales,(((Vitis vinifera)Vitis)Vitaceae)Vitales,(((Theobroma cacao)Theobroma)Malvaceae)Malvales,(((Solanum tuberosum)Solanum)Solanaceae)Solanales)Magnoliopsida,((((Physcomitrium patens)Physcomitrium)Funariaceae)Funariales)Bryopsida)Streptophyta)Viridiplantae,((((((Candida albicans)Candida)Debaryomycetaceae,((Saccharomyces cerevisiae)Saccharomyces)Saccharomycetaceae)Saccharomycetales)Saccharomycetes,((((Schizosaccharomyces pombe)Schizosaccharomyces)Schizosaccharomycetaceae)Schizosaccharomycetales)Schizosaccharomycetes,((((Aspergillus fumigatus,Aspergillus nidulans)Aspergillus)Aspergillaceae)Eurotiales)Eurotiomycetes,((((Botrytis cinerea)Botrytis)Sclerotiniaceae)Helotiales)Leotiomycetes,((((Neurospora crassa)Neurospora)Sordariaceae)Sordariales)Sordariomycetes)Ascomycota,(((((Coprinopsis cinerea)Coprinopsis)Psathyrellaceae)Agaricales)Agaricomycetes,((((Ustilago maydis)Ustilago)Ustilaginaceae)Ustilaginales)Ustilaginomycetes,((((Cryptococcus neoformans)Cryptococcus)Cryptococcaceae)Tremellales)Tremellomycetes)Basidiomycota)Fungi,((((((Thalassiosira pseudonana)Thalassiosira)Thalassiosiraceae)Thalassiosirales)Coscinodiscophyceae,((((Phaeodactylum tricornutum)Phaeodactylum)Phaeodactylaceae)Naviculales)Bacillariophyceae)Bacillariophyta)BacillariophytP+,((((((Tetrahymena thermophila)Tetrahymena)Tetrahymenidae)Hymenostomatida,(((Paramecium tetraurelia)Paramecium)Parameciidae)Peniculida)Oligohymenophorea)Ciliophora)CiliophorP+,((((((Plasmodium falciparum)Plasmodium)Plasmodiidae)Haemosporida)Aconoidasida,((((Toxoplasma gondii)Toxoplasma)Sarcocystidae)Eucoccidiorida)Conoidasida)Apicomplexa)ApicomplexP+,((((((Dictyostelium discoideum)Dictyostelium)Dictyosteliaceae)Dictyosteliales)Eumycetozoa)Evosea)EvoseP+,((((((Giardia lamblia)Giardia)Hexamitidae)Diplomonadida)DiplomonadidO+)Fornicata)FornicatP+,((((((Cyanidioschyzon merolae)Cyanidioschyzon)Cyanidiaceae)Cyanidiales)Bangiophyceae)Rhodophyta)RhodophytP+,((((((Trypanosoma brucei)Trypanosoma,(Leishmania major)Leishmania)Trypanosomatidae)Trypanosomatida)Kinetoplastea)Euglenozoa)EuglenozoP+)Eukaryota,(((((((Ignicoccus hospitalis)Ignicoccus)Desulfurococcaceae)Desulfurococcales,(((Sulfolobus solfataricus)Sulfolobus,(Sulfolobus acidocaldarius)Sulfolobus)Sulfolobaceae)Sulfolobales)Thermoprotei)Crenarchaeota,(((((Nitrosopumilus maritimus)Nitrosopumilus)Nitrosopumilaceae)Nitrosopumilales)NitrosopumilaleO+)Thaumarchaeota,(((((Methanococcus maripaludis)Methanococcus)Methanococcaceae,((Methanocaldococcus jannaschii)Methanocaldococcus)Methanocaldococcaceae)Methanococcales)Methanococci,((((Pyrococcus furiosus)Pyrococcus,(Thermococcus kodakarensis)Thermococcus)Thermococcaceae)Thermococcales)Thermococci,((((Haloferax volcanii)Haloferax)Haloferacaceae)Haloferacales,(((Halobacterium salinarum)Halobacterium)Halobacteriaceae)Halobacteriales)Halobacteria,((((Archaeoglobus fulgidus)Archaeoglobus)Archaeog...
"""
# (p.ej. pega aquí todo tu Newick completo)

# Output file
OUTFILE = "cleaned_newick.nwk"

# --- processing: remove internal node labels that come after a closing parenthesis
import re

s = NEWICK.strip()

# Pattern: a ')' followed by some non-delimiter characters (the internal label),
# only when immediately followed by a comma, another ), or semicolon.
pattern = re.compile(r'\)([^,);\s]+)(?=[,);\s])')

# Apply repeatedly (some files can have nested labels)
prev = None
count = 0
while prev != s:
    prev = s
    s, nsub = pattern.subn(')', s)
    count += nsub

# Also remove stray tokens like multiple spaces produced
s = re.sub(r'\s+', ' ', s)

# Final cleanup: remove space before comma/paren/semicolon
s = re.sub(r'\s+([,);\]])', r'\1', s)

# Ensure semicolon termination
if not s.strip().endswith(';'):
    s = s.strip() + ';'

# Write out
with open(OUTFILE, 'w', encoding='utf-8') as out:
    out.write(s + "\n")

print(f"Replacements performed: {count}")
print("Cleaned Newick (first 400 chars):")
print(s[:400])
print(f"\nCleaned file written to {OUTFILE}")
