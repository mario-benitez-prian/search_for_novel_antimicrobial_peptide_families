from io import StringIO
from Bio import Phylo

# Ruta al archivo Newick original
input_newick = "iphylo_tree_with_lowbar.txt"
output_newick = "iphylo_tree_with_lowbar_branch_length.txt"

# Longitud fija por rama
FIXED_LENGTH = 0.5

# Leer el árbol
tree = Phylo.read(input_newick, "newick")

# Recorrer nodos y asignar branch_length si no existe
for clade in tree.find_clades():
    # saltar la raíz si lo deseas (root.branch_length puede ser None)
    if clade is tree.root:
        continue
    clade.branch_length = FIXED_LENGTH

# Escribir el árbol modificado
Phylo.write(tree, output_newick, "newick")

print("Árbol modificado guardado en:", output_newick)
