# ===============================================
# Venn diagram for three sequence ID files
# ===============================================


import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from matplotlib_venn import venn2, venn3
from pathlib import Path

"""

# ---- Input files ----
file_a = "../data/novel_seqs/mmseqs_novel_seqs_ids.txt"
file_b = "../data/novel_seqs/psi-cd-hit_novel_seqs_ids.txt"
file_c = "../data/novel_seqs/interpro_novel_seqs_ids.txt"
file_d = "../data/novel_seqs/predicted_seqs_ids.txt"

# ---- Output files ----
results_dir = Path("../results/venn_diagram/")
results_dir.mkdir(parents=True, exist_ok=True)

consensus_file_triple = results_dir / "triple_consensus.txt"
venn_diagram_file_triple = results_dir / "triple_venn_diagram.png"

consensus_file_double = results_dir / "double_consensus.txt"
venn_diagram_file_double = results_dir / "double_venn_diagram.png"

consensus_file_all = results_dir / "all_consensus.txt"
venn_diagram_file_all = results_dir / "all_venn_diagram.png"


# ---- Load sequence IDs ----
def load_ids(filename):
    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())

set_a = load_ids(file_a)
set_b = load_ids(file_b)
set_c = load_ids(file_c)
set_d = load_ids(file_d)

# ---- Compute intersections ----
inter_all = set_a & set_b & set_c

# ---- Save triple intersection ----
with open(consensus_file_triple, "w") as out:
    out.write("\n".join(sorted(inter_all)))

print(f"{len(inter_all)} sequences are shared by all three files.")
print(f"Saved as: {consensus_file_triple}")

#
#
# ---- Plot Venn diagram for three databases----
#
#

plt.figure(figsize=(6, 6))

v = venn3(
    subsets=(set_a, set_b, set_c),
    set_labels=("Novel AMPs (mmseqs)", "Novel AMPs (psi-cd-hit)", "Novel AMPs (interpro)"),
)

# ---- Optional aesthetic improvements ----
for text in v.set_labels:
    text.set_fontsize(14)
for text in v.subset_labels:
    if text:
        text.set_fontsize(12)

# High-resolution output for publication
plt.savefig(venn_diagram_file_triple, dpi=400, bbox_inches='tight')
print(f"Saved as: {venn_diagram_file_triple}")

plt.show()



#
#
# ---- Plot Venn diagram for psi and mmseqs databases----
#
#

# ---- Compute intersections ----
inter_all = set_a & set_b

# ---- Save triple intersection ----
with open(consensus_file_double, "w") as out:
    out.write("\n".join(sorted(inter_all)))

print(f"{len(inter_all)} sequences are shared by all three files.")
print(f"Saved as: {consensus_file_double}")

plt.figure(figsize=(6, 6))

v = venn2(
    subsets=(set_a, set_b),
    set_labels=("Novel AMPs (mmseqs)", "Novel AMPs (psi-cd-hit)")
)

# ---- Optional aesthetic improvements ----
for text in v.set_labels:
    text.set_fontsize(14)
for text in v.subset_labels:
    if text:
        text.set_fontsize(12)

# High-resolution output for publication
plt.savefig(venn_diagram_file_double, dpi=400, bbox_inches='tight')
print(f"Saved as: {venn_diagram_file_double}")

plt.show()"""


import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from upsetplot import UpSet, from_contents
import sys

# Configuración de backend
try:
    matplotlib.use("TkAgg")
except:
    matplotlib.use("Agg")

# ---- Rutas (Ajusta según tu entorno) ----
# Asegúrate de que estas rutas sean correctas en tu máquina
file_a = "../data/novel_seqs/mmseqs_novel_seqs_ids.txt"
file_b = "../data/novel_seqs/psi-cd-hit_novel_seqs_ids.txt"
file_c = "../data/novel_seqs/interpro_novel_seqs_ids.txt"
file_d = "../data/novel_seqs/predicted_seqs_ids.txt"

results_dir = Path("../results/venn_diagram/")
results_dir.mkdir(parents=True, exist_ok=True)
output_file = results_dir / "upset_relative_predicted_clean.png"

# ---- Funciones de carga con chequeo ----
def load_ids(filename, label):
    path = Path(filename)
    if not path.exists():
        print(f"ERROR: No se encuentra el archivo para {label}: {path.resolve()}")
        # Creamos sets dummy para que el script corra si faltan archivos
        return set()
    with open(path) as f:
        ids = set(line.strip() for line in f if line.strip())
    print(f"Cargado {label}: {len(ids)} secuencias")
    return ids

print("--- Cargando datos ---")
set_mmseqs = load_ids(file_a, "MMseqs")
set_psicd = load_ids(file_b, "Psi-CD-Hit")
set_interpro = load_ids(file_c, "Interpro")
set_predicted_total = load_ids(file_d, "PREDICTED (Superset)")

# --- PASO CLAVE 1: Definir el denominador ---
TOTAL_N = len(set_predicted_total)
if TOTAL_N == 0:
    print("\nError Crítico: El archivo 'Predicted' está vacío o no se cargó.")
    sys.exit(1)

print(f"\n>>> TOTAL BASE PARA PORCENTAJES (N): {TOTAL_N} <<<")

# --- PASO CLAVE 2: Preparar datos SIN incluir el superset 'Predicted' visualmente ---
# Solo graficamos los subconjuntos, pero usaremos TOTAL_N para los cálculos.
contents = {
    "MMseqs": set_mmseqs,
    "Psi-CD-Hit": set_psicd,
    "Interpro": set_interpro,
    # NOTA: NO incluimos "Predicted" aquí para evitar la fila redundante
}

# Verificación de seguridad: ¿Hay cosas en los subconjuntos que NO estén en Predicted?
union_subsets = set_mmseqs | set_psicd | set_interpro
not_in_predicted = union_subsets - set_predicted_total
if len(not_in_predicted) > 0:
    print(f"ADVERTENCIA: Hay {len(not_in_predicted)} secuencias en los subconjuntos que NO están en el archivo 'Predicted'.")
    print("Los porcentajes podrían ser ligeramente engañosos si 'Predicted' no es un superconjunto perfecto.")

data = from_contents(contents)

# ---- Configuración del Gráfico ----
fig = plt.figure(figsize=(12, 8))

# Instanciamos UpSet SIN mostrar etiquetas automáticas
upset = UpSet(
    data,
    subset_size='count',
    show_counts=False,      
    show_percentages=False, 
    sort_by='cardinality',
    min_subset_size=10 # Filtramos intersecciones muy pequeñas para limpiar el gráfico
)

plot_result = upset.plot(fig=fig)

# ---- PASO CLAVE 3: Etiquetado manual usando TOTAL_N ----
ax_bars = plot_result['intersections']

# Iteramos sobre cada barra (rectángulo) del gráfico
for rect in ax_bars.patches:
    height = rect.get_height() # La cantidad absoluta en esa intersección
    
    # Calculamos el porcentaje respecto al TOTAL_N (Predicted)
    percentage = (height / TOTAL_N) * 100
        
    # Texto a mostrar: Cantidad absoluta + (Porcentaje relativo al total predicted)
    label_text = f"{int(height)}\n({percentage:.1f}%)"
    
    ax_bars.text(
        rect.get_x() + rect.get_width() / 2,
        height + (height * 0.02), # Un poco por encima de la barra
        label_text,
        ha='center', va='bottom', fontsize=10, color='black', fontweight='bold'
    )

# Ajustes finales de estilo
ax_bars.set_ylabel("Tamaño de la Intersección")
# Aumentar el límite Y para que quepan los textos
ax_bars.set_ylim(0, ax_bars.get_ylim()[1] * 1.2) 

plt.suptitle(f"Subgrupos de secuencias y su % respecto al Total Predicted (N={TOTAL_N})", 
             fontsize=15, y=0.99, fontweight='bold')

# Añadir una nota explicativa al pie
fig.text(0.5, 0.01, 
         f"*Nota: Los porcentajes indican qué fracción del total de 'Predicted' ({TOTAL_N}) representa cada barra.*", 
         ha='center', fontsize=9, style='italic')

plt.savefig(output_file, dpi=400, bbox_inches='tight')
print(f"\nGráfico mejorado guardado en: {output_file}")

plt.show()


















import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

# Configuración de backend
try:
    matplotlib.use("TkAgg")
except:
    matplotlib.use("Agg")

# ---- 1. Carga de Datos ----
file_a = "../data/novel_seqs/mmseqs_novel_seqs_ids.txt"
file_b = "../data/novel_seqs/psi-cd-hit_novel_seqs_ids.txt"
file_c = "../data/novel_seqs/interpro_novel_seqs_ids.txt"
file_d = "../data/novel_seqs/predicted_seqs_ids.txt"

results_dir = Path("../results/venn_diagram/")
results_dir.mkdir(parents=True, exist_ok=True)
output_file = results_dir / "custom_upset_with_predicted_row.png"

def load_ids(filename, label):
    path = Path(filename)
    if not path.exists(): 
        print(f"Aviso: {label} no existe.")
        return set()
    with open(path) as f: return set(line.strip() for line in f if line.strip())

mmseqs = load_ids(file_a, "MMseqs")
psi = load_ids(file_b, "Psi")
interpro = load_ids(file_c, "Interpro")
predicted = load_ids(file_d, "Predicted")

TOTAL_PREDICTED = len(predicted)
if TOTAL_PREDICTED == 0: TOTAL_PREDICTED = 1 

# ---- 2. Configuración de Filas y Columnas ----

# Nombres de las filas de la matriz (De ABAJO hacia ARRIBA)
# Ponemos Predicted abajo del todo como la "base"
row_names = ["Predicted", "Interpro", "Psi-CD-Hit", "MMseqs"]

# Definimos las columnas. 
# IMPORTANTE: Como todos son subconjuntos de Predicted, 
# añadimos "Predicted" a la lista 'active' de TODOS.

data_config = [
    # --- BARRA DE REFERENCIA (100%) ---
    {
        "label": "Predicted\n(Universe)", 
        "active": ["Predicted"], 
        "data": predicted,
        "color": "grey" # La pintamos gris para diferenciarla
    },
    
    # --- COLUMNAS INDIVIDUALES (TOTALES) ---
    {
        "label": "MMseqs\n(Total)",   
        "active": ["MMseqs", "Predicted"],    
        "data": mmseqs,
        "color": "black"
    },
    {
        "label": "Psi-CD\n(Total)",   
        "active": ["Psi-CD-Hit", "Predicted"],       
        "data": psi,
        "color": "black"
    },
    {
        "label": "Interpro\n(Total)", 
        "active": ["Interpro", "Predicted"],  
        "data": interpro,
        "color": "black"
    },
    
    # --- INTERSECCIONES DOBLES ---
    {
        "label": "MM & Psi",        
        "active": ["MMseqs", "Psi-CD-Hit", "Predicted"],      
        "data": mmseqs & psi,
        "color": "black"
    },
    {
        "label": "MM & Inter",      
        "active": ["MMseqs", "Interpro", "Predicted"], 
        "data": mmseqs & interpro,
        "color": "black"
    },
    {
        "label": "Psi & Inter",     
        "active": ["Psi-CD-Hit", "Interpro", "Predicted"],    
        "data": psi & interpro,
        "color": "black"
    },
    
    # --- INTERSECCIÓN TRIPLE ---
    {
        "label": "All Three",       
        "active": ["MMseqs", "Psi-CD-Hit", "Interpro", "Predicted"], 
        "data": mmseqs & psi & interpro,
        "color": "black"
    },
]

# ---- 3. Preparación de Valores ----
counts = []
percentages = []
matrix_states = [] 
bar_colors = []

for item in data_config:
    # Intersección segura con el universo
    overlap = item["data"] & predicted
    count = len(overlap)
    pct = (count / TOTAL_PREDICTED) * 100
    
    counts.append(count)
    percentages.append(pct)
    bar_colors.append(item.get("color", "black"))
    
    # Determinar estado de la matriz basado en row_names
    actives = item["active"]
    # Creamos una lista de True/False correspondiente a row_names
    col_state = [name in actives for name in row_names]
    matrix_states.append(col_state)

# ---- 4. Construcción del Gráfico ----
fig = plt.figure(figsize=(14, 9))
# Ajustamos ratios: el gráfico de barras un poco más alto
gs = gridspec.GridSpec(2, 1, height_ratios=[3.5, 1.5]) 

# --- A. Gráfico de Barras (Panel Superior) ---
ax_bar = plt.subplot(gs[0])

x_pos = range(len(data_config))
bars = ax_bar.bar(x_pos, percentages, color=bar_colors, width=0.5)

# Etiquetas encima de las barras
for bar, count, pct in zip(bars, counts, percentages):
    height = bar.get_height()
    ax_bar.text(
        bar.get_x() + bar.get_width()/2, 
        height + (max(percentages)*0.02), 
        f"{count}\n({pct:.1f}%)", 
        ha='center', va='bottom', fontsize=10, fontweight='bold'
    )

ax_bar.set_ylabel(f"% Coverage of Predicted (N={TOTAL_PREDICTED})", fontsize=12)
ax_bar.set_title("Cobertura de Herramientas sobre el conjunto Predicted", fontsize=15, pad=20)
ax_bar.set_ylim(0, max(percentages) * 1.2) 
ax_bar.set_xticks([]) 
ax_bar.set_xlim(-0.6, len(data_config) - 0.4)
ax_bar.spines['top'].set_visible(False)
ax_bar.spines['right'].set_visible(False)
ax_bar.spines['bottom'].set_visible(False)

# --- B. Matriz de Puntos (Panel Inferior) ---
ax_matrix = plt.subplot(gs[1], sharex=ax_bar)

# Dibujar filas (líneas horizontales de guía)
y_positions = range(len(row_names))
for i in y_positions:
    if i % 2 == 0: # Bandas de color alternas para legibilidad
        ax_matrix.axhspan(i - 0.5, i + 0.5, color='#f4f4f4', zorder=0)

# Iterar columnas
for x, state in enumerate(matrix_states):
    # Encontrar índices activos para dibujar la línea vertical
    active_indices = [i for i, is_active in enumerate(state) if is_active]
    
    if len(active_indices) > 1:
        # Línea vertical conectando el punto más bajo con el más alto
        ax_matrix.vlines(x, min(active_indices), max(active_indices), color='black', linewidth=2, zorder=1)
    
    # Dibujar los puntos
    for y, is_active in enumerate(state):
        if is_active:
            color = 'black'
            size = 140
        else:
            color = '#e0e0e0' # Gris muy claro
            size = 50
            
        ax_matrix.scatter(x, y, c=color, s=size, zorder=2)

# Configurar etiquetas de filas
ax_matrix.set_yticks(y_positions)
ax_matrix.set_yticklabels(row_names, fontsize=12, fontweight='medium')
ax_matrix.set_ylim(-0.5, len(row_names) - 0.5)
ax_matrix.set_xlim(-0.6, len(data_config) - 0.4)

# Limpieza visual final
ax_matrix.spines['top'].set_visible(False)
ax_matrix.spines['bottom'].set_visible(False)
ax_matrix.spines['right'].set_visible(False)
ax_matrix.spines['left'].set_visible(False)
ax_matrix.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax_matrix.tick_params(axis='y', length=0) 

plt.tight_layout()
plt.subplots_adjust(hspace=0.0) # Unir visualmente los paneles

plt.savefig(output_file, dpi=400, bbox_inches='tight')
print(f"Gráfico generado en: {output_file}")

plt.show()