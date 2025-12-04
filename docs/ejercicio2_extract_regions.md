# Ejercicio 2 — Extraer regiones upstream/downstream (versión mejorada)

## **Parte A — Implementación base**

Escribe un programa `extract_regions.py` que:

### **Uso**
```
python extract_regions.py --gff genes.gff --fasta genome.fasta --up 500 --down 50 --out regions.fna
```

### **Objetivo**

Para cada `gene` en el archivo GFF:

1. Obtener región **upstream** de tamaño `--up`.
2. Obtener región **downstream** de tamaño `--down`.
3. Considerar **correctamente el strand**:
   - strand `+`: upstream = antes; downstream = después.
   - strand `-`: upstream = después; downstream = antes.
4. **Incluir coordenadas genómicas reales** en los encabezados FASTA.

---

# 🧬 **Formato esperado del encabezado FASTA**

Cada secuencia debe incluir:

- el ID del gene  
- si es upstream o downstream  
- coordenadas genómicas exactas (`start-end`)  
- strand del gen  
- coordenadas originales del gene

### Ejemplo para strand positivo:

```
>araC_upstream 12345-12844 strand=+ original_gene_coords=12845-13500
ATGCTAGCTAGCTAGGCTAGCTACGTAC
>araC_downstream 13501-13550 strand=+ original_gene_coords=12845-13500
TTTTTTGCGCGATTAACCCTT
```

### Ejemplo para strand negativo:

```
>geneX_upstream 8801-9300 strand=- original_gene_coords=8000-8799
ACGTGCGTACGTACGTAGCG
>geneX_downstream 7500-7999 strand=- original_gene_coords=8000-8799
TTAAGCCGTTTTTTGCGGTA
```

---

# **Parte B — Extensión**

Agregar argumento opcional:

```
--type T
```

Ejemplo:

```
python extract_regions.py --gff genes.gff --fasta genome.fasta --up 200 --type CDS
```

Solo procesa features cuyo tipo (columna 3) sea **T**.

---

# Requisitos técnicos

- uso obligatorio de `argparse`
- manejo de archivos FASTA y GFF
- slicing correcto según strand
- validación de límites (no permitir índices < 0)
- uso de funciones, docstrings, comprensión de listas
- Pruebas con asserts + documento de pruebas, o bien con pytest.
