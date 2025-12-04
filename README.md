# Extracción de Regiones Upstream/Downstream

Script para extraer regiones upstream y downstream de genes en secuencias genómicas a partir de archivos FASTA y GFF.

## Instalación

No requiere dependencias externas. Solo Python 3.7+

## Uso General

```bash
python src/extraccion_up_and_down.py <fasta> <gff> <output> [opciones]
```

### Argumentos Requeridos

- `<fasta>`: Archivo FASTA con secuencias genómicas
- `<gff>`: Archivo GFF con anotaciones de genes
- `<output>`: Archivo FASTA de salida con regiones extraídas

## Todas las Formas de Ejecutar el Código

### 1. **Forma Básica (Valores por defecto)**
Extrae 1000 pb upstream y 1000 pb downstream para todos los genes:

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta
```

**Parámetros por defecto:**
- upstream: 1000 bp
- downstream: 1000 bp
- type: gene

---

### 2. **Especificar Longitudes Diferentes (Upstream y Downstream Separados)**

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --upstream 2000 --downstream 1500
```

---

### 3. **Usar Longitud Única para Ambos Lados (--bp)**

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --bp 5000
```

**Nota:** Si se especifica `--bp`, los valores de `--upstream` y `--downstream` son ignorados.

---

### 4. **Procesar Solo Upstream**
Establecer downstream a 0:

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --upstream 1000 --downstream 0
```

---

### 5. **Procesar Solo Downstream**
Establecer upstream a 0:

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --upstream 0 --downstream 1000
```

---

### 6. **Filtrar por Tipo de Feature**

Procesar solo features de tipo `exon`:

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --type exon
```

Procesar solo features de tipo `mRNA`:

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --type mRNA
```

Procesar solo features de tipo `CDS`:

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --type CDS
```

**Nota:** Por defecto se procesan features de tipo `gene`.

---

### 7. **Combinar Todas las Opciones**

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --upstream 3000 --downstream 2000 --type exon
```

---

### 8. **Usar --bp con --type**

```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --bp 5000 --type mRNA
```

---

### 9. **Obtener Ayuda**

Ver todos los parámetros disponibles:

```bash
python src/extraccion_up_and_down.py --help
```

Salida:
```
usage: extraccion_up_and_down.py [-h] [--bp BP] [--upstream UPSTREAM] [--downstream DOWNSTREAM] [--type FEATURE_TYPE] fasta gff output

Extraer regiones upstream/downstream de genes en secuencias genómicas.

positional arguments:
  fasta                 Archivo FASTA con secuencias genómicas.
  gff                   Archivo GFF con anotaciones de genes.
  output                Archivo FASTA de salida con regiones extraídas.

optional arguments:
  -h, --help            show this help message and exit
  --bp BP               Single length for upstream and downstream.
  --upstream UPSTREAM   Longitud upstream (por defecto: 1000).
  --downstream DOWNSTREAM
                        Longitud downstream (por defecto: 1000).
  --type FEATURE_TYPE   Solo procesar features cuyo tipo sea FEATURE_TYPE (por defecto: 'gene').
```

---

## Ejemplos Prácticos

### Ejemplo 1: Análisis de promotores (upstream solamente)
```bash
python src/extraccion_up_and_down.py chromosomes.fasta annotations.gff promoters.fasta --upstream 2000 --downstream 0
```

### Ejemplo 2: Análisis de regiones completas con distintos tipos de features
```bash
python src/extraccion_up_and_down.py data/genome.fa data/genes.gff results/regions.fa --bp 1000 --type CDS
```

### Ejemplo 3: Análisis asimétrico (más upstream que downstream)
```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff output.fasta --upstream 5000 --downstream 500
```

### Ejemplo 4: Regiones amplias para análisis de elementos regulatorios
```bash
python src/extraccion_up_and_down.py genome.fasta genes.gff regulatory_regions.fasta --bp 10000
```

---

## Características

- ✅ Maneja hebras positivas (+) y negativas (-)
- ✅ Evita coordenadas fuera de rango
- ✅ Soporta secuencias multilínea en FASTA
- ✅ Filtra features por tipo dinámicamente
- ✅ Genera encabezados informativos con coordenadas
- ✅ Valida que los seqid del GFF existan en el FASTA
- ✅ Salida FASTA con 60 caracteres por línea

---

## Formato de Salida

Cada secuencia extraída genera un encabezado con la siguiente estructura:

```
>gene_id_region coords strand=strand original_gene_coords=start-end
```

**Ejemplo:**
```
>AT1G01010_upstream 87500-89500 strand=+ original_gene_coords=89501-91000
ATGCATGCATGC...
>AT1G01010_downstream 91001-93001 strand=+ original_gene_coords=89501-91000
GCTAGCTAGCTA...
```

---

## Notas Importantes

1. **Coordenadas GFF**: Se asumen coordenadas 1-based (estándar en GFF)
2. **Secuencias complementarias**: Para la hebra negativa (-), las coordenadas se invierten automáticamente respecto al gen
3. **Seqid no encontrado**: Si un seqid del GFF no existe en el FASTA, el script muestra una advertencia y continúa con el siguiente gen
4. **Archivos grandes**: El script carga todo el FASTA en memoria; para genomas muy grandes, considerar procesamiento por chunks

---

## Validación

Antes de ejecutar, asegúrate de que:
- Los archivos FASTA y GFF existen y son accesibles
- El formato de ambos archivos es válido
- El directorio de salida tiene permisos de escritura
- Tienes espacio suficiente en disco para el archivo de salida
