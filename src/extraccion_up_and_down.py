"""Extraer regiones upstream/downstream de genes en secuencias genómicas.

Fase A: lógica central para analizar un GFF (genes), cargar un genoma FASTA,
calcular coordenadas upstream y downstream considerando la hebra, extraer secuencias y escribir una salida FASTA única.
Este archivo proporciona:
- load_fasta(path) -> dict {seqid: sequence}
- parse_gff_genes(path) -> lista de diccionarios de genes con claves:
    id, seqid, start, end, strand 
- compute_up_down(gene, up_len, down_len, chrom_len) -> dict con coords para up/down
- extract_subseq(chrom_seqs, seqid, start, end) -> cadena de secuencia
- write_regions_fasta(output_path, records) -> escribe fasta con encabezados personalizados
Author: Diana Barrientos con ayuda de Pylia
"""

from typing import Dict, List, Tuple
import re


def load_fasta(path: str) -> Dict[str, str]:#función para cargar archivo fasta
    """
    cargar un archivo FASTA en un diccionario {seqid: sequence}
    """
    seqs: Dict[str, str] = {}# permitirá gaurdar tipo "gene1": ["ATGC", "TTAA"]
    current_id = None #todavía no se gaurdarn como id
    with open(path, "r") as fh:#abre el archivo fasta en modo lectura
        for line in fh:#recorre cada línea del archivo
            line = line.rstrip("\n") #quitar salto de línea evitando errores
            #al unir secuencias multilínea
            if not line:#si no tiene nada la linea continúa, o sea si está vacía
                continue
            if line.startswith(">"):# es un encabezado que se tonme como id
                current_id = line[1:].split()[0]#saca el id sin el ">" 
                #separa por espacios y toma la primera parte
                seqs[current_id] = []#entrada del diccionario con id y lista vacía
            else:
                if current_id is None:#cuando no es gen
                    raise ValueError("FASTA file malformed: sequence data before header")#lanza error
                seqs[current_id].append(line.strip())#agrega linea de secuencia a la lista
    #une las listas de secuencias en cadenas y convierte a mayúsculas
    return {k: "".join(v).upper() for k, v in seqs.items()}#diccionario ya formado


def parse_gff_genes(path: str, feature_type_filter: str = "gene") -> List[Dict[str, str]]:
    """
    Analiza un archivo GFF y extrae features del tipo `feature_type_filter` en una lista de diccionarios con claves:
    id, seqid, start, end, strand
    Por defecto procesa 'gene' para mantener comportamiento previo.
    """
    genes: List[Dict[str, str]] = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line:
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            # ahora comparamos contra el filtro dinámico en lugar de "gene" fijo
            if feature_type != feature_type_filter:
                continue
            match = re.search(r"ID=([^;]+)", attributes)
            gene_id = match.group(1) if match else "unknown"
            gene_info = {
                "id": gene_id,
                "seqid": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand
            }
            genes.append(gene_info)
    return genes

def compute_up_down(#función para calcular coordenadas upstream y downstream
         start: int, end: int, strand: str, up_len: int, down_len: int, chrom_len: int #info de la sequencia upstream y downstream
         #tal cual dice donde empieza el gen y donde termina, la hebra, longitud en bases upstream y downstream y longitud del cromosoma
) -> Dict[str, Tuple[int, int]]:#guarda las coordenadas en un diccionario en tupla incio y fin para que no se modifiquen
    """
    Calcula las coordenadas upstream y downstream considerando la hebra.
    Devuelve un diccionario con claves 'upstream' y 'downstream' y valores como tuplas (start, end).
    """
    if strand == "+": #es hebra positiva?
        up_start = max(1, start - up_len)#busca el inicio upstream o sea a la izquierda del gen, 
        #el intero tiene q ser  mayor a 1, fuerza a que no sea menor a 1
        up_end = start - 1#hace que el final upstream sea justo antes del inicio del gen
            #entonces es hebra positiva, busca al lado izquierdo del gen, justo una base antes del gen
        down_start = end + 1#downstream en hebra + justo una base después del gen
        down_end = min(chrom_len, end + down_len)#el final downstream no puede exceder la longitud del cromosoma
    elif strand == "-":#entonces es hebra negativa, busca al lado derecho del gen
        up_start = end + 1#upstream en hebra - justo una base después del gen
        up_end = min(chrom_len, end + up_len)#sin pasar la longitud del cromosoma, el final upstream es después del gen
        down_start = max(1, start - down_len)#retrocede hacia la izquierda del gen sin pasar de 1
        down_end = start - 1#antes del inicio del gen
    else:
        raise ValueError(f"Invalid strand: {strand}")#lanza error si la hebra no es + o -
    
    return {#diccionario con las coordenadas
        "upstream": (up_start, up_end),#coordenadas de inicio y fin upstream como tupla
        "downstream": (down_start, down_end)#coordenadas de inicio y fin downstream como tupla
    }


def extract_subseq( #función para extraer subsecuencia de una secuencia cromosómica
        chrom_seqs: Dict[str, str],#diccionario con secuencias cromosómicas
        seqid: str,#id de la secuencia
        start: int,#inicio de la subsecuencia
        end: int#fin de la subsecuencia 
) -> str:#devuelve la subsecuencia como cadena
    """
    Extrae una subsecuencia de chrom_seqs dada seqid, start y end (1-based, inclusivo).
    """
    if (start, end) == (0, -1):# caso especial para secuencia vacía
        return ""#devuelve cadena vacía
    seq = chrom_seqs[seqid]  # may raise KeyError cuando seqid no existe
    zero_start = max(0, start - 1)# convertir a 0 el índice del inicio
    zero_end = min(len(seq), end) #calcula limite superior sin pasarse del largo de la secuencia
    return seq[zero_start:zero_end]#extrae y devuelve la subsecuencia

def format_region_header(gene_id: str, region_name: str, coords: Tuple[int, int], strand: str, orig_start: int, orig_end: int) -> str:
    #función para formatear encabezado de región
    start, end = coords#desempaqueta las coordenadas en inicio y fin
    if start == 0 and end == -1:# caso especial para secuencia vacía
        coords_str = "empty"#si es secuencia vacía lo indica
    else:
        coords_str = f"{start}-{end}"#formatea las coordenadas como cadena
    return f">{gene_id}_{region_name} {coords_str} strand={strand} original_gene_coords={orig_start}-{orig_end}"
#regresa el encabezado formateado con id del gen, nombre de la región, coordenadas, hebra y coordenadas originales del gen


def write_regions_fasta(output_path: str, records: List[Tuple[str, str]]) -> None:
    """
    Escribe un archivo FASTA con las secuencias extraídas.
    records es una lista de tuplas (header, sequence).
    """
    with open(output_path, "w") as fh:#abre el archivo de salida en modo escritura
        for header, sequence in records:#recorre cada tupla de encabezado y secuencia
            fh.write(f"{header}\n")#escribe el encabezado en el archivo
            #escribe la secuencia en líneas de 60 caracteres
            for i in range(0, len(sequence), 60):#recorre la secuencia en pasos de 60
                fh.write(f"{sequence[i:i+60]}\n")#escribe la secuencia en bloques de 60 caracteres por línea  

#PARTE B

def parse_args():#función para parsear argumentos de línea de comandos
    import argparse
    parser = argparse.ArgumentParser(
        description="Extraer regiones upstream/downstream de genes en secuencias genómicas."
    )
    parser.add_argument("fasta", help="Archivo FASTA con secuencias genómicas.")
    parser.add_argument("gff", help="Archivo GFF con anotaciones de genes.")
    parser.add_argument("output", help="Archivo FASTA de salida con regiones extraídas.")
    parser.add_argument("--bp", type=int,
                   help="Single length for upstream and downstream.")
    parser.add_argument("--upstream", type=int, default=1000, help="Longitud upstream (por defecto: 1000).")
    parser.add_argument("--downstream", type=int, default=1000, help="Longitud downstream (por defecto: 1000).")
    parser.add_argument("--type", dest="feature_type", default="gene",
                        help="Solo procesar features cuyo tipo (columna 3) sea T (por defecto: 'gene').")
    return parser.parse_args()

def resolve_lengths(args) -> Tuple[int, int]:
    #función para resolver longitudes upstream y downstream
    if args.bp is not None:#si se especifica --bp
        up_len = down_len = args.bp#ambas longitudes son iguales a bp
    else:
        up_len = args.upstream#si no, usa los valores individuales
        down_len = args.downstream
    return up_len, down_len

def main():
    args = parse_args()#parsea los argumentos de línea de comandos
    up_len, down_len = resolve_lengths(args)#resuelve las longitudes upstream y downstream

    chrom_seqs = load_fasta(args.fasta)#carga las secuencias cromosómicas desde el archivo FASTA
    # pasar el filtro de tipo de feature al parser del GFF
    genes = parse_gff_genes(args.gff, feature_type_filter=args.feature_type)#analiza el archivo GFF y extrae los genes

    records: List[Tuple[str, str]] = []#lista para guardar las tuplas de encabezado y secuencia

    for gene in genes:#recorre cada gen
        seqid = gene["seqid"]#obtiene el id de la secuencia

        # Si el seqid no está en el FASTA, avisar y saltar el gen en lugar de fallar.
        if seqid not in chrom_seqs:
            print(f"Warning: seqid '{seqid}' for gene '{gene['id']}' not found in FASTA. Skipping.", flush=True)
            continue

        chrom_len = len(chrom_seqs[seqid])#longitud del cromosoma correspondiente
        coords = compute_up_down( #calcula las coordenadas upstream y downstream
            start=gene["start"],#inicio del gen
            end=gene["end"],#fin del gen
            strand=gene["strand"],#hebra del gen
            up_len=up_len,# longitud upstream
            down_len=down_len,# longitud downstream
            chrom_len=chrom_len# longitud del cromosoma
        )
        #empezamos a extraer las secuencias upstream y downstream
        #extrae y guarda la secuencia upstream
        up_seq = extract_subseq(#extrae la subsecuencia upstream
            chrom_seqs,#diccionario de secuencias cromosómicas
            seqid,#id de la secuencia
            coords["upstream"][0],#inicio upstream
            coords["upstream"][1]#fin upstream
        )
        up_header = format_region_header(#formatea el encabezado upstream
            gene_id=gene["id"],#id del gen
            region_name="upstream",#nombre de la región
            coords=coords["upstream"],#coordenadas upstream
            strand=gene["strand"],#hebra del gen
            orig_start=gene["start"],#  inicio original del gen
            orig_end=gene["end"]# fin original del gen
        )
        records.append((up_header, up_seq))#agrega la tupla a la lista

        #extrae y guarda la secuencia downstream
        down_seq = extract_subseq(#extrae la subsecuencia downstream
            chrom_seqs,#diccionario de secuencias cromosómicas
            seqid,#id de la secuencia
            coords["downstream"][0],#inicio downstream
            coords["downstream"][1]#fin downstream
        )
        down_header = format_region_header(#formatea el encabezado downstream
            gene_id=gene["id"],#id del gen
            region_name="downstream",#nombre de la región
            coords=coords["downstream"],#coordenadas downstream
            strand=gene["strand"],#hebra del gen
            orig_start=gene["start"],# inicio original del gen
            orig_end=gene["end"]# fin original del gen
        )
        records.append((down_header, down_seq))#agrega la tupla a la lista

    write_regions_fasta(args.output, records)#escribe las regiones extraídas en el archivo FASTA de salida
    print(f"Regiones extraídas escritas en {args.output}")#mensaje de confirmación

# Añadir guard para permitir ejecución directa desde la línea de comandos
if __name__ == "__main__":
    main()