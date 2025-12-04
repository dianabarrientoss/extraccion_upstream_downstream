"""
Extraer regiones upstream/downstream de genes en secuencias genómicas.

Fase A: lógica central para analizar un GFF (genes), cargar un genoma FASTA,
calcular coordenadas upstream y downstream considerando la hebra, extraer secuencias y escribir una salida FASTA única.

Este archivo proporciona:
- load_fasta(path) -> dict {seqid: sequence}
- parse_gff_genes(path) -> lista de diccionarios de genes
- compute_up_down(start, end, strand, up_len, down_len, chrom_len)
- extract_subseq(chrom_seqs, seqid, start, end)
- format_region_header(...)
- write_regions_fasta(output_path, records)

Author: Diana Barrientos
"""

from typing import Dict, List, Tuple


def load_fasta(path: str) -> Dict[str, str]:
    """
    TODO: Cargar un archivo FASTA en un diccionario {seqid: sequence}
    """
    pass


def parse_gff_genes(path: str) -> List[Dict[str, str]]:
    """
    TODO: Analizar un archivo GFF y extraer genes en una lista de diccionarios con:
    id, seqid, start, end, strand
    """
    pass


def compute_up_down(
    start: int,
    end: int,
    strand: str,
    up_len: int,
    down_len: int,
    chrom_len: int
) -> Dict[str, Tuple[int, int]]:
    """
    TODO: Calcular coordenadas upstream y downstream considerando la hebra.
    """
    pass


def extract_subseq(
    chrom_seqs: Dict[str, str],
    seqid: str,
    start: int,
    end: int
) -> str:
    """
    TODO: Extraer una subsecuencia de chrom_seqs dada seqid, start y end.
    """
    pass


def format_region_header(
    gene_id: str,
    region_name: str,
    coords: Tuple[int, int],
    strand: str,
    orig_start: int,
    orig_end: int
) -> str:
    """
    TODO: Formatear el encabezado FASTA de cada región.
    """
    pass


def write_regions_fasta(output_path: str, records: List[Tuple[str, str]]) -> None:
    """
    TODO: Escribir un archivo FASTA con las secuencias extraídas.
    """
    pass


# =========================
# PARTE B - CLI
# =========================

def parse_args():
    """
    TODO: Parsear argumentos desde línea de comandos.
    """
    pass


def resolve_lengths(args) -> Tuple[int, int]:
    """
    TODO: Resolver longitudes upstream y downstream a partir de los argumentos.
    """
    pass


def main():
    """
    TODO: Función principal que coordina todo el flujo del programa.
    type
    """
    pass


if __name__ == "__main__":
    main()