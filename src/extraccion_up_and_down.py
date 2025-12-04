## Extraccion de secuencias upstream y downstream
import os

def load_fasta(file_path):
    # CARGAR EL FASTA EN FORMA DE DICCIONARIO para tenerlo de tipo
    # gene -> secuencia
    sequences = {}
    with open(file_path, 'r') as f:#abrr en lectura
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):# nueva secuencia
                if header:
                    sequences[header] = ''.join(seq_lines)# guardar la secuencia anterior
                header = line[1:]  # Remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            sequences[header] = ''.join(seq_lines)
    return sequences
def parse_gff_genes(path:str)->dict:
    #ANALIZAR EL GFF PARA EXTRAER GENES EN LA LISTA DE DICCIONARIOS
    genes = {}
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')# dividir por tabulador
            if len(parts) < 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            if feature_type == 'gene':
                attr_dict = {}
                for attr in attributes.split(';'):
                    key_value = attr.split('=')
                    if len(key_value) == 2:
                        key, value = key_value
                        attr_dict[key] = value
                gene_id = attr_dict.get('ID', None)
                if gene_id:
                    genes[gene_id] = {
                        'seqid': seqid,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand
                    }
    return genes