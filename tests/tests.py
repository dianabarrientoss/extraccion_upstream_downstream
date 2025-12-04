"""
Tests para el módulo extraccion_up_and_down.py
Cubre todas las funciones principales con assertions.
"""

import unittest
import tempfile
import os
import sys

# Importar las funciones a probar
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from extraccion_up_and_down import (
    load_fasta,
    parse_gff_genes,
    compute_up_down,
    extract_subseq,
    format_region_header,
    write_regions_fasta,
    resolve_lengths
)


class TestLoadFasta(unittest.TestCase):
    """Tests para la función load_fasta"""
    
    def test_load_fasta_single_sequence(self):
        """Test carga de un FASTA con una sola secuencia"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\n")
            f.write("ATGCATGC\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = load_fasta(temp_path)
            assert "chr1" in result, "seqid 'chr1' no se encuentra en el diccionario"
            assert result["chr1"] == "ATGCATGC", f"Secuencia esperada 'ATGCATGC', obtuvo '{result['chr1']}'"
        finally:
            os.unlink(temp_path)
    
    def test_load_fasta_multiline_sequence(self):
        """Test carga de secuencias multilínea"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGC\n")
            f.write("ATGC\n")
            f.write("TTAA\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = load_fasta(temp_path)
            assert result["seq1"] == "ATGCATGCTTAA", f"Secuencias no se unieron correctamente: {result['seq1']}"
        finally:
            os.unlink(temp_path)
    
    def test_load_fasta_multiple_sequences(self):
        """Test carga de múltiples secuencias"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\nATGC\n")
            f.write(">chr2\nTTAA\n")
            f.write(">chr3\nGGCC\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = load_fasta(temp_path)
            assert len(result) == 3, f"Se esperaban 3 secuencias, obtuvo {len(result)}"
            assert result["chr1"] == "ATGC", f"chr1 incorrecto: {result['chr1']}"
            assert result["chr2"] == "TTAA", f"chr2 incorrecto: {result['chr2']}"
            assert result["chr3"] == "GGCC", f"chr3 incorrecto: {result['chr3']}"
        finally:
            os.unlink(temp_path)
    
    def test_load_fasta_lowercase_conversion(self):
        """Test conversión a mayúsculas"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\n")
            f.write("atgc\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = load_fasta(temp_path)
            assert result["chr1"] == "ATGC", f"No se convirtió a mayúsculas: {result['chr1']}"
        finally:
            os.unlink(temp_path)
    
    def test_load_fasta_header_with_spaces(self):
        """Test encabezado con espacios (solo toma primer parte)"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1 description here\n")
            f.write("ATGC\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = load_fasta(temp_path)
            assert "chr1" in result, "seqid 'chr1' no encontrado"
            assert "chr1 description here" not in result, "No debería tener descripción en el id"
        finally:
            os.unlink(temp_path)


class TestParseGffGenes(unittest.TestCase):
    """Tests para la función parse_gff_genes"""
    
    def test_parse_gff_single_gene(self):
        """Test parseo de un gen"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = parse_gff_genes(temp_path, feature_type_filter="gene")
            assert len(result) == 1, f"Se esperaba 1 gen, obtuvo {len(result)}"
            assert result[0]["id"] == "gene1", f"ID incorrecto: {result[0]['id']}"
            assert result[0]["seqid"] == "chr1", f"seqid incorrecto: {result[0]['seqid']}"
            assert result[0]["start"] == 1000, f"start incorrecto: {result[0]['start']}"
            assert result[0]["end"] == 2000, f"end incorrecto: {result[0]['end']}"
            assert result[0]["strand"] == "+", f"strand incorrecto: {result[0]['strand']}"
        finally:
            os.unlink(temp_path)
    
    def test_parse_gff_multiple_genes(self):
        """Test parseo de múltiples genes"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n")
            f.write("chr1\t.\tgene\t3000\t4000\t.\t-\t.\tID=gene2\n")
            f.write("chr2\t.\tgene\t500\t1500\t.\t+\t.\tID=gene3\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = parse_gff_genes(temp_path, feature_type_filter="gene")
            assert len(result) == 3, f"Se esperaban 3 genes, obtuvo {len(result)}"
            assert result[0]["id"] == "gene1", f"Primer gen ID incorrecto"
            assert result[1]["strand"] == "-", f"Segundo gen strand incorrecto"
            assert result[2]["seqid"] == "chr2", f"Tercer gen seqid incorrecto"
        finally:
            os.unlink(temp_path)
    
    def test_parse_gff_filter_by_type(self):
        """Test filtrado por tipo de feature"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n")
            f.write("chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=exon1\n")
            f.write("chr1\t.\tCDS\t1100\t1400\t.\t+\t.\tID=cds1\n")
            f.flush()
            temp_path = f.name
        
        try:
            result_genes = parse_gff_genes(temp_path, feature_type_filter="gene")
            result_exons = parse_gff_genes(temp_path, feature_type_filter="exon")
            result_cds = parse_gff_genes(temp_path, feature_type_filter="CDS")
            
            assert len(result_genes) == 1, f"Se esperaba 1 gene, obtuvo {len(result_genes)}"
            assert len(result_exons) == 1, f"Se esperaba 1 exon, obtuvo {len(result_exons)}"
            assert len(result_cds) == 1, f"Se esperaba 1 CDS, obtuvo {len(result_cds)}"
            assert result_genes[0]["id"] == "gene1", "ID de gene incorrecto"
            assert result_exons[0]["id"] == "exon1", "ID de exon incorrecto"
            assert result_cds[0]["id"] == "cds1", "ID de CDS incorrecto"
        finally:
            os.unlink(temp_path)
    
    def test_parse_gff_with_comments(self):
        """Test ignora comentarios"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("# This is a comment\n")
            f.write("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n")
            f.write("# Another comment\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = parse_gff_genes(temp_path, feature_type_filter="gene")
            assert len(result) == 1, f"Los comentarios no se ignoraron correctamente: {len(result)}"
        finally:
            os.unlink(temp_path)
    
    def test_parse_gff_id_extraction(self):
        """Test extracción correcta del ID con atributos múltiples"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=AT1G01010;Name=NAC001\n")
            f.flush()
            temp_path = f.name
        
        try:
            result = parse_gff_genes(temp_path, feature_type_filter="gene")
            assert result[0]["id"] == "AT1G01010", f"ID no extraído correctamente: {result[0]['id']}"
        finally:
            os.unlink(temp_path)


class TestComputeUpDown(unittest.TestCase):
    """Tests para la función compute_up_down"""
    
    def test_compute_upstream_positive_strand(self):
        """Test cálculo upstream para hebra positiva"""
        result = compute_up_down(start=1000, end=2000, strand="+", 
                                 up_len=500, down_len=500, chrom_len=10000)
        assert result["upstream"] == (500, 999), f"Upstream incorrecto: {result['upstream']}"
    
    def test_compute_downstream_positive_strand(self):
        """Test cálculo downstream para hebra positiva"""
        result = compute_up_down(start=1000, end=2000, strand="+", 
                                 up_len=500, down_len=500, chrom_len=10000)
        assert result["downstream"] == (2001, 2500), f"Downstream incorrecto: {result['downstream']}"
    
    def test_compute_upstream_negative_strand(self):
        """Test cálculo upstream para hebra negativa"""
        result = compute_up_down(start=1000, end=2000, strand="-", 
                                 up_len=500, down_len=500, chrom_len=10000)
        assert result["upstream"] == (2001, 2500), f"Upstream negativo incorrecto: {result['upstream']}"
    
    def test_compute_downstream_negative_strand(self):
        """Test cálculo downstream para hebra negativa"""
        result = compute_up_down(start=1000, end=2000, strand="-", 
                                 up_len=500, down_len=500, chrom_len=10000)
        assert result["downstream"] == (500, 999), f"Downstream negativo incorrecto: {result['downstream']}"
    
    def test_compute_upstream_at_chromosome_start(self):
        """Test upstream no excede inicio del cromosoma"""
        result = compute_up_down(start=100, end=200, strand="+", 
                                 up_len=500, down_len=100, chrom_len=10000)
        assert result["upstream"][0] == 1, f"Upstream inicio no respeta límite inferior: {result['upstream'][0]}"
    
    def test_compute_downstream_at_chromosome_end(self):
        """Test downstream no excede fin del cromosoma"""
        result = compute_up_down(start=9800, end=9900, strand="+", 
                                 up_len=100, down_len=500, chrom_len=10000)
        assert result["downstream"][1] == 10000, f"Downstream fin no respeta límite superior: {result['downstream'][1]}"
    
    def test_compute_invalid_strand(self):
        """Test lanza error con hebra inválida"""
        try:
            compute_up_down(start=1000, end=2000, strand="X", 
                           up_len=500, down_len=500, chrom_len=10000)
            assert False, "Debería lanzar ValueError para hebra inválida"
        except ValueError as e:
            assert "Invalid strand" in str(e), f"Mensaje de error incorrecto: {e}"
    
    def test_compute_zero_lengths(self):
        """Test con longitudes cero"""
        result = compute_up_down(start=1000, end=2000, strand="+", 
                                 up_len=0, down_len=0, chrom_len=10000)
        assert result["upstream"] == (1000, 999), f"Upstream con len 0 incorrecto: {result['upstream']}"
        assert result["downstream"] == (2001, 2000), f"Downstream con len 0 incorrecto: {result['downstream']}"


class TestExtractSubseq(unittest.TestCase):
    """Tests para la función extract_subseq"""
    
    def test_extract_valid_region(self):
        """Test extracción de región válida"""
        seqs = {"chr1": "ATGCATGCTTAA"}
        result = extract_subseq(seqs, "chr1", 1, 4)
        assert result == "ATGC", f"Extracción inválida: '{result}' != 'ATGC'"
    
    def test_extract_full_sequence(self):
        """Test extracción de secuencia completa"""
        seqs = {"chr1": "ATGC"}
        result = extract_subseq(seqs, "chr1", 1, 4)
        assert result == "ATGC", f"Extracción completa inválida: '{result}'"
    
    def test_extract_single_base(self):
        """Test extracción de una sola base"""
        seqs = {"chr1": "ATGCATGC"}
        result = extract_subseq(seqs, "chr1", 3, 3)
        assert result == "G", f"Extracción de una base inválida: '{result}' != 'G'"
    
    def test_extract_beyond_sequence_end(self):
        """Test extracción que excede el fin de la secuencia"""
        seqs = {"chr1": "ATGC"}
        result = extract_subseq(seqs, "chr1", 2, 100)
        assert result == "TGC", f"Extracción truncada inválida: '{result}' != 'TGC'"
    
    def test_extract_before_sequence_start(self):
        """Test extracción que comienza antes del inicio"""
        seqs = {"chr1": "ATGC"}
        result = extract_subseq(seqs, "chr1", 0, 2)
        assert result == "AT", f"Extracción desde 0 inválida: '{result}' != 'AT'"
    
    def test_extract_empty_region(self):
        """Test extracción de región vacía"""
        seqs = {"chr1": "ATGC"}
        result = extract_subseq(seqs, "chr1", 0, -1)
        assert result == "", f"Región vacía debería retornar cadena vacía, obtuvo: '{result}'"
    
    def test_extract_invalid_seqid(self):
        """Test lanza error con seqid no encontrado"""
        seqs = {"chr1": "ATGC"}
        try:
            extract_subseq(seqs, "chr2", 1, 4)
            assert False, "Debería lanzar KeyError para seqid no encontrado"
        except KeyError:
            pass  # Es lo esperado


class TestFormatRegionHeader(unittest.TestCase):
    """Tests para la función format_region_header"""
    
    def test_format_header_basic(self):
        """Test formato básico de encabezado"""
        header = format_region_header("gene1", "upstream", (100, 200), "+", 500, 1000)
        assert "gene1_upstream" in header, f"ID y región no en encabezado: {header}"
        assert "100-200" in header, f"Coordenadas no en encabezado: {header}"
        assert "strand=+" in header, f"Strand no en encabezado: {header}"
        assert "original_gene_coords=500-1000" in header, f"Coords originales no en encabezado: {header}"
    
    def test_format_header_negative_strand(self):
        """Test formato con hebra negativa"""
        header = format_region_header("gene2", "downstream", (50, 150), "-", 200, 300)
        assert "gene2_downstream" in header, "ID y región incorrecto"
        assert "strand=-" in header, "Strand negativo incorrecto"
    
    def test_format_header_empty_region(self):
        """Test formato para región vacía"""
        header = format_region_header("gene3", "upstream", (0, -1), "+", 100, 200)
        assert "empty" in header, "Región vacía no indicada con 'empty'"


class TestWriteRegionsFasta(unittest.TestCase):
    """Tests para la función write_regions_fasta"""
    
    def test_write_single_record(self):
        """Test escritura de un registro"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_path = f.name
        
        try:
            records = [
                (">gene1_upstream 100-200 strand=+ original_gene_coords=500-1000", "ATGCATGC")
            ]
            write_regions_fasta(temp_path, records)
            
            with open(temp_path, 'r') as f:
                content = f.read()
            
            assert ">gene1_upstream" in content, "Encabezado no escrito"
            assert "ATGCATGC" in content, "Secuencia no escrita"
        finally:
            os.unlink(temp_path)
    
    def test_write_multiple_records(self):
        """Test escritura de múltiples registros"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_path = f.name
        
        try:
            records = [
                (">gene1_upstream 100-200 strand=+", "ATGC"),
                (">gene1_downstream 300-400 strand=+", "TTAA")
            ]
            write_regions_fasta(temp_path, records)
            
            with open(temp_path, 'r') as f:
                content = f.read()
            
            assert ">gene1_upstream" in content, "Primer encabezado no escrito"
            assert ">gene1_downstream" in content, "Segundo encabezado no escrito"
        finally:
            os.unlink(temp_path)
    
    def test_write_long_sequence_wrapped(self):
        """Test que secuencias largas se dividen en líneas de 60 caracteres"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            temp_path = f.name
        
        try:
            long_seq = "A" * 150
            records = [(">gene1_upstream", long_seq)]
            write_regions_fasta(temp_path, records)
            
            with open(temp_path, 'r') as f:
                lines = f.readlines()
            
            assert lines[0].strip() == ">gene1_upstream", "Encabezado incorrecto"
            assert len(lines[1].strip()) == 60, f"Primera línea debe tener 60 chars, tiene {len(lines[1].strip())}"
            assert len(lines[2].strip()) == 60, f"Segunda línea debe tener 60 chars, tiene {len(lines[2].strip())}"
            assert len(lines[3].strip()) == 30, f"Tercera línea debe tener 30 chars, tiene {len(lines[3].strip())}"
        finally:
            os.unlink(temp_path)


class TestResolveLengths(unittest.TestCase):
    """Tests para la función resolve_lengths"""
    
    def test_resolve_with_bp(self):
        """Test resolución con --bp"""
        class Args:
            bp = 2000
            upstream = 1000
            downstream = 1000
        
        up_len, down_len = resolve_lengths(Args())
        assert up_len == 2000, f"Upstream debe ser 2000, obtuvo {up_len}"
        assert down_len == 2000, f"Downstream debe ser 2000, obtuvo {down_len}"
    
    def test_resolve_without_bp(self):
        """Test resolución sin --bp usa valores individuales"""
        class Args:
            bp = None
            upstream = 1500
            downstream = 2500
        
        up_len, down_len = resolve_lengths(Args())
        assert up_len == 1500, f"Upstream debe ser 1500, obtuvo {up_len}"
        assert down_len == 2500, f"Downstream debe ser 2500, obtuvo {down_len}"
    
    def test_resolve_bp_overrides_individual(self):
        """Test --bp prevalece sobre valores individuales"""
        class Args:
            bp = 3000
            upstream = 1000
            downstream = 1000
        
        up_len, down_len = resolve_lengths(Args())
        assert up_len == 3000, f"Upstream debe ser 3000, obtuvo {up_len}"
        assert down_len == 3000, f"Downstream debe ser 3000, obtuvo {down_len}"
        assert up_len != 1000, "bp no debería ser ignorado"


class TestIntegration(unittest.TestCase):
    """Tests de integración completos"""
    
    def test_full_workflow_positive_strand(self):
        """Test flujo completo para hebra positiva"""
        # Crear archivos temporales
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\n")
            f.write("A" * 500 + "ATGC" * 250 + "T" * 500 + "\n")
            fasta_path = f.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("chr1\t.\tgene\t501\t1500\t.\t+\t.\tID=gene1\n")
            gff_path = f.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            output_path = f.name
        
        try:
            # Cargar datos
            seqs = load_fasta(fasta_path)
            genes = parse_gff_genes(gff_path, feature_type_filter="gene")
            
            assert len(genes) == 1, f"Se esperaba 1 gen, obtuvo {len(genes)}"
            assert genes[0]["strand"] == "+", f"Strand debe ser +, obtuvo {genes[0]['strand']}"
            
            # Procesar gen
            gene = genes[0]
            chrom_len = len(seqs["chr1"])
            coords = compute_up_down(gene["start"], gene["end"], gene["strand"], 
                                     500, 500, chrom_len)
            
            up_seq = extract_subseq(seqs, "chr1", coords["upstream"][0], coords["upstream"][1])
            down_seq = extract_subseq(seqs, "chr1", coords["downstream"][0], coords["downstream"][1])
            
            assert len(up_seq) > 0, "Upstream debería tener contenido"
            assert len(down_seq) > 0, "Downstream debería tener contenido"
        finally:
            os.unlink(fasta_path)
            os.unlink(gff_path)
            os.unlink(output_path)
    
    def test_full_workflow_negative_strand(self):
        """Test flujo completo para hebra negativa"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">chr1\n")
            f.write("ATGC" * 300 + "\n")
            fasta_path = f.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.gff', delete=False) as f:
            f.write("chr1\t.\tgene\t800\t900\t.\t-\t.\tID=gene_neg\n")
            gff_path = f.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            output_path = f.name
        
        try:
            seqs = load_fasta(fasta_path)
            genes = parse_gff_genes(gff_path, feature_type_filter="gene")
            
            assert len(genes) == 1, "Debe haber 1 gen"
            assert genes[0]["strand"] == "-", "Strand debe ser negativo"
            
            gene = genes[0]
            chrom_len = len(seqs["chr1"])
            coords = compute_up_down(gene["start"], gene["end"], gene["strand"], 
                                     200, 200, chrom_len)
            
            # Para hebra negativa, upstream está después del gen
            assert coords["upstream"][0] > gene["end"], "Para hebra -, upstream debe estar después"
        finally:
            os.unlink(fasta_path)
            os.unlink(gff_path)
            os.unlink(output_path)


if __name__ == "__main__":
    # Ejecutar todos los tests con verbosidad
    unittest.main(verbosity=2)
