def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def import_dna_file(file_path):
    """Importa sequência de ADN e devolve dicionário cuja chave é o ID do cromosoma e o item é a string correspondente à sequência"""
    chr_dict = {}
    current_chr = None
    seq_buffer = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_chr and seq_buffer:
                    chr_dict[current_chr] = ''.join(seq_buffer)
                current_chr = line.strip().split(" ")[0]
                seq_buffer = []
            else:
                seq_buffer.append(line)

        if current_chr and seq_buffer:
            chr_dict[current_chr] = ''.join(seq_buffer)

    return chr_dict