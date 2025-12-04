def reverse_complement(seq):
    """ Returns the reverse complement of a DNA sequence. """
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def import_dna_file(file_path):
    """ Imports DNA sequence and returns dict in which the key is the chromosome ID
    and the value is the string that correspodns to the sequence."""
    
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