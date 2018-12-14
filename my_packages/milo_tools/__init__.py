
def reverse_transcribe(seq):
    """reverse transcribes a dna sequence (does not convert any non-atcg/ATCG characters)"""
    watson_crick = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join([watson_crick.setdefault(c, c) for c in seq[::-1]])

def FourLineFastq(handle):
    """
    Reads 4 lines from the file and returns 1, 2, and 4 with newlines stripped
    The only check for fastq format is that line 3 starts with '+'
    """
    while True:
        line = handle.readline()     
        if not line:
            # end of file
            break       
        title = line.rstrip()
        seq = handle.readline().rstrip()
        jnk_line = handle.readline()
        if jnk_line[0] != "+":
            print(title, seq, jnk_line)
            raise ValueError("Looks like this isn’t a strictly 4-line fastq file")
        qual = handle.readline().rstrip()
        yield (title, seq, qual)