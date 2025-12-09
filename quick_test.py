from Bio import SeqIO

dna = 0
prot = 0

for record in SeqIO.parse("test_mixed.fasta", "fasta"):
    seq = str(record.seq)
    # crude classifier just for sanity check, not final rules
    if set(seq) <= set("ACGTUNacgtun"):
        dna += 1
    else:
        prot += 1
    print(record.id, len(seq))

print("DNA:", dna, "Protein:", prot)
