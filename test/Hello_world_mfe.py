import RNA

# The RNA sequence
seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"

# compute minimum free energy (MFE) and corresponding structure
(ss, mfe) = RNA.fold(seq)

# print output
print("{}\n{} [ {:6.2f} ]".format(seq, ss, mfe))


# The RNA sequence alignment
sequences = [
    "CUGCCUCACAACGUUUGUGCCUCAGUUACCCGUAGAUGUAGUGAGGGU",
    "CUGCCUCACAACAUUUGUGCCUCAGUUACUCAUAGAUGUAGUGAGGGU",
    "---CUCGACACCACU---GCCUCGGUUACCCAUCGGUGCAGUGCGGGU",
]

# # compute the consensus sequence
cons = RNA.consensus(sequences)

# predict Minmum Free Energy and corresponding secondary structure
(ss, mfe) = RNA.alifold(sequences)

# print output
print("{}\n{} [ {:6.2f} ]".format(cons, ss, mfe))


# The RNA sequence
seq = "GAGUAGUGGAACCAGGCUAUGUUUGUGACUCGCAGACUAACA"

# create a new model details structure
md = RNA.md()

# change temperature and dangle model
md.temperature = 20.0  # 20 Deg Celcius
md.dangles = 1  # Dangle Model 1

# create a fold compound
fc = RNA.fold_compound(seq, md)

# predict Minmum Free Energy and corresponding secondary structure
(ss, mfe) = fc.mfe()

# print sequence, structure and MFE
print("{}\n{} [ {:6.2f} ]".format(seq, ss, mfe))
