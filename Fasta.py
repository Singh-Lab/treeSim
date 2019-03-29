import RetrieveOrthogroup as RO

class FastaFile:
    def __init__(self, fasta = []):
        """
        Reads input from a fasta file. If unspecified, then an object with empty
        names, sequences, and speciesIDs lists is returned.

        Args:
            fasta  (list): A list representation of a fasta file (with headers untouched)
        """
        if fasta == []:
            self.names, self.sequences, self.speciesIDs = [], [], []
        else:
            self.names, self.sequences, self.speciesIDs = RO.fastaToSeqs(fasta)

    def initFromFasta(self, geneID, level, species=-1):
        """
        Given a gene ID and an orthogroup level to work at, pulls the containing orthogroup
        from orthoDB and parses it.
        """
        if species == -1:
            species = level
        oid = RO.findOrthogroupID(geneID, level, species)
        fasta = RO.getFasta(oid, level, species)
        self.names, self.sequences, self.speciesIDs = RO.fastaToSeqs(fasta)

    def delete(self, index):
        """Deletes the ith entry from the object"""
        self.names.pop(index)
        self.sequences.pop(index)
        self.speciesIDs.pop(index)