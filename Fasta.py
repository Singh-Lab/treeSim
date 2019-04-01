import RetrieveOrthogroup as RO
import TreeUtils

class FastaFile:
    """
    A FastaFile instance reads a fasta file and stores each sequence seen in the file,
    as well as the pubGeneID of each sequence and the species ID of the species it came
    from, and the domains in each sequence.
    """

    def __init__(self, fasta = [], hmmfile = ""):
        """
        Reads input from a fasta file. If unspecified, then an object with empty
        names, sequences, and speciesIDs lists is returned.

        Args:
            fasta (list): A list representation of a fasta file (with headers untouched)
        """
        if fasta == []:
            self.names, self.sequences, self.speciesIDs = [], [], []
        else:
            self.names, self.sequences, self.speciesIDs = RO.fastaToSeqs(fasta)

            if hmmfile != "":
                self.domains, self.domainStarts, self.domainEnds = [], [], []
                for i in range(len(self.sequences)):
                    a,b,c = TreeUtils.findDomains(self.sequences[i], hmmfile)
                    self.domainStarts.append(a)
                    self.domainEnds.append(b)
                    self.domains.append(c)

    #Basically another init, but can't override init :(
    def initFromGeneID(self, geneID, level, hmmfile):
        """
        Given a gene ID and an orthogroup level to work at, pulls the containing orthogroup
        from orthoDB and parses it.

        Args:
            geneID  (str ): The gene ID whose orthogroup to search for
            level   (str ): The NCBI Taxonomy ID level at which to search 
            hmmfile (str ): path to hmm model of domains in sequence
        """
        species = level
        oid = RO.findOrthogroupID(geneID, level, species)
        fasta = RO.getFasta(oid, level, species)
        self.names, self.sequences, self.speciesIDs = RO.fastaToSeqs(fasta)

        self.domains, self.domainStarts, self.domainEnds = [], [], []
        for i in range(len(self.sequences)):
            a,b,c = TreeUtils.findDomains(self.sequences[i], hmmfile)
            self.domainStarts.append(a)
            self.domainEnds.append(b)
            self.domains.append(c)

    #Getters (Setting is not allowed)
    def length(self):
        return len(self.sequences)

    def get(self, index):
        """Returns the ith entry from the object"""
        return (self.names[index], self.sequences[index], self.speciesIDs[index])
    def getName(self, index):
        return self.names[index]
        
    def getSequence(self, index):
        return self.sequences[index]
        
    def getSpeciesID(self, index):
        return self.speciesIDs[index]

    def getDomains(self, index):
        return self.domains[index]

    def getDomainStarts(self, index):
        return self.domainStarts[index]

    def getDomainEnds(self, index):
        return self.domainEnds[index]

    def getAllDomains(self):
        return self.domains

    def getAllSequences(self):
        return self.sequences

    def getAllSpeciesIDs(self):
        return self.speciesIDs

    def delete(self, index):
        """Deletes the ith entry from the object"""
        self.names.pop(index)
        self.sequences.pop(index)
        self.speciesIDs.pop(index)