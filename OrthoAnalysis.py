#Contains functions for manipulating and analyzing orthogroup files
import os
from TreeUtils import findDomains

#Checks if the input string is a valid zf-C2H2 domain
def isValid(domain):
	valid = len(domain) == 23 and domain[2] == "C" and domain[5] == "C"
	valid &= domain[18] == "H" and domain[22] == "H"
	return valid

def groupDomains(names, sequences, hmmfile):
	"""
	Takes a list of input sequences and returns a list of domain strings for each.
	Leaves an empty string at position i of the jth list if the jth sequence does 
	not have a copy of domain i. This aligns all existing domains and makes it clear
	which domains are present in which sequence

	Example (domains marked as xxx):

	sequences = ["AAxxxAAxxxAAAAAAA",
				 "AAAAAAAxxxAAxxxAA",
				 "AAxxxAAAAAAAxxxAA"]


	grouped = [[dom1, dom2, ''  ],
			  [''  , dom2, dom3],
			  [dom1, ''  , dom3]]

	Args:
		sequences (list): A list of sequences 
		hmmfile   (str ): The name of the hmmfile containing the desired domain model

	Returns:
		grouped   (list): A list of lists of all domain sequences from each domain
		domNames  (list): A list of lists of domain names for each domain in each sequence
	"""

	domStarts = [sorted(findDomains(i, hmmfile)[0]) for i in sequences]
	domNames = []
	allStarts = sorted(list(set.union(*[set(i) for i in domStarts])))
	grouped = []
	for i in range(len(domStarts)):
		domains = ['' for _ in range(len(allStarts))]
		dnames = ['' for _ in range(len(allStarts))]
		for start in domStarts[i]:
			domains[allStarts.index(start)] = sequences[i][start: start+23]
			dnames[allStarts.index(start)] = names[i] + "_" + str(start)
		grouped.append(domains)
		domNames.append(dnames)
	return grouped, domNames