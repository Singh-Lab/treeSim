from ConfigParser import ConfigParser as CP

#------------------------- PRELIMINARIES ------------------------#
#TODO: Change this so that the 4 hardcoded ones mean the right thing
#pylint: disable=no-member
pam = [i[:-1].split()[1:] for i in list(open(CP.PAM,'r'))][1:] 
aa = list(open(CP.PAM,'r'))[0].split()

b2i = {}
for i in range(len(aa)):
    b2i[aa[i]] = i

b2i['-'] = b2i['X']
b2i['B'] = b2i['D']
b2i['J'] = b2i['L']
b2i['Z'] = b2i['E']

important_positions = [0,1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21,22]

#----------------------- END PRELIMINARIES -----------------------#

def domainSim(a,b):
    """
    Returns the normalized similarity between two zinc finger domain instances. 
    Both must be length 23. Uses the PAM matrix
    """

    a = a.upper()
    b = b.upper()

    def simHelper(a,b):
        similarity = 0
        #Only works for zinc fingers
        for i in important_positions:
            x = b2i[a[i]]
            y = b2i[b[i]]
            similarity += int(pam[x][y]) #/ entropy[i]
        return similarity

    return simHelper(a, b) / float(max(simHelper(a,a),simHelper(b,b)))

def sequenceSim(a,b):
    """
    Computes pairwise similarity between two general AA sequences 
    a and b using the PAM matrix.
    """

    def simHelper(a,b):
        similarity = 0
        for i in range(len(a)):
            x = b2i[a[i]]
            y = b2i[b[i]]
            similarity += int(pam[x][y]) #TODO: Why is this an int?
        return similarity
    
    return simHelper(a,b) / float(max(simHelper(a,a), simHelper(b,b)))