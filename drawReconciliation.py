from ete3 import Tree
from Tkinter import *
import math
import argparse

"""
t = Tree("(A,(B,(C,D)));")
t2 = Tree("(((A,B),(C,D)),E);")

host = Tree("((A,B)D,C)E;", format=1)
#guest = Tree("((a,b)e,(c,d)f)g;", format=1)
guest = Tree("((((a,b)f,c)g,d)h,e)i;", format=1)

nodemap = {} #guest -> host
simplemap = {'a':"A", 'b':"A", 'c':"B", 'd':"B",  'e':"C", 'f':"A", 'g':"D", 'h':"D", 'i':"E"}

for key in simplemap:
    nodemap[guest&key] = host&simplemap[key]
"""

def readMapping(host, guest, mapping):
    m = {}
    for name in mapping:
        m[guest&name] = host&mapping[name]
    return m

def levelOrder(t):
    levels = []
    level = [leaf for leaf in t]
    seen = set(level)
    while len(level) > 0:
        levels.append(level)
        level = list(set([node.up for node in level if node.up != None]))
        for i in range(len(level)-1, -1, -1):
            node = level[i]
            if node.children != [] and (node.children[0] not in seen or node.children[1] not in seen):
                level.pop(i).name
        for node in level:
            seen.add(node)
    return levels

class Main:
    def __init__(self, host, guest, nodemap):

        window = Tk()
        window.title("Reconciliation Viewer")
        self.width = 1800
        self.height = 1200
        self.canvas = Canvas(window, width = self.width, height = self.height, bg = "white")
        self.canvas.pack()

        self.host = host
        self.guest = guest
        self.nodemap = nodemap

        self.generatePoints()
        self.drawSpeciesTree2()
        self.generateGuestPoints()
        self.drawGuestTree()

        window.mainloop()

    def drawBox(self, point, name=''):
        x, y = point
        size = 5
        self.canvas.create_rectangle(x-size, y-size, x+size, y+size, outline="#000000", fill="#1f1", width=2)
        self.canvas.create_text((x, y-size-7), text=name)

    def generatePoints(self):

        levels = levelOrder(self.host)
        points = []

        for j in range(len(levels)):

            level = levels[j]
            ycoord = int(self.height / float(len(levels) + 1) * (j+1))
            y = self.height -  ycoord

            if j == 0:
                xcoords = [int(float(self.width) / (len(level) + 1) * i) for i in range(1, len(level) + 1)]
                for i in range(len(level)):
                    node = level[i]
                    x = xcoords[i]
                    node.add_feature("coord", (x, y))
                    points.append((x,y))
                    text_coord = (x, y+20)
                    self.canvas.create_text(text_coord, text=node.name, font=("purisa", 18))

            else:
                for i in range(len(level)):
                    node = level[i]
                    x = (node.children[0].coord[0] + node.children[1].coord[0]) / 2
                    node.add_feature("coord", (x, y))
                    points.append((x,y))

    #Tries to minimize overlap between branches
    def drawSpeciesTree(self, pipeWidth = 75, color = "#d4d4d4"):

        def eulerTour(t):
            tour.append((t, 1))
            if t.children != []:
                eulerTour(t.children[0])
                tour.append((t, 2))
                eulerTour(t.children[1])
                tour.append((t, 3))

        tour = []
        eulerTour(self.host)

        points = []

        for pair in tour:
            node, passnum = pair
            x, y = node.coord
            if passnum == 1:
                points.append(x - pipeWidth)
                points.append(y)
                if node.children == []:
                    points.append(x + pipeWidth)
                    points.append(y)
            elif passnum == 2:
                points.append(x)
                points.append(y + pipeWidth)
            elif passnum == 3:
                points.append(x + pipeWidth)
                points.append(y)

        self.canvas.create_polygon(points, fill = color)

    #Looks more natural, but more overlap between branches
    def drawSpeciesTree2(self, pipeWidth = 75, color = "#d4d4d4"):

        for node in self.host.traverse():
            if node.children == []:
                continue

            c1 = node.coord
            c2 = node.children[0].coord
            points = [c1[0] - pipeWidth, c1[1], c2[0] - pipeWidth, c2[1], c2[0] + pipeWidth, c2[1], c1[0] + pipeWidth, c1[1]]
            self.canvas.create_polygon(points, fill = color)

            c1 = node.coord
            c2 = node.children[1].coord
            points = [c1[0] - pipeWidth, c1[1], c2[0] - pipeWidth, c2[1], c2[0] + pipeWidth, c2[1], c1[0] + pipeWidth, c1[1]]
            self.canvas.create_polygon(points, fill = color)

        #Stub above the root
        c2 = self.host.coord
        c1 = list(self.host.coord)
        c1[1] -= 100
        points = [c1[0] - pipeWidth, c1[1], c2[0] - pipeWidth, c2[1], c2[0] + pipeWidth, c2[1], c1[0] + pipeWidth, c1[1]]
        self.canvas.create_polygon(points, fill = color)

    def drawGuestTree(self):
        for node in self.guest.traverse():
            if node.children != []:
                for child in node.children:
                    x,y = node.coord
                    a,b = child.coord
                    self.canvas.create_line(x,y,a,b)

        for node in self.guest.traverse():
            self.drawBox(node.coord, node.name)
  
    def generateGuestPoints(self, pipeWidth = 75):

        #Add in loss nodes
        for leaf in self.guest:
            node = leaf
            while node != self.guest:
                host_me = self.nodemap[node]
                host_parent = self.nodemap[node.up]
                if len(node.children) == 2:
                    lchild = self.nodemap[node.children[0]] == host_me and self.nodemap[node.children[1]] != host_me
                    rchild = self.nodemap[node.children[1]] == host_me and self.nodemap[node.children[0]] != host_me
                if len(node.children) == 2 and (lchild or rchild):
                    if self.nodemap[node.children[0]] == host_me:
                        tofix = node.children[1]
                    else:
                        tofix = node.children[0]
                    temp = Tree()
                    temp.name = "L_" + node.name
                    nodemap[temp] = host_me
                    temp.up = tofix.up
                    temp.children = [tofix]
                    tofix.up = temp
                    if tofix == node.children[0]:
                        node.children[0] = temp
                    else:
                        node.children[1] = temp
                if host_me != host_parent and host_me.up != host_parent:
                    #Add loss nodes in
                    dist = host_parent.get_distance(host_me, topology_only=True)
                    guest_parent = node.up
                    curr = node
                    for i in range(int(dist)):
                        temp = Tree()
                        temp.name = "L_" + str(i) + "_" + guest_parent.name
                        nodemap[temp] = nodemap[curr].up
                        temp.up = curr.up
                        temp.children = [curr]
                        curr.up = temp
                        if curr == guest_parent.children[0]:
                            guest_parent.children[0] = temp
                        else:
                            guest_parent.children[1] = temp
                        curr = temp
                    guest_parent = node
                else:
                    node = node.up

        #Add levels 
        for node in self.guest.traverse():
            node.add_feature('level', -1)

        for leaf in self.guest:
            node = leaf
            node.level = 0
            currmap = self.nodemap[node]
            currlevel = 0
            node = node.up
            while node != None:
                mymap = self.nodemap[node]
                if mymap == currmap:
                    node.level = max(node.level, currlevel + 1)
                else:
                    node.level = max(node.level, 0)
                currlevel = node.level
                currmap = mymap
                node = node.up

        #How many points at each level of a node in the host tree?
        rmap = {} #map of host -> guest
        for key in self.nodemap:
            rkey = self.nodemap[key]
            if rkey in rmap:
                rmap[rkey].append(key)
            else:
                rmap[rkey] = [key]

        hostlevels = {} # hostnode -> levelcounts
        usedlevels = {} # same as hostlevels, but will count how many of each level have been used so far
        for key in rmap:
            nodes = rmap[key]

            maxlevel = 0
            for node in nodes:
                maxlevel = max(maxlevel, node.level)

            levelsizes = [0 for _ in range(maxlevel+1)]
            for node in nodes:
                levelsizes[node.level] += 1

            hostlevels[key] = levelsizes
            usedlevels[key] = [0 for _ in range(maxlevel+1)]

        #Generate Points - this only works for generateSpeciesTree2
        for node in self.guest.traverse():
            hostnode, level = self.nodemap[node], node.level
            used = usedlevels[hostnode][level]
            maxlevel = len(usedlevels[hostnode])
            usedlevels[hostnode][level] += 1

            bottom = hostnode.coord
            if hostnode == self.host:
                top = list(self.host.coord)
                top[1] -= 100
            else:
                top = hostnode.up.coord

            ydiff = bottom[1] - top[1]
            yused = ydiff * level / maxlevel
            y = bottom[1] - yused

            xlow, xhigh = bottom[0], top[0]
            xmid = int(xlow + (xhigh - xlow) * (yused / float(ydiff)))
            xused = int(pipeWidth * 2 * (used + 1) / (float(hostlevels[hostnode][level]) + 1))
            x = xmid + xused - pipeWidth
            node.add_feature('coord', (x,y))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('host', type=str, help='The input host tree in newick format')
	parser.add_argument('guest', type=str, help='The input guest tree in newick format')
	parser.add_argument('mapping', type=str, help='Path to txt file containing guest->host mapping')

	args = parser.parse_args()

	host = Tree(args.host, format=1)
	guest = Tree(args.guest, format=1)

	nodemap = {}
	mapfile = open(args.mapping)
	for line in mapfile:
			gname, hname = line.strip().split('\t')
			nodemap[guest&gname] = host&hname

	Main(host, guest, nodemap)



"""
host = Tree('genes.stree')
guest = Tree('0.nt.raxml.treefix.tree')

#Add names
i = 0
for node in host.traverse():
    if node.name == '':
        node.name = str(i)
        i += 1

i = 0
for node in guest.traverse():
    if node.name == '':
        node.name = str(i)
        i += 1

nmap = list(open('nodemap.txt'))
temp = {}
for line in nmap:
    line = line.split('\t')
    temp[line[0]] = line[1].strip()
nodemap = readMapping(host, guest, temp)
Main(host, guest, nodemap)
"""

















