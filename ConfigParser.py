from sys import exit

class ConfigParser(object):
    pass

f = list(open('config.txt'))
for line in f:
    if line[0] == "#" or len(line) < 3:
        continue
    line = line.split()
    if len(line) != 3 or line[1] != '=':
        print "Config File Parsing Error"
        exit(1)
    setattr(ConfigParser, line[0], line[-1])