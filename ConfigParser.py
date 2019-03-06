class ConfigParser(object):
    pass

f = list(open('config.txt'))
for line in f:
    if line[0] == "#" or len(line) < 3:
        continue
    line = line.split()
    assert len(line) == 3
    assert line[1] == '='
    setattr(ConfigParser, line[0], line[-1])