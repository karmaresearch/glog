import sys

tgfile = sys.argv[1]

edbnodes = set()
nodes = set()
nedges = 0
maxdepth = 0
curdepth = 0
for line in open(tgfile, 'rt'):
    nedges += 1
    line = line[:-1]
    tkns = line.split('\t')
    tknid = int(tkns[0])
    nodes.add(tknid)
    tknid = int(tkns[2][5:])
    nodes.add(tknid)
    t1 = tkns[1]
    if t1.startswith("EDB"):
        tknid = t1[4:]
        edbnodes.add(int(tknid))
        if curdepth > maxdepth:
            maxdepth = curdepth
        curdepth = 1
    else:
        curdepth += 1
print("n. nodes", len(nodes))
print("n. edb nodes", len(edbnodes))
print("n. edges", nedges)
print("max depth", maxdepth)
