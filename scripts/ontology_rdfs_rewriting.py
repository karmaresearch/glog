import sys
import gzip

onto_file = sys.argv[1]
rules_file = None
if len(sys.argv) > 2:
    rules_file = sys.argv[2]

classes = {}
invClasses = {}
properties = {}
invProperties = {}

if rules_file is not None:
    for line in open(rules_file, 'rt'):
        line = line[:-1]
        tkns = line.split(' :- ')
        predicateName = tkns[0][0:tkns[0].find('(')]
        body = tkns[1]
        if 'TE' in body:
            tkns = body.split(',')
            obj = tkns[2]
            if '<' in obj:
                # unary
                assert(tkns[1] == 'rdf:type')
                obj = obj[:-1]
                classes[predicateName] = obj
                invClasses[obj] = predicateName
            else:
                properties[predicateName] = tkns[1]
                invProperties[tkns[1]] = predicateName

# Read all the statements which have the predicate:
# subClassOf
# subPropertyOf
# range
# domain
# type
if onto_file.endswith(".gz"):
    fin = gzip.open(onto_file, "rt")
else:
    fin = open(onto_file, "rt")

def addClass(c):
    if c not in invClasses:
        idx = c.rfind("#")
        if idx == -1:
            idx = c.rfind("/")
        relname = c[idx+1:-1]
        relname = relname.upper()
        relname = relname.replace("_","")
        relname = relname.replace(":","")
        relname = relname.replace(".","")
        relname = relname.replace("(","")
        relname = relname.replace(")","")
        # Check that relname is unique
        while relname in classes:
            relname += "_1"
        classes[relname] = c
        invClasses[c] = relname

def addProperty(p):
    if p not in invProperties:
        idx = p.rfind("#")
        if idx == -1:
            idx = p.rfind("/")
        relname = p[idx+1:-1]
        relname = relname.upper()
        relname = relname.replace(":","")
        relname = relname.replace("_","")
        relname = relname.replace("(","")
        relname = relname.replace(")","")
        # Check that relname is unique
        while relname in properties:
            relname += "_1"
        properties[relname] = p
        invProperties[p] = relname

subClasses = {}
subProperties = {}
superProperties = {}
ranges = {}
domains = {}

def addSubClass(a,b):
    a_id = invClasses[a]
    b_id = invClasses[b]
    if a_id not in subClasses:
        subClasses[a_id] = []
    subClasses[a_id].append(b_id)
    if 'http://www.w3.org/2000/01/rdf-schema#' in a or 'http://www.w3.org/2000/01/rdf-schema#' in b:
        raise("Rewriting not possible!")

def addSubProperty(a,b):
    a_id = invProperties[a]
    b_id = invProperties[b]
    if a_id not in subProperties:
        subProperties[a_id] = []
    subProperties[a_id].append(b_id)
    if 'http://www.w3.org/2000/01/rdf-schema#' in a or 'http://www.w3.org/2000/01/rdf-schema#' in b:
        raise("Rewriting not possible!")
    if b_id not in superProperties:
        superProperties[b_id] = []
    superProperties[b_id].append(a_id)

def addRange(p,c):
    p_id = invProperties[p]
    c_id = invClasses[c]
    if p_id not in ranges:
        ranges[p_id] = []
    ranges[p_id].append(c_id)

def addDomain(p,c):
    p_id = invProperties[p]
    c_id = invClasses[c]
    if p_id not in domains:
        domains[p_id] = []
    domains[p_id].append(c_id)

def getAllSubproperties(p, subprops=set()):
    if p in superProperties:
        sp = superProperties[p]
        for subproperty in sp:
            if subproperty not in subprops:
                subprops.add(subproperty)
                subprops = getAllSubproperties(subproperty, subprops)
    else:
        subprops.add(p)
    return subprops


count = 0
for line in fin:
    if count > 0 and count % 1000000 == 0:
        print(count)
    count += 1
    try:
        line = line[:-1]
        tkns = line.split(' ')
        subj = tkns[0]
        pred = tkns[1]
        obj = ''
        for i in range(2, len(tkns)):
            obj += tkns[i]
        obj = obj[:-1]
        if pred == '<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>':
            addClass(obj)
        if pred == '<http://www.w3.org/2000/01/rdf-schema#subClassOf>':
            addClass(subj)
            addClass(obj)
            addSubClass(subj,obj)
        if pred == '<http://www.w3.org/2000/01/rdf-schema#range>':
            addClass(obj)
            addProperty(subj)
            addRange(subj, obj)
        if pred == '<http://www.w3.org/2000/01/rdf-schema#domain>':
            addClass(obj)
            addProperty(subj)
            addDomain(subj, obj)
        if pred == '<http://www.w3.org/2000/01/rdf-schema#subPropertyOf>':
            addProperty(subj)
            addProperty(obj)
            addSubProperty(subj,obj)
    except Exception as e:
        print(e)
        print("Ignored line", line)

EDBpred = 'TE'
# Write rules that get the data from the EDB layer
for k, v in classes.items():
    print("{}(X) :- {}(X,rdf:type,{})".format(k,EDBpred,v))

for k, v in properties.items():
    print("{}(X,Y) :- {}(X,{},Y)".format(k,EDBpred,v))

for k, values in subClasses.items():
    for v in values:
        print("{}(X) :- {}(X)".format(v,k))

for k, values in subProperties.items():
    for v in values:
        print("{}(X,Y) :- {}(X,Y)".format(v,k))

for k, values in domains.items():
    for v in values:
        # Get all subproperties
        sp = set()
        sp = getAllSubproperties(k, subprops=sp)
        for subprop in sp:
            print("{}(X) :- {}(X,Y)".format(v,subprop))

for k, values in ranges.items():
    for v in values:
        # Get all subproperties
        sp = set()
        sp = getAllSubproperties(k, subprops=sp)
        for subprop in sp:
            print("{}(Y) :- {}(X,Y)".format(v,subprop))



