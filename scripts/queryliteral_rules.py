import requests
import sys
import json
import pprint

address = sys.argv[1]
predicate = sys.argv[2]
rules = sys.argv[3]
gbchase = "true"
rewriteProgram = "true"
limit = "10"

data = { "predicate" : predicate, "rules" : rules, "gbchase" : gbchase, 
        "rewriteProgram" : rewriteProgram, "limit" : limit 
        }
r = requests.post(address, data=data)
print(r.status_code, r.reason)
out = json.loads(r.text)
pprint.pprint(out)
