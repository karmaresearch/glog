#include <glog/gbquerier.h>

JSON GBQuerier::getDerivationTree(size_t nodeId, size_t factId) {
    JSON out;
    if (nodeId >= g.getNNodes() ||
            factId >= g.getNodeSize(nodeId)) {
        return out;
    }
    exportNode(out, nodeId, factId);
    return out;
}

JSON GBQuerier::getDerivationTree(
        std::shared_ptr<const TGSegment> data,
        size_t nodeId,
        size_t factId,
        PredId_t predId,
        size_t ruleIdx,
        size_t step,
        const std::vector<size_t> &incomingEdges)
{
    JSON out;
    exportNode(out, nodeId, factId, predId, data, ruleIdx, step, incomingEdges);
    return out;
}

Literal GBQuerier::getFact(PredId_t predId, std::shared_ptr<const
        TGSegment> data, size_t factId)
{
    auto pred = p.getPredicate(predId);
    auto card = pred.getCardinality();
    assert(data->getNColumns() == card);
    auto row = data->getRow(factId);
    VTuple t(card);
    for(size_t i = 0; i < card; ++i) {
        t.set(VTerm(0, row[i]), i);
    }
    Literal l(pred, t);
    return l;
}

std::string GBQuerier::getTupleIDs(Literal &l)
{
    std::string out = "[";
    for(size_t i = 0; i < l.getTupleSize(); ++i) {
        auto t = l.getTermAtPos(i);
        out += std::to_string(t.getValue());
        if (i < l.getTupleSize() - 1)
            out += ",";
    }
    out += "]";
    return out;
}

void GBQuerier::exportEDBNode(JSON &out, Literal &l, size_t factId)
{
    auto itr = this->l.getIterator(l);
    size_t i = 0;
    while (itr->hasNext()) {
        itr->next();
        if (i == factId) {
            auto predId = itr->getPredicateID();
            auto tuple = l.getTuple();
            for(size_t j = 0; j < tuple.getSize(); ++j) {
                tuple.set(VTerm(0, itr->getElementAt(j)), j);
            }
            auto fl = Literal(l.getPredicate(), tuple);
            auto sfl = fl.toprettystring(&p, &this->l);
            out.put("rule", "none");
            out.put("ruleIdx", "none");
            out.put("fact", sfl);
            out.put("nodeId", "none");
            out.put("step", "none");
            out.put("factId", factId);
            out.put("tupleIds", getTupleIDs(fl));

            break;
        }
        i++;
    }
}

void __convertStringTupleIDsIntoNumbers(
        std::string &in,
        std::vector<Term_t> &out) {
    //String of the form [x1, x2, ...]
    std::string n = "";
    for(size_t i = 0; i < in.size(); ++i) {
        switch (in[i]) {
            case '[':
            case ']':
                if (n.size() > 0) {
                    out.push_back(stoul(n));
                }
                break;
            case ',':
                if (n.size() > 0) {
                    out.push_back(stoul(n));
                }
                n = "";
                break;
            default:
                n += in[i];
        };
    }
}

bool GBQuerier::checkSoundnessDerivationTree(JSON &root)
{
    bool response = true;
    //Check soundness node
    auto sRuleIdx = root.get("ruleIdx");
    if (sRuleIdx == "none") {
        return true;
    }

    auto ruleIdx = stoul(sRuleIdx);
    if (ruleIdx == ~0ul) {
        return true;
    }
    auto sTupleIDs = root.get("tupleIds");
    std::vector<Term_t> tupleIDs;
    __convertStringTupleIDsIntoNumbers(sTupleIDs, tupleIDs);

    auto &rule = p.getRule(ruleIdx);
    std::map<size_t, Term_t> mappingVars2Consts;
    assert(rule.getHeads().size() == 1);
    auto litHead = rule.getHeads()[0];
    assert(tupleIDs.size() == litHead.getTupleSize());
    for(size_t i = 0; i < tupleIDs.size(); ++i) {
        auto t = litHead.getTermAtPos(i);
        if (t.isVariable()) {
            auto varId = t.getId();
            if (mappingVars2Consts.count(varId)) {
                if (tupleIDs[i] != mappingVars2Consts[varId]) {
                    throw 10;
                }
            } else {
                mappingVars2Consts.insert(std::make_pair(varId, tupleIDs[i]));
            }
        } else {
            if (tupleIDs[i] != t.getValue()) {
                throw 10;
            }
        }
    }

    if (root.containsChild("parents")) {
        auto bodyAtoms = rule.getBody();
        auto parents = root.getChild("parents");
        auto ps = parents.getListChildren();
        for (size_t i = 0; i < ps.size(); ++i) {
            auto p = ps[i];
            auto bodyAtom = bodyAtoms[i];
            sTupleIDs = p.get("tupleIds");
            tupleIDs.clear();
            __convertStringTupleIDsIntoNumbers(sTupleIDs, tupleIDs);
            assert(tupleIDs.size() == bodyAtom.getTupleSize());
            for(size_t i = 0; i < tupleIDs.size(); ++i) {
                auto t = bodyAtom.getTermAtPos(i);
                if (t.isVariable()) {
                    auto varId = t.getId();
                    if (mappingVars2Consts.count(varId)) {
                        if (tupleIDs[i] != mappingVars2Consts[varId]) {
                            throw 10;
                        }
                    } else {
                        mappingVars2Consts.insert(std::make_pair(varId, tupleIDs[i]));
                    }
                } else {
                    if (tupleIDs[i] != t.getValue()) {
                        throw 10;
                    }
                }
            }
            response = response & checkSoundnessDerivationTree(p);
        }
    } else {
        //Parents should be there, otherwise how did the rule fire?
        throw 10;
    }
    return response;
}

void GBQuerier::exportNode(JSON &out, size_t nodeId, size_t factId)
{
    auto data = g.getNodeData(nodeId);
    size_t ruleIdx = g.getNodeRuleIdx(nodeId);
    size_t step = g.getNodeStep(nodeId);
    auto nodePred = g.getNodePredicate(nodeId);
    auto incomingEdges = g.getNodeIncomingEdges(nodeId);
    return exportNode(out, nodeId, factId, nodePred, data, ruleIdx, step,
            incomingEdges);
}

void GBQuerier::exportNode(JSON &out,
        size_t nodeId, size_t factId,
        PredId_t nodePred,
        std::shared_ptr<const TGSegment> data,
        size_t ruleIdx,
        size_t step,
        const std::vector<size_t> &incomingEdges)
{
    std::string sRule = "";
    if (ruleIdx != ~0ul)
        sRule = p.getRule(ruleIdx).toprettystring(&p, &l);
    out.put("rule", sRule);
    out.put("ruleIdx", ruleIdx);
    out.put("step", step);
    out.put("nodeId", nodeId);
    out.put("factId", factId);

    //Construct the fact
    auto f = getFact(nodePred, data, factId);
    out.put("fact", f.toprettystring(&p, &l));
    out.put("tupleIds", getTupleIDs(f));

    //Add the parents
    JSON parents;
    if (ruleIdx != ~0ul) {
        auto ie = incomingEdges;
        auto row = data->getRow(factId);
        size_t begin = data->getNColumns() + 1; //Skip the node
        size_t end = data->getNColumns() + data->getNOffsetColumns();
        size_t j = 0;
        auto bodyLiterals = p.getRule(ruleIdx).getBody();
        for(size_t i = begin; i < end; ++i) {
            JSON parentNode;
            auto bodyLiteral = bodyLiterals[i - begin];
            auto offset = row[i];
            if (bodyLiteral.getPredicate().getType() == EDB) {
                exportEDBNode(parentNode, bodyLiteral, offset);
                if (j < ie.size() && ie[j] == ~0ul) {
                    j++;
                }
            } else {
                auto nodeId = ie[j];
                exportNode(parentNode, nodeId, offset);
                j++;
            }
            parents.push_back(parentNode);
        }
    }
    out.add_child("parents", parents);
}

std::vector<std::string> GBQuerier::getListPredicates() const {
    std::vector<std::string> out;
    auto predIds = g.getPredicateIDs();
    for(auto &predId : predIds) {
        std::string label = p.getPredicateName(predId);
        out.push_back(label);
    }
    return out;
}

JSON GBQuerier::getNodeDetailsWithPredicate(std::string predName) const
{
    JSON out;
    auto pred = p.getPredicate(predName);
    auto nodeIds = g.getNodeIDsWithPredicate(pred.getId());
    for(auto nodeId : nodeIds) {
        JSON nodeDetails;
        nodeDetails.put("id", nodeId);
        nodeDetails.put("ruleIdx", g.getNodeRuleIdx(nodeId));
        nodeDetails.put("n_facts", g.getNodeSize(nodeId));
        nodeDetails.put("step", g.getNodeStep(nodeId));
        out.push_back(nodeDetails);
    }
    return out;
}

JSON GBQuerier::getNodeFacts(size_t nodeId) const
{
    JSON out;
    auto data = g.getNodeData(nodeId);
    auto card = data->getNColumns();
    auto itr = data->iterator();
    while (itr->hasNext()) {
        itr->next();
        JSON row;
        for(size_t i = 0; i < card; ++i) {
            auto t = itr->get(i);
            auto str = l.getDictText(t);
            row.push_back(str);
        }
        out.push_back(row);
    }
    return out;
}


std::string GBQuerier::getTermText(Term_t t) const
{
    return l.getDictText(t);
}
