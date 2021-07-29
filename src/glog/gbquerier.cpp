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

std::vector<Term_t> GBQuerier::getLeavesInDerivationTree(
        size_t nodeId,
        size_t factId,
        std::vector<Literal> &out)
{
    auto data = g.getNodeData(nodeId);
    size_t ruleIdx = g.getNodeRuleIdx(nodeId);
    size_t step = g.getNodeStep(nodeId);
    auto nodePred = g.getNodePredicate(nodeId);
    auto &incomingEdges = g.getNodeIncomingEdges(nodeId);
    getLeaves(nodeId, factId, nodePred, data, ruleIdx, step, incomingEdges, out);
    return data->getRow(factId);
}

void GBQuerier::getMappings(const Literal &l,
        const std::vector<uint64_t> &row,
        std::vector<std::pair<Term_t, Term_t>> &mappings)
{
    for(size_t i = 0; i < l.getTupleSize(); ++i)
    {
        auto t = l.getTermAtPos(i);
        if (t.isVariable())
        {
            auto varId = t.getId();
            bool exists = false;
            for (size_t j = 0; j < mappings.size(); ++j)
            {
                if (mappings[j].first == varId)
                {
                    assert(row.size() > i);
                    assert(row[i] == mappings[j].second);
                    exists = true;
                    break;
                }
            }
            if (!exists)
            {
                assert(row.size() > i);
                mappings.push_back(std::make_pair(varId, row[i]));
            }
        }
    }
}

Literal GBQuerier::ground(const Literal &l, const std::vector<
        std::pair<Term_t, Term_t>> &mappings,
        bool &ok)
{
    ok = true;
    auto tuple = l.getTuple();
    for(size_t j = 0; j < tuple.getSize(); ++j)
    {
        auto term = tuple.get(j);
        if (term.isVariable())
        {
            bool found = false;
            for(size_t i = 0; i < mappings.size(); ++i)
            {
                if (mappings[i].first == term.getId())
                {
                    found = true;
                    tuple.set(VTerm(0, mappings[i].second), j);
                    break;
                }
            }
            if (!found)
            {
                ok = false;
            }
        }
    }
    return Literal(l.getPredicate(), tuple);
}

void GBQuerier::getLeaves(
        size_t nodeId, size_t factId,
        PredId_t nodePred,
        std::shared_ptr<const TGSegment> data,
        size_t ruleIdx,
        size_t step,
        const std::vector<size_t> &incomingEdges,
        std::vector<Literal> &out)
{
    if (ruleIdx != ~0ul) {
        auto ie = incomingEdges;
        auto row = data->getRow(factId);
        size_t begin = data->getNColumns() + 1; //Skip the node
        size_t end = data->getNColumns() + data->getNOffsetColumns();
        size_t j = 0;
        auto bodyLiterals = p.getRule(ruleIdx).getBody();

        std::vector<std::pair<Term_t, Term_t>> mappings;
        getMappings(p.getRule(ruleIdx).getHeads()[0], row, mappings);

        for(size_t i = begin; i < end; ++i) {
            auto bodyLiteral = bodyLiterals[i - begin];
            auto offset = row[i];
            if (bodyLiteral.getPredicate().getType() == EDB) {
                bool isFullyGrounded = true;
                auto groundedAtom = ground(bodyLiteral, mappings, isFullyGrounded);
                if (isFullyGrounded)
                {
                    out.push_back(groundedAtom);
                } else {
                    exportEDBNode(groundedAtom, offset, out);
                }
                if (j < ie.size() && ie[j] == ~0ul) {
                    j++;
                }
            } else {
                auto nodeId = ie[j];
                assert(nodeId != ~0ul);
                std::vector<Term_t> row =
                    getLeavesInDerivationTree(nodeId, offset, out);
                getMappings(bodyLiteral, row, mappings);
                j++;
            }
        }
    } else {
        //TODO: I'm not sure what to do here, this happens with magic atoms
        throw 10;
    }
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
    this->l.releaseIterator(itr);
}

void GBQuerier::exportEDBNode(Literal &l, size_t factId, std::vector<Literal> &out)
{
    LOG(ERRORL) << "Slow method, must be fixed!";
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
            out.push_back(fl);
            break;
        }
        i++;
    }
    this->l.releaseIterator(itr);
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
        std::map<Var_t, Term_t> mappings;
        for(size_t i = begin; i < end; ++i) {
            JSON parentNode;
            auto bodyLiteral = bodyLiterals[i - begin];
            auto offset = row[i];
            if (bodyLiteral.getPredicate().getType() == EDB) {
                if (!l.isQueryAllowed(bodyLiteral)) {
                    //Copy all the existing variables into the literal
                    auto tuple = bodyLiteral.getTuple();
                    for(size_t j = 0; j < tuple.getSize(); ++j) {
                        auto t = tuple.get(j);
                        if (t.isVariable() && mappings.count(t.getId())) {
                            tuple.set(VTerm(0, mappings[t.getId()]), j);
                        }
                    }
                    Literal newLit = Literal(bodyLiteral.getPredicate(), tuple);
                    exportEDBNode(parentNode, newLit, offset);
                } else {
                    exportEDBNode(parentNode, bodyLiteral, offset);
                }
                if (j < ie.size() && ie[j] == ~0ul) {
                    j++;
                }
            } else {
                auto nodeId = ie[j];
                exportNode(parentNode, nodeId, offset);
                j++;
            }

            //Code to copy the current mappings from variables to ground terms
            auto sTupleIDs = parentNode.get("tupleIds");
            std::vector<Term_t> tupleIDs;
            __convertStringTupleIDsIntoNumbers(sTupleIDs, tupleIDs);
            for(size_t i = 0; i < bodyLiteral.getTupleSize(); ++i) {
                auto t = bodyLiteral.getTermAtPos(i);
                if (t.isVariable()) {
                    mappings[t.getId()] = tupleIDs[i];
                }
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
