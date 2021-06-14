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

Literal GBQuerier::getFact(PredId_t predId, std::shared_ptr<const
        TGSegment> data, size_t factId) {
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

std::string GBQuerier::getTupleIDs(Literal &l) {
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

void GBQuerier::exportEDBNode(JSON &out, Literal &l, size_t factId) {
    auto itr = this->l.getIterator(l);
    size_t i = 0;
    while (itr->hasNext()) {
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
        itr->next();
        i++;
    }
}

void GBQuerier::exportNode(JSON &out, size_t nodeId, size_t factId) {
    auto data = g.getNodeData(nodeId);
    size_t ruleIdx = g.getNodeRuleIdx(nodeId);

    std::string sRule = "";
    if (ruleIdx != ~0ul)
        sRule = p.getRule(ruleIdx).toprettystring(&p, &l);
    out.put("rule", sRule);
    out.put("ruleIdx", ruleIdx);
    out.put("step", g.getNodeStep(nodeId));
    out.put("nodeId", nodeId);
    out.put("factId", factId);

    //Construct the fact
    auto f = getFact(g.getNodePredicate(nodeId), data, factId);
    out.put("fact", f.toprettystring(&p, &l));
    out.put("tupleIds", getTupleIDs(f));

    //Add the parents
    JSON parents;
    if (ruleIdx != ~0ul) {
        auto ie = g.getNodeIncomingEdges(nodeId);
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

JSON GBQuerier::getNodeDetailsWithPredicate(std::string predName) const {
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

JSON GBQuerier::getNodeFacts(size_t nodeId) const {
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
