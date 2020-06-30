#include <vlog/gbchase/gbgraph.h>

bool GBGraph::isRedundant_checkTypeAtoms(const std::vector<Literal> &atoms) {
    for(size_t i = 1; i < atoms.size(); ++i) {
        auto typ = atoms[i].getPredicate().getType();
        auto prevtyp = atoms[i-1].getPredicate().getType();
        if (typ != prevtyp) {
            return false;
        }
    }
    return true;
}

bool GBGraph::isRedundant(size_t ruleIdx,
        std::vector<size_t> &bodyNodeIdxs, bool &retainFree) {
    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();

    //Get the rule
    const Rule &rule = allRules[ruleIdx];

    //Perform some checks
    if (rule.getHeads().size() != 1) {
        LOG(WARNL) << "Query containment does not work with multiple head"
            " atoms";
        return false;
    }

    const auto &body = rule.getBody();
    //Either all IDB or EDB
    if (!isRedundant_checkTypeAtoms(body)) {
        LOG(WARNL) << "The rule mixes EDB and IDB body atoms. Query "
            "containment does not support it";
        return false;
    }

    const auto &h = rule.getFirstHead();
    const auto predId = h.getPredicate().getId();
    if (!areNodesWithPredicate(predId)) {
        return false;
    }

#ifdef DEBUG
    LOG(INFOL) << "Original rule " << rule.tostring(program, layer);
#endif

    //Create a conjunctive query for the node we would like to add
    std::vector<Literal> outputQueryBody;
    std::vector<size_t> rangesOutputQueryBody;
    const auto outputQueryHead = createQueryFromNode(outputQueryBody,
            rangesOutputQueryBody, rule, bodyNodeIdxs, false);

#ifdef DEBUG
    std::string query = "";
    for(auto &l : outputQueryBody) {
        query += " " + l.tostring(program, layer);
    }
    LOG(INFOL) << "Checking redundacy for QUERY (H) " << outputQueryHead->
        tostring(program, layer) << " " << query;
#endif

    retainFree = true; //by default, I assume that checking for duplicates
    //will not be necessary
    for(auto &nodeId : getNodeIDsWithPredicate(predId)) {
        std::vector<Substitution> allSubs;
        //First check the head
        const auto &headNodeLiteral = getNodeHeadQuery(nodeId);
        auto nsubsHead = Literal::getSubstitutionsA2B(allSubs,
                headNodeLiteral,
                *outputQueryHead.get());
        if (nsubsHead == -1) {
            continue;
        }

        const auto &bodyNodeLiterals = getNodeBodyQuery(nodeId);

#ifdef DEBUG
        std::string match = "";
        for(auto &l : bodyNodeLiterals) {
            match += " " + l.tostring(program, layer);
        }
        LOG(INFOL) << "COMPARING against " << headNodeLiteral.tostring(
                program, layer) << match;
#endif

        //First perform query containment
        bool hr = true;
        for(const auto &bodyNodeLiteral : bodyNodeLiterals) {
            bool found = false;
            for(const auto &targetLiteral : outputQueryBody) {
                std::vector<Substitution> subs;
                auto nsubs = Literal::getSubstitutionsA2B(subs, bodyNodeLiteral,
                        targetLiteral);
                if (nsubs != -1) {
                    //Now I need to check that the subs are compatible with the
                    //other ones
                    bool isCompatible = true;
                    for (const auto &s : subs) {
                        //Does the variable already exist?
                        for(const auto &s2 : allSubs) {
                            if (s.origin == s2.origin && s.destination !=
                                    s2.destination) {
                                isCompatible = false;
                                break;
                            }
                        }
                    }
                    if (isCompatible) {
                        for(const auto &s : subs) {
                            allSubs.push_back(s);
                        }
                        found = true;
                        break;
                    }
                }
            }

            if (!found) {
                hr = false;
                break;
            }
        }
        if (hr) {
#ifdef DEBUG
            std::string match = "";
            for(auto &l : bodyNodeLiterals) {
                match += " " + l.tostring(program, layer);
            }
            LOG(INFOL) << "Found MATCH for " << headNodeLiteral.tostring(
                    program, layer) << match;
#endif
            std::chrono::duration<double, std::milli> dur =
                std::chrono::steady_clock::now() - start;
            //LOG(WARNL) << "2) Rule " << ruleIdx << " DUR: " << dur.count();
            durationQueryContain += dur;
            retainFree = true;
            return true;
        }

        bool rt = false;
        if (isRedundant_checkEquivalenceEDBAtoms(
                    rt,
                    bodyNodeIdxs,
                    h,
                    body,
                    outputQueryHead.get(),
                    outputQueryBody,
                    rangesOutputQueryBody,
                    nodeId)) {
            retainFree = true;
            std::chrono::duration<double, std::milli> dur =
                std::chrono::steady_clock::now() - start;
            //LOG(WARNL) << "1) Rule " << ruleIdx << " DUR: " << dur.count();
            durationQueryContain += dur;

            return true;
        }
        if (!rt) {
            retainFree = false;
        }
    }
    std::chrono::duration<double, std::milli> dur =
        std::chrono::steady_clock::now() - start;
    //LOG(WARNL) << "3) Rule " << ruleIdx << " DUR: " << dur.count();
    durationQueryContain += dur;
    return false;
}

struct __coord {
    size_t bodyAtomIdx;
    int posVarInLiteral;
    int posVar;
    size_t card;
    size_t nhits;
    size_t probedhits;
    __coord() : card(0), nhits(0), probedhits(0) {}
    bool operator < (const __coord& s2) const {
        size_t cost1 = nhits / probedhits;
        size_t cost2 = nhits / probedhits;
        if (cost1 < cost2) {
            return true;
        } else if (cost1 == cost2) {
            return card < s2.card;
        } else {
            return false;
        }
    }
};

bool GBGraph::isRedundant_checkEquivalenceEDBAtoms_one(
        bool &retainFree,
        std::vector<size_t> &bodyNodeIdxs,
        const Literal &originalRuleHead,
        const std::vector<Literal> &originalRuleBody,
        const Literal *rewrittenRuleHead,
        const std::vector<Literal> &rewrittenRuleBody,
        const std::vector<size_t> &rangeRewrittenRuleBody,
        const size_t nodeId) {

    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();

    const uint32_t vId = originalRuleHead.getTermAtPos(0).getId();
    std::vector<__coord> bodyAtomsWithHeadVar;
    for(size_t i = 0; i < originalRuleBody.size(); ++i) {
        auto &l = originalRuleBody[i];
        int varIdx = 0;
        for(size_t j = 0; j < l.getTupleSize(); ++j) {
            if (l.getTermAtPos(j).isVariable()) {
                if (l.getTermAtPos(j).getId() == vId) {
                    __coord c;
                    c.bodyAtomIdx = i;
                    c.posVarInLiteral = j;
                    c.posVar = varIdx;
                    bodyAtomsWithHeadVar.push_back(c);
                }
                varIdx++;
            }
        }
    }

    //Existing tuples
    auto nodeData = getNodeData(nodeId);

    //If we can select more body atoms, then we pick the one with the
    //smallest cardinality and highest hit ratio
    std::vector<Term_t> termsToLookup;
    for(int i = 0; i < bodyAtomsWithHeadVar.size(); ++i) {
        termsToLookup.clear();
        auto &p = bodyAtomsWithHeadVar[i];
        auto nodeIdx = bodyNodeIdxs[p.bodyAtomIdx];
        p.card = getNodeSize(nodeIdx);

        //The numbers are too small to try it
        if (p.card < 1000) {
            retainFree = false;
            //std::chrono::duration<double, std::milli> dur =
            //    std::chrono::steady_clock::now() - start;
            //LOG(WARNL) << "ONE 1 DUR: " << dur.count();
            return false;
        }

        //Try to join up to 20 elements
        auto itr = getNodeData(nodeIdx)->iterator();
        size_t maxCount = 20;
        size_t currentCount = 0;
        while (itr->hasNext() && currentCount++ < maxCount) {
            itr->next();
            Term_t valueToLookup = itr->get(p.posVarInLiteral);
            termsToLookup.push_back(valueToLookup);
        }
        p.nhits = nodeData->countHits(termsToLookup, 0);
        p.probedhits = termsToLookup.size();
    }
    std::sort(bodyAtomsWithHeadVar.begin(), bodyAtomsWithHeadVar.end());
    size_t selectedBodyAtomIdx = bodyAtomsWithHeadVar[0].bodyAtomIdx;
    size_t selectedPos = bodyAtomsWithHeadVar[0].posVar;
    size_t selectedPosInLiteral = bodyAtomsWithHeadVar[0].posVarInLiteral;

    if (bodyAtomsWithHeadVar[0].nhits < 18) {
        retainFree = false;
        return false;
    }

    start = std::chrono::steady_clock::now();

    //New tuples
    auto newNodeData = getNodeData(bodyNodeIdxs[selectedBodyAtomIdx]);

    //Here I check a EDB column (which we should process) against an existing
    //node which is not a EDB predicate
    std::vector<Term_t> retainedTerms;
    if (newNodeData->hasColumnarBackend() && nodeData->hasColumnarBackend()) {
        isRedundant_checkEquivalenceEDBAtoms_one_edb_edb(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0);
    } else if (nodeData->hasColumnarBackend()) {
        isRedundant_checkEquivalenceEDBAtoms_one_mem_edb(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0);
    } else if (newNodeData->hasColumnarBackend()) {
        isRedundant_checkEquivalenceEDBAtoms_one_edb_mem(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0);
    } else {
        isRedundant_checkEquivalenceEDBAtoms_one_mem_mem(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0);
    }

    if (retainedTerms.empty()) {
        retainFree = true;
        return true;
    } else {
        retainFree = false;

        //I should create a new temporary node here
        //1: Get the original body atom that was rewritten in the EDB
        //just considered
        int64_t idxBodyAtom = -1;
        for(size_t i = 0; i < originalRuleBody.size(); ++i) {
            auto &b = originalRuleBody[i];
            if (b.getTupleSize() == 1) {
                auto t = b.getTermAtPos(0);
                if (t.isVariable() && t.getId() == vId) {
                    idxBodyAtom = i;
                    break;
                }
            }
        }

        if (idxBodyAtom != -1) {
            const Literal &b = originalRuleBody[idxBodyAtom]; //This is
            //the atom the we should consider for the replacement.
            assert(idxBodyAtom < bodyNodeIdxs.size());
            size_t nodeToReplace = bodyNodeIdxs[idxBodyAtom];
            assert(nodeToReplace != ~0ul);
            //2: Create a temporary node with only the facts that
            //lead to new derivations
            std::shared_ptr<const TGSegment> d(new
                    UnaryWithConstProvTGSegment(
                        retainedTerms,
                        nodeToReplace, true, 0));
            auto newNodeId = addTmpNode(
                    b.getPredicate().getId(), d);

            //3: Use the temporary node instead
            bodyNodeIdxs[idxBodyAtom] = newNodeId;
            retainFree = true;
        }
    }

    return false;
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_edb_edb(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld) {

    //New
    auto newColNode = ((TGSegmentLegacy*)newSeg.get())->getColumn(posNew);
    assert(newColNode->isEDB());
    std::vector<uint8_t> posInL1;
    posInL1.push_back(((EDBColumn*)newColNode.get())->posColumnInLiteral());

    //Old
    auto oldColNode = ((TGSegmentLegacy*)oldSeg.get())->getColumn(0);
    assert(oldColNode->isEDB());
    std::vector<uint8_t> posInL2;
    posInL2.push_back(((EDBColumn*)oldColNode.get())->posColumnInLiteral());

    std::shared_ptr<Column> retainedColumn;
    auto retainedValues = layer->checkNewIn(
            ((EDBColumn*)newColNode.get())->getLiteral(),
            posInL1,
            ((EDBColumn*)oldColNode.get())->getLiteral(),
            posInL2);
    retainedColumn = retainedValues[0];
    assert(retainedColumn->isBackedByVector());
    ((InmemoryColumn*)retainedColumn.get())->swap(out);
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_mem_edb(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld) {

    //Old (edb)
    auto oldColNode = ((TGSegmentLegacy*)oldSeg.get())->getColumn(0);
    assert(oldColNode->isEDB());
    auto retainedValues = layer->checkNewIn(
            newSeg,
            posNew,
            ((EDBColumn*)oldColNode.get())->getLiteral(),
            ((EDBColumn*)oldColNode.get())->posColumnInLiteral());
    retainedValues.swap(out);
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_edb_mem(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld) {

    //Old (edb)
    auto newColNode = ((TGSegmentLegacy*)newSeg.get())->getColumn(posNew);
    assert(newColNode->isEDB());
    auto retainedValues = layer->checkNewIn(
            ((EDBColumn*)newColNode.get())->getLiteral(),
            ((EDBColumn*)newColNode.get())->posColumnInLiteral(),
            newSeg);
    retainedValues.swap(out);
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_mem_mem(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld) {
    auto itrOld = oldSeg->iterator();
    if (posNew != 0) {
        std::vector<uint8_t> pos;
        pos.push_back(posNew);
        newSeg = newSeg->sortBy(pos);
    }
    auto itrNew = newSeg->iterator();

    Term_t vold = ~0ul;
    Term_t vnew = ~0ul;

    size_t processedTerms = 0;
    size_t existingTerms = 0;

    while (true) {
        if (vold == ~0ul) {
            if (itrOld->hasNext()) {
                itrOld->next();
                vold = itrOld->get(0);
            } else {
                break;
            }
        }
        if (vnew == ~0ul) {
            if (itrNew->hasNext()) {
                itrNew->next();
                vnew = itrNew->get(posNew);
            } else {
                break;
            }
        }
        if (vold == vnew) {
            existingTerms++;
            processedTerms++;
            vnew = ~0ul;
        } else if (vold < vnew) {
            vold = ~0ul;
        } else {
            out.push_back(vnew);
            processedTerms++;
            vnew = ~0ul;
        }
    }
    if (vnew != ~0ul)
        out.push_back(vnew);
    while (itrNew->hasNext()) {
        itrNew->next();
        out.push_back(itrNew->get(posNew));
    }
    std::sort(out.begin(), out.end());
    auto newend = std::unique(out.begin(), out.end());
    out.erase(newend, out.end());
}

bool GBGraph::isRedundant_checkEquivalenceEDBAtoms_two(
        bool &retainFree,
        std::vector<size_t> &bodyNodeIdxs,
        const Literal &originalRuleHead,
        const std::vector<Literal> &originalRuleBody,
        const Literal *rewrittenRuleHead,
        const std::vector<Literal> &rewrittenRuleBody,
        const std::vector<size_t> &rangeRewrittenRuleBody,
        const size_t nodeId) {
    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    auto t1 = rewrittenRuleHead->getTermAtPos(0);
    auto t2 = rewrittenRuleHead->getTermAtPos(1);
    assert(t1.isVariable());
    assert(t2.isVariable());
    const uint32_t vId1 = t1.getId();
    const uint32_t vId2 = t2.getId();
    const uint32_t originalVId1 = originalRuleHead.getTermAtPos(0).getId();
    const uint32_t originalVId2 = originalRuleHead.getTermAtPos(1).getId();

    std::vector<
        std::pair<size_t, std::pair<uint8_t, uint8_t>>> bodyAtomsWithHeadVar;
    for(size_t i = 0; i < rewrittenRuleBody.size(); ++i) {
        auto &l = rewrittenRuleBody[i];
        int pos1 = -1;
        int pos2 = -1;
        for(size_t j = 0; j < l.getTupleSize(); ++j) {
            if (l.getTermAtPos(j).isVariable() &&
                    l.getTermAtPos(j).getId() == vId1) {
                pos1 = j;
            } else if (l.getTermAtPos(j).isVariable() &&
                    l.getTermAtPos(j).getId() == vId2) {
                pos2 = j;
            }
        }
        if (pos1 != -1 && pos2 != -1) {
            bodyAtomsWithHeadVar.push_back(
                    std::make_pair(i,std::make_pair(pos1, pos2)));
        }
    }
    if (bodyAtomsWithHeadVar.size() == 0)
        return false;

    //If we can select more body atoms, then we pick the one with the
    //smallest cardinality
    size_t selectedBodyAtomIdx = 0;
    size_t selectedPos1 = 0;
    size_t selectedPos2 = 0;
    size_t minCard = ~0ul;
    //I select the atom with the smallest cardinality
    for(int i = 0; i < bodyAtomsWithHeadVar.size(); ++i) {
        auto p = bodyAtomsWithHeadVar[i];
        auto c = layer->getCardinality(rewrittenRuleBody[p.first]);
        if (c < minCard) {
            minCard = c;
            selectedBodyAtomIdx = p.first;
            selectedPos1 = p.second.first;
            selectedPos2 = p.second.second;
        }
    }

    //Do a join between the bodyAtom and the node to see whether it's
    //redundant
    const auto &bl = rewrittenRuleBody[selectedBodyAtomIdx];
    auto nodeData = getNodeData(nodeId);
    assert(bl.getPredicate().getType() == EDB);

    //Check
    auto itrOld = nodeData->iterator();
    if (!itrOld->hasNext()) {
        LOG(ERRORL) << "Cannot be empty";
    }
    itrOld->next();
    Term_t vold1 = itrOld->get(0);
    Term_t vold2 = itrOld->get(1);

    auto itrNew = layer->getIterator(bl);
    itrNew->next();
    Term_t vnew1 = itrNew->getElementAt(selectedPos1);
    Term_t vnew2 = itrNew->getElementAt(selectedPos2);
    while (true) {
        if (vold1 == vnew1 && vold2 == vnew2) {
            if (itrNew->hasNext()) {
                itrNew->next();
                vnew1 = itrNew->getElementAt(selectedPos1);
                vnew2 = itrNew->getElementAt(selectedPos2);
            } else {
                layer->releaseIterator(itrNew);
                return true; // is redundant
            }
        } else if (vold1 < vnew1 || (vold1 == vnew1 && vold2 < vnew2)) {
            if (itrOld->hasNext()) {
                itrOld->next();
                vold1 = itrOld->get(0);
                vold2 = itrOld->get(1);
            } else {
                break;
            }
        } else {
            //One element in vnew cannot be found. break;
            break;
        }
    }
    layer->releaseIterator(itrNew);
    return false;
}

bool GBGraph::isRedundant_checkEquivalenceEDBAtoms(
        bool &retainFree,
        std::vector<size_t> &bodyNodeIdxs,
        const Literal &originalRuleHead,
        const std::vector<Literal> &originalRuleBody,
        const Literal *rewrittenRuleHead,
        const std::vector<Literal> &rewrittenRuleBody,
        const std::vector<size_t> &rangeRewrittenRuleBody,
        const size_t nodeId) {
    if (rewrittenRuleHead->getTupleSize() == 1) {
        return isRedundant_checkEquivalenceEDBAtoms_one(
                retainFree,
                bodyNodeIdxs,
                originalRuleHead,
                originalRuleBody,
                rewrittenRuleHead,
                rewrittenRuleBody,
                rangeRewrittenRuleBody,
                nodeId);
    } else if (rewrittenRuleHead->getTupleSize() == 2) {
        return isRedundant_checkEquivalenceEDBAtoms_two(
                retainFree,
                bodyNodeIdxs,
                originalRuleHead,
                originalRuleBody,
                rewrittenRuleHead,
                rewrittenRuleBody,
                rangeRewrittenRuleBody,
                nodeId);
    }
    retainFree = false;
    return false;
}
