#include <glog/gbgraph.h>
#include <glog/gblegacysegment.h>
#include <glog/gbsegmentinserter.h>

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
        std::vector<size_t> &bodyNodeIdxs,
        bool edbCheck,
        bool &retainFree) {
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
    const bool isHeadUnary = h.getTupleSize() == 1;

#ifdef DEBUG
    LOG(DEBUGL) << "Original rule " << rule.tostring(program, layer);
#endif

    //Create a conjunctive query for the node we would like to add
    std::vector<Literal> outputQueryBody;
    std::vector<size_t> rangesOutputQueryBody;
    uint32_t outputCounter = counterFreshVarsQueryCont;
    const auto outputQueryHead = GBGraph::GBGraph_Node::createQueryFromNode(
            outputCounter,
            outputQueryBody,
            rangesOutputQueryBody,
            rule,
            bodyNodeIdxs,
            *this);

#ifdef DEBUG
    std::string query = "";
    for(auto &l : outputQueryBody) {
        query += " " + l.tostring(program, layer);
    }
    LOG(DEBUGL) << "Checking redundacy for QUERY (H) " << outputQueryHead->
        tostring(program, layer) << " " << query;
#endif

    retainFree = true; //by default, I assume that checking for duplicates
    //will not be necessary
    for(auto &nodeId : getNodeIDsWithPredicate(predId)) {
        assert(!isTmpNode(nodeId));
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
        LOG(DEBUGL) << "COMPARING against " << headNodeLiteral.tostring(
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
            LOG(DEBUGL) << "Found MATCH for " << headNodeLiteral.tostring(
                    program, layer) << match;
#endif
            std::chrono::duration<double, std::milli> dur =
                std::chrono::steady_clock::now() - start;
            if (isHeadUnary)
                durationQueryContain1 += dur;
            else
                durationQueryContain2 += dur;
            durationQueryContain += dur;
            retainFree = true;
            return true;
        }

        bool rt = false;
        std::chrono::steady_clock::time_point startE =
            std::chrono::steady_clock::now();
        bool responseEDBCheck = edbCheck &&
            isRedundant_checkEquivalenceEDBAtoms(
                    rt,
                    bodyNodeIdxs,
                    h,
                    body,
                    outputQueryHead.get(),
                    outputQueryBody,
                    rangesOutputQueryBody,
                    nodeId);
        std::chrono::duration<double, std::milli> durE =
            std::chrono::steady_clock::now() - startE;
        durationEDBCheck += durE;

        if (responseEDBCheck) {
            retainFree = true;
            std::chrono::duration<double, std::milli> dur =
                std::chrono::steady_clock::now() - start;
            if (isHeadUnary)
                durationQueryContain1 += dur;
            else
                durationQueryContain2 += dur;
            durationQueryContain += dur;
            return true;
        }
        if (!rt) {
            retainFree = false;
        }
    }
    std::chrono::duration<double, std::milli> dur =
        std::chrono::steady_clock::now() - start;
    if (isHeadUnary)
        durationQueryContain1 += dur;
    else
        durationQueryContain2 += dur;
    durationQueryContain += dur;
    return false;
}

struct __coord {
    size_t bodyAtomIdx;
    int posVarInLiteral;
    int posVarInLiteral2;
    int posVar; //Not used
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
    bool suitableForReplacement = false;
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
                    if (l.getTupleSize() == 1) {
                        suitableForReplacement = true;
                    }
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
        std::sort(termsToLookup.begin(), termsToLookup.end());
        p.nhits = nodeData->countHits(termsToLookup, 0);
        p.probedhits = termsToLookup.size();
    }
    std::sort(bodyAtomsWithHeadVar.begin(), bodyAtomsWithHeadVar.end());
    size_t selectedBodyAtomIdx = bodyAtomsWithHeadVar[0].bodyAtomIdx;
    size_t selectedPosInLiteral = bodyAtomsWithHeadVar[0].posVarInLiteral;

    if (bodyAtomsWithHeadVar[0].nhits < 18 || (!suitableForReplacement &&
                bodyAtomsWithHeadVar[0].nhits !=
                bodyAtomsWithHeadVar[0].probedhits)) {
        retainFree = false;
        return false;
    }

    start = std::chrono::steady_clock::now();

    //New tuples
    auto newNodeData = getNodeData(bodyNodeIdxs[selectedBodyAtomIdx]);
    bool isExistingColumnEDB = nodeData->hasColumnarBackend() &&
        ((TGSegmentLegacy*)nodeData.get())->getColumn(0)->isEDB();
    bool isNewColumnEDB = newNodeData->hasColumnarBackend() &&
        ((TGSegmentLegacy*)newNodeData.get())->getColumn(selectedPosInLiteral)->isEDB();
    bool fullProv = provenanceType == FULLPROV;

    //Here I check a EDB column (which we should process) against an existing
    //node which is not a EDB predicate
    //I don't need to store the node
    std::unique_ptr<GBSegmentInserter> retainedTermsIns;
    std::vector<Term_t> retainedTerms;
    if (!fullProv && isNewColumnEDB && isExistingColumnEDB) {
        isRedundant_checkEquivalenceEDBAtoms_one_edb_edb(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0,
                !suitableForReplacement);
    } else if (!fullProv && isExistingColumnEDB) {
        isRedundant_checkEquivalenceEDBAtoms_one_mem_edb(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0,
                !suitableForReplacement);
    } else if (!fullProv && isNewColumnEDB) {
        isRedundant_checkEquivalenceEDBAtoms_one_edb_mem(retainedTerms,
                newNodeData, selectedPosInLiteral, nodeData, 0,
                !suitableForReplacement);
    } else {
        if (!fullProv) {
            isRedundant_checkEquivalenceEDBAtoms_one_mem_mem(retainedTerms,
                    newNodeData, selectedPosInLiteral, nodeData, 0,
                    !suitableForReplacement);
        } else {
            retainedTermsIns = GBSegmentInserter::getInserter(
                    newNodeData->getNColumns() + 2, 2, false);
            isRedundant_checkEquivalenceEDBAtoms_one_mem_mem(retainedTermsIns,
                    newNodeData, selectedPosInLiteral, nodeData, 0,
                    !suitableForReplacement);
        }
    }

    if ((!fullProv && retainedTerms.empty()) || (fullProv && retainedTermsIns->isEmpty())) {
        retainFree = true;
        return true;
    } else if (!suitableForReplacement) {
        retainFree = false;
        return false;
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

        if (idxBodyAtom != -1 && !isTmpNode(bodyNodeIdxs[idxBodyAtom])) {
            const Literal &b = originalRuleBody[idxBodyAtom]; //This is
            //the atom the we should consider for the replacement.
            assert(idxBodyAtom < bodyNodeIdxs.size());
            size_t nodeToReplace = bodyNodeIdxs[idxBodyAtom];
            assert(!isTmpNode(nodeToReplace));
            assert(nodeToReplace != ~0ul);
            //2: Create a temporary node with only the facts that
            //lead to new derivations
            std::shared_ptr<const TGSegment> d;
            if (!fullProv) {
                d = std::shared_ptr<const TGSegment>(new UnaryWithConstProvTGSegment(
                            retainedTerms, nodeToReplace, true, 0));
            } else {
                d = retainedTermsIns->getSegment(nodeToReplace, true, 0,
                        newNodeData->getProvenanceType(), 2);
            }
            auto newNodeId = addTmpNode(b.getPredicate().getId(), d);

            //3: Use the temporary node instead
            bodyNodeIdxs[idxBodyAtom] = newNodeId;
            retainFree = true;
        } else if (idxBodyAtom == -1) {
            LOG(ERRORL) << "Should not occur";
            throw 10;
        }
    }

    return false;
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_edb_edb(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld,
        bool stopAfterFirst) {

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
            posInL2,
            stopAfterFirst);
    retainedColumn = retainedValues[0];
    assert(retainedColumn->isBackedByVector());
    ((InmemoryColumn*)retainedColumn.get())->swap(out);
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_mem_edb(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld,
        bool stopAfterFirst) {

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
        int posOld,
        bool stopAfterFirst) {

    //Old (edb)
    auto newColNode = ((TGSegmentLegacy*)newSeg.get())->getColumn(posNew);
    assert(newColNode->isEDB());
    auto retainedValues = layer->checkNewIn(
            ((EDBColumn*)newColNode.get())->getLiteral(),
            ((EDBColumn*)newColNode.get())->posColumnInLiteral(),
            oldSeg);
    retainedValues.swap(out);
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_mem_mem(
        std::vector<Term_t> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld,
        bool stopAfterFirst) {
    assert(newSeg->getNColumns() == 1); //I suspect that there is a bug if newSeg has more than a column. In that case, we should store in out also the values of the other columns in out. This assert checks for the condition when this bug might become active. This bug may not exist if later there are other checks that somehow rewrite the position of the variables in the original node into variables in the rewritten ones.
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
            if (stopAfterFirst)
                break;
        }
    }
    if (vnew != ~0ul)
        out.push_back(vnew);
    while (itrNew->hasNext()) {
        itrNew->next();
        out.push_back(itrNew->get(posNew));
        if (stopAfterFirst) {
            break;
        }
    }
    std::sort(out.begin(), out.end());
    auto newend = std::unique(out.begin(), out.end());
    out.erase(newend, out.end());
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_one_mem_mem(
        std::unique_ptr<GBSegmentInserter> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld,
        bool stopAfterFirst) {
    assert(provenanceType == FULLPROV); //This method is implemented only for the FULLPROV scenario
    auto itrOld = oldSeg->iterator();
    const size_t ncols = newSeg->getNColumns();
    std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(new Term_t[ncols + 2]);
    //Node + offset to original node
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
    int64_t rowIdNew = -1; //When I store the offset, I store the link to the
    //row in the original node

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
                rowIdNew++;
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
            for(size_t i = 0; i < ncols; ++i)
                row[i] = itrNew->get(i);
            row[ncols] = itrNew->getNodeId();
            row[ncols+1] = rowIdNew;
            out->add(row.get());
            processedTerms++;
            vnew = ~0ul;
            if (stopAfterFirst)
                break;
        }
    }
    if (vnew != ~0ul) {
        for(size_t i = 0; i < ncols; ++i)
            row[i] = itrNew->get(i);
        row[ncols] = itrNew->getNodeId();
        row[ncols+1] = rowIdNew;
        out->add(row.get());
    }
    while (itrNew->hasNext()) {
        itrNew->next();
        rowIdNew++;
        for(size_t i = 0; i < ncols; ++i)
            row[i] = itrNew->get(i);
        row[ncols] = itrNew->getNodeId();
        row[ncols+1] = rowIdNew;
        out->add(row.get());
        if (stopAfterFirst) {
            break;
        }
    }
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_two_mem_mem(
        std::vector<std::pair<Term_t,Term_t>> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew1,
        int posNew2,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld1,
        int posOld2) {
    if (posNew1 != 0) {
        std::vector<uint8_t> pos;
        pos.push_back(posNew1);
        pos.push_back(posNew2);
        newSeg = newSeg->sortBy(pos);
    }
    auto itrNew = newSeg->iterator();
    auto itrOld = newSeg->iterator();

    Term_t vold1 = ~0ul;
    Term_t vold2 = ~0ul;
    Term_t vnew1 = ~0ul;
    Term_t vnew2 = ~0ul;
    while (true) {
        if (vold1 == ~0ul) {
            if (itrOld->hasNext()) {
                itrOld->next();
                vold1 = itrOld->get(0);
                vold2 = itrOld->get(1);
            } else {
                break;
            }
        }
        if (vnew1 == ~0ul) {
            if (itrNew->hasNext()) {
                itrNew->next();
                vnew1 = itrNew->get(posNew1);
                vnew2 = itrNew->get(posNew2);
            } else {
                break;
            }
        }
        if (vold1 == vnew1 && vold2 == vnew2) {
            vnew1 = vnew2 = ~0ul;
        } else if (vold1 < vnew1 || (vold1 == vnew1 && vold2 < vnew2)) {
            vold1 = vold2 = ~0ul;
        } else {
            out.push_back(std::make_pair(vnew1, vnew2));
            vnew1 = vnew2 = ~0ul;
        }
    }
    if (vnew1 != ~0ul)
        out.push_back(std::make_pair(vnew1, vnew2));
    while (itrNew->hasNext()) {
        itrNew->next();
        vnew1 = itrNew->get(posNew1);
        vnew2 = itrNew->get(posNew2);
        out.push_back(std::make_pair(vnew1, vnew2));
    }
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_two_mem_edb(
        std::vector<std::pair<Term_t,Term_t>> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew1,
        int posNew2,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld1,
        int posOld2) {

    //Old
    auto oldColNode1 = ((TGSegmentLegacy*)oldSeg.get())->getColumn(0);
    assert(oldColNode1->isEDB());
    auto oldColNode2 = ((TGSegmentLegacy*)oldSeg.get())->getColumn(1);
    assert(oldColNode2->isEDB());

    //Check the columns come from the same literal
    const Literal &l3 = ((EDBColumn*)oldColNode1.get())->getLiteral();
    const Literal &l4 = ((EDBColumn*)oldColNode2.get())->getLiteral();
    std::vector<Substitution> subs;
    if (!l3.sameVarSequenceAs(l3) ||
            l3.subsumes(subs, l3, l4) == -1) {
        //The columns come from different literals. This is not supported
        throw 10;
    }

    int posInL2_1 = ((EDBColumn*)oldColNode1.get())->posColumnInLiteral();
    int posInL2_2 = ((EDBColumn*)oldColNode2.get())->posColumnInLiteral();
    auto retainedValues = layer->checkNewIn(
            newSeg, posNew1, posNew2,
            ((EDBColumn*)oldColNode2.get())->getLiteral(),
            posInL2_1, posInL2_2);
    out.swap(retainedValues);
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_two_edb_mem(
        std::vector<std::pair<Term_t,Term_t>> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew1,
        int posNew2,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld1,
        int posOld2) {
    LOG(WARNL) << "Not implemented";
    throw 10;
}

void GBGraph::isRedundant_checkEquivalenceEDBAtoms_two_edb_edb(
        std::vector<std::pair<Term_t,Term_t>> &out,
        std::shared_ptr<const TGSegment> newSeg,
        int posNew1,
        int posNew2,
        std::shared_ptr<const TGSegment> oldSeg,
        int posOld1,
        int posOld2) {
    //New
    assert(newSeg->getNColumns() >= 2);
    auto newColNode1 = ((TGSegmentLegacy*)newSeg.get())->getColumn(posNew1);
    assert(newColNode1->isEDB());
    auto newColNode2 = ((TGSegmentLegacy*)newSeg.get())->getColumn(posNew2);
    assert(newColNode2->isEDB());

    //Check the columns come from the same literal
    const Literal &l1 = ((EDBColumn*)newColNode1.get())->getLiteral();
    const Literal &l2 = ((EDBColumn*)newColNode2.get())->getLiteral();
    std::vector<Substitution> subs;
    if (!l1.sameVarSequenceAs(l2) ||
            l1.subsumes(subs, l1, l2) == -1) {
        //The columns come from different literals. This is not supported
        throw 10;
    }

    std::vector<uint8_t> posInL1;
    posInL1.push_back(((EDBColumn*)newColNode1.get())->posColumnInLiteral());
    posInL1.push_back(((EDBColumn*)newColNode2.get())->posColumnInLiteral());

    //Old
    auto oldColNode1 = ((TGSegmentLegacy*)oldSeg.get())->getColumn(0);
    assert(oldColNode1->isEDB());
    auto oldColNode2 = ((TGSegmentLegacy*)oldSeg.get())->getColumn(1);
    assert(oldColNode2->isEDB());

    //Check the columns come from the same literal
    const Literal &l3 = ((EDBColumn*)oldColNode1.get())->getLiteral();
    const Literal &l4 = ((EDBColumn*)oldColNode2.get())->getLiteral();
    subs.clear();
    if (!l3.sameVarSequenceAs(l3) ||
            l3.subsumes(subs, l3, l4) == -1) {
        //The columns come from different literals. This is not supported
        throw 10;
    }

    std::vector<uint8_t> posInL2;
    posInL2.push_back(((EDBColumn*)oldColNode1.get())->posColumnInLiteral());
    posInL2.push_back(((EDBColumn*)oldColNode2.get())->posColumnInLiteral());

    std::shared_ptr<Column> retainedColumn1, retainedColumn2;
    auto retainedValues = layer->checkNewIn(
            ((EDBColumn*)newColNode1.get())->getLiteral(),
            posInL1,
            ((EDBColumn*)oldColNode2.get())->getLiteral(),
            posInL2);
    retainedColumn1 = retainedValues[0];
    retainedColumn2 = retainedValues[1];
    assert(retainedColumn1->isBackedByVector());
    assert(retainedColumn2->isBackedByVector());

    const auto &v1 = retainedColumn1->getVectorRef();
    const auto &v2 = retainedColumn2->getVectorRef();
    assert(v1.size() == v2.size());
    for(size_t i = 0; i < v1.size(); ++i) {
        out.push_back(std::make_pair(v1[i], v2[i]));
    }
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

    //Existing tuples
    auto nodeData = getNodeData(nodeId);
    //Ignore if node is too small
    if (nodeData->getNRows() < 1000) {
        retainFree = false;
        return false;
    }

    const uint32_t vId1 = originalRuleHead.getTermAtPos(0).getId();
    const uint32_t vId2 = originalRuleHead.getTermAtPos(1).getId();

    std::vector<__coord> bodyAtomsWithHeadVar;
    for(size_t i = 0; i < originalRuleBody.size(); ++i) {
        auto &l = originalRuleBody[i];
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
            __coord c;
            c.bodyAtomIdx = i;
            c.posVarInLiteral = pos1;
            c.posVarInLiteral2 = pos2;
            bodyAtomsWithHeadVar.push_back(c);
        }
    }
    if (bodyAtomsWithHeadVar.size() == 0)
        return false;

    //If we can select more body atoms, then we pick the one with the
    //smallest cardinality and highest hit ratio
    std::vector<std::pair<Term_t,Term_t>> termsToLookup;
    for(int i = 0; i < bodyAtomsWithHeadVar.size(); ++i) {
        termsToLookup.clear();
        auto &p = bodyAtomsWithHeadVar[i];
        auto nodeIdx = bodyNodeIdxs[p.bodyAtomIdx];
        p.card = getNodeSize(nodeIdx);

        //The numbers are too small to try it
        if (p.card < 1000) {
            retainFree = false;
            return false;
        }

        //std::chrono::steady_clock::time_point starth =
        //    std::chrono::steady_clock::now();
        //Try to join up to 20 elements
        auto itr = getNodeIterator(nodeIdx);
        size_t maxCount = 20;
        size_t currentCount = 0;
        while (itr->hasNext() && currentCount++ < maxCount) {
            itr->next();
            Term_t valueToLookup1 = itr->get(p.posVarInLiteral);
            Term_t valueToLookup2 = itr->get(p.posVarInLiteral2);
            termsToLookup.push_back(std::make_pair(valueToLookup1,
                        valueToLookup2));
        }
        std::sort(termsToLookup.begin(), termsToLookup.end());
        p.nhits = nodeData->countHits(termsToLookup, 0, 1);
        p.probedhits = termsToLookup.size();
    }
    std::sort(bodyAtomsWithHeadVar.begin(), bodyAtomsWithHeadVar.end());
    size_t selectedBodyAtomIdx = bodyAtomsWithHeadVar[0].bodyAtomIdx;
    size_t selectedPos1 = bodyAtomsWithHeadVar[0].posVarInLiteral;
    size_t selectedPos2 = bodyAtomsWithHeadVar[0].posVarInLiteral2;

    if (bodyAtomsWithHeadVar[0].nhits < 18) {
        //Too many terms are new (no hits). Avoid to do the check
        retainFree = false;
        return false;
    }

    //New
    auto newNodeData = getNodeData(bodyNodeIdxs[selectedBodyAtomIdx]);
    bool isExistingColumnEDB = nodeData->hasColumnarBackend() &&
        ((TGSegmentLegacy*)nodeData.get())->getColumn(0)->isEDB() &&
        ((TGSegmentLegacy*)nodeData.get())->getColumn(1)->isEDB();
    bool isNewColumnEDB = newNodeData->hasColumnarBackend() &&
        ((TGSegmentLegacy*)newNodeData.get())->getColumn(selectedPos1)->isEDB() &&
        ((TGSegmentLegacy*)newNodeData.get())->getColumn(selectedPos2)->isEDB();

    std::vector<std::pair<Term_t, Term_t>> retainedTerms;
    if (isNewColumnEDB && isExistingColumnEDB) {
        isRedundant_checkEquivalenceEDBAtoms_two_edb_edb(retainedTerms,
                newNodeData, selectedPos1, selectedPos2, nodeData, 0, 1);
    } else if (isExistingColumnEDB) {
        isRedundant_checkEquivalenceEDBAtoms_two_mem_edb(retainedTerms,
                newNodeData, selectedPos1, selectedPos2, nodeData, 0, 1);
    } else if (isNewColumnEDB) {
        isRedundant_checkEquivalenceEDBAtoms_two_edb_mem(retainedTerms,
                newNodeData, selectedPos1, selectedPos2, nodeData, 0, 1);
    } else {
        isRedundant_checkEquivalenceEDBAtoms_two_mem_mem(retainedTerms,
                newNodeData, selectedPos1, selectedPos2, nodeData, 0, 1);

    }

    //std::chrono::duration<double, std::milli> dur =
    //    std::chrono::steady_clock::now() - start;

    if (retainedTerms.empty()) {
        retainFree = true;
        return true;
    } else {
        retainFree = false;
        return false;
    }

    //If I decide to create a temporary node like with unary predicates and fullProv is activated,
    //then I must store also the offsets; hence retainedTerms is not sufficient.
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
