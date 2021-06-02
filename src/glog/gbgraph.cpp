#include <glog/gbgraph.h>
#include <glog/gbquerier.h>
#include <glog/gbruleexecutor.h>
#include <glog/gbcompositesegment.h>

#include <vlog/support.h>


const Literal &GBGraph::GBGraph_Node::getQueryHead(GBGraph &g)
{
    if (!queryCreated) {
        createQueryFromNode(g);
    }
    return *(queryHead.get());
}

const std::vector<Literal> &GBGraph::GBGraph_Node::getQueryBody(GBGraph &g)
{
    if (!queryCreated) {
        createQueryFromNode(g);
    }
    return queryBody;
}

const std::vector<size_t> &GBGraph::GBGraph_Node::getIncomingEdges(bool c) const
{
#ifdef DEBUG
    if (c) {
        for(auto n : incomingEdges)
            if (n == ~0ul)
                throw 10;
    }
#endif
    return incomingEdges;
}

const std::vector<PredId_t> GBGraph::getPredicateIDs() const
{
    std::vector<PredId_t> out;
    for(auto &p : pred2Nodes)
        out.push_back(p.first);
    return out;
}

void GBGraph::GBGraph_Node::createQueryFromNode(GBGraph &g)
{
    auto newHead = GBGraph_Node::createQueryFromNode(
            g.counterFreshVarsQueryCont,
            queryBody,
            rangeQueryBody,
            g.allRules[ruleIdx],
            incomingEdges,
            g);
    this->queryHead = std::move(newHead);
    queryCreated = true;
}

std::unique_ptr<Literal> GBGraph::GBGraph_Node::createQueryFromNode(
        uint32_t &counterFreshVars,
        std::vector<Literal> &outputQueryBody,
        std::vector<size_t> &rangeOutputQueryBody,
        const Rule &rule,
        const std::vector<size_t> &incomingEdges,
        GBGraph &g)
{
    assert(rule.getHeads().size() == 1);
    //Rewrite the rule with fresh variables
    counterFreshVars++;
    Rule newRule = rule.rewriteWithFreshVars(counterFreshVars);
    auto &newBody = newRule.getBody();

    //const Rule &newRule = rule;
    //auto &newBody = rule.getBody();

    outputQueryBody.clear();
    int idxIncomingEdge = 0;
    for(int i = 0; i < newBody.size(); ++i) {
        if (i > 0) {
            rangeOutputQueryBody.push_back(outputQueryBody.size());
        }
        auto &l = newBody[i];
        if (l.getPredicate().getType() == EDB) {
            outputQueryBody.push_back(l);
        } else {
            assert(incomingEdges.size() > 0);
            size_t incEdge = incomingEdges[idxIncomingEdge++];
            //Get the head atom associated to the node
            const auto &litIncEdge = g.getNodeHeadQuery(incEdge);
            //Compute the MGU
            std::vector<Substitution> subs;
            auto nsubs = Literal::getSubstitutionsA2B(subs, litIncEdge, l);
            assert(nsubs != -1);

            //First replace all remaining variables with fresh ones
            const auto &bodyIncEdge = g.getNodeBodyQuery(incEdge);
            std::set<uint32_t> av;
            for(auto &litBody : bodyIncEdge) {
                for(size_t i = 0; i < litBody.getTupleSize(); ++i) {
                    auto t = litBody.getTermAtPos(i);
                    if (t.isVariable())
                        av.insert(t.getId());
                }
            }
            for(auto v : av) {
                bool found = false;
                for(auto &s : subs) {
                    if (s.origin == v) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    subs.push_back(Substitution(v, VTerm(
                                    counterFreshVars++, 0)));
                }
            }

            for(auto &litBody : bodyIncEdge) {
                outputQueryBody.push_back(litBody.substitutes(subs));
            }
        }
    }
    assert(idxIncomingEdge == incomingEdges.size());
    return std::unique_ptr<Literal>(new Literal(newRule.getFirstHead()));
}

void GBGraph::addNode(PredId_t predId,
        size_t step,
        std::vector<std::vector<std::string>> &facts) {

    //Convert the facts
    size_t nExtraColumns = 0;
    if (shouldTrackProvenance()) {
        nExtraColumns++;
        if (provenanceType == FULLPROV)
            nExtraColumns++;
    }
    size_t cardPred = program->getPredicate(predId).getCardinality();
    size_t rowSize = cardPred + nExtraColumns;
    auto ins = GBSegmentInserter::getInserter(rowSize, nExtraColumns, false);
    std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(
            new Term_t[rowSize]);

    auto nodeId = getNNodes();
    row[cardPred] = nodeId;
    for(size_t i = 1; i < nExtraColumns; ++i)
        row[cardPred + i] = ~0ul;
    for(auto &fact : facts) {
        assert(fact.size() == cardPred);
        for(size_t j = 0; j < fact.size(); ++j) {
            Term_t id;
            layer->getOrAddDictNumber(fact[j].c_str(), fact[j].size(), id);
            row[j] = id;
        }
        ins->add(row.get());
    }
    auto data = ins->getSegment(nodeId, false, 0, getSegProvenanceType(false),
            nExtraColumns);
    auto sortedData = data->sort();

    //Add the node
    if (!shouldTrackProvenance()) {
        addNodeNoProv(predId, ~0ul, step, sortedData);
    } else {
        std::vector<size_t> incomingEdges;
        addNodeProv(predId, ~0ul, step, sortedData, incomingEdges);
    }
}

void GBGraph::addNodeNoProv(PredId_t predId,
        size_t ruleIdx,
        size_t step,
        std::shared_ptr<const TGSegment> data) {
    if (shouldTrackProvenance()) {
        LOG(ERRORL) << "This method should not be called if provenance is"
            " activated";
        throw 10;
    }

    if (queryContEnabled) {
        LOG(ERRORL) << "Query containment is available only if provenance "
            "is activated";
        throw 10;
    }

#ifdef DEBUG
    if (data->getName() == "CompositeTGSegment") {
        LOG(ERRORL) << "CompositeTGSegment is meant to be used only during the"
            " execution of a single rule. It should not be added to the graph";
        throw 10;
    }
#endif

    auto nodeId = getNNodes();
    nodes.emplace_back();
    GBGraph_Node &outputNode = nodes.back();
    outputNode.predid = predId;
    outputNode.ruleIdx = ruleIdx;
    outputNode.step = step;
    outputNode.setData(data);

    pred2Nodes[predId].push_back(nodeId);
    LOG(DEBUGL) << "Added node ID " << nodeId << " with # facts=" <<
        data->getNRows();
}

void GBGraph::addNodesProv(PredId_t predId,
        size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> seg,
        const std::vector<std::shared_ptr<Column>> &provenance,
        bool filterFromProvenance) {

#ifdef DEBUG
    if (provenanceType == ProvenanceType::NOPROV) {
        LOG(ERRORL) << "This method should be called only if type provenance is"
            " not set to NOPROV";
        throw 10;
    }
#endif

    if (provenance.size() == 0) {
        //Single EDB or IDB atom
        if (!seg->isNodeConstant()) {
            //Must split the nodes
            auto resortedSeg = seg->sortByProv();

            std::vector<size_t> provNodes;
            auto chunks = resortedSeg->sliceByNodes(getNNodes(),
                    provNodes);
            assert(chunks.size() == provNodes.size());
            std::vector<size_t> currentNodeList(1);
            for (size_t i = 0; i < chunks.size(); ++i) {
                auto c = chunks[i];
                currentNodeList[0] = provNodes[i];
                if (currentNodeList[0] == ~0ul) {
                    assert(chunks.size() == 1);
                    addNodeProv(predId, ruleIdx, step, c,
                            std::vector<size_t>());
                } else {
                    if (filterFromProvenance) {
                        auto retainedData = retainFromDerivationTree(predId,
                                c, currentNodeList);
                        if (retainedData->getNRows() > 0) {
                            addNodeProv(predId, ruleIdx, step, retainedData,
                                    currentNodeList);
                        }
                    } else {
                        addNodeProv(predId, ruleIdx, step, c, currentNodeList);
                    }
                }
            }
        } else {
            std::vector<size_t> provnodes;
            if (seg->getNodeId() == ~0ul) {
                //EDB body atom
#if DEBUG
                assert(allRules[ruleIdx].getBody().size() == 1 &&
                        allRules[ruleIdx].getBody()[0].getPredicate().getType() ==
                        EDB);
#endif
            } else {
                provnodes.push_back(seg->getNodeId());
            }
            auto nodeId = getNNodes();
            auto dataToAdd = seg->slice(nodeId, 0, seg->getNRows());
            if (filterFromProvenance) {
                auto retainedData = retainFromDerivationTree(predId, dataToAdd,
                        provnodes);
                if (retainedData->getNRows() > 0) {
                    addNodeProv(predId, ruleIdx, step, retainedData, provnodes);
                }
            } else {
                addNodeProv(predId, ruleIdx, step, dataToAdd, provnodes);
            }
        }
    } else {
        const auto nnodes = (provenance.size() + 2) / 2;
        const auto nrows = seg->getNRows();
        std::vector<size_t> provnodes(nrows * nnodes);

        auto itr = seg->iterator();
        size_t i = 0;
        while (itr->hasNext()) {
            itr->next();
            size_t provRowIdx = itr->getNodeId();
            for(int j = nnodes - 1; j >= 0; j--) {
                if (j == 0) {
                    provnodes[i * nnodes] = provenance[0]->getValue(provRowIdx);
                } else {
                    provnodes[i * nnodes + j] = provenance[(j - 1)*2 + 1]->
                        getValue(provRowIdx);
                    if (j > 1) {
                        auto tmp = provenance[(j - 1) * 2]->getValue(provRowIdx);
                        if (tmp == ~0ul) {
                            provRowIdx = 0;
                        } else {
                            provRowIdx = tmp;
                        }
                    }
                }
            }
            i++;
        }
        assert(nrows == i);

        //For each tuple, now I know the sequence of nodes that derived them.
        //I re-sort the nodes depending on the sequence of nodes
        std::vector<size_t> providxs;
        auto resortedSeg = seg->sortByProv(nnodes, providxs, provnodes);
        size_t startidx = 0;
        std::vector<size_t> currentNodeList(nnodes);
        for(size_t i = 0; i < nrows; ++i) {
            bool hasChanged = i == 0;
            for(size_t j = 0; j < nnodes && !hasChanged; ++j) {
                const size_t m = providxs[i] * nnodes + j;
                if (currentNodeList[j] != provnodes[m]) {
                    hasChanged = true;
                }
            }
            if (hasChanged) {
                if (startidx < i) {
                    //Create a new node
                    auto nodeId = getNNodes();
                    auto dataToAdd = resortedSeg->slice(nodeId, startidx, i);
                    if (filterFromProvenance) {
                        auto retainedData = retainFromDerivationTree(predId,
                                dataToAdd, currentNodeList);
                        if (retainedData->getNRows() > 0) {
                            addNodeProv(predId, ruleIdx, step, retainedData,
                                    currentNodeList);
                        }
                    } else {
                        addNodeProv(predId, ruleIdx, step, dataToAdd,
                                currentNodeList);
                    }
                }
                startidx = i;
                for(size_t j = 0; j < nnodes; ++j) {
                    const size_t m = providxs[i] * nnodes + j;
                    currentNodeList[j] = provnodes[m];
                }
            }
        }
        //Copy the last segment
        if (startidx < nrows) {
            auto nodeId = getNNodes();
            auto dataToAdd = resortedSeg->slice(nodeId, startidx, nrows);
            if (filterFromProvenance) {
                auto retainedData = retainFromDerivationTree(predId, dataToAdd,
                        currentNodeList);
                if (retainedData->getNRows() > 0) {
                    addNodeProv(predId, ruleIdx, step, retainedData,
                            currentNodeList);
                }
            } else {
                addNodeProv(predId, ruleIdx,
                        step, dataToAdd, currentNodeList);
            }
        }
    }
}

void GBGraph::addNodeProv(PredId_t predid,
        size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> data,
        const std::vector<size_t> &incomingEdges) {

#ifdef DEBUG
    if (provenanceType == ProvenanceType::NOPROV) {
        LOG(ERRORL) << "This method should be called only if type provenance is"
            " not set to NOPROV";
        throw 10;
    }
    if (provenanceType == FULLPROV && data->getProvenanceType() != SEG_FULLPROV) {
        LOG(ERRORL) << "The node does not have a good provenance type";
        throw 10;
    }
#endif

    //Check that the segment is sorted
#ifdef DEBUG
    auto ncols = data->getNColumns();
    size_t noffcols = 0;
    if (data->getNOffsetColumns() > 1)
        noffcols = data->getNOffsetColumns() - 1;
    auto itr = data->iterator();
    assert(itr->hasNext());
    itr->next();
    auto values = std::unique_ptr<Term_t[]>(new Term_t[ncols]);
    std::unique_ptr<Term_t[]> valuesOff;
    if (noffcols > 0)
        valuesOff = std::unique_ptr<Term_t[]>(new Term_t[noffcols]);

    for(int i = 0; i < ncols; ++i) {
        values[i] = itr->get(i);
    }
    for(size_t i = 0; i < noffcols; ++i) {
        valuesOff[i] = itr->getProvenanceOffset(i);
    }
    while (itr->hasNext()) {
        itr->next();
        bool equal = true;
        for(int i = 0; i < ncols; ++i) {
            if (itr->get(i) > values[i]) {
                equal = false;
                break;
            } else if (itr->get(i) < values[i]) {
                LOG(ERRORL) << "Segment not sorted";
                throw 10;
            }
        }
        if (equal) {
            if (!duplAllowed) {
                LOG(ERRORL) << "duplicate values are not allowed";
                throw 10;
            } else {
                //The offset should be different
                equal = noffcols > 0;
                for(size_t i = 0; i < noffcols; ++i) {
                    if (itr->getProvenanceOffset(i) != valuesOff[i]) {
                        equal = false;
                        break;
                    }
                }
                if (equal) {
                    LOG(ERRORL) << "The offset should be different";
                    throw 10;
                }
            }
        }
        for(int i = 0; i < ncols; ++i) {
            values[i] = itr->get(i);
        }
        for(size_t i = 0; i < noffcols; ++i) {
            valuesOff[i] = itr->getProvenanceOffset(i);
        }
    }
#endif

#ifdef DEBUG
    //ruleIdx can be ~0ul if the node is manually added (this happens, e.g., during MS)
    if (ruleIdx != ~0ul) {
        auto rule = program->getRule(ruleIdx);
        auto firstHead = rule.getHeads()[0];
        auto card = firstHead.getTupleSize();
        assert(data->getNColumns() == card);
        if (provenanceType == FULLPROV) {
            assert(data->getProvenanceType() == SEG_FULLPROV);
            auto n = rule.getBody().size();
            assert(data->getNOffsetColumns() == 1 + n);

            size_t idxIncomingEdge = 0;
            for(int idxBodyAtom = 0; idxBodyAtom < rule.getBody().size();
                    idxBodyAtom++) {
                auto b = rule.getBody()[idxBodyAtom];
                if (b.getPredicate().getType() != EDB) {
                    assert(idxIncomingEdge < incomingEdges.size());
                    auto node = incomingEdges[idxIncomingEdge];
                    auto card = getNodeSize(node);
                    auto itr = data->iterator();
                    while (itr->hasNext()) {
                        itr->next();
                        auto off = itr->getProvenanceOffset(idxIncomingEdge);
                        assert(off < card);
                    }
                    idxIncomingEdge += 1;
                }
            }
        }
    }
#endif

    auto nodeId = getNNodes();
    nodes.emplace_back();
    GBGraph_Node &outputNode = nodes.back();
    outputNode.predid = predid;
    outputNode.ruleIdx = ruleIdx;
    outputNode.step = step;
    outputNode.setData(data);
    assert(nodeId == data->getNodeId());

#ifdef DEBUG
    if (data->getName() == "CompositeTGSegment") {
        LOG(ERRORL) << "CompositeTGSegment is meant to be used only during the"
            " execution of a single rule. It should not be added to the graph";
        throw 10;
    }

    //The check below should never be triggered because the cliques
    //are never entirely novel. Thus, duplicate elimination will create a
    //new data structure.
    if (data->hasColumnarBackend()) {
        //Check if one of the columns come from a EDB layout that changes.
        //Issue a warning if that happens
        auto s = (TGSegmentLegacy*)data.get();
        for(size_t i = 0; i < data->getNColumns(); ++i) {
            auto c = s->getColumn(i);
            if (c->isEDB()) {
                auto ec = (EDBColumn*)c.get();
                const auto &layer = ec->getEDBLayer();
                const auto &query = ec->getLiteral();
                const auto table = layer.getEDBTable(query.getPredicate().getId());
                if (table->canChange()) {
                    throw 10;
                }
            }
        }
    }
#endif

    for(auto n : incomingEdges)
        if (isTmpNode(n) && n != ~0ul) {
            throw 10;
        }
    outputNode.setIncomingEdges(incomingEdges);
    pred2Nodes[predid].push_back(nodeId);
    LOG(DEBUGL) << "Added node ID " << nodeId << " with # facts=" <<
        data->getNRows();
}

uint64_t GBGraph::addTmpNode(PredId_t predId,
        std::shared_ptr<const TGSegment> data) {
    mapTmpNodes.insert(std::make_pair(
                counterTmpNodes, GBGraph_Node()));
    GBGraph_Node &n = mapTmpNodes[counterTmpNodes];
    n.predid = predId;
    n.ruleIdx = ~0ul;
    n.step = ~0ul;
    n.setData(data);
    return counterTmpNodes++;
}

void GBGraph::addNodeToBeRetained(PredId_t predId,
        std::shared_ptr<const TGSegment> data,
        std::vector<std::shared_ptr<Column>> &nodes,
        size_t ruleIdx,
        size_t step) {

#ifdef DEBUG
    auto ncols = data->getNColumns();
    auto itr = data->iterator();
    assert(itr->hasNext());
    itr->next();
    auto values = std::unique_ptr<Term_t[]>(new Term_t[ncols]);
    for(int i = 0; i < ncols; ++i) {
        values[i] = itr->get(i);
    }
    while (itr->hasNext()) {
        itr->next();
        for(int i = 0; i < ncols; ++i) {
            if (itr->get(i) > values[i]) {
                break;
            } else if (itr->get(i) < values[i]) {
                LOG(ERRORL) << "Segment not sorted";
                throw 10;
            }
        }
        for(int i = 0; i < ncols; ++i) {
            values[i] = itr->get(i);
        }
    }
#endif

    if (!mapPredTmpNodes.count(predId)) {
        mapPredTmpNodes.insert(std::make_pair(predId,
                    std::vector<GBGraph_TmpPredNode>()));
    }
    GBGraph_TmpPredNode n;
    n.data = data;
    n.nodes = nodes;
    n.ruleIdx = ruleIdx;
    n.step = step;
    mapPredTmpNodes[predId].push_back(n);
}

void GBGraph::replaceEqualTerms(
        size_t ruleIdx,
        size_t step,
        std::shared_ptr<const TGSegment> data) {
    std::vector<std::pair<Term_t, Term_t>> termsToReplace;
    //Fill termsToReplace
    auto itr = data->iterator();
    assert(data->getNColumns() == 2);
    while (itr->hasNext()) {
        itr->next();
        auto v1 = itr->get(0);
        auto v2 = itr->get(1);
        if (v1 == v2)
            continue;
        if (v1 < v2)
            termsToReplace.push_back(std::make_pair(v1, v2));
        else
            termsToReplace.push_back(std::make_pair(v2, v1));
    }
    if (termsToReplace.empty())
        return;
    std::sort(termsToReplace.begin(), termsToReplace.end());
    auto it = std::unique (termsToReplace.begin(), termsToReplace.end());
    termsToReplace.resize(std::distance(termsToReplace.begin(),it));

    //Create a map with all the elements to substitute
    FinalEGDTermMap map;
    map.set_empty_key((Term_t) -1);
    for(auto &pair : termsToReplace) {
        uint64_t key = pair.first;
        uint64_t value = pair.second;
        if (key == value)
            continue;
        if (((key & RULEVARMASK) == 0) && ((value & RULEVARMASK) == 0)) {
            LOG(ERRORL) << "Due to UNA, the chase does not exist (" <<
                key << "," << value << ")";
            throw 10;
        }
        assert(key < value);
        if (!map.count(value)) {
            map.insert(std::make_pair(value, key));
        } else {
            auto prevKey = map[value];
            bool prevKeySmallerThanKey = prevKey < key;
            if (!prevKeySmallerThanKey) {
                map[value] = key;
            }
        }
    }
    while (true) {
        bool replacedEntry = false;
        for(auto &pair : map) {
            if (map.count(pair.second)) {
                //Replace the value with the current one
                auto &v = map[pair.second];
                pair.second = v;
                replacedEntry = true;
            }
        }
        if (!replacedEntry) {
            break;
        }
    }

    //Consider all the nodes one-by-one and do the replacement
    assert(map.size() > 0);
    for(auto pair : pred2Nodes) {
        PredId_t predid = pair.first;
        auto &nodeIDs = pair.second;
        assert(!nodes.empty());
        const auto card = getNodeData(nodeIDs[0])->getNColumns();
        const int extraColumns = shouldTrackProvenance() ? 1 : 0;
        const auto nfields = card + extraColumns;
        std::unique_ptr<GBSegmentInserter> rewrittenTuples =
            GBSegmentInserter::getInserter(nfields, extraColumns, false);
        std::unique_ptr<Term_t[]> row = std::unique_ptr<Term_t[]>(
                new Term_t[nfields]);
        for(auto &nodeId : nodeIDs) {
            auto data = getNodeData(nodeId);
            size_t countAffectedTuples = 0;
            size_t countUnaffectedTuples = 0;
            std::unique_ptr<GBSegmentInserter> oldTuples;
            auto itr = data->iterator();
            while (itr->hasNext()) {
                itr->next();
                bool found = false;
                for(int i = 0; i < card; ++i) {
                    auto v = itr->get(i);
                    if (map.count(v)) {
                        found = true;
                        row[i] = map[v];
                    } else {
                        row[i] = v;
                    }
                }
                if (found) {
                    if (shouldTrackProvenance()) {
                        row[card] = ~0ul;
                    }
                    rewrittenTuples->add(row.get());
                    countAffectedTuples++;
                    if (countUnaffectedTuples > 0 &&
                            oldTuples.get() == NULL) {
                        //Copy all previous tuples
                        oldTuples = GBSegmentInserter::getInserter(nfields,
                                extraColumns, false);
                        auto itr2 = data->iterator();
                        for(size_t i = 0; i < countUnaffectedTuples; ++i) {
                            itr2->next();
                            for(int i = 0; i < card; ++i) {
                                row[i] = itr2->get(i);
                            }
                            if (shouldTrackProvenance()) {
                                row[card] = itr2->getNodeId();
                            }
                            oldTuples->add(row.get());
                        }
                    }
                } else {
                    if (shouldTrackProvenance()) {
                        row[card] = itr->getNodeId();
                    }
                    if (countAffectedTuples > 0 &&
                            oldTuples.get() == NULL) {
                        oldTuples = GBSegmentInserter::getInserter(nfields,
                                extraColumns, false);
                    }
                    if (oldTuples.get() != NULL) {
                        oldTuples->add(row.get());
                    }
                    countUnaffectedTuples++;
                }
            }
            if (oldTuples.get() != NULL) {
                auto tuples = oldTuples->getSegment(nodes[nodeId].step,
                        true, 0, getSegProvenanceType());
                nodes[nodeId].setData(tuples);
            } else {
                if (rewrittenTuples->getNRows() == data->getNRows()) {
                    oldTuples = GBSegmentInserter::getInserter(nfields,
                            extraColumns, false);
                    auto tuples = oldTuples->getSegment(nodes[nodeId].step,
                            true, 0, getSegProvenanceType());
                    nodes[nodeId].setData(tuples);
                } else {
                    assert(countUnaffectedTuples == data->getNRows());
                }
            }
        }
        //Create a new node with the replaced tuples
        if (!rewrittenTuples->isEmpty()) {
            //Invalidate the cache?
            if (cacheRetainEnabled && cacheRetain.count(predid)) {
                cacheRetain.erase(cacheRetain.find(predid));
            }
            auto tuples = rewrittenTuples->getSegment(~0ul,
                    false,
                    0,
                    getSegProvenanceType());
            tuples = tuples->sort()->unique();
            //Retain
            auto retainedTuples = retain(predid, tuples);
            bool nonempty = !(retainedTuples == NULL ||
                    retainedTuples->isEmpty());
            if (nonempty) {
                //Add new nodes
                if (shouldTrackProvenance()) {
                    auto nodeId = getNNodes();
                    auto dataToAdd = tuples->slice(nodeId, 0,
                            tuples->getNRows());
                    std::vector<size_t> nodes;
                    //Merge nodes lose the provenance
                    addNodeProv(predid, ruleIdx, step, dataToAdd,
                            nodes);
                } else {
                    //Add a single node
                    addNodeNoProv(predid, ruleIdx, step, retainedTuples);
                }
            }
        }
    }
}

size_t GBGraph::getNEdges() const {
    size_t out = 0;
    for (const auto &n : nodes) {
        out += n.getIncomingEdges(false).size();
    }
    return out;
}

size_t GBGraph::getNFacts() const {
    size_t out = 0;
    for (const auto &n : nodes) {
        out += n.getData()->getNRows();
    }
    return out;
}

std::shared_ptr<GBQuerier> GBGraph::getQuerier() const {
    return std::shared_ptr<GBQuerier>(new GBQuerier(*this, *program, *layer));
}

SegProvenanceType GBGraph::getSegProvenanceType(bool multipleNodes) const {
    if (provenanceType == ProvenanceType::NOPROV) {
        return SegProvenanceType::SEG_NOPROV;
    } else if (provenanceType == ProvenanceType::NODEPROV) {
        if (multipleNodes) {
            return SegProvenanceType::SEG_DIFFNODES;
        } else {
            return SegProvenanceType::SEG_SAMENODE;
        }
    } else if (provenanceType == ProvenanceType::FULLPROV) {
        return SegProvenanceType::SEG_FULLPROV;
    } else {
        throw 10;
    }
}
