#include <glog/gbgraph.h>
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
        const std::vector<std::shared_ptr<Column>> &provenance) {

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
                    addNodeProv(predId, ruleIdx, step, c,
                            currentNodeList);
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
            addNodeProv(predId, ruleIdx, step, dataToAdd, provnodes);
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
                    addNodeProv(predId, ruleIdx, step, dataToAdd,
                            currentNodeList);
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
            addNodeProv(predId, ruleIdx,
                    step, dataToAdd, currentNodeList);
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

#ifdef DEBUG
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
        if (isTmpNode(n))
            throw 10;
    outputNode.incomingEdges = incomingEdges;
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

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<std::shared_ptr<Column>> &derivationNodes) {
    //Special case for unary relations
    if (provenanceType != FULLPROV && existuples->getNColumns() == 1) {
        return retainVsNodeFast_one(existuples,
                newtuples, derivationNodes);
    } else if (provenanceType != FULLPROV && existuples->getNColumns() == 2) {
        return retainVsNodeFast_two(existuples, newtuples, derivationNodes);
    } else {
        return retainVsNodeFast_generic(existuples, newtuples, derivationNodes);
    }
}

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast_one(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<std::shared_ptr<Column>> &derivationNodes) {
    if (newtuples->hasColumnarBackend()) {
        ColumnWriter writer;
        bool allNew = true;
        auto newColumn = ((TGSegmentLegacy*)newtuples.get())->getColumn(0);
        if (existuples->hasColumnarBackend() && derivationNodes.size() == 0) {
            auto existColumn = ((TGSegmentLegacy*)existuples.get())
                ->getColumn(0);
            if (newColumn->isEDB() && existColumn->isEDB()) {
                //Get literal and pos join
                EDBLayer &layer = ((EDBColumn*)newColumn.get())->getEDBLayer();
                const Literal &l1 = ((EDBColumn*)newColumn.get())->getLiteral();
                uint8_t pos1 = ((EDBColumn*)newColumn.get())
                    ->posColumnInLiteral();
                std::vector<uint8_t> posColumns1;
                posColumns1.push_back(pos1);

                const Literal l2 = ((EDBColumn*)existColumn.get())->getLiteral();
                uint8_t pos2 = ((EDBColumn*)existColumn.get())
                    ->posColumnInLiteral();
                std::vector<uint8_t> posColumns2;
                posColumns2.push_back(pos2);

                std::vector<std::shared_ptr<Column>> columns = layer.
                    checkNewIn(l1, posColumns1, l2,
                            posColumns2);

                auto retainedColumn = columns[0];
                if (retainedColumn->isEmpty()) {
                    return std::shared_ptr<const TGSegment>();
                } else {
                    assert(retainedColumn->isBackedByVector());
                    std::vector<Term_t> tuples;
                    ((InmemoryColumn*)retainedColumn.get())->swap(tuples);
                    if (shouldTrackProvenance()) {
                        return std::shared_ptr<const TGSegment>(
                                new UnaryWithConstProvTGSegment(tuples,
                                    newtuples->getNodeId(),
                                    true, 0));
                    } else {
                        return std::shared_ptr<const TGSegment>(
                                new UnaryTGSegment(tuples,
                                    newtuples->getNodeId(),
                                    true, 0));
                    }
                }
            } else {
                bool allNew = Column::antijoin(newColumn, existColumn, writer);
                if (allNew) {
                    return newtuples;
                } else {
                    if (writer.isEmpty()) {
                        return std::shared_ptr<const TGSegment>();
                    } else {
                        std::vector<Term_t> tuples;
                        tuples.swap(writer.getValues());
                        if (shouldTrackProvenance()) {
                            return std::shared_ptr<const TGSegment>(
                                    new UnaryWithConstProvTGSegment(tuples,
                                        newtuples->getNodeId(),
                                        true, 0));
                        } else {
                            return std::shared_ptr<const TGSegment>(
                                    new UnaryTGSegment(tuples,
                                        newtuples->getNodeId(),
                                        true, 0));
                        }
                    }
                }
            }
        } else {
            //TODO antijoin between a column and a TGUnarySegment
            return retainVsNodeFast_generic(existuples, newtuples, derivationNodes);
        }
    } else {
        return retainVsNodeFast_generic(existuples, newtuples, derivationNodes);
    }
}

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast_two(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<std::shared_ptr<Column>> &derivationNodes) {


    if (newtuples->hasColumnarBackend() && derivationNodes.size() == 0) {
        auto newColumn1 = ((TGSegmentLegacy*)newtuples.get())->getColumn(0);
        auto newColumn2 = ((TGSegmentLegacy*)newtuples.get())->getColumn(1);
        if (newColumn1->isEDB() && newColumn2->isEDB()) {
            ColumnWriter writer;
            bool allNew = true;
            //Get literal and pos join
            EDBLayer &layer = ((EDBColumn*)newColumn1.get())->getEDBLayer();
            const Literal &l1 = ((EDBColumn*)newColumn1.get())->getLiteral();
            uint8_t pos1 = ((EDBColumn*)newColumn1.get())
                ->posColumnInLiteral();
            const Literal &l2 = ((EDBColumn*)newColumn2.get())->getLiteral();
            uint8_t pos2 = ((EDBColumn*)newColumn2.get())
                ->posColumnInLiteral();
            std::vector<Substitution> subs;
            if (!l1.sameVarSequenceAs(l2) ||
                    l1.subsumes(subs, l1, l2) == -1) {
                //The columns come from different literals. This is not yet
                //supported
                throw 10;
            }
            std::vector<uint8_t> posColumnsNew;
            posColumnsNew.push_back(pos1);
            posColumnsNew.push_back(pos2);

            if (existuples->hasColumnarBackend()) {
                auto existColumn1 = ((TGSegmentLegacy*)existuples.get())
                    ->getColumn(0);
                auto existColumn2 = ((TGSegmentLegacy*)existuples.get())
                    ->getColumn(1);

                if (existColumn1->isEDB() && existColumn2->isEDB()) {
                    const Literal l3 = ((EDBColumn*)existColumn1.get())->getLiteral();
                    uint8_t pos3 = ((EDBColumn*)existColumn1.get())
                        ->posColumnInLiteral();
                    const Literal l4 = ((EDBColumn*)existColumn2.get())->getLiteral();
                    uint8_t pos4 = ((EDBColumn*)existColumn2.get())
                        ->posColumnInLiteral();
                    if (!l3.sameVarSequenceAs(l4) ||
                            l1.subsumes(subs, l3, l4) == -1) {
                        //The columns come from different literals. This is not yet
                        //supported
                        throw 10;
                    }
                    std::vector<uint8_t> posColumnsExs;
                    posColumnsExs.push_back(pos3);
                    posColumnsExs.push_back(pos4);

                    std::vector<std::shared_ptr<Column>> retainedColumns = layer.
                        checkNewIn(l1, posColumnsNew, l3, posColumnsExs);
                    assert(retainedColumns.size() == 2);

                    auto retainedColumn1 = retainedColumns[0];
                    auto retainedColumn2 = retainedColumns[1];
                    if (retainedColumn1->isEmpty()) {
                        return std::shared_ptr<const TGSegment>();
                    } else {
                        assert(retainedColumn1->isBackedByVector());
                        assert(retainedColumn2->isBackedByVector());
                        std::vector<std::pair<Term_t, Term_t>> tuples;

                        const auto &v1 = retainedColumn1->getVectorRef();
                        const auto &v2 = retainedColumn2->getVectorRef();
                        for(size_t i = 0; i < v1.size(); ++i) {
                            tuples.push_back(std::make_pair(v1[i], v2[i]));
                        }

                        if (shouldTrackProvenance()) {
                            return std::shared_ptr<const TGSegment>(
                                    new BinaryWithConstProvTGSegment(tuples,
                                        newtuples->getNodeId(),
                                        true, 0));
                        } else {
                            return std::shared_ptr<const TGSegment>(
                                    new BinaryTGSegment(tuples,
                                        newtuples->getNodeId(),
                                        true, 0));
                        }
                    }
                }
            } else {
                assert(existuples->getNColumns() == 2);
                assert(existuples->getProvenanceType() != SEG_DIFFNODES);
                const std::vector<std::pair<Term_t, Term_t>> &t =
                    existuples->getProvenanceType() == SEG_NOPROV ?
                    ((BinaryTGSegment*)existuples.get())->getTuples() :
                    ((BinaryWithConstProvTGSegment*)existuples.get())->getTuples();

                std::vector<std::pair<Term_t, Term_t>> retained = layer.
                    checkNewIn(l1, posColumnsNew, t);
                if (retained.empty()) {
                    return std::shared_ptr<const TGSegment>();
                } else {
                    if (shouldTrackProvenance()) {
                        return std::shared_ptr<const TGSegment>(
                                new BinaryWithConstProvTGSegment(retained,
                                    newtuples->getNodeId(),
                                    true, 0));
                    } else {
                        return std::shared_ptr<const TGSegment>(
                                new BinaryTGSegment(retained,
                                    newtuples->getNodeId(),
                                    true, 0));
                    }
                }
            }
        }
    }
    return retainVsNodeFast_generic(existuples, newtuples, derivationNodes);
}

/*void GBGraph::filterOutDerivationNodes(std::vector<size_t> &idsToFilter,
  std::vector<std::shared_ptr<Column>> &derivationNodes) {
//I filter only the last two columns
assert(derivationNodes.size() >= 2);
std::vector<ColumnWriter> writers(2);
std::vector<std::unique_ptr<ColumnReader>> readers;
for(size_t i = derivationNodes.size() - 2; i < derivationNodes.size(); ++i) {
readers.push_back(derivationNodes[i]->getReader());
}
size_t nrows = derivationNodes.back()->size();
size_t currentIdx = 0;
for(size_t i = 0; i < nrows; ++i) {
for(size_t j = 0; j < readers.size(); ++j) {
assert(readers[j]->hasNext());
}
if(i != idsToFilter[currentIdx]) {
for(size_t j = 0; j < readers.size(); ++j) {
auto v = readers[j]->next();
writers[j].add(v);
}
} else {
for(size_t j = 0; j < readers.size(); ++j) {
readers[j]->next();
}
currentIdx++;
}
}
derivationNodes[derivationNodes.size() - 2] = writers[0].getColumn();
derivationNodes[derivationNodes.size() - 1] = writers[1].getColumn();
}*/

/*void GBGraph::shuffleDerivationNodes(std::vector<size_t> &idsToShuffle,
  std::vector<std::shared_ptr<Column>> &derivationNodes) {
//I filter only the last two columns
assert(derivationNodes.size() >= 2);
auto c1 = derivationNodes[derivationNodes.size() - 2];
auto c2 = derivationNodes[derivationNodes.size() - 1];

std::vector<ColumnWriter> writers(2);
for(size_t i = 0; i < idsToShuffle.size(); ++i) {
auto idx = idsToShuffle[i];
auto v1 = c1->getValue(idx);
writers[0].add(v1);
auto v2 = c2->getValue(idx);
writers[1].add(v2);
}
derivationNodes[derivationNodes.size() - 2] = writers[0].getColumn();
derivationNodes[derivationNodes.size() - 1] = writers[1].getColumn();
}*/

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast_generic(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<std::shared_ptr<Column>> &derivationNodes) {

    std::unique_ptr<GBSegmentInserter> inserter;
    const uint8_t ncols = newtuples->getNColumns();
    const size_t extracols = newtuples->getNOffsetColumns();
    Term_t row[ncols + extracols];


    //Do outer join
    auto leftItr = existuples->iterator();
    auto rightItr = newtuples->iterator();
    bool moveLeftItr = true;
    bool moveRightItr = true;
    bool activeRightValue = false;
    size_t countNew = 0;
    bool isFiltered = false;
    size_t startCopyingIdx = 0;
    int64_t countRight = -1;
    while (true) {
        if (moveRightItr) {
            if (rightItr->hasNext()) {
                rightItr->next();
                activeRightValue = true;
                countRight++;
            } else {
                activeRightValue = false;
                break;
            }
            moveRightItr = false;
        }

        if (moveLeftItr) {
            if (leftItr->hasNext()) {
                leftItr->next();
            } else {
                break;
            }
            moveLeftItr = false;
        }

        //Compare the iterators
        int res = TGSegmentItr::cmp(leftItr.get(), rightItr.get());
        if (res < 0) {
            moveLeftItr = true;
        } else if (res > 0) {
            moveRightItr = true;
            if (isFiltered) {
                for(int i = 0; i < ncols; ++i) {
                    row[i] = rightItr->get(i);
                }
                if (extracols > 0) {
                    row[ncols] = rightItr->getNodeId();
                    assert(row[ncols] != ~0ul);
                    for(size_t i = 1; i < extracols; ++i) {
                        row[ncols + i] = rightItr->getProvenanceOffset(i-1);
                    }
                }
                inserter->add(row);
            } else {
                countNew++;
            }
        } else {
            moveLeftItr = moveRightItr = true;
            activeRightValue = false;
            if (!isFiltered && countNew == 0)
                startCopyingIdx++;

            //The tuple must be filtered out
            if (!isFiltered && countNew > 0) {
                inserter = GBSegmentInserter::getInserter(ncols + extracols,
                        extracols, false);
                //Copy all the previous new tuples in the right iterator
                size_t i = 0;
                auto itrTmp = newtuples->iterator();
                while (i < (startCopyingIdx + countNew) && itrTmp->hasNext()) {
                    itrTmp->next();
                    if (i >= startCopyingIdx) {
                        for(int i = 0; i < ncols; ++i) {
                            row[i] = itrTmp->get(i);
                        }
                        if (extracols > 0) {
                            row[ncols] = itrTmp->getNodeId();
                            assert(row[ncols] != ~0ul);
                            for(size_t i = 1; i < extracols; ++i) {
                                row[ncols + i] = itrTmp->getProvenanceOffset(i-1);
                            }
                        }
                        inserter->add(row);
                    }
                    i++;
                }
                isFiltered = true;
            }
        }
    }

    if (isFiltered) {
        if (activeRightValue) {
            for(int i = 0; i < ncols; ++i) {
                row[i] = rightItr->get(i);
            }
            if (extracols > 0) {
                row[ncols] = rightItr->getNodeId();
                for(size_t i = 1; i < extracols; ++i) {
                    row[ncols + i] = rightItr->getProvenanceOffset(i-1);
                }
            }
            inserter->add(row);
        }
        while (rightItr->hasNext()) {
            rightItr->next();
            for(int i = 0; i < ncols; ++i) {
                row[i] = rightItr->get(i);
            }
            if (extracols > 0) {
                row[ncols] = rightItr->getNodeId();
                for(size_t i = 1; i < extracols; ++i) {
                    row[ncols + i] = rightItr->getProvenanceOffset(i-1);
                }
            }
            inserter->add(row);
        }

        return inserter->getSegment(newtuples->getNodeId(), true, 0,
                getSegProvenanceType(), extracols);
    } else {
        if (countNew > 0 || activeRightValue) {
            if (startCopyingIdx == 0) {
                //They are all new ...
                return newtuples;
            } else {
                //Remove the initial duplicates
                inserter = GBSegmentInserter::getInserter(ncols + extracols,
                        extracols, false);
                //Copy all the previous new tuples in the right iterator
                size_t i = 0;
                auto itrTmp = newtuples->iterator();
                while (itrTmp->hasNext()) {
                    itrTmp->next();
                    if (i >= startCopyingIdx) {
                        for(int i = 0; i < ncols; ++i) {
                            row[i] = itrTmp->get(i);
                        }
                        if (extracols > 0) {
                            row[ncols] = itrTmp->getNodeId();
                            for(size_t i = 1; i < extracols; ++i) {
                                row[ncols + i] = itrTmp->getProvenanceOffset(i-1);
                            }
                        }
                        inserter->add(row);
                    }
                    i++;
                }

                return inserter->getSegment(newtuples->getNodeId(),
                        true, 0, getSegProvenanceType(), extracols);
            }
        } else {
            //They are all duplicates. Leave the derivationNodes unchanged, anyway
            //the node won't be added
            return std::shared_ptr<const TGSegment>();
        }
    }
}

std::shared_ptr<const TGSegment> GBGraph::retain(
        PredId_t p,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<std::shared_ptr<Column>> &derivationNodes) {
    //Check that it is sorted
#ifdef DEBUG
    auto ncols = newtuples->getNColumns();
    auto itr = newtuples->iterator();
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

    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    if (!areNodesWithPredicate(p)) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto dur = end - start;
        durationRetain += dur;
        return newtuples;
    }
    auto &nodeIdxs = getNodeIDsWithPredicate(p);
    if (cacheRetainEnabled &&
            nodeIdxs.size() > 1) {
        //Merge and sort the segments
        if (!cacheRetain.count(p) || cacheRetain[p].nnodes < nodeIdxs.size()) {
            //Get the arity
            if (newtuples->getNColumns() == 1) {
                std::vector<Term_t> tuples;
                size_t startNodeIdx = 0;
                if (cacheRetain.count(p)) {
                    cacheRetain[p].seg->appendTo(0, tuples);
                    startNodeIdx = cacheRetain[p].nnodes;
                }
                for(size_t i = startNodeIdx; i < nodeIdxs.size(); ++i) {
                    getNodeData(nodeIdxs[i])->appendTo(0, tuples);
                }
                std::sort(tuples.begin(), tuples.end());
                auto seg = std::shared_ptr<const TGSegment>(
                        new UnaryTGSegment(tuples, ~0ul, true, 0));
                CacheRetainEntry entry;
                entry.nnodes = nodeIdxs.size();
                entry.seg = seg;
                cacheRetain[p] = entry;
            } else if (newtuples->getNColumns() == 2) {
                //std::chrono::steady_clock::time_point startCache =
                //    std::chrono::steady_clock::now();

                std::vector<std::pair<Term_t, Term_t>> tuples;
                size_t startNodeIdx = 0;
                if (cacheRetain.count(p)) {
                    cacheRetain[p].seg->appendTo(0, 1, tuples);
                    startNodeIdx = cacheRetain[p].nnodes;
                }
                for(size_t i = startNodeIdx; i < nodeIdxs.size(); ++i) {
                    getNodeData(nodeIdxs[i])->appendTo(0, 1, tuples);
                }
                std::sort(tuples.begin(), tuples.end());
#ifdef DEBUG
                auto e = std::unique(tuples.begin(), tuples.end());
                auto d = std::distance(e, tuples.end());
                if (d > 0) {
                    LOG(ERRORL) << "Duplicates should not occur here!";
                    throw 10;
                }
#endif
                auto seg = std::shared_ptr<const TGSegment>(
                        new BinaryTGSegment(tuples, ~0ul, true, 0));
                CacheRetainEntry entry;
                entry.nnodes = nodeIdxs.size();
                entry.seg = seg;
                cacheRetain[p] = entry;

            } else if (newtuples->getNColumns() > 2) {
                size_t ncolumns = newtuples->getNColumns();
                std::vector<int> posFields;
                for(int i = 0; i < ncolumns; ++i)
                    posFields.push_back(i);
                std::vector<std::vector<Term_t>> tuples(ncolumns);
                size_t startNodeIdx = 0;
                if (cacheRetain.count(p)) {
                    cacheRetain[p].seg->appendTo(posFields, tuples);
                    startNodeIdx = cacheRetain[p].nnodes;
                }
                for(size_t i = startNodeIdx; i < nodeIdxs.size(); ++i) {
                    getNodeData(nodeIdxs[i])->appendTo(posFields, tuples);
                }
                size_t nrows = tuples[0].size();
                //Create columns from the content of tuples
                std::vector<std::shared_ptr<Column>> columns;
                for(int i = 0; i < ncolumns; ++i) {
                    columns.push_back(std::shared_ptr<Column>(
                                new InmemoryColumn(tuples[i], true)));
                }
                auto seg = std::shared_ptr<const TGSegment>(
                        new TGSegmentLegacy(columns, nrows));
                seg = seg->sort();
                CacheRetainEntry entry;
                entry.nnodes = nodeIdxs.size();
                entry.seg = seg;
                cacheRetain[p] = entry;
            } else {
                LOG(ERRORL) << "Retain with arity = 0 is not supported";
                throw 10;
            }
        }
        auto existingTuples = cacheRetain[p].seg;
        newtuples = retainVsNodeFast(existingTuples, newtuples,
                derivationNodes);
        if (newtuples == NULL || newtuples->isEmpty()) {
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            auto dur = end - start;
            durationRetain += dur;
            return std::shared_ptr<const TGSegment>();
        }
    } else {
        for(auto &nodeIdx : nodeIdxs) {
            auto nodeData = getNodeData(nodeIdx);
            newtuples = retainVsNodeFast(nodeData, newtuples,
                    derivationNodes);
            if (newtuples == NULL || newtuples->isEmpty()) {
                std::chrono::steady_clock::time_point end =
                    std::chrono::steady_clock::now();
                auto dur = end - start;
                durationRetain += dur;
                return std::shared_ptr<const TGSegment>();
            }
        }
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    auto dur = end - start;
    durationRetain += dur;
    return newtuples;
}

std::shared_ptr<const TGSegment> GBGraph::mergeNodes_special_unary1(
        std::shared_ptr<const TGSegment> seg,
        const std::vector<size_t> &nodeIdxs,
        const std::vector<int> &copyVarPos, bool lazyMode,
        bool replaceOffsets) const {
    assert(nodeIdxs.size() > 0);
    auto ncols = copyVarPos.size();
    bool project = ncols > 0 && ncols < getNodeData(nodeIdxs[0])->getNColumns();
    bool shouldSortAndUnique = (project ||
            (copyVarPos.size() > 0 && copyVarPos[0] != 0));

    std::vector<std::shared_ptr<Column>> projectedColumns;
    std::vector<std::shared_ptr<Column>> outColumns;
    seg->projectTo(copyVarPos, projectedColumns);
    assert(shouldTrackProvenance() || projectedColumns.size() == 1);
    outColumns.push_back(projectedColumns[0]);
    if (shouldSortAndUnique && provenanceType != FULLPROV) {
        if (outColumns.size() == 1) {
            auto sortedCol = outColumns[0]->sort();
            outColumns[0] = sortedCol->unique();
        } else {
            LOG(ERRORL) << "Should never happen";
            throw 10;
        }
    }
    size_t nrows = outColumns[0]->size();

    if (shouldTrackProvenance()) {
        assert(projectedColumns[1]->isConstant());
        bool isSorted = true;
        size_t nprovcolumns = seg->getNOffsetColumns();
        if (provenanceType != FULLPROV) {
            //Add one column with the provenance
            outColumns.push_back(std::shared_ptr<Column>(
                        new CompressedColumn(seg->getNodeId(),
                            nrows)));
        } else {
            if (replaceOffsets) {
                outColumns.push_back(std::shared_ptr<Column>(
                            new CompressedColumn(seg->getNodeId(),
                                nrows)));
                outColumns.push_back(std::shared_ptr<Column>(
                            new CompressedColumn(0,
                                nrows, 1)));
                nprovcolumns = 2;
            } else {
                for (size_t i = 0; i < nprovcolumns; ++i) {
                    outColumns.push_back(projectedColumns[1 + i]);
                }
            }
            isSorted = copyVarPos[0] == 0;
        }

        return std::shared_ptr<const TGSegment>(
                new TGSegmentLegacy(
                    outColumns,
                    nrows,
                    isSorted,
                    0,
                    getSegProvenanceType(nodeIdxs.size() > 1),
                    nprovcolumns));
    } else {
        return std::shared_ptr<const TGSegment>(
                new TGSegmentLegacy(
                    outColumns,
                    nrows,
                    true,
                    0,
                    getSegProvenanceType(nodeIdxs.size() > 1)));
    }
}

std::shared_ptr<const TGSegment> GBGraph::mergeNodes_special_unary2(
        const std::vector<size_t> &nodeIdxs,
        const std::vector<int> &copyVarPos, bool lazyMode,
        bool replaceOffsets) const {
    auto ncols = copyVarPos.size();
    bool project = ncols > 0 && ncols < getNodeData(nodeIdxs[0])->getNColumns();
    bool shouldSortAndUnique = (project ||
            (copyVarPos.size() > 0 && copyVarPos[0] != 0));

    if (shouldTrackProvenance()) {
        if (provenanceType == FULLPROV) {
            size_t maxOffset = 0;
            if (replaceOffsets) {
                maxOffset = 2;
            } else {
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    auto n = getNodeData(idbBodyAtomIdx)->getNOffsetColumns();
                    if (n > maxOffset) {
                        maxOffset = n;
                    }
                }
            }
            assert(maxOffset > 0);
            if (maxOffset == 1) {
                throw 10;
            } else if (maxOffset == 2) {
                std::vector<UnWithFullProv> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    size_t start = tuples.size();
                    getNodeData(idbBodyAtomIdx)->appendTo(
                            copyVarPos[0], tuples);
                    if (replaceOffsets) {
                        for(size_t i = start; i < tuples.size(); ++i) {
                            tuples[i].prov = i - start;
                        }
                    }
                }

                return std::shared_ptr<const TGSegment>(
                        new UnaryWithFullProvTGSegment(tuples, ~0ul,
                            false, 0));
            } else {
                throw 10;
            }
        } else {
            assert(provenanceType != FULLPROV);
            std::vector<std::pair<Term_t,Term_t>> tuples;
            for(auto idbBodyAtomIdx : nodeIdxs) {
                getNodeData(idbBodyAtomIdx)->appendTo(copyVarPos[0], tuples);
            }
            if (shouldSortAndUnique) {
                std::sort(tuples.begin(), tuples.end());
                auto itr = std::unique(tuples.begin(), tuples.end());
                tuples.erase(itr, tuples.end());
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithProvTGSegment(tuples, ~0ul, true, 0));
            } else {
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithProvTGSegment(tuples, ~0ul, false, 0));
            }
        }
    } else {
        std::vector<Term_t> tuples;
        for(auto idbBodyAtomIdx : nodeIdxs)
            getNodeData(idbBodyAtomIdx)->appendTo(copyVarPos[0], tuples);
        if (shouldSortAndUnique) {
            std::sort(tuples.begin(), tuples.end());
            auto itr = std::unique(tuples.begin(), tuples.end());
            tuples.erase(itr, tuples.end());
            return std::shared_ptr<const TGSegment>(
                    new UnaryTGSegment(tuples, ~0ul, true, 0));
        } else {
            return std::shared_ptr<const TGSegment>(
                    new UnaryTGSegment(tuples, ~0ul, false, 0));
        }
    }
}

std::shared_ptr<const TGSegment> GBGraph::mergeNodes(
        const std::vector<size_t> &nodeIdxs,
        const std::vector<int> &copyVarPos, bool lazyMode,
        bool replaceOffsets) const {
    assert(nodeIdxs.size() > 0);
    auto ncols = copyVarPos.size();
    bool project = ncols > 0 && ncols < getNodeData(nodeIdxs[0])->getNColumns();
    bool shouldSortAndUnique = (project ||
            (copyVarPos.size() > 0 && copyVarPos[0] != 0));

    if (copyVarPos.size() == 1) {
        if (nodeIdxs.size() == 1) {
            if (!project && (provenanceType != GBGraph::ProvenanceType::FULLPROV
                        || !replaceOffsets)) {
                return getNodeData(nodeIdxs[0]);
            }
            auto seg = getNodeData(nodeIdxs[0]);
            if (seg->hasColumnarBackend()) {
                return mergeNodes_special_unary1(seg,
                        nodeIdxs,
                        copyVarPos,
                        lazyMode,
                        replaceOffsets);
            }
        }

        if (lazyMode) {
            return std::shared_ptr<const TGSegment>(
                    new CompositeTGSegment(*this, nodeIdxs, copyVarPos,
                        false, 0, getSegProvenanceType(nodeIdxs.size() > 1), false, false,
                        replaceOffsets));
        }

        return mergeNodes_special_unary2(
                nodeIdxs,
                copyVarPos,
                lazyMode,
                replaceOffsets);

    } else if (copyVarPos.size() == 2) {
        if (nodeIdxs.size() == 1 && !project &&
                copyVarPos[0] == 0 && copyVarPos[1] == 1 &&
                (provenanceType != GBGraph::ProvenanceType::FULLPROV || !replaceOffsets)) {
            return getNodeData(nodeIdxs[0]);
        } else {
            if (lazyMode) {
                return std::shared_ptr<const TGSegment>(
                        new CompositeTGSegment(*this, nodeIdxs, copyVarPos,
                            false, 0, getSegProvenanceType(), false, false,
                            replaceOffsets));
            }

            if (shouldTrackProvenance()) {
                if (provenanceType == FULLPROV) {
                    size_t maxOffset = 0;
                    if (replaceOffsets) {
                        maxOffset = 2;
                    } else {
                        for(auto idbBodyAtomIdx : nodeIdxs) {
                            auto n = getNodeData(idbBodyAtomIdx)->getNOffsetColumns();
                            if (n > maxOffset) {
                                maxOffset = n;
                            }
                        }
                    }
                    assert(maxOffset > 0);
                    if (maxOffset == 1) {
                        throw 10;
                    } else if (maxOffset == 2) {
                        std::vector<BinWithFullProv> tuples;
                        for(auto idbBodyAtomIdx : nodeIdxs) {
                            auto start = tuples.size();
                            getNodeData(idbBodyAtomIdx)->appendTo(
                                    copyVarPos[0],
                                    copyVarPos[1], tuples);
                            if (replaceOffsets) {
                                for(size_t i = start; i < tuples.size(); ++i) {
                                    tuples[i].prov = i - start;
                                }
                            }
                        }
                        return std::shared_ptr<const TGSegment>(
                                new BinaryWithFullProvTGSegment(tuples,
                                    nodeIdxs.size() == 1 ? nodeIdxs[0] : ~0ul,
                                    false, 0));
                    } else {
                        throw 10;
                    }
                } else {
                    std::vector<BinWithProv> tuples;
                    for(auto idbBodyAtomIdx : nodeIdxs) {
                        getNodeData(idbBodyAtomIdx)->appendTo(
                                copyVarPos[0], copyVarPos[1], tuples);
                    }
                    if (shouldSortAndUnique) {
                        std::sort(tuples.begin(), tuples.end());
                        auto itr = std::unique(tuples.begin(), tuples.end());
                        tuples.erase(itr, tuples.end());
                        return std::shared_ptr<const TGSegment>(
                                new BinaryWithProvTGSegment(tuples, ~0ul, true, 0));
                    } else {
                        return std::shared_ptr<const TGSegment>(
                                new BinaryWithProvTGSegment(tuples, ~0ul, false, 0));
                    }
                }
            } else {
                std::vector<std::pair<Term_t,Term_t>> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    getNodeData(idbBodyAtomIdx)->appendTo(
                            copyVarPos[0], copyVarPos[1], tuples);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                    return std::shared_ptr<const TGSegment>(
                            new BinaryTGSegment(tuples, ~0ul, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new BinaryTGSegment(tuples, ~0ul, false, 0));
                }
            }
        }
    } else {
        //Could be 0 or > 2
        assert(copyVarPos.size() == 0 || copyVarPos.size() > 2);
        size_t ncolumns = copyVarPos.size();
        size_t maxOffset = 0;
        if (shouldTrackProvenance()) {
            ncolumns += 1; //Node
            if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
                if (replaceOffsets) {
                    ncolumns += 1;
                } else {
                    if (replaceOffsets) {
                        maxOffset = 2;
                    } else {
                        for(auto idbBodyAtomIdx : nodeIdxs) {
                            auto n = getNodeData(idbBodyAtomIdx)->
                                getNOffsetColumns() - 1;
                            if (n > maxOffset) {
                                maxOffset = n;
                            }
                        }
                    }
                    ncolumns += maxOffset;
                }
            }
        }
        std::vector<std::vector<Term_t>> tuples(ncolumns);
        for(auto i : nodeIdxs) {
            auto data = getNodeData(i);
            if (shouldTrackProvenance()) {
                if (!replaceOffsets &&
                        provenanceType == GBGraph::ProvenanceType::FULLPROV) {
                    assert(data->getNOffsetColumns() == maxOffset + 1);
                    //The other case is not implemented
                }
                data->appendTo(copyVarPos, tuples, true);
            } else {
                data->appendTo(copyVarPos, tuples, false);
            }
        }
        size_t nrows = tuples[0].size();
        //Create columns from the content of tuples
        std::vector<std::shared_ptr<Column>> columns;
        for(int i = 0; i < copyVarPos.size(); ++i) {
            columns.push_back(std::shared_ptr<Column>(
                        new InmemoryColumn(tuples[i], true)));
        }
        std::shared_ptr<const TGSegment> seg;
        if (shouldTrackProvenance()) {
            //Add the column with the node IDs
            columns.push_back(std::shared_ptr<Column>(
                        new InmemoryColumn(tuples.back(), true)));
            if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
                throw 10;
                //If replaceOffset is true, then I should have another column with
                //the offsets, otherwise there should be more columns with the
                //existing offsets
            }
        }
        seg = std::shared_ptr<const TGSegment>(
                new TGSegmentLegacy(columns, nrows, false, 0,
                    getSegProvenanceType(nodeIdxs.size() > 1)));
        if (shouldSortAndUnique) {
            auto sortedSeg = seg->sort();
            return sortedSeg->unique();
        } else {
            return seg;
        }
    }
}

bool retainAndAddFromTmpNodes_uniq(
        const std::pair<Term_t, Term_t> &a,
        const std::pair<Term_t, Term_t> &b) {
    return a.first == b.first;
}

void GBGraph::retainAndAddFromTmpNodes(PredId_t predId) {
    if (!mapPredTmpNodes.count(predId))
        return;

    const auto &nodes = mapPredTmpNodes[predId];
    auto card = program->getPredicateCard(predId);
    assert(nodes.size() > 0);
    if (nodes.size() >  (1 << 24)) {
        //The last three bytes are used for a
        //counter
        throw 10;
    };
    
    const size_t n_columns = nodes[0].data->getNColumns();
    size_t n_offsetcolumns = 1;
    for (auto &node : nodes) {
        assert(node.data->getNColumns() == n_columns);
        if (node.data->getNOffsetColumns() > n_offsetcolumns) {
            n_offsetcolumns = node.data->getNOffsetColumns();
        }
        assert(node.data->getNRows() < 0xFFFFFFFFFF); //Otherwise, the offset no longer works
    }
    auto inserter = GBSegmentInserter::getInserter(n_columns + n_offsetcolumns, n_offsetcolumns, true);
    auto row = std::unique_ptr<Term_t[]>(new Term_t[n_columns + n_offsetcolumns]);
    size_t counter = 0;
    for (auto &node : nodes) {
        auto &d = node.data;
        auto itr = d->iterator();
        const auto max_offset_node = node.data->getNOffsetColumns();
        assert(max_offset_node <= n_offsetcolumns);
        assert(node.data->getNodeId() < 0xFFFFFFFFFFl || node.data->getNodeId() == ~0ul);        
        while (itr->hasNext()) {
            itr->next();
            for(int i = 0; i < n_columns; ++i) {
                row[i] = itr->get(i);
            }
            if (provenanceType != NOPROV) {
                if (itr->getNodeId() == ~0ul)
                    row[n_columns] = counter + 0xFFFFFFFFFFl;
                else {
                    row[n_columns] = counter + itr->getNodeId();
                }
                for(int i = 1; i < max_offset_node; ++i) {
                    row[n_columns + i] = itr->getProvenanceOffset(i-1);
                }
            } else {
                row[n_columns] = counter;
            }
            inserter->add(row.get());
        }
        counter += (size_t)1 << 40;
    }
    
    auto seg = inserter->getSegment(~0ul, false, 0, getSegProvenanceType(true), n_offsetcolumns);
    auto sortedSeg = seg->sort();
    auto uniqueSeg = sortedSeg->unique();
    auto retainedSeg = retain(predId, uniqueSeg);
    if (retainedSeg.get() == NULL)
        return;
    auto toBeAddedSeg = retainedSeg->sortByProv();
    auto itr = toBeAddedSeg->iterator();
    auto node = nodes.begin();
    size_t beginSegment = 0;
    size_t endSegment = (size_t) 1 << 40;
    size_t counterSlice = 0;
    size_t startSlice = 0;
    
    while (itr->hasNext()) {
        itr->next();
        auto counter = itr->getNodeId();
        if (counter >= endSegment) {
            assert(node != nodes.end());
            if (startSlice < counterSlice) {
                auto slicedNode = toBeAddedSeg->slice(startSlice, counterSlice);
                auto rewrittenNode = retainAndAddFromTmpNodes_rewriteNode(slicedNode, *node);
                addNodesProv(predId, node->ruleIdx, node->step, rewrittenNode, node->nodes);
             }
            startSlice = counterSlice;
            
            //Advance to a node that has a limit
            while (counter >= endSegment) {
                if (node >= nodes.end() - 1) {
                    throw 10; //This should not happen!
                }
                node++;
                beginSegment += (size_t)1 << 40;
                endSegment += (size_t)1 << 40;
            }
        }
        counterSlice++;
    }
    assert(node != nodes.end());
    auto slicedNode = toBeAddedSeg->slice(startSlice, counterSlice);
    auto rewrittenNode = retainAndAddFromTmpNodes_rewriteNode(slicedNode, *node);
    addNodesProv(predId, node->ruleIdx, node->step, rewrittenNode, node->nodes);
}

std::shared_ptr<const TGSegment> GBGraph::retainAndAddFromTmpNodes_rewriteNode(
                                        std::shared_ptr<const TGSegment> toBeRewrittenNode,
                                        const GBGraph_TmpPredNode &node) {
    auto n_columns = toBeRewrittenNode->getNColumns();
    auto n_offsetcolumns = node.data->getNOffsetColumns(); //It's important that we pick the number of offset columns from the old node
    auto inserter = GBSegmentInserter::getInserter(n_columns + n_offsetcolumns, n_offsetcolumns, true);
    auto row = std::unique_ptr<Term_t[]>(new Term_t[n_columns + n_offsetcolumns]);
    auto itr = toBeRewrittenNode->iterator();
    
    size_t currentNode = ~0ul;
    bool firstRow = true;
    bool constantNode = true;
    while (itr->hasNext()) {
        itr->next();
        for(int i = 0; i < n_columns; ++i) {
            row[i] = itr->get(i);
        }
        if (provenanceType != NOPROV) {
            auto nodeId = itr->getNodeId() & 0xFFFFFFFFFFl; //Removing the initial prefix
            if (nodeId == 0xFFFFFFFFFFl)
                nodeId = ~0ul;
            row[n_columns] = nodeId;
            for(int i = 1; i < n_offsetcolumns; ++i) {
                row[n_columns + i] = itr->getProvenanceOffset(i - 1);
            }
            if (firstRow) {
                currentNode = nodeId;
                firstRow = false;
            } else if (currentNode != nodeId) {
                constantNode = false;
            }
        }
        inserter->add(row.get());
    }
    if (!constantNode) {
        currentNode = ~0ul;
    } else {
        assert(node.data->getProvenanceType() == SEG_FULLPROV || node.data->getProvenanceType() == SEG_NOPROV);
    }
    return inserter->getSegment(currentNode, true, 0, node.data->getProvenanceType(), n_offsetcolumns);
}

/*void GBGraph::retainAndAddFromTmpNodes_add_unary(
        bool storeNode,
        std::vector<Term_t> &t1,
        std::vector<std::pair<Term_t, Term_t>> &t2,
        PredId_t predId,
        const GBGraph_TmpPredNode &node,
        size_t beginSegment) {
    if (!storeNode && !t1.empty()) {
        auto segToStore = std::shared_ptr<TGSegment>(
                new UnaryWithConstProvTGSegment(t1, node.data->getNodeId(),
                    true, 0));
        addNodesProv(predId, node.ruleIdx, node.step, segToStore, node.nodes);
        t1 = std::vector<Term_t>();
    }
    if (storeNode && !t2.empty()) {
        //Rewrite the counters so that they will be equivalent
        //to the ones in the original
        for(auto &t : t2) {
            assert(t.second >= beginSegment);
            t.second = t.second - beginSegment;
            if (t.second >= ((size_t)1 << 40))
                throw 10;
        }
        auto segToStore = std::shared_ptr<TGSegment>(
                new UnaryWithProvTGSegment(t2, node.data->getNodeId(), true, 0));
        addNodesProv(predId, node.ruleIdx, node.step, segToStore, node.nodes);
        t2 = std::vector<std::pair<Term_t, Term_t>>();
    }
}*/

uint64_t GBGraph::mergeNodesWithPredicateIntoOne(PredId_t predId) {
    if (!areNodesWithPredicate(predId)) {
        return 0;
    }
    std::vector<size_t> sortedNodeIDs = getNodeIDsWithPredicate(predId);
    //If there is only one node, then removing duplicates might not be
    //necessary
    if (sortedNodeIDs.size() == 1) {
        return getNodeSize(sortedNodeIDs[0]);
    }

    std::sort(sortedNodeIDs.begin(), sortedNodeIDs.end());
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    //Merge the tuples into a single segment
    auto data = getNodeData(sortedNodeIDs[0]);
    auto card = data->getNColumns();
    std::vector<int> copyVarPos;
    for(size_t i = 0; i < card; ++i)
        copyVarPos.push_back(i);
    auto tuples = mergeNodes(sortedNodeIDs, copyVarPos);
    std::chrono::duration<double> dur = std::chrono::system_clock::now() - start;
    LOG(DEBUGL) << "Merged " << sortedNodeIDs.size() << " in a container with"
        << "  " << tuples->getNRows() << " in " << dur.count() * 1000 << "ms";

    start = std::chrono::system_clock::now();
    tuples = tuples->sort();
    dur = std::chrono::system_clock::now() - start;
    LOG(DEBUGL) << "Sorted " << tuples->getNRows() << " tuples in " <<
        dur.count() * 1000 << "ms";

    start = std::chrono::system_clock::now();
    tuples = tuples->unique();
    dur = std::chrono::system_clock::now() - start;
    LOG(DEBUGL) << "Retained " << tuples->getNRows() << " unique tuples in " <<
        dur.count() * 1000 << "ms";

    //Remove the content of the pre-existing nodes
    size_t lastStep = 0;
    for(auto &nodeId : sortedNodeIDs) {
        nodes[nodeId].setData(tuples->slice(nodeId, 0, 0));
        assert(getNodeSize(nodeId) == 0);
        if (nodes[nodeId].step > lastStep)
            lastStep = nodes[nodeId].step;
    }

    //Create a new node
    if (shouldTrackProvenance()) {
        std::vector<size_t> nodes;
        addNodeProv(predId, ~0ul, lastStep, tuples, nodes);
    } else {
        addNodeNoProv(predId, ~0ul, lastStep, tuples);
    }
    return tuples->getNRows();
}

size_t GBGraph::getNEdges() const {
    size_t out = 0;
    for (const auto &n : nodes) {
        out += n.incomingEdges.size();
    }
    return out;
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
