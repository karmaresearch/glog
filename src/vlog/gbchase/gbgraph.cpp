#include <vlog/gbchase/gbgraph.h>
#include <vlog/gbchase/gbruleexecutor.h>
#include <vlog/gbchase/gbcompositesegment.h>

#include <vlog/support.h>

std::unique_ptr<Literal> GBGraph::createQueryFromNode(
        std::vector<Literal> &outputQueryBody,
        std::vector<size_t> &rangeOutputQueryBody,
        const Rule &rule,
        const std::vector<size_t> &incomingEdges,
        bool incrementCounter) {
    assert(rule.getHeads().size() == 1);

    //Rewrite the rule with fresh variables
    auto oldCounter = counterFreshVarsQueryCont;
    counterFreshVarsQueryCont++;
    Rule newRule = rule.rewriteWithFreshVars(counterFreshVarsQueryCont);
    auto &newBody = newRule.getBody();

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
            const auto &litIncEdge = getNodeHeadQuery(incEdge);
            //Compute the MGU
            std::vector<Substitution> subs;
            auto nsubs = Literal::getSubstitutionsA2B(subs, litIncEdge, l);
            assert(nsubs != -1);

            //First replace all remaining variables with fresh ones
            const auto &bodyIncEdge = getNodeBodyQuery(incEdge);
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
                                    counterFreshVarsQueryCont++, 0)));
                }
            }

            for(auto &litBody : bodyIncEdge) {
                outputQueryBody.push_back(litBody.substitutes(subs));
            }
        }
    }
    assert(idxIncomingEdge == incomingEdges.size());

    if (!incrementCounter)
        counterFreshVarsQueryCont = oldCounter;


    return std::unique_ptr<Literal>(new Literal(newRule.getFirstHead()));
}

void GBGraph::addNodeNoProv(PredId_t predId,
        size_t ruleIdx,
        size_t step,
        std::shared_ptr<const TGSegment> data) {
    if (trackProvenance) {
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

void GBGraph::addNodeProv(PredId_t predid,
        size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> data,
        const std::vector<size_t> &incomingEdges) {
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
#endif

    for(auto n : incomingEdges)
        if (isTmpNode(n))
            throw 10;

    if (queryContEnabled) {
        assert(allRules != NULL && program != NULL && layer != NULL);
        outputNode.incomingEdges = incomingEdges;
        //Create a query and associate it to the node
        auto queryHead = createQueryFromNode(outputNode.queryBody,
                outputNode.rangeQueryBody,
                allRules[ruleIdx],
                incomingEdges);

#ifdef DEBUG
        LOG(INFOL) << "Node " << nodeId << " created with rule " <<
            allRules[ruleIdx].tostring(program, layer);
        LOG(INFOL) << "     QH: " << queryHead->tostring(program, layer);
        std::string nodes = "";
        for(auto i : incomingEdges) {
            nodes += " " + std::to_string(i);
        }
        LOG(INFOL) << "     NO:" << nodes;
        std::string qb = "";
        for(auto &l : outputNode.queryBody) {
            qb += " " + l.tostring(program, layer);
        }
        LOG(INFOL) << "     QB:" << qb;
#endif

        outputNode.queryHead = std::move(queryHead);
        assert(outputNode.queryHead.get() != NULL);
    }

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
        const auto nfields = trackProvenance ? card + 1 : card;
        std::unique_ptr<GBSegmentInserter> rewrittenTuples =
            GBSegmentInserter::getInserter(nfields);
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
                    if (trackProvenance) {
                        row[card] = ~0ul;
                    }
                    rewrittenTuples->add(row.get());
                    countAffectedTuples++;
                    if (countUnaffectedTuples > 0 &&
                            oldTuples.get() == NULL) {
                        //Copy all previous tuples
                        oldTuples = GBSegmentInserter::getInserter(nfields);
                        auto itr2 = data->iterator();
                        for(size_t i = 0; i < countUnaffectedTuples; ++i) {
                            itr2->next();
                            for(int i = 0; i < card; ++i) {
                                row[i] = itr2->get(i);
                            }
                            if (trackProvenance) {
                                row[card] = itr2->getNodeId();
                            }
                            oldTuples->add(row.get());
                        }
                    }
                } else {
                    if (trackProvenance) {
                        row[card] = itr->getNodeId();
                    }
                    if (countAffectedTuples > 0 &&
                            oldTuples.get() == NULL) {
                        oldTuples = GBSegmentInserter::getInserter(nfields);
                    }
                    if (oldTuples.get() != NULL) {
                        oldTuples->add(row.get());
                    }
                    countUnaffectedTuples++;
                }
            }
            if (oldTuples.get() != NULL) {
                auto tuples = oldTuples->getSegment(nodes[nodeId].step,
                        true, 0, trackProvenance);
                nodes[nodeId].setData(tuples);
            } else {
                if (rewrittenTuples->getNRows() == data->getNRows()) {
                    oldTuples = GBSegmentInserter::getInserter(nfields);
                    auto tuples = oldTuples->getSegment(nodes[nodeId].step,
                            true, 0, trackProvenance);
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
                    trackProvenance);
            tuples = tuples->sort()->unique();
            //Retain
            auto retainedTuples = retain(predid, tuples);
            bool nonempty = !(retainedTuples == NULL ||
                    retainedTuples->isEmpty());
            if (nonempty) {
                //Add new nodes
                if (trackProvenance) {
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
        std::shared_ptr<const TGSegment> newtuples) {
    //Special case for unary relations
    if (existuples->getNColumns() == 1) {
        return retainVsNodeFast_one(existuples,
                newtuples);
    } else if (existuples->getNColumns() == 2) {
        return retainVsNodeFast_two(existuples, newtuples);
    } else {
        return retainVsNodeFast_generic(existuples, newtuples);
    }
}

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast_one(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples) {
    if (newtuples->hasColumnarBackend()) {
        ColumnWriter writer;
        bool allNew = true;
        auto newColumn = ((TGSegmentLegacy*)newtuples.get())->getColumn(0);
        if (existuples->hasColumnarBackend()) {
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
                    if (trackProvenance) {
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
                        if (trackProvenance) {
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
            return retainVsNodeFast_generic(existuples, newtuples);
        }
    } else {
        return retainVsNodeFast_generic(existuples, newtuples);
    }
}

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast_two(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples) {
    if (newtuples->hasColumnarBackend()) {
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

                        if (trackProvenance) {
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
                assert(existuples->getProvenanceType() != 2);
                const std::vector<std::pair<Term_t, Term_t>> &t =
                    existuples->getProvenanceType() == 0 ?
                    ((BinaryTGSegment*)existuples.get())->getTuples() :
                    ((BinaryWithConstProvTGSegment*)existuples.get())->getTuples();

                std::vector<std::pair<Term_t, Term_t>> retained = layer.
                    checkNewIn(l1, posColumnsNew, t);
                if (retained.empty()) {
                    return std::shared_ptr<const TGSegment>();
                } else {
                    if (trackProvenance) {
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
    return retainVsNodeFast_generic(existuples, newtuples);
}


std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast_generic(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples) {

    std::unique_ptr<GBSegmentInserter> inserter;
    const uint8_t ncols = newtuples->getNColumns();
    const uint8_t extracol = trackProvenance &&
        newtuples->getProvenanceType() == 2 ? 1 : 0;
    Term_t row[ncols + extracol];
    const bool copyNode = newtuples->getProvenanceType() == 2;

    //Do outer join
    auto leftItr = existuples->iterator();
    auto rightItr = newtuples->iterator();
    bool moveLeftItr = true;
    bool moveRightItr = true;
    bool activeRightValue = false;
    size_t countNew = 0;
    bool isFiltered = false;
    size_t startCopyingIdx = 0;
    while (true) {
        if (moveRightItr) {
            if (rightItr->hasNext()) {
                rightItr->next();
                activeRightValue = true;
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
                if (copyNode) {
                    row[ncols] = rightItr->getNodeId();
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

            //The tuple must be filtered
            if (!isFiltered && countNew > 0) {
                inserter = GBSegmentInserter::getInserter(ncols + extracol);
                //Copy all the previous new tuples in the right iterator
                size_t i = 0;
                auto itrTmp = newtuples->iterator();
                while (i < (startCopyingIdx + countNew) && itrTmp->hasNext()) {
                    itrTmp->next();
                    if (i >= startCopyingIdx) {
                        for(int i = 0; i < ncols; ++i) {
                            row[i] = itrTmp->get(i);
                        }
                        if (copyNode) {
                            row[ncols] = rightItr->getNodeId();
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
            if (copyNode) {
                row[ncols] = rightItr->getNodeId();
            }
            inserter->add(row);
        }
        while (rightItr->hasNext()) {
            rightItr->next();
            for(int i = 0; i < ncols; ++i) {
                row[i] = rightItr->get(i);
            }
            if (copyNode) {
                row[ncols] = rightItr->getNodeId();
            }
            inserter->add(row);
        }
        return inserter->getSegment(newtuples->getNodeId(), true, 0,
                trackProvenance, copyNode);
    } else {
        if (countNew > 0 || activeRightValue) {
            if (startCopyingIdx == 0) {
                //They are all new ...
                return newtuples;
            } else {
                //Remove the initial duplicates
                inserter = GBSegmentInserter::getInserter(ncols + extracol);
                //Copy all the previous new tuples in the right iterator
                size_t i = 0;
                auto itrTmp = newtuples->iterator();
                while (itrTmp->hasNext()) {
                    itrTmp->next();
                    if (i >= startCopyingIdx) {
                        for(int i = 0; i < ncols; ++i) {
                            row[i] = itrTmp->get(i);
                        }
                        if (copyNode) {
                            row[ncols] = rightItr->getNodeId();
                        }
                        inserter->add(row);
                    }
                    i++;
                }
                return inserter->getSegment(newtuples->getNodeId(),
                        true, 0, trackProvenance, copyNode);
            }
        } else {
            //They are all duplicates
            return std::shared_ptr<const TGSegment>();
        }
    }
}

std::shared_ptr<const TGSegment> GBGraph::retain(
        PredId_t p,
        std::shared_ptr<const TGSegment> newtuples) {
    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    if (!areNodesWithPredicate(p)) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto dur = end - start;
        durationRetain += dur;
        return newtuples;
    }
    auto &nodeIdxs = getNodeIDsWithPredicate(p);
    if (cacheRetainEnabled && nodeIdxs.size() > 1) {
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
                auto seg = std::shared_ptr<const TGSegment>(
                        new BinaryTGSegment(tuples, ~0ul, true, 0));
                CacheRetainEntry entry;
                entry.nnodes = nodeIdxs.size();
                entry.seg = seg;
                cacheRetain[p] = entry;

                //std::chrono::steady_clock::time_point end =
                //    std::chrono::steady_clock::now();
                //std::chrono::duration<double, std::milli> dur = end - start;
                //LOG(INFOL) << "RETAIN: " << dur.count();

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
        newtuples = retainVsNodeFast(existingTuples, newtuples);
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
            newtuples = retainVsNodeFast(nodeData, newtuples);
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

std::shared_ptr<const TGSegment> GBGraph::mergeNodes(
        const std::vector<size_t> &nodeIdxs,
        const std::vector<int> &copyVarPos, bool lazyMode) const {

    assert(nodeIdxs.size() > 0);
    auto ncols = copyVarPos.size();
    bool project = ncols < getNodeData(nodeIdxs[0])->getNColumns();
    bool shouldSortAndUnique = project || copyVarPos[0] != 0;

    if (copyVarPos.size() == 1) {
        //Special cases
        if (nodeIdxs.size() == 1) {
            if (!project) {
                return getNodeData(nodeIdxs[0]);
            }

            /*if (lazyMode) {
                return std::shared_ptr<const TGSegment>(
                        new CompositeTGSegment(*this, nodeIdxs, copyVarPos,
                            false, 0, trackProvenance));
            }*/

            auto seg = getNodeData(nodeIdxs[0]);
            if (seg->hasColumnarBackend()) {
                std::vector<std::shared_ptr<Column>> projectedColumns;
                std::vector<std::shared_ptr<Column>> outColumns;
                seg->projectTo(copyVarPos, projectedColumns);
                assert(trackProvenance || projectedColumns.size() == 1);
                if (shouldSortAndUnique) {
                    auto col = projectedColumns[0]->sort();
                    outColumns.push_back(col->unique());
                } else {
                    outColumns.push_back(projectedColumns[0]);
                }
                size_t nrows = outColumns[0]->size();

                if (trackProvenance) {
                    //Add one column with the provenance
                    assert(projectedColumns[1]->isConstant());
                    outColumns.push_back(std::shared_ptr<Column>(
                                new CompressedColumn(seg->getNodeId(), nrows)));

                    return std::shared_ptr<const TGSegment>(
                            new TGSegmentLegacy(
                                outColumns,
                                nrows,
                                true,
                                0,
                                true));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new TGSegmentLegacy(
                                outColumns,
                                nrows,
                                true,
                                0,
                                false));
                }
            }
        }

        if (lazyMode) {
            return std::shared_ptr<const TGSegment>(
                    new CompositeTGSegment(*this, nodeIdxs, copyVarPos,
                        false, 0, trackProvenance));
        }

        if (trackProvenance) {
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
    } else if (copyVarPos.size() == 2) {
        if (nodeIdxs.size() == 1 && !project &&
                copyVarPos[0] == 0 && copyVarPos[1] == 1) {
            return getNodeData(nodeIdxs[0]);
        } else {

            if (lazyMode) {
                return std::shared_ptr<const TGSegment>(
                        new CompositeTGSegment(*this, nodeIdxs, copyVarPos,
                            false, 0, trackProvenance));
            }

            if (trackProvenance) {
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
        assert(copyVarPos.size() > 2);
        if (trackProvenance) {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        } else {
            std::vector<std::vector<Term_t>> tuples(copyVarPos.size());
            for(auto i : nodeIdxs) {
                auto data = getNodeData(i);
                data->appendTo(copyVarPos, tuples);
            }
            size_t nrows = tuples[0].size();
            //Create columns from the content of tuples
            std::vector<std::shared_ptr<Column>> columns;
            for(int i = 0; i < copyVarPos.size(); ++i) {
                columns.push_back(std::shared_ptr<Column>(
                            new InmemoryColumn(tuples[i], true)));
            }
            auto seg = std::shared_ptr<const TGSegment>(
                    new TGSegmentLegacy(columns, nrows));
            if (shouldSortAndUnique) {
                auto sortedSeg = seg->sort();
                return sortedSeg->unique();
            } else {
                return seg;
            }
        }
    }
}

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
    if (trackProvenance) {
        std::vector<size_t> nodes;
        addNodeProv(predId, ~0ul, lastStep, tuples, nodes);
    } else {
        addNodeNoProv(predId, ~0ul, lastStep, tuples);
    }
    return tuples->getNRows();
}
