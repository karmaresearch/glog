#include <glog/gbgraph.h>
#include <glog/gbsegmentinserter.h>

void GBGraph::retainFromDerivationTree_getNodes(size_t nodeId,
        size_t offsetNodeId,
        PredId_t predId,
        std::vector<GBGraph_PathDerivationTree> currentPath,
        std::vector<std::vector<GBGraph_PathDerivationTree>> &out) {
    assert(nodeId < getNNodes());
    auto r = getNodeRuleIdx(nodeId);
    if (r == ~0ul) {
        //Node is a node that has been manually added
        return;
    }
    GBGraph_PathDerivationTree t;
    t.offset = offsetNodeId;
    t.nodeId = nodeId;
    currentPath.push_back(t);
    Rule &rule = allRules[r];
    auto &heads = rule.getHeads();
    assert(heads.size() == 1);
    auto id = heads[0].getPredicate().getId();
    if (id == predId) {
        currentPath.back().tocheck = true;
        out.push_back(currentPath);
    }
    auto &inc = getNodeIncomingEdges(nodeId);
    for(size_t i = 0; i < inc.size(); ++i) {
        auto pNodeId = inc[i];
        if (pNodeId != ~0ul) {
            retainFromDerivationTree_getNodes(pNodeId, i, predId, currentPath,
                    out);
        }
    }
}

bool GBGraph::GBGraph_PathDerivationTree::pathSorter(
        const std::vector<GBGraph_PathDerivationTree> &a,
        const std::vector<GBGraph_PathDerivationTree> &b) {
    for(size_t i = 0; i < a.size() && i < b.size(); ++i)
    {
        if (a[i].nodeId != b[i].nodeId) {
            return a[i].nodeId > b[i].nodeId;
        } else {
            if (a[i].offset != b[i].offset) {
                return a[i].offset > b[i].offset;
            }
        }
    }
    return a.size() > b.size();
}

void GBGraph::GBGraph_PathDerivationTree::removePrefixes(
        std::vector<std::vector<GBGraph_PathDerivationTree>> &paths)
{
    std::sort(paths.begin(), paths.end(), pathSorter);
    for(int64_t j = paths.size() - 1; j >= 1; j--) {
        auto &cp = paths[j];
        auto &pp = paths[j-1];
        bool isPrefix = true;
        for(size_t m = 0; m < cp.size(); ++m) {
            if (m >= pp.size()) {
                isPrefix = false;
                break;
            }
            if (cp[m].nodeId != pp[m].nodeId || cp[m].offset != pp[m].offset) {
                isPrefix = false;
                break;
            }
        }
        if (isPrefix && pp.size() > cp.size()) {
            paths.erase(paths.begin() + j);
        }
    }
}

std::shared_ptr<const TGSegment> GBGraph::retainFromDerivationTree(
        PredId_t p,
        size_t ruleIdx,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<size_t> derivationNodes) {

    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();

    if (!derivationNodes.empty() &&
            derivationNodes.size() != newtuples->getNOffsetColumns() - 1) {
        throw 10;
    }

    std::vector<std::vector<GBGraph_PathDerivationTree>> nodesToCheckAgainst;
    for(size_t i = 0; i < derivationNodes.size(); ++i) {
        size_t node = derivationNodes[i];
        if (node == ~0ul) {
            continue;
        } else {
            //Search for nodes with predid=p
            retainFromDerivationTree_getNodes(node, i, p,
                    std::vector<GBGraph_PathDerivationTree>(),
                    nodesToCheckAgainst);
        }
    }

    //Sort the nodes and retain only the longest paths
    GBGraph_PathDerivationTree::removePrefixes(nodesToCheckAgainst);

    std::vector<bool> toBeDeleted(newtuples->getNRows());

    for(auto &pathNodeToCheckAgainst : nodesToCheckAgainst) {
        auto offsetColumn = pathNodeToCheckAgainst.front().offset;

        //Create offsetMap
        std::vector<
            std::tuple<size_t, size_t,std::vector<size_t>::iterator>> pathOffsetMap;
        std::vector<size_t> rowStats(newtuples->getNRows());
        size_t idx = 0;
        auto itr = newtuples->iterator();
        while (itr->hasNext()) {
            itr->next();
            if (!toBeDeleted[idx]) {
                rowStats[idx] = itr->getNProofs();
                pathOffsetMap.push_back(std::make_tuple(
                            idx, itr->getProvenanceOffset(0, offsetColumn),
                            rowStats.begin() + idx));
                for(size_t j = 1; j < itr->getNProofs(); ++j) {
                    pathOffsetMap.push_back(std::make_tuple(
                                idx, itr->getProvenanceOffset(j, offsetColumn),
                                rowStats.begin() + idx));
                }
            }
            idx++;
        }
        if (pathOffsetMap.empty()) {
            break;
        }

        for(size_t j = 0; j < pathNodeToCheckAgainst.size(); ++j) {
            auto n = pathNodeToCheckAgainst[j].nodeId;
            if (pathNodeToCheckAgainst[j].tocheck) {
                auto existData = getNodeData(n);
                size_t card = existData->getNColumns();
                for(int64_t m = pathOffsetMap.size() - 1; m >= 0; m--)
                {
                    const auto &p = pathOffsetMap[m];
                    auto r = std::get<1>(p);
                    bool equal = true;
                    for(size_t z = 0; z < card; ++z) {
                        if (existData->getValueAtRow(r, z) !=
                                newtuples->getValueAtRow(std::get<0>(p), z)) {
                            equal = false;
                            break;
                        }
                    }
                    if (equal) {
                        auto &counter = std::get<2>(p);
                        assert((*counter) > 0);
                        (*counter)--;
                        if (*counter == 0) {
                            toBeDeleted[std::get<0>(p)] = true;
                        }
                        pathOffsetMap.erase(pathOffsetMap.begin() + m);
                    }
                }
            } else {
                //Take the offset specified in the next node
                size_t off = pathNodeToCheckAgainst[j+1].offset;
                auto d = getNodeData(n);
                for(int64_t m = pathOffsetMap.size() - 1; m >= 0; m--) {
                    auto oldOff = std::get<1>(pathOffsetMap[m]);
                    auto value = d->getOffsetAtRow(oldOff, 0, off);
                    std::get<1>(pathOffsetMap[m]) = value;
                    for (size_t j = 1; j < d->getNProofsAtRow(oldOff); ++j) {
                        std::get<2>(pathOffsetMap[m])++;
                        value = d->getOffsetAtRow(oldOff, j, off);
                        pathOffsetMap.push_back(
                                std::make_tuple(std::get<0>(pathOffsetMap[m]),
                                    value,
                                    std::get<2>(pathOffsetMap[m])));
                    }
                }
            }
        }
    }

    std::vector<size_t> newRowIdxs;
    for(size_t j = 0; j < toBeDeleted.size(); ++j)
        if (!toBeDeleted[j])
            newRowIdxs.push_back(j);
    std::shared_ptr<const TGSegment> out;
    if (newRowIdxs.size() == newtuples->getNRows()) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto dur = end - start;
        durationRetain += dur;
        out = newtuples;
    } else {
        /*const uint8_t ncols = newtuples->getNColumns();
        const size_t extracols = newtuples->getNOffsetColumns();
        auto inserter = GBSegmentInserter::getInserter(ncols + extracols,
                extracols, false);
        Term_t row[ncols + extracols];
        size_t rowIdx = 0;
        size_t curDuplRowIdx = 0;
        auto itr = newtuples->iterator();
        while (itr->hasNext()) {
            itr->next();
            if (curDuplRowIdx < duplicateRowIdxs.size() &&
                    rowIdx == duplicateRowIdxs[curDuplRowIdx]) {
                //skip it
                curDuplRowIdx++;
            } else {
                for(size_t i = 0; i < ncols; ++i) {
                    row[i] = itr->get(i);
                }
                row[ncols] = itr->getNodeId();
                for(size_t i = 1; i < extracols; ++i) {
                    row[ncols + i] = itr->getProvenanceOffset(0, i-1);
                }

                //add the row
                inserter->add(row);
            }
            rowIdx++;
        }*/

        if (newRowIdxs.empty()) {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            auto dur = end - start;
            durationRetain += dur;
            out = std::shared_ptr<const TGSegment>();
        } else {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            auto dur = end - start;
            durationRetain += dur;
            out = newtuples->shuffle(newRowIdxs);
        }
    }
    return out;
}

std::shared_ptr<const TGSegment> GBGraph::retain(
        PredId_t p,
        std::shared_ptr<const TGSegment> newtuples,
        std::vector<std::shared_ptr<Column>> &derivationNodes,
        bool inputWithDuplicates) {
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
    if (inputWithDuplicates) {
        newtuples = newtuples->unique();
    }

    std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    if (!areNodesWithPredicate(p)) {
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        auto dur = end - start;
        durationRetain += dur;
        return newtuples;
    }
    auto &nodeIdxs = getNodeIDsWithPredicate(p);
    assert(nodeIdxs.size() > 0);
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
            newtuples = retainVsNodeFast(nodeData, newtuples, derivationNodes);
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
            assert(itr->getNProofs() == 1);
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
                    row[n_columns + i] = itr->getProvenanceOffset(0, i-1);
                }
            } else {
                row[n_columns] = counter;
            }
            inserter->add(row.get());
        }
        counter += (size_t)1 << 40;
    }

    auto seg = inserter->getSegment(~0ul, false, 0,
            getSegProvenanceType(true), n_offsetcolumns);
    auto sortedSeg = seg->sort();
    auto uniqueSeg = sortedSeg->unique();
    auto retainedSeg = retain(predId, uniqueSeg);
    if (retainedSeg.get() == NULL)
        return;
    std::shared_ptr<const TGSegment> toBeAddedSeg;
    if (!retainedSeg->isNodeConstant()) {
        toBeAddedSeg = retainedSeg->sortByProv();
    } else {
        toBeAddedSeg = retainedSeg;
    }
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
        assert(itr->getNProofs() == 1);
        for(int i = 0; i < n_columns; ++i) {
            row[i] = itr->get(i);
        }
        if (provenanceType != NOPROV) {
            auto nodeId = itr->getNodeId() & 0xFFFFFFFFFFl; //Removing the initial prefix
            if (nodeId == 0xFFFFFFFFFFl)
                nodeId = ~0ul;
            row[n_columns] = nodeId;
            for(int i = 1; i < n_offsetcolumns; ++i) {
                row[n_columns + i] = itr->getProvenanceOffset(0, i - 1);
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
    }
    return inserter->getSegment(currentNode, true, 0,
            node.data->getProvenanceType(), n_offsetcolumns);
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
                    assert(rightItr->getNProofs() == 1);
                    for(size_t i = 1; i < extracols; ++i) {
                        row[ncols + i] = rightItr->getProvenanceOffset(0, i-1);
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
                            assert(itrTmp->getNProofs() == 1);
                            for(size_t i = 1; i < extracols; ++i) {
                                row[ncols + i] = itrTmp->getProvenanceOffset(0, i-1);
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
                assert(rightItr->getNProofs() == 1);
                row[ncols] = rightItr->getNodeId();
                for(size_t i = 1; i < extracols; ++i) {
                    row[ncols + i] = rightItr->getProvenanceOffset(0, i-1);
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
                assert(rightItr->getNProofs() == 1);
                row[ncols] = rightItr->getNodeId();
                for(size_t i = 1; i < extracols; ++i) {
                    row[ncols + i] = rightItr->getProvenanceOffset(0, i-1);
                }
            }
            inserter->add(row);
        }

        return inserter->getSegment(newtuples->getNodeId(), true, 0,
                newtuples->getProvenanceType(), extracols);
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
                            assert(itrTmp->getNProofs() == 1);
                            row[ncols] = itrTmp->getNodeId();
                            for(size_t i = 1; i < extracols; ++i) {
                                row[ncols + i] = itrTmp->getProvenanceOffset(0, i-1);
                            }
                        }
                        inserter->add(row);
                    }
                    i++;
                }

                return inserter->getSegment(newtuples->getNodeId(),
                        true, 0, newtuples->getProvenanceType(), extracols);
            }
        } else {
            //They are all duplicates. Leave the derivationNodes unchanged, anyway
            //the node won't be added
            return std::shared_ptr<const TGSegment>();
        }
    }
}
