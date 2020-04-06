#include <vlog/gbchase/gbgraph.h>
#include <vlog/gbchase/gbruleexecutor.h>

#include <vlog/support.h>

void GBGraph::addNode(PredId_t predid, size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> data) {
    auto nodeId = getNNodes();
    nodes.emplace_back();
    GBGraph_Node &outputNode = nodes.back();
    outputNode.predid = predid;
    outputNode.ruleIdx = ruleIdx;
    outputNode.step = step;
    outputNode.data = data;
    pred2Nodes[predid].push_back(nodeId);
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
        //std::unique_ptr<SegmentInserter> rewrittenTuples = std::unique_ptr<
        //    SegmentInserter>(new SegmentInserter(nfields));
        std::unique_ptr<GBSegmentInserter> rewrittenTuples = GBSegmentInserter::getInserter(nfields);
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
                        //oldTuples = std::unique_ptr<SegmentInserter>(
                        //        new SegmentInserter(nfields));
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
                        //oldTuples = std::unique_ptr<SegmentInserter>(
                        //        new SegmentInserter(nfields));
                        oldTuples = GBSegmentInserter::getInserter(nfields);
                    }
                    if (oldTuples.get() != NULL) {
                        oldTuples->add(row.get());
                    }
                    countUnaffectedTuples++;
                }
            }
            if (oldTuples.get() != NULL) {
                //std::shared_ptr<const Segment> seg = oldTuples->getSegment();
                //auto tuples = GBRuleExecutor::fromSeg2TGSeg(
                //        seg , nodes[nodeId].step, true, 0, trackProvenance);
                auto tuples = oldTuples->getSegment(nodes[nodeId].step, true, 0, trackProvenance);
                nodes[nodeId].data = tuples;
            } else {
                if (rewrittenTuples->getNRows() == data->getNRows()) {
                    LOG(ERRORL) << "The case where all the tuples are replaced"
                        " is not (yet) supported";
                    throw 10;
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
            //TODO
            //std::shared_ptr<const Segment> seg = rewrittenTuples->getSegment();
            //auto tuples = GBRuleExecutor::fromSeg2TGSeg(
            //        seg, ~0ul, false, 0, trackProvenance);
            auto tuples = rewrittenTuples->getSegment(~0ul,
                    false,
                    0,
                    trackProvenance);
            tuples = tuples->sort()->unique();
            //Retain
            auto retainedTuples = retain(predid, tuples);
            bool nonempty = !(retainedTuples == NULL || retainedTuples->isEmpty());
            if (nonempty) {
                //Add new nodes
                if (trackProvenance) {
                    auto nodeId = getNNodes();
                    auto dataToAdd = tuples->slice(nodeId, 0, tuples->getNRows());
                    addNode(predid, ruleIdx, step, dataToAdd);
                } else {
                    //Add a single node
                    addNode(predid, ruleIdx, step, retainedTuples);
                }
            }
        }
    }
}

std::shared_ptr<const TGSegment> GBGraph::retainVsNodeFast(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples) {
    std::unique_ptr<GBSegmentInserter> inserter;
    const uint8_t ncols = newtuples->getNColumns();
    const uint8_t extracol = trackProvenance ? 1 : 0;
    Term_t row[ncols + extracol];

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
            inserter->add(row);
        }
        while (rightItr->hasNext()) {
            rightItr->next();
            for(int i = 0; i < ncols; ++i) {
                row[i] = rightItr->get(i);
            }
            inserter->add(row);
        }
        return inserter->getSegment(0, true, 0, trackProvenance);
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
                        inserter->add(row);
                    }
                    i++;
                }
                return inserter->getSegment(0, true, 0, trackProvenance);
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
                if (cacheRetain.count(p)) {
                    cacheRetain[p].seg->appendTo(0, tuples);
                }
                for(size_t i = cacheRetain[p].nnodes; i < nodeIdxs.size(); ++i) {
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
                std::vector<std::pair<Term_t, Term_t>> tuples;
                if (cacheRetain.count(p)) {
                    cacheRetain[p].seg->appendTo(0, 1, tuples);
                }
                for(size_t i = cacheRetain[p].nnodes; i < nodeIdxs.size(); ++i) {
                    getNodeData(nodeIdxs[i])->appendTo(0, 1, tuples);
                }
                std::sort(tuples.begin(), tuples.end());
                auto seg = std::shared_ptr<const TGSegment>(
                        new BinaryTGSegment(tuples, ~0ul, true, 0));
                CacheRetainEntry entry;
                entry.nnodes = nodeIdxs.size();
                entry.seg = seg;
                cacheRetain[p] = entry;
            } else {
                LOG(ERRORL) << "Retain with arity > 2 is not supported";
                throw 10;
            }
        }
        auto existingTuples = cacheRetain[p].seg;
        newtuples = retainVsNodeFast(existingTuples, newtuples);
        if (newtuples == NULL || newtuples->isEmpty()) {
            std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
            auto dur = end - start;
            durationRetain += dur;
            return std::shared_ptr<const TGSegment>();
        }
    } else {
        for(auto &nodeIdx : nodeIdxs) {
            auto nodeData = getNodeData(nodeIdx);
            newtuples = retainVsNodeFast(nodeData, newtuples);
            if (newtuples == NULL || newtuples->isEmpty()) {
                std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
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
        std::vector<int> &copyVarPos) const {
    if (nodeIdxs.size() == 1) {
        size_t idbBodyAtomIdx = nodeIdxs[0];
        auto data  = getNodeData(idbBodyAtomIdx);

        auto ncols = copyVarPos.size(); //ncol=arity of the relation
        bool project = ncols < data->getNColumns();
        if (copyVarPos.size() == 1) {
            if (project) {
                std::vector<Term_t> tuples;
                data->appendTo(copyVarPos[0], tuples);
                std::sort(tuples.begin(), tuples.end());
                auto itr = std::unique(tuples.begin(), tuples.end());
                tuples.erase(itr, tuples.end());
                if (trackProvenance) {
                    return std::shared_ptr<const TGSegment>(
                            new UnaryWithConstProvTGSegment(tuples, idbBodyAtomIdx, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new UnaryTGSegment(tuples, idbBodyAtomIdx, true, 0));
                }
            } else {
                return data; //Can be reused
            }
        } else if (copyVarPos.size() == 2) {
            if (copyVarPos[0] == 0 && copyVarPos[1] == 1) {
                return data; //Can be reused
            } else {
                //Swapped relation
                std::vector<std::pair<Term_t, Term_t>> tuples;
                data->appendTo(1, 0, tuples);
                if (trackProvenance) {
                    return std::shared_ptr<const TGSegment>(
                            new BinaryWithConstProvTGSegment(tuples, idbBodyAtomIdx, false, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new BinaryTGSegment(tuples, idbBodyAtomIdx, false, 0));
                }
            }
        } else {
            if (project) {
                LOG(ERRORL) << "Not implemented";
                throw 10;
            } else {
                return data;
            }
        }
    } else {
        auto ncols = copyVarPos.size(); //ncol=arity of the relation
        bool shouldSortAndUnique = ncols < getNodeData(nodeIdxs[0])->getNColumns() || copyVarPos[0] != 0;
        if (ncols == 1) { //If it has only one column ...
            if (trackProvenance) {
                std::vector<std::pair<Term_t,Term_t>> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    getNodeData(idbBodyAtomIdx)->appendTo(copyVarPos[0], tuples);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                    return std::shared_ptr<const TGSegment>(new UnaryWithProvTGSegment(tuples, ~0ul, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(new UnaryWithProvTGSegment(tuples, ~0ul, false, 0));
                }
            } else {
                std::vector<Term_t> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs)
                    getNodeData(idbBodyAtomIdx)->appendTo(copyVarPos[0], tuples);
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                    return std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, ~0ul, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, ~0ul, false, 0));
                }
            }
        } else if (ncols == 2) { //If it has only two columns ...
            if (trackProvenance) {
                std::vector<BinWithProv> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    getNodeData(idbBodyAtomIdx)->appendTo(copyVarPos[0], copyVarPos[1], tuples);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                    return std::shared_ptr<const TGSegment>(new BinaryWithProvTGSegment(tuples, ~0ul, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(new BinaryWithProvTGSegment(tuples, ~0ul, false, 0));
                }
            } else {
                std::vector<std::pair<Term_t,Term_t>> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    getNodeData(idbBodyAtomIdx)->appendTo(copyVarPos[0], copyVarPos[1], tuples);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                    return std::shared_ptr<const TGSegment>(new BinaryTGSegment(tuples, ~0ul, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(new BinaryTGSegment(tuples, ~0ul, false, 0));
                }
            }
        } else {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }
    }
}

uint64_t GBGraph::removeDuplicatesFromNodes(const std::vector<size_t> &nodeIDs) {
    std::vector<size_t> sortedNodeIDs = nodeIDs;
    std::sort(sortedNodeIDs.begin(), sortedNodeIDs.end());

    uint64_t removedTuples = 0;

    //Merge the tuples into a single segment
    assert(nodeIDs.size() > 0);
    auto data = getNodeData(sortedNodeIDs[0]);
    auto card = data->getNColumns();
    std::vector<int> copyVarPos;
    for(size_t i = 0; i < card; ++i)
        copyVarPos.push_back(i);
    auto tuples = mergeNodes(sortedNodeIDs, copyVarPos);
    tuples = tuples->sort();
    tuples = tuples->unique();

    //TODO: Retain only the new tuples

    //Replace the content of the nodes
    tuples = tuples->sortByProv();
    auto itr = tuples->iterator();
    size_t begin = 0;
    size_t end = 0;
    size_t prevNode = ~0ul;
    auto currentNodeID = sortedNodeIDs.begin();
    while (itr->hasNext()) {
        if (itr->getNodeId() != prevNode) {
            if (begin < end) {
                while (prevNode > *currentNodeID) {
                    LOG(ERRORL) << "Some nodes must be completely emptied. This"
                        " operation is not currently supported";
                    throw 10;
                }
                currentNodeID++;
                nodes[prevNode].data = tuples->slice(prevNode, begin, end);
            }
            prevNode = itr->getNodeId();
        }
        end++;
    }
    if (begin < end) {
        while (prevNode > *currentNodeID) {
            LOG(ERRORL) << "Some nodes must be completely emptied. This"
                " operation is not currently supported";
            throw 10;
        }
        currentNodeID++;
        nodes[prevNode].data = tuples->slice(prevNode, begin, end);
    }

    return removedTuples;
}
