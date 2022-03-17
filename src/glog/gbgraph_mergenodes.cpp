#include <glog/gbgraph.h>
#include <glog/gblegacysegment.h>
#include <glog/gbcompositesegment.h>

std::shared_ptr<const TGSegment> GBGraph::mergeNodes_special_unary1(
        std::shared_ptr<const TGSegment> seg,
        const std::vector<size_t> &nodeIdxs,
        const std::vector<int> &copyVarPos, bool lazyMode,
        bool replaceOffsets,
        bool removeDuplicates) const {
    assert(nodeIdxs.size() == 1);
    assert(!replaceOffsets || !isTmpNode(nodeIdxs[0])); //This procedure should not be called if there is a temporary nodes
    //because temporary nodes do not create nodes with columnar layouts. If this condition is false, then we must fix the code
    //below to do the replacement only if the node is not temporary
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
            if (removeDuplicates)
                outColumns[0] = sortedCol->unique();
            else
                outColumns[0] = sortedCol;
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
        bool replaceOffsets,
        bool removeDuplicates) const {
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
                bool isSorted = replaceOffsets;
                Term_t prevElement = 0;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    size_t start = tuples.size();
                    getNodeData(idbBodyAtomIdx)->appendTo(
                            copyVarPos[0], tuples);
                    if (replaceOffsets) {
                        if (!isTmpNode(idbBodyAtomIdx)) {
                            for(size_t i = start; i < tuples.size(); ++i) {
                                tuples[i].prov = i - start;
                                if (isSorted && tuples[i].first < prevElement) {
                                    isSorted = false;
                                } else {
                                    prevElement = tuples[i].first;
                                }
                            }
                        }
                    }
                }

                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    return std::shared_ptr<const TGSegment>(
                            new UnaryWithFullProvTGSegment(tuples, ~0ul,
                                true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new UnaryWithFullProvTGSegment(tuples, ~0ul,
                                isSorted, 0));
                }
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
                if (removeDuplicates) {
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                }
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
            if (removeDuplicates) {
                auto itr = std::unique(tuples.begin(), tuples.end());
                tuples.erase(itr, tuples.end());
            }
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
        const std::vector<Term_t> &filterConstants,
        const std::vector<int> &copyVarPos, bool lazyMode,
        bool replaceOffsets, bool removeDuplicates) const {

    if (lazyMode && provenanceType == NOPROV) {
        LOG(DEBUGL) << "Lazymode is deactivated with NOPROV";
        lazyMode = false;
    }

    assert(nodeIdxs.size() > 0);
    auto ncols = copyVarPos.size();
    bool project = ncols > 0 && ncols < getNodeData(nodeIdxs[0])->getNColumns();
    bool shouldSortAndUnique = (project ||
            (copyVarPos.size() > 0 && copyVarPos[0] != 0));

    if (filterConstants.empty() && copyVarPos.size() == 1) {
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
                        replaceOffsets,
                        removeDuplicates);
            }
        }

        if (lazyMode) {
            return std::shared_ptr<const TGSegment>(
                    new CompositeTGSegment(*this, nodeIdxs, copyVarPos,
                        false, 0, getSegProvenanceType(nodeIdxs.size() > 1),
                        false, false, replaceOffsets));
        }

        return mergeNodes_special_unary2(
                nodeIdxs,
                copyVarPos,
                lazyMode,
                replaceOffsets,
                removeDuplicates);

    } else if (filterConstants.empty() && copyVarPos.size() == 2) {
        //Check that node of the nodes is temporary. If they are,
        //then I must fix the replaceOffset procedure, which should not apply
        //for temporaty nodes
        for(auto ni : nodeIdxs) {
            if (isTmpNode(ni))
                throw 10;
        }

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
                        if (removeDuplicates) {
                            auto itr = std::unique(tuples.begin(), tuples.end());
                            tuples.erase(itr, tuples.end());
                        }
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
                    if (removeDuplicates) {
                        auto itr = std::unique(tuples.begin(), tuples.end());
                        tuples.erase(itr, tuples.end());
                    }
                    return std::shared_ptr<const TGSegment>(
                            new BinaryTGSegment(tuples, ~0ul, true, 0));
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new BinaryTGSegment(tuples, ~0ul, false, 0));
                }
            }
        }
    } else {
        return mergeNodes_general(nodeIdxs, filterConstants, copyVarPos,
                lazyMode, replaceOffsets, shouldSortAndUnique,
                shouldSortAndUnique && removeDuplicates);
    }
}

std::shared_ptr<const TGSegment> GBGraph::mergeNodes_general(
        const std::vector<size_t> &nodeIdxs,
        const std::vector<Term_t> &filterConstants,
        const std::vector<int> &copyVarPos,
        bool lazyMode,
        bool replaceOffsets,
        bool shouldSort,
        bool shouldRemoveDuplicates) const {
    //Check that node of the nodes is temporary. If they are, then I must
    //fix the replaceOffset procedure, which should not apply for temporaty nodes
    for(auto ni : nodeIdxs) {
        if (isTmpNode(ni))
            throw 10;
    }

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
    if (filterConstants.empty()) {
        for(auto i : nodeIdxs) {
            auto data = getNodeData(i);
            if (shouldTrackProvenance()) {
                if (!replaceOffsets &&
                        provenanceType == GBGraph::ProvenanceType::FULLPROV) {
                    assert(data->getNOffsetColumns() == maxOffset + 1);
                    //The other case is not implemented
                }
                auto start = tuples.back().size();
                data->appendTo(copyVarPos, tuples, true);
                if (replaceOffsets) {
                    if (data->getProvenanceType() == SEG_FULLPROV)
                    {
                        //Otherwise, the last column contains something else
                        for(size_t i = start; i < tuples.back().size(); ++i) {
                            tuples.back()[i] = i - start;
                        }
                    } else if (data->getProvenanceType() == SEG_SAMENODE) {
                        //Do nothing, I don't have to record the offset
                    } else {
                        //I suspect it's the same as for SEG_SAMENODE, but I'm not sure
                        throw 10; //not implemented
                    }
                }
            } else {
                data->appendTo(copyVarPos, tuples, false);
            }
        }
    } else {
        //Filter rows based on the constants
        for(auto i : nodeIdxs) {
            auto data = getNodeData(i);
            auto itr = data->iterator();
            size_t idxRow = 0;
            while (itr->hasNext()) {
                itr->next();
                if (itr->getNProofs() != 1 && !replaceOffsets)
                    throw 10; //Otherwise, it should be expanded
                bool ok = true;
                for(size_t j = 0; j < filterConstants.size(); ++j) {
                    if (filterConstants[j] != ~0ul) {
                        if (itr->get(j) != filterConstants[j]) {
                            ok = false;
                            break;
                        }
                    }
                }
                if (ok) {
                    for(size_t j = 0; j < copyVarPos.size(); ++j) {
                        tuples[j].push_back(itr->get(copyVarPos[j]));
                    }
                    if (shouldTrackProvenance()) {
                        tuples[copyVarPos.size()].push_back(itr->getNodeId());
                        if (provenanceType == FULLPROV) {
                            if (replaceOffsets) {
                                assert(tuples.size() == copyVarPos.size() + 2);
                                tuples.back().push_back(idxRow);
                            } else {
                                //Copy the remaining offsets
                                for(size_t j = 1; j < maxOffset; ++j) {
                                    tuples[copyVarPos.size() + j].push_back(
                                            itr->getProvenanceOffset(0, j-1));
                                }
                            }
                        }
                    }
                }
                idxRow++;
            }
        }
    }
    size_t nrows = tuples[0].size();
    if (nrows == 0)
    {
        return std::shared_ptr<const TGSegment>();
    } else {
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
                        new InmemoryColumn(tuples[copyVarPos.size()], true)));
            if (provenanceType == GBGraph::ProvenanceType::FULLPROV) {
                if (replaceOffsets) {
                    columns.push_back(std::shared_ptr<Column>(
                                new InmemoryColumn(tuples.back(), true)));
                } else {
                    throw 10; //Not implemented. Here we should copy all the offset columns
                }
            }
        }
        seg = std::shared_ptr<const TGSegment>(
                new TGSegmentLegacy(columns, nrows, false,
                    0, getSegProvenanceType(nodeIdxs.size() > 1),
                    columns.size() - copyVarPos.size()));
        if (shouldSort) {
            auto sortedSeg = seg->sort();
            if (shouldRemoveDuplicates)
                return sortedSeg->unique();
            else
                return sortedSeg;
        } else {
            return seg;
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
    if (shouldTrackProvenance()) {
        std::vector<size_t> nodes;
        addNodeProv(predId, ~0ul, lastStep, tuples, nodes);
    } else {
        addNodeNoProv(predId, ~0ul, lastStep, tuples);
    }
    return tuples->getNRows();
}
