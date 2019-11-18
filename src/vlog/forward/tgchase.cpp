#include <vlog/tgchase.h>
#include <vlog/tgsegmentcache.h>

TGChase::TGChase(EDBLayer &layer, Program *program) : layer(layer),
    program(program),
    durationMergeSort(0),
    durationJoin(0),
    durationRetain(0),
    durationCreateHead(0),
    durationFirst(0),
    trackProvenance(false)
{
}

Program *TGChase::getProgram() {
    return program;
}

EDBLayer &TGChase::getEDBLayer() {
    return layer;
}

std::shared_ptr<const TGSegment> fromSeg2TGSeg(std::shared_ptr<const Segment> seg, size_t nodeId, bool isSorted, uint8_t sortedField) {
    auto ncols = seg->getNColumns();
    if (ncols == 1) {
        const std::vector<Term_t> &orig = seg->getColumn(0)->getVectorRef();
        std::vector<Term_t> tuples(orig);
        return std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, nodeId, isSorted, sortedField));
    } else if (ncols == 2) {
        auto &col1 = seg->getColumn(0)->getVectorRef();
        auto &col2 = seg->getColumn(1)->getVectorRef();
        size_t nrows = col1.size();
        std::vector<std::pair<Term_t, Term_t>> out(nrows);
        for(size_t i = 0; i < nrows; ++i) {
            out[i].first = col1[i];
            out[i].second = col2[i];
        }
        return std::shared_ptr<const TGSegment>(new BinaryTGSegment(out, nodeId, isSorted, sortedField));
    } else {
        std::vector<std::shared_ptr<Column>> columns;
        for(int i = 0; i < ncols; ++i) {
            columns.push_back(seg->getColumn(i));
        }
        return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns, seg->getNRows(), isSorted, sortedField));
    }
}

std::shared_ptr<const Segment> fromTGSeg2Seg(std::shared_ptr<const TGSegment> seg) {
    auto nrows = seg->getNRows();
    auto ncols = seg->getNColumns();
    std::vector<std::vector<Term_t>> matrix(ncols);
    auto itr = seg->iterator();
    while (itr->hasNext()) {
        itr->next();
        for(int i = 0; i < ncols; ++i) {
            matrix[i].push_back(itr->get(i));
        }
    }

    std::vector<std::shared_ptr<Column>> columns;
    for(int i = 0; i < ncols; ++i) {
        columns.push_back(std::shared_ptr<Column>(new InmemoryColumn(matrix[i])));
    }
    return std::shared_ptr<const Segment>(new
            Segment(columns.size(), columns));
}


void TGChase::run() {
    initRun();
    size_t nnodes = 0;
    size_t step = 0;
    rules = program->getAllRules();

    do {
        step++;
        LOG(INFOL) << "Step " << step;
        currentIteration = step;
        nnodes = nodes.size();

        std::vector<TGChase_SuperNode> newnodes;
        for (size_t ruleIdx = 0; ruleIdx < rules.size(); ++ruleIdx) {
            auto &rule = rules[ruleIdx];
            std::vector<std::vector<size_t>> nodesForRule;
            bool empty = false;
            for (auto &bodyAtom : rule.getBody()) {
                Predicate pred = bodyAtom.getPredicate();
                if (pred.getType() != EDB) {
                    if (!pred2Nodes.count(pred.getId())) {
                        empty = true;
                        break;
                    }
                    nodesForRule.push_back(pred2Nodes[pred.getId()]);
                }
            }
            if (empty) continue;

            if (rule.getNIDBPredicates() == 0) {
                //It's a rule with only EDB body atoms. Create a single node
                //I only execute these rules in the first step
                if (step == 1) {
                    newnodes.emplace_back();
                    TGChase_SuperNode &newnode = newnodes.back();
                    newnode.ruleIdx = ruleIdx;
                    newnode.step = step;
                    newnode.incomingEdges = std::vector<std::vector<size_t>>();
                }
            } else {
                size_t prevstep = 0;
                if (step > 1) prevstep = step - 1;
                for(int pivot = 0; pivot < nodesForRule.size(); ++pivot) {
                    //First consider only combinations where at least one node
                    //is in the previous level
                    std::vector<std::vector<size_t>> acceptableNodes;
                    bool acceptableNodesEmpty = false;
                    for(int j = 0; j < nodesForRule.size(); ++j) {
                        std::vector<size_t> selection;
                        if (j < pivot) {
                            for(auto &nodeId : nodesForRule[j]) {
                                const auto &node = nodes[nodeId];
                                if (node.step < prevstep) {
                                    selection.push_back(nodeId);
                                } else if (node.step >= prevstep) {
                                    break;
                                }
                            }
                        } else if (j == pivot) {
                            for(auto &nodeId : nodesForRule[j]) {
                                const auto &node = nodes[nodeId];
                                if (node.step == prevstep) {
                                    selection.push_back(nodeId);
                                }
                            }
                        } else {
                            selection = nodesForRule[j];
                        }
                        if (selection.empty()) {
                            acceptableNodesEmpty = true;
                            break;
                        }
                        acceptableNodes.push_back(selection);
                    }
                    if (acceptableNodesEmpty) continue;

                    newnodes.emplace_back();
                    TGChase_SuperNode &newnode = newnodes.back();
                    newnode.ruleIdx = ruleIdx;
                    newnode.step = step;
                    newnode.incomingEdges = acceptableNodes;
                }
            }
        }

        //Execute the rule associated to the node
        auto nnodes = nodes.size();

        auto nodesToProcess = newnodes.size();
        LOG(INFOL) << "Nodes to process " << nodesToProcess;
        for(size_t idxNode = 0; idxNode < nodesToProcess; ++idxNode) {
            executeRule(newnodes[idxNode]);
        }

        if (nnodes != nodes.size()) {
            auto derivedTuples = 0;
            for(size_t idx = nnodes; idx < nodes.size(); ++idx) {
                auto nrows = nodes[idx].data->getNRows();
                derivedTuples += nrows;
            }
            LOG(INFOL) << "Derived Tuples: " << derivedTuples;
        }
    } while (nnodes != nodes.size());

    SegmentCache::getInstance().clear();
    LOG(INFOL) << "Time first (ms): " << durationFirst.count();
    LOG(INFOL) << "Time mergesort (ms): " << durationMergeSort.count();
    LOG(INFOL) << "Time joins (ms): " << durationJoin.count();
    LOG(INFOL) << "Time head (ms): " << durationCreateHead.count();
    LOG(INFOL) << "Time retain (ms): " << durationRetain.count();
    stopRun();
}

void TGChase::computeVarPos(std::vector<size_t> &leftVars,
        int bodyAtomIdx,
        const std::vector<Literal> &bodyAtoms,
        const Literal &head,
        std::pair<int, int> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight) {
    std::set<int> futureOccurrences;
    for(size_t i = bodyAtomIdx + 1; i < bodyAtoms.size(); ++i) {
        auto &lit = bodyAtoms[i];
        auto vars = lit.getAllVars();
        for(auto v : vars)
            futureOccurrences.insert(v);
    }
    auto headVars = head.getAllVars();
    for(auto v : headVars)
        futureOccurrences.insert(v);

    auto &rightBodyAtom = bodyAtoms[bodyAtomIdx];
    auto rightVars = rightBodyAtom.getAllVars();
    std::set<int> addedVars;
    if (bodyAtomIdx > 0) {
        //First check whether variables in the left are used in the next body
        //atoms or in the head
        for(size_t i = 0; i < leftVars.size(); ++i) {
            if (futureOccurrences.count(leftVars[i])) {
                copyVarPosLeft.push_back(i);
                addedVars.insert(leftVars[i]);
            }
        }
        //Search for the join variable
        bool found = false;
        for(int i = 0; i < rightVars.size(); ++i) {
            for(int j = 0; j < leftVars.size(); ++j) {
                if (rightVars[i] == leftVars[j]) {
                    if (found) {
                        LOG(ERRORL) << "Only one variable can be joined";
                        throw 10;
                    }
                    joinVarPos.first = j;
                    joinVarPos.second = i;
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            LOG(ERRORL) << "Could not find join variable";
            throw 10;
        }
    }

    //Copy variables from the right atom
    for(size_t i = 0; i < rightVars.size(); ++i) {
        if (futureOccurrences.count(rightVars[i]) && !addedVars.count(rightVars[i])) {
            copyVarPosRight.push_back(i);
            addedVars.insert(rightVars[i]);
        }
    }
}

std::shared_ptr<const TGSegment> TGChase::mergeNodes(
        std::vector<size_t> &nodeIdxs,
        std::vector<int> &copyVarPos) {
    if (nodeIdxs.size() == 1) {
        size_t idbBodyAtomIdx = nodeIdxs[0];
        auto bodyNode = nodes[idbBodyAtomIdx];
        return processFirstAtom_IDB(bodyNode.data, copyVarPos, idbBodyAtomIdx);
    } else {
        auto ncols = copyVarPos.size(); //ncol=arity of the relation
        bool shouldSortAndUnique = ncols < nodes[nodeIdxs[0]].data->getNColumns();
        if (ncols == 1) {
            if (trackProvenance) {
                std::vector<std::pair<Term_t,Term_t>> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    nodes[idbBodyAtomIdx].data->appendToWithProv(copyVarPos[0], tuples, idbBodyAtomIdx);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                }
                return std::shared_ptr<const TGSegment>(new UnaryWithProvTGSegment(tuples, ~0ul));

            } else {
                std::vector<Term_t> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs)
                    nodes[idbBodyAtomIdx].data->appendTo(copyVarPos[0], tuples);
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                }
                return std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, ~0ul, false, 0));
            }
        } else if (ncols == 2) {
            if (trackProvenance) {
                std::vector<BinWithProv> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    nodes[idbBodyAtomIdx].data->appendToWithProv(copyVarPos[0], copyVarPos[1], tuples, idbBodyAtomIdx);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                }
                return std::shared_ptr<const TGSegment>(new BinaryWithProvTGSegment(tuples, ~0ul));
            } else {
                std::vector<std::pair<Term_t,Term_t>> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    nodes[idbBodyAtomIdx].data->appendTo(copyVarPos[0], copyVarPos[1], tuples);
                }
                if (shouldSortAndUnique) {
                    std::sort(tuples.begin(), tuples.end());
                    auto itr = std::unique(tuples.begin(), tuples.end());
                    tuples.erase(itr, tuples.end());
                }
                return std::shared_ptr<const TGSegment>(new BinaryTGSegment(tuples, ~0ul, false, 0));
            }
        } else {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }
    }
}

std::shared_ptr<const TGSegment> TGChase::processFirstAtom_IDB(
        std::shared_ptr<const TGSegment> &input,
        std::vector<int> &copyVarPos,
        size_t nodeId) {
    auto ncols = copyVarPos.size(); //ncol=arity of the relation
    bool project = ncols < input->getNColumns();
    if (copyVarPos.size() == 1) {
        if (project) {
            std::vector<Term_t> tuples;
            input->appendTo(copyVarPos[0], tuples);
            std::sort(tuples.begin(), tuples.end());
            auto itr = std::unique(tuples.begin(), tuples.end());
            tuples.erase(itr, tuples.end());
            if (trackProvenance) {
                return std::shared_ptr<const TGSegment>(new UnaryWithConstProvTGSegment(tuples, nodeId));
            } else {
                return std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, nodeId, true, 0));
            }
        } else {
            return input;
        }
    } else if (copyVarPos.size() == 2) {
        if (copyVarPos[0] == 0 && copyVarPos[1] == 1) {
            return input; //Can be reused
        } else {
            std::vector<std::pair<Term_t, Term_t>> tuples;
            input->appendTo(1, 0, tuples);
            if (trackProvenance) {
                return std::shared_ptr<const TGSegment>(new BinaryWithConstProvTGSegment(tuples, nodeId));
            } else {
                return std::shared_ptr<const TGSegment>(new BinaryTGSegment(tuples, nodeId, false, 0));
            }
        }
    } else {
        LOG(ERRORL) << "Not implemented";
        throw 10;
    }
}

std::shared_ptr<const TGSegment> TGChase::processFirstAtom_EDB(
        const Literal &atom,
        std::vector<int> &copyVarPos) {
    PredId_t p = atom.getPredicate().getId();
    if (!edbTables.count(p)) {
        auto edbTable = layer.getEDBTable(p);
        edbTables[p] = edbTable;
    }
    //Get the columns
    auto &table = edbTables[p];
    if (!table->useSegments()) {
        LOG(ERRORL) << "EDB table not supported";
        throw 10;
    }
    if (copyVarPos.size() != atom.getNVars()) {
        LOG(ERRORL) << "EDB table not supported";
        throw 10;
    }
    auto seg = table->getSegment();
    std::vector<std::shared_ptr<Column>> columns;
    for(int i = 0; i < copyVarPos.size(); ++i) {
        int pos = copyVarPos[i];
        auto col = seg->getColumn(pos);
        columns.push_back(col);
    }

    auto nrows = columns[0]->size();
    return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns, nrows, true));
}

int TGChase::cmp(std::unique_ptr<TGSegmentItr> &inputLeft,
        std::unique_ptr<TGSegmentItr> &inputRight,
        std::pair<int, int> &joinVarPos) {
    const auto valLeft = inputLeft->get(joinVarPos.first);
    const auto valRight = inputRight->get(joinVarPos.second);
    if (valLeft < valRight)
        return -1;
    else if (valLeft > valRight)
        return 1;
    else
        return 0;
}

int TGChase::cmp(std::unique_ptr<TGSegmentItr> &inputLeft,
        std::unique_ptr<TGSegmentItr> &inputRight) {
    auto n = inputLeft->getNFields();
    for(size_t i = 0; i < n; ++i) {
        auto vl = inputLeft->get(i);
        auto vr = inputRight->get(i);
        if (vl < vr) {
            return -1;
        } else if (vl > vr) {
            return 1;
        }
    }
    return 0;
}

void TGChase::join(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::vector<size_t> &bodyNodeIdxs,
        std::pair<int, int> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<SegmentInserter> &output) {
    std::shared_ptr<const TGSegment> inputRight;
    if (bodyNodeIdxs.size() == 1) {
        size_t idbBodyAtomIdx = bodyNodeIdxs[0];
        auto bodyNode = nodes[idbBodyAtomIdx];
        inputRight = bodyNode.data;
    } else {
        auto ncols = nodes[bodyNodeIdxs[0]].data->getNColumns();
        std::vector<int> projectedPos;
        for(int i = 0; i < ncols; ++i)
            projectedPos.push_back(i);
        inputRight = mergeNodes(bodyNodeIdxs, projectedPos);
    }
    mergejoin(inputLeft,
            nodesLeft,
            inputRight,
            bodyNodeIdxs,
            joinVarPos,
            copyVarPosLeft,
            copyVarPosRight,
            output);
}

void TGChase::mergejoin(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::shared_ptr<const TGSegment> inputRight,
        const std::vector<size_t> &nodesRight,
        std::pair<int, int> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<SegmentInserter> &output) {
    std::chrono::system_clock::time_point startL = std::chrono::system_clock::now();

    //Sort the left segment by the join variable
    std::chrono::system_clock::time_point startS = std::chrono::system_clock::now();
    std::unique_ptr<TGSegmentItr> itrLeft;
    if (!inputLeft->isSortedBy(joinVarPos.first)) {
        if (nodesLeft.size() > 0) {
            SegmentCache &c = SegmentCache::getInstance();
            if (!c.contains(nodesLeft, joinVarPos.first)) {
                inputLeft = inputLeft->sortBy(joinVarPos.first);
                c.insert(nodesLeft, joinVarPos.first, inputLeft);
            } else {
                inputLeft = c.get(nodesLeft, joinVarPos.first);
            }
        } else {
            inputLeft = inputLeft->sortBy(joinVarPos.first);
        }
    }
    itrLeft = inputLeft->iterator();
    //Sort the right segment by the join variable
    std::unique_ptr<TGSegmentItr> itrRight;
    if (!inputRight->isSortedBy(joinVarPos.second)) {
        if (nodesRight.size() > 0) {
            SegmentCache &c = SegmentCache::getInstance();
            if (!c.contains(nodesRight, joinVarPos.second)) {
                inputRight = inputRight->sortBy(joinVarPos.second);
                c.insert(nodesRight, joinVarPos.second, inputRight);
            } else {
                inputRight = c.get(nodesRight, joinVarPos.second);
            }
        } else {
            inputRight = inputRight->sortBy(joinVarPos.second);
        }
    }
    itrRight = inputRight->iterator();
    durationMergeSort += std::chrono::system_clock::now() - startS;

    //Do the merge join
    if (itrLeft->hasNext()) {
        itrLeft->next();
    } else {
        return;
    }

    if (itrRight->hasNext()) {
        itrRight->next();
    } else {
        return;
    }

#if DEBUG
    size_t total = 0;
    size_t max = 65536;
#endif

    long countLeft = -1;
    size_t currentKey = 0;
    auto sizerow = copyVarPosLeft.size() + copyVarPosRight.size();
    Term_t currentrow[sizerow + 2];
    int res = cmp(itrLeft, itrRight, joinVarPos);
    while (true) {
        //Are they matching?
        while (res < 0 && itrLeft->hasNext()) {
            itrLeft->next();
            res = cmp(itrLeft, itrRight, joinVarPos);
        }

        if (res < 0) //The first iterator is finished
            break;

        while (res > 0 && itrRight->hasNext()) {
            itrRight->next();
            res = cmp(itrLeft, itrRight, joinVarPos);
        }

        if (res > 0) { //The second iterator is finished
            break;
        }

        if (res == 0) {
            if (countLeft == -1) {
                itrLeft->mark();
                countLeft = 1;
                currentKey = itrLeft->get(joinVarPos.first);
                while (itrLeft->hasNext()) {
                    itrLeft->next();
                    auto newKey = itrLeft->get(joinVarPos.first);
                    if (newKey != currentKey) {
                        break;
                    }
                    countLeft++;
                }
            }
#if DEBUG
            total += countLeft;
            while (total >= max) {
                LOG(TRACEL) << "Count = " << countLeft << ", total = " << total;
                max = max + max;
            }
#endif

            itrLeft->reset();
            //Move the left iterator countLeft times and emit tuples
            for(int idx = 0; idx < copyVarPosRight.size(); ++idx) {
                auto rightPos = copyVarPosRight[idx];
                currentrow[copyVarPosLeft.size() + idx] = itrRight->get(rightPos);
            }
            auto c = 0;
            bool leftActive = true;
            while (c < countLeft) {
                for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
                    auto leftPos = copyVarPosLeft[idx];
                    auto el = itrLeft->get(leftPos);
                    currentrow[idx] = el;
                }
                output->addRow(currentrow);
                if (itrLeft->hasNext()) {
                    itrLeft->next();
                } else {
                    leftActive = false;
                }
                c++;
            }

            //Move right
            if (!itrRight->hasNext()) {
                break;
            } else {
                itrRight->next();
            }
            //Does the right row have not the same value as before?
            auto newKey = itrRight->get(joinVarPos.second);
            if (newKey != currentKey) {
                countLeft = -1;
                //LeftItr is already pointing to the next element ...
                if (!leftActive) {
                    break;
                }
                res = cmp(itrLeft, itrRight, joinVarPos);
            }
        }
    }
#if DEBUG
    LOG(TRACEL) << "Total = " << total;
    std::chrono::duration<double> secL = std::chrono::system_clock::now() - startL;
    LOG(TRACEL) << "merge_join: time : " << secL.count() * 1000;
#endif
}

std::shared_ptr<const TGSegment> TGChase::retainVsNodeFast(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples) {
    std::unique_ptr<SegmentInserter> inserter;
    const uint8_t ncols = newtuples->getNColumns();
    Term_t row[ncols];

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
        int res = cmp(leftItr, rightItr);
        if (res < 0) {
            moveLeftItr = true;
        } else if (res > 0) {
            moveRightItr = true;
            if (isFiltered) {
                for(int i = 0; i < ncols; ++i) {
                    row[i] = rightItr->get(i);
                }
                inserter->addRow(row);
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
                inserter = std::unique_ptr<SegmentInserter>(
                        new SegmentInserter(rightItr->getNFields()));
                //Copy all the previous new tuples in the right iterator
                size_t i = 0;
                auto itrTmp = newtuples->iterator();
                while (i < (startCopyingIdx + countNew) && itrTmp->hasNext()) {
                    itrTmp->next();
                    if (i >= startCopyingIdx) {
                        for(int i = 0; i < ncols; ++i) {
                            row[i] = itrTmp->get(i);
                        }
                        inserter->addRow(row);
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
            inserter->addRow(row);
        }
        while (rightItr->hasNext()) {
            rightItr->next();
            for(int i = 0; i < ncols; ++i) {
                row[i] = rightItr->get(i);
            }
            inserter->addRow(row);
        }
        return fromSeg2TGSeg(inserter->getSegment(), 0, true, 0); //TODO
    } else {
        if (countNew > 0 || activeRightValue) {
            if (startCopyingIdx == 0) {
                //They are all new ...
                return newtuples;
            } else {
                //Remove the initial duplicates
                inserter = std::unique_ptr<SegmentInserter>(
                        new SegmentInserter(rightItr->getNFields()));
                //Copy all the previous new tuples in the right iterator
                size_t i = 0;
                auto itrTmp = newtuples->iterator();
                while (itrTmp->hasNext()) {
                    itrTmp->next();
                    if (i >= startCopyingIdx) {
                        for(int i = 0; i < ncols; ++i) {
                            row[i] = itrTmp->get(i);
                        }
                        inserter->addRow(row);
                    }
                    i++;
                }
                return fromSeg2TGSeg(inserter->getSegment(), 0, true, 0); //TODO
            }
        } else {
            //They are all duplicates
            return std::shared_ptr<const TGSegment>();
        }
    }
}

std::shared_ptr<const TGSegment> TGChase::retain(
        PredId_t p,
        std::shared_ptr<const TGSegment> newtuples) {
    if (!pred2Nodes.count(p)) {
        return newtuples;
    }
    auto &nodeIdxs = pred2Nodes[p];
    for(auto &nodeIdx : nodeIdxs) {
        auto node = nodes[nodeIdx];
        newtuples = retainVsNodeFast(node.data, newtuples);

        if (newtuples == NULL || newtuples->isEmpty()) {
            return std::shared_ptr<const TGSegment>();
        } else {
        }
    }
    return newtuples;
}

size_t TGChase::getNDerivedFacts() {
    size_t nderived = 0;
    for(auto &node : nodes) {
        nderived += node.data->getNRows();
    }
    return nderived;
}

size_t TGChase::getNnodes() {
    return nodes.size();
}

void TGChase::shouldSortDelDupls(const Literal &head,
        const std::vector<Literal> &bodyAtoms,
        const std::vector<std::vector<size_t>> &bodyNodes,
        bool &shouldSort,
        bool &shouldDelDupl) {
    if (bodyAtoms.size() == 1) {
        auto &bodyAtom = bodyAtoms[0];
        auto th = head.getTuple();
        auto tb = bodyAtom.getTuple();
        int sortedFields = 0;
        for(int i = 0; i < th.getSize() && i < tb.getSize(); ++i) {
            if (th.get(i).getId() != tb.get(i).getId()) {
                break;
            }
            sortedFields++;
        }
        shouldSort = !(sortedFields == th.getSize()) || (bodyNodes.size() > 0 && bodyNodes[0].size() > 1);
        shouldDelDupl = th.getSize() < tb.getSize();
    } else {
        shouldSort = shouldDelDupl = true;
    }
}

void TGChase::postprocessJoin(
        std::shared_ptr<const TGSegment> &intermediateResults,
        std::vector<std::shared_ptr<Column>> &intermediateResultsNodes,
        bool replace) {
    //Take out the last two columns
    /*auto ncolumns = intermediateResults->getNColumns();
      if (ncolumns < 2) {
      LOG(ERRORL) << "Cannot happen";
      throw 10;
      }
      std::vector<std::shared_ptr<Column>> columns;
      for(int i = 0; i < ncolumns - 2; ++i) {
      columns.push_back(intermediateResults->getColumn(i));
      }
      if (replace) {
    //Add one extra column
    CompressedColumnBlock b(0,1,columns[0]->size());
    std::vector<CompressedColumnBlock> blocks;
    blocks.push_back(b);
    columns.push_back(std::shared_ptr<Column>(
    new CompressedColumn(blocks, columns[0]->size())));
    }

    //Save the columns with the mappings to the nodes
    auto col1 = intermediateResults->getColumn(ncolumns - 2);
    auto col2 = intermediateResults->getColumn(ncolumns - 1);
    intermediateResultsNodes.push_back(col1);
    intermediateResultsNodes.push_back(col2);

    //Create new intermediate results
    intermediateResults = std::shared_ptr<const Segment>(new Segment(
    columns.size(), columns));*/
    //TODO
}

bool TGChase::executeRule(TGChase_SuperNode &node) {
    auto &bodyNodes = node.incomingEdges;
    Rule &rule = rules[node.ruleIdx];
#ifdef WEBINTERFACE
    currentRule = rule.tostring();
    currentPredicate = rule.getFirstHead().getPredicate().getId();
#endif

    //LOG(INFOL) << "Executing rule " << rule.tostring(program, &layer) <<
    //    " " << rule.getFirstHead().getPredicate().getId() << " " << node.ruleIdx;

    //Perform the joins and populate the head
    auto &bodyAtoms = rule.getBody();
    //Maybe Rearrange the body atoms? Don't forget to also re-arrange the body
    //nodes

    std::vector<size_t> varsIntermediate;
    std::shared_ptr<const TGSegment> intermediateResults;
    //The following data structure is used only if trackProvenance=true
    std::vector<std::shared_ptr<Column>> intermediateResultsNodes;
    size_t currentBodyNode = 0;
    bool firstBodyAtomIsIDB = false;

    for(size_t i = 0; i < bodyAtoms.size(); ++i) {
        const Literal &currentBodyAtom = bodyAtoms[i];
        VTuple currentVars = currentBodyAtom.getTuple();

        //Which are the positions (left,right) of the variable to join?
        std::pair<int, int> joinVarPos;
        //Which variables should be copied from the existing collection?
        std::vector<int> copyVarPosLeft;
        //Which variables should be copied from the new body atom?
        std::vector<int> copyVarPosRight;
        //Compute the value of the vars
        computeVarPos(varsIntermediate, i, bodyAtoms, rule.getFirstHead(),
                joinVarPos, copyVarPosLeft, copyVarPosRight);

        //Get the node of the current bodyAtom (if it's IDB)
        auto &bodyAtom = bodyAtoms[i];
        std::vector<size_t> newVarsIntermediateResults;

        if (i == 0) {
            //Process first body atom, no join
            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            if (bodyAtom.getPredicate().getType() == EDB) {
                intermediateResults = processFirstAtom_EDB(bodyAtom,
                        copyVarPosRight);
            } else {
                firstBodyAtomIsIDB = true;
                intermediateResults = mergeNodes(
                        bodyNodes[currentBodyNode], copyVarPosRight);
                currentBodyNode++;
            }
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            durationFirst += end - start;
        } else {
            const uint8_t extraColumns = trackProvenance ? 2 : 0;
            std::unique_ptr<SegmentInserter> newIntermediateResults
                (new SegmentInserter(copyVarPosLeft.size() +
                                     copyVarPosRight.size() + extraColumns));

            //Perform merge join between the intermediate results and
            //the new collection
            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            join(intermediateResults,
                    i == 1 && firstBodyAtomIsIDB ? bodyNodes[0] : noBodyNodes,
                    bodyNodes[currentBodyNode],
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    newIntermediateResults);
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            durationJoin += end - start;

            intermediateResults = fromSeg2TGSeg(newIntermediateResults->getSegment(), 0, false, 0); //TODO
            if (trackProvenance) {
                //Process the output of nodes
                postprocessJoin(intermediateResults, intermediateResultsNodes,
                        i < bodyAtoms.size() - 1);
            }

            //Update the list of variables from the left atom
            for(auto varIdx = 0; varIdx < copyVarPosLeft.size(); ++varIdx) {
                auto var = copyVarPosLeft[varIdx];
                newVarsIntermediateResults.push_back(varsIntermediate[var]);
            }
            currentBodyNode++;
        }
        //Update the list of variables from the right atom
        for(auto varIdx = 0; varIdx < copyVarPosRight.size(); ++varIdx) {
            auto varPos = copyVarPosRight[varIdx];
            newVarsIntermediateResults.push_back(
                    currentVars.get(varPos).getId());
        }

        varsIntermediate = newVarsIntermediateResults;
        if (intermediateResults->isEmpty()) {
            intermediateResults = std::shared_ptr<const TGSegment>();
            break;
        }
    }

    //Filter out the derivations produced by the rule
    auto nonempty = !(intermediateResults == NULL || intermediateResults->isEmpty());
    if (nonempty) {
        //Compute the head
        std::chrono::steady_clock::time_point start =
            std::chrono::steady_clock::now();
        Literal head = rule.getFirstHead();
        bool shouldSort = true, shouldDelDupl = true;
        shouldSortDelDupls(head, bodyAtoms, bodyNodes, shouldSort, shouldDelDupl);
        intermediateResults = projectHead(head,
                varsIntermediate, intermediateResults,
                shouldSort,
                shouldDelDupl);
        std::chrono::steady_clock::time_point end =
            std::chrono::steady_clock::now();
        auto dur = end - start;
        durationCreateHead += dur;

        start = std::chrono::steady_clock::now();
        auto retainedTuples = retain(currentPredicate, intermediateResults);
        end = std::chrono::steady_clock::now();
        dur = end - start;
        durationRetain += dur;

        nonempty = !(retainedTuples == NULL || retainedTuples->isEmpty());
        if (nonempty) {
            auto nodeId = nodes.size();
            nodes.emplace_back();
            TGChase_Node &outputNode = nodes.back();
            outputNode.ruleIdx = node.ruleIdx;
            outputNode.step = node.step;
            outputNode.data = retainedTuples;
            //Index the non-empty node
            pred2Nodes[currentPredicate].push_back(nodeId);
        }
    }
    return nonempty;
}

std::shared_ptr<const TGSegment> TGChase::projectHead(const Literal &head,
        std::vector<size_t> &vars,
        std::shared_ptr<const TGSegment> intermediateResults,
        bool shouldSort,
        bool shouldDelDupl) {
    //Project the columns to instantiate the head
    const auto &tupleHead = head.getTuple();
    assert(tupleHead.getSize() == vars.size());
    assert(intermediateResults->getNColumns() == vars.size());
    if (tupleHead.getSize() > 1) {
        if (tupleHead.getSize() == 2) {
            if (vars[0] != tupleHead.get(0).getId()) {
                //I have to swap the variables
                intermediateResults = intermediateResults->swap();
            }
        } else {
            LOG(ERRORL) << "Not implemented";
            throw 10;
        }
    }
    if (shouldSort) {
        if (shouldDelDupl) {
            return intermediateResults->sort()->unique();
        } else {
            return intermediateResults->sort();
        }
    } else {
        if (shouldDelDupl) {
            return intermediateResults->unique();
        } else {
            return intermediateResults;
        }
    }
}

size_t TGChase::getSizeTable(const PredId_t predid) const {
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

FCIterator TGChase::getTableItr(const PredId_t predid) {
    LOG(ERRORL) << "Method not implemented";
    throw 10;
}

FCTable *TGChase::getTable(const PredId_t predid) {
    //Create a FCTable obj with the nodes in predid
    uint8_t card = program->getPredicateCard(predid);
    VTuple tuple(card);
    Literal query(program->getPredicate(predid), tuple);
    FCTable *t = new FCTable(NULL, card);
    if (pred2Nodes.count(predid)) {
        auto &nodeIDs = pred2Nodes[predid];
        for (auto &nodeId : nodeIDs) {
            auto &node = nodes[nodeId];
            auto data = fromTGSeg2Seg(node.data);
            std::shared_ptr<const FCInternalTable> table =
                std::shared_ptr<const FCInternalTable>(
                        new InmemoryFCInternalTable(card, nodeId, true, data));
            FCBlock block(nodeId, table, query, 0, NULL, 0, true);
            t->addBlock(block);
        }
    }
    return t;
}

size_t TGChase::getCurrentIteration() {
    return currentIteration;
}

#ifdef WEBINTERFACE
std::string TGChase::getCurrentRule() {
    return currentRule;
}

PredId_t TGChase::getCurrentPredicate() {
    return currentPredicate;
}
#endif
