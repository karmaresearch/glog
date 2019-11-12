#include <vlog/tgchase.h>

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

std::shared_ptr<const Segment> TGChase::processFirstAtom_IDB(
        std::vector<size_t> &nodeIdxs,
        std::vector<int> &copyVarPos) {
    if (nodeIdxs.size() == 1) {
        size_t idbBodyAtomIdx = nodeIdxs[0];
        auto bodyNode = nodes[idbBodyAtomIdx];
        return processFirstAtom_IDB(bodyNode.data, copyVarPos, idbBodyAtomIdx);
    } else {
        const uint8_t extraCol = trackProvenance ? 1 : 0;
        SegmentInserter ins(copyVarPos.size() + extraCol);
        for(auto idbBodyAtomIdx : nodeIdxs) {
            auto bodyNode = nodes[idbBodyAtomIdx];
            for(int i = 0; i < copyVarPos.size(); ++i) {
                int pos = copyVarPos[i];
                auto col = bodyNode.data->getColumn(pos);
                ins.addColumn(i, col, false);
            }
            if (trackProvenance) {
                //Add an extra column with the node ID
                ins.addColumn(copyVarPos.size(), std::shared_ptr<Column>(
                            new CompressedColumn(idbBodyAtomIdx,
                                bodyNode.data->getNRows())), true);
            }
        }
        return ins.getSortedAndUniqueSegment();
    }
}

std::shared_ptr<const Segment> TGChase::processFirstAtom_IDB(
        std::shared_ptr<const Segment> &input,
        std::vector<int> &copyVarPos,
        size_t nodeId) {
    std::vector<std::shared_ptr<Column>> columns;
    for(int i = 0; i < copyVarPos.size(); ++i) {
        int pos = copyVarPos[i];
        auto col = input->getColumn(pos);
        columns.push_back(col);
    }
    if (trackProvenance) {
        auto nrows = columns[0]->size();
        columns.push_back(std::shared_ptr<Column>(
                    new CompressedColumn(nodeId, nrows)));
    }
    return std::shared_ptr<const Segment>(
            new Segment(columns.size(), columns));
}

std::shared_ptr<const Segment> TGChase::processFirstAtom_EDB(
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
    auto seg = table->getSegment();
    std::vector<std::shared_ptr<Column>> columns;
    for(int i = 0; i < copyVarPos.size(); ++i) {
        int pos = copyVarPos[i];
        auto col = seg->getColumn(pos);
        columns.push_back(col);
    }
    if (trackProvenance) {
        //Add a column with a fictional node
        auto nrows = columns[0]->size();
        columns.push_back(std::shared_ptr<Column>(new CompressedColumn(~0lu, nrows)));
    }
    return std::shared_ptr<const Segment>(
            new Segment(columns.size(), columns));
}

int TGChase::cmp(std::unique_ptr<SegmentIterator> &inputLeft,
        std::unique_ptr<SegmentIterator> &inputRight,
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

int TGChase::cmp(std::unique_ptr<SegmentIterator> &inputLeft,
        std::unique_ptr<SegmentIterator> &inputRight) {
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

std::shared_ptr<const Segment> TGChase::concatenate(
        std::vector<size_t> &bodyNodeIdxs) {
    std::unique_ptr<SegmentInserter> ins;
    int nfields = -1;
    for(size_t i = 0; i < bodyNodeIdxs.size(); ++i) {
        auto idbBodyAtomIdx = bodyNodeIdxs[i];
        auto &node = nodes[idbBodyAtomIdx].data;
        if (i == 0) {
            nfields = node->getNColumns();
            const uint8_t extraCol = trackProvenance ? 1 : 0;
            ins = std::unique_ptr<SegmentInserter>(
                    new SegmentInserter(nfields + extraCol));
        }
        for(int j = 0; j < nfields; ++j) {
            ins->addColumn(j, node->getColumn(j), false);
        }
        if (trackProvenance) {
            //Add an extra column with the node ID
            ins->addColumn(nfields, std::shared_ptr<Column>(
                        new CompressedColumn(idbBodyAtomIdx,
                            node->getNRows())), true);
        }
    }
    return ins->getSegment();
}

void TGChase::join(
        std::shared_ptr<const Segment> inputLeft,
        std::vector<size_t> &bodyNodeIdxs,
        std::pair<int, int> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<SegmentInserter> &output) {
    std::shared_ptr<const Segment> inputRight;
    if (bodyNodeIdxs.size() == 1) {
        size_t idbBodyAtomIdx = bodyNodeIdxs[0];
        auto bodyNode = nodes[idbBodyAtomIdx];
        inputRight = bodyNode.data;
        if (trackProvenance) {
            //Add extra column
            std::vector<std::shared_ptr<Column>> columns;
            for(int j = 0; j < inputRight->getNColumns(); ++j) {
                columns.push_back(inputRight->getColumn(j));
            }
            columns.push_back(std::shared_ptr<Column>(new CompressedColumn(
                            idbBodyAtomIdx, inputRight->getNRows())));
            inputRight = std::shared_ptr<const Segment>(new Segment(
                        columns.size(), columns));
        }
    } else {
        inputRight = concatenate(bodyNodeIdxs);
    }
    mergejoin(inputLeft, inputRight, joinVarPos, copyVarPosLeft, copyVarPosRight,
            output);
}

void TGChase::mergejoin(
        std::shared_ptr<const Segment> inputLeft,
        std::shared_ptr<const Segment> inputRight,
        std::pair<int, int> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<SegmentInserter> &output) {
    std::chrono::system_clock::time_point startL = std::chrono::system_clock::now();

    //Sort the left segment by the join variable
    std::chrono::system_clock::time_point startS = std::chrono::system_clock::now();
    auto sortedInputLeft = inputLeft->sortByField(joinVarPos.first);
    //Sort the right segment by the join variable
    auto sortedInputRight = inputRight->sortByField(joinVarPos.second);
    durationMergeSort += std::chrono::system_clock::now() - startS;
    auto itrLeft = sortedInputLeft->iterator();
    auto itrRight = sortedInputRight->iterator();

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

std::shared_ptr<const Segment> TGChase::retainVsNode(
        std::shared_ptr<const Segment> existuples,
        std::shared_ptr<const Segment> newtuples) {

    //Do outer join
    auto leftItr = existuples->iterator();
    auto rightItr = newtuples->iterator();
    bool activeRightValue = false;
    bool moveLeftItr = true;
    bool moveRightItr = true;
    std::unique_ptr<SegmentInserter> inserter = std::unique_ptr<SegmentInserter>(
            new SegmentInserter(rightItr->getNFields()));
    while (true) {
        if (moveLeftItr) {
            if (leftItr->hasNext()) {
                leftItr->next();
            } else {
                break;
            }
            moveLeftItr = false;
        }

        if (moveRightItr) {
            if (rightItr->hasNext()) {
                activeRightValue = true;
                rightItr->next();
            } else {
                activeRightValue = false;
                break;
            }
            moveRightItr = false;
        }

        //Compare the iterators
        int res = cmp(leftItr, rightItr);
        if (res < 0) {
            moveLeftItr = true;
        } else if (res > 0) {
            moveRightItr = true;
            inserter->addRow(*rightItr.get());
        } else {
            moveLeftItr = moveRightItr = true;
            activeRightValue = false;
        }
    }

    if (activeRightValue)
        inserter->addRow(*rightItr.get());
    while (rightItr->hasNext()) {
        rightItr->next();
        inserter->addRow(*rightItr.get());
    }
    return inserter->getSegment();
}

std::shared_ptr<const Segment> TGChase::retainVsNodeFast(
        std::shared_ptr<const Segment> existuples,
        std::shared_ptr<const Segment> newtuples) {
    std::unique_ptr<SegmentInserter> inserter;

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
                inserter->addRow(*rightItr.get());
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
                    if (i >= startCopyingIdx)
                        inserter->addRow(*itrTmp.get());
                    i++;
                }
                isFiltered = true;
            }
        }
    }

    if (isFiltered) {
        if (activeRightValue)
            inserter->addRow(*rightItr.get());
        while (rightItr->hasNext()) {
            rightItr->next();
            inserter->addRow(*rightItr.get());
        }
        return inserter->getSegment();
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
                    if (i >= startCopyingIdx)
                        inserter->addRow(*itrTmp.get());
                    i++;
                }
                return inserter->getSegment();
            }
        } else {
            //They are all duplicates
            return std::shared_ptr<const Segment>();
        }
    }
}

std::shared_ptr<const Segment> TGChase::retain(
        PredId_t p,
        std::shared_ptr<const Segment> newtuples) {
    if (!pred2Nodes.count(p)) {
        return newtuples;
    }
    auto &nodeIdxs = pred2Nodes[p];
    for(auto &nodeIdx : nodeIdxs) {
        auto node = nodes[nodeIdx];
        newtuples = retainVsNodeFast(node.data, newtuples);

        if (newtuples == NULL || newtuples->isEmpty()) {
            return std::shared_ptr<const Segment>();
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
        shouldSort = !(sortedFields == th.getSize());
        shouldDelDupl = th.getSize() < tb.getSize();
    } else {
        shouldSort = shouldDelDupl = true;
    }
}

void TGChase::postprocessJoin(
        std::shared_ptr<const Segment> &intermediateResults,
        std::vector<std::shared_ptr<Column>> &intermediateResultsNodes,
        bool replace) {
    //Take out the last two columns
    auto ncolumns = intermediateResults->getNColumns();
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
                columns.size(), columns));
}

bool TGChase::executeRule(TGChase_SuperNode &node) {
    auto &bodyNodes = node.incomingEdges;
    Rule &rule = rules[node.ruleIdx];
#ifdef WEBINTERFACE
    currentRule = rule.tostring();
    currentPredicate = rule.getFirstHead().getPredicate().getId();
#endif

    //    LOG(DEBUGL) << "Executing rule " << rule.tostring(program, &layer) <<
    //        " " << rule.getFirstHead().getPredicate().getId() << " " << node.ruleIdx;

    //Perform the joins and populate the head
    auto &bodyAtoms = rule.getBody();
    //Maybe Rearrange the body atoms? Don't forget to also re-arrange the body
    //nodes

    std::vector<size_t> varsIntermediate;
    std::shared_ptr<const Segment> intermediateResults;
    //Used only if trackProvenance=true
    std::vector<std::shared_ptr<Column>> intermediateResultsNodes;
    size_t currentBodyNode = 0;

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
                intermediateResults = processFirstAtom_IDB(
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
                    bodyNodes[currentBodyNode],
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    newIntermediateResults);
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            durationJoin += end - start;

            intermediateResults = newIntermediateResults->getSegment();
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
            intermediateResults = std::shared_ptr<const Segment>();
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
        shouldSortDelDupls(head, bodyAtoms, shouldSort, shouldDelDupl);
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

std::shared_ptr<const Segment> TGChase::projectHead(const Literal &head,
        std::vector<size_t> &vars,
        std::shared_ptr<const Segment> intermediateResults,
        bool shouldSort,
        bool shouldDelDupl) {
    //Project the columns to instantiate the head
    std::vector<std::shared_ptr<Column>> columns;
    const auto &tupleHead = head.getTuple();
    for(int i = 0; i < tupleHead.getSize(); ++i) {
        for(int j = 0; j < vars.size(); ++j) {
            if (vars[j] == tupleHead.get(i).getId()) {
                columns.push_back(intermediateResults->getColumn(j));
                break;
            }
        }
    }
    std::shared_ptr<const Segment> segment(new Segment(columns.size(),
                columns));

    if (shouldSort) {
        if (shouldDelDupl) {
            return SegmentInserter::unique(segment->sortBy(NULL));
        } else {
            return segment->sortBy(NULL);
        }
    } else {
        if (shouldDelDupl) {
            return SegmentInserter::unique(segment);
        } else {
            return segment;
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
            std::shared_ptr<const FCInternalTable> table =
                std::shared_ptr<const FCInternalTable>(
                        new InmemoryFCInternalTable(card, nodeId, true, node.data));
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