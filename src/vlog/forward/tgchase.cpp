#include <vlog/tgchase.h>

TGChase::TGChase(EDBLayer &layer, Program *program) : layer(layer),
    program(program),
    durationJoin(0),
    durationRetain(0),
    durationCreateHead(0),
    durationFirst(0)
{
}

Program *TGChase::getProgram() {
    return program;
}

EDBLayer &TGChase::getEDBLayer() {
    return layer;
}

void TGChase::recursiveCreateNode(
        const size_t step,
        const size_t ruleIdx,
        std::vector<std::vector<size_t>> &input,
        std::vector<size_t> &currentRow,
        const size_t columnIdx) {
    size_t rowIdx = 0;
    while (rowIdx < input[columnIdx].size()) {
        currentRow[columnIdx] = input[columnIdx][rowIdx];
        if (columnIdx < input.size() - 1) {
            recursiveCreateNode(step, ruleIdx, input, currentRow,
                    columnIdx + 1);
        } else {
            //Create a new node
            nodes.emplace_back();
            TGChase_Node &newnode = nodes.back();
            newnode.step = step;
            newnode.ruleIdx = ruleIdx;
            newnode.incomingEdges = currentRow;
        }
        rowIdx++;
    }
}

void TGChase::run() {
    initRun();
    size_t nnodes = 0;
    size_t step = 0;
    rules = program->getAllRules();

    do {
        step++;
        currentIteration = step;
        nnodes = nodes.size();
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
            if (empty)
                continue;

            if (rule.getNIDBPredicates() == 0) {
                //It's a rule with only EDB body atoms. Create a single node
                //I only execute these rules in the first step
                if (step == 1) {
                    nodes.emplace_back();
                    TGChase_Node &newnode = nodes.back();
                    newnode.step = step;
                    newnode.ruleIdx = ruleIdx;
                    newnode.incomingEdges = std::vector<size_t>();
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

                    //Create new nodes with all possible combinations
                    std::vector<size_t> bodyNodes;
                    bodyNodes.resize(acceptableNodes.size());
                    recursiveCreateNode(step, ruleIdx, acceptableNodes,
                            bodyNodes, 0);
                }
            }
        }

        //TODO: Prune nodes that lead only to duplicate derivations

        //Execute the rule associated to the node
        size_t idxNode = nnodes;
        while (idxNode < nodes.size()) {
            bool res = executeRule(idxNode);
            if (!res) {
                nodes.erase(nodes.begin() + idxNode);
            } else {
                idxNode++;
            }
        }
    } while (nnodes != nodes.size());

    LOG(INFOL) << "Time first (ms): " << durationFirst.count();
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
                addedVars.insert(rightVars[i]);
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
        std::shared_ptr<const Segment> &input,
        std::vector<int> &copyVarPos) {
    std::vector<std::shared_ptr<Column>> columns;
    for(int i = 0; i < copyVarPos.size(); ++i) {
        int pos = copyVarPos[i];
        auto col = input->getColumn(pos);
        columns.push_back(col);
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

void TGChase::mergejoin(
        std::shared_ptr<const Segment> inputLeft,
        std::shared_ptr<const Segment> inputRight,
        std::pair<int, int> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<SegmentInserter> &output) {
    std::chrono::system_clock::time_point startL = std::chrono::system_clock::now();

    //Sort the left segment by the join variable
    auto sortedInputLeft = inputLeft->sortByField(joinVarPos.first);
    //Sort the right segment by the join variable
    auto sortedInputRight = inputRight->sortByField(joinVarPos.second);
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
    Term_t currentrow[sizerow];
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
    std::unique_ptr<SegmentInserter> inserter;

    bool isFiltered = false;
    size_t startCopyingIdx = 0;

    //Do outer join
    auto leftItr = existuples->iterator();
    auto rightItr = newtuples->iterator();
    bool moveLeftItr = true;
    bool moveRightItr = true;
    size_t countNew = 0;
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
                rightItr->next();
            } else {
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
            if (isFiltered) {
                inserter->addRow(*rightItr.get());
            } else {
                countNew++;
            }
        } else {
            moveLeftItr = moveRightItr = true;
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

    if (!isFiltered && countNew > 0 && startCopyingIdx > 0) {
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
        isFiltered = true;
    }

    if (isFiltered) {
        if (!moveRightItr)
            inserter->addRow(*rightItr.get());
        while (rightItr->hasNext()) {
            rightItr->next();
            inserter->addRow(*rightItr.get());
        }
        return inserter->getSegment();
    } else {
        if (countNew == 0 && moveRightItr) {
            return std::shared_ptr<const Segment>();
        } else {
            //They are all new
            return newtuples;
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
        newtuples = retainVsNode(node.data, newtuples);

        if (newtuples == NULL || newtuples->isEmpty())
            return std::shared_ptr<const Segment>();
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

bool TGChase::executeRule(size_t nodeId) {
    auto &node = nodes[nodeId];
    auto &bodyNodes = node.incomingEdges;
    Rule &rule = rules[node.ruleIdx];
#ifdef WEBINTERFACE
    currentRule = rule.tostring();
    currentPredicate = rule.getFirstHead().getPredicate().getId();
#endif

    LOG(DEBUGL) << "Executing rule " << rule.tostring(program, &layer) << " nodeId=" << nodeId << " " << rule.getFirstHead().getPredicate().getId();


    //Perform the joins and populate the head
    auto &bodyAtoms = rule.getBody();
    //Maybe Rearrange the body atoms? Don't forget to also re-arrange the body nodes

    std::vector<size_t> varsIntermediate;
    std::shared_ptr<const Segment> intermediateResults;
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
                size_t idbBodyAtomIdx = bodyNodes[currentBodyNode];
                auto bodyNode = nodes[idbBodyAtomIdx];
                intermediateResults = processFirstAtom_IDB(bodyNode.data,
                        copyVarPosRight);
                currentBodyNode++;
            }
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            durationFirst += end - start;
        } else {
            size_t idbBodyAtomIdx = bodyNodes[currentBodyNode];
            auto bodyNode = nodes[idbBodyAtomIdx];
            std::unique_ptr<SegmentInserter> newIntermediateResults
                (new SegmentInserter(copyVarPosLeft.size() + copyVarPosRight.size()));

            //Perform merge join between the intermediate results and
            //the new collection
            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();
            mergejoin(intermediateResults,
                    bodyNode.data,
                    joinVarPos,
                    copyVarPosLeft,
                    copyVarPosRight,
                    newIntermediateResults);
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            durationJoin += end - start;

            intermediateResults = newIntermediateResults->getSegment();

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
            node.data = retainedTuples;
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
