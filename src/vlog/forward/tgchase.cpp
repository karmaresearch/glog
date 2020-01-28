#include <vlog/tgchase.h>
#include <vlog/tgsegmentcache.h>

TGChase::TGChase(EDBLayer &layer, Program *program, bool useCacheRetain) :
    layer(layer),
    program(program),
    durationMergeSort(0),
    durationJoin(0),
    durationRetain(0),
    durationCreateHead(0),
    durationFirst(0),
    trackProvenance(true),
    cacheRetainEnabled(useCacheRetain)
{
    if (! program->stratify(stratification, nStratificationClasses)) {
        LOG(ERRORL) << "Program could not be stratified";
        throw std::runtime_error("Program could not be stratified");
    }
    LOG(DEBUGL) << "nStratificationClasses = " << nStratificationClasses;
}

Program *TGChase::getProgram() {
    return program;
}

EDBLayer &TGChase::getEDBLayer() {
    return layer;
}

std::shared_ptr<const TGSegment> TGChase::fromSeg2TGSeg(
        std::shared_ptr<const Segment> seg,
        size_t nodeId, bool isSorted, uint8_t sortedField) {
    auto ncols = seg->getNColumns();
    if (ncols == 1) {
        const std::vector<Term_t> &orig = seg->getColumn(0)->getVectorRef();
        std::vector<Term_t> tuples(orig);
        return std::shared_ptr<const TGSegment>(
                new UnaryTGSegment(tuples, nodeId, isSorted, sortedField));
    } else if (ncols == 2) {
        auto &col1 = seg->getColumn(0)->getVectorRef();
        size_t nrows = col1.size();
        if (trackProvenance) {
            bool constantNodeVal = seg->getColumn(1)->isConstant();
            if (constantNodeVal) {
                std::vector<Term_t> tuples(col1);
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithConstProvTGSegment(
                            tuples, nodeId, isSorted, sortedField));
            } else {
                auto itrProv = seg->getColumn(1)->getReader();
                std::vector<std::pair<Term_t, Term_t>> out(nrows);
                for(size_t i = 0; i < nrows; ++i) {
                    out[i].first = col1[i];
                    if (!itrProv->hasNext()) {
                        throw 10;
                    }
                    out[i].second = itrProv->next();
                }
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithProvTGSegment(out, ~0ul, isSorted, sortedField));
            }
        } else {
            auto &col2 = seg->getColumn(1)->getVectorRef();
            std::vector<std::pair<Term_t, Term_t>> out(nrows);
            for(size_t i = 0; i < nrows; ++i) {
                out[i].first = col1[i];
                out[i].second = col2[i];
            }
            return std::shared_ptr<const TGSegment>(
                    new BinaryTGSegment(out, nodeId, isSorted, sortedField));
        }
    } else if (ncols == 3 && trackProvenance) {
        auto &col1 = seg->getColumn(0)->getVectorRef();
        auto &col2 = seg->getColumn(1)->getVectorRef();
        bool constantNodeVal = seg->getColumn(2)->isConstant();
        size_t nrows = col1.size();
        if (constantNodeVal) {
            std::vector<std::pair<Term_t, Term_t>> out(nrows);
            for(size_t i = 0; i < nrows; ++i) {
                out[i].first = col1[i];
                out[i].second = col2[i];
            }
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithConstProvTGSegment(
                        out, seg->getColumn(2)->first(), isSorted, sortedField));
        } else {
            std::vector<BinWithProv> out(nrows);
            auto itrNode = seg->getColumn(2)->getReader();
            for(size_t i = 0; i < nrows; ++i) {
                out[i].first = col1[i];
                out[i].second = col2[i];
                if (!itrNode->hasNext()) {
                    throw 10;
                }
                out[i].node = itrNode->next();
            }
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithProvTGSegment(
                        out, ~0ul, isSorted, sortedField));
        }
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

bool lowerStrat(Rule &rule, int currentStrat, std::vector<int> &stratification) {
    for (auto &bodyAtom : rule.getBody()) {
        Predicate pred = bodyAtom.getPredicate();
        if (pred.getType() != EDB && stratification[pred.getId()] >= currentStrat) {
            return false;
        }
    }
    return true;
}

void TGChase::run() {
    initRun();
    size_t nnodes = 0;
    size_t step = 0;
    rules = program->getAllRules();

    for (int currentStrat = 0; currentStrat < nStratificationClasses; currentStrat++) {
        LOG(INFOL) << "Strat: " << currentStrat;
        size_t saved_step = step;
        do {
            step++;
            LOG(INFOL) << "Step " << step;
            currentIteration = step;
            nnodes = nodes.size();

            std::vector<TGChase_SuperNode> newnodes;
            for (size_t ruleIdx = 0; ruleIdx < rules.size(); ++ruleIdx) {
                auto &rule = rules[ruleIdx];
                PredId_t id = rule.getFirstHead().getPredicate().getId();
                if (nStratificationClasses > 1 && stratification[id] != currentStrat) {
                    LOG(DEBUGL) << "Skipping rule, wrong strat: " << rule.tostring();
                    continue;
                }
                LOG(DEBUGL) << "Considering rule " << rule.tostring();
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
                if (empty) {
                    LOG(DEBUGL) << "Empty; continuing";
                    continue;
                }

                if (rule.getNIDBPredicates() == 0) {
                    //It's a rule with only EDB body atoms.
                    //Create a single node. I only execute these rules in the first iteration
                    //of the current strat (which should be the first strat, I think).
                    if (step == saved_step + 1) {
                        newnodes.emplace_back();
                        TGChase_SuperNode &newnode = newnodes.back();
                        newnode.ruleIdx = ruleIdx;
                        newnode.step = step;
                        newnode.incomingEdges = std::vector<std::vector<size_t>>();
                        LOG(DEBUGL) << "Pushing node for rule " << rule.tostring();
                    } else {
                        LOG(DEBUGL) << "Skipping EDB rule";
                    }
                } else if (currentStrat > 0 && lowerStrat(rule, currentStrat, stratification)) {
                    // All IDBs of the body are of a lower strat, so we only execute
                    // this rule in the first iteration of the current strat.
                    if (step == saved_step + 1) {
                        for(int j = 0; j < nodesForRule.size(); ++j) {
                            if (! nodesForRule[j].empty()) {
                                LOG(DEBUGL) << "Pushing node for lowerStrat rule " << rule.tostring();
                                newnodes.emplace_back();
                                TGChase_SuperNode &newnode = newnodes.back();
                                newnode.ruleIdx = ruleIdx;
                                newnode.step = step;
                                newnode.incomingEdges = std::vector<std::vector<size_t>>();
                                newnode.incomingEdges.push_back(nodesForRule[j]);
                            }
                        }
                    } else {
                        LOG(DEBUGL) << "Skipping lowerStrat rule";
                    }
                } else {
                    size_t prevstep = 0;
                    if (step > 1) prevstep = step - 1;
                    // Why not just prevstep = step - 1;??? step is at least 1 here. --Ceriel
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
                                    } else {
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

                        if (acceptableNodesEmpty) {
                            continue;
                        }

                        LOG(DEBUGL) << "Pushing node for rule " << rule.tostring();
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
    }

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
        std::vector<std::pair<int, int>> &joinVarPos,
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
                    joinVarPos.push_back(std::pair<int,int>(j, i));
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
        bool shouldSortAndUnique = ncols < nodes[nodeIdxs[0]].data->getNColumns() || copyVarPos[0] != 0;
        if (ncols == 1) { //If it has only one column ...
            if (trackProvenance) {
                std::vector<std::pair<Term_t,Term_t>> tuples;
                for(auto idbBodyAtomIdx : nodeIdxs) {
                    nodes[idbBodyAtomIdx].data->appendTo(copyVarPos[0], tuples);
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
                    nodes[idbBodyAtomIdx].data->appendTo(copyVarPos[0], tuples);
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
                    nodes[idbBodyAtomIdx].data->appendTo(copyVarPos[0], copyVarPos[1], tuples);
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
                    nodes[idbBodyAtomIdx].data->appendTo(copyVarPos[0], copyVarPos[1], tuples);
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
                return std::shared_ptr<const TGSegment>(new UnaryWithConstProvTGSegment(tuples, nodeId, true, 0));
            } else {
                return std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, nodeId, true, 0));
            }
        } else {
            return input; //Can be reused
        }
    } else if (copyVarPos.size() == 2) {
        if (copyVarPos[0] == 0 && copyVarPos[1] == 1) {
            return input; //Can be reused
        } else {
            //Swapped relation
            std::vector<std::pair<Term_t, Term_t>> tuples;
            input->appendTo(1, 0, tuples);
            if (trackProvenance) {
                return std::shared_ptr<const TGSegment>(new BinaryWithConstProvTGSegment(tuples, nodeId, false, 0));
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
    if (trackProvenance) {
        CompressedColumnBlock b(~0ul, 0, nrows);
        std::vector<CompressedColumnBlock> blocks;
        blocks.push_back(b);
        columns.push_back(std::shared_ptr<Column>(
                    new CompressedColumn(blocks, nrows)));
    }

    return std::shared_ptr<const TGSegment>(new TGSegmentLegacy(columns, nrows, true, 0, trackProvenance));
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
        std::unique_ptr<TGSegmentItr> &inputRight,
        std::vector<std::pair<int, int>> &joinVarPos) {
    for (int i = 0; i < joinVarPos.size(); i++) {
        const auto valLeft = inputLeft->get(joinVarPos[i].first);
        const auto valRight = inputRight->get(joinVarPos[i].second);
        if (valLeft < valRight)
            return -1;
        else if (valLeft > valRight)
            return 1;
    }
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

void TGChase::leftjoin(
        std::shared_ptr<const TGSegment> inputLeft,
        std::vector<size_t> &bodyNodeIdxs,
        std::vector<std::pair<int, int>> &joinVarPos,
        std::vector<int> &copyVarPosLeft,
        std::unique_ptr<SegmentInserter> &output) {
    std::shared_ptr<const TGSegment> inputRight;
    LOG(DEBUGL) << "bodyNodeIdxs.size() = " << bodyNodeIdxs.size();
    for (int i = 0; i < bodyNodeIdxs.size(); i++) {
        LOG(DEBUGL) << "bodyNodeIdxs[" << i << "] = " << bodyNodeIdxs[i];
    }
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
    LOG(DEBUGL) << "inputRight cols = " << inputRight->getNColumns() << ", size = " << inputRight->getNRows();

    std::vector<uint8_t> fields1;
    std::vector<uint8_t> fields2;

    for (uint32_t i = 0; i < joinVarPos.size(); ++i) {
        fields1.push_back(joinVarPos[i].first);
        fields2.push_back(joinVarPos[i].second);
    }
    auto sortedInputLeft = inputLeft->sortBy(fields1);
    auto sortedInputRight = inputRight->sortBy(fields2);
    auto itrLeft = sortedInputLeft->iterator();
    auto itrRight = sortedInputRight->iterator();

    bool leftActive = false;
    bool rightActive = false;
    if (itrLeft->hasNext()) {
        itrLeft->next();
        leftActive = true;
    } else {
        return;
    }
    if (itrRight->hasNext()) {
        itrRight->next();
        rightActive = true;
    }

    auto sizerow = copyVarPosLeft.size();
    Term_t currentrow[sizerow + 2];

    while (leftActive && rightActive) {
        int res = cmp(itrLeft, itrRight, joinVarPos);
        if (res <= 0) {
            if (res < 0) {
                for(int idx = 0; idx < copyVarPosLeft.size(); ++idx) {
                    auto leftPos = copyVarPosLeft[idx];
                    auto el = itrLeft->get(leftPos);
                    currentrow[idx] = el;
                }
                output->addRow(currentrow);
            }
            if (itrLeft->hasNext()) {
                itrLeft->next();
            } else {
                leftActive = false;
            }
        } else {
            if (itrRight->hasNext()) {
                itrRight->next();
            } else {
                rightActive = false;
            }
        }
    }
    while (leftActive) {
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
    }
}

void TGChase::join(
        std::shared_ptr<const TGSegment> inputLeft,
        const std::vector<size_t> &nodesLeft,
        std::vector<size_t> &bodyNodeIdxs,
        std::vector<std::pair<int, int>> &joinVarPos,
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
        std::vector<std::pair<int, int>> &joinVarsPos,
        std::vector<int> &copyVarPosLeft,
        std::vector<int> &copyVarPosRight,
        std::unique_ptr<SegmentInserter> &output) {
    std::chrono::system_clock::time_point startL = std::chrono::system_clock::now();

    std::vector<uint8_t> fields1;
    std::vector<uint8_t> fields2;
    for (uint32_t i = 0; i < joinVarsPos.size(); ++i) {
        fields1.push_back(joinVarsPos[i].first);
        fields2.push_back(joinVarsPos[i].second);
    }
    std::chrono::system_clock::time_point startS = std::chrono::system_clock::now();

    //Sort the left segment by the join variable
    if (!inputLeft->isSortedBy(fields1)) {
        if (nodesLeft.size() > 0) {
            SegmentCache &c = SegmentCache::getInstance();
            if (!c.contains(nodesLeft, fields1)) {
                inputLeft = inputLeft->sortBy(fields1);
                c.insert(nodesLeft, fields1, inputLeft);
            } else {
                inputLeft = c.get(nodesLeft, fields1);
            }
        } else {
            inputLeft = inputLeft->sortBy(fields1);
        }
    }
    std::unique_ptr<TGSegmentItr> itrLeft = inputLeft->iterator();
    //Sort the right segment by the join variable
    if (!inputRight->isSortedBy(fields2)) {
        if (nodesRight.size() > 0) {
            SegmentCache &c = SegmentCache::getInstance();
            if (!c.contains(nodesRight, fields2)) {
                inputRight = inputRight->sortBy(fields2);
                c.insert(nodesRight, fields2, inputRight);
            } else {
                inputRight = c.get(nodesRight, fields2);
            }
        } else {
            inputRight = inputRight->sortBy(fields2);
        }
    }
    std::unique_ptr<TGSegmentItr> itrRight = inputRight->iterator();
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
    std::vector<Term_t> currentKey;
    auto sizerow = copyVarPosLeft.size() + copyVarPosRight.size();
    Term_t currentrow[sizerow + 2];
    for(int i = 0; i < sizerow + 2; ++i) currentrow[i] = 0;
    int res = cmp(itrLeft, itrRight, joinVarsPos);
    while (true) {
        //Are they matching?
        while (res < 0 && itrLeft->hasNext()) {
            itrLeft->next();
            res = cmp(itrLeft, itrRight, joinVarsPos);
        }

        if (res < 0) //The first iterator is finished
            break;

        while (res > 0 && itrRight->hasNext()) {
            itrRight->next();
            res = cmp(itrLeft, itrRight, joinVarsPos);
        }

        if (res > 0) { //The second iterator is finished
            break;
        }

        if (res == 0) {
            if (countLeft == -1) {
                currentKey.clear();
                itrLeft->mark();
                countLeft = 1;
                for (int i = 0; i < fields1.size(); i++) {
                    currentKey.push_back(itrLeft->get(fields1[i]));
                }
                while (itrLeft->hasNext()) {
                    itrLeft->next();
                    bool equal = true;
                    for (int i = 0; i < fields1.size(); i++) {
                        auto k = itrLeft->get(fields1[i]);
                        if (k != currentKey[i]) {
                            equal = false;
                            break;
                        }
                    }
                    if (! equal) {
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
                if (trackProvenance) {
                    currentrow[sizerow] = itrLeft->getNodeId();
                    currentrow[sizerow + 1] = itrRight->getNodeId();
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
            bool equal = true;
            for (int i = 0; i < fields2.size(); i++) {
                auto newKey = itrRight->get(fields2[i]);
                if (newKey != currentKey[i]) {
                    equal = false;
                    break;
                }
            }
            if (! equal) {
                countLeft = -1;
                //LeftItr is already pointing to the next element ...
                if (!leftActive) {
                    break;
                }
                res = cmp(itrLeft, itrRight, joinVarsPos);
            }
        }
    }
#if DEBUG
    LOG(TRACEL) << "Total = " << total;
    std::chrono::duration<double> secL = std::chrono::system_clock::now() - startS;
    LOG(TRACEL) << "merge_join: time : " << secL.count() * 1000;
#endif
}

std::shared_ptr<const TGSegment> TGChase::retainVsNodeFast(
        std::shared_ptr<const TGSegment> existuples,
        std::shared_ptr<const TGSegment> newtuples) {
    std::unique_ptr<SegmentInserter> inserter;
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
                        new SegmentInserter(ncols + extracol));
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
                        new SegmentInserter(ncols + extracol));
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
                    nodes[nodeIdxs[i]].data->appendTo(0, tuples);
                }
                std::sort(tuples.begin(), tuples.end());
                auto seg = std::shared_ptr<const TGSegment>(new UnaryTGSegment(tuples, ~0ul, true, 0));
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
                    nodes[nodeIdxs[i]].data->appendTo(0, 1, tuples);
                }
                std::sort(tuples.begin(), tuples.end());
                auto seg = std::shared_ptr<const TGSegment>(new BinaryTGSegment(tuples, ~0ul, true, 0));
                CacheRetainEntry entry;
                entry.nnodes = nodeIdxs.size();
                entry.seg = seg;
                cacheRetain[p] = entry;
            } else {
                LOG(ERRORL) << "Not supported";
                throw 10;
            }
        }
        auto existingTuples = cacheRetain[p].seg;
        newtuples = retainVsNodeFast(existingTuples, newtuples);
        if (newtuples == NULL || newtuples->isEmpty()) {
            return std::shared_ptr<const TGSegment>();
        }
    } else {
        for(auto &nodeIdx : nodeIdxs) {
            auto node = nodes[nodeIdx];
            newtuples = retainVsNodeFast(node.data, newtuples);
            if (newtuples == NULL || newtuples->isEmpty()) {
                return std::shared_ptr<const TGSegment>();
            }
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

std::shared_ptr<const Segment> TGChase::postprocessJoin(
        std::shared_ptr<const Segment> &intermediateResults,
        std::vector<std::shared_ptr<Column>> &intermediateResultsNodes) {
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
    //For now, always add one extra column
    CompressedColumnBlock b(0, 1, columns[0]->size());
    std::vector<CompressedColumnBlock> blocks;
    blocks.push_back(b);
    columns.push_back(std::shared_ptr<Column>(
                new CompressedColumn(blocks, columns[0]->size())));

    //Save the columns with the mappings to the nodes
    auto col1 = intermediateResults->getColumn(ncolumns - 2);
    auto col2 = intermediateResults->getColumn(ncolumns - 1);
    intermediateResultsNodes.push_back(col1);
    intermediateResultsNodes.push_back(col2);

    //Create new intermediate results
    return std::shared_ptr<const Segment>(new Segment(columns.size(), columns));
}

bool TGChase::executeRule(TGChase_SuperNode &node) {
    auto &bodyNodes = node.incomingEdges;
    Rule &rule = rules[node.ruleIdx];
#ifdef WEBINTERFACE
    currentRule = rule.tostring();
#endif
    currentPredicate = rule.getFirstHead().getPredicate().getId();

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
        std::vector<std::pair<int, int>> joinVarPos;
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

            //If empty then stop
            if (intermediateResults->isEmpty()) {
                intermediateResults = std::shared_ptr<const TGSegment>();
                break;
            }
        } else {
            const uint8_t extraColumns = trackProvenance ? 2 : 0;
            std::unique_ptr<SegmentInserter> newIntermediateResults
                (new SegmentInserter(copyVarPosLeft.size() +
                                     copyVarPosRight.size() + extraColumns));

            std::chrono::steady_clock::time_point start =
                std::chrono::steady_clock::now();

            if (currentBodyAtom.isNegated()) {
                // Negated atoms should not introduce new variables.
                assert(copyVarPosRight.size() == 0);
                leftjoin(intermediateResults,
                        bodyNodes[currentBodyNode],
                        joinVarPos,
                        copyVarPosLeft,
                        newIntermediateResults);
            } else {
                //Perform a join between the intermediate results and
                //the new collection
                if (joinVarPos.size() == 0) {
                    LOG(ERRORL) << "Join variables required";
                    throw 10;
                }
                join(intermediateResults,
                        i == 1 && firstBodyAtomIsIDB ? bodyNodes[0] : noBodyNodes,
                        bodyNodes[currentBodyNode],
                        joinVarPos,
                        copyVarPosLeft,
                        copyVarPosRight,
                        newIntermediateResults);
            }
            std::chrono::steady_clock::time_point end =
                std::chrono::steady_clock::now();
            durationJoin += end - start;

            //If empty then stop
            if (newIntermediateResults->isEmpty()) {
                intermediateResults = std::shared_ptr<const TGSegment>();
                break;
            }

            std::shared_ptr<const Segment> seg = newIntermediateResults->getSegment();
            if (trackProvenance) {
                //Process the output of nodes
                seg = postprocessJoin(seg, intermediateResultsNodes);
            }
            intermediateResults = fromSeg2TGSeg(seg , ~0ul, false, 0);
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
            if (trackProvenance) {
                createNewNodesWithProv(node.ruleIdx, node.step,
                        retainedTuples, intermediateResultsNodes);
            } else {
                //Add a single node
                auto nodeId = nodes.size();
                nodes.emplace_back();
                TGChase_Node &outputNode = nodes.back();
                outputNode.ruleIdx = node.ruleIdx;
                outputNode.step = node.step;
                outputNode.data = retainedTuples;
                pred2Nodes[currentPredicate].push_back(nodeId);
            }
        }
    }
    return nonempty;
}

void TGChase::createNewNodesWithProv(size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> seg,
        std::vector<std::shared_ptr<Column>> &provenance) {
    if (provenance.size() == 0) {
        //There was no join. Must replace the nodeID with a new one
        //Note at this point calling slice will create a new vector
        auto nodeId = nodes.size();
        nodes.emplace_back();
        TGChase_Node &outputNode = nodes.back();
        outputNode.ruleIdx = ruleIdx;
        outputNode.step = step;
        outputNode.data = seg->slice(nodeId, 0, seg->getNRows());
        pred2Nodes[currentPredicate].push_back(nodeId);
    } else {
        const auto nnodes = (provenance.size() + 2) / 2;
        const auto nrows = seg->getNRows();
        std::vector<size_t> provnodes(nrows * nnodes);
        for(size_t i = 0; i < nrows; ++i) {
            size_t provRowIdx = i;
            for(int j = nnodes - 1; j >= 0; j--) {
                if (j == 0) {
                    provnodes[i * nnodes] = provenance[0]->getValue(provRowIdx);
                } else {
                    provnodes[i * nnodes + j] = provenance[(j - 1)*2 + 1]->getValue(provRowIdx);
                    if (j > 1)
                        provRowIdx = provenance[(j-1)*2]->getValue(provRowIdx);
                }
            }
        }
        //For each tuple, now I know the sequence of nodes that derived them.
        //I re-sort the nodes depending on the sequence of nodes
        std::vector<size_t> providxs;
        auto resortedSeg = seg->sortByProv(nnodes, providxs, provnodes);
        size_t startidx = 0;
        std::vector<size_t> currentNodeList(nnodes);
        for(size_t i = 0; i < nrows; ++i) {
            bool hasChanged = i == 0;
            for(size_t j = 0; j < nnodes && !hasChanged; ++j) {
                const size_t m = i * nnodes + j;
                if (currentNodeList[j] != provnodes[m]) {
                    hasChanged = true;
                }
            }
            if (hasChanged) {
                if (startidx < i) {
                    //Create a new node
                    auto nodeId = nodes.size();
                    nodes.emplace_back();
                    TGChase_Node &outputNode = nodes.back();
                    outputNode.ruleIdx = ruleIdx;
                    outputNode.step = step;
                    outputNode.data = resortedSeg->slice(nodeId, startidx, i);
                    pred2Nodes[currentPredicate].push_back(nodeId);
                }
                startidx = i;
                for(size_t j = 0; j < nnodes; ++j) {
                    const size_t m = i * nnodes + j;
                    currentNodeList[j] = provnodes[m];
                }
            }
        }
        //Copy the last segment
        if (startidx < nrows) {
            auto nodeId = nodes.size();
            nodes.emplace_back();
            TGChase_Node &outputNode = nodes.back();
            outputNode.ruleIdx = ruleIdx;
            outputNode.step = step;
            outputNode.data = resortedSeg->slice(nodeId, startidx, nrows);
            pred2Nodes[currentPredicate].push_back(nodeId);
        }
    }
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
