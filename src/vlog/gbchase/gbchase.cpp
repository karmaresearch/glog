#include <vlog/gbchase/gbchase.h>
#include <vlog/gbchase/gbsegmentcache.h>

GBChase::GBChase(EDBLayer &layer, Program *program, bool useCacheRetain,
        bool trackProvenance,
        bool filterQueryCont) :
    layer(layer),
    program(program),
    trackProvenance(trackProvenance),
    filterQueryCont(filterQueryCont),
    g(trackProvenance, useCacheRetain, filterQueryCont),
    executor(trackProvenance, g, layer, program) {
        if (!program->stratify(stratification, nStratificationClasses)) {
            LOG(ERRORL) << "Program could not be stratified";
            throw std::runtime_error("Program could not be stratified");
        }
        LOG(DEBUGL) << "nStratificationClasses = " << nStratificationClasses;
        rules = program->getAllRules();
        if (trackProvenance) {
            g.setRulesProgramLayer(rules.data(), program, &layer);
        }
    }

Program *GBChase::getProgram() {
    return program;
}

EDBLayer &GBChase::getEDBLayer() {
    return layer;
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

bool lowerStrat(const Rule &rule, int currentStrat,
        const std::vector<int> &stratification) {
    for (auto &bodyAtom : rule.getBody()) {
        Predicate pred = bodyAtom.getPredicate();
        if (pred.getType() != EDB && stratification[pred.getId()] >= currentStrat) {
            return false;
        }
    }
    return true;
}

void GBChase::prepareRuleExecutionPlans(
        const size_t &ruleIdx,
        const size_t prevstep,
        const size_t step,
        std::vector<GBRuleInput> &newnodes) {
    //The first parameter in the pair records whether the associated
    //atom is negated. In this case, we do not apply semi-naive
    //evaluation
    const auto &rule = rules[ruleIdx];
    std::vector<std::pair<bool,std::vector<size_t>>> nodesForRule;
    for (auto &bodyAtom : rule.getBody()) {
        Predicate pred = bodyAtom.getPredicate();
        if (pred.getType() != EDB) {
            bool negated = bodyAtom.isNegated();
            nodesForRule.push_back(std::make_pair(negated,
                        g.getNodeIDsWithPredicate(pred.getId())));
        }
    }

    if (nodesForRule.size() == 0) {
        //The rule contains only EDB predicates
        assert(rule.getNIDBPredicates() == 0);
        newnodes.emplace_back();
        GBRuleInput &newnode = newnodes.back();
        newnode.ruleIdx = ruleIdx;
        newnode.step = step;
        newnode.incomingEdges = std::vector<std::vector<size_t>>();
    } else {
        for(int pivot = 0; pivot < nodesForRule.size(); ++pivot) {
            if (nodesForRule[pivot].first) { //Is negated, skip
                continue;
            }
            //First consider only combinations where at least one node
            //is derived in the previous step
            std::vector<std::vector<size_t>> acceptableNodes;
            bool acceptableNodesEmpty = false;
            for(int j = 0; j < nodesForRule.size(); ++j) {
                std::vector<size_t> selection;
                bool isNegated = nodesForRule[j].first;
                auto &nodes = nodesForRule[j].second;
                if (j < pivot && !isNegated) {
                    for(auto &nodeId : nodes) {
                        auto nodeStep = g.getNodeStep(nodeId);
                        if (nodeStep < prevstep) {
                            selection.push_back(nodeId);
                        } else {
                            break;
                        }
                    }
                } else if (j == pivot && !isNegated) {
                    for(auto &nodeId : nodes) {
                        auto nodeStep = g.getNodeStep(nodeId);
                        if (nodeStep >= prevstep) {
                            selection.push_back(nodeId);
                        }
                    }
                } else {
                    selection = nodes;
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

            if (filterQueryCont) {
                size_t nCombinations = 1;
                for(size_t i = 0; i < acceptableNodes.size(); ++i) {
                    nCombinations = nCombinations * acceptableNodes[i].size();
                }
                if (acceptableNodes.size() > 1 && nCombinations > 10) {
                    //The combinations are too many to test...
                    newnodes.emplace_back();
                    GBRuleInput &newnode = newnodes.back();
                    newnode.ruleIdx = ruleIdx;
                    newnode.step = step;
                    newnode.incomingEdges = acceptableNodes;
                } else {
                    size_t filteredCombinations = 0;
                    size_t redundantFreeCombinations = 0;
                    std::vector<int> processedCombs(acceptableNodes.size());
                    std::vector<size_t> currentComb(acceptableNodes.size());
                    for(size_t i = 0; i < acceptableNodes.size(); ++i) {
                        processedCombs[i] = -1;
                    }
                    size_t currentIdx = 0;

                    //Remember temporary nodes
                    std::map<size_t, std::vector<size_t>> replacements;

                    //Check all combinations
                    while (true) {
                        if (processedCombs[currentIdx] ==
                                acceptableNodes[currentIdx].size() - 1) {
                            //Go back one level
                            if (currentIdx == 0) {
                                break;
                            } else {
                                processedCombs[currentIdx] = -1;
                                currentIdx--;
                            }
                        } else {
                            processedCombs[currentIdx]++;
                            auto idx = processedCombs[currentIdx];
                            currentComb[currentIdx] = acceptableNodes[currentIdx][idx];
                            //Am I at the last? Then check
                            if (currentIdx == acceptableNodes.size() - 1) {
                                //Do the check!
                                bool rFree = false;
                                if (g.isRedundant(ruleIdx, currentComb,
                                            rFree)) {
                                    filteredCombinations++;
                                } else {
                                    //Take the tmp nodes into account
                                    for(int i = 0; i < currentComb.size(); ++i) {
                                        auto nodeId = currentComb[i];
                                        if (g.isTmpNode(nodeId)) {
                                            auto originalNode = acceptableNodes[
                                                i][processedCombs[i]];
                                            if (!replacements.count(originalNode)) {
                                                replacements.insert(
                                                        std::make_pair(originalNode,
                                                            std::vector<size_t>()));
                                            }
                                            replacements[originalNode].push_back(
                                                    nodeId);
                                            currentComb[i] = originalNode;
                                        }
                                    }
                                }
                                if (rFree) {
                                    redundantFreeCombinations++;
                                }
                            } else {
                                currentIdx++;
                            }
                        }
                    }

                    if (filteredCombinations == nCombinations) {
                        //All combinations are redundant, skip
                    } else {
                        if (filteredCombinations != 0) {
                            //LOG(WARNL) << "FIXME: (gbchase) Not (yet) implemented"
                            //    << filteredCombinations << " " << nCombinations <<
                            //    " rule " << ruleIdx;
                        }
                        //Use the temporary nodes
                        if (!replacements.empty()) {
                            std::map<size_t, size_t> finalReplacementMap;
                            for(auto &p : replacements) {
                                assert(p.second.size() > 0);
                                if (p.second.size() == 1) {
                                    finalReplacementMap.insert(std::make_pair(
                                                p.first, p.second[0]));
                                } else {
                                    size_t idx = 0;
                                    size_t card = ~0ul;
                                    for(size_t i = 0; i < p.second.size(); ++i) {
                                        auto nodeId = p.second[i];
                                        if (g.getNodeSize(nodeId) < card) {
                                            idx = i;
                                            card = g.getNodeSize(nodeId);
                                        }
                                    }
                                    finalReplacementMap.insert(std::make_pair(
                                                p.first, p.second[idx]));
                                }
                            }
                            //Do the replacement
                            for(size_t i = 0; i < acceptableNodes.size(); ++i) {
                                for(size_t j = 0; j < acceptableNodes[i].size();
                                        ++j) {
                                    auto n = acceptableNodes[i][j];
                                    if (finalReplacementMap.count(n)) {
                                        acceptableNodes[i][j] =
                                            finalReplacementMap[n];
                                    }
                                }
                            }
                        }

                        newnodes.emplace_back();
                        GBRuleInput &newnode = newnodes.back();
                        newnode.ruleIdx = ruleIdx;
                        newnode.step = step;
                        newnode.incomingEdges = acceptableNodes;
                        if (redundantFreeCombinations == nCombinations) {
                            newnode.retainFree = true;
                        } else {
                            newnode.retainFree = false;
                        }
                    }
                }
            } else {
                newnodes.emplace_back();
                GBRuleInput &newnode = newnodes.back();
                newnode.ruleIdx = ruleIdx;
                newnode.step = step;
                newnode.incomingEdges = acceptableNodes;
            }
        }
    }
}

std::pair<bool, size_t> GBChase::determineAdmissibleRule(
        const size_t &ruleIdx,
        const size_t stratumLevel,
        const size_t stepStratum,
        const size_t step) const {
    auto &rule = rules[ruleIdx];
    bool empty = false;
    for (auto &bodyAtom : rule.getBody()) {
        Predicate pred = bodyAtom.getPredicate();
        if (pred.getType() != EDB) {
            bool negated = bodyAtom.isNegated();
            if (!g.areNodesWithPredicate(pred.getId()) && !negated) {
                empty = true;
                break;
            }
        }
    }
    if (empty) {
        return std::make_pair(false, 0);
    }

    size_t prevstep = step - 1;
    if (rule.getNIDBPredicates() == 0) {
        //It's a rule with only EDB body atoms.
        //I only execute these rules in the first iteration
        //of the current strat (which should be the first strat,
        //I think).
        if (step != stepStratum + 1 || stratumLevel != 0) {
            return std::make_pair(false, 0);
        }
        prevstep = 0;
    } else if (stratumLevel > 0 &&
            lowerStrat(rule, stratumLevel, stratification)) {
        // All IDBs of the body are of a lower strat, so we only execute
        // this rule in the first iteration of the current strat.
        if (step != stepStratum + 1) {
            return std::make_pair(false, 0);
        }
        prevstep = 0;
    }
    return std::make_pair(true, prevstep);
}

void GBChase::determineAdmissibleRules(
        const std::vector<size_t> &ruleIdxs,
        const size_t stratumLevel,
        const size_t stepStratum,
        const size_t step,
        std::vector<std::pair<size_t,size_t>> &admissibleRules) const {
    for (size_t ruleIdx : ruleIdxs) {
        auto response = determineAdmissibleRule(ruleIdx, stratumLevel,
                stepStratum, step);
        if (response.first) {
            size_t prevstep = response.second;
            admissibleRules.push_back(std::make_pair(ruleIdx, prevstep));
        }
    }
}

size_t GBChase::executeRulesInStratum(
        const std::vector<size_t> &ruleIdxs,
        const size_t stratumLevel,
        const size_t stepStratum,
        size_t &step) {

    //Identify the rules that can be executed
    std::vector<std::pair<size_t,size_t>> admissibleRules;
    determineAdmissibleRules(
            ruleIdxs,
            stratumLevel,
            stepStratum,
            step,
            admissibleRules);

    //Create different execution plans using semi-naive evaluation
    std::vector<GBRuleInput> newnodes;
    for(auto p : admissibleRules) {
        auto ruleIdx = p.first;
        size_t prevstep = p.second;
        prepareRuleExecutionPlans(ruleIdx, prevstep, step, newnodes);
    }

    //Execute the rule associated to the node
    auto nnodes = g.getNNodes();
    auto nodesToProcess = newnodes.size();
    LOG(INFOL) << "Nodes to process " << nodesToProcess;
    for(size_t idxNode = 0; idxNode < nodesToProcess; ++idxNode) {
        executeRule(newnodes[idxNode]);
    }

    if (nnodes != g.getNNodes()) {
        auto derivedTuples = 0;
        for(size_t idx = nnodes; idx < g.getNNodes(); ++idx) {
            auto nrows = g.getNodeSize(idx);
            derivedTuples += nrows;
        }
        LOG(INFOL) << "Derived Tuples: " << derivedTuples;
        return derivedTuples;
    } else {
        return 0;
    }
}

void GBChase::run() {
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    initRun();
    size_t nnodes = 0;
    size_t step = 0;

    for (int currentStrat = 0;
            currentStrat < nStratificationClasses;
            currentStrat++) {
        LOG(INFOL) << "Strat: " << currentStrat;
        size_t stepStratum = step;
        std::vector<size_t> rulesInStratum;
        for (size_t ruleIdx = 0; ruleIdx < rules.size(); ++ruleIdx) {
            auto &rule = rules[ruleIdx];
            PredId_t id = rule.getFirstHead().getPredicate().getId();
            if (nStratificationClasses > 1 &&
                    stratification[id] != currentStrat) {
                continue;
            }
            rulesInStratum.push_back(ruleIdx);
        }

        do {
            step++;
            LOG(INFOL) << "Step " << step;
            currentIteration = step;
            g.cleanTmpNodes();
            nnodes = g.getNNodes();
            size_t derivedTuples = executeRulesInStratum(rulesInStratum,
                    currentStrat, stepStratum, step);
            //g.printStats();
        } while (g.getNNodes() != nnodes);
        g.cleanTmpNodes();
    }

    std::chrono::duration<double, std::milli> dur =
        std::chrono::system_clock::now() - start;
    LOG(INFOL) << "Runtime chase: " << dur.count();

    SegmentCache::getInstance().clear();
    executor.printStats();
    g.printStats();
    stopRun();
}

size_t GBChase::getNDerivedFacts() {
    return g.getNDerivedFacts();
}

size_t GBChase::getNnodes() {
    return g.getNNodes();
}

bool GBChase::executeRule(GBRuleInput &node, bool cleanDuplicates) {
    auto &bodyNodes = node.incomingEdges;
    Rule &rule = rules[node.ruleIdx];
#ifdef WEBINTERFACE
    currentRule = rule.tostring();
#endif

#ifdef DEBUG
    LOG(INFOL) << "Executing rule " << node.ruleIdx << " "
        << rule.tostring(program, &layer);
#endif

    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();

    size_t nders = 0;
    size_t nders_un = 0;

    auto &heads = rule.getHeads();
    int headIdx = 0;
    currentPredicate = heads[headIdx].getPredicate().getId();
    auto outputsRule = executor.executeRule(rule, node);
    bool nonempty = false;

    std::chrono::system_clock::time_point starth =
        std::chrono::system_clock::now();
    for (GBRuleOutput &outputRule : outputsRule) {
        currentPredicate = heads[headIdx++].getPredicate().getId();
        auto derivations = outputRule.segment;
        nders_un += derivations->getNRows();
        auto derivationNodes = outputRule.nodes;
        const bool shouldCleanDuplicates = cleanDuplicates &&
            !outputRule.uniqueTuples;
        nonempty |= !(derivations == NULL || derivations->isEmpty());
        if (nonempty) {
            //Keep only the new derivations
            std::shared_ptr<const TGSegment> retainedTuples;
            size_t oldNodeId = derivations->getNodeId();
            if (shouldCleanDuplicates && !node.retainFree) {
                retainedTuples = g.retain(currentPredicate, derivations);
            } else {
                retainedTuples = derivations;
            }
            nonempty = !(retainedTuples == NULL || retainedTuples->isEmpty());
            if (nonempty) {
                nders += retainedTuples->getNRows();
                if (rule.isEGD()) {
                    g.replaceEqualTerms(node.ruleIdx, node.step, retainedTuples);
                } else {
                    //Add new nodes
                    if (trackProvenance) {
                        createNewNodesWithProv(node.ruleIdx, node.step,
                                retainedTuples, derivationNodes);
                    } else {
                        //Add a single node
                        g.addNodeNoProv(currentPredicate, node.ruleIdx,
                                node.step, retainedTuples);
                    }
                }
            }
        }
    }

    if (shouldStoreStats()) {
        std::chrono::duration<double, std::milli> retainRuntime =
            std::chrono::system_clock::now() - starth;
        std::chrono::duration<double, std::milli> totalRuntime =
            std::chrono::system_clock::now() - start;

        StatsRule stats;
        stats.step = node.step;
        stats.idRule = node.ruleIdx;
        stats.nderivations_final = nders;
        stats.nderivations_unfiltered = nders_un;
        stats.timems = totalRuntime.count();
        stats.timems_first = executor.getDuration(DurationType::DUR_FIRST).count();
        stats.timems_merge = executor.getDuration(DurationType::DUR_MERGE).count();
        stats.timems_join = executor.getDuration(DurationType::DUR_JOIN).count();
        stats.timems_createhead = executor.getDuration(DurationType::DUR_HEAD).count();
        stats.timems_retain = retainRuntime.count();
        stats.nbdyatoms = executor.getStat(StatType::N_BDY_ATOMS);
        saveStatistics(stats);

    }

    return nonempty;
}

void GBChase::createNewNodesWithProv(size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> seg,
        std::vector<std::shared_ptr<Column>> &provenance) {
    if (provenance.size() == 0) {
        //Single EDB or IDB atom
        if (seg->getProvenanceType() == 2) {
            //Must split the nodes
            auto resortedSeg = seg->sortByProv();

            std::vector<size_t> provNodes;
            auto chunks = resortedSeg->sliceByNodes(g.getNNodes(),
                    provNodes);
            assert(chunks.size() == provNodes.size());
            std::vector<size_t> currentNodeList(1);
            for (size_t i = 0; i < chunks.size(); ++i) {
                auto c = chunks[i];
                currentNodeList[0] = provNodes[i];
                assert(currentNodeList[0] != ~0ul);
                g.addNodeProv(currentPredicate, ruleIdx, step, c,
                        currentNodeList);
            }

            /*size_t startidx = 0;
              auto itr = resortedSeg->iterator();
              size_t i = 0;
              size_t currentNode = ~0ul;
              std::vector<size_t> currentNodeList(1);
              while (itr->hasNext()) {
              itr->next();
              bool hasChanged = i == 0 || currentNode != itr->getNodeId();
              if (hasChanged) {
              if (startidx < i) {
            //Create a new node
            currentNodeList[0] = currentNode;
            auto nodeId = g.getNNodes();
            auto dataToAdd = resortedSeg->slice(nodeId, startidx, i);
            g.addNodeProv(currentPredicate, ruleIdx, step, dataToAdd,
            currentNodeList);
            }
            startidx = i;
            currentNode = itr->getNodeId();
            }
            i++;
            }
            //Copy the last segment
            if (startidx < i) {
            currentNodeList[0] = currentNode;
            auto nodeId = g.getNNodes();
            auto dataToAdd = resortedSeg->slice(nodeId, startidx, i);
            g.addNodeProv(currentPredicate, ruleIdx,
            step, dataToAdd, currentNodeList);
            }*/
        } else {
            std::vector<size_t> provnodes;
            if (seg->getNodeId() == ~0ul) {
                //EDB body atom
#if DEBUG
                assert(rules[ruleIdx].getBody().size() == 1 &&
                        rules[ruleIdx].getBody()[0].getPredicate().getType() ==
                        EDB);
#endif
            } else {
                provnodes.push_back(seg->getNodeId());
            }
            auto nodeId = g.getNNodes();
            auto dataToAdd = seg->slice(nodeId, 0, seg->getNRows());
            g.addNodeProv(currentPredicate, ruleIdx, step, dataToAdd, provnodes);
        }
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
                const size_t m = providxs[i] * nnodes + j;
                if (currentNodeList[j] != provnodes[m]) {
                    hasChanged = true;
                }
            }
            if (hasChanged) {
                if (startidx < i) {
                    //Create a new node
                    auto nodeId = g.getNNodes();
                    auto dataToAdd = resortedSeg->slice(nodeId, startidx, i);
                    g.addNodeProv(currentPredicate, ruleIdx, step, dataToAdd,
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
            auto nodeId = g.getNNodes();
            auto dataToAdd = resortedSeg->slice(nodeId, startidx, nrows);
            g.addNodeProv(currentPredicate, ruleIdx,
                    step, dataToAdd, currentNodeList);
        }
    }
}

size_t GBChase::getSizeTable(const PredId_t predid) const {
    if (g.areNodesWithPredicate(predid)) {
        const auto &nodeIDs = g.getNodeIDsWithPredicate(predid);
        size_t size = 0;
        for(auto nodeID : nodeIDs) {
            auto data = g.getNodeSize(nodeID);
            size += data;
        }
        return size;
    }
    return 0;
}

FCIterator GBChase::getTableItr(const PredId_t predid) {
    FCTable *t = getTable(predid);
    cacheFCTables.insert(std::make_pair(predid, std::shared_ptr<FCTable>(t)));
    return t->read(0);
}

FCTable *GBChase::getTable(const PredId_t predid) {
    //Create a FCTable obj with the nodes in predid
    uint8_t card = program->getPredicateCard(predid);
    VTuple tuple(card);
    Literal query(program->getPredicate(predid), tuple);
    FCTable *t = new FCTable(NULL, card);
    if (g.areNodesWithPredicate(predid)) {
        auto &nodeIDs = g.getNodeIDsWithPredicate(predid);
        for (auto &nodeId : nodeIDs) {
            auto nodeData = g.getNodeData(nodeId);
            auto data = fromTGSeg2Seg(nodeData);
            std::shared_ptr<const FCInternalTable> table =
                std::shared_ptr<const FCInternalTable>(
                        new InmemoryFCInternalTable(card, nodeId, true, data));
            FCBlock block(nodeId, table, query, 0, NULL, 0, true);
            t->addBlock(block);
        }
    }
    return t;
}

size_t GBChase::getCurrentIteration() {
    return currentIteration;
}

#ifdef WEBINTERFACE
std::string GBChase::getCurrentRule() {
    return currentRule;
}

PredId_t GBChase::getCurrentPredicate() {
    return currentPredicate;
}
#endif
