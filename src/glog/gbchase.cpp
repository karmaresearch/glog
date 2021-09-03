#include <glog/gbchase.h>
#include <glog/gbsegmentcache.h>
#include <glog/gbquerier.h>

#include <unordered_set>

GBChase::GBChase(EDBLayer &layer, Program *program, bool useCacheRetain,
        GBGraph::ProvenanceType provenanceType,
        bool filterQueryCont,
        bool edbCheck,
        bool rewriteCliques,
        bool duplAllowed) :
    provenanceType(provenanceType),
    filterQueryCont(filterQueryCont),
    edbCheck(edbCheck),
    program(program),
    layer(layer),
    g(provenanceType, useCacheRetain, filterQueryCont, duplAllowed),
    triggers(0),
    durationRuleExec(0),
    currentIteration(0),
    startStep(0),
    maxStep(~0ul),
    lastStep(0),
    durationPreparation(0)
{
    LOG(INFOL) << "Query cont=" << filterQueryCont <<
        " EDB check=" << edbCheck;
    if (rewriteCliques) {
        //Rewrite rules of the form P(X,Y),P(Y,Z)->P(X,Z) and P(X,Y)->P(Y,X)
        program->rewriteCliques();
    }

    if (!program->stratify(stratification, nStratificationClasses)) {
        LOG(ERRORL) << "Program could not be stratified";
        throw std::runtime_error("Program could not be stratified");
    }
    LOG(DEBUGL) << "nStratificationClasses = " << nStratificationClasses;

    rules = program->getAllRules();
    g.setRulesProgramLayer(rules.data(), program, &layer);
    executor = std::unique_ptr<
        GBRuleExecutor>(new GBRuleExecutor(g, layer, program));
}

Program *GBChase::getProgram() {
    return program;
}

EDBLayer &GBChase::getEDBLayer() {
    return layer;
}

GBGraph &GBChase::getGBGraph() {
    return g;
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

void GBChase::prepareRuleExecutionPlans_queryContainment(
        std::vector<GBRuleInput> &newnodes,
        std::vector<std::vector<size_t>> &acceptableNodes,
        const size_t ruleIdx,
        const size_t step)
{
    size_t nCombinations = 1;
    for(size_t i = 0; i < acceptableNodes.size(); ++i) {
        nCombinations = nCombinations * acceptableNodes[i].size();
    }
    if (nCombinations > 10) {
        LOG(DEBUGL) << "Skipping checking redundancies " << acceptableNodes.size() <<
            " " << nCombinations << " ...";
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
                                edbCheck, rFree)) {
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
            newnode.retainFree = false;
        }
    }
}

void GBChase::prepareRuleExecutionPlans_SNE(
        const std::vector<std::pair<bool, std::vector<size_t>>> nodesForRule,
        size_t ruleIdx,
        size_t step,
        size_t prevstep,
        std::vector<GBRuleInput> &newnodes)
{
    for(int pivot = 0; pivot < nodesForRule.size(); ++pivot) {
        if (nodesForRule[pivot].first) { //Is negated, skip
            continue;
        }
        //First consider only combinations where at least one node
        //is derived in the previous step (semi-naive evaluation)
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
            //The implementation of query containment does not work if the
            //rule contains constants in the body
            bool canCheck = true;
            for (auto &b : rules[ruleIdx].getBody()) {
                if (b.getNConstants() > 0) {
                    canCheck = false;
                    break;
                }
            }
            if (canCheck) {
                prepareRuleExecutionPlans_queryContainment(newnodes,
                        acceptableNodes,
                        ruleIdx,
                        step);
            } else {
                newnodes.emplace_back();
                GBRuleInput &newnode = newnodes.back();
                newnode.ruleIdx = ruleIdx;
                newnode.step = step;
                newnode.incomingEdges = acceptableNodes;
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
        if (rule.isTransitive()) {
            bool isOptimizationApplicable = true;
            assert(nodesForRule.size() == 2);
            for(size_t i = 0; i < nodesForRule.size(); ++i) {
                if (nodesForRule[i].first) { //is negated
                    isOptimizationApplicable = false;
                    break;
                }
            }
            if (!isOptimizationApplicable) {
                prepareRuleExecutionPlans_SNE(nodesForRule,
                        ruleIdx, step, prevstep, newnodes);
            } else {
                //Get the nodes not produced by the transitive rule
                std::vector<size_t> leftSidePrevStep;
                std::vector<size_t> leftSideBeforePrevStep;
                for (auto &nodeId : nodesForRule[0].second) {
                    auto nodeRuleIdx = g.getNodeRuleIdx(nodeId);
                    if (nodeRuleIdx != ruleIdx) {
                        auto nodeStep = g.getNodeStep(nodeId);
                        if (nodeStep == prevstep) {
                            leftSidePrevStep.push_back(nodeId);
                        } else {
                            leftSideBeforePrevStep.push_back(nodeId);
                        }
                    }
                }

                if (leftSideBeforePrevStep.size() > 0) {
                    std::vector<std::vector<size_t>> acceptableNodes;
                    acceptableNodes.resize(nodesForRule.size());
                    acceptableNodes[0] = leftSideBeforePrevStep;
                    for (auto &nodeId : nodesForRule[1].second) {
                        if (g.getNodeStep(nodeId) == prevstep) {
                            acceptableNodes[1].push_back(nodeId);
                        }
                    }
                    if (!acceptableNodes[0].empty() &&
                            !acceptableNodes[1].empty()) {
                        newnodes.emplace_back();
                        GBRuleInput &newnode = newnodes.back();
                        newnode.ruleIdx = ruleIdx;
                        newnode.step = step;
                        newnode.incomingEdges = acceptableNodes;
                    }
                }
                if (leftSidePrevStep.size() > 0) {
                    std::vector<std::vector<size_t>> acceptableNodes;
                    acceptableNodes.resize(nodesForRule.size());
                    acceptableNodes[0] = leftSidePrevStep;
                    for (auto &nodeId : nodesForRule[1].second) {
                        acceptableNodes[1].push_back(nodeId);
                    }
                    if (!acceptableNodes[0].empty() &&
                            !acceptableNodes[1].empty()) {
                        newnodes.emplace_back();
                        GBRuleInput &newnode = newnodes.back();
                        newnode.ruleIdx = ruleIdx;
                        newnode.step = step;
                        newnode.incomingEdges = acceptableNodes;
                    }
                }
            }
        } else {
            prepareRuleExecutionPlans_SNE(nodesForRule,
                    ruleIdx, step, prevstep, newnodes);
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
    bool edbCanChange = false;
    for (auto &bodyAtom : rule.getBody()) {
        Predicate pred = bodyAtom.getPredicate();
        if (pred.getType() != EDB) {
            bool negated = bodyAtom.isNegated();
            if (!g.areNodesWithPredicate(pred.getId()) && !negated) {
                empty = true;
                break;
            }
        } else {
            if (layer.canChange(pred.getId())) {
                edbCanChange = true;
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
        //I make an exception if the EDB predicates can change
        if (!edbCanChange && (step != stepStratum + 1 || stratumLevel != 0)) {
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

    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
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
    durationPreparation +=  std::chrono::system_clock::now() - start;

    //Execute the rule associated to the node
    auto nnodes = g.getNNodes();
    auto nodesToProcess = newnodes.size();

    LOG(INFOL) << "Nodes to process " << nodesToProcess;
    for(size_t idxNode = 0; idxNode < nodesToProcess; ++idxNode) {
        executeRule(newnodes[idxNode]);
    }

    //Retain all the predicates that should be cleaned at the end
    for (auto &predId : predToBeRetainedEndStep) {
        g.retainAndAddFromTmpNodes(predId);
    }

    if (nnodes != g.getNNodes()) {
        auto derivedTuples = 0;
        for(size_t idx = nnodes; idx < g.getNNodes(); ++idx) {
            auto nrows = g.getNodeSize(idx);
            derivedTuples += nrows;
        }
        LOG(INFOL) << "Derived Tuples: " << derivedTuples;
        LOG(INFOL) << "#nodes: " << g.getNNodes();
        return derivedTuples;
    } else {
        return 0;
    }
}

void GBChase::prepareRun(size_t startStep, size_t maxStep) {
    this->startStep = startStep;
    this->maxStep = maxStep;
}

void GBChase::run() {
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    initRun();
    size_t nnodes = 0;
    size_t step = startStep;

    //Mark the predicates that should be cleaned at the end
    for (auto &predId : program->getAllPredicateIDs()) {
        if (program->getNRulesByPredicate(predId) > 100) {
            predToBeRetainedEndStep.insert(predId);
        }
    }

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
            if (step == maxStep)
                break;
            layer.setContext(&g, step);
            LOG(INFOL) << "Step " << step;
            currentIteration = step;
            g.cleanTmpNodes();
            nnodes = g.getNNodes();
            size_t derivedTuples = executeRulesInStratum(rulesInStratum,
                    currentStrat, stepStratum, step);
        } while (g.getNNodes() != nnodes);
        g.cleanTmpNodes();
    }

    lastStep = step;

    std::chrono::duration<double, std::milli> dur =
        std::chrono::system_clock::now() - start;
    LOG(INFOL) << "Runtime chase: " << dur.count();

    layer.clearContext();
    SegmentCache::getInstance().clear();
    executor->printStats();
    LOG(INFOL) << "(GBChase) Time stratum preparation (ms): " << durationPreparation.count();
    LOG(INFOL) << "(GBChase) Time rule exec (ms): " << durationRuleExec.count();
    //LOG(INFOL) << "(GBChase) Time debug (ms): " << durationDebug.count();
    g.printStats();
    stopRun();

    /*std::cout << "Testing derivation trees..." << std::endl;
      GBQuerier q = GBQuerier(g, *program, layer);
      auto out = q.getDerivationTree(0, 0);*/
}

size_t GBChase::getNDerivedFacts() {
    return g.getNDerivedFacts();
}

size_t GBChase::getNnodes() {
    return g.getNNodes();
}

size_t GBChase::getNedges() {
    return g.getNEdges();
}

size_t GBChase::getNTriggers() {
    return triggers;
}

size_t GBChase::getNSteps() {
    return lastStep;
}

bool GBChase::shouldRetainAtEnd(PredId_t pred) {
    bool outcome = predToBeRetainedEndStep.count(pred);
    return outcome;
}

bool GBChase::executeRule(GBRuleInput &node, bool cleanDuplicates) {
    auto &bodyNodes = node.incomingEdges;
    Rule &rule = rules[node.ruleIdx];
#ifdef WEBINTERFACE
    currentRule = rule.tostring();
#endif

//#ifdef DEBUG
    LOG(DEBUGL) << "Execute rule " << node.ruleIdx << " " <<
        rule.tostring(program, &layer);
//#endif

    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();

    size_t nders = 0;
    size_t nders_un = 0;
    auto &heads = rule.getHeads();
    int headIdx = 0;
    currentPredicate = heads[headIdx].getPredicate().getId();
    auto outputsRule = executor->executeRule(rule, node);
    bool nonempty = false;
    std::chrono::duration<double, std::milli> execRuntime =
        std::chrono::system_clock::now() - start;
    durationRuleExec += execRuntime;

    std::chrono::system_clock::time_point starth =
        std::chrono::system_clock::now();
    for (GBRuleOutput &outputRule : outputsRule) {
        currentPredicate = heads[headIdx++].getPredicate().getId();
        auto derivations = outputRule.segment;
        nders_un += derivations->getNRows();
        triggers += derivations->getNRows();
        auto derivationNodes = outputRule.nodes;
        const bool shouldCleanDuplicates = cleanDuplicates &&
            !outputRule.uniqueTuples;
        nonempty |= !(derivations == NULL || derivations->isEmpty());

        if (nonempty) {
            if (shouldRetainAtEnd(currentPredicate) && !rule.isEGD()) {
                g.addNodeToBeRetained(currentPredicate, derivations,
                        derivationNodes, node.ruleIdx, node.step);
            } else {
                //Keep only the new derivations
                std::shared_ptr<const TGSegment> retainedTuples;
                if (shouldCleanDuplicates && !node.retainFree) {
                    retainedTuples = g.retain(currentPredicate, derivations,
                            derivationNodes);
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
                        if (shouldTrackProvenance()) {
                            g.addNodesProv(currentPredicate, node.ruleIdx,
                                    node.step, retainedTuples, derivationNodes);
                        } else {
                            //Add a single node
                            g.addNodeNoProv(currentPredicate, node.ruleIdx,
                                    node.step, retainedTuples);
                        }
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
        stats.timems_first = executor->getDuration(DurationType::DUR_FIRST).count();
        stats.timems_merge = executor->getDuration(DurationType::DUR_MERGE).count();
        stats.timems_join = executor->getDuration(DurationType::DUR_JOIN).count();
        stats.timems_createhead = executor->getDuration(DurationType::DUR_HEAD).count();
        stats.timems_retain = retainRuntime.count();
        stats.nbdyatoms = executor->getStat(StatType::N_BDY_ATOMS);
        saveStatistics(stats);
    }

    return nonempty;
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
            std::shared_ptr<const FCInternalTable> table;
            if (nodeData->getNColumns() > 0) {
                auto data = fromTGSeg2Seg(nodeData);
                table = std::shared_ptr<const FCInternalTable>(
                        new InmemoryFCInternalTable(card, nodeId, true, data));
            } else {
                table = std::shared_ptr<const FCInternalTable>(new SingletonTable(nodeId));
            }
            FCBlock block(nodeId, table, query, 0, NULL, 0, true);
            t->addBlock(block);
        }
    }
    return t;
}

size_t GBChase::getCurrentIteration() {
    return currentIteration;
}


std::vector<GBRuleOutput> GBChase::executeRule(size_t ruleIdx)
{
    Rule rule = program->getRule(ruleIdx);
    GBRuleInput in;
    in.ruleIdx = ruleIdx;
    in.step = ~0ul;
    in.retainFree = true;
    //Get all the nodes
    auto &ruleBody = rule.getBody();
    for(auto literalBody : ruleBody) {
        if (literalBody.getPredicate().getType() != EDB) {
            auto predId = literalBody.getPredicate().getId();
            auto nnodes = g.getNodeIDsWithPredicate(predId);
            if (nnodes.empty()) {
                return std::vector<GBRuleOutput>();
            }
            in.incomingEdges.push_back(nnodes);
        }
    }
    auto output = executor->executeRule(rule, in);
    return output;
}

#ifdef WEBINTERFACE
std::string GBChase::getCurrentRule() {
    return currentRule;
}

PredId_t GBChase::getCurrentPredicate() {
    return currentPredicate;
}
#endif
