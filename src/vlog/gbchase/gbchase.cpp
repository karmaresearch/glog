#include <vlog/gbchase.h>
#include <vlog/tgsegmentcache.h>

GBChase::GBChase(EDBLayer &layer, Program *program, bool useCacheRetain) :
    layer(layer),
    program(program),
    trackProvenance(true),
    g(trackProvenance, useCacheRetain),
    executor(trackProvenance, g, layer) {
        if (!program->stratify(stratification, nStratificationClasses)) {
            LOG(ERRORL) << "Program could not be stratified";
            throw std::runtime_error("Program could not be stratified");
        }
        LOG(DEBUGL) << "nStratificationClasses = " << nStratificationClasses;
        /*for (auto &r : program->getAllRules()) {
           if (r.getHeads().size() > 1) {
                LOG(ERRORL) << "GBChase does not support the execution of rules with more atoms in the head. set --rerwriteMultihead to true to rewrite these rulese";
                throw std::runtime_error("Program could not be stratified");
            }
        }*/
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

bool lowerStrat(Rule &rule, int currentStrat, std::vector<int> &stratification) {
    for (auto &bodyAtom : rule.getBody()) {
        Predicate pred = bodyAtom.getPredicate();
        if (pred.getType() != EDB && stratification[pred.getId()] >= currentStrat) {
            return false;
        }
    }
    return true;
}

void GBChase::run() {
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
            nnodes = g.getNNodes();

            std::vector<GBRuleInput> newnodes;
            for (size_t ruleIdx = 0; ruleIdx < rules.size(); ++ruleIdx) {
                auto &rule = rules[ruleIdx];
                PredId_t id = rule.getFirstHead().getPredicate().getId();
                if (nStratificationClasses > 1 && stratification[id] != currentStrat) {
                    LOG(DEBUGL) << "Skipping rule, wrong strat: " << rule.tostring();
                    continue;
                }
                LOG(DEBUGL) << "Considering rule " << rule.tostring();
                //The first parameter in the pair records whether the associated
                //atom is negated. In this case, we do not apply semi-naive
                //evaluation
                std::vector<std::pair<bool,std::vector<size_t>>> nodesForRule;
                bool empty = false;
                for (auto &bodyAtom : rule.getBody()) {
                    Predicate pred = bodyAtom.getPredicate();
                    if (pred.getType() != EDB) {
                        bool negated = bodyAtom.isNegated();
                        if (!g.areNodesWithPredicate(pred.getId()) && !negated) {
                            empty = true;
                            break;
                        }
                        nodesForRule.push_back(std::make_pair(negated,
                                    g.getNodeIDsWithPredicate(pred.getId())));
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
                        GBRuleInput &newnode = newnodes.back();
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
                        /*for(int j = 0; j < nodesForRule.size(); ++j) {
                          if (!nodesForRule[j].second.empty()) {
                          LOG(DEBUGL) << "Pushing node for lowerStrat rule " << rule.tostring();
                          newnodes.emplace_back();
                          GBRuleInput &newnode = newnodes.back();
                          newnode.ruleIdx = ruleIdx;
                          newnode.step = step;
                          newnode.incomingEdges = std::vector<std::vector<size_t>>();
                          newnode.incomingEdges.push_back(nodesForRule[j]);
                          }
                          }*/
                        //Jacopo: I think the code above is wrong... put a message
                        //to make sure we check it
                        LOG(ERRORL) << "This case is not implemented. FIXME";
                        throw 10;
                    } else {
                        LOG(DEBUGL) << "Skipping lowerStrat rule";
                    }
                } else {
                    size_t prevstep = step - 1;
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
                                    if (nodeStep == prevstep) {
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

                        LOG(DEBUGL) << "Pushing node for rule " << rule.tostring();
                        newnodes.emplace_back();
                        GBRuleInput &newnode = newnodes.back();
                        newnode.ruleIdx = ruleIdx;
                        newnode.step = step;
                        newnode.incomingEdges = acceptableNodes;
                    }
                }
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
                    auto nrows = g.getNodeData(idx)->getNRows();
                    derivedTuples += nrows;
                }
                LOG(INFOL) << "Derived Tuples: " << derivedTuples;
            }
        } while (nnodes != g.getNNodes());
    }

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

    auto &heads = rule.getHeads();
    int headIdx = 0;
    currentPredicate = heads[headIdx].getPredicate().getId();
    auto outputsRule = executor.executeRule(rule, node);
    bool nonempty = false;
    for (GBRuleOutput &outputRule : outputsRule) {
        currentPredicate = heads[headIdx++].getPredicate().getId();
        auto derivations = outputRule.segment;
        auto derivationNodes = outputRule.nodes;
        const bool shouldCleanDuplicates = cleanDuplicates && !outputRule.uniqueTuples;
        nonempty |= !(derivations == NULL || derivations->isEmpty());
        if (nonempty) {
            //Keep only the new derivations
            std::shared_ptr<const TGSegment> retainedTuples;
            if (shouldCleanDuplicates) {
                retainedTuples = g.retain(currentPredicate, derivations);
            } else {
                retainedTuples = derivations;
            }
            nonempty = !(retainedTuples == NULL || retainedTuples->isEmpty());
            if (nonempty) {
                if (trackProvenance) {
                    createNewNodesWithProv(node.ruleIdx, node.step,
                            retainedTuples, derivationNodes);
                } else {
                    //Add a single node
                    g.addNode(currentPredicate, node.ruleIdx,
                            node.step, retainedTuples);
                }
            }
        }
    }
    return nonempty;
}

void GBChase::createNewNodesWithProv(size_t ruleIdx, size_t step,
        std::shared_ptr<const TGSegment> seg,
        std::vector<std::shared_ptr<Column>> &provenance) {
    if (provenance.size() == 0) {
        //There was no join. Must replace the nodeID with a new one
        //Note at this point calling slice will create a new vector
        auto nodeId = g.getNNodes();
        auto dataToAdd = seg->slice(nodeId, 0, seg->getNRows());
        g.addNode(currentPredicate, ruleIdx, step, dataToAdd);
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
                    auto nodeId = g.getNNodes();
                    auto dataToAdd = resortedSeg->slice(nodeId, startidx, i);
                    g.addNode(currentPredicate, ruleIdx, step, dataToAdd);
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
            auto nodeId = g.getNNodes();
            auto dataToAdd = resortedSeg->slice(nodeId, startidx, nrows);
            g.addNode(currentPredicate, ruleIdx, step, dataToAdd);
        }
    }
}

size_t GBChase::getSizeTable(const PredId_t predid) const {
    if (g.areNodesWithPredicate(predid)) {
        const auto &nodeIDs = g.getNodeIDsWithPredicate(predid);
        size_t size = 0;
        for(auto nodeID : nodeIDs) {
            auto data = g.getNodeData(nodeID);
            size += data->getNRows();

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
