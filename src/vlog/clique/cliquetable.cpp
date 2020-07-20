#include <vlog/clique/cliquetable.h>

CliqueTable::CliqueTable(PredId_t predid, PredId_t targetPredicate) :
    predid(predid), targetPredicate(targetPredicate), componentIDCounter(0),
    recompute(true) {
    }

void CliqueTable::query(QSQQuery *query, TupleTable *outputTable,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t CliqueTable::estimateCardinality(const Literal &query) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

size_t CliqueTable::getCardinality(const Literal &query) {
    computeConnectedComponents();
    size_t card = 0;
    for (auto &p : components) {
        card += p.second.size() * p.second.size();
    }
    return card;
}

size_t CliqueTable::getCardinalityColumn(const Literal &query,
        uint8_t posColumn) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

bool CliqueTable::isEmpty(const Literal &query,
        std::vector<uint8_t> *posToFilter,
        std::vector<Term_t> *valuesToFilter) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *CliqueTable::getIterator(const Literal &query) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

EDBIterator *CliqueTable::getSortedIterator(const Literal &query,
        const std::vector<uint8_t> &fields) {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

void CliqueTable::releaseIterator(EDBIterator *itr) {
    delete itr;
}

bool CliqueTable::getDictNumber(const char *text, const size_t sizeText,
        uint64_t &id) {
    return false;
}

bool CliqueTable::getDictText(const uint64_t id, char *text) {
    return false;
}

bool CliqueTable::getDictText(const uint64_t id, std::string &text) {
    return false;
}

uint64_t CliqueTable::getNTerms() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint64_t CliqueTable::getSize() {
    LOG(ERRORL) << "Not implemented";
    throw 10;
}

uint8_t CliqueTable::getArity() const {
    return 2;
}

void CliqueTable::setContext(GBGraph *g, size_t step) {
    this->g = g;
    this->step = step;
    recompute = true;
}

void CliqueTable::clearContext() {
    this->g = NULL;
    this->step = 0;
    components.clear();
    term2component.clear();
    componentIDCounter = 0;
    recompute = true;
}

CliqueIterator *CliqueTable::iterator() {
    CliqueIterator *itr = new CliqueIterator(predid,
            components,
            term2component.begin(),
            term2component.end());
    return itr;
}

std::vector<std::pair<Term_t, Term_t>> CliqueTable::checkNewIn(
        const Literal &l1,
        std::vector<uint8_t> &posInL1,
        const std::vector<std::pair<Term_t, Term_t>> &existing) {
    assert(!recompute);

    //Return all the new tuples that are not in existing
    assert(posInL1.size() == 2 && posInL1[0] == 0 && posInL1[1] == 1);
    assert(l1.getTermAtPos(0).isVariable() && l1.getTermAtPos(1).isVariable());

    auto itr = iterator();
    std::vector<std::pair<Term_t, Term_t>> out;
    //Term_t prevnewt[2] = { 0, 0};
    Term_t newt[2];
    newt[0] = newt[1] = ~0ul;

    Term_t ext[2];
    ext[0] = ext[1] = ~0ul;

    size_t i = 0;
    while (true) {
        if (newt[0] == ~0ul) {
            if (itr->hasNext()) {
                itr->next();
                newt[0] = itr->getElementAt(0);
                newt[1] = itr->getElementAt(1);

                /*if (prevnewt[0] > newt[0] ||
                        (prevnewt[0] == newt[0] && prevnewt[1] > newt[1])) {
                    throw 10;
                } else if (prevnewt[0] == newt[0] &&
                        prevnewt[1] == newt[1]) {
                    throw 10;
                }
                prevnewt[0] = newt[0];
                prevnewt[1] = newt[1];*/
            } else {
                break;
            }
        }
        if (ext[0] == ~0ul) {
            if (i == existing.size()) {
                break;
            } else {
                ext[0] = existing[i].first;
                ext[1] = existing[i].second;
                i++;
            }
        }
        if (newt[0] < ext[0] || (newt[0] == ext[0] &&
                    newt[1] < ext[1])) {
            //new
            out.push_back(std::make_pair(newt[0], newt[1]));
            newt[0] = newt[1] = ~0ul;
        } else if (newt[0] == ext[0] && newt[1] == ext[1]) {
            //existing
            newt[0] = newt[1] = ~0ul;
            ext[0] = ext[1] = ~0ul;
        } else {
            ext[0] = ext[1] = ~0ul;
        }
    }

    if (newt[0] != ~0ul) {
        out.push_back(std::make_pair(newt[0], newt[1]));
        while (itr->hasNext()) {
            itr->next();
            out.push_back(std::make_pair(itr->getElementAt(0),
                        itr->getElementAt(1)));
        }
    }
    //LOG(WARNL) << "NEW TERMS for pred " << predid << " " << out.size();

    releaseIterator(itr);
    return out;
}

void CliqueTable::computeConnectedComponents() {
    if (!recompute)
        return;

    //Take the content of all nodes produced in the previous step
    //Check each pair, if they are in the same component, do nothing. Otherwise
    //merge different components into a bigger one
    if (g->areNodesWithPredicate(targetPredicate)) {
        const auto &nodes = g->getNodeIDsWithPredicate(targetPredicate);
        assert(nodes.size() > 0);
        assert(step > 0);
        LOG(DEBUGL) << "Creating components for " << nodes.size() << " nodes";
        std::set<size_t> modifiedComponents;
        for(int64_t i = nodes.size() - 1; i >= 0; i--) {
            auto nodeId = nodes[i];
            if (g->getNodeStep(nodeId) < step - 1) {
                break;
            }
            auto itr = g->getNodeIterator(nodeId);
            while (itr->hasNext()) {
                assert(itr->getNFields() == 2);
                itr->next();
                auto term1 = itr->get(0);
                auto term2 = itr->get(1);
                if (!term2component.count(term1) &&
                        !term2component.count(term2)) {
                    //Create a new component
                    term2component.insert(std::make_pair(term1,
                                componentIDCounter));
                    term2component.insert(std::make_pair(term2,
                                componentIDCounter));
                    components.insert(std::make_pair(componentIDCounter,
                                std::vector<Term_t>()));
                    components[componentIDCounter].push_back(term1);
                    //if (term1 != term2)
                        components[componentIDCounter].push_back(term2);
                    componentIDCounter++;
                } else if (!term2component.count(term2)) {
                    //Add term2 to the component of term1
                    auto cId = term2component[term1];
                    components[cId].push_back(term2);
                    term2component.insert(std::make_pair(term2, cId));
                    modifiedComponents.insert(cId);
                } else if (!term2component.count(term1)) {
                    //Add term1 to the component of term2
                    auto cId = term2component[term2];
                    components[cId].push_back(term1);
                    term2component.insert(std::make_pair(term1, cId));
                    modifiedComponents.insert(cId);
                } else {
                    auto cId1 = term2component[term1];
                    auto cId2 = term2component[term2];
                    if (cId1 != cId2) {
                        //Must merge two components
                        auto &c1 = components[cId1];
                        auto &c2 = components[cId2];
                        for (auto c : c2) {
                            c1.push_back(c);
                            term2component[c] = cId1;
                        }
                        components.erase(cId2);
                        modifiedComponents.insert(cId1);
                    }
                }
            }
        }
        for(auto cId : modifiedComponents) {
            auto &v = components[cId];
            size_t oldSize = v.size();
            std::sort(v.begin(), v.end());
            auto e = std::unique(v.begin(), v.end());
            v.erase(e, v.end());
        }
        LOG(DEBUGL) << predid << "N. components " << components.size();
        LOG(DEBUGL) << predid << "term2component " << term2component.size();
    }
    recompute = false;
}
