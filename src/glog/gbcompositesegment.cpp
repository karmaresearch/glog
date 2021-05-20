#include <glog/gbcompositesegment.h>

std::shared_ptr<const TGSegment> CompositeTGSegment::merge() const {
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    std::shared_ptr<const TGSegment> t = g.mergeNodes(nodes, copyVarPos,
            false, replaceOffsets);
    std::chrono::system_clock::time_point end =
        std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> dur = end - start;
    return t;
}

size_t CompositeTGSegment::getNRows() const {
    assert(nodes.size() > 0);
    assert(copyVarPos.size() == 1 || copyVarPos[0] != copyVarPos[1]);
    if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns()) {
        if (nodes.size() == 1) {
            auto nodeData = g.getNodeData(nodes[0]);
            if (copyVarPos.size() == 2) {
                std::vector<std::shared_ptr<Column>> out;
                std::vector<int> columns;
                columns.push_back(copyVarPos[0]);
                columns.push_back(copyVarPos[1]);
                nodeData->projectTo(columns, out);
                if (out[0]->isEDB() && out[1]->isEDB()) {
                    EDBColumn *c1 = (EDBColumn*) out[0].get();
                    EDBColumn *c2 = (EDBColumn*) out[1].get();
                    EDBLayer &layer = c1->getEDBLayer();
                    const Literal &l1 = c1->getLiteral();
                    const Literal &l2 = c2->getLiteral();
                    auto card = layer.getCardinality(l1);
                    return card;
                }
            }
        }
        LOG(ERRORL) << "Case not yet implemented";
        throw 10;
    } else {
        size_t nrows = 0;
        for(auto n : nodes)
            nrows += g.getNodeSize(n);
        return nrows;
    }
}

std::unique_ptr<TGSegmentItr> CompositeTGSegment::iterator(
        std::shared_ptr<const TGSegment> selfref) const {
    auto mergedSegment = merge();
    if (sortBeforeAccess) {
        mergedSegment = mergedSegment->sort();
    }
    if (removeDuplBeforeAccess) {
        mergedSegment = mergedSegment->unique();
    }
    return mergedSegment->iterator(mergedSegment);
}

bool CompositeTGSegment::isSortedBy(std::vector<uint8_t> &fields) const {
    return false;
}

std::shared_ptr<const TGSegment> CompositeTGSegment::sort() const {
    assert(nodes.size() > 0);
    assert(copyVarPos.size() == 1 || copyVarPos[0] != copyVarPos[1]);
    /*bool seq = true;
      for(int i = 0; i < copyVarPos.size(); ++i) {
      if (copyVarPos[i] != i) {
      seq = false;
      break;
      }
      }

      if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns() ||
      !seq) {*/
    auto mergedSegment = merge();
    if (mergedSegment->isSorted())
        return mergedSegment;
    else
        return mergedSegment->sort();
    /*} else {
      throw 10; //Not sure that this code below is correct. Are we sure
    //the segment is always sorted?
    return std::shared_ptr<const TGSegment>(new CompositeTGSegment(g,
    nodes,
    copyVarPos,
    f_isSorted,
    sortedField,
    provenanceType,
    !isSorted(),
    false,
    replaceOffsets));
    }*/
}

void CompositeTGSegment::argsort(std::vector<size_t> &indices) const {
    auto mergedSegment = merge();
    mergedSegment->argsort(indices);
}

void CompositeTGSegment::argunique(std::vector<size_t> &indices) const {
    auto mergedSegment = merge();
    mergedSegment->argunique(indices);
}

std::shared_ptr<TGSegment> CompositeTGSegment::sortBy(
        std::vector<uint8_t> &fields) const {
    auto mergedSegment = merge();
    return mergedSegment->sortBy(fields);
}

std::shared_ptr<const TGSegment> CompositeTGSegment::sortByProv(size_t ncols,
        std::vector<size_t> &idxs,
        std::vector<size_t> &nodes) const {
    auto mergedSegment = merge();
    return mergedSegment->sortByProv(ncols, idxs, nodes);
}

std::shared_ptr<const TGSegment> CompositeTGSegment::sortByProv() const {
    std::vector<size_t> sortedNodes = nodes;
    std::sort(sortedNodes.begin(), sortedNodes.end());
    return std::shared_ptr<const TGSegment>(new CompositeTGSegment(g,
                sortedNodes,
                copyVarPos,
                f_isSorted,
                sortedField,
                provenanceType,
                sortBeforeAccess,
                removeDuplBeforeAccess,
                replaceOffsets));
}

std::shared_ptr<const TGSegment> CompositeTGSegment::unique() const {
    assert(nodes.size() > 0);
    assert(copyVarPos.size() == 1 || copyVarPos[0] != copyVarPos[1]);
    if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns()) {
        auto mergedSegment = merge();
        return mergedSegment->unique();
    } else {
        return std::shared_ptr<const TGSegment>(new CompositeTGSegment(g,
                    nodes,
                    copyVarPos,
                    f_isSorted,
                    sortedField,
                    provenanceType,
                    sortBeforeAccess,
                    true,
                    replaceOffsets));
    }
}

void CompositeTGSegment::projectTo(const std::vector<int> &posFields,
        std::vector<std::shared_ptr<Column>> &out) const {
    LOG(ERRORL) << "projectTo() NOT IMPLEMENTED";
    throw 10;
}

std::shared_ptr<TGSegment> CompositeTGSegment::swap() const {
    if (copyVarPos.size() != 2) {
        LOG(ERRORL) << "Not implemented";
        throw 10;
    }
    std::vector<int> swappedPos;
    swappedPos.push_back(copyVarPos[1]);
    swappedPos.push_back(copyVarPos[0]);
    return std::shared_ptr<TGSegment>(new CompositeTGSegment(g, nodes,
                swappedPos, false, 0, provenanceType,
                sortBeforeAccess,
                removeDuplBeforeAccess,
                replaceOffsets));
}

std::vector<std::shared_ptr<const TGSegment>>
CompositeTGSegment::sliceByNodes(size_t startNodeIdx,
        std::vector<size_t> &provNodes) const {
    std::vector<std::shared_ptr<const TGSegment>> out;
    bool needMerging = false;
    for(auto n : nodes) {
        if (copyVarPos.size() != g.getNodeData(n)->getNColumns()) {
            needMerging = true;
        }
    }

    if (needMerging) //not supported yet
        throw 10;

    if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns()) {
        throw 10;
    } else {
        for (auto n : nodes) {
            std::vector<size_t> listNodes;
            listNodes.push_back(n);
            std::shared_ptr<const TGSegment> t =
                g.mergeNodes(listNodes, copyVarPos, replaceOffsets);
            //I cannot copy n in provNodes because it could be temporary
            auto nid = t->getNodeId();
            if (nid != n && !g.isTmpNode(n)) {
                //It can be that GBGraph has reshaped t, thus t->getNodeId()
                //is no longer equal to 'n'. In this case, I use the original n
                nid = n;
            }
            assert(nid != ~0ul && !g.isTmpNode(nid));
            provNodes.push_back(nid);
            out.push_back(t->slice(startNodeIdx++, 0, t->getNRows()));
        }
    }
    return out;
}

void CompositeTGSegment::appendTo(uint8_t colPos,
        std::vector<Term_t> &out) const {
    int posToAppend = copyVarPos[colPos];
    for(auto n : nodes) {
        g.getNodeData(n)->appendTo(posToAppend, out);
    }
}

void CompositeTGSegment::appendTo(uint8_t colPos,
        std::vector<std::pair<Term_t, Term_t>> &out) const {
    int posToAppend = copyVarPos[colPos];
    for(auto n : nodes) {
        g.getNodeData(n)->appendTo(posToAppend, out);
    }
}

void CompositeTGSegment::appendTo(uint8_t colPos1, uint8_t colPos2,
        std::vector<BinWithProv> &out) const {
    assert(!sortBeforeAccess);
    assert(!removeDuplBeforeAccess);
    int p1 = copyVarPos[colPos1];
    int p2 = copyVarPos[colPos2];
    for(auto n : nodes) {
        g.getNodeData(n)->appendTo(p1, p2, out);
    }
}
