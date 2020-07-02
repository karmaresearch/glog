#include <vlog/gbchase/gbcompositesegment.h>

std::shared_ptr<const TGSegment> CompositeTGSegment::merge() const {
    std::chrono::system_clock::time_point start =
        std::chrono::system_clock::now();
    std::shared_ptr<const TGSegment> t = g.mergeNodes(nodes, copyVarPos);
    std::chrono::system_clock::time_point end =
        std::chrono::system_clock::now();
    std::chrono::duration<double, std::milli> dur = end - start;
    return t;
}

size_t CompositeTGSegment::getNRows() const {
    assert(nodes.size() > 0);
    assert(copyVarPos.size() == 1 || copyVarPos[0] != copyVarPos[1]);
    if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns()) {
        LOG(ERRORL) << "Case not supported";
        throw 10;
    }

    size_t nrows = 0;
    for(auto n : nodes)
        nrows += g.getNodeSize(n);
    return nrows;
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
    if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns()) {
        auto mergedSegment = merge();
        if (mergedSegment->isSorted())
            return mergedSegment;
        else
            return mergedSegment->sort();
    } else {
        return std::shared_ptr<const TGSegment>(new CompositeTGSegment(g,
                    nodes,
                    copyVarPos,
                    f_isSorted,
                    sortedField,
                    trackProvenance,
                    !isSorted()));
    }
}

std::shared_ptr<TGSegment> CompositeTGSegment::sortBy(
        std::vector<uint8_t> &fields) const {
    auto mergedSegment = merge();
    return mergedSegment->sortBy(fields);
}

std::shared_ptr<TGSegment> CompositeTGSegment::sortByProv(size_t ncols,
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
                trackProvenance,
                sortBeforeAccess,
                removeDuplBeforeAccess));
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
                    trackProvenance,
                    sortBeforeAccess,
                    true));
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
                swappedPos, false, 0, trackProvenance,
                sortBeforeAccess,
                removeDuplBeforeAccess));
}

std::vector<std::shared_ptr<const TGSegment>>
CompositeTGSegment::sliceByNodes(size_t startNodeIdx,
        std::vector<size_t> &provNodes) const {
    std::vector<std::shared_ptr<const TGSegment>> out;
    if (copyVarPos.size() != g.getNodeData(nodes[0])->getNColumns()) {
        throw 10;
    } else {
        for (auto n : nodes) {
            provNodes.push_back(n);
            std::vector<size_t> listNodes;
            listNodes.push_back(n);
            std::shared_ptr<const TGSegment> t =
                g.mergeNodes(listNodes, copyVarPos);
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
