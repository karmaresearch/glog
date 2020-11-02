#ifndef _GB_SEGMENT_INSERTER_H
#define _GB_SEGMENT_INSERTER_H

#include <vlog/gbchase/gbsegment.h>
#include <vlog/gbchase/gblegacysegment.h>

#include <vlog/term.h>
#include <vlog/column.h>
#include <vlog/segment.h>

#include <cstddef>
#include <memory>

#define THRESHOLD_CHECK_DUPLICATES 32*1000000
#include <google/dense_hash_set>

typedef google::dense_hash_set<Term_t> GBSegmentInserterEntities;

class GBSegmentInserter {
    private:
        std::vector<BuiltinFunction> fns;

    protected:
        virtual void addRow(Term_t *row) = 0;

    public:
        static std::unique_ptr<GBSegmentInserter> getInserter(size_t card,
                size_t nodeColumns, //n. columns which contain nodes
                bool delDupl);

        void addBuiltinFunctions(std::vector<BuiltinFunction> &fns) {
            this->fns = fns;
        }

        size_t getNBuiltinFunctions() const {
            return fns.size();
        }

        virtual void add(Term_t *row);

        virtual void swap(std::vector<Term_t> &t) {
            throw 10;
        }

        virtual void swap(std::vector<std::pair<Term_t, Term_t>> &t) {
            throw 10;
        }

        virtual GBSegmentInserterEntities getEntitiesAddedSoFar(int pos) {
            throw 10;
        }

        virtual bool isEmpty() const = 0;

        virtual size_t getNRows() const = 0;

        virtual std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance,
                bool lastColumnIsNode = false) = 0;

        virtual void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) = 0;

        virtual ~GBSegmentInserter() {}
};

template<class K, typename V, typename H=std::hash<V>>
class GBSegmentInserterImpl : public GBSegmentInserter {
    protected:
        const size_t card;

        bool shouldRemoveDuplicates;
        size_t checkDuplicatesAfter;
        bool useDuplicateMap;
        google::dense_hash_set<V, H> novelTuples;
        virtual bool isInMap(Term_t *row) = 0;
        virtual void populateMap() = 0;

        size_t processedRecords;

        K tuples;

    public:
        GBSegmentInserterImpl(size_t card,
                bool shouldRemoveDuplicates) : card(card),
        shouldRemoveDuplicates(shouldRemoveDuplicates),
        checkDuplicatesAfter(THRESHOLD_CHECK_DUPLICATES) {
            processedRecords = 0;
            useDuplicateMap = false;
        }

        bool isEmpty() const {
            return tuples.empty();
        }

        void add(Term_t *row) {
            if (shouldRemoveDuplicates) {
                processedRecords++;
                if (processedRecords % 10000000 == 0)
                    LOG(DEBUGL) << "Processed records: " << processedRecords;
                if (useDuplicateMap && isInMap(row)) {
                    return;
                }

                if (tuples.size() > checkDuplicatesAfter) {
                    //Sort and remove duplicates
                    std::sort(tuples.begin(), tuples.end());
                    auto e = std::unique(tuples.begin(), tuples.end());
                    //Are there duplicates?
                    LOG(DEBUGL) << "Removed tuples " <<
                        std::distance(tuples.end(), e);
                    if (e != tuples.end()) {
                        size_t oldSize = tuples.size();
                        tuples.erase(e, tuples.end());
                        size_t newSize = tuples.size();
                        size_t diff = oldSize - newSize;
                        LOG(DEBUGL) << "Diff: " << diff << " " <<
                            0.8 * THRESHOLD_CHECK_DUPLICATES;
                        if (!useDuplicateMap &&
                                diff > 0.8 * THRESHOLD_CHECK_DUPLICATES) {
                            LOG(DEBUGL) << "Decided to use the map with " <<
                                tuples.size() << " elements to filter out "
                                "duplicates immediately";
                            useDuplicateMap = true;
                            populateMap();
                        }
                        checkDuplicatesAfter = tuples.size() +
                            THRESHOLD_CHECK_DUPLICATES;
                    } else {
                        shouldRemoveDuplicates = false;
                        LOG(DEBUGL) << "No longer check!";
                    }
                }
            }
            GBSegmentInserter::add(row);
        }

        void swap(K &t) {
            tuples.swap(t);
        }

        virtual size_t getNRows() const {
            return tuples.size();
        }

        void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) {
            LOG(ERRORL) << "Should not occur";
            throw 10;
        }
};

class GBSegmentInserterUnary :
    public GBSegmentInserterImpl<std::vector<Term_t>,Term_t>
{
    protected:
        void addRow(Term_t *row) {
            tuples.push_back(row[0]);
        }

        bool isInMap(Term_t *row) {
            return novelTuples.count(row[0]);
        }

        void populateMap() {
            for(auto &v : tuples) {
                novelTuples.insert(v);
            }
        }

    public:
        GBSegmentInserterUnary(bool delDupl) :
            GBSegmentInserterImpl(1, delDupl)
    {
        novelTuples.set_empty_key(~0ul);
    }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance,
                bool lastColumnIsNode = false) {
            if (trackProvenance) {
                assert(lastColumnIsNode == false); //Otherwise
                //the arity is zero
                return std::shared_ptr<const TGSegment>(
                        new UnaryWithConstProvTGSegment(tuples, nodeId, isSorted,
                            sortedField));

            } else {
                return std::shared_ptr<const TGSegment>(
                        new UnaryTGSegment(tuples, nodeId, isSorted,
                            sortedField));
            }
        }
};

class GBSegmentInserterBinary : public GBSegmentInserterImpl<
                                std::vector<std::pair<Term_t, Term_t>>,
                                std::pair<Term_t, Term_t>>
{
    private:
        bool secondFieldConstant;

    protected:
        void addRow(Term_t *row) {
            if (secondFieldConstant && !tuples.empty() &&
                    tuples.back().second != row[1]) {
                secondFieldConstant = false;
            }
            tuples.push_back(std::make_pair(row[0], row[1]));
        }

        GBSegmentInserterEntities getEntitiesAddedSoFar(int pos) {
            GBSegmentInserterEntities e;
            e.set_empty_key(~0ul);
            if (pos == 0) {
                for(auto t : tuples)
                    e.insert(t.first);
            } else if (pos == 1) {
                for(auto t : tuples)
                    e.insert(t.second);
            } else {
                throw 10;
            }
            return e;
        }

        bool isInMap(Term_t *row) {
            return novelTuples.count(std::make_pair(row[0], row[1]));
        }

        void populateMap() {
            for(auto &v : tuples) {
                novelTuples.insert(v);
            }
        }

    public:
        GBSegmentInserterBinary(bool delDupl) : GBSegmentInserterImpl(1, delDupl),
        secondFieldConstant(true) {
            novelTuples.set_empty_key(std::make_pair(~0ul, ~0ul));
        }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance,
                bool lastColumnIsNode = false) {
            if (trackProvenance) {
                //lsatColumnIsNode indicates whether the second column contains
                //nodes instead of terms
                if (lastColumnIsNode) {
                    if (secondFieldConstant) {
                        std::vector<Term_t> t;
                        for(auto tuple : tuples) {
                            t.push_back(tuple.first);
                        }
                        return std::shared_ptr<const TGSegment>(
                                new UnaryWithConstProvTGSegment(
                                    t, tuples.back().second, isSorted,
                                    sortedField));
                    } else {
                        return std::shared_ptr<const TGSegment>(
                                new UnaryWithProvTGSegment(tuples, ~0ul,
                                    isSorted,
                                    sortedField));
                    }
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new BinaryWithConstProvTGSegment(tuples, nodeId,
                                isSorted,
                                sortedField));
                }
            } else {
                return std::shared_ptr<const TGSegment>(
                        new BinaryTGSegment(tuples, nodeId,
                            isSorted, sortedField));
            }
        }
};

//Use only for inserts where duplicates should be removed
struct BinWithDoubleProv {
    size_t first, second, node1, node2;

    bool operator <(const BinWithDoubleProv &rhs) const {
        return first < rhs.first ||
            (first == rhs.first && second < rhs.second);
    }

    bool operator==(const BinWithDoubleProv& rhs) const {
        return first == rhs.first && second == rhs.second;
    }
};

class GBSegmentInserterBinaryWithDoubleProv : public GBSegmentInserterImpl<
                                              std::vector<BinWithDoubleProv>,
                                              std::pair<Term_t, Term_t>>
{
    private:
        bool node1Constant;
        bool node2Constant;

    protected:
        void addRow(Term_t *row) {
            if (node1Constant && !tuples.empty() &&
                    tuples.back().node1 != row[2]) {
                node1Constant = false;
            }
            if (node2Constant && !tuples.empty() &&
                    tuples.back().node2 != row[3]) {
                node2Constant = false;
            }

            BinWithDoubleProv p;
            p.first = row[0];
            p.second = row[1];
            p.node1 = row[2];
            p.node2 = row[3];
            tuples.push_back(p);
        }

        bool isInMap(Term_t *row) {
            return novelTuples.count(std::make_pair(row[0], row[1]));
        }

        GBSegmentInserterEntities getEntitiesAddedSoFar(int pos) {
            GBSegmentInserterEntities e;
            e.set_empty_key(~0ul);
            if (pos == 0) {
                for(auto t : tuples)
                    e.insert(t.first);
            } else if (pos == 1) {
                for(auto t : tuples)
                    e.insert(t.second);
            } else {
                throw 10;
            }
            return e;
        }

        void populateMap() {
            novelTuples.set_empty_key(std::make_pair(~0ul, ~0ul));
            for(auto &v : tuples) {
                novelTuples.insert(std::make_pair(v.first, v.second));
            }
        }

    public:
        GBSegmentInserterBinaryWithDoubleProv(bool delDupl) :
            GBSegmentInserterImpl(2, delDupl),
            node1Constant(true),
            node2Constant(true) {
                if (!delDupl) {
                    throw 10;
                }
                novelTuples.set_empty_key(std::make_pair(~0ul, ~0ul));
            }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance,
                bool lastColumnIsNode = false) {
            if (!trackProvenance || !lastColumnIsNode) {
                //lsatColumnIsNode indicates whether the last column contains
                //nodes instead of terms. With this data structure it should
                //always be equal to true
                throw 10;
            }
            std::vector<BinWithProv> out(tuples.size());
            for(size_t i = 0; i < tuples.size(); ++i) {
                out[i].first = tuples[i].first;
                out[i].second = tuples[i].second;
                out[i].node = i;
            }
            return std::shared_ptr<const TGSegment>(
                    new BinaryWithProvTGSegment(out, ~0ul,
                        isSorted, sortedField));

        }

        void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) {
            if (node1Constant) {
                intermediateResultsNodes.push_back(
                        std::shared_ptr<Column>(
                            new CompressedColumn(tuples.back().node1,
                                tuples.size())));
            } else {
                ColumnWriter w;
                for(auto &v : tuples) {
                    w.add(v.node1);
                }
                intermediateResultsNodes.push_back(
                        w.getColumn());
            }

            if (node2Constant) {
                intermediateResultsNodes.push_back(
                        std::shared_ptr<Column>(
                            new CompressedColumn(tuples.back().node2,
                                tuples.size())));
            } else {
                ColumnWriter w;
                for(auto &v : tuples) {
                    w.add(v.node2);
                }
                intermediateResultsNodes.push_back(
                        w.getColumn());
            }

        }
};

class GBSegmentInserterNAry : public GBSegmentInserter
{
    private:
        std::vector<ColumnWriter> writers;
        std::vector<std::shared_ptr<Column>> columns;
        const size_t card;
        size_t addedRows;
        bool isFinal;

    protected:
        void addRow(Term_t *row) {
            for(size_t i = 0; i < card; ++i) {
                writers[i].add(row[i]);
            }
            addedRows++;
        }

    public:
        GBSegmentInserterNAry(size_t card) :
            writers(card),
            card(card),
            addedRows(0), isFinal(false) { }

        bool isEmpty() const  {
            return addedRows == 0;
        }

        size_t getNRows() const {
            return addedRows;
        }

        std::shared_ptr<const TGSegment> getSegment(size_t nodeId,
                bool isSorted,
                uint8_t sortedField,
                bool trackProvenance,
                bool lastColumnIsNode = false) {

            if (!isFinal) {
                for(int i = 0; i < card; ++i) {
                    columns.push_back(writers[i].getColumn());
                }
                isFinal = true;
            }

            auto ncols = columns.size();
            if (ncols == 2) {
                //In this case, we had the original vector with three
                //columns. One for terms, the other two for nodes. Then,
                //the last two columns have been replaced. Thus,
                //lastColumnIsNode == true and trackProvenance == true
                assert(lastColumnIsNode == true);
                auto &col1 = columns[0]->getVectorRef();
                size_t nrows = col1.size();
                if (trackProvenance) {
                    bool constantNodeVal = columns[1]->isConstant();
                    if (constantNodeVal) {
                        std::vector<Term_t> tuples(col1);
                        return std::shared_ptr<const TGSegment>(
                                new UnaryWithConstProvTGSegment(
                                    tuples, nodeId, isSorted, sortedField));
                    } else {
                        auto itrProv = columns[1]->getReader();
                        std::vector<std::pair<Term_t, Term_t>> out(nrows);
                        for(size_t i = 0; i < nrows; ++i) {
                            out[i].first = col1[i];
                            if (!itrProv->hasNext()) {
                                throw 10;
                            }
                            out[i].second = itrProv->next();
                        }
                        return std::shared_ptr<const TGSegment>(
                                new UnaryWithProvTGSegment(out, ~0ul,
                                    isSorted, sortedField));
                    }
                } else {
                    //This case should not happen
                    throw 10;
                }
            } else if (ncols == 3 && trackProvenance && lastColumnIsNode) {
                //Another special case. We had two columns with terms
                //and two columns with nodes. The last two were replaced by
                //a single one
                auto &col1 = columns[0]->getVectorRef();
                auto &col2 = columns[1]->getVectorRef();
                bool constantNodeVal = columns[2]->isConstant();
                size_t nrows = col1.size();
                if (constantNodeVal) {
                    std::vector<std::pair<Term_t, Term_t>> out(nrows);
                    for(size_t i = 0; i < nrows; ++i) {
                        out[i].first = col1[i];
                        out[i].second = col2[i];
                    }
                    return std::shared_ptr<const TGSegment>(
                            new BinaryWithConstProvTGSegment(
                                out, columns[2]->first(), isSorted, sortedField));
                } else {
                    std::vector<BinWithProv> out(nrows);
                    auto itrNode = columns[2]->getReader();
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
                if (trackProvenance) {
                    if (lastColumnIsNode) {
                        return std::shared_ptr<const TGSegment>(
                                new TGSegmentLegacy(columns, addedRows, isSorted,
                                    sortedField, true));
                    } else {
                        //In this case, I must add a special column with a
                        //constant node because all the columns in ``columns''
                        //contain terms
                        auto newcolumns = columns;
                        newcolumns.push_back(std::shared_ptr<Column>(
                                    new CompressedColumn(nodeId, addedRows)));
                        return std::shared_ptr<const TGSegment>(
                                new TGSegmentLegacy(newcolumns, addedRows,
                                    isSorted, sortedField, true));
                    }
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new TGSegmentLegacy(columns, addedRows, isSorted,
                                sortedField, false));
                }
            }
        }

        void postprocessJoin(std::vector<std::shared_ptr<Column>>
                &intermediateResultsNodes) {
            //Take out the last two columns
            assert(card > 2);
            for(int i = 0; i < card - 2; ++i) {
                columns.push_back(writers[i].getColumn());
            }
            //For now, always add one extra column
            CompressedColumnBlock b(0, 1, addedRows);
            std::vector<CompressedColumnBlock> blocks;
            blocks.push_back(b);
            columns.push_back(std::shared_ptr<Column>(
                        new CompressedColumn(blocks, addedRows)));

            //Save the columns with the mappings to the nodes
            auto col1 = writers[card - 2].getColumn();
            auto col2 = writers[card - 1].getColumn();
            intermediateResultsNodes.push_back(col1);
            intermediateResultsNodes.push_back(col2);
            isFinal = true;
        }
};

#endif
