#ifndef _GB_SEGMENT_INSERTER_H
#define _GB_SEGMENT_INSERTER_H

#include <vlog/gbchase/gbsegment.h>

#include <vlog/term.h>
#include <vlog/column.h>
#include <vlog/segment.h>

#include <cstddef>
#include <memory>

class GBSegmentInserter {
    private:
        std::vector<BuiltinFunction> fns;

    protected:
        virtual void addRow(Term_t *row) = 0;

    public:
        static std::unique_ptr<GBSegmentInserter> getInserter(size_t card);

        void addBuiltinFunctions(std::vector<BuiltinFunction> &fns) {
            this->fns = fns;
        }

        size_t getNBuiltinFunctions() const {
            return fns.size();
        }

        void add(Term_t *row);

        virtual void swap(std::vector<Term_t> &t) {
            throw 10;
        }

        virtual void swap(std::vector<std::pair<Term_t, Term_t>> &t) {
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

template<class K>
class GBSegmentInserterImpl : public GBSegmentInserter {
    protected:
        const size_t card;
        K tuples;

    public:
        GBSegmentInserterImpl(size_t card) : card(card) {
        }

        bool isEmpty() const {
            return tuples.empty();
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

class GBSegmentInserterUnary : public GBSegmentInserterImpl<std::vector<Term_t>>
{
    public:
        GBSegmentInserterUnary() : GBSegmentInserterImpl(1) { }

        void addRow(Term_t *row) {
            tuples.push_back(row[0]);
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
                                std::vector<std::pair<Term_t, Term_t>>>
{
    private:
        bool secondFieldConstant;

    public:
        GBSegmentInserterBinary() : GBSegmentInserterImpl(1),
        secondFieldConstant(true) { }

        void addRow(Term_t *row) {
            if (secondFieldConstant && !tuples.empty() &&
                    tuples.back().second != row[1]) {
                secondFieldConstant = false;
            }
            tuples.push_back(std::make_pair(row[0], row[1]));
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
                                    t, nodeId, isSorted, sortedField));
                    } else {
                        return std::shared_ptr<const TGSegment>(
                                new UnaryWithProvTGSegment(tuples, ~0ul,
                                    isSorted,
                                    sortedField));
                    }
                } else {
                    return std::shared_ptr<const TGSegment>(
                            new BinaryWithConstProvTGSegment(tuples, ~0ul,
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

class GBSegmentInserterNAry : public GBSegmentInserter
{
    private:
        std::vector<ColumnWriter> writers;
        std::vector<std::shared_ptr<Column>> columns;
        const size_t card;
        size_t addedRows;
        bool isFinal;

    public:
        GBSegmentInserterNAry(size_t card) : writers(card),
        card(card),
        addedRows(0), isFinal(false) { }

        bool isEmpty() const  {
            return addedRows == 0;
        }

        size_t getNRows() const {
            return addedRows;
        }

        void addRow(Term_t *row) {
            for(size_t i = 0; i < card; ++i) {
                writers[i].add(row[i]);
            }
            addedRows++;
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
                    /*auto &col2 = columns[1]->getVectorRef();
                      std::vector<std::pair<Term_t, Term_t>> out(nrows);
                      for(size_t i = 0; i < nrows; ++i) {
                      out[i].first = col1[i];
                      out[i].second = col2[i];
                      }
                      return std::shared_ptr<const TGSegment>(
                      new BinaryTGSegment(out, nodeId, isSorted,
                      sortedField));
                      */
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
                        columns.push_back(std::shared_ptr<Column>(
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
