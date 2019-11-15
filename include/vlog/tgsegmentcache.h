#ifndef _TGSEGMENTCACHE_H
#define _TGSEGMENTCACHE_H

#include <vlog/concepts.h>
#include <vlog/tgsegmentitr.h>

#include <map>
#include <vector>

class SegmentCache {
    private:
        static SegmentCache instance;

        std::map<const Term_t *, std::unique_ptr<std::vector<Term_t>>> cacheUnique;
        std::map<const std::pair<Term_t,Term_t> *, std::unique_ptr<std::vector<std::pair<Term_t, Term_t>>>> cacheDouble;
        std::map<const std::pair<Term_t,Term_t> *, std::unique_ptr<std::vector<std::pair<Term_t, Term_t>>>> cacheDoubleInv;
        std::map<const BinWithProv *, std::unique_ptr<std::vector<BinWithProv>>> cacheBinWithProv;
        std::map<const BinWithProv *, std::unique_ptr<std::vector<BinWithProv>>> cacheBinWithProvInv;

    public:
        static SegmentCache &getInstance() {
            return instance;
        }

        bool contains(const Term_t *key) {
            return cacheUnique.count(key);
        }

        bool contains(const std::pair<Term_t,Term_t> *key) {
            return cacheDouble.count(key);
        }

        bool contains(const std::pair<Term_t,Term_t> *key, const uint8_t field) {
            if (field == 0)
                return cacheDouble.count(key);
            else
                return cacheDoubleInv.count(key);
        }

        bool contains(const BinWithProv *key, const uint8_t field) {
            if (field == 0)
                return cacheBinWithProv.count(key);
            else
                return cacheBinWithProvInv.count(key);
        }

        void insert(const Term_t *key, std::vector<Term_t> &value) {
            cacheUnique[key] = std::unique_ptr<std::vector<Term_t>>(new std::vector<Term_t>(value));
        }

        void insert(const std::pair<Term_t, Term_t> *key, std::vector<
                std::pair<Term_t, Term_t>> &value) {
            cacheDouble[key] = std::unique_ptr<std::vector<std::pair<Term_t, Term_t>>>(
                    new std::vector<std::pair<Term_t, Term_t>>(value));
        }

        void insert(const std::pair<Term_t, Term_t> *key, const uint8_t field,
                std::vector<std::pair<Term_t, Term_t>> &value) {
            if (field == 0) {
                cacheDouble[key] = std::unique_ptr<std::vector<std::pair<Term_t, Term_t>>>(
                        new std::vector<std::pair<Term_t, Term_t>>(value));
            } else {
                cacheDoubleInv[key] = std::unique_ptr<std::vector<std::pair<Term_t, Term_t>>>(
                        new std::vector<std::pair<Term_t, Term_t>>(value));
            }
        }

        void insert(const BinWithProv *key, const uint8_t field,
                std::vector<BinWithProv> &value) {
            if (field == 0) {
                cacheBinWithProv[key] = std::unique_ptr<std::vector<BinWithProv>>(
                        new std::vector<BinWithProv>(value));
            } else {
                cacheBinWithProvInv[key] = std::unique_ptr<std::vector<BinWithProv>>(
                        new std::vector<BinWithProv>(value));
            }
        }

        const std::vector<Term_t> &get(const Term_t *key) {
            auto &v = cacheUnique[key];
            return *v.get();
        }

        const std::vector<std::pair<Term_t, Term_t>> &get(
                const std::pair<Term_t,Term_t> *key) {
            auto &v = cacheDouble[key];
            return *v.get();
        }

        const std::vector<std::pair<Term_t, Term_t>> &get(
                const std::pair<Term_t,Term_t> *key, const uint8_t field) {
            if (field == 0) {
                auto &v = cacheDouble[key];
                return *v.get();
            } else {
                auto &v = cacheDoubleInv[key];
                return *v.get();
            }
        }

        const std::vector<BinWithProv> &get(
                const BinWithProv *key, const uint8_t field) {
            if (field == 0) {
                auto &v = cacheBinWithProv[key];
                return *v.get();
            } else {
                auto &v = cacheBinWithProvInv[key];
                return *v.get();
            }
        }
};

#endif
