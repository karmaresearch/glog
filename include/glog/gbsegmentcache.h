#ifndef _TGSEGMENTCACHE_H
#define _TGSEGMENTCACHE_H

#include <vlog/concepts.h>

#include <vlog/gbchase/gbsegmentitr.h>

#include <map>
#include <vector>

class CacheEntry {
    private:
        std::vector<size_t> nodes;

    public:
        CacheEntry(const std::vector<size_t> &n) : nodes(n) {}

        bool operator <(const CacheEntry& rhs) const {
            if (nodes.size() != rhs.nodes.size()) {
                return nodes.size() < rhs.nodes.size();
            } else {
                for(int i = 0; i < nodes.size(); ++i) {
                    if (nodes[i] != rhs.nodes[i]) {
                        return nodes[i] < rhs.nodes[i];
                    }
                }
            }
            return false;
        }
};

class TGSegment;
class SegmentCache {
    private:
        static SegmentCache instance;

        std::map<CacheEntry, std::shared_ptr<const TGSegment>> cacheVar0;
        std::map<CacheEntry, std::shared_ptr<const TGSegment>> cacheVar1;

    public:
        static SegmentCache &getInstance() {
            return instance;
        }

        bool contains(const std::vector<size_t> &key, const std::vector<uint8_t> &fields) const {
            if (fields.size() != 1)
                return false;

            auto field = fields[0];
            CacheEntry k(key);
            if (field == 0)
                return cacheVar0.count(k);
            else if (field == 1)
                return cacheVar1.count(k);
            else {
                LOG(INFOL) << "Not supported";
                throw 10;
            }
        }

        void insert(const std::vector<size_t> &key, const std::vector<uint8_t> &fields,
                std::shared_ptr<const TGSegment> value) {
            assert(fields.size() == 1);
            auto field = fields[0];
            CacheEntry k(key);
            if (field == 0) {
                cacheVar0[k] = value;
            } else {
                assert(field == 1);
                cacheVar1[k] = value;
            }
        }

        const std::shared_ptr<const TGSegment> get(
                const std::vector<size_t> &key,
                const std::vector<uint8_t> &fields) {
            assert(fields.size() == 1);
            auto field = fields[0];
            CacheEntry k(key);
            if (field == 0) {
                return cacheVar0[k];
            } else {
                assert(field == 1);
                return cacheVar1[k];
            }
        }

        void clear() {
            cacheVar0.clear();
            cacheVar1.clear();
        }
};

#endif
