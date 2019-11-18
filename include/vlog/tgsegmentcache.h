#ifndef _TGSEGMENTCACHE_H
#define _TGSEGMENTCACHE_H

#include <vlog/concepts.h>
#include <vlog/tgsegmentitr.h>

#include <map>
#include <vector>

class TGSegment;
class SegmentCache {
    private:
        static SegmentCache instance;

        std::map<long, std::shared_ptr<const TGSegment>> cacheVar0;
        std::map<long, std::shared_ptr<const TGSegment>> cacheVar1;

    public:
        static SegmentCache &getInstance() {
            return instance;
        }

        static long hash(std::vector<size_t> &nodes);

        bool contains(long key, const uint8_t field) {
            if (field == 0)
                return cacheVar0.count(key);
            else if (field == 1)
                return cacheVar1.count(key);
            else {
                LOG(INFOL) << "Not supported";
                throw 10;
            }
        }

        void insert(long key, const uint8_t field,
                std::shared_ptr<const TGSegment> value) {
            if (field == 0) {
                cacheVar0[key] = value;
            } else {
                assert(field == 1);
                cacheVar1[key] = value;
            }
        }

        const std::shared_ptr<const TGSegment> get(long key,
                const uint8_t field) {
            if (field == 0) {
                return cacheVar0[key];
            } else {
                assert(field == 1);
                return cacheVar1[key];
            }
        }

        void clear() {
            cacheVar0.clear();
            cacheVar1.clear();
        }
};

#endif
