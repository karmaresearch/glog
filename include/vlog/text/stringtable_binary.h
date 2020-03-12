#ifndef _STRING_TABLE_BINARY_H
#define _STRING_TABLE_BINARY_H

#include <vlog/text/stringtable.h>

class StringTableBinary : public StringTable {
    private:
        std::unique_ptr<char[]> buffer1;
        std::unique_ptr<char[]> buffer2;

    protected:
        bool execFunction(const uint64_t t1, const uint64_t t2);

    public:
        StringTableBinary(PredId_t predid,
                EDBLayer *layer,
                std::string fname);

        uint8_t getArity() const {
            return 2;
        }
};

#endif
