#ifndef _STRING_TABLE_UNARY_H
#define _STRING_TABLE_UNARY_H

#include <vlog/text/stringtable.h>

class StringTableUnary : public StringTable {
    private:
        std::unique_ptr<char[]> buffer1;
        const std::string param;

    protected:
        bool execFunction(const uint64_t t1);

    public:
        StringTableUnary(PredId_t predid,
                EDBLayer *layer,
                std::string fname,
                std::string param = "");

        uint8_t getArity() const {
            return 1;
        }
};

#endif
