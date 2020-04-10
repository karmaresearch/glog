#ifndef _STRING_TABLE_BINARY_H
#define _STRING_TABLE_BINARY_H

#include <vlog/text/stringtable.h>

class StringTableBinary : public StringTable {
    private:
        std::unique_ptr<char[]> buffer1;
        std::unique_ptr<char[]> buffer2;
        int param1;

    protected:
        bool execFunction(const uint64_t t1, const uint64_t t2);

        bool builtinFunction(Term_t *t, uint8_t *pos);

    public:
        StringTableBinary(PredId_t predid,
                EDBLayer *layer,
                std::string fname,
                std::string param1);

        uint8_t getArity() const {
            return 2;
        }

        BuiltinFunction getBuiltinFunction() {
            BuiltinFunction fn;
            fn.fn = std::bind(&StringTableBinary::builtinFunction,
                    this,
                    std::placeholders::_1, std::placeholders::_2);
            return fn;

        }
};

#endif
