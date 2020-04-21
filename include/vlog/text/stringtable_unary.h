#ifndef _STRING_TABLE_UNARY_H
#define _STRING_TABLE_UNARY_H

#include <vlog/text/stringtable.h>

class StringTableUnary : public StringTable {
    private:
        std::unique_ptr<char[]> buffer1;
        const std::string param;

        int param1_int;
        char param1_char;

    protected:
        bool execFunction(const uint64_t t1);

        bool builtinFunction(Term_t *t, uint8_t *pos) {
            return execFunction(t[pos[0]]);
        }

    public:
        StringTableUnary(PredId_t predid,
                EDBLayer *layer,
                std::string fname,
                std::string param = "");

        uint8_t getArity() const {
            return 1;
        }

        BuiltinFunction getBuiltinFunction() {
            BuiltinFunction fn;
            fn.fn = std::bind(&StringTableUnary::builtinFunction,
                    this,
                    std::placeholders::_1, std::placeholders::_2);
            return fn;
        }
};

#endif
