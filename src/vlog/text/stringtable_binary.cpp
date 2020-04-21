#include <vlog/text/stringtable_binary.h>

#include <vlog/edb.h>

#include <cstring>
#include <chrono>
#include <algorithm>
#include <string>

StringTableBinary::StringTableBinary(PredId_t predid,
        EDBLayer *layer,
        std::string fname,
        std::string param1) : StringTable(predid, layer, fname) {
    buffer1 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
    buffer2 = std::unique_ptr<char[]>(new char[MAX_TERM_SIZE]);
    this->param1 = 0;
    if (param1 != "") {
        this->param1 = stoi(param1);
    }
}

//Function copied from:
//https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/
//Levenshtein_distance#C++
//
unsigned int leven_distance(const std::string& s1, const std::string& s2)
{
    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<std::vector<unsigned int>> d(len1 + 1,
            std::vector<unsigned int>(len2 + 1));

    d[0][0] = 0;
    for(unsigned int i = 1; i <= len1; ++i) d[i][0] = i;
    for(unsigned int i = 1; i <= len2; ++i) d[0][i] = i;

    for(unsigned int i = 1; i <= len1; ++i)
        for(unsigned int j = 1; j <= len2; ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1,
                    d[i][j - 1] + 1,
                    d[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ?
                            0 : 1) });
    return d[len1][len2];
}

bool StringTableBinary::execFunction(const uint64_t t1, const uint64_t t2) {
    //std::chrono::steady_clock::time_point begin =
    //std::chrono::steady_clock::now();

    bool found1 = layer->getDictText(t1, buffer1.get());
    bool found2 = layer->getDictText(t2, buffer2.get());
    if (!found1 || !found2) {
        return false;
    }

    bool outcome = false;
    if (fname == "containedIn") {
        auto resp = std::string(buffer2.get()).find(std::string(buffer1.get()));
        outcome = resp != std::string::npos;
    } else if (fname == "equal") {
        auto s1 = std::string(buffer1.get());
        auto s2 = std::string(buffer2.get());
        outcome = s1 == s2;
        //outcome = strcmp(buffer1.get(), buffer2.get()) == 0;
    } else if (fname == "levenshtein") {
        auto s1 = std::string(buffer1.get());
        std::transform(s1.begin(), s1.end(),s1.begin(), ::toupper);
        auto s2 = std::string(buffer2.get());
        std::transform(s2.begin(), s2.end(),s2.begin(), ::toupper);
        auto dis = leven_distance(s1, s2);
        outcome = dis <= param1;
    } else {
        LOG(ERRORL) << "(StringTableBinary) Function " << fname << " is unknown";
        throw 10;
    }

    //std::chrono::steady_clock::time_point end =
    //std::chrono::steady_clock::now();
    //auto micros = std::chrono::duration_cast<std::chrono::nanoseconds>(end
    //- begin).count();
    //std::cout << "Runtime (nanos): " << micros << std::endl;

    return outcome;
}

bool StringTableBinary::builtinFunction(Term_t *t, uint8_t *pos) {
    return execFunction(t[pos[0]], t[pos[1]]);
}
