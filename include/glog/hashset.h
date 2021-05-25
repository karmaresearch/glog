#ifndef hashset_h
#define hashset_h

#include <vlog/concepts.h>
#include <google/dense_hash_set>

#include <inttypes.h>

#define HASHSET_SIZEBLOCK 1*1000*1000
class HashSet {
private:
    struct Hasher {
        size_t card;
        Hasher(size_t card) : card(card) {}
        size_t operator() (Term_t *x) const {
            if (x == NULL)
                return 0;
            
            size_t result = 0;
            for (size_t i = 0; i < card; i++) {
                result = (result + (324723947 + x[i])) ^93485734985;
            }
            return (size_t) result;
        }
    };
    
    struct EqualChecker {
        size_t card;
        EqualChecker(size_t card) : card(card) {}
        bool operator()(Term_t* s1, Term_t* s2) const {
            if (s1 == s2) {
                return true;
            } else if (s1 == NULL && s2 != NULL)
                return false;
            else if (s1 != NULL && s2 == NULL)
                return false;
            else {
                for(size_t i = 0; i < card; ++i)
                    if (s1[i] != s2[i])
                        return false;
                return true;
            }
        }
    };
    
    const size_t card;
    Hasher h;
    EqualChecker e;
    google::dense_hash_set<Term_t *, Hasher, EqualChecker> elements;
    std::vector<std::unique_ptr<Term_t[]>> blocks;
    uint32_t blockCounter;
    Term_t *currentblock;
    
public:
    HashSet(size_t card, size_t estSize) : card(card), h(card), e(card), elements(estSize, h, e) {
        elements.set_empty_key(NULL);
        blockCounter = 0;
        currentblock = NULL;
    }
    
    bool isIn(Term_t *term) {
        return elements.count(term);
    }
    
    void add(Term_t *term) {
        if (!currentblock || blockCounter >= HASHSET_SIZEBLOCK) {
            //Create a new block
            std::unique_ptr<Term_t[]> n =
                std::unique_ptr<Term_t[]>(new Term_t[card * HASHSET_SIZEBLOCK]);
            currentblock = n.get();
            blocks.push_back(std::move(n));
            blockCounter = 0;
        }
        for(uint8_t i = 0; i < card; ++i) {
            currentblock[i] = term[i];
        }
        elements.insert(currentblock);
        currentblock += card;
        blockCounter++;
    }
};

#endif
