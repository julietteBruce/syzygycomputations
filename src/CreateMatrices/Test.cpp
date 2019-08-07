#include <catch2/catch.hpp>
#include "ToricVariety.h"
#include "Basis.h"
#include <vector>
#include <algorithm>

using namespace std;

bool isNonDecreasing(const vector<int>& v, size_t start, size_t stop){
    for(size_t i = start+1;i<stop;i++){
        if(v[i]<v[i-1]){
            return false;
        }
    }
    return true;
}

int sum(const vector<int>& v){
    int ret = 0;
    for(int i : v)
        ret += i;
    return ret;
}

TEST_CASE("Pn multidegrees", "[Pn]"){
    Pn p(3);
    SECTION("multidegree count"){
        REQUIRE(p.multidegrees({2},false).size() == 10);
        REQUIRE(p.multidegrees({2},true).size() == 2);
        REQUIRE(p.multidegrees({3},false).size() == 20);
        REQUIRE(p.multidegrees({3},true).size() == 3);
    }
    SECTION("multidegree sizes"){
        auto m = p.multidegrees({3},false);
        REQUIRE(all_of(m.begin(),m.end(),[&p](vector<int> v){return v.size()==p.mdegSize;}));
    }
    SECTION("multidegree sums"){
        auto m = p.multidegrees({3},false);
        REQUIRE(all_of(m.begin(),m.end(),[&p](vector<int> v){return sum(v)==3;}));
    }
    SECTION("deduped multidegrees are non-decreasing"){
        auto m = p.multidegrees({5},true);
        REQUIRE(all_of(m.begin(),m.end(),
                       [](const vector<int>& v){return isNonDecreasing(v,0,v.size());}));
    }
}

TEST_CASE("Product multidegrees", "[Product]"){
    Product p(vector<ToricVariety*>({new Pn(3),new Pn(3)}));
    SECTION("multidegree count"){
        REQUIRE(p.multidegrees({2,2},false).size() == 100);
        REQUIRE(p.multidegrees({2,2},true).size() == 4);
        REQUIRE(p.multidegrees({3,2},false).size() == 200);
        REQUIRE(p.multidegrees({3,2},true).size() == 6);
    }
    SECTION("multidegree sizes"){
        auto m = p.multidegrees({3,3},false);
        REQUIRE(all_of(m.begin(),m.end(),[&p](vector<int> v){return v.size()==p.mdegSize;}));
        REQUIRE(p.mdegSize==8);
        REQUIRE(p.degreeSize==2);
    }
    SECTION("multidegree sums"){
        auto m = p.multidegrees({3,3},false);
        //This is an imperfect test
        REQUIRE(all_of(m.begin(),m.end(),[&p](vector<int> v){return sum(v)==6;}));
    }
    SECTION("deduped multidegrees are non-decreasing"){
        auto m = p.multidegrees({3,3},true);
        REQUIRE(all_of(m.begin(),m.begin(),
                       [](const vector<int>& v){
                           return isNonDecreasing(v,0,4) &&
                               isNonDecreasing(v,4,8);}));
    }
}

TEST_CASE("P1Bundle multidegrees",""){
    P1Bundle p(new Pn(1),{1});
    SECTION("multidegree count"){
        REQUIRE(p.multidegrees({2,2},false).size() == 6);
        REQUIRE(p.multidegrees({2,2},true).size() == 4);
        REQUIRE(p.multidegrees({3,2},false).size() == 9);
        REQUIRE(p.multidegrees({3,2},true).size() == 5);
    }
    SECTION("multidegree sizes"){
        auto m = p.multidegrees({4,2},false);
        REQUIRE(all_of(m.begin(),m.end(),[&p](vector<int> v){return v.size()==p.mdegSize;}));
        REQUIRE(p.mdegSize==4);
        REQUIRE(p.degreeSize==2);
    }
}

TEST_CASE("WedgeBasis","[WedgeBasis]"){
    Product tv({new Pn (3),new Pn (3)});
    SECTION("size"){
        WedgeBasis basis = createToricWedgeBasis(tv,{3,2},0);
        REQUIRE(basis.getMonomials()==tv.multidegrees({3,2},false));
        REQUIRE(basis.size()==1);
    }
    SECTION("rank/unrank"){
        WedgeBasis basis = createToricWedgeBasis(tv,{3,2},2);
        REQUIRE(basis.size()>=1000);
        long long r = GENERATE(0,take(100,random(0,1000)));
        REQUIRE(basis.rank(basis.unrank(r))==r);
    }
    SECTION("iter"){
        WedgeBasis basis = createToricWedgeBasis(tv,{3,2},2);
        auto iter = basis.getIter();
        long long i=0;
        do{
            REQUIRE(iter.getCurr()==basis.unrank(i));
            i++;
        }while(iter.next());
        REQUIRE(i==basis.size());
        iter.reset();
        i = 0;
        do{
            REQUIRE(iter.getCurr()==basis.unrank(i));
            i++;
        }while(iter.next());
    }
}
