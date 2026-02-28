#include <iostream>
#include <string>
#include <random>
#include <array>

int main()
{
    constexpr int dim = 3;
    int L = 6;
    int lin_l = 3;
    std::array<size_t, dim + 1> plengths;
    
    plengths.front() = 1;
    std::fill(plengths.begin() + 1, plengths.end(), L);
    std::partial_sum(plengths.begin(), plengths.end(), plengths.begin(),
                        std::multiplies<>{});

    for (int i = 0; i < plengths.size(); i++)
        std::cout << plengths[i] << " ";
    std::cout << std::endl;
        
    //std::vector<int> ti = {0,2,4,12,14,16,24,26,28,72,74,76,84,86,88,
    //                    96,98,100,144,146,148,156,158,160,168,170,172 };
    std::vector<int> ti = {0,3,18,21,108,111,126,129};
    //std::vector<int> ti = {0,4,8,48,52,56,96,100,104,576,580,584,1152};
    for (int j = ti.size()-1; j >= 0; j--) {
        int uidx =ti[j];
        
        auto decr = [&](int i, auto&& decr) {
                if(i == 0) {uidx -= lin_l; return;}
                auto upper = plengths[i];
                auto lower = plengths[i-1];
                size_t shift = i-1 > 0 ? (lin_l-1)*lower : 0;
                if (uidx % (lin_l* upper) == 0) {
                    //std::cout << shift << "   ";
                    uidx += -(lin_l-1)*upper - lin_l - shift;
                    //uidx = (uidx/upper - lin_l-1) * upper;
                    return;
                }
                decr(i - 1, decr);
        };
        decr(plengths.size()-2,decr);
        /*
        auto incr = [&](int i, auto&& incr) {
            if(i == 0) {
                idx += lin_l; return;
            }  
            int sub = i-1 > 0 ? plengths[i-1] - lin_l*plengths[i-2] : 0;
            if ((idx + lin_l*plengths[i-1] - sub) % plengths[i] == 0) {
                idx += lin_l*plengths[i-1];
                idx = (idx/plengths[i] + lin_l-1) * plengths[i];
                return;
            }
            incr(i - 1, incr);
        };
        incr(plengths.size()-2,incr);
        */
        std::cout << ti[j] << " -> " << uidx << std::endl;
    }
    /*
    for (int j = 0; j < ti.size(); j++) {
        int idx = ti[j];
        if ((idx + lin_l*plengths[1] - (plengths[1] - lin_l*plengths[0]))
            % plengths[2] == 0) {
            idx += lin_l*plengths[1];
            idx = (idx/plengths[2] + lin_l-1)* plengths[2];
        }
        else if ((idx + lin_l*plengths[0]) % plengths[1] == 0) {
            idx += lin_l*plengths[0];
            idx = (idx/plengths[1] + lin_l-1)* plengths[1];
        }
        else {
            idx += lin_l;
        }
        std::cout << ti[j] << " -> " << idx << std::endl;
    }
    */
    return 0;
}
