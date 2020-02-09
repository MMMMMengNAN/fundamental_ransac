#include "fundamental.ransac.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
int main()
{
    //加载归一化之后的匹配对
    correspondances2D2D corr_all;
    std::ifstream in(".\\correspondences.txt");
    
    std::string line,word;
    int n_line = 0;
    while(std::getline(in,line))
    {
        std::stringstream stream(line);
        if(n_line == 0)
        {
            int n_corrs = 0;
            stream>>n_corrs;
            corr_all.resize(n_corrs);
            n_line++;
            continue;
        }
        if(n_line > 0)
        {
            stream>>corr_all[n_line - 1].p1[0] >> corr_all[n_line - 1].p1[1]
                    >>corr_all[n_line - 1].p2[0] >> corr_all[n_line - 1].p2[1];
        }
        n_line++;
    }

    std::cout<<corr_all.size()<<std::endl;
    FundamentalMatrix F = calc_fundamental_ransac(corr_all);
    std::cin.get();
}