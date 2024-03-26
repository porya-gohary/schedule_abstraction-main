#ifndef NP_STATISTICS_HPP
#define NP_STATISTICS_HPP
#include <iostream>
#include <ostream>
#include <string>

namespace NP {
  struct StatCollect {
    int count[50];
    int total;
    int div;
    int max;
    int rettrue;
    int retfalse;
    bool printit;
    const char *label;

    StatCollect(const char *lab)
      : total(0), div(1), max(0),rettrue(0),retfalse(0),printit(false),label(lab)
    {
      for (int i=0; i<50; i++) count[i]=0;
    }

    void tick(int value)
    {
      if (value>max) {
	printit=true;
	max = value;
      }
      total++;
      if (total%500==0) printit=true;
      if (value<50) count[value]++;
    }
    void tick(bool retvalue)
    {
      if (retvalue)  rettrue++;
      else retfalse++;
    }
    void print()
    {
      if (printit) {
	std::cerr << label << " (" << total << ", " << max << ") : "<< count[0];
	for (int i=1; i<50 && i<=max; i++)
	  std::cerr << ", " << count[i];
	std::cerr << "\n";
	std::cerr << label << " returns true: " << rettrue << ", false: "<< retfalse << "\n";
      }
      printit=false;
    }
  };
}
#endif
