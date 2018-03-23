#include <iostream>
#include "../src/incomplete_lq_builder.h"

int main(int argc, char **argv) {
  if(argc!=1) {
    std::cerr << "Parameters must be mesh file name\n";
    return 1;
  }
  return 0;
}
