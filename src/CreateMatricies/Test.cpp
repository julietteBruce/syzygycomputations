
#include "Combinatorics.h"
#include <iostream>

using namespace std;

int main(){
  auto a=createProductIntegerVectors(1,1,4,5);
  for(int i=0;i<a.size();i++){
    for(int j=0;j<a[i].size();j++){
      cout << a[i][j] << " ";
    }
    cout << endl;
  }
  return 0;
}

