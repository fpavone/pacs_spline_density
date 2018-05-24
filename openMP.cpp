#include <iostream>
#include <vector>
#include <omp.h>


class OBJ
{

public:
  //OBJ(int nn): n(nn) {};
int n = 2;
  void foo(int & a, int & out){
    //std::cout << n*a << std::endl;
    std::cout << "n = " << n << " thread = " << omp_get_thread_num() << '\n';
    n += omp_get_thread_num();
    out = n*a;
  }
};

int main()
{
  std::vector<int> data{1,2,3,4,5,6};
  std::vector<int> out(6);
  OBJ dens;

  #pragma omp parallel private(dens.n)
  {
    #pragma omp for
      for(int i = 0; i < data.size(); i++){
        dens.foo(data[i],out[i]);
        //std::cout << "THREAD: " << omp_get_thread_num() << std::endl;
      }
        //std::cout << i << " THREAD:" << omp_get_thread_num() << std::endl;
  }

  for(auto it:out)
    std::cout << it << '\n';

  return 0;
}
