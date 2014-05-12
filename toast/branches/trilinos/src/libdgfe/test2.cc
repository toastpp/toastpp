#include <iostream>
#include <vector>
#include<algorithm>

using namespace std;

main()
{
   vector< vector<int> > vI2Matrix;    // Declare two dimensional array
   vector<int> A, B;
   vector< vector<int> >::iterator iter_ii;
   vector<int>::iterator                 iter_jj;

   A.push_back(10);
   A.push_back(20);
   A.push_back(30);
   B.push_back(100);
   B.push_back(200);
   B.push_back(300);

   vI2Matrix.push_back(A);
   vI2Matrix.push_back(B);

   cout << endl << "Using Iterator:" << endl;

   for(iter_ii=vI2Matrix.begin(); iter_ii!=vI2Matrix.end(); iter_ii++)
   {
      for(iter_jj=(*iter_ii).begin(); iter_jj!=(*iter_ii).end(); iter_jj++)
      {
         cout << *iter_jj << endl;
      }
   }

   int *test; //[7] = {50, 70, 20, 10, 60, 30, 20 };
   test = (int *) malloc(7*sizeof(int));
   test[0] = 50; test[1] = 70; test[2] = 20; test[3] = 10; 
   test[4] = 60; test[5] = 30; test[6] = 20; 
   int n=7;

   std::sort(test, test+n);
   printf("%d %d %d %d %d %d %d\n", test[0], test[1], test[2], test[3], test[4], test[5], test[6]);
}
