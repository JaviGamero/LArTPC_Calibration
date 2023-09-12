// Libraries
#include <string>
#include <iostream>
#include <fstream>
#include <unistd.h>

// namespaces
using namespace std::filesystem;
using namespace std;

void test()
{
    int array1[] = {1,2,3,4,5,6,7,8,9,10};
    int array2[] = {9,8,7,6,5,4,3,2,1,0};

    int *array2d[2]; 

    array2d[0] = array1;
    array2d[1] = array2;

    cout << array2d[0][0] << endl;
    cout << array2d[1][2] << endl;

    return ;
}