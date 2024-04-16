/* ### Численные методы решения линейного одномерного уравнения переноса ###
Автор: Климов О.Д. ФН2-61Б, GitHub: @gerc.klim
 */



#include <iostream>
#include "LTE.h"



int main() {

    /* Тест 1 */
    LTE test1;
    test1.a = 1;
    test1.T = 10;
    test1.l1 = 0;
    test1.l2 = 10;
    test1.h = 0.1;
    test1.tau = 0.1;
    test1.set_U0([&](double x) { return (x - 0) / (10 - 0) ;});
    test1.info();


    test1.SolveLD2e("test1_LD2e.txt");


    std::cout << "Complete!" << std::endl;
    return 0;
}
