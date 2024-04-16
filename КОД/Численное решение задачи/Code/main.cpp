/* ### Численные методы решения линейного одномерного уравнения переноса ###
Автор: Климов О.Д. ФН2-61Б, GitHub: @gerc.klim
 */



#include <iostream>
#include "LTE.h"


/* Тест 1: Левый треугольник */
void test1() {
    LTE test1;
    test1.a = 1;
    test1.T = 1;

    test1.l1 = -1.0;
    test1.l2 = 1.0;

    test1.h = 0.1;
    test1.tau = 0.1;
    test1.set_U0([&](double x) { return ((x - test1.l1) / (test1.l2 - test1.l1)) ;});
    test1.gamma = test1.a * test1.tau / test1.h;
    //test1.info();

    test1.SolveLD2e("./Output/test1/test1_LD2e.txt");

    test1.SolveLax("./Output/test1/test1_Lax.txt");

    test1.SolveLaxWen("./Output/test1/test1_LaxWen.txt");

    test1.SolveLD2e("./Output/test1/test1_LD2e.txt");

    test1.SolveLD2e("./Output/test1/test1_LD2e.txt");

    test1.SolveLD2e("./Output/test1/test1_LD2e.txt");
}


/* Тест 2: Прямоугольник */
void test2() {
    /* Тест 2: Прямоугольник */
    LTE test2;
    test2.a = 1;
    test2.T = 1;

    test2.l1 = -1.0;
    test2.l2 = 1.0;

    test2.h = 0.1;
    test2.tau = 0.1;
    test2.set_U0([&](double x) { return 1.0 ;});
    test2.gamma = test2.a * test2.tau / test2.h;
    //test2.info();

    test2.SolveLD2e("./Output/test2/test2_LD2e.txt");

    test2.SolveLax("./Output/test2/test2_Lax.txt");

    test2.SolveLaxWen("./Output/test2/test2_LaxWen.txt");

    test2.SolveLD2e("./Output/test2/test2_LD2e.txt");

    test2.SolveLD2e("./Output/test2/test2_LD2e.txt");

    test2.SolveLD2e("./Output/test2/test2_LD2e.txt");
}


/* Тест 3: Косинус */
void test3(){
    /* Тест 3: Косинус */
    LTE test3;
    test3.a = 1;
    test3.T = 1;

    test3.l1 = -1.0;
    test3.l2 = 1.0;

    test3.h = 0.1;
    test3.tau = 0.1;
    test3.set_U0([&](double x) { return (0.5 * (1 - cos( ((2 * M_PI * ( x - test3.l1)) / (test3.l2 - test3.l1)))));});
    test3.gamma = test3.a * test3.tau / test3.h;
    //test3.info();

    test3.SolveLD2e("./Output/test3/test3_LD2e.txt");

    test3.SolveLax("./Output/test3/test3_Lax.txt");

    test3.SolveLaxWen("./Output/test3/test3_LaxWen.txt");

    test3.SolveLD2e("./Output/test3/test3_LD2e.txt");

    test3.SolveLD2e("./Output/test3/test3_LD2e.txt");

    test3.SolveLD2e("./Output/test3/test3_LD2e.txt");
}


/* Тест 4: Зуб */
void test4() {
    LTE test4;
    test4.a = 1;
    test4.T = 1;

    test4.l1 = -1.0;
    test4.l2 = 1.0;

    test4.h = 0.1;
    test4.tau = 0.1;
    double l11 = -0.1, l22 = 0.1;
    test4.set_U0([&](double x) {
        if ((x >= test4.l1) and (x < l11)) {
            return (-2 / 3 * (l11 - test4.l1) * (x - test4.l1) + 1);
        } else if ((x >= l11) and (x <= l22)) {
            return 1./3;
        } else if ((x > l22) and (x < test4.l2)) {
            return (2/3 * (test4.l2 - l22) * (x - test4.l2) + 1);
        } else {
            return 0.0;
        }
    });
    test4.gamma = test4.a * test4.tau / test4.h;
    //test4.info();


    test4.SolveLD2e("./Output/test4/test4_LD2e.txt");

    test4.SolveLax("./Output/test4/test4_Lax.txt");

    test4.SolveLaxWen("./Output/test4/test4_LaxWen.txt");

    test4.SolveLD2e("./Output/test4/test4_LD2e.txt");

    test4.SolveLD2e("./Output/test4/test4_LD2e.txt");

    test4.SolveLD2e("./Output/test4/test4_LD2e.txt");

}


/* Тест 5: М */
void test5() {
    LTE test5;
    test5.a = 1;
    test5.T = 1;

    test5.l1 = -1.0;
    test5.l2 = 1.0;

    test5.h = 0.1;
    test5.tau = 0.1;

    double l12 = 0.;
    test5.set_U0([&](double x) {
        if ((x >= test5.l1) and (x < l12)) {
            return (-2 / 3 * (l12 - test5.l1) * (x - test5.l1) + 1);
        } else if ((x >= l12) and (x < test5.l2)) {
            return (2/3 * (test5.l2 - l12) * (x - test5.l2) + 1);
        } else {
            return 0.0;
        }
    });
    test5.gamma = test5.a * test5.tau / test5.h;
    //test5.info();

    test5.SolveLD2e("./Output/test5/test5_LD2e.txt");

    test5.SolveLax("./Output/test5/test5_Lax.txt");

    test5.SolveLaxWen("./Output/test5/test5_LaxWen.txt");

    test5.SolveLD2e("./Output/test5/test5_LD2e.txt");

    test5.SolveLD2e("./Output/test5/test5_LD2e.txt");

    test5.SolveLD2e("./Output/test5/test5_LD2e.txt");

}


/* Тест 6: Правый треугольник */
void test6() {
    LTE test6;
    test6.a = 1;
    test6.T = 1;

    test6.l1 = -1.0;
    test6.l2 = 1.0;

    test6.h = 0.1;
    test6.tau = 0.1;
    test6.set_U0([&](double x) { return ((test6.l2 - x) / (test6.l2 - test6.l1)) ;});
    test6.gamma = test6.a * test6.tau / test6.h;
    //test6.info();

    test6.SolveLD2e("./Output/test6/test6_LD2e.txt");

    test6.SolveLax("./Output/test6/test6_Lax.txt");

    test6.SolveLaxWen("./Output/test6/test6_LaxWen.txt");

    test6.SolveLD2e("./Output/test6/test6_LD2e.txt");

    test6.SolveLD2e("./Output/test6/test6_LD2e.txt");

    test6.SolveLD2e("./Output/test6/test6_LD2e.txt");
}



int main() {

    test1();
    std::cout << "Complete!!!" << std::endl;
    return 0;
}
