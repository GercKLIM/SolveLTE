/* ### Численные методы решения линейного одномерного уравнения переноса ###
Автор: Климов О.Д. ФН2-61Б, GitHub: @gerc.klim
 */



#include <iostream>
#include "LTE.h"


/* Тест 1: Левый треугольник */
void test1(double h = 0.1, double tau = 0.1, double L1 = -10., double L2 = 10., double T = 10., double l1 = -1, double l2 = 1) {
    LTE test1;
    test1.a = 1;
    test1.T = T;

    test1.l1 = L1;
    test1.l2 = L2;

    test1.h = h;
    test1.tau = tau;

    test1.set_U0([&](double x) { return (x - l1) / (l2 - l1);});

    test1.set_GU_left([&](double t) { return 1./2. * (1. - t + L1);});

    test1.gamma = test1.a * test1.tau / test1.h;
    //test1.info();


    // Явная схема с левой разностью 1-го пор.
    test1.SolveLD2e("./Output/test1/test1_LD2e.txt");

    // Неявная схема с левой разностью 1-го пор.
    test1.SolveLD2i("./Output/test1/test1_LD2i.txt");

    // Схема Лакса
    test1.SolveLax("./Output/test1/test1_Lax.txt");

    // Схема Лакса-Вендрофа
    test1.SolveLaxWen("./Output/test1/test1_LaxWen.txt");

    // Явная схема с левой разностью 2-го пор.
    test1.SolveLD3e("./Output/test1/test1_LD3e.txt");

    // Неявная схема с левой разностью 2-го пор.
    test1.SolveLD3i("./Output/test1/test1_LD3i.txt");
}

/* Тест 2: Правый треугольник */
void test2(double h = 0.1, double tau = 0.1, double L1 = -10., double L2 = 10., double T = 10., double l1 = -1., double l2 = 1.) {
    LTE test2;
    test2.a = 1;
    test2.T = T;

    test2.l1 = L1;
    test2.l2 = L2;

    test2.h = h;
    test2.tau = tau;
    test2.set_U0([&](double x) { return ((l2 - x) / (l2 - l1)) ;});
    test2.set_GU_left([&](double t) { return 1./2. * (1. + t - L1);});
    test2.gamma = test2.a * test2.tau / test2.h;
    //test6.info();

    test2.SolveLD2e("./Output/test2/test2_LD2e.txt");

    test2.SolveLax("./Output/test2/test2_Lax.txt");

    test2.SolveLaxWen("./Output/test2/test2_LaxWen.txt");

    test2.SolveLD3e("./Output/test2/test2_LD3e.txt");

    test2.SolveLD2i("./Output/test2/test2_LD2i.txt");

    test2.SolveLD3i("./Output/test2/test2_LD3i.txt");
}


/* Тест 3: Прямоугольник */
void test3(double h = 0.1, double tau = 0.1, double L1 = -10., double L2 = 10., double T = 10., double l = 1.) {
    /* Тест 3: Прямоугольник */
    LTE test3;
    test3.a = 1;
    test3.T = T;

    test3.l1 = L1;
    test3.l2 = L2;

    test3.h = h;
    test3.tau = tau;
    test3.set_U0([&](double x) { return l ;});
    test3.set_GU_left([&](double t) { return l;});
    test3.gamma = test3.a * test3.tau / test3.h;
    //test3.info();

    test3.SolveLD2e("./Output/test3/test3_LD2e.txt");

    test3.SolveLax("./Output/test3/test3_Lax.txt");

    test3.SolveLaxWen("./Output/test3/test3_LaxWen.txt");

    test3.SolveLD3e("./Output/test3/test3_LD3e.txt");

    test3.SolveLD2i("./Output/test3/test3_LD2i.txt");

    test3.SolveLD3i("./Output/test3/test3_LD3i.txt");
}


/* Тест 4: Косинус */
void test4(double h = 0.1, double tau = 0.1, double L1 = -10., double L2 = 10., double T = 10., double l11 = -1., double l22 = 1.){
    /* Тест 4: Косинус */
    LTE test4;
    test4.a = 1;
    test4.T = T;

    test4.l1 = L1;
    test4.l2 = L2;

    test4.h = h;
    test4.tau = tau;
    test4.set_U0([&](double x) { return (1./3. * (1 - cos( ((2 * M_PI * ( x - l11)) / (l22 - l11)))));});
    test4.set_GU_left([&](double t) { return 1./3. * (1 - cos(M_PI * (1 - t + L1)));});
    test4.gamma = test4.a * test4.tau / test4.h;
    //test3.info();

    test4.SolveLD2e("./Output/test4/test4_LD2e.txt");

    test4.SolveLax("./Output/test4/test4_Lax.txt");

    test4.SolveLaxWen("./Output/test4/test4_LaxWen.txt");

    test4.SolveLD3e("./Output/test4/test4_LD3e.txt");

    test4.SolveLD2i("./Output/test4/test4_LD2i.txt");

    test4.SolveLD3i("./Output/test4/test4_LD3i.txt");
}


/* Тест 5: Зуб */
void test5(double h = 0.1, double tau = 0.1, double L1 = -10., double L2 = 10., double T = 10., double l1 = -1, double l2 = 1, double l11 = -0.5, double l22 = 0.5) {
    LTE test5;
    test5.a = 1;
    test5.T = T;

    test5.l1 = L1;
    test5.l2 = L2;

    test5.h = h;
    test5.tau = tau;;
    test5.set_U0([&](double x) {
        if ((x >= l1) and (x < l11)) {
            return (-2./ 3. * (l11 - l1) * (x - l1) + 1.);
        } else if ((x >= l11) and (x <= l22)) {
            return 1./3.;
        } else if ((x > l22) and (x <= l2)) {
            return (2./3. * (l2 - l22) * (x - l2) + 1.);
        } else {
            return 0.;
        }
    });

    test5.set_GU_left([&](double x) {
        if ((x >= l1) and (x < l11)) {
            return (-2./ 3. * (l11 - l1) * (x - l1) + 1.);
        } else if ((x >= l11) and (x <= l22)) {
            return 1./3.;
        } else if ((x > l22) and (x <= l2)) {
            return (2./3. * (l2 - l22) * (x - l2) + 1.);
        } else {
            return 0.;
        }
    });
    test5.gamma = test5.a * test5.tau / test5.h;
    //test4.info();


    test5.SolveLD2e("./Output/test5/test5_LD2e.txt");

    test5.SolveLax("./Output/test5/test5_Lax.txt");

    test5.SolveLaxWen("./Output/test5/test5_LaxWen.txt");

    test5.SolveLD3e("./Output/test5/test5_LD3e.txt");

    test5.SolveLD2i("./Output/test5/test5_LD2i.txt");

    test5.SolveLD3i("./Output/test5/test5_LD3i.txt");

}


/* Тест 6: М */
//void test6(double h = 0.1, double tau = 0.1, double L1 = -10., double L2 = 10., double T = 10., double l1 = -1., double l2 = 1., double l12 = 0.0) {
//    LTE test6;
//    test6.a = 1;
//    test6.T = T;
//
//    test6.l1 = L1;
//    test6.l2 = L2;
//
//    test6.h = h;
//    test6.tau = tau;
//
//    test6.set_U0([&](double x) {
//        if ((x >= l1) and (x < l12)) {
//            return (-2. / 3 * (l12 - l1) * (x - l1) + 1);
//        } else if ((x >= l12) and (x < l2)) {
//            return (2./3 * (l2 - l12) * (x - l2) + 1);
//        } else {
//            return 0.0;
//        }
//    });
//    test6.gamma = test6.a * test6.tau / test6.h;
//    //test5.info();
//
//    test6.SolveLD2e("./Output/test5/test5_LD2e.txt");
//
//    test6.SolveLax("./Output/test5/test5_Lax.txt");
//
//    test6.SolveLaxWen("./Output/test5/test5_LaxWen.txt");
//
//    test6.SolveLD3e("./Output/test5/test5_LD3e.txt");
//
//    test6.SolveLD2i("./Output/test5/test5_LD2i.txt");
//
//    test6.SolveLD3i("./Output/test5/test5_LD3i.txt");
//
//}



void make_test_dataset(double h, double tau){

    double L1 = -10., L2 = 10.;
    double T = 1.;

    // Тест 1
    LTE test1;
    test1.a = 1;
    test1.T = T;

    test1.l1 = L1;
    test1.l2 = L2;

    test1.h = h;
    test1.tau = tau;

    test1.set_U0([&](double x) { return (x + 1) / 2;});
    //test1.set_U0([&](double x) { return (x - l1) / (l2 - l1);});
    test1.gamma = test1.a * test1.tau / test1.h;

    // Тест 2
    LTE test2;
    test2.a = 1;
    test2.T = T;

    test2.l1 = L1;
    test2.l2 = L2;

    test2.h = h;
    test2.tau = tau;


    test2.set_U0([&](double x) { return ((1. - x) / (1. + 1.)) ;});
    test2.gamma = test2.a * test2.tau / test2.h;

    // Тест 3
    LTE test3;
    test3.a = 1;
    test3.T = T;

    test3.l1 = L1;
    test3.l2 = L2;
    test3.h = h;
    test3.tau = tau;

    test3.set_U0([&](double x) { return 2./3. ;});
    test3.gamma = test3.a * test3.tau / test3.h;

    // Тест 4
    LTE test4;
    test4.a = 1;
    test4.T = T;

    test4.l1 = L1;
    test4.l2 = L2;

    test4.h = h;
    test4.tau = tau;

    test4.set_U0([&](double x) { return (1./3. * (1 - cos( ((2 * M_PI * ( x + 1.)) / (1. + 1.)))));});
    test4.gamma = test4.a * test4.tau / test4.h;




    LTE test5;
    test5.a = 1;
    test5.T = T;

    test5.l1 = L1;
    test5.l2 = L2;

    test5.h = h;
    test5.tau = tau;;

    test5.set_U0([&](double x) {
        if ((x >= -0.2) and (x < -0.1)) {
            return (-2./ 3. * (-0.1 + 0.2) * (x + 0.2) + 1.);
        } else if ((x >= -0.1) and (x <= 0.1)) {
            return 1./3.;
        } else if ((x > 0.1) and (x <= 0.2)) {
            return (2./3. * (0.2 - 0.1) * (x - 0.2) + 1.);
        } else {
            return 0.;
        }
    });
    test5.gamma = test5.a * test5.tau / test5.h;


    test1.SolveLD2e("./Output/datasets/test1_" + std::to_string(tau) + "_" + std::to_string(h) + ".txt");
    test2.SolveLD2e("./Output/datasets/test2_" + std::to_string(tau) + "_" + std::to_string(h) + ".txt");
    test3.SolveLD2e("./Output/datasets/test3_" + std::to_string(tau) + "_" + std::to_string(h) + ".txt");
    test4.SolveLD2e("./Output/datasets/test4_" + std::to_string(tau) + "_" + std::to_string(h) + ".txt");
    test5.SolveLD2e("./Output/datasets/test5_" + std::to_string(tau) + "_" + std::to_string(h) + ".txt");

}


int main() {

    /* Тест 1. Левый треугольник */
    test1(0.1, 0.01, -100., 100., 10, -1., 1.);

    /* Тест 2. Правый треугольник */
    //test2(0.1, 0.1, -10., 10., 10., -1., 1.);

    /* Тест 3. Прямоугольник */
    //test3(0.1, 0.1, -10., 10., 10., 1.);

    /* Тест 4. Косинус */
    //test4(0.1, 0.1, -10., 10., 10., -1., 1.);

    /* Тест 5. Зуб */
    //test5(0.1, 0.1, -10., 10., 10., -3., 3., -1., 1.);

    //make_test_dataset(0.1, 0.1);
    //make_test_dataset(0.1, 0.05);
    //make_test_dataset(0.05, 0.05);
    //make_test_dataset(0.01, 0.01);

    std::cout << "Complete!!!" << std::endl;
    return 0;
}
