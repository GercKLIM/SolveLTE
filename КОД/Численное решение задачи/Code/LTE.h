#pragma once
#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <fstream>
#include "./Libs/algebra.h"

/* Класс Линейных Уравнений Переноса */
class LTE {
    /* Уравнение вида:
     * du/dt + a * du/dx = 0,
     * u(x, 0) = U0(x).
     */

public:
    /* Параметры класса */

    double a = 1;       // Параметр a
    double t0 = 0.0;    // Начальное время
    double T = 1;       // Конец временного учатска
    double h = 0.1;     // Шаг по пространству
    double tau = 0.1;   // Шаг по времени
    double l1 = -1, l2 = 1; // Область
    double gamma;

    // Функция задания U0
    void set_U0(std::function<double(double)> func) {
        U0 = func;
        U0_is_set = true;
    }

    double u0(double x){
        return U0(x);
    }
    // Функция для вывода информации об объекте класса в консоль
    void info(){
        std::cout << "LTE object information:" << std::endl;
        std::cout << "a = " << a << std::endl;
        std::cout << "T = " << T << std::endl;
        std::cout << "l = (" << l1 << ", " << l2 << ")" << std::endl;
        std::cout << "h = " << h << std::endl;
        std::cout << "tau = " << tau << std::endl;
        std::cout << "gamma = " << /*a * tau / h*/ gamma << std::endl;
        std::cout << "U0 is " << ((U0_is_set) ? "set" : "NOT set") << std::endl;
    }

    /* Заголовки методов численного решения */

    // Явная схема с левой разностью на 2 точках
    bool SolveLD2e(std::string filename);

    std::vector<double> SolveTriDiagonal(const std::vector<std::vector<double>>& A, const std::vector<double>& B);
    // Невная схема с левой разностью на 2 точках
    bool SolveLD2i(std::string filename);

    // Явная схема с левой разностью на 3 точках
    bool SolveLD3e(std::string filename);

    // Неявная схема с левой разностью на 3 точках
    bool SolveLD3i(std::string filename);

    // Схема Лакса
    bool SolveLax(std::string filename);

    // Схема Лакса-Вендрофа
    bool SolveLaxWen(std::string filename);


private:

    std::function<double(double)> U0;

    /* СОстояние заданности параметров */
//    bool a_is_set = false;
//    bool T_is_set = false;
//    bool h_is_set = false;
//    bool tau_is_set = false;
    bool U0_is_set = false;
};
