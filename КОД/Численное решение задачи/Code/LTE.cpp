#include "LTE.h"

/* Реализация методов решения LTE*/


// Явная схема с левой разностью на 2 точках
bool LTE::SolveLD2e(std::string filename){

    // Открытие файла
    std::ofstream file(filename);
    if (file.is_open()) {

        double gamma = 0.5;

        // Временной слой
        double y_t_next;
        double y_t_now;
        double y_t_prev;

        // Пространственный слой
        double y_x_next;
        double y_x_now;
        double y_x_prev;




        int Nt = static_cast<int>((l2 - l1) / tau);
        int Nx = static_cast<int>( (l2 - l1) / h);

        std::vector<double> U_prev(Nx);
        std::vector<double> U_now(Nx);

        // Аппрокисмация начального условия
        for (int x = 1; x < Nx; x += 1){
            U_prev[x] = U0(x * h);
        }
        U_now = U_prev;

        // Запись в файл первой строчки с шагом по пространству
        file << "0 ";
        for (int x = 1; x < Nx; x += 1){
            file << x * h << " ";
        }
        file << std::endl;

        for (int time = 1; time < Nt; time += 1) {

            for (int x = 1; x < Nx; x += 1){
                U_now[x] = (1 - gamma)  * U_prev[x] + gamma * U_prev[x - h];
            }
            U_prev = U_now;

            // Запись в файл
            file << time * tau << ' ';
            for (int x = 1; x < Nx; x += 1){
                file << U_now[x] << " ";
            }
            file << std::endl;
        }



        file.close();
        return true;
    } else {
        std::cout << "LOG: DON'T open file: " << filename << std::endl;
        return false;
    }
}

// Невная схема с левой разностью на 2 точках
bool LTE::SolveLD2i(std::string filename){
    return false;
}

// Явная схема с левой разностью на 3 точках
bool LTE::SolveLD3e(std::string filename){
    return false;
}

// Неявная схема с левой разностью на 3 точках
bool LTE::SolveLD3i(std::string filename){
    return false;
}

// Схема Лакса
bool LTE::SolveLax(std::string filename){
    return false;
}

// Схема Лакса-Вендрофа
bool LTE::SolveLaxWen(std::string filename){
    return false;
}