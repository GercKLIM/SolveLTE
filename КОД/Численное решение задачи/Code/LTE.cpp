#include "LTE.h"

/* Реализация методов решения LTE*/


// Явная схема с левой разностью на 2 точках
bool LTE::SolveLD2e(std::string filename){

    // Открытие файла
    std::ofstream file(filename);
    if (file.is_open()) {

        // Создаем сетки
        int Nx = static_cast<int>(abs(l2 - l1) / h);
        int Nt = static_cast<int>(abs(T - t0) / tau);
        std::vector<double> web_h(Nx);
        std::vector<double> web_tau(Nt);

        for (int x = 0; x <= Nx; x++) {
            web_h[x] = l1 + x * h;
        }
        for (int time = 0; time <= Nt; time++) {
            web_tau[time] = t0 + time * tau;
        }


        // Запись в файл первой строчки с шагом по пространству
        file << "0" << " ";
        for (int x = 0; x <= Nx; x++){
            file << web_h[x] << " ";
        }
        file << std::endl;


        std::vector<double> U_prev(Nx);
        std::vector<double> U_now(Nx);

        // Аппроксимация начального условия
        for (int x = 0; x <= Nx; x ++){
            U_prev[x] = U0(web_h[x]);
        }
        U_now = U_prev;


        // Запись первой итерации <-> начального условия
        file << t0 << " ";
        for (int x = 0; x <= Nx; x ++){
            file << U_prev[x] << " ";
        }
        file << std::endl;

        std::cout << "LOG: Method LD2e complete!" << std::endl;
        for (int time = 1; time <= Nt; time++) {

            for (int x = 1; x <= Nx; x++){
                U_now[x] = (1.0 - gamma)  * U_prev[x] + gamma * U_prev[x - 1];
            }
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x < Nx; x ++){
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

    // Открытие файла
    std::ofstream file(filename);
    if (file.is_open()) {

        // Создаем сетки
        int Nx = static_cast<int>(abs(l2 - l1) / h);
        int Nt = static_cast<int>(abs(T - t0) / tau);
        std::vector<double> web_h(Nx);
        std::vector<double> web_tau(Nt);

        for (int x = 0; x <= Nx; x++) {
            web_h[x] = l1 + x * h;
        }
        for (int time = 0; time <= Nt; time++) {
            web_tau[time] = t0 + time * tau;
        }


        // Запись в файл первой строчки с шагом по пространству
        file << "0" << " ";
        for (int x = 0; x <= Nx; x++){
            file << web_h[x] << " ";
        }
        file << std::endl;


        std::vector<double> U_prev(Nx);
        std::vector<double> U_now(Nx);

        // Аппроксимация начального условия
        for (int x = 0; x <= Nx; x ++){
            U_prev[x] = U0(web_h[x]);
        }
        U_now = U_prev;


        // Запись первой итерации <-> начального условия
        file << t0 << " ";
        for (int x = 0; x <= Nx; x ++){
            file << U_prev[x] << " ";
        }
        file << std::endl;

        std::cout << "LOG: Method Lax complete!" << std::endl;
        for (int time = 1; time <= Nt; time++) {

            for (int x = 1; x <= Nx; x++){

                // Cхема
                U_now[x] = ((U_prev[x+1] + U_prev[x-1]) - gamma * (U_prev[x+1] - U_prev[x-1])) / 2;
                //U_now[x] = (1.0 - gamma)  * U_prev[x] + gamma * U_prev[x - 1];


            }
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x < Nx; x ++){
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

// Схема Лакса-Вендрофа
bool LTE::SolveLaxWen(std::string filename){

    // Открытие файла
    std::ofstream file(filename);
    if (file.is_open()) {

        // Создаем сетки
        int Nx = static_cast<int>(abs(l2 - l1) / h);
        int Nt = static_cast<int>(abs(T - t0) / tau);
        std::vector<double> web_h(Nx);
        std::vector<double> web_tau(Nt);

        for (int x = 0; x <= Nx; x++) {
            web_h[x] = l1 + x * h;
        }
        for (int time = 0; time <= Nt; time++) {
            web_tau[time] = t0 + time * tau;
        }


        // Запись в файл первой строчки с шагом по пространству
        file << "0" << " ";
        for (int x = 0; x <= Nx; x++){
            file << web_h[x] << " ";
        }
        file << std::endl;


        std::vector<double> U_prev(Nx);
        std::vector<double> U_now(Nx);
        double F_p = 0, F_l = 0;

        // Аппроксимация начального условия
        for (int x = 0; x <= Nx; x ++){
            U_prev[x] = U0(web_h[x]);
        }
        U_now = U_prev;


        // Запись первой итерации <-> начального условия
        file << t0 << " ";
        for (int x = 0; x <= Nx; x ++){
            file << U_prev[x] << " ";
        }
        file << std::endl;

        std::cout << "LOG: Method LaxWen complete!" << std::endl;
        for (int time = 1; time <= Nt; time++) {

            for (int x = 1; x <= Nx; x++){

                // Cхема
                F_p = ((U_prev[x+1] + U_prev[x]) - gamma * (U_prev[x+1] - U_prev[x])) / 2;
                F_l = ((U_prev[x] + U_prev[x-1]) - gamma * (U_prev[x] - U_prev[x-1])) / 2;
                U_now[x] = U_prev[x] - gamma * (F_p - F_l);
                //U_now[x] = (1.0 - gamma)  * U_prev[x] + gamma * U_prev[x - 1];


            }
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x < Nx; x ++){
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