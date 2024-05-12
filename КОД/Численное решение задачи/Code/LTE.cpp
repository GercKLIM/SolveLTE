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
        for (int x = 0; x <= Nx; x++){
            U_prev[x] = U0(web_h[x]);
        }
        U_now = U_prev;

        // Запись первой итерации <-> начального условия
        file << t0 << " ";
        for (int x = 0; x <= Nx; x++){
            file << U_prev[x] << " ";
        }
        file << std::endl;

        // Итерационный процесс
        for (int time = 1; time <= Nt; time++) {
            for (int x = 1; x <= Nx; x++) {
                U_now[x] = (1.0 - gamma) * U_prev[x] + gamma * U_prev[x-1];
            }
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG: Method LD2e complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG: DON'T open file: " << filename << std::endl;
        return false;
    }
}

// Функция для решения трехдиагональной системы уравнений
std::vector<double> LTE::SolveTriDiagonal(const std::vector<std::vector<double>>& A, const std::vector<double>& B) {
    int n = B.size();
    std::vector<double> x(n);
    std::vector<double> alpha(n), beta(n);

    alpha[1] = A[0][1] / A[0][0];
    beta[1] = B[0] / A[0][0];

    for (int i = 1; i < n - 1; ++i) {
        double m = 1.0 / (A[i][i] - A[i][i - 1] * alpha[i]);
        alpha[i + 1] = A[i][i + 1] * m;
        beta[i + 1] = (B[i] - A[i][i - 1] * beta[i]) * m;
    }

    x[n - 1] = (B[n - 1] - A[n - 1][n - 2] * beta[n - 1]) / (A[n - 1][n - 1] - A[n - 1][n - 2] * alpha[n - 1]);
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
    }

    return x;
}

// Неявная схема с левой разностью на 2 точках
bool LTE::SolveLD2i(std::string filename){

    // Открытие файла
    std::ofstream file(filename);
    if (file.is_open()) {

        // Создаем сетки
        int Nx = static_cast<int>(abs(l2 - l1) / h);
        int Nt = static_cast<int>(abs(T - t0) / tau);

        std::vector<double> web_h(Nx + 1);
        std::vector<double> web_tau(Nt + 1);

        for (int x = 0; x <= Nx; x++) {
            web_h[x] = l1 + x * h;
        }

        for (int time = 0; time <= Nt; time++) {
            web_tau[time] = t0 + time * tau;
        }

        // Запись в файл первой строчки с шагом по пространству
        file << "0" << " ";
        for (int x = 0; x <= Nx; x++) {
            file << web_h[x] << " ";
        }
        file << std::endl;

        std::vector<double> U_prev(Nx + 1);
        std::vector<double> U_now(Nx + 1);

        // Аппроксимация начального условия
        for (int x = 0; x <= Nx; x++) {
            U_prev[x] = U0(web_h[x]);
        }
        U_now = U_prev;

        // Запись первой итерации <-> начального условия
        file << t0 << " ";
        for (int x = 0; x <= Nx; x++) {
            file << U_prev[x] << " ";
        }
        file << std::endl;

        // Итерационный процесс
        for (int time = 1; time <= Nt; time++) {

            // Формирование матрицы и вектора для неявной схемы
            std::vector<std::vector<double>> A(Nx + 1, std::vector<double>(Nx + 1));
            std::vector<double> B(Nx + 1);

            for (int x = 0; x <= Nx; x++) {
                double a = -gamma / h;
                double b = 1.0 + gamma / h;
                if (x == 0) {
                    A[x][x] = b;
                    A[x][x + 1] = a;
                    B[x] = U_prev[x];
                } else if (x == Nx) {
                    A[x][x] = b;
                    A[x][x - 1] = a;
                    B[x] = U_prev[x];
                } else {
                    A[x][x - 1] = a;
                    A[x][x] = b;
                    A[x][x + 1] = a;
                    B[x] = U_prev[x];
                }
            }

            // Решение системы уравнений
            std::vector<double> U_next = SolveTriDiagonal(A, B);

            // Обновление значений для следующего временного шага
            U_prev = U_next;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_next[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG: Implicit method complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG: Couldn't open file: " << filename << std::endl;
        return false;
    }
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
        std::cout << "LOG: Method Lax complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
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
        std::cout << "LOG: Method LaxWen complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG: DON'T open file: " << filename << std::endl;
        return false;
    }
}