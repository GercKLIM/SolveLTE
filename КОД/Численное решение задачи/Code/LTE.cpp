#include "LTE.h"


/* РЕАЛИЗАЦИЯ МЕТОДОВ РЕШЕНИЯ LTE*/

/* Cхема LD2e - Явная с левой разностью по 2-м точкам */
bool LTE::SolveLD2e(std::string filename) {

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

            // Применение ГУ
            U_now[0] = GU(web_h[0], web_tau[time]);

            for (int x = 1; x <= Nx; x++) {
                /* схема LD2e: y^ = (1 - gamma) * y + gamma * y_-1 */
                U_now[x] = (1.0 - gamma) * U_prev[x] + gamma * U_prev[x-1];

            }

            // Обновление вектора значений
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG[info]: Method LD2e complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG[err]: DON'T open file: " << filename << std::endl;
        return false;
    }
}

/* Cхема LD2i - Неявная с левой разностью по 2-м точкам */
bool LTE::SolveLD2i(std::string filename) {

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

            // Применение ГУ
            U_now[0] = GU(web_h[0], web_tau[time]);

            for (int x = 1; x <= Nx; x++) {
                /* схема LD2i: y^ = gamma / (1 - gamma) * y^_-1 + 1 / (1 + gamma) * y */
                U_now[x] = (gamma / (gamma + 1.)) * U_now[x-1] + (1. / (1. + gamma)) * U_prev[x];

            }

            // Обновление вектора значений
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG[info]: Method LD2i complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG[err]: DON'T open file: " << filename << std::endl;
        return false;
    }
}

/* Старая версия LD2i c cупер точностью */
//bool LTE::SolveLD2i(std::string filename) {
//
//    std::ofstream file(filename);
//    if (file.is_open()) {
//
//        int Nx = static_cast<int>(abs(l2 - l1) / h);
//        int Nt = static_cast<int>(abs(T - t0) / tau);
//        std::vector<double> web_h(Nx + 1);
//        std::vector<double> web_tau(Nt + 1);
//
//        for (int x = 0; x <= Nx; x++) {
//            web_h[x] = l1 + x * h;
//        }
//
//        for (int time = 0; time <= Nt; time++) {
//            web_tau[time] = t0 + time * tau;
//        }
//
//        file << "0 ";
//        for (int x = 0; x <= Nx; x++) {
//            file << web_h[x] << " ";
//        }
//        file << std::endl;
//
//        std::vector<double> U_prev(Nx + 1);
//        std::vector<double> U_now(Nx + 1);
//
//        for (int x = 0; x <= Nx; x++) {
//            U_prev[x] = U0(web_h[x]);
//        }
//        U_now = U_prev;
//
//        file << t0 << " ";
//        for (int x = 0; x <= Nx; x++) {
//            file << U_prev[x] << " ";
//        }
//        file << std::endl;
//
//        for (int time = 1; time <= Nt; time++) {
//            std::vector<std::vector<double>> A(Nx + 1, std::vector<double>(Nx + 1, 0.0));
//            std::vector<double> B(Nx + 1);
//
//            for (int x = 1; x < Nx; x++) {
//                A[x][x - 1] = -gamma;
//                A[x][x] = 1 + gamma;
//                B[x] = U_prev[x];
//            }
//
//            // Граничное условие слева, задаваемое функцией U_left
//            A[0][0] = 1.0;
//            B[0] = GU_left(web_tau[time]);
//
//            // Правое граничное условие периодическое
//            A[Nx][Nx] = 1.0;
//            B[Nx] = U_prev[Nx];
//
//            std::vector<double> U_next = SolveTriDiagonal(A, B);
//
//            U_prev = U_next;
//
//            file << web_tau[time] << " ";
//            for (int x = 0; x <= Nx; x++) {
//                file << U_next[x] << " ";
//            }
//            file << std::endl;
//        }
//
//        file.close();
//        std::cout << "LOG: Implicit method complete! Gamma = " << gamma << " Result saved on: " << filename
//                  << std::endl;
//        return true;
//    } else {
//        std::cout << "LOG: Couldn't open file: " << filename << std::endl;
//        return false;
//    }
//}
//std::vector<double> LTE::SolveTriDiagonal(const std::vector<std::vector<double>>& A, const std::vector<double>& B) {
//    int n = B.size();
//    std::vector<double> x(n);
//    std::vector<double> alpha(n), beta(n);
//
//    alpha[0] = A[0][1] / A[0][0];
//    beta[0] = B[0] / A[0][0];
//
//    for (int i = 1; i < n - 1; ++i) {
//        double m = 1.0 / (A[i][i] - A[i][i - 1] * alpha[i - 1]);
//        alpha[i] = A[i][i + 1] * m;
//        beta[i] = (B[i] - A[i][i - 1] * beta[i - 1]) * m;
//    }
//
//    double m = 1.0 / (A[n - 1][n - 1] - A[n - 1][n - 2] * alpha[n - 2]);
//    beta[n - 1] = (B[n - 1] - A[n - 1][n - 2] * beta[n - 2]) * m;
//
//    x[n - 1] = beta[n - 1];
//    for (int i = n - 2; i >= 0; --i) {
//        x[i] = alpha[i] * x[i + 1] + beta[i];
//    }
//
//    return x;
//}



/* Cхема LD3e - Явная с левой разностью по 3-м точкам */
bool LTE::SolveLD3e(std::string filename) {

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

            // Применение ГУ
            U_now[0] = GU(web_h[0], web_tau[time]);

            // Получение второй точки верхнего слоя с помощью LD2i
            U_now[1] = (gamma / (gamma + 1.)) * U_now[0] + (1. / (1. + gamma)) * U_prev[1];

            for (int x = 1; x <= Nx; x++) {
                /* схема LD3e: y^ = (1 - 1.5 * gamma) y + gamma * (2 * y_-1 - 0.5 * y_-2) */
                U_now[x] =  (1. - 1.5 * gamma) * U_prev[x] + gamma * (2. * U_prev[x-1] - 0.5 * U_prev[x-2]);

            }

            // Обновление вектора значений
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG[info]: Method LD3e complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG[err]: DON'T open file: " << filename << std::endl;
        return false;
    }
}


/* Cхема LD3i - Неявная с левой разностью по 3-м точкам */
bool LTE::SolveLD3i(std::string filename) {

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

            // Применение ГУ
            U_now[0] = GU(web_h[0], web_tau[time]);

            // Получение второй точки верхнего слоя с помощью LD2i
            U_now[1] = (gamma / (gamma + 1.)) * U_now[0] + (1. / (1. + gamma)) * U_prev[1];

            for (int x = 2; x <= Nx; x++) {
                /* схема LD3i: y^ = (2 / (2 + 3 * gamma)) y + (gamma / (2 + 3 * gamma)) * (4 * y^_-1 - y^_-2) */
                U_now[x] =  (2. / (2. + 3. * gamma)) * U_prev[x] + (gamma / (2. + 3. * gamma)) * (4. * U_now[x-1] - U_now[x-2]);

            }

            // Обновление вектора значений
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG[info]: Method LD3i complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG[err]: DON'T open file: " << filename << std::endl;
        return false;
    }
}

/* Решение СЛАУ методом прогонки для LD3i */
//std::vector<double> LTE::TridiagonalSolver(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d) {
//    int n = b.size();
//    std::vector<double> p(n), q(n), x(n);
//
//    p[0] = -c[0] / b[0];
//    q[0] = d[0] / b[0];
//
//    for (int i = 1; i < n; ++i) {
//        double denom = b[i] + a[i] * p[i - 1];
//        p[i] = -c[i] / denom;
//        q[i] = (d[i] - a[i] * q[i - 1]) / denom;
//    }
//
//    x[n - 1] = q[n - 1];
//    for (int i = n - 2; i >= 0; --i) {
//        x[i] = p[i] * x[i + 1] + q[i];
//    }
//
//    return x;
//}
//bool LTE::SolveLD3i(std::string filename) {
//    std::ofstream file(filename);
//    if (file.is_open()) {
//        int Nx = static_cast<int>(std::abs(l2 - l1) / h);
//        int Nt = static_cast<int>(std::abs(T - t0) / tau);
//
//        std::vector<double> web_h(Nx + 1);
//        std::vector<double> web_tau(Nt + 1);
//
//        for (int x = 0; x <= Nx; x++) {
//            web_h[x] = l1 + x * h;
//        }
//
//        for (int time = 0; time <= Nt; time++) {
//            web_tau[time] = t0 + time * tau;
//        }
//
//        file << "0 ";
//        for (int x = 0; x <= Nx; x++) {
//            file << web_h[x] << " ";
//        }
//        file << std::endl;
//
//        std::vector<double> U_prev(Nx + 1);
//        std::vector<double> U_now(Nx + 1);
//
//        for (int x = 0; x <= Nx; x++) {
//            U_prev[x] = U0(web_h[x]);
//        }
//        U_now = U_prev;
//
//        file << t0 << " ";
//        for (int x = 0; x <= Nx; x++) {
//            file << U_prev[x] << " ";
//        }
//        file << std::endl;
//
//        std::vector<double> a(Nx + 1, 0);
//        std::vector<double> b(Nx + 1, 0);
//        std::vector<double> c(Nx + 1, 0);
//        std::vector<double> d(Nx + 1, 0);
//
//        for (int time = 1; time <= Nt; time++) {
//            for (int x = 0; x <= Nx; x++) {
//                d[x] = U_prev[x];
//            }
//
//            b[0] = 1.0;
//            c[0] = 0.0;
//            d[0] = GU(web_h[0], web_tau[time]);
//
//            for (int x = 1; x < Nx; x++) {
//                a[x] = -1.0 / h;
//                b[x] = 1.0 / h + 1.0 / tau;
//                c[x] = -1.0 / tau;
//                d[x] = U_prev[x] / tau;
//            }
//
//            b[Nx] = 1.0;
//            a[Nx] = 0.0;
//            d[Nx] = GU(web_h[Nx], web_tau[time]);
//
//            U_now = TridiagonalSolver(a, b, c, d);
//
//            U_prev = U_now;
//
//            file << web_tau[time] << " ";
//            for (int x = 0; x <= Nx; x++) {
//                file << U_now[x] << " ";
//            }
//            file << std::endl;
//        }
//
//        file.close();
//        std::cout << "LOG[info]: Method LD3i complete! Result saved on: " << filename << std::endl;
//        return true;
//    } else {
//        std::cout << "LOG[err]: Couldn't open file: " << filename << std::endl;
//        return false;
//    }
//}



/* Cхема Lax - Схема Лакса */
bool LTE::SolveLax(std::string filename) {

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

            // Применение ГУ
            U_now[0] = GU(web_h[0], web_tau[time]);
            U_now[Nx] = GU(web_h[Nx], web_tau[time]);

            for (int x = 1; x <= (Nx - 1); x++) {
                /* схема Lax: y^ = ((y_+1 + y_-1) - gamma * (y_+1 - y_-1)) / 2 */
                U_now[x] = ((U_prev[x+1] + U_prev[x-1]) - gamma * (U_prev[x+1] - U_prev[x-1])) / 2.;
            }

            // Обновление вектора значений
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG[info]: Method Lax complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG[err]: DON'T open file: " << filename << std::endl;
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
        double F_p = 0, F_l = 0;

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

            // Применение ГУ
            U_now[0] = GU(web_h[0], web_tau[time]);
            U_now[Nx] = GU(web_h[Nx], web_tau[time]);

            for (int x = 1; x <= (Nx - 1); x++) {
                /* схема LaxWen: y^ = y - gamma * (F_p - F_l),
                 *               F_p = ((y_+1 + y) - gamma * (y_+1 - y)) / 2,
                 *               F_l = ((y + y_-1) - gamma * (y - y_-1)) / 2.
                 * */
                F_p = ((U_prev[x+1] + U_prev[x]) - gamma * (U_prev[x+1] - U_prev[x])) / 2;
                F_l = ((U_prev[x] + U_prev[x-1]) - gamma * (U_prev[x] - U_prev[x-1])) / 2;
                U_now[x] = U_prev[x] - gamma * (F_p - F_l);
            }

            // Обновление вектора значений
            U_prev = U_now;

            // Запись в файл
            file << web_tau[time] << " ";
            for (int x = 0; x <= Nx; x++) {
                file << U_now[x] << " ";
            }
            file << std::endl;
        }

        file.close();
        std::cout << "LOG[info]: Method LaxWen complete! Gamma = " << gamma << " Result saved on: " << filename << std::endl;
        return true;
    } else {
        std::cout << "LOG[err]: DON'T open file: " << filename << std::endl;
        return false;
    }
}

