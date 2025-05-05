#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define TOCHNOST 1e-8       // точность вычислений
#define MAX_SHAGOV 1000     // максимальное число итераций
#define NAHALO_POISKA -10.0 // начало диапазона поиска
#define KONEC_POISKA 10.0   // конец диапазона поиска
#define SHAG_POISKA 0.1     // начальный шаг поиска

double vichislit_funkciyu(double x) {
    return 3 * x + x * x * x - pow(4, x) + 1;
}

double vichislit_proizvodnuyu(double x) {
    return 3 + 3 * x * x - log(4) * pow(4, x);
}

double vichislit_vtoruyu_proizvodnuyu(double x) {
    return 6 * x - log(4) * log(4) * pow(4, x);
}

double iteriruyushaya_funkciya(double x, int interval_num) {
    if (interval_num == 0) {
        return cbrt(pow(4, x) - 3 * x - 1);
    }
    else if (interval_num == 1) {
        return log(x * x * x + 3 * x + 1) / log(4);
    }
    else {
        return x - vichislit_funkciyu(x) / vichislit_proizvodnuyu(x); // альтернативный вариант
    }
}

double proizvodnaya_iteriruyushaya_funkciya(double x, int interval_num) {
    if (interval_num == 0) {
        double numerator = log(4) * pow(4, x) - 3;
        double denominator = 3 * pow(cbrt(pow(4, x) - 3 * x - 1), 2);
        return numerator / denominator;
    }
    else if (interval_num == 1) {
        double numerator = 3 * x * x + 3;
        double denominator = log(4) * (x * x * x + 3 * x + 1);
        return numerator / denominator;
    }
    else {
        return 1 - (vichislit_proizvodnuyu(x) * vichislit_proizvodnuyu(x) -
            vichislit_funkciyu(x) * vichislit_vtoruyu_proizvodnuyu(x)) /
            (vichislit_proizvodnuyu(x) * vichislit_proizvodnuyu(x));
    }
}

void metod_polovinnogo_deleniya(double a, double b, double* koren, int* shagi) {
    *shagi = 0;
    double fa = vichislit_funkciyu(a);
    double fb = vichislit_funkciyu(b);

    if (fa * fb >= 0) {
        *koren = NAN;
        return;
    }

    while ((b - a) > TOCHNOST && *shagi < MAX_SHAGOV) {
        (*shagi)++;
        double c = (a + b) / 2;
        double fc = vichislit_funkciyu(c);

        if (fc == 0 || fabs(fc) < TOCHNOST) {
            *koren = c;
            return;
        }

        if (fc * fa < 0) {
            b = c;
            fb = fc;
        }
        else {
            a = c;
            fa = fc;
        }
    }
    *koren = (a + b) / 2;
}

void metod_hord(double a, double b, double* koren, int* shagi) {
    *shagi = 0;
    double fa = vichislit_funkciyu(a);
    double fb = vichislit_funkciyu(b);

    if (fa * fb >= 0) {
        *koren = NAN;
        return;
    }

    while (*shagi < MAX_SHAGOV) {
        (*shagi)++;
        double x_nov = a - fa * (b - a) / (fb - fa);
        double fx_nov = vichislit_funkciyu(x_nov);

        if (fabs(fx_nov) < TOCHNOST) {
            *koren = x_nov;
            return;
        }

        if (fx_nov * fa < 0) {
            b = x_nov;
            fb = fx_nov;
        }
        else {
            a = x_nov;
            fa = fx_nov;
        }
    }
    *koren = a - fa * (b - a) / (fb - fa);
}

void metod_newtona(double a, double b, double* koren, int* shagi) {
    *shagi = 0;
    double x_pred = (a + b) / 2; // начальное приближение

    // проверка условий сходимости
    double fa = vichislit_funkciyu(a);
    double fb = vichislit_funkciyu(b);
    double fpa = vichislit_proizvodnuyu(a);
    double fpb = vichislit_proizvodnuyu(b);
    double fppa = vichislit_vtoruyu_proizvodnuyu(a);
    double fppb = vichislit_vtoruyu_proizvodnuyu(b);

    if (fa * fb >= 0 || fabs(fpa) < TOCHNOST || fabs(fpb) < TOCHNOST) {
        *koren = NAN;
        return;
    }

    while (*shagi < MAX_SHAGOV) {
        (*shagi)++;
        double fx = vichislit_funkciyu(x_pred);
        double fpx = vichislit_proizvodnuyu(x_pred);

        if (fabs(fpx) < TOCHNOST) {
            *koren = NAN;
            return;
        }

        double x_nov = x_pred - fx / fpx;

        if (fabs(x_nov - x_pred) < TOCHNOST || fabs(fx) < TOCHNOST) {
            *koren = x_nov;
            return;
        }

        x_pred = x_nov;
    }
    *koren = x_pred;
}

void metod_prostyh_iteraciy(double a, double b, double* koren, int* shagi, int interval_num) {
    *shagi = 0;
    double x_pred = (a + b) / 2;

    bool converges = true;
    double test_points[] = { a, (a + b) / 2, b };
    for (int i = 0; i < 3; i++) {
        double g_prime = fabs(proizvodnaya_iteriruyushaya_funkciya(test_points[i], interval_num));
        if (g_prime >= 1) {
            converges = false;
            break;
        }
    }

    if (!converges) {
        interval_num = 2;
        for (int i = 0; i < 3; i++) {
            double g_prime = fabs(proizvodnaya_iteriruyushaya_funkciya(test_points[i], interval_num));
            if (g_prime >= 1) {
                *koren = NAN;
                return;
            }
        }
    }

    while (*shagi < MAX_SHAGOV) {
        (*shagi)++;
        double x_nov = iteriruyushaya_funkciya(x_pred, interval_num);

        if (fabs(x_nov - x_pred) < TOCHNOST) {
            *koren = x_nov;
            return;
        }

        x_pred = x_nov;
    }
    *koren = x_pred;
}

void nayti_vse_intervaly(double intervaly[][2], int* kolvo) {
    *kolvo = 0;
    double x_tek = NAHALO_POISKA;
    double shag = SHAG_POISKA;

    while (x_tek < KONEC_POISKA) {
        double x_sled = x_tek + shag;
        double f_tek = vichislit_funkciyu(x_tek);
        double f_sled = vichislit_funkciyu(x_sled);

        if (f_tek * f_sled < 0) {
            intervaly[*kolvo][0] = x_tek;
            intervaly[*kolvo][1] = x_sled;
            (*kolvo)++;
        }
        else if (fabs(f_tek) < 10 * TOCHNOST) {
            intervaly[*kolvo][0] = x_tek - shag / 2;
            intervaly[*kolvo][1] = x_tek + shag / 2;
            (*kolvo)++;
        }

        double proizv = vichislit_proizvodnuyu(x_tek);
        shag = fmin(SHAG_POISKA, 0.1 / fmax(fabs(proizv), 1.0));
        x_tek += shag;
    }
}

void metod_sekushih(double a, double b, double* koren, int* shagi) {
    *shagi = 0;
    double x0 = a;
    double x1 = b;
    double fx0 = vichislit_funkciyu(x0);
    double fx1 = vichislit_funkciyu(x1);

    if (fabs(fx0) < fabs(fx1)) {
        double temp = x0; x0 = x1; x1 = temp;
        temp = fx0; fx0 = fx1; fx1 = temp;
    }

    while (*shagi < MAX_SHAGOV) {
        (*shagi)++;

        if (fabs(fx1) < TOCHNOST || fabs(x1 - x0) < TOCHNOST) {
            *koren = x1;
            return;
        }

        double df = (fx1 - fx0) / (x1 - x0);
        if (fabs(df) < TOCHNOST) {
            *koren = NAN;
            return;
        }

        double x_nov = x1 - fx1 / df;
        fx0 = fx1;
        x0 = x1;
        x1 = x_nov;
        fx1 = vichislit_funkciyu(x1);
    }
    *koren = x1;
}

void metod_hord_i_kasatelnyh(double a, double b, double* koren, int* shagi) {
    *shagi = 0;
    double fa = vichislit_funkciyu(a);
    double fb = vichislit_funkciyu(b);
    double fpa, fpb;

    if (fa * fb >= 0) {
        *koren = NAN;
        return;
    }

    fpa = vichislit_proizvodnuyu(a);
    fpb = vichislit_proizvodnuyu(b);

    double x_hord, x_kas;

    while (*shagi < MAX_SHAGOV && fabs(b - a) > TOCHNOST) {
        (*shagi)++;
        x_hord = a - fa * (b - a) / (fb - fa);
        x_kas = b - fb / fpb;

        a = x_hord;
        b = x_kas;

        fa = vichislit_funkciyu(a);
        fb = vichislit_funkciyu(b);
        fpb = vichislit_proizvodnuyu(b);

        if (fabs(fa) < TOCHNOST) {
            *koren = a;
            return;
        }
        if (fabs(fb) < TOCHNOST) {
            *koren = b;
            return;
        }
    }
    *koren = (a + b) / 2;
}

int main() {
    system("chcp 1251");
    double intervaly[20][2];
    int kolvo_intervalov;

    nayti_vse_intervaly(intervaly, &kolvo_intervalov);

    printf("Найдено интервалов с корнями: %d\n", kolvo_intervalov);
    for (int i = 0; i < kolvo_intervalov; i++) {
        printf("%2d. [%7.4f, %7.4f]\n", i + 1, intervaly[i][0], intervaly[i][1]);
    }
    while(1){
    int metod;
    printf("\nВыберите метод:\n");
    printf("1 - Метод половинного деления\n");
    printf("2 - Метод хорд\n");
    printf("3 - Метод Ньютона (касательных)\n");
    printf("4 - Метод простых итераций\n");
    printf("5 - Метод секущих\n");
    printf("6 - Комбинированный метод хорд и касательных\n");
    printf("Ваш выбор: ");
    scanf_s("%d", &metod);
    printf("\nРезультаты поиска корней:\n");
    printf("----------------------------------------------------------------------\n");
    printf("| №  |     Интервал     |      Корень      |   f(корень)  | Итерации |\n");
    printf("----------------------------------------------------------------------\n");
    for (int i = 0; i < kolvo_intervalov; i++) {
        double koren;
        int iteracii;

        switch (metod) {
        case 1:
            metod_polovinnogo_deleniya(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
            break;
        case 2:
            metod_hord(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
            break;
        case 3:
            metod_newtona(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
            break;
        case 4:
            metod_prostyh_iteraciy(intervaly[i][0], intervaly[i][1], &koren, &iteracii,
                (intervaly[i][0] + intervaly[i][1]) < 0 ? 0 : 1);
            break;
        case 5:
            metod_sekushih(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
            break;
        case 6:
            metod_hord_i_kasatelnyh(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
            break;
        default:
            return 1;
        }
        if (!isnan(koren)) {
            printf("| %2d | [%6.3f, %6.3f] | %16.8f | %12.2e | %8d |\n",
                i + 1, intervaly[i][0], intervaly[i][1],
                koren, vichislit_funkciyu(koren), iteracii);
        }
    }
    printf("----------------------------------------------------------------------\n");
    }
    return 0;
}