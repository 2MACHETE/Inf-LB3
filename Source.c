#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define TOCHNOST 1e-8       // точность вычислений
#define MAX_SHAGOV 1000        // максимальное число итераций
#define NAHALO_POISKA -10.0  // начало диапазона поиска
#define KONEC_POISKA 10.0   // конец диапазона поиска
#define SHAG_POISKA 0.1     // начальный шаг поиска

double vichislit_funkciyu(double x) {
    return 3 * x + x * x * x - pow(4, x) + 1;
}
double vichislit_proizvodnuyu(double x) {
    return 3 + 3 * x * x - log(4) * pow(4, x);
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
            x_tek = x_sled;
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

int main() {
    system("chcp 1251");
    double intervaly[20][2];
    int kolvo_intervalov;

    nayti_vse_intervaly(intervaly, &kolvo_intervalov);

    printf("Найдено интервалов с корнями: %d\n", kolvo_intervalov);
    for (int i = 0; i < kolvo_intervalov; i++) {
        printf("%2d. [%7.4f, %7.4f]\n", i + 1, intervaly[i][0], intervaly[i][1]);
    }

    int metod;
    printf("\nВыберите метод:\n");
    printf("1 - Метод половинного деления\n");
    printf("2 - Метод хорд\n");
    printf("Ваш выбор: ");
    scanf_s("%d", &metod);

    printf("\nРезультаты поиска корней:\n");
    printf("----------------------------------------------------------------------\n");
    printf("| №  |     Интервал     |      Корень      |   f(корень)  | Итерации |\n");
    printf("----------------------------------------------------------------------\n");

    for (int i = 0; i < kolvo_intervalov; i++) {
        double koren;
        int iteracii;

        if (metod == 1) {
            metod_polovinnogo_deleniya(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
        }
        else {
            metod_hord(intervaly[i][0], intervaly[i][1], &koren, &iteracii);
        }

        if (!isnan(koren)) {
            printf("| %2d | [%6.3f, %6.3f] | %16.8f | %12.2e | %8d |\n",
                i + 1, intervaly[i][0], intervaly[i][1],
                koren, vichislit_funkciyu(koren), iteracii);
        }
    }
    printf("----------------------------------------------------------------------\n");

    return 0;
}
