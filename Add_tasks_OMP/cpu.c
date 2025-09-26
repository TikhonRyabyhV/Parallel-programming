#include <stdio.h>
#include <omp.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

int main() {
    int omp_procs = omp_get_num_procs();
    printf("omp_get_num_procs(): %d\n", omp_procs);

#ifdef _WIN32
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    printf("Logical processors (Windows): %d\n", si.dwNumberOfProcessors);
#else
    long logical = sysconf(_SC_NPROCESSORS_ONLN);
    printf("Logical processors (sysconf): %ld\n", logical);
#endif

    // Для физических ядер (пример для Linux; требует парсинга /proc/cpuinfo)
    // Здесь упрощённо: предполагаем, что OMP возвращает logical
    printf("Ожидаемые физические ядра: (требует ручной проверки через lscpu или /proc/cpuinfo)\n");

    return 0;
}
