import matplotlib.pyplot as plt
import re
import pandas as pd
import numpy as np

def parse_data(filename="output_data.txt"):
    """
    Парсит табличный вывод C-программы и возвращает DataFrame с данными.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    if not lines:
        print("Файл данных пуст.")
        return pd.DataFrame()

    header = lines[0].strip().split()
    data = []
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) == len(header):
            data.append([
                int(parts[0]),      # N_SIZE
                int(parts[1]),      # THREADS
                int(parts[2]),      # THRESHOLD
                float(parts[3]),    # SEQ_TIME
                float(parts[4]),    # PAR_TIME
                float(parts[5]),    # SPEEDUP
                float(parts[6]),    # EFFICIENCY
                bool(int(parts[7])) # IS_CORRECT
            ])
    
    df = pd.DataFrame(data, columns=header)
    return df

def plot_graphs(df):

    if df.empty:
        print("DataFrame пуст, графики не будут построены.")
        return

    plt.style.use('seaborn-v0_8-darkgrid')
    
    # --- Фиксируем N для всех графиков (берем максимальный N из данных) ---
    max_N = df['N_SIZE'].max()
    df_fixed_N = df[df['N_SIZE'] == max_N]

    if df_fixed_N.empty:
        print(f"Нет данных для N={max_N} для построения графиков.")
        return

    print(f"Графики строятся для фиксированного размера массива N = {max_N}")

    # --- График 1: Время выполнения параллельной версии vs. Порог (кривые для каждого числа потоков) ---
    plt.figure(figsize=(12, 6))
    
    # Сначала данные для последовательной версии (THRESHOLD=0, THREADS=0)
    seq_data_fixed_N = df_fixed_N[df_fixed_N['THREADS'] == 0].iloc[0]
    seq_time = seq_data_fixed_N['SEQ_TIME']
    plt.axhline(y=seq_time, color='gray', linestyle='--', label='Последовательное время')

    # Затем данные для параллельной версии для каждого числа потоков
    for num_threads in range(8,9):# sorted(df_fixed_N['THREADS'].unique()):
        if num_threads == 0: continue # Пропускаем последовательную версию
        
        subset = df_fixed_N[(df_fixed_N['THREADS'] == num_threads) & (df_fixed_N['THREADS'] > 0)]
        if not subset.empty:
            plt.plot(subset['THRESHOLD'], subset['PAR_TIME'], marker='o', label=f'{num_threads} потоков')
    
    plt.xlabel('Порог последовательной обработки (Threshold)')
    plt.ylabel('Время выполнения (секунды)')
    plt.title(f'Время выполнения параллельной версии vs. Порог (N={max_N})')
    plt.ylim(0.3, 0.7)
    plt.xscale('log', base=2)
    #plt.xticks(df_fixed_N['THRESHOLD'].unique(), labels=[str(t) for t in df_fixed_N['THRESHOLD'].unique()])
    plt.legend()
    plt.grid(True, which="both", ls="--", c='0.7')
    plt.tight_layout()
    plt.savefig(f'timeVSthreshold.png')
    plt.show()

    # --- График 2: Ускорение vs. Число потоков (кривые для каждого значения Threshold) ---
    plt.figure(figsize=(12, 6))
    
    # Добавляем линию идеального ускорения
    max_threads_val = df_fixed_N['THREADS'].max()
    if max_threads_val > 1:
        ideal_speedup_x = np.array([t for t in df_fixed_N['THREADS'].unique() if t > 0])
        ideal_speedup_y = ideal_speedup_x
        plt.plot(ideal_speedup_x, ideal_speedup_y, linestyle=':', color='red', label='Идеальное ускорение')

    for threshold_val in 16, 64, 256, 1024: # sorted(df_fixed_N['THRESHOLD'].unique()):
        if threshold_val == 0: continue # Пропускаем THRESHOLD=0, если он был для последовательной версии
        
        subset = df_fixed_N[(df_fixed_N['THRESHOLD'] == threshold_val) & (df_fixed_N['THREADS'] > 0)]
        if not subset.empty:
            # Сортируем по количеству потоков для правильного отображения линии
            subset = subset.sort_values(by='THREADS')
            plt.plot(subset['THREADS'], subset['SPEEDUP'], marker='o', label=f'Threshold={threshold_val}')
    
    plt.xlabel('Число потоков')
    plt.ylabel('Ускорение (Speedup)')
    plt.title(f'Ускорение vs. Число потоков (N={max_N})')
    plt.xticks(df_fixed_N['THREADS'].unique())
    plt.ylim(1.0, 5.0)
    plt.legend()
    plt.grid(True, which="both", ls="--", c='0.7')
    plt.tight_layout()
    plt.savefig(f'speedVSthreads.png')
    plt.show()

    # --- График 3: Эффективность vs. Число потоков (кривые для каждого значения Threshold) ---
    plt.figure(figsize=(12, 6))
    
    # Добавляем линию идеальной эффективности (всегда 1.0)
    plt.axhline(y=1.0, color='gray', linestyle=':', label='Идеальная эффективность')

    for threshold_val in 16, 64, 256, 1024: # sorted(df_fixed_N['THRESHOLD'].unique()):
        if threshold_val == 0: continue # Пропускаем THRESHOLD=0
        
        subset = df_fixed_N[(df_fixed_N['THRESHOLD'] == threshold_val) & (df_fixed_N['THREADS'] > 0)]
        if not subset.empty:
            subset = subset.sort_values(by='THREADS')
            plt.plot(subset['THREADS'], subset['EFFICIENCY'], marker='o', label=f'Threshold={threshold_val}')
    
    plt.xlabel('Число потоков')
    plt.ylabel('Эффективность (Efficiency)')
    plt.title(f'Эффективность vs. Число потоков (N={max_N})')
    plt.xticks(df_fixed_N['THREADS'].unique())
    plt.legend()
    plt.grid(True, which="both", ls="--", c='0.7')
    plt.tight_layout()
    plt.savefig(f'effVSthreads.png')
    plt.show()

if __name__ == "__main__":
    df = parse_data()
    print("Полученные данные:")
    print(df.to_string())
    
    # Строим графики на всех данных
    plot_graphs(df)
