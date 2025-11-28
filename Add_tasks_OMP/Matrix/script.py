import subprocess
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from math import log

# Настройки
EXECUTABLE = "matrix_mul_exe" if os.name == 'nt' else "./matrix_mul_exe"
C_SOURCE = "matrix_mul.c"

def compile_c_code():
    if not os.path.exists(EXECUTABLE) or (os.path.exists(C_SOURCE) and os.path.getmtime(C_SOURCE) > os.path.getmtime(EXECUTABLE)):
        print(f"Compilation {C_SOURCE}...")
        compile_cmd = [
            "gcc", "-fopenmp", "-mavx2", "-mfma", 
            C_SOURCE, "-o", "matrix_mul_exe", "-lm"
        ]
        try:
            subprocess.check_call(compile_cmd)
            print("Successful compilation.")
        except subprocess.CalledProcessError:
            print("Error due compilation.")
            sys.exit(1)

def run_benchmark():
    print(f"Benchmarks are running")
    try:
        result = subprocess.check_output([EXECUTABLE], universal_newlines=True)
    except subprocess.CalledProcessError as e:
        print(f"Execution error: {e}")
        sys.exit(1)

    data = {}
    headers = []
    
    for line in result.strip().split('\n'):
        parts = line.split()
        if not parts: continue
        if "Size" in parts:
            headers = parts
            for h in headers: data[h] = []
            continue
            
        if parts[0].isdigit() and headers:
            data["Size"].append(int(parts[0]))
            for i in range(1, len(headers)):
                val = parts[i]
                if val == "Skip" or val == "(FAIL)":
                    data[headers[i]].append(None)
                else:
                    try:
                        data[headers[i]].append(float(val))
                    except:
                        data[headers[i]].append(None)
    return data

def get_slope(N, times):
    valid_points = []
    for n, t in zip(N, times):
        if t is not None and t > 0:
            valid_points.append((log(n), log(t)))
    
    if len(valid_points) < 2: return 0
    
    if len(valid_points) > 3:
        valid_points = valid_points[-3:]
        
    X = np.array([p[0] for p in valid_points])
    Y = np.array([p[1] for p in valid_points])
    k, _ = np.polyfit(X, Y, 1)
    return k

def plot_group(title, filename, x_data, y_datasets, y_labels):
    plt.figure(figsize=(10, 6))
    plt.title(title)
    plt.xlabel("Размер матрицы N (log scale)")
    plt.ylabel("Время, с (log scale)")
    
    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.grid(True, which="both", ls="-", alpha=0.2)
    
    colors = ['r', 'g', 'b', 'm', 'k', 'c', 'y']
    markers = ['o', 's', '^', 'D', '*', 'v', 'x']
    
    for i, (y_data, label) in enumerate(zip(y_datasets, y_labels)):
        valid_x = [x for x, y in zip(x_data, y_data) if y is not None]
        valid_y = [y for x, y in zip(x_data, y_data) if y is not None]
        
        if not valid_x: continue
            
        k = get_slope(valid_x, valid_y)
        label_with_k = f"{label} (k={k:.2f})"
        
        plt.plot(valid_x, valid_y, label=label_with_k, 
                 color=colors[i % len(colors)], 
                 marker=markers[i % len(markers)])

    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Сохранен: {filename}")

def process_and_plot(data):
    N = data["Size"]
    
    # 1. Onethread opt.
    plot_group(
        "Однопоточные оптимизации (Single Thread)", 
        "bench_basic.png",
        N,
        [data.get("Naive_1"), data.get("Transp_1"), data.get("Block_1"), data.get("Strass_1"), data.get("SIMD_1")],
        ["Naive", "Transpose", "Blocked", "Strassen", "SIMD (AVX2)"]
    )


    # 2. Multithread opt.
    plot_group(
        "Многопоточные версии (8 Threads)", 
        "bench_multithreaded.png",
        N,
        [data.get("Naive_8"), data.get("Transp_8"), data.get("Block_8"), data.get("Strass_8"), data.get("SIMD_8")],
        ["Naive + p8", "Transp + p8", "Block + p8", "Strassen + p8", "SIMD + p8"]
    )
    
    # 3. The best opt.
    plot_group(
        "Итоговое сравнение", 
        "bench_final.png",
        N,
        [data.get("Naive_1"), data.get("Strass_8"), data.get("Fast_8")],
        ["Naive (Baseline)", "Strassen + p8", "Fast (Combined + p8)"]
    )

if __name__ == "__main__":
    compile_c_code()
    data = run_benchmark()
    if data:
        process_and_plot(data)
    else:
        print("No data for plots.")

