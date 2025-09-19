"""
CUDA-ускоренная версия для сравнения биологических последовательностей
Требует: cupy-cuda11x или cupy-cuda12x
Установка: pip install cupy-cuda11x  # или cupy-cuda12x для CUDA 12
"""

import numpy as np
try:
    import cupy as cp
    CUDA_AVAILABLE = True
    print("CUDA доступна, используем GPU ускорение")
except ImportError:
    import numpy as cp
    CUDA_AVAILABLE = False
    print("CUDA недоступна, используем CPU версию")

from parser import simple_parse

def scalar_compare_cuda(arr1, arr2):
    """CUDA-ускоренная версия scalar_compare"""
    A = 12.5
    arr1 = cp.asarray(arr1) - A
    arr2 = cp.asarray(arr2) - A

    # Векторизованные скалярные произведения на GPU
    sc11 = cp.dot(arr1, arr1)
    sc12 = cp.dot(arr1, arr2)
    sc22 = cp.dot(arr2, arr2)

    # Косинус угла между векторами
    return sc12 / (cp.sqrt(sc11) * cp.sqrt(sc22))

def compress_cuda(string, compress_size, binom_coeffs):
    """CUDA-ускоренная версия compress"""
    if compress_size <= 1:
        return string

    string = cp.asarray(string, dtype=cp.int32)
    result_len = len(string) - compress_size + 1

    # Создаем матрицу всех сегментов для векторизованной обработки
    segments = cp.zeros((result_len, compress_size), dtype=cp.int32)
    for i in range(result_len):
        segments[i] = string[i:i+compress_size]

    # Векторизованное матричное умножение на GPU
    result = cp.dot(segments, binom_coeffs) % 26
    return result

def haos_compress_binom_gen_cuda(size):
    """Генератор биномиальных коэффициентов для CUDA"""
    # Предвычисляем на CPU с защитой от переполнения
    binom = np.ones(size, dtype=np.int64)

    # Вычисляем биномиальные коэффициенты с применением модулю на каждом шаге
    for i in range(1, size):
        # Применяем модуль сразу для предотвращения переполнения
        numerator = (binom[i-1] * (size - i)) % (26 * 1000000) 
        binom[i] = numerator // i

    # Финальное приведение к модулю 26
    binom = binom % 26
    binom = np.where(binom > 12, binom - 26, binom)

    # Переносим на GPU
    return cp.asarray(binom, dtype=cp.int32)

def string2arr_cuda(string):
    """Конвертирует строку в GPU массив"""
    arr = np.array([ord(c) - ord('A') for c in string], dtype=np.int32)
    return cp.asarray(arr)

def batch_compare_cuda(strings_batch, test_arr, compress_funcs_cache):
    """Векторизованное CUDA сравнение"""
    batch_size = len(strings_batch)
    results = cp.zeros(batch_size)

    for i, string_arr in enumerate(strings_batch):
        # Определяем порядок строк для сжатия
        if len(string_arr) < len(test_arr):
            longer_arr, shorter_arr = test_arr, string_arr
        else:
            longer_arr, shorter_arr = string_arr, test_arr

        delta = len(longer_arr) - len(shorter_arr) + 1

        # Используем кэшированные коэффициенты
        if delta not in compress_funcs_cache:
            compress_funcs_cache[delta] = haos_compress_binom_gen_cuda(delta)

        binom_coeffs = compress_funcs_cache[delta]

        # CUDA-ускоренное сжатие и сравнение
        compressed = compress_cuda(longer_arr, delta, binom_coeffs)
        results[i] = scalar_compare_cuda(compressed, shorter_arr)

    return results

def main_cuda(filename, test_string):
    """Главная CUDA-ускоренная функция"""
    print("CUDA-ускоренная версия сравнения последовательностей")
    print(f"Используем: {'GPU (CUDA)' if CUDA_AVAILABLE else 'CPU (NumPy fallback)'}")

    # Парсинг файла (на CPU)
    print("Парсинг файла...")
    strings, names = simple_parse(filename)
    strings = strings[:50000]
    print(f"Загружено {len(strings)} последовательностей")

    # Конвертация в GPU массивы
    print("Конвертация данных для GPU...")
    test_arr = string2arr_cuda(test_string)

    # Конвертируем все строки и переносим на GPU батчами
    string_arrays = []
    conversion_batch_size = 5000

    for i in range(0, len(strings), conversion_batch_size):
        batch_end = min(i + conversion_batch_size, len(strings))
        batch_strings = [string2arr_cuda(s) for s in strings[i:batch_end]]
        string_arrays.extend(batch_strings)
        if (i // conversion_batch_size) % 10 == 0:
            print(f"Конвертировано {batch_end}/{len(strings)} строк", end="\r")

    print(f"\nВсе данные загружены в {'GPU' if CUDA_AVAILABLE else 'CPU'}")

    # Кэш для биномиальных коэффициентов
    compress_funcs_cache = {}

    # Обработка батчами на GPU
    processing_batch_size = 2000 if CUDA_AVAILABLE else 1000
    total_batches = (len(string_arrays) + processing_batch_size - 1) // processing_batch_size

    print("Начинаем GPU обработку...")
    all_results = []

    for batch_idx in range(total_batches):
        start_idx = batch_idx * processing_batch_size
        end_idx = min(start_idx + processing_batch_size, len(string_arrays))

        if batch_idx % 5 == 0:
            print(f"Батч {batch_idx + 1}/{total_batches} ({start_idx}/{len(string_arrays)})", end="\r")

        batch_strings = string_arrays[start_idx:end_idx]
        batch_results = batch_compare_cuda(batch_strings, test_arr, compress_funcs_cache)

        # Перенос результатов обратно на CPU для накопления
        if CUDA_AVAILABLE:
            batch_results = cp.asnumpy(batch_results)

        all_results.extend(batch_results)

    print(f"\nОбработка завершена!")
    print(f"Использовано {len(compress_funcs_cache)} уникальных размеров сжатия")

    # Финальная сортировка на CPU
    indexed_list = list(enumerate(all_results))
    sorted_list = sorted(indexed_list, key=lambda x: x[1], reverse=True)

    # Вывод результатов
    print("\n" + "="*60)
    print("ТОП-10 НАИБОЛЕЕ ПОХОЖИХ ПОСЛЕДОВАТЕЛЬНОСТЕЙ:")
    print("="*60)

    for i in range(min(10, len(sorted_list))):
        idx, score = sorted_list[i]
        print(f"{i+1:2d}. Сходство: {score:.6f} | Индекс: {idx}")
        print(f"    {names[idx]}")
        print("-" * 60)

if __name__ == "__main__":
    # Проверка на GPU память при запуске
    if CUDA_AVAILABLE:
        mempool = cp.get_default_memory_pool()
        print(f"Доступная GPU память: {mempool.total_bytes() / 1024**2:.1f} MB")

    main_cuda("uniprot_sprot.fasta",
             "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")