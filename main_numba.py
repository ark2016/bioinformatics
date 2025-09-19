import numpy as np
import numba as nb
from numba import jit, prange, njit
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import multiprocessing as mp
import functools
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

@njit(cache=True, parallel=True)
def compute_binom_coeffs_numba(size):
    binom = np.ones(size, dtype=np.int64)
    for i in range(1, size):
        binom[i] = (binom[i-1] * (size - i)) // i

    binom = binom % 26
    for i in prange(size):
        if binom[i] > 12:
            binom[i] -= 26

    return binom.astype(np.int32)

@njit(cache=True, parallel=True)
def compress_numba(string, compress_size, binom):
    if compress_size <= 1:
        return string

    result_len = len(string) - compress_size + 1
    result = np.zeros(result_len, dtype=np.int32)

    for i in prange(result_len):
        sum_val = 0
        for j in range(compress_size):
            sum_val += string[i + j] * binom[j]
        result[i] = sum_val % 26

    return result

@njit(cache=True)
def scalar_compare_numba(arr1, arr2):
    A = 12.5

    if len(arr1) != len(arr2):
        min_len = min(len(arr1), len(arr2))
        arr1 = arr1[:min_len]
        arr2 = arr2[:min_len]

    arr1_c = arr1.astype(np.float32) - A
    arr2_c = arr2.astype(np.float32) - A

    sc11 = np.dot(arr1_c, arr1_c)
    sc12 = np.dot(arr1_c, arr2_c)
    sc22 = np.dot(arr2_c, arr2_c)

    return sc12 / (np.sqrt(sc11) * np.sqrt(sc22))

@njit(cache=True)
def compare_single_numba(stra, strb, binom_cache):
    if len(stra) < len(strb):
        stra, strb = strb, stra

    delta = len(stra) - len(strb) + 1

    if delta > 1:
        if delta in binom_cache:
            binom = binom_cache[delta]
        else:
            binom = compute_binom_coeffs_numba(delta)
            binom_cache[delta] = binom

        stra = compress_numba(stra, delta, binom)

    return scalar_compare_numba(stra, strb)

class OptimizedGenomeComparator:
    def __init__(self):
        self.binom_cache = {}
        self._precompute_common_sizes()

    def _precompute_common_sizes(self):
        common_sizes = range(2, 100)
        for size in common_sizes:
            self.binom_cache[size] = compute_binom_coeffs_numba(size)

    def get_binom_coeffs(self, size):
        if size not in self.binom_cache:
            self.binom_cache[size] = compute_binom_coeffs_numba(size)
        return self.binom_cache[size]

    def compare_pair(self, stra, strb):
        if len(stra) < len(strb):
            stra, strb = strb, stra

        delta = len(stra) - len(strb) + 1

        if delta > 1:
            binom = self.get_binom_coeffs(delta)
            stra = compress_numba(stra, delta, binom)

        return scalar_compare_numba(stra, strb)
    
    def batch_compare_parallel(self, strings_batch, test_arr, n_workers=None):
        if n_workers is None:
            n_workers = mp.cpu_count()

        tasks = [(i, s, test_arr) for i, s in enumerate(strings_batch)]

        def process_task(task):
            idx, string_arr, test = task
            return idx, self.compare_pair(string_arr, test)

        results = np.zeros(len(strings_batch))

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            for idx, score in executor.map(process_task, tasks):
                results[idx] = score

        return results

class SIMDOptimizedComparator:

    def __init__(self):
        self.binom_cache = {}
        self._setup_vectorized_ops()

    def _setup_vectorized_ops(self):
        for size in range(2, 100):
            self.binom_cache[size] = self._compute_binom_vectorized(size)

    def _compute_binom_vectorized(self, size):
        binom = np.ones(size, dtype=np.int64)

        for i in range(1, size):
            binom[i] = (binom[i-1] * (size - i)) // i

        binom = binom % 26
        mask = binom > 12
        binom[mask] -= 26

        return binom.astype(np.int32)

    def compress_vectorized(self, string, compress_size):
        if compress_size <= 1:
            return string

        if compress_size not in self.binom_cache:
            self.binom_cache[compress_size] = self._compute_binom_vectorized(compress_size)
        binom = self.binom_cache[compress_size]

        from numpy.lib.stride_tricks import as_strided

        result_len = len(string) - compress_size + 1
        stride = string.strides[0]

        windows = as_strided(string,
                            shape=(result_len, compress_size),
                            strides=(stride, stride))

        result = np.dot(windows, binom) % 26

        return result.astype(np.int32)
    
    def batch_compare_vectorized(self, strings_batch, test_arr):
        results = np.zeros(len(strings_batch), dtype=np.float32)

        delta_groups = {}
        for i, string_arr in enumerate(strings_batch):
            if len(string_arr) < len(test_arr):
                comp_str, ref_str = test_arr, string_arr
            else:
                comp_str, ref_str = string_arr, test_arr

            delta = len(comp_str) - len(ref_str) + 1

            if delta not in delta_groups:
                delta_groups[delta] = []
            delta_groups[delta].append((i, comp_str, ref_str))

        for delta, group in delta_groups.items():
            if delta > 1:
                for idx, comp_str, ref_str in group:
                    compressed = self.compress_vectorized(comp_str, delta)
                    results[idx] = self._scalar_compare_vectorized(compressed, ref_str)
            else:
                for idx, comp_str, ref_str in group:
                    results[idx] = self._scalar_compare_vectorized(comp_str, ref_str)

        return results

    def _scalar_compare_vectorized(self, arr1, arr2):
        A = 12.5

        if len(arr1) != len(arr2):
            min_len = min(len(arr1), len(arr2))
            arr1 = arr1[:min_len]
            arr2 = arr2[:min_len]

        arr1_c = arr1.astype(np.float32) - A
        arr2_c = arr2.astype(np.float32) - A

        sc11 = np.einsum('i,i->', arr1_c, arr1_c, optimize=True)
        sc12 = np.einsum('i,i->', arr1_c, arr2_c, optimize=True)
        sc22 = np.einsum('i,i->', arr2_c, arr2_c, optimize=True)

        return sc12 / (np.sqrt(sc11) * np.sqrt(sc22))

class HybridOptimizedComparator:

    def __init__(self, n_threads=None):
        self.n_threads = n_threads or mp.cpu_count()
        self.simd_comp = SIMDOptimizedComparator()
        self.numba_comp = OptimizedGenomeComparator()
    
    def process_batch(self, strings_batch, test_arr, use_simd=True):
        batch_size = len(strings_batch)

        if batch_size < 100:
            return self.simd_comp.batch_compare_vectorized(strings_batch, test_arr)

        elif batch_size < 1000:
            return self._process_numba_threaded(strings_batch, test_arr)

        else:
            return self._process_multiprocess(strings_batch, test_arr)
    
    def _process_numba_threaded(self, strings_batch, test_arr):
        results = np.zeros(len(strings_batch), dtype=np.float32)

        def process_chunk(start_idx, end_idx):
            for i in range(start_idx, end_idx):
                results[i] = self.numba_comp.compare_pair(strings_batch[i], test_arr)

        chunk_size = len(strings_batch) // self.n_threads
        threads = []

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            for t in range(self.n_threads):
                start = t * chunk_size
                end = start + chunk_size if t < self.n_threads - 1 else len(strings_batch)
                threads.append(executor.submit(process_chunk, start, end))

            for thread in threads:
                thread.result()

        return results
    
    def _process_multiprocess(self, strings_batch, test_arr):
        chunk_size = len(strings_batch) // self.n_threads
        chunks = []

        for i in range(self.n_threads):
            start = i * chunk_size
            end = start + chunk_size if i < self.n_threads - 1 else len(strings_batch)
            chunks.append((strings_batch[start:end], test_arr, start))

        with ProcessPoolExecutor(max_workers=self.n_threads) as executor:
            futures = []
            for chunk_data, test, offset in chunks:
                future = executor.submit(self._process_chunk_mp, chunk_data, test, offset)
                futures.append(future)

            results = np.zeros(len(strings_batch), dtype=np.float32)
            for future in futures:
                chunk_results, offset = future.result()
                results[offset:offset+len(chunk_results)] = chunk_results

        return results
    
    def _process_chunk_mp(self, chunk, test_arr, offset):
        comp = SIMDOptimizedComparator()
        results = comp.batch_compare_vectorized(chunk, test_arr)
        return results, offset


class FastGenomeProcessor:

    def __init__(self, use_gpu=False, n_threads=None):
        self.use_gpu = use_gpu and self._check_gpu_available()
        self.n_threads = n_threads or mp.cpu_count()

        if self.use_gpu:
            print("Используется GPU ускорение")
            self.processor = None
        else:
            print(f"Используется CPU оптимизация ({self.n_threads} потоков)")
            self.processor = HybridOptimizedComparator(n_threads=self.n_threads)
    
    def _check_gpu_available(self):
        try:
            import cupy as cp
            return cp.cuda.runtime.getDeviceCount() > 0
        except ImportError:
            return False

    def process_file(self, filename, test_genome):
        from parser import simple_parse

        genomes, names = simple_parse(filename)

        genome_arrays = [self.string2arr(g) for g in genomes]
        test_arr = self.string2arr(test_genome)

        if self.use_gpu:
            results = None
        else:
            results = self.processor.process_batch(genome_arrays, test_arr)

        return results, names
    
    @staticmethod
    def string2arr(string):
        return np.array([ord(c) - ord('A') for c in string], dtype=np.int32)

    def find_similar_genomes(self, filename, test_genome, threshold=0.8):
        results, names = self.process_file(filename, test_genome)

        similar = []
        for i, score in enumerate(results):
            if score >= threshold:
                similar.append((names[i], score))

        similar.sort(key=lambda x: x[1], reverse=True)

        return similar

def compress_fast(string, compress_size, compress_rule=None):
    if isinstance(string, list):
        string = np.array(string, dtype=np.int32)

    comp = SIMDOptimizedComparator()
    return comp.compress_vectorized(string, compress_size)

def batch_compare_fast(strings_batch, test_arr, compress_funcs=None):
    comp = HybridOptimizedComparator()

    if not isinstance(test_arr, np.ndarray):
        test_arr = np.array(test_arr, dtype=np.int32)

    batch_arrays = []
    for s in strings_batch:
        if not isinstance(s, np.ndarray):
            batch_arrays.append(np.array(s, dtype=np.int32))
        else:
            batch_arrays.append(s)

    return comp.process_batch(batch_arrays, test_arr)

