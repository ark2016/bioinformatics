import cupy as cp
import numpy as np

class GPUGenomeCompressor:

    def __init__(self, max_size=1000):
        self.max_size = max_size
        self.binom_cache = {}
        self._precompute_binomials()

    def _precompute_binomials(self):
        for size in range(2, min(100, self.max_size)):
            self.binom_cache[size] = self._compute_binom_coeffs(size)

    def _compute_binom_coeffs(self, size):
        binom = cp.ones(size, dtype=cp.int32)
        for i in range(1, size):
            binom[i] = binom[i-1] * (size - i) // i
        binom = binom % 26
        binom = cp.where(binom > 12, binom - 26, binom)
        return binom

    def compress_gpu(self, string, compress_size):
        if compress_size <= 1:
            return string

        if compress_size not in self.binom_cache:
            self.binom_cache[compress_size] = self._compute_binom_coeffs(compress_size)

        binom = self.binom_cache[compress_size]

        result = cp.convolve(string, binom[::-1], mode='valid') % 26

        return result

    def batch_compress(self, strings, compress_sizes):
        results = []

        size_groups = {}
        for i, (string, size) in enumerate(zip(strings, compress_sizes)):
            if size not in size_groups:
                size_groups[size] = []
            size_groups[size].append((i, string))

        for size, group in size_groups.items():
            if size not in self.binom_cache:
                self.binom_cache[size] = self._compute_binom_coeffs(size)

            binom = self.binom_cache[size]

            for idx, string in group:
                compressed = cp.convolve(string, binom[::-1], mode='valid') % 26
                results.append((idx, compressed))

        results.sort(key=lambda x: x[0])
        return [r[1] for r in results]

class GPUGenomeComparator:

    def __init__(self):
        self.compressor = GPUGenomeCompressor()
        self.A = 12.5

    def scalar_compare_gpu(self, arr1, arr2):
        if len(arr1) != len(arr2):
            min_len = min(len(arr1), len(arr2))
            arr1 = arr1[:min_len]
            arr2 = arr2[:min_len]

        arr1_centered = arr1 - self.A
        arr2_centered = arr2 - self.A

        sc11 = cp.dot(arr1_centered, arr1_centered)
        sc12 = cp.dot(arr1_centered, arr2_centered)
        sc22 = cp.dot(arr2_centered, arr2_centered)

        return float(sc12 / (cp.sqrt(sc11) * cp.sqrt(sc22)))

    def compare_single(self, stra, strb):
        if len(stra) < len(strb):
            stra, strb = strb, stra

        delta = len(stra) - len(strb) + 1

        if delta > 1:
            stra = self.compressor.compress_gpu(stra, delta)

        return self.scalar_compare_gpu(stra, strb)
    
    def batch_compare_gpu(self, strings_batch, test_string):
        test_gpu = cp.asarray(test_string)
        batch_gpu = [cp.asarray(s) for s in strings_batch]

        results = cp.zeros(len(strings_batch), dtype=cp.float32)

        size_groups = {}
        for i, string_gpu in enumerate(batch_gpu):
            if len(string_gpu) < len(test_gpu):
                string_to_compress = test_gpu
                reference = string_gpu
                swap = True
            else:
                string_to_compress = string_gpu
                reference = test_gpu
                swap = False

            delta = len(string_to_compress) - len(reference) + 1

            if delta not in size_groups:
                size_groups[delta] = []
            size_groups[delta].append((i, string_to_compress, reference, swap))

        for delta, group in size_groups.items():
            if delta > 1:
                strings_to_compress = [item[1] if not item[3] else item[2] for item in group]
                compressed = self.compressor.batch_compress(strings_to_compress, [delta] * len(group))

                for (idx, _, reference, swap), comp in zip(group, compressed):
                    if swap:
                        results[idx] = self.scalar_compare_gpu(reference, comp)
                    else:
                        results[idx] = self.scalar_compare_gpu(comp, reference)
            else:
                for idx, str_to_cmp, reference, _ in group:
                    results[idx] = self.scalar_compare_gpu(str_to_cmp, reference)

        return cp.asnumpy(results)

def string2arr_gpu(string):
    arr = np.array([ord(c) - ord('A') for c in string], dtype=np.int32)
    return cp.asarray(arr)

def batch_process_genomes(genome_strings, test_genome, batch_size=1000):
    comparator = GPUGenomeComparator()
    test_arr = np.array([ord(c) - ord('A') for c in test_genome], dtype=np.int32)

    all_results = []

    for i in range(0, len(genome_strings), batch_size):
        batch = genome_strings[i:i+batch_size]
        batch_arrays = [np.array([ord(c) - ord('A') for c in s], dtype=np.int32)
                       for s in batch]

        batch_results = comparator.batch_compare_gpu(batch_arrays, test_arr)
        all_results.append(batch_results)

    return np.concatenate(all_results)

class MultiStreamGPUComparator:

    def __init__(self, n_streams=4):
        self.n_streams = n_streams
        self.streams = [cp.cuda.Stream() for _ in range(n_streams)]
        self.comparators = [GPUGenomeComparator() for _ in range(n_streams)]

    def process_parallel(self, genome_strings, test_genome):
        test_arr = np.array([ord(c) - ord('A') for c in test_genome], dtype=np.int32)

        chunk_size = len(genome_strings) // self.n_streams
        chunks = []
        for i in range(self.n_streams):
            start_idx = i * chunk_size
            end_idx = start_idx + chunk_size if i < self.n_streams - 1 else len(genome_strings)
            chunks.append(genome_strings[start_idx:end_idx])

        results = []

        for i, (chunk, stream, comparator) in enumerate(zip(chunks, self.streams, self.comparators)):
            with stream:
                chunk_arrays = [np.array([ord(c) - ord('A') for c in s], dtype=np.int32)
                              for s in chunk]
                chunk_results = comparator.batch_compare_gpu(chunk_arrays, test_arr)
                results.append(chunk_results)

        for stream in self.streams:
            stream.synchronize()

        return np.concatenate(results)

def benchmark_comparison():
    import time

    np.random.seed(42)
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    def generate_random_genome(length):
        return ''.join(np.random.choice(list(alphabet), length))

    n_genomes = 10000
    genome_lengths = np.random.randint(100, 500, n_genomes)
    genomes = [generate_random_genome(l) for l in genome_lengths]
    test_genome = generate_random_genome(300)

    print(f"Тестирование на {n_genomes} геномах...")

    start = time.time()
    gpu_results = batch_process_genomes(genomes, test_genome, batch_size=1000)
    gpu_time = time.time() - start
    print(f"GPU время: {gpu_time:.3f} сек")

    start = time.time()
    multi_comparator = MultiStreamGPUComparator(n_streams=4)
    multi_results = multi_comparator.process_parallel(genomes, test_genome)
    multi_time = time.time() - start
    print(f"Multi-stream GPU время: {multi_time:.3f} сек")

    return gpu_results, multi_results

if __name__ == "__main__":
    if cp.cuda.runtime.getDeviceCount() > 0:
        print(f"Обнаружено GPU устройств: {cp.cuda.runtime.getDeviceCount()}")
        device = cp.cuda.Device()
        print(f"Используется: {device.name}")
        print(f"Память: {device.mem_info[1] / 1024**3:.2f} GB")

        benchmark_comparison()
    else:
        print("GPU не обнаружен. Установите CUDA и CuPy для GPU-ускорения.")