import numpy as np
from parser import simple_parse
from compress import compress, haos_compress_binom_gen
from scalar_mul import scalar_compare

def compare(stra, strb):
	if len(stra) < len(strb):
		stra, strb = strb, stra

	delta = len(stra)-len(strb)+1

	stra = compress(stra, delta, haos_compress_binom_gen(delta))
	return scalar_compare(stra, strb)

def string2arr(string):
	"""Конвертирует строку в массив чисел (A=0, B=1, ...)"""
	return np.array([ord(c) - ord('A') for c in string], dtype=np.int32)

def batch_compare(strings_batch, test_arr, compress_funcs):
	"""Векторизованное сравнение батча строк с тестовой строкой"""
	batch_size = len(strings_batch)
	results = np.zeros(batch_size)

	for i, string_arr in enumerate(strings_batch):
		if len(string_arr) < len(test_arr):
			string_arr, test_arr_local = test_arr, string_arr
		else:
			string_arr, test_arr_local = string_arr, test_arr

		delta = len(string_arr) - len(test_arr_local) + 1

		if delta in compress_funcs:
			compress_func = compress_funcs[delta]
		else:
			compress_func = haos_compress_binom_gen(delta)
			compress_funcs[delta] = compress_func

		compressed = compress(string_arr, delta, compress_func)
		results[i] = scalar_compare(compressed, test_arr_local)

	return results
