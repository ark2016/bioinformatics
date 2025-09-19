
import numpy as np

def compress(string, compress_size, compress_rule):
	"""
	Производит сжатие строки с помощью клеточного автомата(КА) на определенную длинну.
	gets:
		string - строка для сжатия
		compress_size - размер сжатия
		compress_rule - правила КА(пример, функция haos_compress)
	returns: новая строка размера len(string)-compress_size
	"""
	if compress_size <= 1:
		return string

	string = np.asarray(string, dtype=np.int32)
	result_len = len(string) - compress_size + 1
	result = np.zeros(result_len, dtype=np.int32)

	# Векторизованная обработка КА
	for i in range(result_len):
		result[i] = compress_rule(string[i: i+compress_size])

	return result.tolist()


def haos_compress(vector):
	"""Тривиальная функция сжатия"""
	return np.sum(vector) % 26

def haos_compress_binom_gen(size):
	"""Генератор биномиальной функции сжатия с предвычисленными коэффициентами"""
	# Предвычисляем биномиальные коэффициенты
	binom = np.ones(size, dtype=np.int64)
	for i in range(1, size):
		binom[i] = binom[i-1] * (size - i) // i

	# Применяем модуло и центрирование
	binom = binom % 26
	binom = np.where(binom > 12, binom - 26, binom)

	def haos_compress_binom(vector):
		vector = np.asarray(vector)
		return np.dot(vector, binom) % 26

	return haos_compress_binom
	
	
if __name__ == "__main__": # микротесты
	SIZE = 30
	test = "ABCDSIDUFSIDUFHISDUHFISDUFNSIUDFNSDIUFHISDUFJI"
	
	def microtest(test_str, test_size):
		tmp1 = compress(test, test_size, haos_compress_binom_gen(test_size))
	
		tmp2 = test
		for i in range(test_size-1):
			tmp2 = compress(tmp2, 2, haos_compress)
		return tmp1, tmp2
	
	errors = 0
	for i in range(2, SIZE):
		a, b = microtest(test, i)
		if a != b:
			print("error", i, a, b)
			errors+=1
	if errors == 0:
		print("test success")
	else:
		print("found", errors, "errors")
	
