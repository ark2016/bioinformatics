

def compress(string, compress_size, compress_rule):
	"""
	Производит сжатие строки с помощью клеточного автомата(КА) на определенную длинну.
	gets:
		string - строка для сжатия
		compress_size - размер сжатия
		compress_rule - правила КА(пример, функция haos_compress)
	returns: новая строка размера len(string)-compress_size
	"""
	# TODO: нет проверки на безопасность данных

	# Преобразуем строку символов в массив чисел
	buff = list(map(lambda x: ord(x) - ord('A'), string))
	
	# Обрабатываем КА
	for i in range(len(string) - compress_size+1):
		buff[i] = compress_rule(buff[i: i+compress_size])
	
	# Преобразуем массив чисел в строку
	return "".join(list(map(lambda x: chr(x+ord('A')), buff[:-compress_size+1])))


def haos_compress(vector): # тривиальная функция сжатия
	return sum(vector) % 26	
	
def haos_compress_binom_gen(size):
	binom = [1]*size
	for i in range(1, size):
		binom[i] = binom[i-1]*(size-i)//i
	
	binom = list(map(lambda x: x%26-26 if x%26 > 13 else x%26, binom))
	
	def haos_compress_binom(vector):
		res = 0
		for i in range(size):
			res += vector[i]*binom[i]
		return res%26
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
		print("found", erors, "errors")
	
