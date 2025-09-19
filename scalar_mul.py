

def scalar1(arr1, arr2):
	"""
	Производит тривиальное скалярное умножение двух векторов
	"""
	s = 0
	for i in range(min(len(arr1), len(arr2))):
		s += arr1[i]*arr2[i]
	return s
	
def scalar2(arr1,arr2):
	"""
	Производит скалярное умножение, но вектора имеют 0 в середине гиперкуба.
	Нулевых векторов при таком подходе нет
	"""
	A = 12.5
	s = 0
	for i in range(min(len(arr1), len(arr2))):
		s += (arr1[i] - A)*(arr2[i] - A)
	return s
	
def scalar_compare(str1, str2):
	"""
	Производим сравнение двух строк ОДИНАКОВОЙ ДЛИНЫ на их похожесть. 
	returns: cos \in [-1,1]
	0 - ортогональны, 1 - одинаковые, -1 - противостоящие
	"""
	
	# Преобразуем в массивы 
	arr1 = list(map(lambda x: ord(x) - ord('A'), str1))
	arr2 = list(map(lambda x: ord(x) - ord('A'), str2))
	
	# Производим скалярные умножения
	sc11 = scalar2(arr1, arr1)
	sc12 = scalar2(arr1, arr2)
	sc22 = scalar2(arr2, arr2)
	#print(sc11, sc12, sc22)
	
	# Вычисляем косинус угла между векторами
	return sc12/(sc11**0.5 * sc22**0.5)
