
import numpy as np

def scalar1(arr1, arr2):
	"""
	Производит тривиальное скалярное умножение двух векторов
	"""
	arr1 = np.asarray(arr1)
	arr2 = np.asarray(arr2)
	min_len = min(len(arr1), len(arr2))
	return np.dot(arr1[:min_len], arr2[:min_len])

def scalar2(arr1, arr2):
	"""
	Производит скалярное умножение, но вектора имеют 0 в середине гиперкуба.
	Нулевых векторов при таком подходе нет
	"""
	A = 12.5
	arr1 = np.asarray(arr1)
	arr2 = np.asarray(arr2)
	min_len = min(len(arr1), len(arr2))
	arr1_centered = arr1[:min_len] - A
	arr2_centered = arr2[:min_len] - A
	return np.dot(arr1_centered, arr2_centered)

def scalar_compare(arr1, arr2):
	"""
	Производим сравнение двух строк ОДИНАКОВОЙ ДЛИНЫ на их похожесть.
	returns: cos \in [-1,1]
	0 - ортогональны, 1 - одинаковые, -1 - противостоящие
	"""
	A = 12.5
	arr1 = np.asarray(arr1) - A
	arr2 = np.asarray(arr2) - A

	sc11 = np.dot(arr1, arr1)
	sc12 = np.dot(arr1, arr2)
	sc22 = np.dot(arr2, arr2)

	return sc12 / (np.sqrt(sc11) * np.sqrt(sc22))
