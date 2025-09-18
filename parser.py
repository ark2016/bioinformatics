
def simple_parse(filename):
	"""
	Простой парсер файлов с геномами
	gets: имя файла или путь к нему
	returns: два массива геномы и подписи к ним
	"""
	#TODO: парсить подписи к геномам
	names = []
	ret = []
	with open(filename, "r") as f:
		buff = ""
		for line in f:
			l = line.strip()
			if l[0] == ">": # special data
				names.append(l)
				if len(buff) != 0:
					ret.append(buff)
					buff = ""
			else:
				buff += l
	return ret, names
			


if __name__ == "__main__":
    print(simple_parse("uniprot-choline_esterase_reviewed-yes.fasta")[:2])
