from parser import simple_parse
from compress import compress, haos_compress_binom_gen
from scalar_mul import scalar_compare


def compare(stra, strb):
	if len(stra) < len(strb):
		stra, strb = strb, stra
	
	delta = len(stra)-len(strb)+1
	
	stra = compress(stra, delta, haos_compress_binom_gen(delta))
	#print(stra, strb)
	return scalar_compare(stra, strb)

def main(filename, test_string):
	# Парсим файл
	strings, names = simple_parse(filename)
	
	# Определяем разницу
	"""
	mindelta = abs(len(test_string)-len(strings[0]))
	maxdelta = mindelta
	for i in range(1, len(strings)):
		delta = abs(len(test_string)-len(strings[1]))
		if delta > maxdelta:
			maxdelta = delta
		if delta < mindelta:
			mindelta = delta
	"""
	result = [0]*len(strings)
	# основной цикл
	for index in range(len(strings)):
		if index%10 == 0:
			print(index, "/", len(strings), end="\r")
		result[index] = compare(strings[index], test_string)
	
	
	indexed_list = [(index, value) for index, value in enumerate(result)]
	# Сортируем по значению в убывающем порядке
	sorted_list = sorted(indexed_list, key=lambda x: x[1], reverse=True)    
	for i in range(10):
		print(i+1, sorted_list[i], names[sorted_list[i][0]])


if __name__=="__main__":
	main("uniprot_sprot.fasta", "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")
	#main("uniprot-choline_esterase_reviewed-yes.fasta", "MHSKVTIICIRFLFWFLLLCMLIGKSHTEDDIIIATKNGKVRGMNLTVFGGTVTAFLGIPYAQPPLGRLRFKKPQSLTKWSDIWNATKYANSCCQNIDQSFPGFHGSEMWNPNTDLSEDCLYLNVWIPAPKPKNATVLIWIYGGGFQTGTSSLHVYDGKFLARVERVIVVSMNYRVGALGFLALPGNPEAPGNMGLFDQQLALQWVQKNIAAFGGNPKSVTLFGESAGAASVSLHLLSPGSHSLFTRAILQSGSFNAPWAVTSLYEARNRTLNLAKLTGCSRENETEIIKCLRNKDPQEILLNEAFVVPYGTPLSVNFGPTVDGDFLTDMPDILLELGQFKKTQILVGVNKDEGTAFLVYGAPGFSKDNNSIITRKEFQEGLKIFFPGVSEFGKESILFHYTDWVDDQRPENYREALGDVVGDYNFICPALEFTKKFSEWGNNAFFYYFEHRSSKLPWPEWMGVMHGYEIEFVFGLPLERRDNYTKAEEILSRSIVKRWANFAKYGNPNETQNNSTSWPVFKSTEQKYLTLNTESTRIMTKLRAQQCRFWTSFFPKVLEMTGNIDEAEWEWKAGFHRWNNYMMDWKNQFNDYTSKKESCVGL")
	#main("insulin.fasta", "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")
