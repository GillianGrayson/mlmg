input_file_1 = open('GSE87571_matrix1of2.txt', 'r')
input_file_2 = open('GSE87571_matrix1of2.txt', 'r')
output_file = open('average_beta.txt', 'w')

eps = 1e-16

header_1 = input_file_1.readline()
header_2 = input_file_2.readline()

header_1_list = header_1.split('\t')
header_2_list = header_2.split('\t')
record_names_1_raw = header_1_list[1::2]
record_names_2_raw = header_2_list[1::2]
record_names_1 = [int(x[1::]) for x in record_names_1_raw]
record_names_2 = [int(x[1::]) for x in record_names_2_raw]


line_1 = header_1
line_2 = header_2
while line_1 != '' and line_2 != '':
    line_1 = input_file_1.readline()
    line_2 = input_file_2.readline()
    line_1_list = line_1.split('\t')
    line_2_list = line_2.split('\t')
    cpg = line_1_list[0]
    line_1_list = line_1_list[1::]
    line_2_list = line_2_list[1::]

    line_1_first = [float(x) if x != 'NA' else 0 for x in line_1_list[0::2]]
    line_1_second = [float(x) if x != 'NA' else 0 for x in line_1_list[1::2]]
    vals_1_list = [line_1_first[x] if line_1_first[x] > eps else line_1_second[x] for x in range(0, len(line_1_first))]

    line_2_first = [float(x) if x != 'NA' else 0 for x in line_2_list[0::2]]
    line_2_second = [float(x) if x != 'NA' else 0 for x in line_2_list[1::2]]
    vals_2_list = [line_2_first[x] if line_2_first[x] > eps else line_2_second[x] for x in range(0, len(line_2_first))]



    a = 1