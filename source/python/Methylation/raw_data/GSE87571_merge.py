input_file_1 = open('GSE87571_matrix1of2.txt', 'r')
input_file_2 = open('GSE87571_matrix1of2.txt', 'r')
output_file = open('average_beta.txt', 'w')

header_1 = input_file_1.readline()
header_2 = input_file_2.readline()

header_1_list = header_1.split('\t')
record_names_1 = header_1_list[1:2:len(header_1_list)]


line_1 = header_1
line_2 = header_2
while line_1 != '' and line_2 != '':
    line_1 = input_file_1.readline()
    line_2 = input_file_2.readline()