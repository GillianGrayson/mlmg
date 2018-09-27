input_file = open('D:/Work/mlmg/data/GSE87571/average_beta.txt', 'r')
output_file = open('average_beta.txt', 'w')

bad_indexes = [115, 439, 716]

header = input_file.readline()
header_list = header.split('\t')
passed_header_list = [header_list[x] for x in range(0, len(header_list)) if x not in bad_indexes]
output_file.write('\t'.join(map(str, passed_header_list)))

line = header
line_num = 0
while line != '':
    line = input_file.readline()
    line_list = line.split('\t')
    line_list_passed = [line_list[x] for x in range(0, len(line_list)) if x not in bad_indexes]
    output_file.write('\t'.join(map(str, line_list_passed)))

    line_num += 1
    if divmod(line_num, 1000)[1] == 0:
        print(line_num)

input_file.close()
output_file.close()