import sys
if (len(sys.argv)!=2):
	print("Bad argument")
	exit()

else:
	if(sys.argv[1].endswith('.netlist')):
		f = open(sys.argv[1],'r')

#to extract the circuit definition
flag = 0
circuit_def_array = []
lines = f.readlines()
for line in lines:
    if line == '.circuit\n' :
        flag = 1
    if line == '.end\n' or line == ".end":
        flag = 0
    if flag == 1:
        circuit_def_array.append(line)
#to avoid addiinf '.circuit' to the circuit def
circuit_def_array = circuit_def_array[1:]

non_comment_def =[]

for i in circuit_def_array:
    non_comment_def.append(i.split(' #')[0])

stripped_array = []

for i in non_comment_def:
	stripped_array.append(i.rstrip())

list.reverse(stripped_array)

final_array= []

for i in stripped_array:
	list_bla = i.split(" ")
	list.reverse(list_bla)
	final_array.append(list_bla)

final_string = []

for i in final_array:
	string = ' '.join(i)
	final_string.append(string)
for i in final_string:
	print(i)
