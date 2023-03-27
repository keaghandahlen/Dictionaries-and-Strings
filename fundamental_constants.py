# exercise 6.1
#using constants.txt

def read_data(filename):
    infile = open(filename, 'r').readlines()[2:-1]  # we open the file and are in 'r' mode (read mode).
    # we use readlines to start at line 2 and go to the end of the file
    constants = {}  # opens a dictionary that we will add stuff to below
    for line in infile:
        words = line.split()  # creates a word after every space
        constant = float(words[-2])  # where our number is in the constants.txt file
        if len(words[:-2]) == 3:  # some constants have 3 words and some have 2
            # we combine them into one word in the next couple lines
            object = words[0] + ' ' + words[1] + ' ' + words[2]
        elif len(words[:-2]) == 2:
            object = words[0] + ' ' + words[1]
        else:
            object = words[0]
        constants[object] = constant
    return constants


constants = read_data('constants.txt')
for object, constant in constants.items():
    print(f'{object}, {constant}')
