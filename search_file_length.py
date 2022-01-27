import os
from sys import maxsize

# 100MBは2の20乗
file_limit = 2**(10 + 10) * 100
# print(100*1000*1000)
# print(file_limit)

def my_func(dir_path, max_name, max_size):

    l = os.listdir(path = dir_path)

    for i in range(len(l)):
        f = dir_path + '/' + l[i]
        if os.path.isdir(f):
            max_name, max_size = my_func(f,max_name, max_size)
        else:
            size_f = os.path.getsize(f)
            if size_f > max_size:
                max_size = size_f
                max_name = f

            if size_f > file_limit:
                print(f, size_f)
                print('   ↑でかい')
    
    return max_name, max_size

max_name, max_size = my_func('.', "", 0)
print(max_name, max_size)