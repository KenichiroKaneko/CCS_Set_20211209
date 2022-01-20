import os

# 100MBは2の20乗
file_limit = 2**(10 + 10) * 100
print(100*1000*1000)
print(file_limit)

def my_func(dir_path):

    l = os.listdir(path = dir_path)

    for i in range(len(l)):
        f = dir_path + '/' + l[i]
        if os.path.isdir(f):
            my_func(f)
        else:
            print(f, os.path.getsize(f))
            if os.path.getsize(f) > file_limit:
                print('   ↑でかい')

    return

my_func('.')