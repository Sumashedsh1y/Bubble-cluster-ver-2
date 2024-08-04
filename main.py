import subprocess
exe_path = r"Z:\Программирование\Diplom\x64\Release\Diplom.exe"
number = 0
for i in range(45):
    with open(r"Z:\Программирование\Diplom\Diplom\count.txt","r") as file:
        number = int(file.read())
    subprocess.call(exe_path)
    number += 1
    with open(r"Z:\Программирование\Diplom\Diplom\count.txt","w") as file:
        file.write(str(number))