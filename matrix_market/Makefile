
all : write.exe

mmio.o : mmio.c
	gcc -c $^

write.o : write.c
	gcc -c $^

write.exe : write.o mmio.o
	gcc $^ -o $@