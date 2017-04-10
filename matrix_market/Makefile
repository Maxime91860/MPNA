
all : write.exe read.exe

mmio.o : mmio.c
	gcc -c $^

write.o : write.c
	gcc -c $^

write.exe : write.o mmio.o
	gcc $^ -o $@

read.o : read.c
	gcc -c $^

read.exe : read.o mmio.o
	gcc $^ -o $@

clean : 
	rm -rf *.exe *.o 