all: tpsdemo 2d-morph

tpsdemo: main.o
	g++ main.o -lglut -L/usr/lib/nvidia-331 -lGL -lGLU -o tpsdemo

main.o: main.cpp
	g++ -c -O2 main.cpp

2d-morph: 2d-morph.o
	g++ 2d-morph.o -o 2d-morph
	
2d-morph.o: 2d-morph.cpp
	g++ -c -O2 2d-morph.cpp

.PHONEY: clean
clean:
	rm main.o tpsdemo 2d-morph.o 2d-morph
