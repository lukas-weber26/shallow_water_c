ocean: main.c
	cc ./main.c -O2 -o ocean -fopenmp -lglfw -lGL -lX11 -lpthread -lXrandr -lXi -lm -ldl -lGLEW

clean: 
	rm *.o output

