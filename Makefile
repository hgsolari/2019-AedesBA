aedes:	aedes22.c d-alimentacion22.c rk2.c algoritmos.c deriv22.c  rates22.c vm22.c binomial.c f-auxiliares22.c aedes22.h CONSTANTES22 reparto22_b.c
	gcc  -O3  aedes22.c -o aedes22 -lm

promedios: PromAdultos.c PromHuevos.c PromPesoComida.c PromFloats.c PromOviPositivas.c  PromOviPositivas.c PromHuevos_todos.c
	gcc -O3 -o adultos  PromAdultos.c -lm
	gcc -O3 -o huevos   PromHuevos.c -lm
	gcc -O3 -o huevosDiap   PromHuevosDiap.c -lm
	gcc -O3 -o huevosT   PromHuevos_todos.c -lm
	gcc -O3 -o comida   PromPesoComida.c -lm
	gcc -O3 -o ovi    PromOviPositivas.c -lm
	gcc -O3 -o PromTodas PromTodas.c -lm

