aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

aprox_cgm: main.o splines.o points.o aproksymator_na_bazie_cgm.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie_cgm.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

intrp_cgm: main.o splines.o points.o interpolator_cgm.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator_cgm.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h gaus/cgmethod_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

aproksymator_na_bazie_cgm.o: makespl.h points.h gaus/piv_ge_solver.h gaus/cgmethod_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c -o aproksymator_na_bazie_cgm.o -DCGM

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h gaus/cgmethod_solver.h
	$(CC) -I gaus -c interpolator.c

interpolator_cgm.o: makespl.h points.h gaus/piv_ge_solver.h gaus/cgmethod_solver.h
	$(CC) -I gaus -c interpolator.c -o interpolator_cgm.o -DCGM

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta
