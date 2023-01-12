libge: matrix.o pivot.o piv_ge_solver.o cgmethod_solver.o
	ar rvs libge.a $^

matrix.o: matrix.c matrix.h
pivot.o: pivot.c matrix.h
piv_ge_solver.o: piv_ge_solver.c piv_ge_solver.h matrix.h
cgmethod_solver.o: cgmethod_solver.c cgmethod_solver.h

.PHONY: clean

clean:
	-rm *.o libge.a
