#include "matrix.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
	FILE *in;
	if (argc > 1 && (in = fopen (argv[1], "r")) != NULL)
	{
    		matrix_t *m = read_matrix(in);
		printf("Wczytana macierz:\n");
        	write_matrix(m, stdout);

		if (cgm_solver(m) != 0)
		{
			printf("Blad rozwiaazywania ukladu.\n");
			return 1;
		}

		printf("\nUklad po przejsciu funkcji:\n");
		write_matrix(m, stdout);

		printf("\nRozwiazania ukladu:\n");
		for (int i = 0; i < m->rn; i++)
		{
			printf("\t%g\n", m->e[(i * m->cn) + m->cn - 1]);
		}

    	}	
    	return 0;
}
