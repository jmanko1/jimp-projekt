#include "cgmethod_solver.h"
#include <stdlib.h>

int cgmethod_solver(matrix_t * eqs)
{
	if (eqs != NULL)
	{
    		if (cgm_solver(eqs) == 0)
		{
      			return 0;
   	 	}
    		else
		{
      			return 1;
    		}
  	}
  	else
    		return 1;
}
