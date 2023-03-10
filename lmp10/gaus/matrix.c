#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

matrix_t *
make_matrix (int rn, int cn)
{
  matrix_t *new_mat = malloc (sizeof *new_mat);
  if (new_mat == NULL)
    return NULL;
  if ((new_mat->e =
       malloc ((size_t) rn * (size_t) cn * sizeof *new_mat->e)) == NULL) {
    free (new_mat);
    return NULL;
  }
  new_mat->rn = rn;
  new_mat->cn = cn;
  memset (new_mat->e, 0, (size_t) (rn * (size_t) cn * sizeof *new_mat->e));
  return new_mat;
}

void
free_matrix (matrix_t * m)
{
  free (m->e);
  free (m);
}

void
put_entry_matrix (matrix_t * m, int i, int j, double val)
{
  if (i >= 0 && i < m->rn && j >= 0 && j <= m->cn)
    m->e[i * m->cn + j] = val;
}

void
add_to_entry_matrix (matrix_t * m, int i, int j, double val)
{
  if (i >= 0 && i < m->rn && j >= 0 && j <= m->cn)
    m->e[i * m->cn + j] += val;
}

double
get_entry_matrix (matrix_t * m, int i, int j )
{
  if (i >= 0 && i < m->rn && j >= 0 && j <= m->cn)
    return m->e[i * m->cn + j];
	else
		return -999;
}

matrix_t *
read_matrix (FILE * in)
{
  int rn, cn;
  int i, j;
  matrix_t *new_mat;
  if (fscanf (in, "%d %d", &rn, &cn) != 2)
    return NULL;

  if ((new_mat = make_matrix (rn, cn)) == NULL)
    return NULL;
  for (i = 0; i < rn; i++)
    for (j = 0; j < cn; j++)
      if (fscanf (in, "%lf", &new_mat->e[i * cn + j]) != 1) {
        free_matrix (new_mat);
        return NULL;
      }

  return new_mat;
}

void
write_matrix (matrix_t * m, FILE * out)
{
  int i, j;
  if (m == NULL) {
    fprintf (out, "Matrix is NULL\n");
    return;
  }

  fprintf (out, "%d %d\n", m->rn, m->cn);
  for (i = 0; i < m->rn; i++) {
    for (j = 0; j < m->cn - 1; j++)
      fprintf (out, "%8.5f ", m->e[i * m->cn + j]);
    fprintf (out, "%8.5f\n", m->e[i * m->cn + j]);
  }
}

matrix_t *
copy_matrix (matrix_t * s)
{
  matrix_t *d = NULL;
  if (s != NULL)
    d = make_matrix (s->rn, s->cn);
  if (d != NULL) {
    memcpy (d->e, s->e, s->rn * s->cn * sizeof *s->e);
    /* int i;
       for( i= 0; i < s->rn*s->cn; i++ )
       *(d->e+i)= *(s->e+i);
     */
    /* d->rn= s->rn; d->cn= s->cn; done in make_matrix */
  }
  return d;
}

matrix_t *
transpose_matrix (matrix_t * s)
{
  matrix_t *d = NULL;
  if (s != NULL)
    d = make_matrix (s->cn, s->rn);
  if (d != NULL) {
    int i, j;
    for (i = 0; i < s->rn; i++)
      for (j = 0; j < s->cn; j++)
        *(d->e + j * d->cn + i) = *(s->e + i * s->cn + j);
    /* d->rn= s->cn; d->cn= s->rn; done in make_matrix */
  }
  return d;
}

void
xchg_rows (matrix_t * m, int i, int j)
{
  if (m != NULL && i >= 0 && i < m->rn && j >= 0 && j < m->rn) {
    int k;
    double tmp;
    for (k = 0; k < m->cn; k++) {
      tmp = *(m->e + i * m->cn + k);
      *(m->e + i * m->cn + k) = *(m->e + j * m->cn + k);
      *(m->e + j * m->cn + k) = tmp;
    }
  }
}

void
xchg_cols (matrix_t * m, int i, int j)
{
  if (m != NULL && i >= 0 && i < m->cn && j >= 0 && j < m->cn) {
    int k;
    double tmp;
    for (k = 0; k < m->rn; k++) {
      tmp = *(m->e + k * m->cn + i);
      *(m->e + k * m->cn + i) = *(m->e + k * m->cn + j);
      *(m->e + k * m->cn + j) = tmp;
    }
  }
}

matrix_t *
mull_matrix (matrix_t * a, matrix_t * b)
{
  if (a == NULL || b == NULL || a->cn != b->rn)
    return NULL;
  else {
    matrix_t *c = make_matrix (a->rn, b->cn);
    int i, j, k;
    if (c == NULL)
      return NULL;
    for (i = 0; i < c->rn; i++) {
      for (j = 0; j < c->cn; j++) {
        double s = 0;
        for (k = 0; k < a->cn; k++)
          s += *(a->e + i * a->cn + k) * *(b->e + k * b->cn + j);
        *(c->e + i * c->cn + j) = s;
      }
    }
    return c;
  }
}

matrix_t *
ge_matrix (matrix_t * a)
{
  matrix_t *c = copy_matrix (a);
  if (c != NULL) {
    int i, j, k;
    int cn = c->cn;
    int rn = c->rn;
    double *e = c->e;
    for (k = 0; k < rn - 1; k++) {      /* eliminujemy (zerujemy) kolumn?? nr k */
      for (i = k + 1; i < rn; i++) {    /* p??tla po kolejnych
                                           wierszach poni??ej diagonalii k,k */
        double d = *(e + i * cn + k) / *(e + k * cn + k);
        for (j = k; j < cn; j++)
          *(e + i * cn + j) -= d * *(e + k * cn + j);
      }
    }
  }
  return c;
}

int
bs_matrix (matrix_t * a)
{
  if (a != NULL) {
    int r, c, k;
    int cn = a->cn;
    int rn = a->rn;
    double *e = a->e;

    for (k = rn; k < cn; k++) { /* p??tla po prawych stronach */
      for (r = rn - 1; r >= 0; r--) {   /* petla po niewiadomych */
        double rhs = *(e + r * cn + k); /* warto???? prawej strony */
        for (c = rn - 1; c > r; c--)    /* odejmujemy warto??ci ju?? wyznaczonych niewiadomych *
                                           pomno??one przed odpowiednie wsp????czynniki */
          rhs -= *(e + r * cn + c) * *(e + c * cn + k);
        *(e + r * cn + k) = rhs / *(e + r * cn + r);    /* nowa warto???? to prawa strona / el. diagonalny */
      }
    }
    return 0;
  }
  else
    return 1;
}


void add_matrix(matrix_t *a, matrix_t *b, double scalar, matrix_t **out)
// Dodawanie do siebie dwoch macierzy.
{
	for(int i = 0; i < (a->rn * a->cn); i++)
		(*out)->e[i] = a->e[i] + (b->e[i] * scalar);
}


void subtract_matrix(matrix_t *a, matrix_t *b, double scalar, matrix_t **out)
// Odejmomwanie macierzy b od a.
{
	for(int i = 0; i < (a->rn * a->cn); i++)
		(*out)->e[i] = a->e[i] - (b->e[i] * scalar);
}


matrix_t *scalar_matrix(double scalar, matrix_t *a)
// Mnozenie macierzy przez skalar.
{
	matrix_t *b = copy_matrix(a);
	for(int i = 0; i < (a->rn * a->cn); i++)
		b->e[i] *= scalar;
	return b;
}


matrix_t *flip_matrix(matrix_t *a)
{
	matrix_t *new = copy_matrix(a);

	new->rn = a->cn;
	new->cn = a->rn;

	return new;
}


int cgm_solver(matrix_t *m)
// Rozwiazywanie liniowego ukladu rownan metoda gradientow sprzezonych.
{
	int k, i, j, n, idx1, idx2;
	double a, tmpa;

	n = m->rn;

	// Rozdzielamy macierz m na A i b
	matrix_t *A = make_matrix(n, n);
	matrix_t *b = make_matrix(n, 1);

	// Przepisywanie liczb do macierzy A
	idx1 = 0;
	idx2 = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			A->e[idx1] = m->e[idx2];
			idx1++;
			idx2++;
		}
		idx2++;
	}

	// Przepisywanie liczb do macierzy b
	for (i = 0; i < n; i++)
	{
		b->e[i] = m->e[(i * (n + 1)) + n];
	}

	// Zerowanie macierzy rozwiazan
	matrix_t *x = make_matrix(n, 1);
	for (i = 0; i < n; i++)
		x->e[i] = 0.0;

	matrix_t *tmp1 = NULL;
	matrix_t *tmp2 = NULL;
	matrix_t *rT = NULL;
	matrix_t *pT = NULL;

	// Wektor residualny
	matrix_t *r = copy_matrix(b);
	matrix_t *p = copy_matrix(b);
	
	while(1)
	{
		// Zwalnianie macierzy, by uniknac wyciekow.
		if (tmp1 != NULL)
			free_matrix(tmp1);
		if (tmp2 != NULL)
			free_matrix(tmp2);
		if (rT != NULL)
			free_matrix(rT);
		if (pT != NULL)
			free_matrix(pT);
		
		rT = flip_matrix(r);
		
		// Czy konczymy
		tmp1 = mull_matrix(rT, r); 
		if (fabs(tmp1->e[0]) < 0.001)
			break;
		
		// Obliczanie alpha[k]
		free_matrix(tmp1);
		tmp1 = mull_matrix(rT, r); 
		tmpa = tmp1->e[0]; // tu jest licznik ulamka

		pT = flip_matrix(p);
		tmp1 = mull_matrix(pT, A);
		tmp2 = mull_matrix(tmp1, p); // tu jest mianownik ulamka

		a = tmpa / tmp2->e[0];

		// Obliczanie x[k]
		add_matrix(x, p, a, &x);
		
		// Obliczanie r[k+1]
		free_matrix(tmp1);
		tmp1 = mull_matrix(A, p);
		subtract_matrix(r, tmp1, a, &r);

		// Obliczanie beta[k]
		free_matrix(rT); free_matrix(tmp1);
		rT = flip_matrix(r);
		tmp1 = mull_matrix(rT, r); // tu jest licznik ulamka
		
		a = tmp1->e[0] / tmpa; // mianownik bety jest taki sam jak licznik dla alfy
	
		// Obliczanie p[k+1]
		add_matrix(r, p, a, &p);
	}
	
	// Przepisanie odpowiedzi do pierwotnego ukladu
	for (i = 0; i < n; i++)
	{
		m->e[(i * (n + 1)) + n] = x->e[i];
	}

	free_matrix(A);
	free_matrix(b);
	free_matrix(x);
	free_matrix(tmp1);
	free_matrix(rT);
	free(p);
	free(r);

	return 0;
}
