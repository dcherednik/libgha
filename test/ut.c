#include <3rd/fctx/fct.h>

#include <sle.h>

static double eq_matrix_1[3][4] =
	{{ 2,	 1,	-1,	 8},
	{-3,	-1,	 2,	-11},
	{-2,	 1,	 2,	-3}};
static double eq_result_1[3] =
	{2.0, 3.0, -1.0};

static double eq_matrix_2[2][3] =
	{{ 2,	 1,	1},
	{4,	2,	2}};


FCT_BGN()
{
	FCT_SUITE_BGN(simple)
	{
		FCT_TEST_BGN(sle_ok_1)
		{
			double* result = calloc(3, sizeof(double));;
			int i, rv;
			rv = sle_solve(eq_matrix_1[0], 3, result);
			fct_chk_eq_int(rv, 0);
			for (i = 0; i < 3; i++)
				fct_chk_eq_dbl(result[i], eq_result_1[i]);
			free(result);
		}
		FCT_TEST_END();

		FCT_TEST_BGN(sle_not_ok_2)
		{
			double* result = calloc(2, sizeof(double));;
			int rv;
			rv = sle_solve(eq_matrix_2[0], 2, result);
			fct_chk_eq_int(rv, -1);
			free(result);
		}
		FCT_TEST_END();
	}
	FCT_SUITE_END();
}
FCT_END();
