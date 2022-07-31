/**
 *    @file  ops_eig_sol_gcg.c
 *   @brief  ����ֵ����� GCG 
 *
 *  ����ֵ����� GCG
 *
 *  @author  Yu Li, liyu@tjufe.edu.cn
 *
 *       Created:  2020/8/18
 *      Revision:  none
 */

#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include    <math.h>
#include    <memory.h>
#include    <assert.h>
#include    <time.h>

#include    "ops_eig_sol_gcg.h"
#include    "pase_matvec.h"
#include    "pase_app_gcge.h"

#define     DEBUG 0
#define     DEBUG_PASE 0
#define     TIME_GCG 0
#define     PRINT_FIRST_UNCONV 0

double* check;

typedef struct TimeGCG_ {
	double initX_time;
	double checkconv_time;
	double compP_time;
	double compRR_time;
	double rr_matW_time;
	double dsyevx_time;
	double compRV_time;
	double compW_time;
	double compX_time;            
	double linsol_time;
	double time_total;
} TimeGCG;

struct TimeGCG_ time_gcg = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};


static int sizeN, startN, endN;
static int sizeP, startP, endP;
static int sizeW, startW, endW;
static int sizeC, sizeX , sizeV, endX;


static void   **mv_ws[3]; 
static double *dbl_ws; 
static int    *int_ws;
static struct OPS_ *ops_gcg;
static struct GCGSolver_ *gcg_solver;

#if 0
static double tmp_sigma[200];
#endif

/* y = ( A+sigma B ) x 
 * Only for CG (A+sigma*B)y = (lambda+sigma) B x 
 * use z[s:end] as workspace, which is b or p in CG */ 
static void MatDotMultiVecShift(void **x, void **y, 
	int *start, int *end, void **z, int s, struct OPS_ *ops)
{
	void  *A = gcg_solver->A;
	void  *B = gcg_solver->B;

	extern int pase_app_gcge_flag;
	if(pase_app_gcge_flag)
	{
		PASE_Matrix paseA = (PASE_Matrix)A;
		PASE_Matrix paseB = (PASE_Matrix)B;
		A=paseA->A_H;
		B=paseB->A_H;
	}

	double sigma = gcg_solver->sigma;
	ops->MatDotMultiVec(A,x,y,start,end,ops);
//	ops->Printf("sigma=%f\n", sigma);
	if (sigma != 0.0) {
		if (B==NULL) {
#if 1
			ops->MultiVecAxpby(sigma,x,1.0,y,start,end,ops);
#else
			int ncols = end[0]-start[0], col;
			for (col = 0; col < ncols; ++col) {
				ops->MultiVecAxpby(tmp_sigma[col+start[1]],x,1.0,y,start,end,ops);
			}
#endif
		}		
		else {
			//void **z;
			//ops->MultiVecCreateByMat(&z,end[0]-start[0],A,ops);
			int start_tmp[2], end_tmp[2];
			start_tmp[0] = start[0]; end_tmp[0] = end[0];
			start_tmp[1] = s       ; end_tmp[1] = s+end[0]-start[0];
			ops->MatDotMultiVec(B,x,z,start_tmp,end_tmp,ops);
			start_tmp[0] = s       ; end_tmp[0] = s+end[0]-start[0];
			start_tmp[1] = start[1]; end_tmp[1] = end[1];
			ops->MultiVecAxpby(sigma,z,1.0,y,start_tmp,end_tmp,ops);
			//ops->MultiVecDestroy(&z,end[0]-start[0],ops);	
		}
	}
	return;
}

static void InitializeX(void **V, void **ritz_vec, void *B, int nevGiven)
{
#if TIME_GCG
    time_gcg.initX_time -= ops_gcg->GetWtime();
#endif
	int start[2], end[2];	
	start[0] = 0; end[0] = nevGiven;
	start[1] = 0; end[1] = nevGiven;
	ops_gcg->MultiVecAxpby(1.0,ritz_vec,0.0,V,start,end,ops_gcg);	
	ops_gcg->Printf("sizeX = %d, nevGiven = %d, %s\n",
		sizeX,nevGiven,gcg_solver->initX_orth_method);
	/* orth_dbl_ws begin from the end of ss_eval */
	double *orth_dbl_ws = gcg_solver->dbl_ws+gcg_solver->nevMax+2*gcg_solver->block_size;
	if (0 == strcmp("mgs", gcg_solver->initX_orth_method))
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->initX_orth_block_size,
				gcg_solver->initX_orth_max_reorth,
				gcg_solver->initX_orth_zero_tol,
				//ritz_vec,gcg_solver->dbl_ws,ops_gcg);
				ritz_vec,orth_dbl_ws,ops_gcg);
	else if (0 == strcmp("bgs", gcg_solver->initX_orth_method))
			MultiVecOrthSetup_BinaryGramSchmidt(
				gcg_solver->initX_orth_block_size,
				gcg_solver->initX_orth_max_reorth,
				gcg_solver->initX_orth_zero_tol,
				//ritz_vec,gcg_solver->dbl_ws,ops_gcg);
				ritz_vec,orth_dbl_ws,ops_gcg);
	else
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->initX_orth_block_size,
				gcg_solver->initX_orth_max_reorth,
				gcg_solver->initX_orth_zero_tol,
				//ritz_vec,gcg_solver->dbl_ws,ops_gcg);
				ritz_vec,orth_dbl_ws,ops_gcg);

	ops_gcg->MultiVecOrth(V,0,&nevGiven,B,ops_gcg);
	ops_gcg->MultiVecSetRandomValue(V,nevGiven,sizeX,ops_gcg);
	ops_gcg->MultiVecOrth(V,nevGiven,&endX,B,ops_gcg);
	assert(endX==sizeX);

#if TIME_GCG
    time_gcg.initX_time += ops_gcg->GetWtime();
#endif
	return;
}

static void ComputeRitzVec(void **ritz_vec, void **V, double *ss_evec)
{
#if TIME_GCG
    time_gcg.compRV_time -= ops_gcg->GetWtime();
#endif 
	int start[2], end[2]; double *coef;
	start[0] = startN; end[0] = endW;
	start[1] = startN; end[1] = endX;
	coef     = ss_evec;
	ops_gcg->MultiVecLinearComb(V,ritz_vec,0, 
			start,end,coef,sizeV-sizeC,NULL,0,ops_gcg);
#if TIME_GCG
    time_gcg.compRV_time += ops_gcg->GetWtime();
#endif
	return;
}

static int CheckConvergence(void *A, void *B, double *ss_eval, void **ritz_vec, 
	int numCheck, double *tol, int *offset)
{
#if TIME_GCG
    time_gcg.checkconv_time -= ops_gcg->GetWtime();
#endif
	int start[2], end[2], idx; double *inner_prod;
	int nevConv;
	start[0] = startN; end[0] = start[0]+numCheck;
	start[1] = 0     ; end[1] = numCheck;	
	ops_gcg->MatDotMultiVec(A,ritz_vec,mv_ws[0],start,end,ops_gcg);	
	ops_gcg->MatDotMultiVec(B,ritz_vec,mv_ws[1],start,end,ops_gcg);	
	/* lambda Bx */
	ops_gcg->MultiVecLinearComb(NULL,mv_ws[1],0,start,end,
			NULL,0,ss_eval+startN,1,ops_gcg);
	start[0] = 0     ; end[0] = numCheck;
	start[1] = 0     ; end[1] = numCheck;
	/* Ax - lambda Bx */
	ops_gcg->MultiVecAxpby(-1.0,mv_ws[1],1.0,mv_ws[0],start,end,ops_gcg);
	/* ��ʹ�� ss_evec ���� */
	inner_prod = dbl_ws+(sizeV-sizeC)*sizeW;
	ops_gcg->MultiVecInnerProd('D',mv_ws[0],mv_ws[0],0,
			start,end,inner_prod,1,ops_gcg);
	for (idx = 0; idx < numCheck; ++idx) {
		inner_prod[idx] = sqrt(inner_prod[idx]);
	}
	for (idx = 0; idx < numCheck; ++idx) {
		/* ���Բ��� �� ��Բ��� ��ֱ�С�� tol[0] �� tol[1] */
		if (fabs(ss_eval[startN+idx]) > tol[1]) {
			if (inner_prod[idx] > tol[0] || 
					inner_prod[idx] > fabs(ss_eval[startN+idx])*tol[1]) {
#if PRINT_FIRST_UNCONV
				ops_gcg->Printf("GCG: [%d] %6.14e (%6.4e, %6.4e)\n",
						startN+idx,ss_eval[startN+idx],
						inner_prod[idx], inner_prod[idx]/fabs(ss_eval[startN+idx]));
#endif
				break;
			}
		}
		else {
			if (inner_prod[idx] > tol[0]) {
#if PRINT_FIRST_UNCONV
				ops_gcg->Printf("GCG: [%d] %6.14e (%6.4e, %6.4e)\n",
						startN+idx,ss_eval[startN+idx],
						inner_prod[idx], inner_prod[idx]/fabs(ss_eval[startN+idx]));
#endif 
				break;
			}
		}
	}	
	for ( ; idx > 0; --idx) {
		/* ���һ������������ֵ���һ��������������ֵ�����ظ� */
		if ( fabs((ss_eval[startN+idx-1]-ss_eval[startN+idx])/ss_eval[startN+idx-1]) 
		      > gcg_solver->gapMin) {
			break;
		}
	}
	nevConv = sizeC+idx;
	
	/* offset[0] Ϊδ������ĸ���, offset[2n-1] <= idx < offset[2n]
	 * idx �ǲ������ı�� 1 <= n <= offset[0] */
	int state, num_unconv;
	/* 1 1 0 0 1 1 1 1 0 0 1 0 1 0 0 0 0 0 0 */
	offset[0] = 0; state = 1; num_unconv = 0; 
	for (idx = 0; idx < numCheck; ++idx) {
		/* ��һ���ǲ������� */
		if (inner_prod[idx] > tol[0] || 
				inner_prod[idx] > fabs(ss_eval[startN+idx])*tol[1]) {
			/* ��һ���������� */
			if (state) {
				offset[ offset[0]*2+1 ] = startN+idx;
				state = 0;
			}
			++num_unconv;
			if (num_unconv == sizeN) {
				offset[ offset[0]*2+2 ] = startN+idx+1;
				++offset[0];
				break;
			}
		}
		else {
			/* ��һ���ǲ������� */
			if (!state) {
				offset[ offset[0]*2+2 ] = startN+idx;
				++offset[0];
				state = 1;
			}
		}
	}
	if (num_unconv < sizeN) {
		if (state == 1) {
			offset[ offset[0]*2+1 ] = startN+numCheck;	
		}
		offset[ offset[0]*2+2 ] = startN+numCheck+sizeN-num_unconv;
		offset[ offset[0]*2+2 ] = offset[ offset[0]*2+2 ] < endX?
				offset[ offset[0]*2+2 ]:endX;
		assert(offset[ offset[0]*2+1 ]<offset[ offset[0]*2+2 ]);
		++offset[0];	
	}
	
#if TIME_GCG
    time_gcg.checkconv_time += ops_gcg->GetWtime();
#endif
	assert(offset[0]>0); 
	return nevConv;
}

static void ComputeP(void **V, double *ss_evec, int *offset)
{
#if TIME_GCG
    time_gcg.compP_time -= ops_gcg->GetWtime();
#endif
	int length, incx, incy, ldm, block_size;
	int nrows, idx, col, start[2], end[2]; 
	double *source, *destin, *mat, *coef;
	
	/* ���� n ���ֶ�Ӧ���� */
	block_size = 0;
	for (idx = 0; idx < offset[0]; ++idx) {
		length = (sizeV-sizeC)*(offset[idx*2+2]-offset[idx*2+1]);		
		source = ss_evec+(sizeV-sizeC)*(offset[idx*2+1]-sizeC) ; incx = 1;
		destin = ss_evec+(sizeV-sizeC)*(sizeX-sizeC+block_size); incy = 1;
		dcopy(&length,source,&incx,destin,&incy);		
		block_size += offset[idx*2+2]-offset[idx*2+1];
	}
	sizeP = block_size;
	/* ���� np ���� */
	for (idx = 0; idx < offset[0]; ++idx) {
		length = (offset[idx*2+2]-offset[idx*2+1]);
		destin = ss_evec+(sizeV-sizeC)*(sizeX-sizeC)+(offset[idx*2+1]-sizeC);
		for (col = 0; col < sizeP; ++col) {
			memset(destin,0,length*sizeof(double));
			destin += sizeV-sizeC;
		}			
	}
	
	
	/* С��ģ������ */
	mat    = ss_evec; 
	nrows  = sizeV-sizeC; ldm  = sizeV-sizeC ;
	startP = sizeX-sizeC; endP = startP+sizeP;
	double *orth_dbl_ws = ss_evec+ldm*endP;
	/* ss_diag ss_matA ss_evec ʣ�µĿռ� */
	if (0 == strcmp("bqr", gcg_solver->compP_orth_method)) {
		int length_orth_dbl_ws = gcg_solver->length_dbl_ws-(orth_dbl_ws - gcg_solver->dbl_ws); 
		ops_gcg->DenseMatOrth(mat,nrows,ldm,startP,&endP,
			gcg_solver->compP_orth_zero_tol,
			orth_dbl_ws,length_orth_dbl_ws,int_ws);		
	}
	else {
		LAPACKVEC lapack_vec_P, lapack_vec_ws;
		lapack_vec_P.data   = mat;
		lapack_vec_P.ldd    = ldm;
		lapack_vec_P.ncols  = endP;
		lapack_vec_P.nrows  = nrows;
		
		lapack_vec_ws.data  = orth_dbl_ws;
		lapack_vec_ws.ldd   = ldm;
		lapack_vec_ws.ncols = endP-startP;
		lapack_vec_ws.nrows = nrows;
		if (0 == strcmp("mgs", gcg_solver->compP_orth_method))
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->compP_orth_block_size,
				gcg_solver->compP_orth_max_reorth,
				gcg_solver->compP_orth_zero_tol,
				(void*)&lapack_vec_ws,orth_dbl_ws+ldm*(endP-startP),
				ops_gcg->lapack_ops);
		else if (0 == strcmp("bgs", gcg_solver->compP_orth_method))
			MultiVecOrthSetup_BinaryGramSchmidt(
				gcg_solver->compP_orth_block_size,
				gcg_solver->compP_orth_max_reorth,
				gcg_solver->compP_orth_zero_tol,
				(void*)&lapack_vec_ws,orth_dbl_ws+ldm*(endP-startP),
				ops_gcg->lapack_ops);
		else
			MultiVecOrthSetup_ModifiedGramSchmidt(
				gcg_solver->compP_orth_block_size,
				gcg_solver->compP_orth_max_reorth,
				gcg_solver->compP_orth_zero_tol,
				(void*)&lapack_vec_ws,orth_dbl_ws+ldm*(endP-startP),
				ops_gcg->lapack_ops);
					
		ops_gcg->lapack_ops->MultiVecOrth((void*)&lapack_vec_P,
			startP,&endP,NULL,ops_gcg->lapack_ops);		
	}
	startP += sizeC; endP += sizeC; sizeP = endP-startP;
	
	/* ���� P */
	start[0] = startN; end[0] = endW ;
	start[1] = 0     ; end[1] = sizeP;
	coef     = ss_evec+(sizeV-sizeC)*(sizeX-sizeC);
	ops_gcg->MultiVecLinearComb(V,mv_ws[0],0,start,end,
			coef,sizeV-sizeC,NULL,0,ops_gcg);
	start[0] = 0     ; end[0] = sizeP;
	start[1] = startP; end[1] = endP ;
	ops_gcg->MultiVecAxpby(1.0,mv_ws[0],0.0,V,start,end,ops_gcg);
	
#if TIME_GCG
    time_gcg.compP_time += ops_gcg->GetWtime();
#endif	
	return;	
}

static void ComputeX(void **V, void **ritz_vec)
{
#if TIME_GCG
    time_gcg.compX_time -= ops_gcg->GetWtime();
#endif
	int start[2], end[2];
	start[0] = startN; end[0] = endX;
	start[1] = startN; end[1] = endX;
	ops_gcg->MultiVecAxpby(1.0,ritz_vec,0.0,V,start,end,ops_gcg);
#if TIME_GCG
    time_gcg.compX_time += ops_gcg->GetWtime();
#endif
	return;
}

static void ComputeW(void **V, void *A, void *B,
	double *ss_eval, void **ritz_vec, int *offset)
{
#if TIME_GCG
	time_gcg.compW_time -= ops_gcg->GetWtime();
#endif	
	void **b = ritz_vec;
	int start[2], end[2], block_size, length, inc, idx;
	double *destin = dbl_ws;
	
	double sigma = 0.0;
	if (gcg_solver->compW_cg_auto_shift==1) {
		sigma = -ss_eval[sizeC]+((ss_eval[sizeC+1]-ss_eval[sizeC])*0.01);
	}
	gcg_solver->sigma = gcg_solver->compW_cg_shift+sigma; sigma = gcg_solver->sigma;

	/* ��֧�� ͬʱʹ�� �ⲿ����� �� shift */
	assert(gcg_solver->compW_cg_auto_shift==0 || gcg_solver->user_defined_multi_linear_solver == 0);
	/* initialize */
	block_size = 0; startW = endP; inc = 1; 
	for (idx = 0; idx < offset[0]; ++idx) {
		length   = offset[idx*2+2]-offset[idx*2+1];
		/* initialize x */
		start[0] = offset[idx*2+1]  ; end[0] = offset[idx*2+2];
		start[1] = startW+block_size; end[1] = start[1]+length;
		ops_gcg->MultiVecAxpby(1.0,ritz_vec,0.0,V,start,end,ops_gcg);

		/* set b, b = (lambda+sigma) Bx */
		start[0] = offset[idx*2+1]     ; end[0] = offset[idx*2+2];
		start[1] = offset[1]+block_size; end[1] = start[1]+length;
		ops_gcg->MatDotMultiVec(B,V,b,start,end,ops_gcg);

		
		int i;
		/* shift eigenvalues with sigma */
		for (i = start[0]; i < end[0]; ++i) ss_eval[i] += sigma;
		ops_gcg->MultiVecLinearComb(NULL,b,0,start,end,
				NULL,0,ss_eval+start[0],1,ops_gcg);			
		dcopy(&length,ss_eval+start[0],&inc,destin,&inc);
		/* recover eigenvalues */
		for (i = start[0]; i < end[0]; ++i) ss_eval[i] -= sigma;
		destin += length;
		block_size += length;
	}
	endW = startW+block_size;	
#if DEBUG_PASE
	extern bool if_compute_w;
	if_compute_w = true;
#endif	
	/* solve x */
	start[0] = offset[1]; end[0] = start[0]+block_size;
	start[1] = startW   ; end[1] = endW               ;
#if TIME_GCG
    	time_gcg.linsol_time -= ops_gcg->GetWtime();
#endif
	void(*lin_sol)(void*,void**,void**,int*,int*,struct OPS_*);
	void *ws;
	lin_sol = ops_gcg->MultiLinearSolver;
	ws      = ops_gcg->multi_linear_solver_workspace;
	/* b is set to (lambda+sigma) Bx */
	if (gcg_solver->user_defined_multi_linear_solver==2) {
	    ops_gcg->MultiLinearSolver(A,b,V,start,end,ops_gcg);
	}
#if TIME_GCG
    	time_gcg.linsol_time += ops_gcg->GetWtime();
#endif
	if (gcg_solver->user_defined_multi_linear_solver==0||
	    	gcg_solver->user_defined_multi_linear_solver==2) {	        
#if 1
		extern int A_pre;
		/* 20210628 A = sigma B + A */
		if (sigma!=0.0 && B!=NULL && ops_gcg->MatAxpby!=NULL) {
			ops_gcg->MatAxpby(sigma,B,1.0,A,ops_gcg);
			if(A_pre==0){
				MultiLinearSolverSetup_BlockPCG(
					gcg_solver->compW_cg_max_iter,
					gcg_solver->compW_cg_rate,
					gcg_solver->compW_cg_tol,
					gcg_solver->compW_cg_tol_type,
					mv_ws,dbl_ws,int_ws,NULL,NULL,ops_gcg);
			}
			else{
				PASE_Diag_BlockPCG_Setup(
					gcg_solver->compW_cg_max_iter,
					gcg_solver->compW_cg_rate,
					gcg_solver->compW_cg_tol,
					gcg_solver->compW_cg_tol_type,
					mv_ws,dbl_ws,int_ws,NULL,NULL,ops_gcg);
			}
		}
		else {
#endif
			if(A_pre==0){
				MultiLinearSolverSetup_BlockPCG(
					gcg_solver->compW_cg_max_iter,
					gcg_solver->compW_cg_rate,
					gcg_solver->compW_cg_tol,
					gcg_solver->compW_cg_tol_type,
					mv_ws,dbl_ws,int_ws,NULL,MatDotMultiVecShift,ops_gcg);
			}
			else{
				PASE_Diag_BlockPCG_Setup(
					gcg_solver->compW_cg_max_iter,
					gcg_solver->compW_cg_rate,
					gcg_solver->compW_cg_tol,
					gcg_solver->compW_cg_tol_type,
					mv_ws,dbl_ws,int_ws,NULL,MatDotMultiVecShift,ops_gcg);
			}
#if 1
		}
#endif
	}
#if TIME_GCG
    	time_gcg.linsol_time -= ops_gcg->GetWtime();
#endif
	ops_gcg->MultiLinearSolver(A,b,V,start,end,ops_gcg);
#if DEBUG_PASE
	if_compute_w = false;
#endif
#if 1
	/* 20210628 recover A */
	if (sigma!=0.0 && B!=NULL && ops_gcg->MatAxpby!=NULL) {
		/* A = -sigma B + A */
		ops_gcg->MatAxpby(-sigma,B,1.0,A,ops_gcg);
	}
#endif

#if TIME_GCG
    	time_gcg.linsol_time += ops_gcg->GetWtime();
#endif
	ops_gcg->MultiLinearSolver             = lin_sol;
	ops_gcg->multi_linear_solver_workspace = ws;
	/* orth W in V */
	if (0 == strcmp("mgs", gcg_solver->compW_orth_method))
		MultiVecOrthSetup_ModifiedGramSchmidt(
			gcg_solver->compW_orth_block_size,
			gcg_solver->compW_orth_max_reorth,
			gcg_solver->compW_orth_zero_tol,
			mv_ws[0],dbl_ws,ops_gcg);
	else if (0 == strcmp("bgs", gcg_solver->compW_orth_method))
		MultiVecOrthSetup_BinaryGramSchmidt(
			gcg_solver->compW_orth_block_size,
			gcg_solver->compW_orth_max_reorth,
			gcg_solver->compW_orth_zero_tol,
			mv_ws[0],dbl_ws,ops_gcg);
	else
		MultiVecOrthSetup_ModifiedGramSchmidt(
			gcg_solver->compW_orth_block_size,
			gcg_solver->compW_orth_max_reorth,
			gcg_solver->compW_orth_zero_tol,
			mv_ws[0],dbl_ws,ops_gcg);
	
	ops_gcg->MultiVecOrth(V,startW,&endW,B,ops_gcg);

	sizeW = endW-startW;

#if TIME_GCG
    	time_gcg.compW_time += ops_gcg->GetWtime();
#endif	
	return;
}

static void ComputeRayleighRitz(double *ss_matA, double *ss_eval, double *ss_evec, double tol,
		int nevConv, double *ss_diag, void *A, void **V)
{
#if TIME_GCG
    time_gcg.compRR_time -= ops_gcg->GetWtime();
#endif	
	int nrows, ncols, nrowsA, ncolsA, length, incx, incy, idx, start[2], end[2];
	double *source, *destin, alpha;
	if (sizeP>0) {
		/* ���� PtAP ���� */
		nrows  = sizeP      ; ncols  = sizeP      ;
		nrowsA = sizeV-sizeC; ncolsA = sizeV-sizeC;
		/* C = alpha*op(Q)*op(A)*op(P) + beta*C */
		/* dbl_ws: nrows*ncols+nrowA*ncols
		 *       <=(sizeV+sizeP)*sizeP */
		ops_gcg->DenseMatQtAP('L','S',nrowsA,ncolsA,nrows,ncols,
				1.0,ss_evec+(sizeV-sizeC)*(sizeX-sizeC),sizeV-sizeC, /* Q */
				    ss_matA                            ,sizeV-sizeC, /* A */
				    ss_evec+(sizeV-sizeC)*(sizeX-sizeC),sizeV-sizeC, /* P */
				0.0,dbl_ws                             ,nrows      , /* C */
				dbl_ws+nrows*ncols);		
	}	
	
	sizeV  = sizeX+sizeP+sizeW;
	startN = startN + (nevConv-sizeC);
	endN   = endN   + (nevConv-sizeC);
	endN   = (endN<endX)?endN:endX;

	sizeN  = endN - startN; 
	sizeC  = nevConv;

	/* ���� ss_mat ss_evec */
	ss_matA = ss_diag+(sizeV-sizeC);
	ss_evec = ss_matA+(sizeV-sizeC)*(sizeV-sizeC); 

#if TIME_GCG
    time_gcg.rr_matW_time -= ops_gcg->GetWtime();
#endif
	if (sizeW>0) {
		/* ���� VtAW ���� */
		start[0] = startN; end[0] = endW;
		start[1] = startW; end[1] = endW;
		destin = ss_matA+(sizeV-sizeC)*(sizeX+sizeP-sizeC);
		/* (endW-startN)*(endW-startW) �� double 
		 *               (endW-startW) �� ���� */
		ops_gcg->MultiVecQtAP('S','N',V,A,V,0,start,end,destin,sizeV-sizeC,
				mv_ws[0],ops_gcg);
		/* �Գƻ� */
		length = sizeX+sizeP-sizeC;
		source = ss_matA+(sizeV-sizeC)*(sizeX+sizeP-sizeC); incx = 1; 
		destin = ss_matA+(sizeX+sizeP-sizeC); incy = sizeV-sizeC;
		for (idx = 0; idx < sizeW; ++idx) {
			dcopy(&length,source,&incx,destin,&incy);
			source += sizeV-sizeC; destin += 1;
		}
	}
#if TIME_GCG
    time_gcg.rr_matW_time += ops_gcg->GetWtime();
#endif
	
	if (sizeX == sizeV) {		
		int block_size = gcg_solver->block_size;
		destin     = ss_matA;
		length     = sizeX-sizeC;
		block_size = block_size<length?block_size:length;
		start[0] = sizeC; end[0] = sizeX;
		start[1] = sizeC; end[1] = start[1]+block_size;
		while (length) {			
			ops_gcg->MultiVecQtAP('S','N',V,A,V,0,start,end,
				destin,sizeV-sizeC,mv_ws[0],ops_gcg);
			destin    += (sizeV-sizeC)*block_size;
			length    -= block_size;
			block_size = block_size<length?block_size:length;
			start[1] = end[1]; end[1] = start[1]+block_size;
		}	
	}
	else {
		/* ���� X P ����, ���� C ���� */
		length = sizeX+sizeP-sizeC;
		destin = ss_matA;
		for (idx = 0; idx < length; ++idx) {
			memset(destin,0,length*sizeof(double));
			destin += sizeV-sizeC;
		}
		/* ��ֵ X ���ֵĶԽ��� */
		length = sizeX-sizeC;
		source = ss_eval+sizeC; incx = 1              ; 
		destin = ss_matA      ; incy = (sizeV-sizeC)+1;
		dcopy(&length,source,&incx,destin,&incy);
		/* ���� PtAP ����*/
		length = sizeP;
		source = dbl_ws                                           ; incx = 1; 
		destin = ss_matA+(sizeV-sizeC)*(sizeX-sizeC)+(sizeX-sizeC); incy = 1;
		for (idx = 0; idx < length; ++idx) {
			dcopy(&length,source,&incx,destin,&incy);
			source += length; destin += sizeV-sizeC;
		}
	}
	
	/* ��¼�Խ��߲��� */
	length = sizeV-sizeC;
	source = ss_matA; incx = (sizeV-sizeC)+1; 
	destin = ss_diag; incy = 1              ;
	dcopy(&length,source,&incx,destin,&incy);
	
	/* �� ss_matA ���� shift */
	if (gcg_solver->compW_cg_shift != 0.0) {
		alpha = 1.0;
		length = sizeV-sizeC;
		source = &(gcg_solver->compW_cg_shift); incx = 0;
		destin = ss_matA             ; incy = (sizeV-sizeC)+1;
		daxpy(&length,&alpha,source,&incx,destin,&incy);
	}

	/* ����С��ģ����ֵ���� */
	char   JOBZ, RANGE, UPLO; 
	int    LDA, M, LDZ, INFO, N, LWORK, *IWORK, *IFAIL; 
	double ABSTOL, *AA, *W, *Z, *WORK;
	JOBZ   = 'V'        ; RANGE  = 'A'; UPLO  = 'U'        ;
	LDA    = sizeV-sizeC; ABSTOL = tol; LDZ   = sizeV-sizeC; 
	IWORK  = int_ws; INFO   = 0  ;
	/* ���ټ��� C ���� */
	N      = sizeV-sizeC; M = N;
	IFAIL  = int_ws+5*N; 
	AA     = ss_matA;
	W      = ss_eval+sizeC; 
	Z      = ss_evec; 
	WORK   = Z+LDZ*N;
	/* ss_diag ss_matA ss_evec ʣ�µĿռ� */
	LWORK  = gcg_solver->length_dbl_ws-(WORK - gcg_solver->dbl_ws); 


#if OPS_USE_MPI
	/* �� PAS ���� GCG ʱ, ��ʹ�ò�����ô��? 
	 * û��ϵ, PAS ��Ҫ��֤ÿ�����̶����������� 
	 * ͬʱ, �����ķ�������, ��������Ч�ʵ�����
	 * ����Ҫ����, ��֤, ÿ�����̵�����������ȫһ�� */
	int *displs;
	int sendcount, *recvcounts;
	double *recvbuf;
	int IL, IU; int rank, nproc;

	/* ÿ�ж�һ��, ������ֵ��������, ����ͨѶ */
	LDZ  = LDZ+1;
	/* �������������� C �Ĳ��� */
	Z    = ss_evec;	
	/* ���ù����ռ� */ 
	WORK = Z+LDZ*N; LWORK = LWORK-N;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	/* ��������ֵ */
	destin = ss_eval+sizeC;
	length = N;
	/* ÿ������10�� */
	if (gcg_solver->compRR_min_num <= 0) {
	   gcg_solver->compRR_min_num = N/(nproc+2)>10?N/(nproc+2):10;
	}
	displs = malloc((2*nproc+1)*sizeof(int)); /* ������Ҫ 2*nproc+1 */
	if (rank == 0) {
	   SplitDoubleArray(destin,length,nproc,
		 gcg_solver->compRR_min_gap,
		 gcg_solver->compRR_min_num,
		 displs,dbl_ws,int_ws);
	}
	MPI_Bcast(displs,nproc+1,MPI_INT,0,MPI_COMM_WORLD);
 	sendcount  = displs[rank+1]-displs[rank];
	recvcounts = displs+nproc+1;
	for (idx = 0; idx < nproc; ++idx) {
		recvcounts[idx] = displs[idx+1]-displs[idx];
	}
	RANGE = 'I';
	/* 1 <= IL <= IU <= N */
	IL = displs[rank]+1; IU = displs[rank+1]  ;
	M  = IU-IL+1;
	/* ��ͬ���� W Z ��ͬ */
	W += displs[rank]  ; Z += LDZ*displs[rank];	

#if TIME_GCG
    	time_gcg.dsyevx_time -= ops_gcg->GetWtime();
#endif
	if (sendcount > 0) {
		dsyevx(&JOBZ,&RANGE,&UPLO,&N,AA,&LDA,
				NULL,NULL,&IL,&IU,&ABSTOL,&M,
				W,Z,&LDZ,WORK,&LWORK,IWORK,IFAIL,&INFO);
		assert(M==IU-IL+1);
		assert(INFO==0);
	}
#if TIME_GCG
    	time_gcg.dsyevx_time += ops_gcg->GetWtime();
	//ops_gcg->Printf("dsyevx = %.2f\n",time_gcg.dsyevx_time);
#endif
	/* ������õ�������ֵ���Ƶ� Z �����һ�� */
	length  = sendcount;
	source  = W      ; incx    = 1  ;
	destin  = Z+LDZ-1; incy    = LDZ;
	dcopy(&length,source,&incx,destin,&incy);
	recvbuf = ss_evec;
	sendcount *= LDZ;
	for (idx = 0; idx < nproc; ++idx) {
		recvcounts[idx] *= LDZ;
		displs[idx+1]   *= LDZ;
	}
	/* ȫ�ۼ�������, ���ͺͽ��ն����������� */

	MPI_Allgatherv(MPI_IN_PLACE,sendcount,MPI_DOUBLE,
		recvbuf,recvcounts,displs,MPI_DOUBLE,MPI_COMM_WORLD);
	free(displs);
	/* �� Z �����һ�и��Ƹ�����ֵ */
	length = N;
	source = ss_evec+LDZ-1; incx = LDZ;
	destin = ss_eval+sizeC; incy = 1  ;
	dcopy(&length,source,&incx,destin,&incy);
	/* �ƶ��������� */
#if DEBUG
	ops_gcg->Printf("before memmove length = %d\n", length);
#endif
	length = N; destin = ss_evec; source = ss_evec; 
	for (idx = 0; idx < N; ++idx) {
		/* ��֤ source �ڱ�����֮ǰ
		 * ���ص�������ֽڿ����� destin �� */
		memmove(destin,source,length*sizeof(double));
		destin += N; source += LDZ;
	}

#else
#if TIME_GCG
    time_gcg.dsyevx_time -= ops_gcg->GetWtime();
#endif
	/* ��֤ ss_evec ��������һ�� */
	dsyevx(&JOBZ,&RANGE,&UPLO,&N,AA,&LDA,
			NULL,NULL,NULL,NULL,&ABSTOL,&M,
			W,Z,&LDZ,WORK,&LWORK,IWORK,IFAIL,&INFO);
	assert(INFO==0);
#if TIME_GCG
    time_gcg.dsyevx_time += ops_gcg->GetWtime();
#endif
#if DEBUG
	ops_gcg->Printf("dsyevx: N = %d, M = %d\n", N, M);
#endif
	assert(M==N);

#endif

	/* �ָ�ss_matA�Խ��߲��� */
	length = sizeV-sizeC;
	source = ss_diag; incx = 1              ;
	destin = ss_matA; incy = (sizeV-sizeC)+1;
	dcopy(&length,source,&incx,destin,&incy);

	/* �ظ�����ֵ W */
	if (gcg_solver->compW_cg_shift != 0.0) {
		alpha  = -1.0;
		length = sizeV-sizeC;
		source = &(gcg_solver->compW_cg_shift); incx = 0;
		destin = ss_eval+sizeC       ; incy = 1;
		daxpy(&length,&alpha,source,&incx,destin,&incy);
	}
	
#if TIME_GCG
    time_gcg.compRR_time += ops_gcg->GetWtime();
#endif
	return;
}

static void GCG(void *A, void *B, double *eval, void **evec,
		int nevGiven, int *nevConv, struct OPS_ *ops)
{	
	/* offsetW[0] ��ʾ�ж��ٸ���, 
	 * offsetW[1] <= idx < offsetW[2] ��δ�����ı�� */ 
	int *offsetP, *offsetW, *ptr_tmp;
	gcg_solver = (GCGSolver*)ops->eigen_solver_workspace;
	gcg_solver->A = A; gcg_solver->B = B; 
	gcg_solver->nevGiven = nevGiven;
	gcg_solver->nevConv  = *nevConv;	
	ops_gcg = ops;
	int    nevMax, multiMax, block_size, nevInit, nev0, nev;
	int    numIterMax, numIter, numCheck;
	void   **V, **ritz_vec;
	double *ss_matA, *ss_diag, *ss_eval, *ss_evec, *tol;
	int    start[2], end[2], idx; double *coef;

	nevInit    = gcg_solver->nevInit   ;
	nevMax     = gcg_solver->nevMax    ; 
	block_size = gcg_solver->block_size; 
	multiMax   = gcg_solver->multiMax  ; 
	/*  �����ռ���� nevInit nevMax block_size ���� */
	assert(nevInit >= nevGiven);
	assert(nevInit <= nevMax);
	assert(nevInit >= 3*block_size || nevInit==nevMax);
        assert(nevMax  >= *nevConv+block_size);
	assert(nevMax  <= *nevConv+nevInit);
	assert(multiMax<= block_size);
	/* ��ʼ������ sizeX == nevInit ������Ҫ����� sizeX = nevMax ҪС
	 * �����ĺô���, dsyevx_ �Ĺ�ģ��С, �� gcg ��������������, 
	 * ������ֵ������ķǳ���ʱ����Ч�� */

	numIterMax = gcg_solver->numIterMax; tol = gcg_solver->tol;
	/* ȫ�ֱ�����ʼ�� */
	sizeC  = 0    ; sizeN = block_size  ; 
	/* sizeX ��Ҫ���� nevGiven */
	sizeX  = nevInit; sizeP  = 0; sizeW = 0; 
	sizeV  = sizeX+sizeP+sizeW;
	startN = sizeC; endN  = startN+sizeN; endX  = sizeX;
	startP = endX ; endP  = startP+sizeP;
	startW = endP ; endW  = startW+sizeW;
	/* workspace */
	V        = gcg_solver->mv_ws[0]; ritz_vec = evec;
	mv_ws[0] = gcg_solver->mv_ws[1]; mv_ws[1] = gcg_solver->mv_ws[2];
	mv_ws[2] = gcg_solver->mv_ws[3];
	ss_eval  = gcg_solver->dbl_ws; 
	for (idx = 0; idx < (nevMax+2*block_size); ++idx) {
	   ss_eval[idx] = 1.0;
	}
	ss_diag  = ss_eval+(nevMax+2*block_size);
	ss_matA  = ss_diag+(sizeV-sizeC);
	ss_evec  = ss_matA+(sizeV-sizeC)*(sizeV-sizeC); 

	int distance = (nevMax +2*block_size)                   /* ss_eval */ 
	              +(nevInit+2*block_size)                   /* ss_diag */ 
			+(nevInit+2*block_size)*(nevInit+2*block_size)  /* ss_matA */
			+(nevInit+2*block_size)*(nevInit+1*block_size); /* ss_evec */ 
	/* dbl_ws ���� W �Ĳ��� */
	dbl_ws = gcg_solver->dbl_ws+distance;
	gcg_solver->length_dbl_ws = (nevMax+2*block_size)                 /* ss_eval */
	                +2*(nevInit+2*block_size)*(nevInit+2*block_size)  /* ss_matA ss_evec */   
	                +10*(nevInit+2*block_size)                        /* ss_diag WORK */
			        +nevMax*block_size;                               /* for orth */

#if 1
	offsetP = gcg_solver->int_ws;
	offsetW = offsetP + block_size+3; 
	int_ws  = offsetW + block_size+3;
#else	
	int_ws  = gcg_solver->int_ws;
	offsetP = int_ws  + 6*(nevInit+2*block_size);
	offsetW = offsetP + block_size+2;
#endif

#if TIME_GCG
	time_gcg.checkconv_time = 0.0;
	time_gcg.compP_time     = 0.0; 
	time_gcg.compRR_time    = 0.0; 
	time_gcg.compRV_time    = 0.0;
	time_gcg.compW_time     = 0.0;
	time_gcg.compX_time     = 0.0;
	time_gcg.rr_matW_time   = 0.0;
	time_gcg.dsyevx_time    = 0.0;
	time_gcg.initX_time     = 0.0;
	time_gcg.linsol_time    = 0.0;
#endif
	
	/* �� X �������ֵ�� B ������һ�� */
	InitializeX(V,ritz_vec,B,nevGiven);		

	ComputeRayleighRitz(ss_matA,ss_eval,ss_evec,
		gcg_solver->compRR_tol,0,ss_diag,A,V);	

	// for (idx = sizeV; idx < (nevMax+2*block_size); ++idx) {
	//    ss_eval[idx] = ss_eval[sizeV-1];
	// }
	/* ���� ss_mat ss_evec */
	ss_matA = ss_diag+(sizeV-sizeC);
	ss_evec = ss_matA+(sizeV-sizeC)*(sizeV-sizeC);


	ComputeRitzVec(ritz_vec,V,ss_evec);				
	
	*nevConv = (*nevConv)<nevMax?(*nevConv):nevMax;
	/* �û�ϣ�������������Ը��� */
	nev0 = *nevConv; *nevConv = 0; 
	/* ���������ﵽ nev �� P �� W ��������Ϊ X ���� */
	nev  = nevInit<nevMax?2*block_size:nev0;
	nev  = nev<nev0?nev:nev0;
	numIter = 0; /* numIter ȡ��ֵʱ, С�ڵ�����ĵ����������ж������� */
#if PRINT_FIRST_UNCONV
	ops_gcg->Printf("------------------------------\n");
	ops_gcg->Printf("numIter\tnevConv\n",numIter, *nevConv);		
#endif
	do {
		if (numIter <= 0) {
		   numCheck = 0;
		} 
		else {
		   numCheck = (startN+sizeN<endX)?(sizeN):(endX-startN);
		}
		numCheck = numCheck<gcg_solver->check_conv_max_num?numCheck:gcg_solver->check_conv_max_num;
		*nevConv = CheckConvergence(A,B,ss_eval,ritz_vec,numCheck,tol,offsetW);		
#if PRINT_FIRST_UNCONV
		ops_gcg->Printf("%d\t%d\n",numIter, *nevConv);		
#endif
		if (*nevConv >= nev) {
			if (*nevConv >= nev0) {
				break;
			}
			else {
				/* Update sizeX */
				nev   += sizeP+sizeW; 
				nev    = nev<nev0?nev:nev0;
				sizeX += sizeP+sizeW; 
				sizeX  = sizeX<nevMax?sizeX:nevMax;
				/* �� P �� W ����д�� ritz_vec */
				start[0] = startN; end[0] = endW ;
				start[1] = endX  ; end[1] = sizeX; 
				coef     = ss_evec+(sizeV-sizeC)*(endX-sizeC);
				ops_gcg->MultiVecLinearComb(V,ritz_vec,0, 
						start,end,coef,sizeV-sizeC,NULL,0,ops_gcg);
			
				sizeP  = 0; sizeW = 0; sizeV = sizeX;
				startP = endX ; endP = startP;
				startW = endP ; endW = startW; 
				endX   = sizeX; 

				endN   = startN+block_size;
				endN   = endN<endX?endN:endX;
				sizeN  = endN-startN;	

				numIterMax -= numIter; numIter = 0;				
			}
		}
		if (numIter == 0)	{
			sizeP = 0; startP = endX; endP = startP+sizeP;
		}
		else {
	    ComputeP(V,ss_evec,offsetP); /* update sizeP startP endP */
		}

		ComputeX(V,ritz_vec);

		ComputeW(V,A,B,ss_eval,ritz_vec,offsetW); /* update sizeW startW endW */
		ptr_tmp = offsetP; offsetP = offsetW; offsetW = ptr_tmp;
		
	
		/* ������ PtAP ���ֺ��ٸ��� sizeV */
		ComputeRayleighRitz(ss_matA,ss_eval,ss_evec,
			gcg_solver->compRR_tol,*nevConv,ss_diag,A,V); /* update sizeC startN endN sizeN */

		for (idx = sizeV; idx < (nevMax+2*block_size); ++idx) {
		   ss_eval[idx] = ss_eval[sizeV-1];
		}
		ss_matA = ss_diag+(sizeV-sizeC);
		ss_evec = ss_matA+(sizeV-sizeC)*(sizeV-sizeC);

		ComputeRitzVec(ritz_vec,V,ss_evec);
		
		++numIter;
	} while (numIter < numIterMax);
	
	gcg_solver->numIter = numIter+(gcg_solver->numIterMax-numIterMax);
	/* eval evec ���� sizeX �� */
	int inc = 1;
	dcopy(&sizeX,ss_eval,&inc,eval,&inc);
	
#if TIME_GCG
	ops_gcg->Printf("|--GCG----------------------------\n");
	time_gcg.time_total = time_gcg.checkconv_time
		+time_gcg.compP_time
		+time_gcg.compRR_time
		+time_gcg.compRV_time
		+time_gcg.compW_time
		+time_gcg.compX_time
		+time_gcg.initX_time;
	ops_gcg->Printf("|Total Time = %.2f, Avg Time per Iteration = %.2f\n",
		time_gcg.time_total,time_gcg.time_total/gcg_solver->numIter);
	ops_gcg->Printf("|checkconv   compP   compRR   (rr_matW   dsyexv)   compRV   compW   (linsol)   compX   initX\n");
	ops_gcg->Printf("|%.2f\t%.2f\t%.2f\t(%.2f\t%.2f)\t%.2f\t%.2f\t(%.2f)\t%.2f\t%.2f\n",
		time_gcg.checkconv_time,		
		time_gcg.compP_time,		
		time_gcg.compRR_time,
		time_gcg.rr_matW_time,
		time_gcg.dsyevx_time,		
		time_gcg.compRV_time,
		time_gcg.compW_time,
		time_gcg.linsol_time,
		time_gcg.compX_time,
		time_gcg.initX_time);	   	
	ops_gcg->Printf("|%.2f%%\t%.2f%%\t%.2f%%\t(%.2f%%\t%.2f%%)\t%.2f%%\t%.2f%%\t(%.2f%%)\t%.2f%%\t%.2f%%\n",
		time_gcg.checkconv_time/time_gcg.time_total*100,
		time_gcg.compP_time    /time_gcg.time_total*100,
		time_gcg.compRR_time   /time_gcg.time_total*100,
		time_gcg.rr_matW_time  /time_gcg.compRR_time*100,
		time_gcg.dsyevx_time   /time_gcg.compRR_time*100,
		time_gcg.compRV_time   /time_gcg.time_total*100,		
		time_gcg.compW_time    /time_gcg.time_total*100,
		time_gcg.linsol_time   /time_gcg.compW_time*100,
		time_gcg.compX_time    /time_gcg.time_total*100,
		time_gcg.initX_time    /time_gcg.time_total*100);
	ops_gcg->Printf("|--GCG----------------------------\n");
	time_gcg.checkconv_time = 0.0;
	time_gcg.compP_time     = 0.0; 
	time_gcg.compRR_time    = 0.0; 
	time_gcg.compRV_time    = 0.0;
	time_gcg.compW_time     = 0.0;
	time_gcg.compX_time     = 0.0;
	time_gcg.rr_matW_time   = 0.0;
	time_gcg.dsyevx_time    = 0.0;
	time_gcg.initX_time     = 0.0;
	time_gcg.linsol_time    = 0.0;
#endif
	
	return;
}

/* �趨 GCG �Ĺ����ռ� */
void EigenSolverSetup_GCG(
	int    multiMax, double gapMin , 
	int    nevInit , int    nevMax , int block_size, 
	double tol[2]  , int    numIterMax, 
	int user_defined_multi_linear_solver,
	void **mv_ws[4], double *dbl_ws, int *int_ws,	
	struct OPS_ *ops)
{
	static GCGSolver gcg_solver_static = {
		.nevMax     = 3 , .multiMax   = 2, .gapMin = 0.01, 
		.nevInit    = 3 , .nevGiven   = 0,
		.block_size = 1 , .numIterMax = 4, .user_defined_multi_linear_solver = 0,
		.mv_ws      = {}, .dbl_ws  = NULL, .int_ws = NULL,		
		/* �㷨�ڲ����� */		
		.initX_orth_method     = "mgs",
		.initX_orth_block_size = -1   ,
		.initX_orth_max_reorth = 1    ,
		.initX_orth_zero_tol   = 1e-14,
		.check_conv_max_num    = 15   ,	
		.compP_orth_method     = "mgs", 
		.compP_orth_block_size = -1   ,
		.compP_orth_max_reorth = 1    ,
		.compP_orth_zero_tol   = 1e-14,		
		.compW_orth_method     = "mgs",
		.compW_orth_block_size = -1   ,
		.compW_orth_max_reorth = 1    ,
		.compW_orth_zero_tol   = 1e-14,	
		.compW_cg_max_iter     = 40   ,
		.compW_cg_rate         = 1e-2 , 
		.compW_cg_tol          = 1e-8 ,
		.compW_cg_tol_type     = "abs",
		.compW_cg_auto_shift   = 0    ,	
		.compW_cg_shift        = 0.0  ,	
		.compW_cg_order        = 1    ,	
		.compRR_min_gap        = 0.01 ,
		.compRR_min_num        = -1   ,
		.compRR_tol            = 1e-16,
	};
		
	gcg_solver_static.multiMax   = multiMax;
	gcg_solver_static.gapMin     = gapMin;	 
	gcg_solver_static.nevInit    = nevInit;
	gcg_solver_static.nevMax     = nevMax;
	gcg_solver_static.block_size = block_size;
	gcg_solver_static.tol[0]     = tol[0];
	gcg_solver_static.tol[1]     = tol[1];
	gcg_solver_static.numIterMax = numIterMax;
	gcg_solver_static.mv_ws[0]   = mv_ws[0];
	gcg_solver_static.mv_ws[1]   = mv_ws[1];
	gcg_solver_static.mv_ws[2]   = mv_ws[2];
	gcg_solver_static.mv_ws[3]   = mv_ws[3];
	gcg_solver_static.dbl_ws     = dbl_ws;
 	gcg_solver_static.int_ws     = int_ws;
 	
 	gcg_solver_static.compRR_min_gap = gapMin;
 	gcg_solver_static.check_conv_max_num = block_size;
 	gcg_solver_static.user_defined_multi_linear_solver = user_defined_multi_linear_solver;
		
	ops->eigen_solver_workspace = (void *)(&gcg_solver_static);
	ops->EigenSolver            = GCG;
	return;	
}

void EigenSolverCreateWorkspace_GCG(
	int nevInit, int nevMax, int block_size, void *mat,
	void ***mv_ws, double **dbl_ws, int **int_ws, 
	struct OPS_ *ops)
{
	assert(mv_ws!=NULL);
	int sizeV = nevMax+2*block_size; 
	ops->MultiVecCreateByMat(&mv_ws[0],sizeV,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[0],0,sizeV,ops);
	ops->MultiVecCreateByMat(&mv_ws[1],block_size,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[1],0,block_size,ops);
	ops->MultiVecCreateByMat(&mv_ws[2],block_size,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[2],0,block_size,ops);
	ops->MultiVecCreateByMat(&mv_ws[3],block_size,mat,ops);				
	ops->MultiVecSetRandomValue(mv_ws[3],0,block_size,ops);

	/* ���� nevInit ���趨Ҫ�� EigenSolverSetup_GCG �� nevInit һ�� */
	sizeV = nevInit+2*block_size;
	int length_dbl_ws = 2*sizeV*sizeV+10*sizeV
		+(nevMax+2*block_size)+(nevMax)*block_size;
	ops->Printf ( "length_dbl_ws = %d\n", length_dbl_ws );
	int length_int_ws = 6*sizeV+2*(block_size+3);
	ops->Printf ( "length_int_ws = %d\n", length_int_ws );
	if (dbl_ws!=NULL) {
		*dbl_ws = malloc(length_dbl_ws*sizeof(double));
		memset(*dbl_ws,0,length_dbl_ws*sizeof(double));
	} 
	if (int_ws!=NULL) {
		*int_ws = malloc(length_int_ws*sizeof(int));	
		memset(*int_ws,0,length_int_ws*sizeof(int));	
	}	
	return;
}
void EigenSolverDestroyWorkspace_GCG(
	int nevInit, int nevMax, int block_size, void *mat,
	void ***mv_ws, double **dbl_ws, int **int_ws, 
	struct OPS_ *ops)
{
	assert(mv_ws!=NULL);
	ops->MultiVecDestroy(&mv_ws[0],nevMax+2*block_size,ops);
	ops->MultiVecDestroy(&mv_ws[1],block_size,ops);
	ops->MultiVecDestroy(&mv_ws[2],block_size,ops);
	ops->MultiVecDestroy(&mv_ws[3],block_size,ops);
	if (dbl_ws!=NULL) {
		free(*dbl_ws); *dbl_ws = NULL;
	}
	if (int_ws!=NULL) {
		free(*int_ws); *int_ws = NULL;
	}
	return;
}


/* �����趨������Ҫ�� Setup ֮����� */
void EigenSolverSetParameters_GCG(
	int check_conv_max_num,
	const char *initX_orth_method, int initX_orth_block_size, int initX_orth_max_reorth, double initX_orth_zero_tol,
	const char *compP_orth_method, int compP_orth_block_size, int compP_orth_max_reorth, double compP_orth_zero_tol,
	const char *compW_orth_method, int compW_orth_block_size, int compW_orth_max_reorth, double compW_orth_zero_tol,
	int compW_cg_max_iter, double compW_cg_rate, double compW_cg_tol, const char *compW_cg_tol_type, int compW_cg_auto_shift,
	int compRR_min_num, double compRR_min_gap, double compRR_tol, 
	struct OPS_ *ops)
{
	
	struct GCGSolver_ *gcg_solver = (GCGSolver*)ops->eigen_solver_workspace;
	if (check_conv_max_num>0)
		gcg_solver->check_conv_max_num = check_conv_max_num;
	if (initX_orth_method!=NULL)
		strcpy(gcg_solver->initX_orth_method, initX_orth_method);
	if (initX_orth_block_size>0)
		gcg_solver->initX_orth_block_size = initX_orth_block_size;
	if (initX_orth_max_reorth>=0)
		gcg_solver->initX_orth_max_reorth = initX_orth_max_reorth;
	if (initX_orth_zero_tol>0)
		gcg_solver->initX_orth_zero_tol   = initX_orth_zero_tol;
	
	if (compP_orth_method!=NULL)
		strcpy(gcg_solver->compP_orth_method, compP_orth_method);
	if (compP_orth_block_size>0)
		gcg_solver->compP_orth_block_size = compP_orth_block_size;
	if (compP_orth_max_reorth>=0)
		gcg_solver->compP_orth_max_reorth = compP_orth_max_reorth;
	if (compP_orth_zero_tol>0)
		gcg_solver->compP_orth_zero_tol   = compP_orth_zero_tol;	

	if (compW_orth_method!=NULL)
		strcpy(gcg_solver->compW_orth_method, compW_orth_method);
	if (compW_orth_block_size>0)
		gcg_solver->compW_orth_block_size = compW_orth_block_size;
	if (compW_orth_max_reorth>=0)
		gcg_solver->compW_orth_max_reorth = compW_orth_max_reorth;
	if (compW_orth_zero_tol>0)
		gcg_solver->compW_orth_zero_tol   = compW_orth_zero_tol;
	if (compW_cg_max_iter>0)	
		gcg_solver->compW_cg_max_iter = compW_cg_max_iter;
	if (compW_cg_rate>0)
		gcg_solver->compW_cg_rate     = compW_cg_rate;
	if (compW_cg_tol>0)
		gcg_solver->compW_cg_tol      = compW_cg_tol;
	if (compW_cg_tol_type!=NULL)
		strcpy(gcg_solver->compW_cg_tol_type, compW_cg_tol_type);
	gcg_solver->compW_cg_auto_shift       = compW_cg_auto_shift;	
		
	if (compRR_min_gap>0)
		gcg_solver->compRR_min_gap = compRR_min_gap;
	if (compRR_min_num>0)
		gcg_solver->compRR_min_num = compRR_min_num;
	if (compRR_tol>0)
		gcg_solver->compRR_tol     = compRR_tol;
	
	return;	
}

void EigenSolverSetParametersFromCommandLine_GCG(
	int argc, char* argv[], struct OPS_ *ops)
{
	struct GCGSolver_ *gcg_solver = (GCGSolver*)ops->eigen_solver_workspace;

	ops->GetOptionFromCommandLine("-gcge_max_multi"  ,'i',
		&gcg_solver->multiMax  ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_min_gap"    ,'f',
		&gcg_solver->gapMin    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_given_nevec",'i',
		&gcg_solver->nevGiven  ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_max_niter"  ,'i',
		&gcg_solver->numIterMax,argc,argv, ops);	
	ops->GetOptionFromCommandLine("-gcge_abs_tol"    ,'f',
		&gcg_solver->tol[0]    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_rel_tol"    ,'f',
		&gcg_solver->tol[1]    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_user_defined_multi_lin_sol",'i',
		&gcg_solver->user_defined_multi_linear_solver,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_initX_orth_method"    ,'s',
		&gcg_solver->initX_orth_method    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_initX_orth_block_size",'i',
		&gcg_solver->initX_orth_block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_initX_orth_max_reorth",'i',
		&gcg_solver->initX_orth_max_reorth,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_initX_orth_zero_tol"  ,'f',
		&gcg_solver->initX_orth_zero_tol  ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_check_conv_max_num"   ,'i',
		&gcg_solver->check_conv_max_num,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compP_orth_method"    ,'s',
		&gcg_solver->compP_orth_method    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compP_orth_block_size",'i',
		&gcg_solver->compP_orth_block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compP_orth_max_reorth",'i',
		&gcg_solver->compP_orth_max_reorth,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compP_orth_zero_tol"  ,'f',
		&gcg_solver->compP_orth_zero_tol  ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compW_orth_method"    ,'s',
		&gcg_solver->compW_orth_method    ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_orth_block_size",'i',
		&gcg_solver->compW_orth_block_size,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_orth_max_reorth",'i',
		&gcg_solver->compW_orth_max_reorth,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_orth_zero_tol"  ,'f',
		&gcg_solver->compW_orth_zero_tol  ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compW_cg_max_iter"  ,'i',
		&gcg_solver->compW_cg_max_iter  ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_rate"      ,'f',
		&gcg_solver->compW_cg_rate      ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_tol"       ,'f',
		&gcg_solver->compW_cg_tol       ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_tol_type"  ,'s',
		&gcg_solver->compW_cg_tol_type  ,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_auto_shift",'i',
		&gcg_solver->compW_cg_auto_shift,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_shift"     ,'f',
		&gcg_solver->compW_cg_shift,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compW_cg_order"     ,'i',
		&gcg_solver->compW_cg_order     ,argc,argv, ops);
	
	ops->GetOptionFromCommandLine("-gcge_compRR_min_num",'i',
		&gcg_solver->compRR_min_num,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compRR_min_gap",'i',
		&gcg_solver->compRR_min_gap,argc,argv, ops);
	ops->GetOptionFromCommandLine("-gcge_compRR_tol    ",'f',
		&gcg_solver->compRR_tol    ,argc,argv, ops);

	int print_usage = 1;
	ops->GetOptionFromCommandLine("-gcge_print_usage",'i',&print_usage,argc,argv, ops);
	if (print_usage) {
       ops->Printf("\n");
       ops->Printf("Usage: %s [<options>]\n", argv[0]);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_max_multi    <i>: maximum of multiplicity of eigenpairs       %d (default: 6)\n",gcg_solver->multiMax);
       ops->Printf(" -gcge_min_gap      <f>: minimum of gap of eigenvalues relatively    %.2e (default: 1e-2)\n",gcg_solver->gapMin);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_max_niter    <i>: maximum of gcg iterations                   %d (default: 100)\n",gcg_solver->numIterMax);
       ops->Printf(" -gcge_given_nevec  <i>: number of given initial eigenvectors        %d (default: 0)\n",gcg_solver->nevGiven);
       ops->Printf(" -gcge_abs_tol      <f>: absolute convergence tolerance              %.2e (default: 1e-4)\n",gcg_solver->tol[0]);
       ops->Printf(" -gcge_rel_tol      <f>: relative convergence tolerance              %.2e (default: 1e-4)\n",gcg_solver->tol[1]);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_user_defined_multi_lin_sol  <i>: use user-defined multi linear solver  %d (default: 0[1])\n",gcg_solver->user_defined_multi_linear_solver);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_method  <s>: use which kind of orthogonalization for X  %s (default: mgs[bgs])\n",gcg_solver->initX_orth_method);
       ops->Printf(" -gcge_compP_orth_method  <s>: use which kind of orthogonalization for P  %s (default: bqr[bgs|mgs])\n",gcg_solver->compP_orth_method);
       ops->Printf(" -gcge_compW_orth_method  <s>: use which kind of orthogonalization for W  %s (default: mgs[bgs])\n",gcg_solver->compW_orth_method);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_block_size  <i>: size of vectors orthogonalized in one patch for X  %d (default: -1)\n",gcg_solver->initX_orth_block_size);
       ops->Printf(" -gcge_compP_orth_block_size  <i>: size of vectors orthogonalized in one patch for P  %d (default: -1)\n",gcg_solver->compP_orth_block_size);
       ops->Printf(" -gcge_compW_orth_block_size  <i>: size of vectors orthogonalized in one patch for W  %d (default: -1)\n",gcg_solver->compW_orth_block_size);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_zero_tol  <f>: zero tolerance in orthogonal for X  %.2e (default: 1e-16)\n",gcg_solver->initX_orth_zero_tol);
       ops->Printf(" -gcge_compP_orth_zero_tol  <f>: zero tolerance in orthogonal for P  %.2e (default: 1e-16)\n",gcg_solver->compP_orth_zero_tol);
       ops->Printf(" -gcge_compW_orth_zero_tol  <f>: zero tolerance in orthogonal for W  %.2e (default: 1e-16)\n",gcg_solver->compW_orth_zero_tol);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_initX_orth_max_reorth  <i>: maximum reorthogonal times for X  %d (default: 2)\n",gcg_solver->initX_orth_max_reorth);
       ops->Printf(" -gcge_compP_orth_max_reorth  <i>: maximum reorthogonal times for P  %d (default: 2)\n",gcg_solver->compP_orth_max_reorth);
       ops->Printf(" -gcge_compW_orth_max_reorth  <i>: maximum reorthogonal times for W  %d (default: 2)\n",gcg_solver->compW_orth_max_reorth);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_compW_cg_max_iter   <i>: maximum number of cg iteration       %d (default: 30)\n",gcg_solver->compW_cg_max_iter);
       ops->Printf(" -gcge_compW_cg_rate       <f>: descent rate of residual in cg       %.2e (default: 1e-2)\n",gcg_solver->compW_cg_rate);
       ops->Printf(" -gcge_compW_cg_tol        <f>: convergence tolerance in cg          %.2e (default: 1e-8)\n",gcg_solver->compW_cg_tol);
       ops->Printf(" -gcge_compW_cg_tol_type   <s>: type of convergence tolerance in cg  %s (default: abs[rel|user])\n",gcg_solver->compW_cg_tol_type);
       ops->Printf(" -gcge_compW_cg_order      <i>: order of krylov space for W in cg    %d (default: 1[2])\n",gcg_solver->compW_cg_order);
       ops->Printf(" -gcge_compW_cg_auto_shift <i>: shift automatically in cg            %d (default: 0[1])\n",gcg_solver->compW_cg_auto_shift);
       ops->Printf(" -gcge_compW_cg_shift      <f>: shift manually in cg                 %.2e (default: 0.0)\n",gcg_solver->compW_cg_shift);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_compRR_min_num  <i>: minimum number for splitting RR eval  %d (default: 10)\n",gcg_solver->compRR_min_num);
       ops->Printf(" -gcge_compRR_min_gap  <f>: minimum gap for splitting RR eval     %.2e (default: 1e-2)\n",gcg_solver->compRR_min_gap);
       ops->Printf(" -gcge_compRR_tol      <f>: convergence tolerance in RR           %.2e (default: 1e-16)\n",gcg_solver->compRR_tol);
       ops->Printf("---------------------------------------------------------------------------------------------------\n");
       ops->Printf(" -gcge_print_orth_zero  <i>: print the zero index in orthogonal      %d (default: 0[1])\n",1);
       ops->Printf(" -gcge_print_split      <i>: print the split information of RR eval  %d (default: 0[1])\n",0);
       ops->Printf(" -gcge_print_conv       <i>: print convergence in each iteration     %d (default: 1[0])\n",1);
       ops->Printf(" -gcge_print_eval       <i>: print the final eigenvalues             %d (default: 1[0])\n",1);
       ops->Printf(" -gcge_print_evec       <i>: print the final eigenvectors            %d (default: 0[1])\n",0);
       ops->Printf(" -gcge_print_time       <i>: print total time of each part           %d (default: 1[0])\n",1);
       ops->Printf(" -gcge_print_usage      <i>: print usage of gcg eigen solver         %d (default: 1[0])\n",1);
       ops->Printf("--------------------------------------------------------------------------------------------------\n");
       //ops->Printf(" -bpcg_print_res        <i>: print residual per five bpcg iteration  (default: 1[0])\n");
    }
	return;	
}
