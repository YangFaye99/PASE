#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pase_param.h"

void PASE_PARAMETER_Create(PASE_PARAMETER *param, int num_levels, int nev)
{

    *param = (PASE_PARAMETER)malloc(sizeof(PASE_PARAMETER_PRIVATE));

    (*param)->num_levels = num_levels;
    (*param)->initial_level = num_levels - 1;
    (*param)->mg_coarsest_level = num_levels - 1;
    (*param)->mg_finest_level = 0;
    (*param)->aux_coarse_level = num_levels - 1;
    (*param)->aux_fine_level = 0;

    (*param)->max_cycle_count = 10;
    (*param)->max_pre_count = 2;
    (*param)->max_post_count = 2;
    (*param)->max_direct_count = 20;
    (*param)->max_initial_direct_count = 20;

    (*param)->nev = nev;
    (*param)->more_nev = nev;
    (*param)->num_given_eigs = 0;

    (*param)->rtol = 1e-7;
    (*param)->atol = 1e-7;
    (*param)->initial_rtol = 1e-8;
    (*param)->aux_rtol = 1e-8;

    (*param)->multigrid_type = 0;
    (*param)->if_precondition = 1;
    (*param)->print_level = 2;
}

void PASE_PARAMETER_Destroy(PASE_PARAMETER *param)
{
    free(*param);
    *param = NULL;
}

void PASE_PARAMETER_Get_from_command_line(PASE_PARAMETER param, int argc, char *argv[])
{
    int arg_index = 0;

    while (arg_index < argc)
    {
        if (strcmp(argv[arg_index], "-nev") == 0)
        {
            arg_index++;
            param->nev = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-more_nev") == 0)
        {
            arg_index++;
            param->more_nev = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-num_levels") == 0)
        {
            arg_index++;
            param->num_levels = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-print_level") == 0)
        {
            arg_index++;
            param->print_level = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-initial_level") == 0)
        {
            arg_index++;
            param->initial_level = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-aux_coarse_level") == 0)
        {
            arg_index++;
            param->aux_coarse_level = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-aux_fine_level") == 0)
        {
            //要求解的特征值个数
            arg_index++;
            param->aux_fine_level = atoi(argv[arg_index++]);
        }
        else if (strcmp(argv[arg_index], "-if_precondition") == 0)
        {
            //要求解的特征值个数
            arg_index++;
            param->if_precondition = atoi(argv[arg_index++]);
        }
        else
        {
            arg_index++;
        }
    }
}
