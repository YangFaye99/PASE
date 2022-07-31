#BSUB -J job_name
#BSUB -n 108
#BSUB -o %J.submit.out
#BSUB -e %J.submit.err
#BSUB -W 180
#BSUB -q batch
#BSUB -R span[ptile=36]
#BSUB -x

cd $LS_SUBCWD
mpirun ./test0.exe -mat_superlu_dist_rowperm NOROWPERM -mat_superlu_dist_colperm PARMETIS \
                   -mat_superlu_dist_replacetinypivot