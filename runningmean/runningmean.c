#include <stdio.h>
#include <omp.h>
#define MAXNCHANS  10000

void run(unsigned long long int *data, int nsamp, int nchan, unsigned char *outdata, int windowsize){

    float ave[MAXNCHANS];
    unsigned long long int head[MAXNCHANS];
    unsigned long long int tail[MAXNCHANS];
    float sum[MAXNCHANS];
    int j,k;
    float ws = windowsize;
    /*printf("windowsize: %f", ws);*/
///*#pragma omp parallel for default(shared) private(j) shared(sum)*/
    for(j=nchan-1;j>=0;j--){
    /*for(j=nchan/2;j>=nchan/2-1;j--){*/
        sum[j] = 0.;
        for (k=0;k<windowsize;k++){
            sum[j] += data[k*nchan+j];
        }
        ave[j] = sum[j]/ws;
        /*printf("ave[%d] = %f; sum[%d] = %f\n", j, ave[j], j, sum[j]);*/
    }


    for(int i=windowsize; i<nsamp; i++){
#pragma omp parallel for default(shared) private(j) shared(head, tail, ave)
        for(j=nchan-1;j>=0;j--){
        /*for(j=nchan/2;j>=nchan/2-1;j--){*/
            head[j] = data[(i-windowsize)*nchan+j];
            tail[j] = data[i*nchan+j];
            outdata[(i-windowsize)*nchan + nchan-j-1] = (unsigned char) (head[j] > ave[j]);
            ave[j]  += ((float)tail[j] - (float)head[j])/ws;
            /*printf("i %d, j %d, head %d tail %d (tail[j] - head[j])/ws %f, ave[j] %f\n ", i, j, head[j], tail[j],  ((float)tail[j] - (float)head[j])/ws, ave[j]);*/
            /*printf("outdata[] %d;  (head[j] - ave[j]) %f \n",outdata[(i-windowsize)*nchan + nchan-j-1], ((float)head[j] - (float)ave[j]));*/
        }
    }
    for(int i=0;i<windowsize;i++){
#pragma omp parallel for default(shared) private(j) shared(head, tail, ave)
        for(j=nchan-1;j>=0;j--){
        /*for(j=nchan/2;j>=nchan/2-1;j--){*/
            head[j] =  data[(nsamp+i-windowsize)*nchan + j];
            tail[j] =  data[(nchan-1*i)*nchan + j];
            outdata[(i+nsamp-windowsize)*nchan + nchan-j-1] = (unsigned char) (head[j] > ave[j]);
            ave[j]  += ((float)tail[j] - (float)head[j])/ws;
            /*printf("i %d, j %d, head %d tail %d (tail[j] - head[j])/ws %f, ave[j] %f\n ", i, j, head[j], tail[j],  ((float)tail[j] - (float)head[j])/ws, ave[j]);*/
            /*printf("outdata[] %d;  (head[j] - ave[j]) %f \n",outdata[(i+nsamp-windowsize)*nchan + nchan-j-1], ((float)head[j] - (float)ave[j]));*/
        }
    }
}

