#include <stdio.h>
#include <omp.h>
#define MAXNCHANS  10000

void run(unsigned char *data, int nsamp, int nchan, unsigned char *outdata, int windowsize){
    /* Loop through each channel i, using ave[i] to track the running mean and
     * head[i] and tail[i] to track the head and tail. Using the head and tail
     * value of the window to update the runnig mean. sum[i] is used only to
     * calculate the ave[i]
     */

    float ave[MAXNCHANS];
    float head[MAXNCHANS];
    float tail[MAXNCHANS];
    float sum[MAXNCHANS];
    int j,k;
    float ws = windowsize;
    /*printf("windowsize: %f", ws);*/
#pragma omp parallel for default(shared) private(j) shared(sum)
    for(j=nchan-1;j>=0;j--){
    /*for(j=nchan/2;j>=nchan/2-1;j--){*/
        sum[j] = 0.;
        for (k=0;k<windowsize;k++){
            sum[j] += (float)data[2*k*nchan+j] + (float)data[(2*k+1)*nchan+j];
        }
        ave[j] = sum[j]/ws;
        /*printf("ave[%d] = %f; sum[%d] = %f\n", j, ave[j], j, sum[j]);*/
    }


    for(int i=windowsize; i<nsamp; i++){
#pragma omp parallel for default(shared) private(j) shared(head, tail, ave)
        for(j=nchan-1;j>=0;j--){
        /*for(j=nchan/2;j>=nchan/2-1;j--){*/
            k = (i-windowsize);
            head[j] = (float) data[2*k*nchan+j] + (float) data[(2*k+1)*nchan+j];
            tail[j] = (float) data[2*i*nchan+j] + (float) data[(2*i+1)*nchan+j] ;
            outdata[k*nchan + nchan-j-1] = (unsigned char) (head[j] > ave[j]);
            ave[j]  += (tail[j] - head[j])/ws;
            /*printf("i %d, j %d, head %d tail %d (tail[j] - head[j])/ws %f, ave[j] %f\n ", i, j, head[j], tail[j],  ((float)tail[j] - (float)head[j])/ws, ave[j]);*/
            /*printf("outdata[] %d;  (head[j] - ave[j]) %f \n",outdata[(i-windowsize)*nchan + nchan-j-1], ((float)head[j] - (float)ave[j]));*/
        }
    }
    for(int i=0;i<windowsize;i++){
#pragma omp parallel for default(shared) private(j) shared(head, tail, ave)
        for(j=nchan-1;j>=0;j--){
        /*for(j=nchan/2;j>=nchan/2-1;j--){*/
            k = (nsamp+i-windowsize);
            head[j] =  (float) data[2*k*nchan + j] + (float) data[(2*k+1)*nchan + j] ;
            k = (nsamp-i-1);
            tail[j] =  (float) data[2*k*nchan + j] + (float) data[(2*k+1)*nchan + j];
            outdata[(nsamp+i-windowsize)*nchan + nchan-j-1] = (unsigned char) (head[j] > ave[j]);
            ave[j]  += (tail[j] - head[j])/ws;
            /*printf("i %d, j %d, head %d tail %d (tail[j] - head[j])/ws %f, ave[j] %f\n ", i, j, head[j], tail[j],  ((float)tail[j] - (float)head[j])/ws, ave[j]);*/
            /*printf("outdata[] %d;  (head[j] - ave[j]) %f \n",outdata[(i+nsamp-windowsize)*nchan + nchan-j-1], ((float)head[j] - (float)ave[j]));*/
        }
    }
}

