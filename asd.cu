#include <stdio.h> 
#include <string.h> // for ffs
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_profiler_api.h>

/*Debug enabled gives more print statements of bids and how the "Matrix" gets evaluated*/
#define DEBUG 0
#define TRUE 1 
#define FALSE 0
/*Test sets all bids to one, which should give you n=|ITEMS| bids on output*/
#define TEST 1
/*Defines from 0-Range the random will give out*/
#define RANGE 10000
#define ITEMS 25

#define MAX (2 << (ITEMS-1))
#if ITEMS < 8
#define dint uint8_t  
#elif ITEMS < 16
#undef dint
#define dint uint16_t
#elif ITEMS < 32
#undef dint
#define dint uint32_t
#elif ITEMS > 32
#undef dint
#define dint uint64_t
#endif

#define SUBSET(X)((~_conf+(X+1))&_conf)
#define SETSUM(X)(f[setdiff(_conf,X)]+f[X])


static void HandleError( cudaError_t err, const char *file,int line ) {
	if (err != cudaSuccess) {
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ),file, line );
		exit( EXIT_FAILURE );
	}
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))



/*		           0 1 2 3 4 5 6 7 8			*/
dint * bids;
// dint bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
//dint bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
// dint bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6
dint * f;
dint * O;
 
struct _stack {
	dint conf;
	struct _stack * next;
} typedef stack;

struct _locklist {
	unsigned int size;
	unsigned int conf;
	unsigned int * dev_f;
	unsigned int * dev_o;
	struct _locklist * next;
} typedef locklist;

struct _lockstruct {
	unsigned int * dev_lock;
	struct _lockstruct * next;
} typedef lockstruct;


#define setdiff(seta,setb) (seta & ~setb)

inline  dint cardinality( dint seta) {
	return __builtin_popcount(seta);
}

void gen_rand_bids(dint MAXVAL) {
	register dint i = 0;
#if TEST
	for(i = 1; i < MAXVAL;i++) {
		bids[i] = 1;
		O[i] = i;
	}
//	bids[1] =0;
//	bids[2] = 0;
//	bids[32769] = 20;
#else
	for(i = 1; i < MAXVAL;i++) {
		bids[i] = rand() % RANGE;
		O[i] = i;
	}
	for(i = 1; i < MAXVAL;i*=2) {
		bids[i] = rand() % RANGE;

	}
#endif

#if DEBUG
	printf("i =\t");
	for(i = 1; i < MAX;i++) {
		printf("%u\t",i);
	}
	printf("\n");
	printf("val =\t");
	for(i = 1; i < MAX;i++) {
		printf("%u\t",bids[i]);
	}
	printf("\n");
#endif
}

const char *btb(dint y)
{
	int x = y;
	static char b[9];
	b[0] = '\0';
	
	int z;
	for (z = 128; z > 0; z >>= 1)
	{
		strcat(b, ((x & z) == z) ? "1" : "0");
	}	
	return b;
}

inline void printfo() {
#if DEBUG
	dint i;
	printf("i\t");
	for(i =1; i < MAX; i++) {
		printf("%u\t",i);
	}
	printf("\n");
	printf("f[]\t");
	for(i =1; i < MAX; i++) {
		printf("%u\t",f[i]);
	}
	printf("\n");
	printf("O[]\t");
	for(i =1; i < MAX; i++) {
		printf("%u\t",O[i]);
	}
	printf("\n");
#endif
}

/*Sets all bids with one element in it, |n| = 1*/
inline void set_singleton_bid(dint MAXVAL) {
	register  dint i;
	for(i =1; i< MAXVAL; i*=2) {
		f[i] = bids[i];
		if(bids[i] > 0)
			O[i] = i;
	}
}

/* conf a e.g. 1101
 * (~a+i) & a gives a subset of a
 *  i is a integer from 1 to |a|
 *
 * ~1101 = 0010
 * i = 0001
 * (0010+0001)&1101 =
 * 0011&1101 = 0001
 *
 *i = 0011
 * (0010+0011)&1101
 *(0101)&1101 = 0101
 */



void printfo(dint MAXVAL) {
	int i;
	printf("\n");
	for(i = 1; i< MAXVAL; i++)
	{
		printf("Bid[%d]\t%u\tF[%d]\t%u\tO[%d]\t%u\tbin\t%s\n",i,bids[i],i,f[i],i,O[i],btb(i));

	}
}

void parse_wopt(dint MAXVAL) {
	printf("parse maxval = %u\n",MAXVAL);
	//wopt at start contain MAX at wopt[0] which is the combination that goes in bids[wopt[n]]
	stack * root =(stack *) malloc(sizeof(stack));
	stack * sroot = NULL;
	stack * scurr = NULL;
	//DO N	OT REMOVE -1
	root->conf = (MAXVAL)-1;
	int count = 0;
	root->next = NULL;
	stack * curr = root;
	while(curr) {
		dint conf = curr->conf;
		if(conf == 0) {
			printf("EXIT FAILURE\n");
			return;
		}
		/*if something is wrong*/
		if(count > 40) {
			fprintf(stderr,"Something went wrong at line %d in %s\n",__LINE__,__FILE__);
			return;
		}
			printf("curr %u\t",curr->conf);
		if(conf != O[conf]) {

			dint diff = setdiff(conf,O[conf]);
			curr->conf = O[conf];
			stack * tmp = (stack *) malloc(sizeof(stack));
			printf("diff %u\t O[conf] %u f %u\n",diff,O[conf],f[conf]);
			tmp->conf = diff;
			tmp->next = curr;
			root = tmp;
			curr = root;
			count++;
			continue;
		}
		if(sroot == NULL) {
			sroot = curr;
			scurr = curr;
		} else {
			/*set next pointer to the next singelton*/
			scurr->next = curr;
		}

		/*set the current singleton to be the pointer*/
		scurr = curr;
		curr = curr->next;
		/*clear the pinter to avoid infinite loop if cross-referenced*/
		scurr->next = NULL;

	}
	curr = sroot;
	dint tmp = 0;
	while(curr != NULL) {
		printf("conf %u value %u\n",curr->conf,bids[curr->conf]);
		tmp++;
		stack * tmp = curr;
		curr = curr->next;
		free(tmp);
	}
	printf("n = %u\n",tmp);
}

#define I (threadIdx.x + blockDim.x * blockIdx.x)
#define SET_TEST_FETCH(STEP,S1,S2) {					\
		S1 = S2 = 0U;						\
		STEP = SUBSET(ispec);					\
		if(__popc(STEP) <= cardmax && ispec <= maxval) {	\
			S1 = f[setdiff(_conf,STEP)];			\
			S2 = f[STEP];					\
		}							\
		ispec += inc;						\
	}

/* ispec += blockDim.x; change back if not working */

#define COMP_SET(V1,S1,V2,S2) {			\
		if(V1>V2) {			\
			V2 = V1;		\
			S2 = S1;		\
		}				\
	}
  

#define MAXBLOCKSIZE 256U
#define NAGENTS 23
#define NSTREAMS 16 
#define NPERBLOCK 8
#define HALFBLOCK 4
 
__global__ void subsetcomp22(
		 	     unsigned int * __restrict__ f, /*Bid value*/
			     unsigned int * __restrict__ O, /*The move array*/
			     unsigned int * __restrict__ lock,
			     unsigned int _conf, /*The configuration*/
			     unsigned int cardmax, /*cardinality of max allowance*/
			     unsigned int maxval,
			     unsigned int count,
			     unsigned int offset,
			     unsigned int defbid)
{
/*these arrays are shared between all threads in the same block */
	__shared__ unsigned int share[MAXBLOCKSIZE];
	__shared__ unsigned int step[MAXBLOCKSIZE];     

	unsigned int inc = gridDim.x*blockDim.x; //corrected the increment
	unsigned int ispec = I + offset;
//	int val11;
	int i;
	unsigned int val1[HALFBLOCK];//the value for one of the subset sums
	unsigned int val2[HALFBLOCK];//the value for the other subset sums
	unsigned int stept[HALFBLOCK]; // the step array
	step[threadIdx.x] = share[threadIdx.x] = 0U;
	if(ispec <= maxval) {

/*Local for the thread, check all its bid and pick the biggest*/
#pragma unroll 4
		for(i = 0; i < HALFBLOCK; i++) {				
			SET_TEST_FETCH(stept[i],val1[i],val2[i]);
		}
#pragma unroll 4
		for(i = 0; i < HALFBLOCK; i++) {		
			val1[i] += val2[i];
			COMP_SET(val1[i],stept[i],share[threadIdx.x],step[threadIdx.x]);			
			SET_TEST_FETCH(stept[i],val1[i],val2[i]); // pipelined fetch
		}

#pragma unroll 4
		for(i = 0; i < HALFBLOCK; i++) {		
			val1[i] += val2[i];
			COMP_SET(val1[i],stept[i],share[threadIdx.x],step[threadIdx.x]);
			//	SET_TEST_FETCH(stept[i],val1[i],val2[i]);
		}
	}

	ispec = I;
       
	i= blockDim.x >> 1;
	__syncthreads();
/*do max reduction on the shared array for all threads inside the block*/
#pragma unroll
	for (; i>0; i>>=1) {
		if (threadIdx.x < i && (ispec <= maxval)) {
			if(share[threadIdx.x] < share[threadIdx.x + i]) {
				step[threadIdx.x] = step[threadIdx.x+i];
				share[threadIdx.x] = share[threadIdx.x+i];
			}
		}
		__syncthreads();
	}

/*thread 0 will attempt to set to global memory the agreed maximum value inside the block,
* if it is greater than the original bid and the bid in the lock array
*/
	if(threadIdx.x == 0U) {
		i = share[0U];
		if(defbid => i)
			return;
		if(lock[count] < i) {
			if(atomicMax(&(lock[count]),i) < i) {
				O[_conf] = step[0U];
				f[_conf] = i; 
				
			} 
		}
	}
}



#define COMBS(X) ((1 << cardinality(X)) - 1)

int run_test(dint MAXVAL,dint items) {
/*Setup the environment*/
	//dint perm[MAXVAL];
	printfo();
	register unsigned int i, c,count =0;
 	unsigned int *dev_f,*dev_o;

	i = items/2;
	count = 0;

	HANDLE_ERROR(cudaDeviceReset());

	unsigned int * dev_lock1,*dev_lock2,*dev_ptr;
	const	unsigned int devcount = 1024;// count;
	register unsigned int streams = NSTREAMS;
	register unsigned int count2 = 0;
	register unsigned int streamcount = 0;
	register cudaStream_t stream[streams];
	for(int i = 0;i < streams; i++)
		HANDLE_ERROR(cudaStreamCreate(&stream[i]));
	printf("count %u\n",devcount);
	count = 0;
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock1,(10+devcount)*sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock2,(10+devcount)*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_f, MAXVAL*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_o, MAXVAL*sizeof(int)));

 	HANDLE_ERROR(cudaMemcpy(dev_f,bids,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));
 	HANDLE_ERROR(cudaMemcpy(dev_o,O,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));

	HANDLE_ERROR(cudaMemset(dev_lock1,0,devcount*sizeof(int)));
	HANDLE_ERROR(cudaMemset(dev_lock2,0,devcount*sizeof(int)));
	/*2.*/
//	printfo(MAXVAL); printf("before\n");
	dev_ptr = dev_lock1;
	register unsigned int bsize = 0;
	register int blocks;
	int prev =0;
	count2 = 0;
	for(i = 2; i <= items; i++) {
		time_t start,end,t;
		start=clock();
		for(c = (1 << i) -1; c <= MAXVAL;) {

			double tmp = (double) COMBS(c)/NPERBLOCK;
			
			while( bsize < MAXBLOCKSIZE && tmp > bsize) {
				bsize += 32;
			}
			blocks =(int)  ceil((tmp/bsize));
/* #if __CUDA_ARCH__ < 300 */
/* 			int remaindern = blocks - 65535; */
/* 			while( blocks > 65535 ) { */
/* 				bsize += 32; */
/* 				blocks =(int)  ceil((tmp/bsize)); */
/* 			} */
/* 			printf("hello"); */
/* 			//double bsize = BLOCKSIZE; */
/* 			if(remaindern > 0) { */
/* 				blocks =65535; */
/* 				subsetcomp22<<<remaindern,bsize,0,stream[streamcount]>>>(dev_f,dev_o,dev_ptr,c,i/2,tmp,count2,65535*bsize,bids[c]); */
/* 			} */
/* #endif */
			subsetcomp22<<<blocks,bsize,0,stream[streamcount]>>>(dev_f,dev_o,dev_ptr,c,i/2,tmp,count2,0,bids[c]);
			//	printf("blocks %d block size %u stream count %d\n",blocks,bsize,streamcount);
			t = c | (c-1);
			c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c) + 1));
			count++;
			count2++;	
			streamcount++;
			if(streamcount >= streams)
				streamcount = 0;
			if(count2 < devcount)
				continue;
			HANDLE_ERROR(cudaMemset(dev_ptr,0,devcount*sizeof(int)));

			if(dev_ptr == dev_lock1)
				dev_ptr = dev_lock2;
			else
				dev_ptr = dev_lock1;
			count2 = 0;
		}
		end=clock();
		t=(end-start)/CLOCKS_PER_SEC;
		printf("ended card %d blocks\t %d threads/block %u, n kernels %u \t time %lu\n",i,blocks,bsize,count-prev,t);
		prev =	count;
		for (int i = 0; i < streams; ++i)
			HANDLE_ERROR(cudaStreamSynchronize(stream[i]));

		HANDLE_ERROR(cudaDeviceSynchronize());
		printfo();
	}
	for (int i = 0; i < streams; ++i)
	cudaStreamDestroy(stream[i]);

	HANDLE_ERROR(cudaDeviceSynchronize());
	HANDLE_ERROR(cudaMemcpy(f,dev_f,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(O,dev_o,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
	//int i;
	HANDLE_ERROR(cudaFree(dev_f));
	HANDLE_ERROR(cudaFree(dev_o));
	HANDLE_ERROR(cudaFree(dev_lock1));
	HANDLE_ERROR(cudaFree(dev_lock2));

	HANDLE_ERROR(cudaDeviceReset());
//	printfo(MAXVAL);
	return count;
}



int main(void) {
	/*Start n amount of assets*/
	dint from = NAGENTS;
	dint MAXVAL = (2 << (from-1));

	time_t start,end,t;
	O = (dint * ) malloc(sizeof(dint)*(2 << (from-1)));
	bids = (dint * ) malloc(sizeof(dint)*(2 << (from-1)));
	f = bids;


	MAXVAL = (2 << (from-1));
	gen_rand_bids(MAXVAL);
	set_singleton_bid(MAXVAL);
	printf("maxval %u from %u\n",MAXVAL,from);
	start=clock();//predefined  function in c
	int count = run_test(MAXVAL,from);
	end=clock();
	t=(end-start)/CLOCKS_PER_SEC;
	parse_wopt(MAXVAL);
	printf("\nTime taken =%lu for n= %u with count %d average per count %lf\n", (unsigned long) t,from,count,(double)t/count);


	free(O);
	free(f);

	return 0;
}
