#include <stdio.h> 
#include "asd.h"
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


unsigned int pascal[30][30] =  
{{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,3,3,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,4,6,4,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,5,10,10,5,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,6,15,20,15,6,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,7,21,35,35,21,7,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,8,28,56,70,56,28,8,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,9,36,84,126,126,84,36,9,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,10,45,120,210,252,210,120,45,10,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,11,55,165,330,462,462,330,165,55,11,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,12,66,220,495,792,924,792,495,220,66,12,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,13,78,286,715,1287,1716,1716,1287,715,286,78,13,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,14,91,364,1001,2002,3003,3432,3003,2002,1001,364,91,14,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,15,105,455,1365,3003,5005,6435,6435,5005,3003,1365,455,105,15,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,16,120,560,1820,4368,8008,11440,12870,11440,8008,4368,1820,560,120,16,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,17,136,680,2380,6188,12376,19448,24310,24310,19448,12376,6188,2380,680,136,17,1,0,0,0,0,0,0,0,0,0,0,0,0},
 {1,18,153,816,3060,8568,18564,31824,43758,48620,43758,31824,18564,8568,3060,816,153,18,1,0,0,0,0,0,0,0,0,0,0,0},
 {1,19,171,969,3876,11628,27132,50388,75582,92378,92378,75582,50388,27132,11628,3876,969,171,19,1,0,0,0,0,0,0,0,0,0,0},
 {1,20,190,1140,4845,15504,38760,77520,125970,167960,184756,167960,125970,77520,38760,15504,4845,1140,190,20,1,0,0,0,0,0,0,0,0,0},
 {1,21,210,1330,5985,20349,54264,116280,203490,293930,352716,352716,293930,203490,116280,54264,20349,5985,1330,210,21,1,0,0,0,0,0,0,0,0},
 {1,22,231,1540,7315,26334,74613,170544,319770,497420,646646,705432,646646,497420,319770,170544,74613,26334,7315,1540,231,22,1,0,0,0,0,0,0,0},
 {1,23,253,1771,8855,33649,100947,245157,490314,817190,1144066,1352078,1352078,1144066,817190,490314,245157,100947,33649,8855,1771,253,23,1,0,0,0,0,0,0},
 {1,24,276,2024,10626,42504,134596,346104,735471,1307504,1961256,2496144,2704156,2496144,1961256,1307504,735471,346104,134596,42504,10626,2024,276,24,1,0,0,0,0,0},
 {1,25,300,2300,12650,53130,177100,480700,1081575,2042975,3268760,4457400,5200300,5200300,4457400,3268760,2042975,1081575,480700,177100,53130,12650,2300,300,25,1,0,0,0,0},
 {1,26,325,2600,14950,65780,230230,657800,1562275,3124550,5311735,7726160,9657700,10400600,9657700,7726160,5311735,3124550,1562275,657800,230230,65780,14950,2600,325,26,1,0,0,0},
 {1,27,351,2925,17550,80730,296010,888030,2220075,4686825,8436285,13037895,17383860,20058300,20058300,17383860,13037895,8436285,4686825,2220075,888030,296010,80730,17550,2925,351,27,1,0,0},
 {1,28,378,3276,20475,98280,376740,1184040,3108105,6906900,13123110,21474180,30421755,37442160,40116600,37442160,30421755,21474180,13123110,6906900,3108105,1184040,376740,98280,20475,3276,378,28,1,0},
 {1,29,406,3654,23751,118755,475020,1560780,4292145,10015005,20030010,34597290,51895935,67863915,77558760,77558760,67863915,51895935,34597290,20030010,10015005,4292145,1560780,475020,118755,23751,3654,406,29,1}};


/*		           0 1 2 3 4 5 6 7 8			*/
uint32_t * bids,* f, * O;
// dint bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
//dint bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
// dint bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6
 
struct _stack {
	dint conf;
	struct _stack * next;
} typedef stack;



#define setdiff(seta,setb) (seta & ~setb)

inline  dint cardinality( dint seta) {
	return __builtin_popcount(seta);
}
int indexa =0;
void gen_rand_bids(dint MAXVAL) {
	register dint i = 0;
#if TEST
	for(i = 1; i < MAXVAL;i++) {
		bids[i] = 1;
		O[i] = i;
	}
//     unsigned int seed = (unsigned)time ( NULL );
	//   srand(seed);
	indexa++;
	//indexa = rand() % MAXVAL;
     bids[indexa] = 100;
     printf("index %d \n",indexa);
     if(indexa >= MAXVAL) {
	     printf("No error\n");
	     exit(0);
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

/*Sets all bids with one element in it, |n| = 1*/
inline void set_singleton_bid(dint MAXVAL) {
	register  dint i;
	for(i =1; i< MAXVAL; i*=2) {
		f[i] = bids[i];
		if(bids[i] > 0)
			O[i] = i;
	}
}



void printfo(dint MAXVAL) {
	int i;
	printf("\n");
	for(i = 1; i< MAXVAL; i++)
	{
		printf("Bid[%d]\t%u\tF[%d]\t%u\tO[%d]\t%u\tbin\t%s\n",i,bids[i],i,f[i],i,O[i],btb(i));

	}
}

int parse_wopt(dint MAXVAL) {
	//printf("parse maxval = %u\n",MAXVAL);
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
			return 1;
		}
		/*if something is wrong*/
			if(count > 40) {
				fprintf(stderr,"Something went wrong at line %d in %s\n",__LINE__,__FILE__);
				return 1;
			}
		printf("curr %u\t\n",curr->conf);
		if(conf != O[conf]) {

			dint diff = setdiff(conf,O[conf]);
			curr->conf = O[conf];
			stack * tmp = (stack *) malloc(sizeof(stack));
			printf("diff %u\t conf %u\t O[diff] %u\t O[conf] %u\t f %u\n",diff,conf,O[diff],O[conf],f[conf]);
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
	printf("\n");
	curr = sroot;
	dint tmp = 0;
	while(curr != NULL) {
		if(bids[curr->conf]) {
			if(curr->conf == indexa) {
				printf("correct bid\n");
				tmp++;
				return 0;
			}
			printf("conf %u value %u\n",curr->conf,bids[curr->conf]);
			
		}
		stack * tmp = curr;
		curr = curr->next;
		free(tmp);
	}
	if(tmp < 1) {
		printf("something is wrong, no bids\n");
			       
		return 1;

	}
	printf("n = %u\n",tmp);
	return 0;
}

#define I (threadIdx.x + blockDim.x * blockIdx.x)
#define SET_TEST_FETCH(STEP,S1,S2) {					\
									\
		STEP = SUBSET(ispec);					\
		if( ispec < maxval) {					\
			S1 = f[(setdiff(_conf,STEP))];			\
			S2 = f[(STEP)];					\
		} else {						\
			S1 = S2 = 0U;					\
		}							\
		ispec++;						\
	}

/* ispec += blockDim.x; change back if not working */

#define COMP_SET(V1,S1,V2,S2) {			\
		if(V1>V2) {			\
			V2 = V1;		\
			S2 = S1;		\
		}				\
	}
  
    
#define MAXBLOCKSIZE 256U
#define NAGENTS 21  
#define NSTREAMS 8 
#define NPERBLOCK 8
#define HALFBLOCK 8

#define PASCALSIZE 30

static __constant__ unsigned int pascl[PASCALSIZE][PASCALSIZE];

static __constant__ unsigned int cardindex[NAGENTS+1];

__device__ unsigned int get_index(unsigned int set) {
	const unsigned int card = __popc(set);
	const unsigned int cardi = cardindex[card];
	unsigned int sum = 0;
	unsigned int tmp = set;
	int i;
#pragma unroll
	for(i = 1;i<= card; i++) {
		unsigned int fsb = __ffs(tmp) -1;
		sum += pascl[fsb][i];
		tmp &= ~(1 << fsb);
	}
	sum += cardi;
	return sum;
}

__global__ void subsetcomp22(
	/*0*/	uint32_t * __restrict__ f, /*Bid value*/ 
	/*1*/	unsigned int * __restrict__ O, /*The move array*/
	/*2*/	unsigned int * __restrict__ lock,
	/*3*/	unsigned int _conf, /*The configuration*/
	/*5*/	unsigned int maxval,
	/*6*/	unsigned int count,
	/*8*/	unsigned int defbid)
{
/*these arrays are shared between all threads in the same block */
	__shared__ unsigned int share[MAXBLOCKSIZE];
	__shared__ unsigned int step[MAXBLOCKSIZE];     
	unsigned int ispec = NPERBLOCK*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;
	const unsigned int tid = threadIdx.x;
	int i;  
	unsigned int max = 0;
	unsigned int rstep = 0;
	unsigned int val1[NPERBLOCK];//the value for one of the subset sums
	unsigned int val2[NPERBLOCK];//the value for the other subset sums
	unsigned int stept[NPERBLOCK]; // the step array
	if(ispec < maxval) {

/*Local for the thread, check all its bid and pick the biggest*/
#pragma unroll 8
		for(i = 0; i < NPERBLOCK; i++) {
			SET_TEST_FETCH(stept[i],val1[i],val2[i]);	
		}
#pragma unroll 8
		for(i = 0; i < NPERBLOCK; i++) {		
			val1[i] += val2[i];
			COMP_SET(val1[i],stept[i],max,rstep);			
		}

	}
	step[threadIdx.x] = rstep;
	share[threadIdx.x] = max;
	i= blockDim.x >> 1;
	__syncthreads();
/*do max reduction on the shared array for all threads inside the block*/
#pragma unroll
	for (; i>0; i>>=1) {
		if (tid < i /* && (ispec <= maxval) */) {
			if(share[tid] <= share[tid + i]) {
				step[tid] = step[tid+i];
				share[tid] = share[tid+i];
			}
		}
		__syncthreads();
	}

/*thread 0 will attempt to set to global memory the agreed maximum value inside the block,
* if it is greater than the original bid and the bid in the lock array
*/

	if(tid == 0U) {
		i = share[0U];
		if(i == 0)
			return;
		if(defbid >= i)
			return;
		if(lock[count] < i) {
			if(atomicMax(&(lock[count]),i) < i) {
				O[_conf] = step[0U];
				f[_conf] = i; 
				__threadfence();
				return;
			} 
		}
	} else {
		return;
	}
}

int gen_copy_base_index() {
	unsigned int *tmpa =(unsigned int *) malloc(sizeof(unsigned int)*(NAGENTS+1));
	if(tmpa == NULL) {
		fprintf(stderr,"Can not allocate memory at line %s in %s\n",__LINE__,__FILE__);
		exit(1);
	}
	int i;
	int tmp = 0;
	for(i = 1; i <= NAGENTS;i++) {
		tmp += pascal[NAGENTS][i-1];
		tmpa[i] = tmp;
		//printf("tmp = %d for i %d\n",tmp,i);
	       
	}
	HANDLE_ERROR(cudaMemcpyToSymbol(cardindex, tmpa, sizeof(unsigned int)*(NAGENTS+1), 0, cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpyToSymbol(pascl, pascal, sizeof(unsigned int)*(PASCALSIZE*PASCALSIZE), 0, cudaMemcpyHostToDevice));
	return 0;
}

#define COMBS(X) ((1 << cardinality(X)-1) - 1)

int run_test(dint MAXVAL,dint items) {
/*Setup the environment*/
	//dint perm[MAXVAL];

	register unsigned int i, c,count =0;
 	unsigned int *dev_f,*dev_o;

	i = items/2;
	count = 0;

	HANDLE_ERROR(cudaDeviceReset());
//  	cudaDeviceSetCacheConfig( cudaFuncCachePreferL1 );

	unsigned int * dev_lock1,*dev_lock2,*dev_ptr;
	const	unsigned int devcount = 1024;// count;
	register unsigned int streams = NSTREAMS;
	register unsigned int lock_count = 0;
	register unsigned int streamcount = 0;
	register cudaStream_t stream[streams];
	for(int i = 0;i < streams; i++)
		HANDLE_ERROR(cudaStreamCreate(&stream[i]));
//	printf("count %u\n",devcount);
	count = 0;
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock1,(10+devcount)*sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock2,(10+devcount)*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_f, MAXVAL*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_o, MAXVAL*sizeof(int)));

 	HANDLE_ERROR(cudaMemcpy(dev_f,bids,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));
 	HANDLE_ERROR(cudaMemcpy(dev_o,O,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));

	HANDLE_ERROR(cudaMemset(dev_lock1,0,devcount*sizeof(int)));
	HANDLE_ERROR(cudaMemset(dev_lock2,0,devcount*sizeof(int)));
	gen_copy_base_index();
	/*2.*/
	//printfo(MAXVAL); printf("before\n");
	dev_ptr = dev_lock1;
	register unsigned int bsize = 0;
	register int blocks;
	int prev =0;
	lock_count = 0;
	time_t rstart,rend,rt;
	rstart=clock();
	for(i = 2; i <= items; i++) {
		time_t start,end,t;
		
		start=clock();
		int splittings;
		double threads;
		c = c = (1 << i) -1;
		splittings =  COMBS(c);///NPERBLOCK;
		threads = ((double) splittings)/ NPERBLOCK;
		threads = ceil(threads);
		for(; c <= MAXVAL;) {
			while( bsize < MAXBLOCKSIZE && threads > bsize) {
				bsize += 32;
			}
			blocks =(int)  ceil((threads/bsize));

			subsetcomp22<<<blocks,bsize,0,stream[streamcount]>>>(dev_f,dev_o,dev_ptr,c,splittings,lock_count,bids[c]);

			t = c | (c-1);
			c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c) + 1));

			count++;
			lock_count++;	
			streamcount++;
			if(streamcount >= streams)
				streamcount = 0;
			if(lock_count < devcount)
				continue;
			HANDLE_ERROR(cudaMemset(dev_ptr,0,devcount*sizeof(int)));

			if(dev_ptr == dev_lock1)
				dev_ptr = dev_lock2;
			else
				dev_ptr = dev_lock1;
			lock_count = 0;
		}

		for (int t = 0; t < streams; ++t) {
			HANDLE_ERROR(cudaStreamSynchronize(stream[t]));
		}

		HANDLE_ERROR(cudaDeviceSynchronize());
		

		end=clock();
		t=(end-start)/(CLOCKS_PER_SEC/1000);
		printf("ended card %d blocks\t %d threads/block %u, n kernels %u \t time %lu \t splittings %d\n",i,blocks,bsize,count-prev,t,splittings);
		prev =	count;

	}
	for (int i = 0; i < streams; ++i)
	cudaStreamDestroy(stream[i]);

	HANDLE_ERROR(cudaDeviceSynchronize());

	rend=clock();
	rt=(rend-rstart)/(CLOCKS_PER_SEC);
	printf("real time %lu\n",rt);

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
	int ret_val =0;
	int count;
	while(ret_val == 0) {
	
	MAXVAL = (2 << (from-1));
	gen_rand_bids(MAXVAL);
	set_singleton_bid(MAXVAL);
	printf("maxval %u from %u\n",MAXVAL,from);
	start=clock();//predefined  function in c
	 count = run_test(MAXVAL,from);
	end=clock();
	t=(end-start)/CLOCKS_PER_SEC;
	ret_val= parse_wopt(MAXVAL);
	printf("\nTime taken =%lu for n= %u with count %d average per count %lf\n", (unsigned long) t,from,count,(double)t/count);
	}


	free(O);
	free(f);

	return 0;
}
