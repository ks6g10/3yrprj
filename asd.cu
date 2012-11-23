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

#define BLOCKSIZE 512
#define HALFBLOCK 256
#if HALFBLOCK*2 != BLOCKSIZE
#error HALFBLOCK is not set correctly
#endif
#define SUBSET(X)((~_conf+(X+1))&_conf)
#define SETSUM(X)(f[setdiff(_conf,X)]+f[X])
#define I (threadIdx.x + BLOCKSIZE * blockIdx.x)

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
	bids[1] =0;
	bids[2] = 0;
#else
	for(i = 1; i < MAXVAL;i++) {
		bids[i] = rand() % RANGE;
	}
	for(i = 1; i < MAXVAL;i*=2) {
		bids[i] = 1;//rand() % RANGE;
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

/*n 15 t 9 n 16 t 42*/
void max2(dint conf) {
	register dint card = cardinality(conf)/2;
	register dint combinations = 1 << (cardinality(conf)-1);
	register dint max = bids[conf];
	register dint set = conf;
	register dint tmp = 0;
	register dint subset;
	register dint inverse = ~conf;
	register dint i;
	for(i = 1;i<combinations; i++) {
		subset = (inverse+i)&conf;
		if(cardinality(subset) > card)
			continue;
		tmp = f[setdiff(conf,subset)] + f[subset];
		if(max < tmp) {
			max = tmp;
			set = subset;
		}
	}
	f[conf] = max;
	O[conf] = set;
}


//the index




/*
 *
 * 1. gen all combinations of card n
 * 2. for each combination, generate all subset with condition |s| < |c|/2
 * 3. for each subset check if |s| < |c|/2 then compute the sum
 * 4. 
 *
 *
 *
 */
__global__ void setlock(unsigned int * lock) {
		lock[I] = 0;
}
__global__ void subsetcomp22(unsigned int * f, /*Bid value*/
			     unsigned int * O, /*The move array*/
			     unsigned int * lock,
			     unsigned int _conf, /*The configuration*/
			     unsigned int cardmax, /*cardinality of max allowance*/
			     unsigned int maxval,
			     unsigned int count,
			    	unsigned int offset)
{

	__shared__ unsigned int share[BLOCKSIZE];
	__shared__ unsigned int step[BLOCKSIZE];
		/* printf("threadid.x\t%d\tblockDim.x\t%d\tblockIdx.x\t%d\tI\t%u\n", */
	/*        threadIdx.x, */
	/*        blockDim.x, */
	/*        blockIdx.x, */
	/*        conf); */
	unsigned int i = I+offset;
	/*subset var, but also for indexing later on*/
	unsigned int s = SUBSET(i);
//	unsigned int s2 = SUBSET(i);

	unsigned int tid = threadIdx.x;


	step[tid] = share[tid] = 0;
	if(i < maxval) {
		if(__popc(s) <= cardmax ) {
			share[tid] = f[setdiff(_conf,s)] + f[s];
			step[tid] = s;
		}
	}
	s= blockDim.x >> 1;
	__syncthreads();
	for (; s>0; s>>=1) {
		if (tid < s && (i < maxval)) {
			if(share[tid] < share[tid + s]) {
				step[tid] = step[tid+s];
				share[tid] = share[tid+s];
			}
		}
		__syncthreads();
	}

	if(tid == 0) {
		if(lock[count+1] < share[0]) {
			if(atomicMax(&(lock[count]),share[0]) < share[0]) {
				lock[count+1] = share[0];
				O[_conf] = step[0];
				f[_conf] =share[0];
				
			} 
		}
	}
}



#define COMBS(X) ((1 << cardinality(X)) - 1)

int run_test(dint MAXVAL,dint items) {
/*Setup the environment*/
	//dint perm[MAXVAL];
	printfo();
	unsigned int i, c,t,count =0;
//	f = bids;
 	unsigned int *dev_f,*dev_o;

	i = items/2;
	count = 0;
	/* for(c = (1 << i) -1; c <= MAXVAL;) { */
	/* 	count++; */
	/* 	t = c | (c-1); */
	/* 	c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c) + 1)); */
	/* } */
	
	unsigned int * dev_lock;
	unsigned int devcount = 1024;// count;
	unsigned int count2 = 0;
	unsigned int * cpy_lock =(unsigned int *)malloc((devcount+10)*sizeof(int));
	for(i = 0; i< devcount+10; i++)
		cpy_lock[i] = 0;
	printf("count %u\n",devcount);
	

	count = 0;
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock,(10+devcount)*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_f, MAXVAL*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_o, MAXVAL*sizeof(int)));

 	HANDLE_ERROR(cudaMemcpy(dev_f,bids,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));
 	HANDLE_ERROR(cudaMemcpy(dev_o,O,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));
//	HANDLE_ERROR(cudaMemcpy(dev_lock,cpy_lock,(10+devcount)*sizeof(int),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemset(dev_lock,0,devcount*sizeof(int)));
	/*2.*/
//	printfo(MAXVAL); printf("before\n");

	double bsize = 0;
	for(i = 2; i <= items; i++) {
		//	count =0;
		/*Generate all combinations of cardinality i*/
		count2 = 0;// bsize = 0;
		//	c = (1 << i) -1;
		//	printf("blocks %d\n",(COMBS(c)/BLOCKSIZE)+1);
		for(c = (1 << i) -1; c <= MAXVAL;) {

			double tmp = (double) COMBS(c);
			
			while( bsize <= 128 && tmp > bsize) {
				bsize += 32;
			}
			int blocks =(int)  ceil((tmp/bsize));
			int remainder = blocks - 65535;
	//		while( blocks > 65535 ) {
//				bsize += 32;
//				blocks =(int)  ceil((tmp/bsize));
//			}
			//double bsize = BLOCKSIZE;
			if(remainder > 0) {
				blocks =65535;
				subsetcomp22<<<remainder,bsize>>>(dev_f,dev_o,dev_lock,c,i/2,tmp,count2,65535*bsize);
			
			}

			subsetcomp22<<<blocks,bsize>>>(dev_f,dev_o,dev_lock,c,i/2,tmp,count2,0);
		
			t = c | (c-1);
			c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c) + 1));
			count++;
			count2 +=2;	
			if(count2 < devcount)
				continue;
			HANDLE_ERROR(cudaDeviceSynchronize());		
			HANDLE_ERROR(cudaMemset(dev_lock,0,devcount*sizeof(int)));
//			HANDLE_ERROR(cudaMemcpy(dev_lock,cpy_lock,(devcount)*sizeof(int),cudaMemcpyHostToDevice));
			count2 = 0;
		}

		HANDLE_ERROR(cudaDeviceSynchronize());
		HANDLE_ERROR(cudaMemset(dev_lock,0,devcount*sizeof(int)));
//		HANDLE_ERROR(cudaMemcpy(dev_lock,cpy_lock,(10+devcount)*sizeof(int),cudaMemcpyHostToDevice));
      
		printfo();
	}
	HANDLE_ERROR(cudaDeviceSynchronize());
	HANDLE_ERROR(cudaMemcpy(f,dev_f,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaMemcpy(O,dev_o,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
	//int i;
	HANDLE_ERROR(cudaFree(dev_f));
	HANDLE_ERROR(cudaFree(dev_o));
	HANDLE_ERROR(cudaFree(dev_lock));


	HANDLE_ERROR(cudaDeviceReset());
	free(cpy_lock);
//	printfo(MAXVAL);
	//printf("items %u F[%u] = %u\n",items,MAXVAL,f[MAXVAL-1]);
	parse_wopt(MAXVAL);
	return count;
}



int main(void) {
	/*Start n amount of assets*/
	dint from = 23;
	/*End amount of assets , inclusive*/
	dint till = 23;
	dint MAXVAL = (2 << (from-1));


	time_t start,end,t;
	O = (dint * ) malloc(sizeof(dint)*(2 << (till-1)));
	bids = (dint * ) malloc(sizeof(dint)*(2 << (till-1)));
	f = bids;
	/*Run all tests*/
	for(;from <= till;from++) {
		MAXVAL = (2 << (from-1));
		gen_rand_bids(MAXVAL);
		set_singleton_bid(MAXVAL);
		printf("maxval %u from %u\n",MAXVAL,from);
		  start=clock();//predefined  function in c
		  int count = run_test(MAXVAL,from);
		  end=clock();
		  t=(end-start)/CLOCKS_PER_SEC;
		  printf("\nTime taken =%lu for n= %u with count %d average per count %lf\n", (unsigned long) t,from,count,(double)t/count);
/*Reset the arrays*/
		memset(&f,'\0',sizeof(f));
		memset(&O,'\0',sizeof(O));
	}
	free(O);
	free(f);

	return 0;
}








///old stuff

__global__ void setglobal(unsigned int * f, /*Bid value*/
			  unsigned int * O, /*The move array*/
			  unsigned int * tf,
			  unsigned int * to, /*The configuration*/
			  unsigned int conf /*cardinality of max allowance*/
			  ) {

	unsigned int tid = threadIdx.x;
//	while(to[threadIdx.x] == 0);
//	while(tf[threadIdx.x] == 0);
	extern __shared__ unsigned int share[];
	extern __shared__ unsigned int step[];
	share[tid] = tf[tid];
	step[tid]  = to[tid];

	__syncthreads();
	for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
		if (tid < s) {
			if(share[tid] < share[tid + s]) {
				share[tid] = share[tid+s];
				step[tid] = step[tid+s];
			}
		}
		__syncthreads();
	}

	if(tid == 0) {
		f[conf] = share[0];
		O[conf] = step[0];

		__threadfence();
	}
	
}

__global__ void subsetcomp2(unsigned int * f, /*Bid value*/
			   unsigned int * O, /*The move array*/
			   unsigned int cardset, /*The configuration*/
			    unsigned int cardmax, /*cardinality of max allowance*/
			    unsigned int bidperthread, /*how many bids should be looked at per thread*/
			    unsigned int maxval) 
/*The maximum value it can take, e.g tid can not be greater than 129*/
{
	__shared__ unsigned int _conf;
//	__shared__ unsigned int max[128];
//	__shared__ unsigned int step[128];
	unsigned int i,tid;
	unsigned int subset=0,sum = 0;
	unsigned int tmpset = 0,tmpsum = 0;
	/*thread 0 sets up variables*/
	if(threadIdx.x == 0) {
		/*set up first permutaion eg 0011 for cardset 2*/
		_conf = (1 << cardset) -1;
		
		/*generate the conf value*/
		for(i=0;i<blockIdx.x;i++) {
			tid = _conf | (_conf-1);
			_conf = (tid + 1) | (((~tid & -~tid) - 1) >> (__ffs(_conf)));
		}
		/*set the configuration value to shared memory*/
//		_conf = c;
		/*put it also in the global memory*/

		/*make sure that all blocks sees the change, could possibly discard it*/
		__threadfence();
		O[_conf] = _conf;
		/* printf("threadid.x\t%d\tblockDim.x\t%d\tblockIdx.x\t%d\tI\t%u\n", */
		/*        threadIdx.x, */
		/*        blockDim.x, */
		/*        blockIdx.x, */
		/*        _conf); */
	}
//	max[threadIdx.x] = 0;
//	step[threadIdx.x] = 0;
	__syncthreads();		
//	unsigned int comb = (1 << (__popc(_conf)-1));
	
	tid = threadIdx.x;//*2;
	
	for(i = 1; i<= bidperthread && tid < maxval;i++) {
		tmpset = SUBSET(tid);
		if(__popc(tmpset) <= cardmax) {
			tmpsum = f[tmpset]+f[setdiff(_conf,tmpset)];
			if(sum < tmpsum) {
				subset = tmpset;
				sum = tmpsum;
				
			}
			
		}
		tid += blockDim.x;
		//	printf("subset %u sum %u\n",tmpset,f[tmpset]+f[setdiff(_conf,tmpset)]);
	}

	__syncthreads();

//	unsigned int temp = 0;

	for(i = 0; i < blockDim.x;i++) {
		if(threadIdx.x == i) {
			if(f[_conf] < sum) {
				f[_conf] = sum;
				O[_conf] = subset;
			}
		}
		__syncthreads();
	}
}
__global__ void subsetcomp(unsigned int * f, /*Bid value*/
			   unsigned int * O, /*The move array*/
			   unsigned int conf, /*The configuration*/
			   unsigned int cardinality) /*cardinality of max allowance*/
{
	unsigned int max;
	/*tmp_max is a temporary max variable that is not subject to mutex lock*/
//	__shared__ unsigned int tmp_max;
	__shared__ unsigned int tmpstore[192];
	/* printf("threadid.x\t%d\tblockDim.x\t%d\tblockIdx.x\t%d\tI\t%d\n", */
	/*        threadIdx.x, */
	/*        blockDim.x, */
	/*        blockIdx.x, */
	/*        I); */
	unsigned int subset = (~conf+(I+1))&conf;
	if(__popc(subset) <= cardinality)
		tmpstore[I] = f[setdiff(conf,subset)] + f[subset];
	__syncthreads();
//	__threadfence();
	if(threadIdx.x == 0) {
		unsigned int i = 0;
		unsigned int c = 0;
		for(;i < blockDim.x;i++)
		{
			if(tmpstore[i] > max) {
				max = tmpstore[i];
				c = i;
			}
		}
		if(atomicMax(&f[conf],max) < max) {
			subset = (~conf+((c+blockDim.x*blockIdx.x)+1))&conf;
			atomicExch(&O[conf],subset);
		}
	}
}
__global__ void add(unsigned int * p, unsigned int * f, unsigned int * O)
{
	int tid = blockIdx.x;

	unsigned int conf = p[tid];
	unsigned int card = (unsigned int) __popc(conf)/2;
	unsigned int combinations = 1 << (__popc(conf) -1);
	unsigned int max = f[conf];
	unsigned int set = p[tid];

	unsigned int tmp = 0;
	unsigned int subset;
	unsigned int inverse = ~set;
	unsigned int i;
	if(max == 0) {
		printf("hello");
		return;
	}
/**/
	
	for(i = 1;i<combinations; i++) {
		subset = (inverse+i)&conf;
		if(__popc(subset) > card)
			continue;
		tmp = f[setdiff(conf,subset)] + f[subset];
		if(max < tmp) {
			max = tmp;
			set = subset;
		}
	}
	f[conf] = max;
	O[conf] = set;
}
