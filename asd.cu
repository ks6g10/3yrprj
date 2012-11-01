
#include <stdio.h>
#include <string.h> // for ffs
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cuda.h>

/*Debug enabled gives more print statements of bids and how the "Matrix" gets evaluated*/
#define DEBUG 0
/*Test sets all bids to one, which should give you n=|ITEMS| bids on output*/
#define TEST 1
/*Defines from 0-Range the random will give out*/
#define RANGE 20
#define ITEMS 23
#define MAX (2 << (ITEMS-1))
#if ITEMS < 8
#define dint uint8_t
#elif ITEMS < 16
#define dint uint16_t
#elif ITEMS < 32
#define dint uint32_t
#elif ITEMS > 32
#define dint uint64_t
#endif

static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

/*                    1         2        4       8        16*/

char * assets[3] = {"apple","mapple","potato"};
/*		           0 1 2 3 4 5 6 7 8			*/
dint bids[MAX];
// dint bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
//dint bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
// dint bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6
dint wopt[MAX] = {0};
dint f[MAX] = {0};
dint O[MAX] = {0};

inline  dint intersect( dint seta,  dint setb) {
     return (seta & setb);
}

inline  dint _union( dint seta,  dint setb) {
     return (seta | setb);
}

#define setdiff(seta,setb) (seta & ~setb)

inline  dint cardinality( dint seta) {
     return __builtin_popcount(seta);
}

void gen_rand_bids(dint MAXVAL) {
     register dint i = 0;
#if TEST
     for(i = 1; i < MAXVAL;i++) {
	  bids[i] = 1;
     }
#else
     for(i = 1; i < MAX;i++) {
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

/*Reminder of sets
 *f[]
 *O[]
 *bids[]
 */
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



struct _stack {
	dint conf;
	struct _stack * next;
} typedef stack;


void parse_wopt(dint MAXVAL) {
     //wopt at start contain MAX at wopt[0] which is the combination that goes in bids[wopt[n]]
     stack * root =(stack *) malloc(sizeof(stack));
     //DO NOT REMOVE -1
     root->conf = (MAXVAL)-1;
     stack * curr = root;
     while(curr != NULL) {

	     dint conf = curr->conf;
	     if(conf == 0) {
	    	 printf("EXIT FAILURE\n");
	    	 return;
	     }
	     if(conf != O[conf]) {
		     dint diff = setdiff(conf,O[conf]);
		     curr->conf = O[conf];
		     stack * tmp = (stack *) malloc(sizeof(stack));
		     tmp->conf = diff;
		     tmp->next = curr;
		     root = tmp;
		     curr = root;
		     continue;
	     }
	     curr = curr->next;
     }
     curr = root;
     dint tmp = 0;
     while(curr != NULL) {
	  printf("conf %u value %u\n",curr->conf,bids[curr->conf]);
	  tmp++;
	  stack * tmp = curr;
	  curr = curr->next;
	  free(tmp);
     }
     printf("n = %u",tmp);
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















__global__ void gen_subset(void) {



}

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

void run_test(dint MAXVAL,dint items) {
/*Setup the environment*/
     dint perm[MAXVAL];
     gen_rand_bids(MAXVAL);
     set_singleton_bid(MAXVAL);
     printfo();
     dint i, c;

 	dint *dev_f,*dev_o,*dev_p;
 	HANDLE_ERROR(cudaMalloc((void **)&dev_f, MAXVAL*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_o, MAXVAL*sizeof(int)));
 	HANDLE_ERROR(cudaMalloc((void **)&dev_p, MAXVAL*sizeof(int)));
 	//add<<<1,1>>>(2,7,dev_c);

 	HANDLE_ERROR(cudaMemcpy(dev_f,bids,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));
 	HANDLE_ERROR(cudaMemcpy(dev_o,O,MAXVAL*sizeof(int),cudaMemcpyHostToDevice));
 //	printf("2+7 = %d \n",c);
     /*2.*/
 	 int count = 0;
/*     items represent the maximum cardinality a set can have*/
     for(i = 2; i <= items; i++) {
    	 count = 0;
/*1.*/	     for(c = (1 << i) -1; c <= MAXVAL;) {
	    	 perm[count++] = c;
		     //bit hacks "Compute the lexicographically next bit permutation"
		     dint t = c | (c-1);
		     c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c) + 1));
		     //end ref
	     }
	     HANDLE_ERROR(cudaMemcpy(dev_p,perm,count*sizeof(int),cudaMemcpyHostToDevice));
	     add<<<count,1>>>(dev_p,dev_f,dev_o);
	     printfo();
     }
     HANDLE_ERROR(cudaMemcpy(f,dev_f,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
     HANDLE_ERROR(cudaMemcpy(O,dev_o,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
     //int i;
     cudaFree(dev_f);
     cudaFree(dev_o);
     cudaFree(dev_p);

     parse_wopt(MAXVAL);

}



int main(void) {
     /*Start n amount of assets*/
     dint from = 14;
     /*End amount of assets, inclusive*/
     dint till = 14;
     dint MAXVAL = (2 << (from-1));
     if(till > ITEMS) {
	  printf("More than maximum allowed\n");
	  return 1;
     }

     //time_t start,end,t;

     /*Run all tests*/
     for(;from <= till;from++) {
	  MAXVAL = (2 << (from-1));
	  printf("maxval %u\n",MAXVAL);
	//  start=clock();//predefined  function in c
	  run_test(MAXVAL,from);
	//  end=clock();
	//  t=(end-start)/CLOCKS_PER_SEC;
	//  printf("\nTime taken =%lu for n= %u\n", (unsigned long) t,from);
/*Reset the arrays*/
	  memset(&f,'\0',sizeof(f));
	  memset(&O,'\0',sizeof(O));
     }

     return 0;
}

void old_test(void) {
     /*Testing facility*/
     dint i;
     for(i = 1; i < 8; i++) {
	  dint z = intersect(i,i-1);
	  dint x = _union(i,i-1);
	  dint f = setdiff(i,i-1);
	  dint t = cardinality(i);
	  printf("i %u \tinter %u\t union %u\t diff %u\t card %u\n",i,z,x,f,t);
     }
     return;
}


//#include <cuda.h> // helper utility functions
//#include <cuda_runtime.h>
//#include <cuda.h>
//#include <book.h>





__global__ void kernel(void){

}

int m2ain(void)
{

	return 0;
}
