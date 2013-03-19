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

//#define MAX (2 << (ITEMS-1))
#if ITEMS < 8
#define dint uint8_t
#elif ITEMS < 16
#undef dint
#define dint uint16_t
#elif ITEMS < 32
#undef dint
#define dint uint32_t
#endif

static void HandleError( cudaError_t err, const char *file,int line ) {
	if (err != cudaSuccess) {
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ),file, line );
		exit( EXIT_FAILURE );
	}
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


/*		           0 1 2 3 4 5 6 7 8			*/
uint32_t  * O;
unsigned short  * f, * bids;
// dint bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
//dint bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
// dint bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6

struct _stack {
	dint conf;
	struct _stack * next;
} typedef stack;


inline  dint cardinality( dint seta) {
	return __builtin_popcount(seta);
}
int indexa =0;
void gen_rand_bids(dint MAXVAL) {
	register dint i = 0;
#if TEST
	unsigned int seed = (unsigned)time ( NULL );
	srand(seed);
	indexa++;
	for(i = 1; i < MAXVAL;i++) {
		bids[i] = 1;//rand()%10+1;
		O[i] = i;
	}

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


#define SUBSET(Y,X)(((~Y+1)+(X))&Y)

#define setdiff(seta,setb) (seta & ~setb)
//maxreg 32
//256
#define MAXBLOCKSIZE (1024)
//32
#define WARPSIZE (32)
//0
#define MIN_BLOCKS_PER_MP 4
#define NAGENTS (25)
//32  
#define NSTREAMS (16)
//2
#define NPERBLOCK (2)
#define confpwarp (2)
//32
#define CONFPKERNEL ((MAXBLOCKSIZE/32)*confpwarp)
//4
#define parasplittings (8)
#define NPARALLELCONF (4)
#define TIMING (0)

#define COMP(Z) {							\
		if(shared_value[tid][Z] < shared_value[tid+i][Z]) {	\
			shared_value[tid][Z] = shared_value[tid+i][Z];	\
		}							\
	}

#define CHECKPOINT(X) {							\
		stop_time = clock();					\
		if(tid == 0 && blockIdx.x == 0) {			\
			total = stop_time - start_time;			\
			printf(X,stop_time - start_time);		\
		}							\
		start_time =clock();					\
	}

#if (TIMING == 0)
#undef CHECKPOINT
#define CHECKPOINT(X) {}
#endif

#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))


template<int blockSize, int overlastval>
__global__ void
__launch_bounds__(blockSize,MIN_BLOCKS_PER_MP)
	subsetcomp33(
		/*0*/	unsigned short * __restrict__ f, /*Bid value*/
		/*5*/	const unsigned int splittings,		
		const unsigned int lastval,//the value a permutation can not exceed.
		const unsigned int card,
		const unsigned int conf1,
		const unsigned int conf2,
		const unsigned int conf3,
		const unsigned int conf4
)
{
	__shared__ unsigned int conf[confpwarp][(blockSize/32)+1];
	__shared__ unsigned short value[confpwarp][(blockSize/32)+1];
	__shared__  unsigned char shift [(blockSize/32)*confpwarp][NAGENTS];// the shift matrix/array
	const unsigned int tid = threadIdx.x;
	const unsigned int laneId = (tid%32);
	const unsigned int warpId = tid/32;
	const unsigned int specsplittings = (!!(splittings%32))+(splittings/32);
//(laneId < splittings)*(splittings/32)+(laneId < (splittings%32));
	const unsigned int initsplit = tid*specsplittings;//+(laneId >= (splittings%32))*(splittings%32);
	// const unsigned int specsplittings = (splittings/32)+!!(splittings%32);
	// const unsigned int initsplit = tid*specsplittings;
	unsigned int leafsplit[2];
//	unsigned int rootsplit[confpwarp];
	unsigned int rvalue[confpwarp][parasplittings];
	/*Thread 0 of each warp*/
	if(!tid) {
		unsigned int tmp;
		unsigned int conftmp;
		if(blockIdx.x == 0) {
			conftmp = conf[0][warpId] = conf1;
		} else if(blockIdx.x == 1) {
			conftmp = conf[0][warpId] = conf2;
		} else if(blockIdx.x == 2) {
			conftmp = conf[0][warpId] = conf3;
		} else {
			conftmp = conf[0][warpId] = conf4;
		}
		int x;
		tmp = conftmp | (conftmp-1);
		conftmp = (tmp + 1) | (((~tmp & -~tmp) - 1) >> (__ffs(conftmp)));
		conf[1][warpId] = conftmp;
		for(x =1; x < (blockSize/32);x++) {
			tmp = conftmp | (conftmp-1);
			conftmp = (tmp + 1) | (((~tmp & -~tmp) - 1) >> (__ffs(conftmp)));
			conf[0][x] = conftmp;
			// if(conftmp > lastval) {
			// 	conf[1][x] = conftmp;
			// 	continue;
			// }
			tmp = conftmp | (conftmp-1);
			conftmp = (tmp + 1) | (((~tmp & -~tmp) - 1) >> (__ffs(conftmp)));
			conf[1][x] = conftmp;
		}
	}
	__syncthreads();
	if(conf[0][warpId] > lastval) {
		return;
	}
	if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1)) && conf[1][warpId] < lastval) {
		if(laneId == 0) {
			unsigned int index;
			unsigned int count = 0;
			unsigned int conftmp = conf[0][warpId] & conf[1][warpId];

			//if(conf[1][warpId] < lastval) {
			value[1][warpId] = f[conf[1][warpId]];
			//}
			value[0][warpId] = f[conf[0][warpId]];
//#pragma unroll
			while(conftmp) {
				index = __ffs(conftmp) - 1; //find which index is first bit
				conftmp &= ~(1 << index);//set nth bit to 0
				shift[warpId][count] = index;
				shift[warpId+1][count] = index;
				count++;
			}
		}
	} else if(laneId < confpwarp) { //generate the shift arrays
		unsigned int index;
		unsigned int conftmp = conf[laneId][warpId];

		if(conftmp < lastval) {
			value[laneId][warpId] = f[conftmp];
		}
#pragma unroll
		for(int x = 0; x < card;x++) {//could put card in template to unroll
			index = __ffs(conftmp) - 1; //find which index is first bit
			conftmp &= ~(1 << index);//set nth bit to 0
			shift[warpId+laneId][x] = index;
		}

	}
	//__syncthreads();



	if(conf[1][warpId] < lastval) {
		unsigned int index;
		unsigned int splittmp= initsplit;
		leafsplit[0] = leafsplit[1] = 0;
		while(splittmp) {
			index = __ffs(splittmp)-1;
			leafsplit[1] += (1 << shift[warpId+1][index]);
			leafsplit[0] += (1 << shift[warpId][index]);//CHECK
			splittmp &= ~(1 << index);
		}
		if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1))) {
			leafsplit[0] = leafsplit[1] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[1]);
		} else {
			leafsplit[1] = SUBSET(conf[1][warpId],leafsplit[1]);
			leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);
		}
	} else {
		unsigned int index;
		unsigned int splittmp= initsplit;
		leafsplit[0] = 0;
		while(splittmp) {
			index = __ffs(splittmp)-1;
			leafsplit[0] += (1 << shift[warpId][index]);//CHECK
			splittmp &= ~(1 << index);
		}
		leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);
	}
	if(!specsplittings) {
		return;
	}

	int x,y;


	for(x = 0; x < specsplittings;x += parasplittings) {

		if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1)) && (conf[1][warpId] < lastval)) {
#pragma unroll 8
			for(y = 0;y < parasplittings ;y++) {
				rvalue[0][y] = 	rvalue[1][y] = 0;
				if(x+y+initsplit < splittings) {
					int tmp = __popc(leafsplit[0]);
					if((NAGENTS-card) <= tmp && (NAGENTS-card) <= (card-tmp)) {
						rvalue[1][y] = rvalue[0][y] = f[leafsplit[0]];
						rvalue[0][y] += f[setdiff(conf[0][warpId],leafsplit[0])];
						rvalue[1][y] += f[setdiff(conf[1][warpId],leafsplit[1])];
					}
					leafsplit[1] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[1]);
					leafsplit[0] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[0]);
				}
			}

		} else {			
#pragma unroll 8			
			for(y = 0;y < parasplittings ;y++) {
				rvalue[0][y] = rvalue[1][y] = 0;
				if(x+y+initsplit < splittings) {
					int tmp = __popc(leafsplit[0]);
					if((NAGENTS-card) <= tmp && (NAGENTS-card) <= (card-tmp)) {
						rvalue[0][y] = f[leafsplit[0]];
						rvalue[0][y] += f[setdiff(conf[0][warpId],leafsplit[0])];
						if(overlastval) {
						if(conf[1][warpId] < lastval) {
							rvalue[1][y] = f[leafsplit[1]];	
							rvalue[1][y] += f[setdiff(conf[1][warpId],leafsplit[1])];
						}
						} else {
							rvalue[1][y] = f[leafsplit[1]];	
							rvalue[1][y] += f[setdiff(conf[1][warpId],leafsplit[1])];
						}
					}
					leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);
					leafsplit[1] = SUBSET(conf[1][warpId],leafsplit[1]);
					
				}
			}

		}

#pragma unroll	       
		for(y = 1;y < parasplittings;y++) {
			if(x+y+initsplit >= splittings) {
				continue;
			}
			if(rvalue[0][0] < rvalue[0][y]) {
				rvalue[0][0] = rvalue[0][y];
			}
			if(rvalue[1][0] < rvalue[1][y]) {
				rvalue[1][0] = rvalue[1][y];
			}
		}


		// 	if(rvalue[1][0] + rvalue[0][0] >= 100) {

		// 	printf("hello %u\n",rvalue[1][0] + rvalue[0][0]);
		// }
		if(__ballot( ( (rvalue[0][0] > value[0][warpId]) || (rvalue[1][0] > value[1][warpId]) ) ) == 0) {
			continue;
		}


		if(__ballot( ( rvalue[0][0] > value[0][warpId] ) ) ) {
			rvalue[0][1] = __shfl_xor((int)rvalue[0][0],16,32);
			rvalue[0][0] = MAX(rvalue[0][0],rvalue[0][1]);
			rvalue[0][1] = __shfl_xor((int)rvalue[0][0],8,32);
			rvalue[0][0] = MAX(rvalue[0][0],rvalue[0][1]);
			rvalue[0][1] = __shfl_xor((int)rvalue[0][0],4,32);
			rvalue[0][0] = MAX(rvalue[0][0],rvalue[0][1]);
			rvalue[0][1] = __shfl_xor((int)rvalue[0][0],2,32);
			rvalue[0][0] = MAX(rvalue[0][0],rvalue[0][1]);
			rvalue[0][1] = __shfl_xor((int)rvalue[0][0],1,32);
			rvalue[0][0] = MAX(rvalue[0][0],rvalue[0][1]);
		}
		
		if(__ballot( ( rvalue[1][0] > value[1][warpId] ) ) ) {
			if(conf[1][warpId] < lastval) {
			rvalue[1][1] = __shfl_xor((int)rvalue[1][0],16,32);
			rvalue[1][0] = MAX(rvalue[1][0],rvalue[1][1]);
			rvalue[1][1] = __shfl_xor((int)rvalue[1][0],8,32);
			rvalue[1][0] = MAX(rvalue[1][0],rvalue[1][1]);
			rvalue[1][1] = __shfl_xor((int)rvalue[1][0],4,32);
			rvalue[1][0] = MAX(rvalue[1][0],rvalue[1][1]);
			rvalue[1][1] = __shfl_xor((int)rvalue[1][0],2,32);
			rvalue[1][0] = MAX(rvalue[1][0],rvalue[1][1]);
			rvalue[1][1] = __shfl_xor((int)rvalue[1][0],1,32);
			rvalue[1][0] = MAX(rvalue[1][0],rvalue[1][1]);
			}
		}
// #pragma unroll
// 		for(int i = 16;i >=1;i >>=1) {
// 				rvalue[0][1] = __shfl_xor((int)rvalue[0][0],i,32);				
// 				if(rvalue[0][1] > rvalue[0][0]) {
// 					rvalue[0][0] = rvalue[0][1];
// 				}

// 				if(!overlastval) {
// 					rvalue[1][1] = __shfl_xor((int)rvalue[1][0],i,32);
// 					if(rvalue[1][1] > rvalue[1][0]) {
// 						rvalue[1][0] = rvalue[1][1];
// 					}
					
// 				} else {
// 					if(conf[1][warpId] < lastval) {
// 						rvalue[1][1] = __shfl_xor((int)rvalue[1][0],i,32);
// 						if(rvalue[1][1] > rvalue[1][0]) {
// 							rvalue[1][0] = rvalue[1][1];
// 						}
// 					}
// 				}
// 		}

		if(laneId == 0) {
			if(rvalue[0][0] > value[0][warpId]) {
				value[0][warpId] = rvalue[0][0];
			}
			if(rvalue[1][0] > value[1][warpId]) {
				value[1][warpId] = rvalue[1][0];
			}
		}
	}

	if(laneId == 0) {
		if(conf[0][warpId] < lastval) {
			if(value[0][warpId] > f[conf[0][warpId]]) {
				f[conf[0][warpId]] = value[0][warpId];
			}
		}
		if(conf[1][warpId] < lastval) {
			if(value[1][warpId] > f[conf[1][warpId]]) {
				f[conf[1][warpId]] = value[1][warpId];
			}
		}
	}
	return;
}



#define COMBS(X) ((1 << cardinality(X)-1) - 1)

int run_test(unsigned int MAXVAL,dint items) {

	register unsigned int i,c1,count =0;
	unsigned short *dev_bids;

	count = 0;

	HANDLE_ERROR(cudaDeviceReset());
	HANDLE_ERROR(cudaSetDeviceFlags(cudaDeviceScheduleYield));
  	HANDLE_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
//	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitMallocHeapSize,0));
	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitStackSize,0));
	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitMallocHeapSize,0));
	HANDLE_ERROR(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));

	register unsigned int streams = NSTREAMS;
	register unsigned int streamcount = 0;
	register cudaStream_t stream[streams];
	for(i = 0;i < streams; i++)
		HANDLE_ERROR(cudaStreamCreate(&stream[i]));

	count = 0;

 	HANDLE_ERROR(cudaMalloc((void **)&dev_bids, MAXVAL*sizeof(short)));

 	HANDLE_ERROR(cudaMemcpy(dev_bids,bids,MAXVAL*sizeof(short),cudaMemcpyHostToDevice));

	register unsigned int bsize = MAXBLOCKSIZE;
	register int blocks;
	int prev =0;
//	lock_count = 0;
	time_t rstart,rend,rt;
	rstart=clock();
	for(i = 2; i <= items; i++) {
		time_t start,end,t;

		start=clock();
		unsigned int splittings;
		blocks =4;//(int)  ceil((threads/bsize));
		double threads;
		c1 = (1 << i) -1;
		unsigned int c2 = c1;
		unsigned int c3;
		unsigned int ca[blocks];
		unsigned int cacount;
		splittings =  COMBS(c1);///NPERBLOCK;
		threads = ((double) splittings)/ NPERBLOCK;
		threads = ceil(threads);
		
		for(; c1 <= MAXVAL;) {

			cacount = 1;
			ca[0] = c1;			
			for(int x = 0, cacount = 1; x < CONFPKERNEL*blocks;x++) {
				t = c2 | (c2-1);
				c2 = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c2) + 1));
				if(x%CONFPKERNEL == 0 && x > 1) {
					ca[cacount] = c2;
					cacount++;
				}
			}
	
			switch((c2 > MAXVAL)) {
			case 1:
				subsetcomp33 < MAXBLOCKSIZE , 1> <<<blocks,MAXBLOCKSIZE,0,stream[streamcount]>>>(dev_bids,splittings,MAXVAL,i,ca[0],ca[1],ca[2],ca[3]);
				break;
			case 0:
				subsetcomp33 < MAXBLOCKSIZE , 0> <<<blocks,MAXBLOCKSIZE,0,stream[streamcount]>>>(dev_bids,splittings,MAXVAL,i,ca[0],ca[1],ca[2],ca[3]);
				break;

			}
			c1 = c2;


			streamcount++;
			count++;

			if(streamcount >= streams)
				streamcount = 0;
		}
		

		for (int t = 0; t < streams; ++t) {
			HANDLE_ERROR(cudaStreamSynchronize(stream[t]));
		}

		HANDLE_ERROR(cudaDeviceSynchronize());


		end=clock();
		t=(end-start)/(CLOCKS_PER_SEC/1000);
		printf("ended card %d blocks\t %d threads/block %u, n kernels %u \t time %lu \t splittings %d time per kernel %u\n",i,blocks,bsize,count-prev,t,splittings,t/(count-prev));
		prev =	count;

	}
	for (int i = 0; i < streams; ++i)
		cudaStreamDestroy(stream[i]);

	HANDLE_ERROR(cudaDeviceSynchronize());

	rend=clock();
	rt=(rend-rstart)/(CLOCKS_PER_SEC/1000);
	printf("real time %lu\n",rt);

	HANDLE_ERROR(cudaMemcpy(f,dev_bids,MAXVAL*sizeof(short),cudaMemcpyDeviceToHost));
//	HANDLE_ERROR(cudaMemcpy(O,dev_o,MAXVAL*sizeof(int),cudaMemcpyDeviceToHost));
	//int i;
	HANDLE_ERROR(cudaFree(dev_bids));
//	HANDLE_ERROR(cudaFree(dev_o));
	// HANDLE_ERROR(cudaFree(dev_lock1));
	// HANDLE_ERROR(cudaFree(dev_lock2));

	HANDLE_ERROR(cudaDeviceReset());
//	printfo(MAXVAL);
	return count;
}

dint max2(dint conf) {
     register dint card = cardinality(conf)/2;
     register dint combinations = (1 << cardinality(conf)-1)-1;
     register dint max = f[conf];
     register dint tmp = 0;
     register dint subset = 0;
     register const dint inverse = ~conf;
     register dint i;
     for(i = 1;i<=combinations; i++) {
	     subset = ((inverse+1)+subset)&conf;
	     tmp = f[setdiff(conf,subset)] + f[subset];
	     if(max == tmp) {
		     break;
		     //  return subset;
	     }

     }
     return subset;
}

int recur_parse_wopt(dint MAXVAL) {
	stack * root = (stack *) malloc(sizeof(stack));
	stack * sroot = NULL;
	stack * scurr = NULL;
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
//		printf("curr %u\t\n",curr->conf);
		if(f[conf] != bids[conf]) {
			dint proper_subset = max2(conf);
			dint diff = setdiff(conf,proper_subset);
			curr->conf = proper_subset;
			stack * tmp = (stack *) malloc(sizeof(stack));
			//printf("diff %u\t conf %u\t O[diff] %u\t O[conf]\t f %u\n",diff,conf,O[diff],O[conf],f[conf]);
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
			//printf("conf %u value %u\n",curr->conf,bids[curr->conf]);

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

int main(void) {
	/*Start n amount of assets*/
	dint from = NAGENTS;
	dint MAXVAL = (2 << (from-1));

	time_t start,end,t;
	O = (dint * ) malloc(sizeof(dint)*(2 << (from-1)));
	bids = (unsigned short * ) malloc(sizeof(short)*(2 << (from-1)));
	f = (unsigned short * ) malloc(sizeof(short)*(2 << (from-1)));
//	f = bids;
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
		ret_val= recur_parse_wopt(MAXVAL);// parse_wopt(MAXVAL);
		printf("\nTime taken =%lu for n= %u with count %d average per count %lf\n", (unsigned long) t,from,count,(double)t/count);
	}


	free(O);
	free(f);

	return 0;
}


// template<int blockSize,int nparallelconf,int confpkernel,int nperblock,int currblocksize>
// __global__ void
// __launch_bounds__(currblocksize,MIN_BLOCKS_PER_MP)
// 	subsetcomp32(
// 		/*0*/	unsigned short * __restrict__ f, /*Bid value*/
// 		/*1*/	unsigned int * __restrict__ O, /*The move array*/
// 		/*2*/	unsigned int * __restrict__ lock,
// 		/*5*/	unsigned int maxval,
// 		/*6*/	unsigned short count1,
// 		unsigned int conf1,
// 		unsigned int lastval,//the value a permutation can not exceed.
// 		unsigned int card)
// {
// 	//confpkernel = how many configurations the kernel will evaluate
// 	//nparallelconf = how many configurations the kernel will evaluate at the same time
// 	__shared__  unsigned short shared_value[(currblocksize >> 5)+1][nparallelconf];
// 	__shared__  unsigned int shared_conf[(currblocksize >> 5)+1][nparallelconf];
// 	__shared__  unsigned int conf[confpkernel];// the configurations needed for the whole execution
// 	__shared__  unsigned short shift [confpkernel][NAGENTS];// the shift matrix/array
// 	//__shared__ volatile unsigned int tmp [confpkernel][2];
//  	register unsigned int subset_value[2][8];//the value for one of the subset sums
// 	register unsigned int subset_conf[2][8];

// 	register unsigned int count = count1;
// 	register unsigned int const tid = threadIdx.x;
// //	__shared__ unsigned int
// 	register unsigned int ispec = nperblock*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;
// 	register int x,i,z; //counter

// #if (TIMING == 1)
// 	clock_t stop_time, total,start_time;
// 	total = 0;
// 	start_time = clock();
// #endif

// 	if(tid == 0) {
// 		conf[0] = conf1;
// #pragma unroll
// 		for(x = 1;x < confpkernel; x++) { // generate the configurations
// 			conf[x] = conf[x-1];
// 			z = conf[x] | (conf[x]-1);
// 			conf[x] = (z + 1) | (((~z & -~z) - 1) >> (__ffs(conf[x])));
// 		}
// 		CHECKPOINT("l1 %d\n");
// 	}
// 	if(confpkernel > 32) {
// 		__syncthreads();
// 	}

// 	count = count1;
// 	if(tid < confpkernel) { //generate the shift arrays
// 		subset_conf[0][0] = conf[tid]; // re-use registers tmpval
// 		subset_conf[0][1] = 0;// re-use registers index
// #pragma unroll
// 		for(x = 0,i=0; x < card;x++) {//could put card in template to unroll
// 			subset_conf[0][1] = __ffs(subset_conf[0][0]) - 1; //find which index is first bit
// 			subset_conf[0][0] &= ~(1 << subset_conf[0][1]);//set nth bit to 0
// 			shift[tid][x] = subset_conf[0][1];
// 		}
// 	}

// 	CHECKPOINT("l2 %d\n");

// 	__syncthreads();

// #pragma unroll
// 	for(x =0; x < confpkernel; x += nparallelconf) {

// 		if(conf[x] > lastval) {//if the permutation is larger than the full set, c> (1 << NAGENTS)
// 			continue;
// 		}
// //#pragma unroll
// //		for(i = 0; i < nperblock; i++) {
// #pragma unroll
//  		for(z=0; z < nparallelconf;z++) {
//  			subset_value[0][z] = subset_value[1][z] = 0U;
//  		}
// //		}
// 		CHECKPOINT("l3 %d\n");
// 		if(ispec >= maxval) {
// 			goto postfetch;
// 		}
// 			//This for loop initilize the first subset configuration.

// #pragma unroll
// 		for(i = 0; i < nparallelconf; i++) {
// 			unsigned int tmp = ispec;
// 			//unsigned int const tcar = __popc(tmp);
// 			subset_conf[0][i] = 0;
// 			while(tmp) {
// 				unsigned short index = __ffs(tmp)-1;
// 				subset_conf[0][i] += (1 << shift[x+i][index]);//CHECK
// 				tmp &= ~(1 << index);
// 			}
// 		}
// 		CHECKPOINT("l4 %d\n");
// #pragma unroll
// 		for(z=0;z < nparallelconf;z++) {
// 			if(conf[z+x] > lastval) {
// 				continue;
// 			}
// 			subset_conf[0][z] = SUBSET(conf[z+x],subset_conf[0][z]);
// 			subset_value[0][z] = f[(setdiff(conf[z+x],subset_conf[0][z]))] + f[subset_conf[0][z]];
// 			//ispec++;
// 			if((ispec+1) >= maxval) {
// 				continue;
// 			}
// 			subset_conf[1][z] = SUBSET(conf[z+x],subset_conf[0][z]);
// 			subset_value[1][z] = f[(setdiff(conf[z+x],subset_conf[1][z]))] + f[subset_conf[1][z]];
// 		}
// 		CHECKPOINT("l6 %d\n");

// 	postfetch:

// #pragma unroll
// 		for(z = 0; z < nparallelconf;z++) {//warp reduction
// 			if(subset_value[1][z] > subset_value[0][z]) {
// 				subset_value[0][z] = subset_value[1][z];
// 				subset_conf[0][z] = subset_conf[1][z];
// 			}
// #pragma unroll
// 			for(i = 16;i >=1;i >>=1) {
// 				int warp_value = __shfl_xor((int)subset_value[0][z],i,32);
// 				int warp_conf = __shfl_xor((int)subset_conf[0][z],i,32);
// 				if(warp_value > subset_value[0][z]) {
// 					subset_value[0][z] =(unsigned int) warp_value;
// 					subset_conf[0][z] =(unsigned int) warp_conf;
// 				}
// 			}
// 			//tid&(WARPSIZE-1) == tid%WARPSIZE
// 			//Only threads with line id == 0 is allowed to update in the shared memory,
// 			//i.e. the first thread in each warp
// 			if(!(tid&(31))) {
// 				unsigned int index = tid >> 5; // tid >> 5 == tid / 32 which warp it is
// 				shared_value[index][z] = subset_value[0][z];
// 				shared_conf[index][z] = subset_conf[0][z];
// 			}

// 		}
// 		CHECKPOINT("l7 %d\n");

// 		//	CHECKPOINT("l8 %d\n");
// 		if((currblocksize/32) > 1) {
//  		__syncthreads();
// 		}
// 		//how many warps is it, block dimension divided by warp size
// 		//e.g. 256/32 == 256 >> 5
// 		if((currblocksize/32) > 1) {//evaluated by the pre-processor
// 			i = (currblocksize >> 6);//blockDim.x >> 6;
// 			if(tid<i) {//reduction mby move down if you get wrong results gained ~1000 cycles

// #pragma unroll
// 			for(; i > 0; i >>= 1) {
// #pragma unroll
//   					for(z=0; z < nparallelconf;z +=4) {
// 						COMP(z);
//   					}
// 			}
// 			}
// 		}

// 		CHECKPOINT("l9 %d\n");
// 		if((currblocksize) >= blockSize) {//evaluated by the pre-processor
// 			//__syncthreads();
// 			if(tid == 0) {
// #pragma unroll
// 				for(z=0; z < nparallelconf;z++) {
// 					if(f[conf[z+x]] < shared_value[0][z]) {
// 						//	printf("lock val %u shared_val %s\n",lock[count+z] ,shared_value[0][z]);
// 						if(atomicMax(&(lock[count+z]),shared_value[0][z]) < shared_value[0][z]) {
// 							//	O[conf[z+x]] = shared_conf[0][z];
// 							f[conf[z+x]] = shared_value[0][z];
// 						}
// 					}
// 				}
// 			}
// 		} else {
// #pragma unroll
// 			for(z=0; z < nparallelconf;z++) {
// 				if(f[conf[z+x]] < shared_value[0][z]) {
// 					//	O[conf[z+x]] = shared_conf[0][z];
// 					f[conf[z+x]] = shared_value[0][z];

// 				}

// 			}
// 		}

// 		count += nparallelconf;
// 		ispec = NPERBLOCK*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;

// 	}
// }
// template<int blockSize,int nparallelconf,int confpkernel,int nperblock,int currblocksize>
// __global__ void
// __launch_bounds__(currblocksize,MIN_BLOCKS_PER_MP)
// 	subsetcomp35(
// 				/*0*/	unsigned short * __restrict__ f, /*Bid value*/
// 		/*1*/	unsigned int * __restrict__ O, /*The move array*/
// 		/*2*/	unsigned int * __restrict__ lock,
// 		/*5*/	unsigned int maxval,
// 		/*6*/	unsigned short count1,
// 		unsigned int conf1,
// 		unsigned int lastval,//the value a permutation can not exceed.
// 		unsigned int card)
// {
// 	//confpkernel = how many configurations the kernel will evaluate
// 	//nparallelconf = how many configurations the kernel will evaluate at the same time
// 	__shared__  unsigned short shared_value[(currblocksize >> 5)+1][nparallelconf];
// 	__shared__  unsigned int conf[confpkernel];// the configurations needed for the whole execution
// //	__shared__ unsigned short old_values[]
// 	__shared__  unsigned char shift [confpkernel][NAGENTS];// the shift matrix/array
// 	//__shared__ volatile unsigned int tmp [confpkernel][2];
//  	register unsigned int subset_value[nperblock][nparallelconf];//the value for one of the subset sums
// 	register unsigned int subset_conf[2];//[nparallelconf];
// //	register unsigned int count = count1;
// 	register unsigned int const tid = threadIdx.x;
// 	register unsigned int const ispec = nperblock*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;
// 	register int x,i,z; //counter

// #if (TIMING == 1)
// 	clock_t stop_time, total,start_time;
// 	total = 0;
// 	start_time = clock();
// #endif

// 	if(tid == 0) {

// 		subset_conf[0] = conf[0] = conf1;
// #pragma unroll
// 		for(x = 1;x < confpkernel; x++) { // generate the configurations
// 			//conf[x] = conf[x-1];
// 			z = subset_conf[0] | (subset_conf[0]-1);
// 			subset_conf[0] = (z + 1) | (((~z & -~z) - 1) >> (__ffs(subset_conf[0])));
// 			conf[x]= subset_conf[0];
// 		}
// 		CHECKPOINT("l1\t %d\n");
// 	}
// //	__syncthreads();

// #if (TIMING == 1)
// 	start_time = clock();
// #endif
// 	/* Generates the shift array in order in the future to create subsets.
// 	 *
// 	 */

// 	if(tid < confpkernel) { //generate the shift arrays

// 		subset_conf[0] = conf[tid]; // re-use registers tmpval
// 		subset_conf[1] = 0;// re-use registers index
// #pragma unroll
// 		for(x = 0,i=0; x < card;x++) {//could put card in template to unroll
// 			subset_conf[1] = __ffs(subset_conf[0]) - 1; //find which index is first bit
// 			subset_conf[0] &= ~(1 << subset_conf[1]);//set nth bit to 0
// 			shift[tid][x] = subset_conf[1];
// 		}

// 	}

// 	CHECKPOINT("l2\t %d\n");

// 	__syncthreads();


// #pragma unroll
// 	for(x =0; x < confpkernel; x += nparallelconf) {

// #if (TIMING == 1)
// 		start_time = clock();
// #endif
// 		if(conf[x] > lastval) {//if the permutation is larger than the full set, c> (1 << NAGENTS)
// 			continue;
// 		}

// #pragma unroll
//  		for(z=0; z < nparallelconf;z++) {
//  			subset_value[0][z] = subset_value[1][z] = 0U;
//  		}
// 		CHECKPOINT("l3\t %d\n");
// 		if(ispec >= maxval) {
// 			goto postfetch;
// 		}

// 		CHECKPOINT("l3\t %d\n");
// 		unsigned int tmp[2];
// 		tmp[1]= tmp[0] = 0;
// #pragma unroll
// 		for(z=0;z < nparallelconf;z +=2) {

// 			if(conf[z+x] >= lastval) {
// 				continue;
// 			}
// 			unsigned int tmpstore = 0, tmpstore1 = 0;
// 			unsigned int index;
// 			//unsigned int const tcar = __popc(tmp);
// 			unsigned int wshile;
// 			unsigned int icard;

// 			if(conf[z+x+1] < lastval) {
// 				wshile = ispec;
// 				subset_conf[0] = 0;
// 				subset_conf[1] = 0;
// 				while(wshile) {
// 					index = __ffs(wshile)-1;
// 					subset_conf[1] += (1 << shift[x+z+1][index]);
// 					subset_conf[0] += (1 << shift[x+z][index]);//CHECK
// 					wshile &= ~(1 << index);
// 				}
// 				subset_conf[1] = SUBSET(conf[z+x+1],subset_conf[1]);
// 				// icard = __popc(subset_conf[1]);
// 				// if(icard >= (NAGENTS-card)) {
// 				subset_value[0][z+1] = f[subset_conf[1]];
// 				// } else {
// 				// 	subset_value[0][z+1] = 0;
// 				// }
// 			} else {
// 				wshile = ispec;
// 				subset_conf[0] = 0;
// 				while(wshile) {
// 					index = __ffs(wshile)-1;
// 					subset_conf[0] += (1 << shift[x+z][index]);//CHECK
// 					wshile &= ~(1 << index);
// 				}
// 			}
// 			subset_conf[0] = SUBSET(conf[z+x],subset_conf[0]);
// 			// icard = __popc(subset_conf[0]);
// 			// if(icard > (NAGENTS-card)) {
// 			subset_value[0][z] = f[subset_conf[0]];
// 			// } else {
// 			//subset_value[0][z] = 0;
// 			// }
// 			if((ispec+1) < maxval) {
// 				tmp[0] = SUBSET(conf[z+x],subset_conf[0]);
// 				// if(tmp[0] == subset_conf[1]) {\\expensive
// 				// 	subset_value[1][z] = subset_value[0][z];
// 				// } else
// 				if(tmp[0] != tmp[1]) {//dont remove
// 					subset_value[1][z] = f[tmp[0]];
// 				}

// 				if(conf[z+x+1] < lastval) {
// 					tmp[1] = SUBSET(conf[z+x+1],subset_conf[1]);
// 					if(tmp[1] == tmp[0]) {// Do not remove, cost 1 second
// 						subset_value[1][z+1] = subset_value[1][z];
// 					// } else//  if(tmp[1] == subset_conf[0]) {
// 					// 	subset_value[1][z+1] = subset_value[0][z];
// 					} else {
// 						subset_value[1][z+1] = f[tmp[1]];
// 					}
// 				}
// 			}


// 				//	}

// 			subset_conf[0] = (setdiff(conf[z+x],subset_conf[0]));
// 			tmpstore = f[subset_conf[0]];//

// 			//next splitting
// 			if((ispec+1) < maxval) {
// 				tmp[0] = (setdiff(conf[z+x],tmp[0]));
// 				tmpstore1 = f[tmp[0]];
// 			}

// 			//next configuration
// 			if(conf[z+x+1] < lastval) {
// 				subset_conf[1] = (setdiff(conf[z+x+1],subset_conf[1]));
// 				if(z < (nparallelconf-2)){
// 				 	subset_value[1][z+2] = subset_value[1][z+1];
// 				}
// 				if(subset_conf[1] == subset_conf[0]){// dont remove
// 					subset_value[0][z+1] += tmpstore;
// 					//} // else if(subset_conf[1] == tmp[0]){ // expensive
// 				// 	subset_value[0][z+1] += tmpstore1;
// 				} else{
// 					subset_value[0][z+1] += f[subset_conf[1]];
// 				}
// 				//next splitting
// 				if((ispec+1) < maxval) {
// 					tmp[1] = (setdiff(conf[z+x+1],tmp[1]));
// 					// if(tmp[1] == tmp[0]){// expensive
// 					// 	subset_value[1][z+1] += tmpstore1;
// 					// }else
// 					subset_value[1][z+1] += f[tmp[1]];
// 				}
// 			}
// 			subset_value[0][z] += tmpstore;//f[(setdiff(conf[z+x],subset_conf[0]))];//
// 			subset_value[1][z] += tmpstore1;
// 			tmp[1] = (setdiff(conf[z+x+1],tmp[1]));
// 			//	subset_conf[1] = (setdiff(conf[z+x+1],subset_conf[1]));


// 		}
// 		CHECKPOINT("l4\t %d\n");

// 	postfetch:

// #pragma unroll
// 		for(z = 0; z < nparallelconf;z +=2) {//warp reduction
// 			if(conf[z+x] >= lastval) {
// 				continue;
// 			}
// 			if(subset_value[1][z] > subset_value[0][z]) {
// 				subset_value[0][z] = subset_value[1][z];
// 				//subset_conf[0][z] = subset_conf[1][z];
// 			}
// 			if(subset_value[1][z+1] > subset_value[0][z+1]) {
// 				subset_value[0][z+1] = subset_value[1][z+1];
// 				//subset_conf[0][z] = subset_conf[1][z];
// 			}
// #pragma unroll
// 			for(i = 16;i >=1;i >>=1) {
// 				subset_value[1][z] = __shfl_xor((int)subset_value[0][z],i,32);
// 				subset_value[1][z+1] = __shfl_xor((int)subset_value[0][z+1],i,32);
// 				//int warp_conf = __shfl_xor((int)subset_conf[0][z],i,32);
// 				if(subset_value[1][z] > subset_value[0][z]) {
// 					subset_value[0][z] = subset_value[1][z];// (unsigned int) warp_value;
// 					//	subset_conf[0][z] =(unsigned int) warp_conf;
// 				}
// 				if(subset_value[1][z+1] > subset_value[0][z+1]) {
// 					subset_value[0][z+1] = subset_value[1][z+1];
// 				//subset_conf[0][z] = subset_conf[1][z];
// 				}
// 			}
// 			//tid&(WARPSIZE-1) == tid%WARPSIZE
// 			//Only threads with line id == 0 is allowed to update in the shared memory,
// 			//i.e. the first thread in each warp
// 			if(!(tid&(31))) {
// 				unsigned int index = tid >> 5; // tid >> 5 == tid / 32 which warp it is
// 				shared_value[index][z] = subset_value[0][z];
// 				shared_value[index][z+1] = subset_value[0][z+1];
// 				//shared_conf[index][z] = subset_conf[0][z];
// 			}

// 		}
// 		/*WORKING STOP DELETE*/
// 		CHECKPOINT("l5\t %d\n");

// 		//	CHECKPOINT("l8\t %d\n");
// 		if((currblocksize/32) > 1) {
//  		__syncthreads();
// 		}
// 		//how many warps is it, block dimension divided by warp size
// 		//e.g. 256/32 == 256 >> 5
// //		if((currblocksize/32) > 1) {//evaluated by the pre-processor

//  			i = blockDim.x >> 6;
//  			if(tid<i) {//reduction mby move down if you get wrong results gained ~1000 cycles
//  #pragma unroll
//  				for(; i > 0; i >>= 1) {
// //#pragma unroll
// 						for(z=0; z < nparallelconf;z ++) {

// // 						subset_value[0][z] = shared_value[tid][z];
// // //#pragma unroll
// // 						for(i = (currblocksize/32);i >=1;i >>=1) {
// // 							int warp_value = __shfl_down((int)subset_value[0][z],1,32);
// // 							//int warp_conf = __shfl_xor((int)subset_conf[0][z],i,32);
// // 							if(warp_value > subset_value[0][z]) {
// // 								subset_value[0][z] =(unsigned int) warp_value;
// // 								//	subset_conf[0][z] =(unsigned int) warp_conf;
// // 							}
// // 						}
// // 						if(tid == 0) {
// // //							unsigned int index = tid >> 5; // tid >> 5 == tid / 32 which warp it is
// // 							shared_value[0][z] = subset_value[0][z];
// // 							//shared_conf[index][z] = subset_conf[0][z];
// // 						}
// 						COMP(z);
// 						// COMP(0+1);
// 						// COMP(0+2);
// 						// COMP(0+3);

// 							}
// // 				}
//  			}
// 		}

// 		CHECKPOINT("l6\t %d\n");
// 		if((currblocksize) >= blockSize) {//evaluated by the pre-processor
// 			//__syncthreads();
// 			if(tid == 0) {
// #pragma unroll
// 				for(z=0; z < nparallelconf;z++) {
// 					if(conf[z+x] > lastval) {continue;}
// 					if(f[conf[z+x]] < shared_value[0][z]) {
// 						//	printf("lock val %u shared_val %s\n",lock[count+z] ,shared_value[0][z]);
// 						if(atomicMax(&(lock[count1+x+z]),shared_value[0][z]) < shared_value[0][z]) {
// 							//O[conf[z+x]] = shared_conf[0][z];
// 							if(f[conf[z+x]] < shared_value[0][z]) {
// 								f[conf[z+x]] = shared_value[0][z];
// 							}
// 						}
// 					}
// 				}
// 			}
// 		} else {
// #pragma unroll
// 			for(z=0; z < nparallelconf;z++) {
// 				if(conf[z+x] > lastval) {continue;}
// 				if(f[conf[z+x]] < shared_value[0][z]) {
// 					//O[conf[z+x]] = shared_conf[0][z];
// 					f[conf[z+x]] = shared_value[0][z];

// 				}

// 			}
// 		}
// 		CHECKPOINT("l7\t %d\n");
// //		count += nparallelconf;
// 	}
// }
