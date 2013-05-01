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
		//	O[i] = i;
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
		//	O[i] = i;
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
	}
}



void printfo(dint MAXVAL) {
	int i;
	printf("\n");
	for(i = 1; i< MAXVAL; i++)
	{
		//printf("Bid[%d]\t%u\tF[%d]\t%u\tO[%d]\t%u\tbin\t%s\n",i,bids[i],i,f[i],i,O[i],btb(i));

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
#define NAGENTS (20)
//32  
#define NSTREAMS (16)
//2
#define NPERBLOCK (2)
#define confpwarp (2)
//32
#define CONFPKERNEL ((MAXBLOCKSIZE/32)*confpwarp)
//4
#define parasplittings (32)
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
		const unsigned int splittings,
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
	const unsigned int laneId = (tid&31);
	const unsigned int warpId = tid/32;
	const unsigned int specsplittings = (splittings/32)+!!(splittings&31);
//(laneId < splittings)*(splittings/32)+(laneId < (splittings%32));,COMBS(29),(!!(COMBS(29)%32))+(COMBS(29)/32)
	const unsigned int initsplit = laneId*specsplittings;//+(laneId >= (splittings%32))*(splittings%32);
	
	// const unsigned int initsplit = tid*specsplittings;
	unsigned int leafsplit[2];
//	unsigned int rootsplit[confpwarp];
	unsigned int rvalue[confpwarp][2];
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
			if(conftmp > lastval) {
				conf[1][x] = conftmp;
				continue;
			}
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
		unsigned int splittmp = initsplit;
		leafsplit[0] = leafsplit[1] = 0;
		if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1))) {
			while(splittmp) {
			index = __ffs(splittmp)-1;
			leafsplit[1] += (1 << shift[warpId+1][index]);			
			splittmp &= ~(1 << index);
			}
			leafsplit[0] = leafsplit[1] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[1]);
		} else {
			while(splittmp) {
			index = __ffs(splittmp)-1;
			leafsplit[1] += (1 << shift[warpId+1][index]);
			leafsplit[0] += (1 << shift[warpId][index]);//CHECK
			splittmp &= ~(1 << index);
		}
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

	int y;


//	for(x = 0; x < specsplittings;x += parasplittings) {
	rvalue[0][0] =rvalue[0][1] = rvalue[1][1] = 0;
	if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1)) && (conf[1][warpId] < lastval)) {
		
#pragma unroll 4		
		for(y = 0;y < specsplittings ;y +=2) {				
			if(y+initsplit < splittings) {
				
				int tmp = __popc(leafsplit[0]);
				if((NAGENTS-card) <= tmp && (NAGENTS-card+tmp) <= (card)) {
					rvalue[1][1] = rvalue[0][1] = f[leafsplit[0]];
						
					rvalue[0][1] += f[setdiff(conf[0][warpId],leafsplit[0])];						
					rvalue[1][1] += f[setdiff(conf[1][warpId],leafsplit[1])];					
				}
				leafsplit[1] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[1]);
				leafsplit[0] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[0]);
				if(rvalue[0][0] < rvalue[0][1]) {
					rvalue[0][0] = rvalue[0][1];
				}
				if(rvalue[1][0] < rvalue[1][1]) {
					rvalue[1][0] = rvalue[1][1];
				}

			}
			if(y+1+initsplit < splittings) {
//				rvalue[0][1] = rvalue[1][1] = 0;
				int tmp = __popc(leafsplit[0]);
				if((NAGENTS-card) <= tmp && (NAGENTS-card+tmp) <= (card)) {
					rvalue[1][1] = rvalue[0][1] = f[leafsplit[0]];
						
					rvalue[0][1] += f[setdiff(conf[0][warpId],leafsplit[0])];						
					rvalue[1][1] += f[setdiff(conf[1][warpId],leafsplit[1])];					
				}
				leafsplit[1] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[1]);
				leafsplit[0] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[0]);
				if(rvalue[0][0] < rvalue[0][1]) {
					rvalue[0][0] = rvalue[0][1];
				}
				if(rvalue[1][0] < rvalue[1][1]) {
					rvalue[1][0] = rvalue[1][1];
				}

			}

		}

	} else {		
		//	rvalue[0][0] =rvalue[0][1] = rvalue[1][1] = 0;
		if(conf[1][warpId] < lastval) {
#pragma unroll 4		
		for(y = 0;y < specsplittings ;y+=2) {
			if(y+initsplit < splittings) {
				
				int tmp = __popc(leafsplit[0]);
				if((NAGENTS-card) <= tmp && (NAGENTS-card+tmp) <= (card)) {

					rvalue[0][1] = f[leafsplit[0]];
					rvalue[1][1] = f[leafsplit[1]];								
					rvalue[0][1] += f[setdiff(conf[0][warpId],leafsplit[0])];

					rvalue[1][1] += f[setdiff(conf[1][warpId],leafsplit[1])];
		
				}
				leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);
				leafsplit[1] = SUBSET(conf[1][warpId],leafsplit[1]);
				if(rvalue[0][0] < rvalue[0][1]) {
					rvalue[0][0] = rvalue[0][1];
				}
				if(rvalue[1][0] < rvalue[1][1]) {
					rvalue[1][0] = rvalue[1][1];
				}
	
			}
			if(y+1+initsplit < splittings) {
				//rvalue[0][1] = rvalue[1][1] = 0;
				int tmp = __popc(leafsplit[0]);
				if((NAGENTS-card) <= tmp && (NAGENTS-card+tmp) <= (card)) {

					rvalue[0][1] = f[leafsplit[0]];
					rvalue[1][1] = f[leafsplit[1]];								
					rvalue[0][1] += f[setdiff(conf[0][warpId],leafsplit[0])];
					rvalue[1][1] += f[setdiff(conf[1][warpId],leafsplit[1])];					
				}
				leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);
				leafsplit[1] = SUBSET(conf[1][warpId],leafsplit[1]);
				if(rvalue[0][0] < rvalue[0][1]) {
					rvalue[0][0] = rvalue[0][1];
				}
				if(rvalue[1][0] < rvalue[1][1]) {
					rvalue[1][0] = rvalue[1][1];
				}
	
			}
		}
		} else {
#pragma unroll 4
		for(y = 0;y < specsplittings ;y+=2) {
			if(y+initsplit < splittings) {
				
				int tmp = __popc(leafsplit[0]);
				if((NAGENTS-card) <= tmp && (NAGENTS-card+tmp) <= (card)) {

					rvalue[0][1] = f[leafsplit[0]];
						
					rvalue[0][1] += f[setdiff(conf[0][warpId],leafsplit[0])];

				}
				leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);

				if(rvalue[0][0] < rvalue[0][1]) {
					rvalue[0][0] = rvalue[0][1];
				}
			}
			if(y+1+initsplit < splittings) {
				//rvalue[0][1] = rvalue[1][1] = 0;
				int tmp = __popc(leafsplit[0]);
				if((NAGENTS-card) <= tmp && (NAGENTS-card+tmp) <= (card)) {

					rvalue[0][1] = f[leafsplit[0]];
						
					rvalue[0][1] += f[setdiff(conf[0][warpId],leafsplit[0])];

				}
				leafsplit[0] = SUBSET(conf[0][warpId],leafsplit[0]);

				if(rvalue[0][0] < rvalue[0][1]) {
					rvalue[0][0] = rvalue[0][1];
				}
			}
		}
	}

	}

	#// pragma unroll	       
	// 		for(y = 1;y < parasplittings;y++) {
	// 			if(x+y+initsplit >= splittings) {
	// 				continue;
	// 			}

	// 		}


	 // if(__any( ( (rvalue[0][0] > value[0][warpId]) || (rvalue[1][0] > value[1][warpId]) ) ) == 0) {
	 // 	 return;
	 // }
		
	if(__any( ( rvalue[0][0] > value[0][warpId] ) ) ) {
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
		
	if(__any( ( rvalue[1][0] > value[1][warpId] ) ) ) {
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
		
	if(laneId == 0) {
		if(rvalue[0][0] > value[0][warpId]) {
			value[0][warpId] = rvalue[0][0];
		}
		if(rvalue[1][0] > value[1][warpId]) {
			value[1][warpId] = rvalue[1][0];
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



#define COMBS(X) ((1 << (X-1)) - 1)

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
	for(i = 2; i <= NAGENTS; i++) {
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
		splittings = ((1 << (i-1))-1);// COMBS(c1);///NPERBLOCK;
		threads = ((double) splittings)/ NPERBLOCK;
		threads = ceil(threads);
		const int o = 1;
		const int a[3] = {0,1,(2 << 2)};
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
			if((c2 > MAXVAL)) {
						      
		
					subsetcomp33 < MAXBLOCKSIZE , 1> <<<blocks,MAXBLOCKSIZE,0,stream[streamcount]>>>(dev_bids,splittings,MAXVAL,i,ca[0],ca[1],ca[2],ca[3]);
		
			}else{
		
					subsetcomp33 < MAXBLOCKSIZE , 0> <<<blocks,MAXBLOCKSIZE,0,stream[streamcount]>>>(dev_bids,splittings,MAXVAL,i,ca[0],ca[1],ca[2],ca[3]);
			
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
		printf("ended card %d blocks\t %d threads/block %u, n kernels %u \t time %lu \t splittings %d time per kernel %lf\n",i,blocks,bsize,count-prev,t,splittings,(double)t/(count-prev));
		prev =	count;

	}
	for (int i = 0; i < streams; ++i)
		cudaStreamDestroy(stream[i]);

	HANDLE_ERROR(cudaDeviceSynchronize());

	rend=clock();
	rt=(rend-rstart)/(CLOCKS_PER_SEC/1000);
	printf("real time %lu ms\n",rt);

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


int main(void) {
	/*Start n amount of assets*/
	dint from = NAGENTS;
	const unsigned long  MAXVAL = (2 << (from-1));
	time_t start,end,t;
//	O = (dint * ) malloc(sizeof(dint)*(2 << (from-1)));
	bids = (unsigned short * ) malloc(sizeof(short)*(2 << (from-1)));

	f =  (unsigned short * ) malloc(sizeof(short)*(2 << (from-1)));
	
//	return;
//	f = bids;
	int ret_val =1;
	int count;
	while(ret_val == 1) {
		
		//MAXVAL = (2 << (from-1));
		gen_rand_bids(MAXVAL);
		printf("hello\n");
		set_singleton_bid(MAXVAL);

		printf("maxval %u from %u\n",MAXVAL,from);
		start=clock();//predefined  function in c
		count = run_test(MAXVAL,from);
		end=clock();
		t=(end-start)/CLOCKS_PER_SEC;
		ret_val= recur_parse_wopt(MAXVAL);// parse_wopt(MAXVAL);
		printf("\nTime taken =%lu for n= %u with count %d average per count %lf\n", (unsigned long) t,from,count,(double)t/count);
	}
	free(f);

	return 0;
}
