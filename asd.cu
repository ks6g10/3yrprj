#include <stdio.h>
#include <string.h> // for ffs
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_profiler_api.h>


//Not removing useless stuff that may have traces left in the code as I am on my laptop which can not compile cuda code.

//LEGACY START - Could possibly be removed, do not listen to the comments here, it does not work
/*Debug enabled gives more print statements of bids and how the "Matrix" gets evaluated*/
#define DEBUG 0
#define TRUE 1
#define FALSE 0
/*Test sets all bids to one, which should give you n=|ITEMS| bids on output*/
#define TEST 1
/*Defines from 0-Range the random will give out*/
#define RANGE 10000
#define ITEMS 25


#if ITEMS < 8
#define dint uint8_t
#elif ITEMS < 16
#undef dint
#define dint uint16_t
#elif ITEMS < 32
#undef dint
#define dint uint32_t
#endif
//LEGACY END


//Handle errors for cuda
static void HandleError( cudaError_t err, const char *file,int line ) {
	if (err != cudaSuccess) {
		printf( "%s in %s at line %d\n", cudaGetErrorString( err ),file, line );
		exit( EXIT_FAILURE );
	}
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))

//useless
uint32_t  * O;

//the bids and f array, why not make them global eh.
unsigned short  * f, * bids;

//Small stack implementation
struct _stack {
	dint conf;
	struct _stack * next;
} typedef stack;


inline  dint cardinality( dint seta) {
	return __builtin_popcount(seta);
}
int indexa =0; //the joker coalition
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
	bids[indexa] = 100; //set the joker to something high
	printf("index %d \n",indexa);
	if(indexa >= MAXVAL) {
		printf("No error\n");
		exit(0);
	}

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
// //You can touch this
#define MAXBLOCKSIZE (1024)
//32
#define WARPSIZE (32)
//0 //You can touch this
#define MIN_BLOCKS_PER_MP 4


//You can touch this
#define NAGENTS (25)

//32  //You can touch this
#define NSTREAMS (32)

//2 //You should probably not touch this
#define NPERBLOCK (2)

//DO NOT TOUCH THIS, IT IS NOT DYNAMIC
#define confpwarp (2)

//Cant touch this na na na na
#define CONFPKERNEL ((MAXBLOCKSIZE/32)*confpwarp)

//You can touch this
#define parasplittings (8)

//not used any more
#define NPARALLELCONF (4)

//have not inserted any timings in the code
#define TIMING (0)

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
     //coalition structures  indexed [0-1][warpId]
	__shared__ unsigned int conf[confpwarp][(blockSize/32)+1];
	//values
	__shared__ unsigned short value[confpwarp][(blockSize/32)+1];
	//shift array for init split
	__shared__  unsigned char shift [(blockSize/32)*confpwarp][NAGENTS];// the shift matrix/array
	const unsigned int tid = threadIdx.x;
	const unsigned int laneId = (tid%32);
	const unsigned int warpId = tid/32;
	//How many splittings can we do per thread
	const unsigned int specsplittings = (!!(splittings%32))+(splittings/32);
	//the first splitting
	const unsigned int initsplit = tid*specsplittings;
	//Holds splittings
	unsigned int leafsplit[2];
	//values in registers
	unsigned int rvalue[confpwarp][parasplittings];
	/*Thread 0 of each warp*/
	if(!tid) {
		unsigned int tmp;
		unsigned int conftmp;
		//assign right coalition to right block
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
		//genereate the next coalition structures
		tmp = conftmp | (conftmp-1);
		conftmp = (tmp + 1) | (((~tmp & -~tmp) - 1) >> (__ffs(conftmp)));
		conf[1][warpId] = conftmp;
		//genereate the next coalition structures
		for(x =1; x < (blockSize/32);x++) {
			tmp = conftmp | (conftmp-1);
			conftmp = (tmp + 1) | (((~tmp & -~tmp) - 1) >> (__ffs(conftmp)));
			conf[0][x] = conftmp;

			tmp = conftmp | (conftmp-1);
			conftmp = (tmp + 1) | (((~tmp & -~tmp) - 1) >> (__ffs(conftmp)));
			conf[1][x] = conftmp;
		}
	}
	__syncthreads();
	//if coalition greater than what we can calculate, i.e. out of range, e.g. 24th bit set in a 23 sized problem
	if(conf[0][warpId] > lastval) {
		return;
	}

	//the collision detection, if |c & c'| == |c|-1 && |c| == |c'| then we can generate all splittings for both of c and c' using their intersection
	// here we set up the shift array
	if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1)) && conf[1][warpId] < lastval) {
	     //one thread
		if(laneId == 0) {
			unsigned int index;
			unsigned int count = 0;
			unsigned int conftmp = conf[0][warpId] & conf[1][warpId];


			//fetch the coalitions values to shared memory
			value[1][warpId] = f[conf[1][warpId]];

			value[0][warpId] = f[conf[0][warpId]];

			//look at the paper
			while(conftmp) {
				index = __ffs(conftmp) - 1; //find which index is first bit
				conftmp &= ~(1 << index);//set nth bit to 0
				shift[warpId][count] = index;
				shift[warpId+1][count] = index;
				count++;
			}
		}
		//two threads
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


	//Generate the initial splitting, if the second coalition is less than lastval do both else do the single one
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
		//as we base our index from zero, do one run with nextsplit
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

	//fetch splittings
	for(x = 0; x < specsplittings;x += parasplittings) {
	     //if there is a collision
		if((__popc(conf[1][warpId] & conf[0][warpId]) == (card - 1)) && (conf[1][warpId] < lastval)) {
#pragma unroll 8
		     for(y = 0;y < parasplittings ;y++) {
			  //reset value
				rvalue[0][y] = 	rvalue[1][y] = 0;
				if(x+y+initsplit < splittings) {
				     //IDP
					int tmp = __popc(leafsplit[0]);
					if((NAGENTS-card) <= tmp && (NAGENTS-card) <= (card-tmp)) {
					     //fetch the subset in common
						rvalue[1][y] = rvalue[0][y] = f[leafsplit[0]];
						//fetch the set difference
						rvalue[0][y] += f[setdiff(conf[0][warpId],leafsplit[0])];
						rvalue[1][y] += f[setdiff(conf[1][warpId],leafsplit[1])];
					}
					//nextsplit
					leafsplit[1] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[1]);
					leafsplit[0] = SUBSET((conf[1][warpId] & conf[0][warpId]),leafsplit[0]);
				}
			}

		} else {			
#pragma unroll 8			
			for(y = 0;y < parasplittings ;y++) {
				rvalue[0][y] = rvalue[1][y] = 0;
				if(x+y+initsplit < splittings) {
				     //idp
					int tmp = __popc(leafsplit[0]);
					if((NAGENTS-card) <= tmp && (NAGENTS-card) <= (card-tmp)) {
						rvalue[0][y] = f[leafsplit[0]];
						rvalue[0][y] += f[setdiff(conf[0][warpId],leafsplit[0])];
						//template optimisation in order to remove one if statement from kernel launches which does not exceed
						//the maximum coalition structure
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

		//register reduction
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

		//if their values are less than the value in shared memory, continue onto the next loop
		if(__ballot( ( (rvalue[0][0] > value[0][warpId]) || (rvalue[1][0] > value[1][warpId]) ) ) == 0) {
			continue;
		}

		//if rvalue is gfreater than value in shared memory, do warp reduction, this is a unrolled version
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
		//same but different value
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
		//update shared memory
		if(laneId == 0) {
			if(rvalue[0][0] > value[0][warpId]) {
				value[0][warpId] = rvalue[0][0];
			}
			if(rvalue[1][0] > value[1][warpId]) {
				value[1][warpId] = rvalue[1][0];
			}
		}
	}
	//lastly update global memory, no need for atomics as there is no contender
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

	HANDLE_ERROR(cudaFree(dev_bids));

	HANDLE_ERROR(cudaDeviceReset());

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
	dint MAXVAL = (2 << (from-1));

	time_t start,end,t;
//	O = (dint * ) malloc(sizeof(dint)*(2 << (from-1)));
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

