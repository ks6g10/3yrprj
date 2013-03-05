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
//legacy code
#define TRUE 1
#define FALSE 0
/*Test sets all bids to one, which should give you n=|ITEMS| bids on output*/
#define TEST 1
/*Defines from 0-Range the random will give out*/
//legacy code
#define RANGE 10000
#define ITEMS 25

#define MAX (2 << (ITEMS-1))

//legacy code
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


//legacy var
uint32_t  * O;

unsigned short  * f, * bids;

struct _stack {
	dint conf;
	struct _stack * next;
} typedef stack;


inline  dint cardinality( dint seta) {
	return __builtin_popcount(seta);
}
int indexa =0;

//generate the bids, indexa is the controll bid, should always be present
void gen_rand_bids(dint MAXVAL) {
	register dint i = 0;
#if TEST
	unsigned int seed = (unsigned)time ( NULL );
	srand(seed);
	indexa++;
	for(i = 1; i < MAXVAL;i++) {
		bids[i] = rand()%10+1;
		O[i] = i;
	}
	bids[indexa] = 100;
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
//legacy code
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

//legacy code
/*Sets all bids with one element in it, |n| = 1*/
inline void set_singleton_bid(dint MAXVAL) {
	register  dint i;
	for(i =1; i< MAXVAL; i*=2) {
		f[i] = bids[i];
		if(bids[i] > 0)
			O[i] = i;
	}
}


//legacy code
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
#define MAXBLOCKSIZE (256)
//32
#define WARPSIZE (32)
//0
#define MIN_BLOCKS_PER_MP 8
#define NAGENTS (24)
//32 
#define NSTREAMS (32) 
//2
#define NPERBLOCK (2)
//32
#define CONFPKERNEL (32)
//4
#define NPARALLELCONF (4) 
#define TIMING (0)

#define COMP(Z) {							\
		if(shared_value[tid][Z] < shared_value[tid+i][Z]) {	\
			shared_value[tid][Z] = shared_value[tid+i][Z];	\
		}							\
	}

//timing function
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

template<int blockSize,int nparallelconf,int confpkernel,int nperblock,int currblocksize>
__global__ void
__launch_bounds__(currblocksize,MIN_BLOCKS_PER_MP)  
	subsetcomp33(
		/*0*/	unsigned short * __restrict__ f, /*value*/	       
		/*2*/	unsigned int * __restrict__ lock, /*pointer to the atomic lock array*/
		/*5*/	unsigned int maxval, /*the maximum number of subsets possible*/
		/*6*/	unsigned short count1, /*which index in the lock array*/
		unsigned int conf1, /*the intitial coalition structure*/
		unsigned int lastval//the value a coalition permutation can not exceed.
		)
{  
	//shared memory
	__shared__  unsigned short shared_value[(currblocksize >> 5)+1][nparallelconf]; // the values for each warp
	__shared__  unsigned int conf[confpkernel];// the configurations needed for the whole execution
	__shared__  unsigned char shift [confpkernel][NAGENTS];// the shift matrix/array
//registers
 	register unsigned int subset_value[nperblock][nparallelconf];//the value for one of the subset sums
	register unsigned int subset_conf[2];// temporary array to hold the subset coalition
	register unsigned int const tid = threadIdx.x;
	register unsigned int const ispec = nperblock*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;
	register int x,i,z; //counter

#if (TIMING == 1)
	clock_t stop_time, total,start_time;
	total = 0;
	start_time = clock(); 
#endif

	if(tid == 0) {
		subset_conf[0] = conf[0] = conf1;
		/* Generates the coalition structures using the bit tricks
		 * google "Compute the lexicographically next bit permutation"
		 *
		 */
#pragma unroll	
		for(x = 1;x < confpkernel; x++) { 			
			z = subset_conf[0] | (subset_conf[0]-1);
			subset_conf[0] = (z + 1) | (((~z & -~z) - 1) >> (__ffs(subset_conf[0])));
			conf[x]= subset_conf[0];
		}
		CHECKPOINT("l1\t %d\n"); // first checkpoint, see paper
	}
//	__syncthreads();

#if (TIMING == 1)
	start_time = clock(); 
#endif
	/* Generates the shift array in order to in the future to create subsets.
	 * as confpkernel is <= 32, this will be done in one warp, hence no need
	 * for the syncthreads above.
	 */
	
	if(tid < confpkernel) { //generate the shift arrays		
		subset_conf[0] = conf[tid]; // re-use registers tmpval
		unsigned int cnt = 0;
		while(subset_conf[0]) {
			subset_conf[1] = __ffs(subset_conf[0]) - 1; //find which index is first bit
			subset_conf[0] &= ~(1 << subset_conf[1]);//set nth bit to 0
			shift[tid][cnt++] = subset_conf[1];
		}
	}

	CHECKPOINT("l2\t %d\n");

	__syncthreads(); //all other warps except the first one will wait here for the first one to finish

#pragma unroll	
	for(x =0; x < confpkernel; x += nparallelconf) { // fetches the values in batches

#if (TIMING == 1)
		start_time = clock(); 
#endif
		if(conf[x] > lastval) {//if the permutation is larger than the full set, c> (1 << NAGENTS)
			continue;
		}
#pragma unroll	
 		for(z=0; z < nparallelconf;z++) {//reset values
 			subset_value[0][z] = subset_value[1][z] = 0U;			
 		}

		CHECKPOINT("l3\t %d\n");

		if(ispec >= maxval) { // no need to fetch, have greater splitting index than splittings possible
			goto postfetch;
		}

		CHECKPOINT("l3\t %d\n");

		/* Hold onto your hat, this will be a wild ride
		 * Here is where the collision detection happens
		 * Which collitions do I try to find?
		 * The ones that causes collisions and those that make the runtime faster
		 * ergo not all that causes a collision
		 * How did I find which causes collision?
		 * if(coalision_x == coalision_y) print(hello world)
		 * I evaluate two coallisions at the same time, 
		 * will reference them as second and first coalition
		 * and subset means one half of a splitting, and setdiff means the other half
		 */

		unsigned int tmp[2];
		tmp[1]= tmp[0] = 0; //reset
		subset_conf[1] = 0; //reset
#pragma unroll
		for(z=0;z < nparallelconf;z +=2) {

			if(conf[z+x] >= lastval) {
				continue;
			}		 
			unsigned int tmpstore = 0, tmpstore1 = 0;
			unsigned int index;
			//unsigned int const tcar = __popc(tmp);
			unsigned int wshile = ispec;
			subset_conf[0] = 0;
			
			/* Generate the subset (intitialsplit) and fetch the value for the second coalition structure first
			 * Why? It is so much faster than generating and fetching the first
			 * Why? Do not know, but my guess is that it is better to block early one time than to block
			 * twice even though it should be the same time, I guess better scheduling, did not profile that part.
			 * Escept for speed, gained like 5-6 seconds @ 24 agents
			 */
			if(conf[z+x+1] < lastval) {
				wshile = ispec;
				subset_conf[1] = 0;
				while(wshile) {
					index = __ffs(wshile)-1;
					subset_conf[1] += (1 << shift[x+z+1][index]);
					wshile &= ~(1 << index);	
				}
				/*As I start from index 0, need to use makro SUBSET( nextSplit) once before you get the right value*/
				subset_conf[1] = SUBSET(conf[z+x+1],subset_conf[1]); 
				/*Fetch the value*/
				subset_value[0][z+1] = f[subset_conf[1]];
			}
			while(wshile) {
				index = __ffs(wshile)-1;
				subset_conf[0] += (1 << shift[x+z][index]);
				wshile &= ~(1 << index);
				if(wshile) {
					index = __ffs(wshile)-1;
					subset_conf[0] += (1 << shift[x+z][index]);
					wshile &= ~(1 << index);
				}
			}
			subset_conf[0] = SUBSET(conf[z+x],subset_conf[0]);

			/*If there is enough splittings, we can get the next splitting as well
			* just use SUBSET (nextsplit) to generate the next splitting
			* this is for the first coalition structure
			**/
			if((ispec+1) < maxval) {
				tmp[0] = SUBSET(conf[z+x],subset_conf[0]);
				/*The first collision, I have pre-emptively stored the previouse 
				 * subset of the second splitting of the second coalision structure's
				 * value in the register for this second splitting of the first coalision structure
				 * if they are not the same you shal fetch me a new value
				 */
				if(tmp[0] != tmp[1]) {
					subset_value[1][z] = f[tmp[0]];
				}
			}
			/* Fetch the first coalitions first splittings value
			 */
			subset_value[0][z] = f[subset_conf[0]];		
			
			/*The second splitting of the second coalition
			 */
			if((ispec+1) < maxval && conf[z+x+1] < lastval) {
				tmp[1] = SUBSET(conf[z+x+1],subset_conf[1]);
				//Collision
				if(tmp[1] == tmp[0]) {
					subset_value[1][z+1] = subset_value[1][z];
				} else {
					subset_value[1][z+1] = f[tmp[1]];
				}
			}

			subset_conf[0] = (setdiff(conf[z+x],subset_conf[0]));
			tmpstore = f[subset_conf[0]];//

			//the setdiff of the second splitting @ first coalition
			if((ispec+1) < maxval) {
				tmp[0] = (setdiff(conf[z+x],tmp[0]));
				tmpstore1 = f[tmp[0]];
			}
			
			//second coalision
			if(conf[z+x+1] < lastval) {
				subset_conf[1] = (setdiff(conf[z+x+1],subset_conf[1]));
				if(z < (nparallelconf-2)){
				 	subset_value[1][z+2] = subset_value[1][z+1]; // preemptive storage for first Collision
				}
				//Collision
				if(subset_conf[1] == subset_conf[0]){// dont remove
					subset_value[0][z+1] += tmpstore;					
				} else{
					subset_value[0][z+1] += f[subset_conf[1]];
				}
				//next splitting
				if((ispec+1) < maxval) {
					tmp[1] = (setdiff(conf[z+x+1],tmp[1]));
					subset_value[1][z+1] += f[tmp[1]];										
				}
			}
			subset_value[0][z] += tmpstore;
			subset_value[1][z] += tmpstore1;
			tmp[1] = (setdiff(conf[z+x+1],tmp[1])); //get the subset instead of the setdiff in order to do the first Collision

		}
		CHECKPOINT("l4\t %d\n");
					
	postfetch:

#pragma unroll	
		for(z = 0; z < nparallelconf;z++) {
			if(conf[z+x] >= lastval) {
				continue;
			}
			if(subset_value[1][z] > subset_value[0][z]) { // see which one is greater
				subset_value[0][z] = subset_value[1][z];
			}
#pragma unroll	
			for(i = 16;i >=1;i >>=1) {//warp reduction
				int warp_value = __shfl_xor((int)subset_value[0][z],i,32); // exchange values between threads
				if(warp_value > subset_value[0][z]) {
					subset_value[0][z] =(unsigned int) warp_value;
				}
			}
			//tid&(WARPSIZE-1) == tid%WARPSIZE
			//Only threads with line id == 0 is allowed to update in the shared memory,
			//i.e. the first thread in each warp
			if(!(tid&(31))) {
				unsigned int index = tid >> 5; // tid >> 5 == tid / 32 which warp it is
				shared_value[index][z] = subset_value[0][z];
			}

		}

		CHECKPOINT("l5\t %d\n");		

		//template optimisation
		if((currblocksize/32) > 1) {
 		__syncthreads();
		}
		//how many warps is it, block dimension divided by warp size
		//e.g. 256/32 == 256 >> 5
//		if((currblocksize/32) > 1) {//evaluated by the pre-processor
			
		i = blockDim.x >> 6;
		if(tid<i) {
			for(; i > 0; i >>= 1) {//reduction in shared memory
				for(z=0; z < nparallelconf;z +=4) {
					COMP(z); // small macro, lazy to put it back
				}
			}
		}
		
		CHECKPOINT("l6\t %d\n");
		/*
		 *Update the value in f using the atomic lock array
		 */
		if((currblocksize) >= blockSize) {//evaluated by the pre-processor
			//__syncthreads(); // no need as the shared memory reduction is in the same warp
			if(tid == 0) {
#pragma unroll	
				for(z=0; z < nparallelconf;z++) {
					if(conf[z+x] > lastval) {continue;}
					if(f[conf[z+x]] < shared_value[0][z]) {
						if(atomicMax(&(lock[count1+x+z]),shared_value[0][z]) < shared_value[0][z]) {
							if(f[conf[z+x]] < shared_value[0][z]) {
								f[conf[z+x]] = shared_value[0][z];
							}
						}
					}
				}
			}
		} else {
			/*if only one block*/
#pragma unroll
			for(z=0; z < nparallelconf;z++) {
				if(conf[z+x] > lastval) {continue;}
				if(f[conf[z+x]] < shared_value[0][z]) {
					f[conf[z+x]] = shared_value[0][z];
						
				}

			}
		}
		CHECKPOINT("l7\t %d\n");
	}
}


#define COMBS(X) ((1 << cardinality(X)-1) - 1)

int run_test(dint MAXVAL,dint items) {
/*Setup the environment*/
	//dint perm[MAXVAL];

	register unsigned int i,c1,count =0;
	unsigned short *dev_bids;

	count = 0;
 
	HANDLE_ERROR(cudaDeviceReset());
	HANDLE_ERROR(cudaSetDeviceFlags(cudaDeviceScheduleYield)); 
  	HANDLE_ERROR(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitMallocHeapSize,0));
	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitStackSize,0));
	HANDLE_ERROR(cudaDeviceSetLimit(cudaLimitMallocHeapSize,0));
	HANDLE_ERROR(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
	unsigned int * dev_lock1,*dev_lock2,*dev_ptr;
	const	unsigned int devcount = 1024*CONFPKERNEL;// count;
	register unsigned int streams = NSTREAMS;
	register unsigned short lock_count = 0;
	register unsigned int streamcount = 0;
	register cudaStream_t stream[streams];
	for(i = 0;i < streams; i++)
		HANDLE_ERROR(cudaStreamCreate(&stream[i]));

	count = 0;
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock1,(devcount)*sizeof(int)));
	HANDLE_ERROR(cudaMalloc((void **)&dev_lock2,(devcount)*sizeof(int)));

 	HANDLE_ERROR(cudaMalloc((void **)&dev_bids, MAXVAL*sizeof(short)));

 	HANDLE_ERROR(cudaMemcpy(dev_bids,bids,MAXVAL*sizeof(short),cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemset(dev_lock1,0,devcount*sizeof(int)));
	HANDLE_ERROR(cudaMemset(dev_lock2,0,devcount*sizeof(int)));
	/*2.*/
	dev_ptr = dev_lock1;
	register unsigned int bsize = 32;
	register int blocks;
	int prev =0;
	lock_count = 0;
	time_t rstart,rend,rt;
	rstart=clock();
	// start with small coalition structures and go big
	for(i = 2; i <= items; i++) {
		time_t start,end,t;

		start=clock();
		int splittings;

		double threads;
		c1 = (1 << i) -1;
		splittings =  COMBS(c1);///NPERBLOCK;
		threads = ((double) splittings)/ NPERBLOCK;
		threads = ceil(threads);
		
		for(; c1 <= MAXVAL;) {
			while( bsize < MAXBLOCKSIZE && threads > bsize ) {//dynamicly adjust the blocksize, 32 -> 256
				bsize += 32;
			}
			blocks =(int)  ceil((threads/bsize)); // number of blocks
			

			switch(bsize) {
			case 32:
				subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK,32> <<<blocks,32,0,stream[streamcount]>>>(dev_bids,dev_ptr,splittings,lock_count,c1,MAXVAL);
				break;
			case 64:
				subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK,64> <<<blocks,64,0,stream[streamcount]>>>(dev_bids,dev_ptr,splittings,lock_count,c1,MAXVAL);
				break;
			case 128:
				subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK,128> <<<blocks,128,0,stream[streamcount]>>>(dev_bids,dev_ptr,splittings,lock_count,c1,MAXVAL);
				break;
			case 256:
				subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK,256> <<<blocks,256,0,stream[streamcount]>>>(dev_bids,dev_ptr,splittings,lock_count,c1,MAXVAL);
				break;
			// case 512:
			// 	subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK,512> <<<blocks,512,0,stream[streamcount]>>>(dev_bids,dev_ptr,splittings,lock_count,c1,MAXVAL,i);
			// 	break;
			//  case 1024:
			//  	subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK,1024> <<<blocks,1024,0,stream[streamcount]>>>(dev_bids,dev_ptr,splittings,lock_count,c1,MAXVAL,i);
			//  	break;	
				
			}
			//subsetcomp33 < MAXBLOCKSIZE , NPARALLELCONF , CONFPKERNEL , NPERBLOCK > <<<blocks,bsize,0,stream[streamcount]>>>(dev_bids,dev_o,dev_ptr,splittings,lock_count,c1,MAXVAL,i);

			for(int x = 0; x < CONFPKERNEL;x++) {
				t = c1 | (c1-1);
				c1 = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c1) + 1));
			}
			streamcount++;
			count++;
  			
			lock_count += CONFPKERNEL;
			
			/*
			 *Handle the streams and locks
			 *have two lock arrays I switch between
			 * once one is full, switch and issue an memset on the full
			 *
			 */
			
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
	rt=(rend-rstart)/(CLOCKS_PER_SEC/1000);
	printf("real time %lu\n",rt);

	HANDLE_ERROR(cudaMemcpy(f,dev_bids,MAXVAL*sizeof(short),cudaMemcpyDeviceToHost));
	HANDLE_ERROR(cudaFree(dev_bids));
	HANDLE_ERROR(cudaFree(dev_lock1));
	HANDLE_ERROR(cudaFree(dev_lock2));

	HANDLE_ERROR(cudaDeviceReset());
	return count;
}
//max function in DP in order to get the final splittings just like in idp of removing the splittings table of a given C
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
//this one handels the the retrival of the splittings like IDP
int recur_parse_wopt(dint MAXVAL) {
	stack * root = (stack *) malloc(sizeof(stack));
	stack * sroot = NULL;
	stack * scurr = NULL;
	root->conf = (MAXVAL)-1; //just like CS = {A} in DP
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
//legacy code of the above function with the final splittings in an array like DP
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
	while(ret_val == 0) { // remove this to only test one valye, not all

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

//legacy code
template<int blockSize,int nparallelconf,int confpkernel,int nperblock,int currblocksize>
__global__ void
__launch_bounds__(currblocksize,MIN_BLOCKS_PER_MP)  
	subsetcomp32(
		/*0*/	unsigned short * __restrict__ f, /*Bid value*/
		/*1*/	unsigned int * __restrict__ O, /*The move array*/
		/*2*/	unsigned int * __restrict__ lock,
		/*5*/	unsigned int maxval,
		/*6*/	unsigned short count1,
		unsigned int conf1,
		unsigned int lastval,//the value a permutation can not exceed.
		unsigned int card)
{  
	//confpkernel = how many configurations the kernel will evaluate
	//nparallelconf = how many configurations the kernel will evaluate at the same time
	__shared__  unsigned short shared_value[(currblocksize >> 5)+1][nparallelconf];
	__shared__  unsigned int shared_conf[(currblocksize >> 5)+1][nparallelconf];
	__shared__  unsigned int conf[confpkernel];// the configurations needed for the whole execution
	__shared__  unsigned short shift [confpkernel][NAGENTS];// the shift matrix/array
	//__shared__ volatile unsigned int tmp [confpkernel][2];
 	register unsigned int subset_value[2][8];//the value for one of the subset sums
	register unsigned int subset_conf[2][8];
   
	register unsigned int count = count1;
	register unsigned int const tid = threadIdx.x;
//	__shared__ unsigned int
	register unsigned int ispec = nperblock*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;
	register int x,i,z; //counter

#if (TIMING == 1)
	clock_t stop_time, total,start_time;
	total = 0;
	start_time = clock(); 
#endif

	if(tid == 0) {
		conf[0] = conf1;
#pragma unroll	
		for(x = 1;x < confpkernel; x++) { // generate the configurations
			conf[x] = conf[x-1];
			z = conf[x] | (conf[x]-1);
			conf[x] = (z + 1) | (((~z & -~z) - 1) >> (__ffs(conf[x])));
		}
		CHECKPOINT("l1 %d\n");
	}
	if(confpkernel > 32) {
		__syncthreads();
	}
	
	count = count1;
	if(tid < confpkernel) { //generate the shift arrays
		subset_conf[0][0] = conf[tid]; // re-use registers tmpval
		subset_conf[0][1] = 0;// re-use registers index
#pragma unroll		
		for(x = 0,i=0; x < card;x++) {//could put card in template to unroll
			subset_conf[0][1] = __ffs(subset_conf[0][0]) - 1; //find which index is first bit
			subset_conf[0][0] &= ~(1 << subset_conf[0][1]);//set nth bit to 0
			shift[tid][x] = subset_conf[0][1];
		}
	}

	CHECKPOINT("l2 %d\n");

	__syncthreads();

#pragma unroll	
	for(x =0; x < confpkernel; x += nparallelconf) {

		if(conf[x] > lastval) {//if the permutation is larger than the full set, c> (1 << NAGENTS)
			continue;
		}
//#pragma unroll	
//		for(i = 0; i < nperblock; i++) {
#pragma unroll	
 		for(z=0; z < nparallelconf;z++) {
 			subset_value[0][z] = subset_value[1][z] = 0U;			
 		}
//		}
		CHECKPOINT("l3 %d\n");
		if(ispec >= maxval) {
			goto postfetch;
		}
			//This for loop initilize the first subset configuration.

#pragma unroll	
		for(i = 0; i < nparallelconf; i++) { 
			unsigned int tmp = ispec;
			//unsigned int const tcar = __popc(tmp);
			subset_conf[0][i] = 0;
			while(tmp) {
				unsigned short index = __ffs(tmp)-1;
				subset_conf[0][i] += (1 << shift[x+i][index]);//CHECK
				tmp &= ~(1 << index);				
			}
		}
		CHECKPOINT("l4 %d\n");	
#pragma unroll
		for(z=0;z < nparallelconf;z++) {
			if(conf[z+x] > lastval) {
				continue;
			}
			subset_conf[0][z] = SUBSET(conf[z+x],subset_conf[0][z]);
			subset_value[0][z] = f[(setdiff(conf[z+x],subset_conf[0][z]))] + f[subset_conf[0][z]];
			//ispec++;
			if((ispec+1) >= maxval) {
				continue;
			}
			subset_conf[1][z] = SUBSET(conf[z+x],subset_conf[0][z]);
			subset_value[1][z] = f[(setdiff(conf[z+x],subset_conf[1][z]))] + f[subset_conf[1][z]];	
		}
		CHECKPOINT("l6 %d\n");
					
	postfetch:

#pragma unroll	
		for(z = 0; z < nparallelconf;z++) {//warp reduction
			if(subset_value[1][z] > subset_value[0][z]) {
				subset_value[0][z] = subset_value[1][z];
				subset_conf[0][z] = subset_conf[1][z];
			}
#pragma unroll	
			for(i = 16;i >=1;i >>=1) {
				int warp_value = __shfl_xor((int)subset_value[0][z],i,32);
				int warp_conf = __shfl_xor((int)subset_conf[0][z],i,32);
				if(warp_value > subset_value[0][z]) {
					subset_value[0][z] =(unsigned int) warp_value;
					subset_conf[0][z] =(unsigned int) warp_conf;
				}
			}
			//tid&(WARPSIZE-1) == tid%WARPSIZE
			//Only threads with line id == 0 is allowed to update in the shared memory,
			//i.e. the first thread in each warp
			if(!(tid&(31))) {
				unsigned int index = tid >> 5; // tid >> 5 == tid / 32 which warp it is
				shared_value[index][z] = subset_value[0][z];
				shared_conf[index][z] = subset_conf[0][z];
			}

		}
		CHECKPOINT("l7 %d\n");		

		//	CHECKPOINT("l8 %d\n");
		if((currblocksize/32) > 1) {
 		__syncthreads();
		}
		//how many warps is it, block dimension divided by warp size
		//e.g. 256/32 == 256 >> 5
		if((currblocksize/32) > 1) {//evaluated by the pre-processor
			i = (currblocksize >> 6);//blockDim.x >> 6;
			if(tid<i) {//reduction mby move down if you get wrong results gained ~1000 cycles

#pragma unroll	
			for(; i > 0; i >>= 1) {
#pragma unroll	
  					for(z=0; z < nparallelconf;z +=4) {
						COMP(z);
  					}
			}
			}
		}
		
		CHECKPOINT("l9 %d\n");
		if((currblocksize) >= blockSize) {//evaluated by the pre-processor
			//__syncthreads();
			if(tid == 0) {
#pragma unroll	
				for(z=0; z < nparallelconf;z++) {
					if(f[conf[z+x]] < shared_value[0][z]) {
						//	printf("lock val %u shared_val %s\n",lock[count+z] ,shared_value[0][z]);
						if(atomicMax(&(lock[count+z]),shared_value[0][z]) < shared_value[0][z]) {
							//	O[conf[z+x]] = shared_conf[0][z];
							f[conf[z+x]] = shared_value[0][z];
						}
					}
				}
			}
		} else {
#pragma unroll
			for(z=0; z < nparallelconf;z++) {
				if(f[conf[z+x]] < shared_value[0][z]) {
					//	O[conf[z+x]] = shared_conf[0][z];
					f[conf[z+x]] = shared_value[0][z];
						
				}

			}
		}

		count += nparallelconf;
		ispec = NPERBLOCK*(threadIdx.x + blockDim.x * blockIdx.x);//I + offset;

	}
}
