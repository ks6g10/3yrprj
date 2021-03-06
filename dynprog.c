#include <stdio.h>
#include <string.h> // for ffs
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/*Debug enabled gives more print statements of bids and how the "Matrix" gets evaluated*/
#define DEBUG 0
/*Test sets all bids to one, which should give you n=|ITEMS| bids on output*/
#define TEST 1
/*Defines from 0-Range the random will give out*/
#define RANGE 24
#define ITEMS 26
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

dint bids[MAX];
// dint bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
//dint bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
// dint bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6
dint wopt[MAX] = {0};
dint f[MAX] = {0};
dint O[MAX] = {0};
  
inline  dint setdiff( dint seta,  dint setb) {

     return (seta & ~setb);
}

inline  dint cardinality( dint seta) {
     return __builtin_popcount(seta);
}
int indexa =0;
int gen_rand_bids(dint MAXVAL) {
     register dint i = 0;
#if TEST
     for(i = 1; i < MAXVAL;i++) {
	  bids[i] = 1;
	  O[i] = i;
     }
     //unsigned int seed = (unsigned)time ( NULL );
     //srand(seed);
     //indexa = rand() % MAXVAL;
     indexa++;
     bids[indexa] = 100;
     //    printf("-------index %d \n",indexa);
     if(MAXVAL <= indexa) {
	     printf("clean run\n");
	     return 1;
     }
     return 0;

#else
     for(i = 1; i < MAX;i++) {
	  bids[i] = rand() % RANGE;
     }
#endif
}

/*Reminder of sets
 *f[]
 *O[]
 *bids[]	
 */


/*Sets all bids with one element in it, |n| = 1*/
inline void set_singleton_bid(dint MAXVAL) {
     register  dint i;
     for(i =1; i< MAXVAL; i*=2) {
	  f[i] = bids[i];
	  //if(bids[i] > 0)
	       O[i] = i;
     }
}


struct _stack {
	 dint conf;
	struct _stack * next;
} typedef stack;


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
			return;
		}
		/*if something is wrong*/
		if(count > 40) {
			fprintf(stderr,"Something went wrong at line %d in %s\n",__LINE__,__FILE__);
			return;
		}
		//	printf("curr %u\t",curr->conf);
		if(conf != O[conf]) {

			dint diff = setdiff(conf,O[conf]);
			curr->conf = O[conf];
			stack * tmp = (stack *) malloc(sizeof(stack));
			//	printf("diff %u\t O[conf] %u f %u\n",diff,O[conf],f[conf]);
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
//	printf("\n");
	curr = sroot;
	dint tmp = 0;
	unsigned int bool = 0;
	while(curr != NULL) {
		if(bids[curr->conf]) {
				if(curr->conf == indexa) {
					//printf("something is wrong index %d\n",indexa);
					
					return 0;
				}
				printf("conf %u value %u\n",curr->conf,bids[curr->conf]);
			tmp++;
		}
		stack * tmp = curr;
		curr = curr->next;
		free(tmp);
	}
	printf("n = %u\n",tmp);
	return 1;
}

void printfo(dint MAXVAL) {
	int i;
	printf("\n");
	for(i = 1; i< MAXVAL; i++)
	{
		printf("Bid[%d]\t%u\tF[%d]\t%u\tO[%d]\t%u\t\n",i,bids[i],i,f[i],i,O[i]);

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

/*n 15 t 9 n 16 t 42*/

uint64_t no,yes;

void max3(dint conf1,dint conf2,dint idp, unsigned int doidp) {
     register dint card = cardinality(conf)/2;
     register dint combinations = (1 << cardinality(conf)-1)-1;
     register dint max = bids[conf];
     register dint set = conf;
     register dint tmp = 0;
     register dint subset = 0;
     register const dint inverse = ~conf;
     register dint i;
     O[conf] = set;
     int do_collision = (cardinality(conf1 & conf2) ==( cardinality(conf1)-1));
     subset = ((inverse+1)+subset)&(conf1&conf2);
     for(i = 0;i<combinations; i++) {
	     setcard = cardinality(subset);

     }
     if(doidp) {
	     for(i = 1;i<=combinations; i++) {
			  subset = ((inverse+1)+subset)&conf;

			  register unsigned int setcard = cardinality(subset);
			  register unsigned int tmpset = setdiff(conf,subset);
			  register unsigned int tmpcard = cardinality(tmpset);
	  //  if(setcard > card) {
		  //  no++;
	  //	  continue;
	  //  }
			  if(setcard < tmpcard) {
				  setcard = tmpcard;
			  }
			  
			  if(setcard < idp) {
				  continue;
			  }
			  // yes++;	
	  
			  tmp = f[tmpset] + f[subset];
			  if(max < tmp) {
				  max = tmp;
				  set = subset;
			  }
		  }
     }
     else {
	     for(i = 1;i<=combinations; i++) {
		     subset = ((inverse+1)+subset)&conf;
		     tmp = f[setdiff(conf,subset)] + f[subset];
		     
		     if(max < tmp) {
			     max = tmp;
			     set = subset;
		     }
	  
	     }
     }
     f[conf] = max;
     O[conf] = set;
}

void max2(dint conf,dint idp, unsigned int doidp) {
     register dint card = cardinality(conf)/2;
     register dint combinations = (1 << cardinality(conf)-1)-1;
     register dint max = bids[conf];
     register dint set = conf;
     register dint tmp = 0;
     register dint subset = 0;
     register const dint inverse = ~conf;
     register dint i;
     O[conf] = set;
     if(doidp) {
	          for(i = 1;i<=combinations; i++) {
			  subset = ((inverse+1)+subset)&conf;

			  register unsigned int setcard = cardinality(subset);
			  register unsigned int tmpset = setdiff(conf,subset);
			  register unsigned int tmpcard = cardinality(tmpset);
	  //  if(setcard > card) {
		  //  no++;
	  //	  continue;
	  //  }
			  if(setcard < tmpcard) {
				  setcard = tmpcard;
			  }
			  
			  if(setcard < idp) {
				  continue;
			  }
			  // yes++;	
	  
			  tmp = f[tmpset] + f[subset];
			  if(max < tmp) {
				  max = tmp;
				  set = subset;
			  }
		  }
     }
     else {
	     for(i = 1;i<=combinations; i++) {
		     subset = ((inverse+1)+subset)&conf;
		     tmp = f[setdiff(conf,subset)] + f[subset];
		     
		     if(max < tmp) {
			     max = tmp;
			     set = subset;
		     }
	  
	     }
     }
     f[conf] = max;
     O[conf] = set;
}

void run_test(dint MAXVAL,dint items) {
//	     printfo(MAXVAL);
	    printf("before\n");
	    register dint i, c;
	    /*Select which cardinality we inspect*/
	    for(i = 2; i <= items; i++) {
		    /*Generate the set of cardinality i*/
		    for(c = (1 << i) -1; c <= MAXVAL;) {
			    max2(c,items-i,0);			  
			    //bit hacks "Compute the lexicographically next bit permutation"
			    register dint t = c | (c-1);
			    c = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(c) + 1)); 
			    //end ref
		    }
     }
//               printfo(MAXVAL);
      
}


int main(void) {
     /*Start n amount of assets*/
     dint from = 26;
     /*End amount of assets, inclusive*/
     dint till = 26;
     dint MAXVAL = (2 << (from-1));
     if(till > ITEMS) {
	  printf("More than maximum allowed, increase macro ITEMS\n");
	  return 1;
     }
 
     time_t start,end,t;
     int ret_val = 0;
     /*Run all tests*/
     int count = 0;
     for(;from <= till;from++) {
	     printf("n = %u\n",from);
	     indexa = 1;
//      while(ret_val == 0) {
	     no = yes = 0;
	  MAXVAL = (2 << (from-1));
	  
	  if(gen_rand_bids(MAXVAL))
		  break;
	  set_singleton_bid(MAXVAL);
	  start=clock();//predefined  function in c
	  run_test(MAXVAL,from);
	  end=clock();
	  t=(end-start)/CLOCKS_PER_SEC;
	  printf("\nTime taken =%lu for n= %u yes %lu no %lu\n", (unsigned long) t,from,yes,no);
/*Reset the arrays*/
	  ret_val = parse_wopt(MAXVAL);
	  memset(&f,'\0',sizeof(f));
	  memset(&O,'\0',sizeof(O));
//     }
//     }

     }
     return 0;
}
