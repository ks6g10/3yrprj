#include <stdio.h>
#include <strings.h> // for ffs
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#define ITEMS 20
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

struct _bid {
	dint value;
	 dint comb;
} typedef bid;
/*0000 0 0
 *0001 1 1 1
 *0010 1 2 2
 *0011 2 3 2
 *0100 1 4 4
 *0101 2 5 1
 *0110 2 6 3
 *0111 3 7
 *1000 1 8 8
 *1001 2 9 1
 *1010 2 10 2
 *1011 3 11
 *1100 2 12 
 *1101 3 13
 *1110 3 14
 *1111 4 15
 *10000  
 *
 *        111  
 *        9
 *        /|\
 *    110 101 011
 *     10   6   6  
 *     |\ / \ /|      
 *    100 010 001
 *     4   3   20
 *
 *  
 * 
 *
 *
 */

/*                    1         2        4       8        16*/

char * assets[3] = {"apple","mapple","potato"};
/*		           0 1 2 3 4 5 6 7 8			*/
// dint bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
// dint bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
// dint bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6
dint bids[MAX];
 dint wopt[MAX] = {0};
 dint f[MAX] = {0};
 dint O[MAX] = {0};
  
inline  dint intersect( dint seta,  dint setb) {
     return (seta & setb);
}

inline  dint _union( dint seta,  dint setb) {
     return (seta | setb);
}

inline  dint setdiff( dint seta,  dint setb) {

     return (seta & ~setb);
}

inline  dint cardinality( dint seta) {
     return __builtin_popcount(seta);
}

/*Reminder of sets
 *f[]
 *O[]
 *bids[]	
 */

void gen_rand_bids() {
	register dint i = 0;
	printf("i =\t");
	for(i = 1; i < MAX;i++) {
		bids[i] = 1;//rand() % 10;
		//printf("%u\t",i);
	}
	return;
		printf("\n");
		printf("val =\t");
	for(i = 1; i < MAX;i++) {
		printf("%u\t",bids[i]);
	}
	printf("\n");
}

struct _bnode{
     void * ptr[2];// = {NULL};
      dint value;
} typedef node;




/*Bids are arranged like following: bid[n][0] is the set config bid[n][1] is the value*/
node * parse_bids(dint nitems, dint nbids, dint ** bids ) {
     node * root = malloc(sizeof(node));
     dint i,b;
     dint tmp =0;
     node * curr;
     // go through all bids
     for(b = 0; b < nbids;b++) {
	  //always start from root
	  curr = root;
	  tmp = bids[b][0];
	  //go through all bits in config
	  for(i=0;i<nitems;i++) {
	       /* if the bit i is 1*/
	       dint bit = (1 & (tmp>> i));
	       if(curr->ptr[bit] == NULL) {
		    curr->ptr[bit] = malloc(sizeof(node));
	       }
	       curr = curr->ptr[bit];
	       //  continue if there is more bits to parse
	       if(i < nitems-1)
		    continue;
	       //else set the value
	       curr->value = bids[b][1];
	       break;
	  }
     }
}



void printfo() {
	return;
      dint i;
     printf("i\t");
     for(i =1; i < MAX; i++) 
	  printf("%u\t",i);
     printf("\n");
     printf("f\t");
     for(i =1; i < MAX; i++) 
	  printf("%u\t",f[i]);
     printf("\n");
     printf("o\t");
     for(i =1; i < MAX; i++) 
	  printf("%u\t",O[i]);
     printf("\n");
}

inline void set_singleton_bid() {
  register  dint i;
  for(i =1; i< MAX; i*=2) {
       f[i] = bids[i];
       if(bids[i] > 0)
	    O[i] = i;
  }
}

dint max(dint conf) {
     register  dint card = cardinality(conf)/2;
      dint i;
      dint max = 0;
      dint set = 0;
      dint tmp = 0;
     for(i=1;i<MAX;i++) {
	  if(cardinality(i) > card)
	       continue;
	  if(i != (i&conf))
	       continue;
	  tmp = f[setdiff(conf,i)] + f[i];
	  if(max < tmp) {
	       max = tmp;
	       set = i;
	  }
     }
     f[conf] = max;
     return set;
}

struct _bleaf {
     void * ptr[2];// = {NULL};
     union{
	  void * ptr0;// = NULL;
	  void * ptr1;// = NULL;
     };
      dint value;
} typedef leaf;

struct _stack {
	 dint conf;
	struct _stack * next;
} typedef stack;

void parse_wopt(void) {
     //wopt at start contain MAX at wopt[0] which is the combination that goes in bids[wopt[n]]
     stack * root = malloc(sizeof(stack));
     //DO NOT REMOVE -1
     root->conf = (MAX)-1;
     stack * curr = root;
     while(curr != NULL) {
	      dint conf = curr->conf;
	     if(conf != O[conf]) {
		     dint diff = setdiff(conf,O[conf]);
		     curr->conf = O[conf];
		     stack * tmp = malloc(sizeof(stack));
		     tmp->conf = diff;
		     tmp->next = curr;
		     root = tmp;
		     curr = root;
		     continue;
	     }
	     curr = curr->next;
     }
     curr = root;
     while(curr != NULL) {
	     printf("conf %u value %u\n",curr->conf,bids[curr->conf]);
	     curr = curr->next;
     }
}

dint main(void) {
      printf("max %d\n", MAX);
      register dint i;
      /*1.*/
      gen_rand_bids();
      set_singleton_bid();
      //printfo();
      /*2.*/
      for(i = 2; i <MAX; i++) {
	       dint c;
	      for(c = 1; c < MAX; c++) {
		      if(cardinality(c) == i && bids[c] > 0) {
			      //printf("hit card %u \n",i);
			       dint tmpset = max(c);
			      if(f[c] >= bids[c]) {//b
				      O[c] = tmpset;//net t	o set
			      }
			      else {//c	
				      f[c] = bids[c];
				      O[c] = c;
			      }
			      //printfo();
		      }
	      }
      }
      	     	      	     
      parse_wopt();
      //printfo();

      return 0;
      /*Testing facility*/
	 for(i = 1; i < 8; i++) {
		 
		 dint z = intersect(i,i-1);
		 dint x = _union(i,i-1);
		 dint f = setdiff(i,i-1);
		 dint t = cardinality(i);
		 printf("i %u \tinter %u\t union %u\t diff %u\t card %u\n",i,z,x,f,t);
	 }
	 return 0;
} 
