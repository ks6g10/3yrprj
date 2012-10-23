#include <stdio.h>
#include <strings.h> // for ffs
struct _bid {
	unsigned int value;
	unsigned int comb;
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
#define MAX 8
char * assets[3] = {"apple","mapple","potato"};
/*		           0 1 2 3 4 5 6 7 8			*/
//unsigned int bids[MAX] ={0,2,3,6,4,7,6,9}; //conf 5 and 2
//unsigned int bids[MAX] =  {0,2,3,6,4,6,6,9}; //conf 3 and 4
unsigned int bids[MAX] =  {0,20,3,6,4,6,10,9}; //conf 1 and 6
unsigned int wopt[MAX] = {0};
unsigned int f[MAX] = {0};
unsigned int O[MAX] = {0};
  
inline unsigned int intersect(unsigned int seta, unsigned int setb) {
     return (seta & setb);
}

inline unsigned int _union(unsigned int seta, unsigned int setb) {
     return (seta | setb);
}

inline unsigned int setdiff(unsigned int seta, unsigned int setb) {

     return (seta & ~setb);
}

inline unsigned int cardinality(unsigned int seta) {
     return __builtin_popcount(seta);
}

/*Reminder of sets
 *f[]
 *O[]
 *bids[]	
 */

struct _bnode{
     void * ptr[2];// = {NULL};
     unsigned int value;
} typedef node;




/*Bids are arranged like following: bid[n][0] is the set config bid[n][1] is the value*/
node * parse_bids(int nitems, int nbids, int ** bids ) {
     node * root = malloc(sizeof(node));
     int i,b;
     int tmp =0;
     node * curr;
     // go through all bids
     for(b = 0; b < nbids;b++) {
	  //always start from root
	  curr = root;
	  tmp = bids[b][0];
	  //go through all bits in config
	  for(i=0;i<nitems;i++) {
	       /* if the bit i is 1*/
	       int bit = (1 & (tmp>> i));
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
     unsigned int i;
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
  register unsigned int i;
  for(i =1; i< MAX; i*=2) {
       f[i] = bids[i];
       if(bids[i] > 0)
	    O[i] = i;
  }
}

unsigned int max(unsigned conf) {
     register unsigned int card = cardinality(conf)/2;
     unsigned int i;
     unsigned int max = 0;
     unsigned int set = 0;
     unsigned int tmp = 0;
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
     unsigned int value;
} typedef leaf;

struct _stack {
	unsigned int conf;
	struct _stack * next;
} typedef stack;

void parse_wopt(void) {
     unsigned int c;
     //wopt at start contain MAX at wopt[0] which is the combination that goes in bids[wopt[n]]
     stack * root = malloc(sizeof(stack));
     root->conf = MAX-1;
     stack * curr = root;
     while(curr != NULL) {
	     unsigned int conf = curr->conf;
	     if(conf != O[conf]) {
		     unsigned int diff = setdiff(conf,O[conf]);
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

int main(void) {
      register unsigned int i;
      /*1.*/
      set_singleton_bid();
       printfo();
      /*2.*/
      //register unsigned int i;
      for(i = 2; i <MAX; i++) {
	   unsigned int c;
	   for(c = 1; c < MAX; c++) {
		if(cardinality(c) == i && bids[c] > 0) {
		     printf("hit card %u \n",i);
		     unsigned int tmpset = max(c);
		     if(f[c] >= bids[c]) {//b
			     O[c] = tmpset;//net to set
		     }
		     else {//c
			  f[c] = bids[c];
			  O[c] = c;
		     }
		     printfo();
		}
	   }
      }

      parse_wopt();
      //printfo();

      return 0;
      /*Testing facility*/
	 for(i = 1; i < 8; i++) {
		 
		 unsigned     int z = intersect(i,i-1);
		 unsigned     int x = _union(i,i-1);
		 unsigned     int f = setdiff(i,i-1);
		 unsigned int t = cardinality(i);
		 printf("i %d \tinter %u\t union %u\t diff %	u\t card %u\n",i,z,x,f,t);
	 }
	 return 0;
} 
