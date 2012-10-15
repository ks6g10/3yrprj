#include <stdio.h>

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
 */

/*                    1         2        4       8        16*/
#define MAX 8
char * assets[3] = {"apple","mapple","potato"};
/*		      0 1 2 3 4 5 6 7 8			*/
unsigned int bids[MAX] {0,2,3,6,2,3,5,7}; 
unsigned int oset[MAX] {0,1,2,3,4,5,6,7};
unsigned int f[MAX] = {0};
unsigned int O[MAX] = {0};

int[] set_singleton_bid(register int[] _bids) {
	register unsigned int i;
	register unsigned int _f[MAX];
	for(i =1; i< MAX; i*=2) {
		_f[i] = _bids[i]
	}
	return _f;
}
  
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

/*
 *f[]
 *O[]
 *	
 */

 int main(void) {
	 /*1.*/
	 f = set_singleton_bid(bids);
	 O = set_singleton_bid(oset);
	 /*2.*/
	 register unsigned int i;
	 for(i = 2; i <MAX; i++) {




	 }
	 return 0;
	 unsigned int i;
	 for(i = 1; i < 8; i++) {
		 
		 unsigned     int z = intersect(i,i-1);
		 unsigned     int x = _union(i,i-1);
		 unsigned     int f = setdiff(i,i-1);
		 unsigned int t = cardinality(i);
		 printf("i %d \tinter %u\t union %u\t diff %	u\t card %u\n",i,z,x,f,t);
	 }
	 return 0;
} 
