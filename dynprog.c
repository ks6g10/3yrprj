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
 */

/*                    1         2        4       8        16*/
#define MAX 8
char * assets[3] = {"apple","mapple","potato"};
/*		         0 1 2 3 4 5 6 7 8			*/
unsigned int bids[MAX] ={0,2,2,2,2,2,2,2}; 
unsigned int oset[MAX] = {0,1,2,3,4,5,6,7};
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
//     printf("o\t");
     //   for(i =1; i < MAX; i++) 
//	  printf("%u\t",O[i]);
     //   printf("\n");
}

 int main(void) {
      register unsigned int i;
      /*1.*/
      set_singleton_bid(bids,f);
      set_singleton_bid(oset,O);
      printfo();
      /*2.*/
      //register unsigned int i;
      for(i = 2; i <MAX; i++) {
	   unsigned int c;
	   for(c = 1; c < MAX; c++) {
		if(cardinality(c) == i && bids[c] > 0) {
		     printf("hit card %u \n",i);
		     unsigned int tmpset = c & (~c + 1);
		     f[c] = f[setdiff(c,tmpset)] + f[tmpset];//a
		     if(f[c] >= bids[c])//b
			  O[c];//net to set
		     else {//c
			  f[c] = bids[c];
			  O[c] = c;
		     }
		     printfo();
		}
	   }
      }


      printfo();

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
