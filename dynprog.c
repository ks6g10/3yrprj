#include <stdio.h>

/*                  1         2      4       8        16*/
char * assets[5] = {"apple","mapple","pear","potato","asdf"};
unsigned int set1 = 1; //apple pear potato
unsigned int set2 = 6; // mapple pear potato


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


 int main(void) {
      unsigned int i;
      for(i = 1; i < 8; i++) {
      
      unsigned     int z = intersect(i,i-1);
      unsigned     int x = _union(i,i-1);
      unsigned     int f = setdiff(i,i-1);
      unsigned int t = cardinality(i);
      printf("i %d \tinter %u\t union %u\t diff %u\t card %u\n",i,z,x,f,t);
      }
      return 0;
} 
