#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#define WORDSIZE (sizeof(unsigned int)*8)




//typedef unsigned int int_fast32_t;

struct configuration {
	unsigned int * allocation;
	unsigned int * bin;
	unsigned int * id;
	// each bids dummy, if any else is set to 0
	unsigned int * dummies;
	unsigned int * bin_to_offset;
	unsigned int * bin_count;
	float * offer;
	float * score;
	float * average;
	// the allocation dummy	
	unsigned int goods;
	unsigned int words;
	unsigned int bids;
	unsigned int dummy;
	unsigned int singletons;
	unsigned int * curr_allocation;
	unsigned int * allocation_id;
	unsigned int * allocation_dummy;
	unsigned int allocation_id_index; 
	float allocation_value;
};


 unsigned int ints =0;
/* unsigned int * MASK = void; */

/* double h(struct allocation * pi,struct bid2 ** bins,unsigned int * bin_counts) { */
/* 	unsigned int u_bound = 0; */
/* //	unsigned int free[ints]; */
/* 	unsigned int * pi_conf = pi->conf; */
/* 	int x; */
/* //	double * means[goods]; */
/* 	double upper_bound =0; */
/* 	unsigned int tmp; */
/* 	for(x = 0; x < ints; x++) { */
/* 		tmp  = ~pi_conf[x] & MASK[x]; */
/* 		while(tmp) { */
/* 			unsigned int index = __builtin_ffs(tmp) -1; */
/* 			tmp  &= ~(1 << index); */
/* 			upper_bound +=  v(bins[x],pi_conf,bin_counts[x]); */
/* 		} */
/* 	} */
/* 	return u_bound; */
/* } */

/* oid init_mask(void) { */
/* 	MASK = calloc(sizeof(int),ints); */
/* 	int x; */
/* 	unsigned int mod,tmp; */
/* 	for(x =0; x < ints) { */
/* 		if(x < (ints-1)) { */
/* 			MASK[x] = ~0; */
/* 		} else { */
/* 			mod = NGOODS % (INT); */
/* 			tmp = (1 << mod)-1; */
/* 			MASK[x] = tmp; */
/* 		} */
/* 	}	       	 */
/* } */

/* int do_not_intersect(struct bid2 * a_bid, unsigned int * conf) { */
/* 	int x; */
/* 	unsigned int tmp; */

/* 	for(x = 0; x < ints; x++)  { */
/* 		tmp = a_bid->conf[x] & conf[x]; */
/* 		if(tmp) { */
/* 			return 0; */
/* 		} */
/* 	} */
/* 	return 1; */
/* } */

/* double v(struct bid2 * bin,unsigned int * conf,unsigned int bin_count) { */
/* 	if(!bin_count) { */
/* 		return 0.0; */
/* 	} */
/* 	int x; */
/* 	for(x =0; x < bin_count;x++) { */
/* 		if(do_not_intersect( bin[x], conf )) { */
/* 			return  */
/* 		} */
/* 	} */
/* } */

unsigned int next_index(unsigned int a_index) {		
	return __builtin_ffs(a_index) -1;
//	tmp  &= ~(1 << index);
//	upper_bound +=  v(bins[x],pi_conf,bin_counts[x]);
	
}


void print_binary(unsigned int * allocation,unsigned int goods) {
	int x;
	for(x=goods-1;x>=0;x--) {
		printf("%u", !!(*allocation & (1 << x) ) );
	}
	printf("\n");
}




struct configuration * get_configuration(FILE * fp) {
	const char * s_goods = "goods";
	const char * s_bids = "bids";
	const char * s_dummy = "dummy";

	ssize_t read;
	char * line = NULL;
	size_t len = 0;
	unsigned int all = 0;
	struct configuration * ret = malloc(sizeof(struct configuration));
	while ((read = getline(&line, &len, fp)) != -1 && !all) {
		if(line[0] == '%' || line[0] == '\n') {
			continue;
		}
		if(strncmp(line,s_goods,strlen(s_goods)) == 0) {      
			ret->goods = atoi(line+strlen(s_goods)+1);
			printf("Number of goods %u\n",ret->goods);				
		} else if(strncmp(line,s_bids,strlen(s_bids)) == 0) {
			ret->bids = atoi(line+strlen(s_bids)+1);				
			printf("Number of bids %u\n",ret->bids);
		} else if(strncmp(line,s_dummy,strlen(s_dummy)) == 0) {
			ret->dummy = atoi(line+strlen(s_dummy)+1);
			printf("Number of dummy %u\n",ret->dummy);
//			ints = 1+(goods-1)/32;
		}
//			total_goods = goods + dummy;
		all = !!(ret->goods && ret->bids && ret->dummy);	 
	}
	ret->words = 1+((ret->goods-1)/WORDSIZE);
	return ret;
}

unsigned int * get_bincount(FILE * fp,
			    struct configuration * conf,
			    unsigned int * have_singelton) {
	unsigned int * bin_count = malloc(sizeof(int)*conf->goods);
	int x;
	for(x=0;x<conf->goods;x++) {
		bin_count[x] = 0;
	}
	
	char * head, * tail,*line;
	size_t len =0;
	ssize_t read;
	while ((read = getline(&line, &len, fp)) != -1) {
		if(!isdigit(line[0])){ continue;}		
		head = tail = line;
		int tab_count = 0;
		head++;
		while(tab_count < 3) {
			if(*head == '\t') {
				tab_count++;
				if(tab_count <= 2) {
					tail = head;
					head++;
				}
			} else {
				head++;
			}
		}
		int which_bin = strtol(tail,&head,10);
		head++;
		bin_count[which_bin]++;
		int goods_count =1;
		while(*head != '#' && *head != '\0') {
			if(*head == '\t') {
				tail = head;
				goods_count++;
			}
			head++;
		}
		if(goods_count == 1) {			
			have_singelton[which_bin] = 1;
		}
//		printf("Bin %d count %u\n",which_bin,bin_count[which_bin]);
	}
	return bin_count;       
}

void allocate_all_bids(FILE * fp,
		       struct configuration * conf,
		       unsigned int * have_singelton,
		       unsigned int * bin_count) {

	conf->allocation = malloc(sizeof(unsigned int)*conf->bids);
	conf->bin = malloc(sizeof(unsigned int)*conf->bids);
	conf->id = malloc(sizeof(unsigned int)*conf->bids);
	conf->dummies = malloc(sizeof(unsigned int)*conf->bids);
	conf->bin = malloc(sizeof(unsigned int)*conf->bids);
	conf->bin_to_offset = malloc(sizeof(unsigned int)*conf->goods);
	conf->offer = malloc(sizeof(float)*conf->bids);
	conf->average = malloc(sizeof(float)*conf->bids);
	conf->score = malloc(sizeof(float)*conf->goods);
	conf->bin_count = bin_count;
	char * head, * tail,*line; 
	unsigned long total_goods_count[conf->goods];

	unsigned int numbids_count[conf->goods];
	unsigned int goods[conf->goods];	
	size_t len =0;
	ssize_t read;	
	int x;
	for(x = 0; x< conf->goods;x++) {
		conf->bin_to_offset[x] = 0;
	}
	unsigned int bin_index[conf->goods];
	bin_index[0] = 0;

	for(x = 1; x< conf->goods;x++) {
		bin_index[x] = bin_count[x-1]+bin_index[x-1];
	}

	while ((read = getline(&line, &len, fp)) != -1) {
		if(!isdigit(line[0])) {continue;}
		head = tail = line;

		while(*head != '\t' && *head != '\0') {
			head++;
		}
		int id = strtol(tail,&head,10);
		tail = head;
		head++;
		//get offer or value
		while(*head != '\t' && *head != '\0') {
			head++;
		}
		float offer = strtod(tail,&head);
		tail = head;
		head++;
		unsigned int goods_count = 0;		
		unsigned int good = 0;	
		unsigned int dummy = 0;
		unsigned int tmp_allocation[conf->words];
		for(x = 0;x<conf->words;x++) {
			tmp_allocation[x] = 0;
		}
		//reset the temporary goods array, used to determin the score
		goods_count = 0;
		for(x = 0;x<conf->goods;x++) {
			goods[x] = 0;
		}

		while(*head != '#' && *head != '\0') {
			if(*head == '\t') {
				good = strtol(tail,&head,10);
				
				//sscanf(tail,"\t%u\t",tmp2);
				if(good < conf->goods) {
					tmp_allocation[(good/WORDSIZE)] += (1 << good);
					if(id == 13) {
					printf("good %u ",tmp_allocation[0]);
				}
					goods[goods_count] = good;	
				} else {
					dummy = good;
				}
				tail = head;
				goods_count++;
			}
			head++;
			
		}
		int index = bin_index[goods[0]];
		conf->average[index] = offer/goods_count;
		for(x=0;x<goods_count;x++) {
			total_goods_count[goods[x]] += goods_count;
			numbids_count[goods[x]]++;
		}
		conf->offer[index] = offer;
		conf->bin[index] = goods[0];
		conf->dummies[index] = dummy;
		conf->id[index] = id;
		for(x = 0;x<conf->words;x++) {
			conf->allocation[index*conf->words+x] = tmp_allocation[x];
		}
//		printf("id %d bin %u count %u value %.3lf\n",bid_count,bin_for_bid,tmp_count[bin_for_bid],0);
		bin_index[goods[0]]++;
	}
	int singleton_count = conf->bids-conf->singletons;	
	for(x =0;x< conf->goods;x++) {
		if(!have_singelton[x]) {
			int index = bin_index[x];
			int y;
			for(y = 0;y<conf->words;y++) {
				conf->allocation[index*conf->words+y] = 0;
			}
			int word_index = x/WORDSIZE;
			int bit_index = x % WORDSIZE;			
			conf->allocation[index*conf->words+word_index] = (1<<bit_index);
			conf->average[index] = 0;
			conf->bin[index] = x;
			conf->offer[index] = 0;
			conf->dummies[index] = 0;
			conf->id[index] = singleton_count;
			total_goods_count[x] +=1; //add one more to the score stat
			numbids_count[x] +=1; // also add one more to the number of bids to the score stat
			singleton_count++; // next singleton bid will have an consecutive bid id
		}
		double score =0;
		double avg;
		if(numbids_count[x]) {
		avg = ((double)total_goods_count[x])/((double)numbids_count[x]);
		score = ((double)numbids_count[x])/avg;
		}
		conf->score[x] = score;
	}
	return ;	
}


//wordsize in bits
#define WORD (32)
#define BIN (0)
#define INDEX (1)

int is_conflict(struct configuration * conf,
		unsigned int index_to_allocate) {
	int x = 0;
	int ret = 0;
	int id = index_to_allocate;
	//check if the allocation is conflicting
	for(x=0;x < conf->words;x++)  {
		ret |= (conf->curr_allocation[x] & conf->allocation[id*conf->words+x]);		
	}
	int size = conf->allocation_id_index;
	unsigned int dummy = conf->dummies[index_to_allocate];
	//special case if there is no dummy good
	if(dummy ==0) {return ret;}
	
	//check if the dummy bids are conflicting	
	for(x=0;x<size;x++) {
		if(dummy == conf->allocation_dummy[x]) {
			ret = 1;
			break;
		}
	}
	return ret;
	

}

void print_debug(struct configuration * conf,
		 unsigned int (* allocation_count), 
		 unsigned int allocation_id_index,
		 unsigned int low_order_good,
		 unsigned int (* bin_count)) {
	return;
	printf("low order good %u allocation_id_index %u\n",low_order_good,allocation_id_index);
	int x;
	printf("allocation:\n");
	for(x=0;x < conf->goods;x++) {
		printf("%u/%u\t",allocation_count[x],bin_count[x]);
	}
	printf("\n\n");


}


#define _BID(X) (bins[X].bids[count[X]])
#define BID (_BID(good))
int dealloc_id(struct configuration * conf) {
	conf->allocation_id_index--;
	unsigned int allocation_id_index = conf->allocation_id_index;
	unsigned int index_to_deallocate = conf->allocation_id[allocation_id_index];
	int x;
	for(x =0; x < conf->words;x++) {
		int index = index_to_deallocate*conf->words+x;
		conf->curr_allocation[x] &= ~(conf->allocation[index]);
	}
	conf->allocation_value -= conf->offer[index_to_deallocate];
	conf->allocation_dummy[conf->allocation_id_index] = 0;//conf->dummies[index_to_deallocate];
//	int good = 
	return conf->bin[index_to_deallocate];
//	conf->bin_count[good]++;

}

void alloc_id(unsigned int index_to_allocate,struct configuration * conf) {

	int x;
	for(x =0; x < conf->words;x++) {
		int index = index_to_allocate*conf->words+x;
		conf->curr_allocation[x] |= conf->allocation[index];
	}
	conf->allocation_value += conf->offer[index_to_allocate];
	conf->allocation_id[conf->allocation_id_index] = index_to_allocate;
	conf->allocation_dummy[conf->allocation_id_index] = conf->dummies[index_to_allocate];
	
	conf->allocation_id_index++;
}




void calc_best(struct configuration * conf) {
	const unsigned int goods = conf->goods;
	const unsigned int words = conf->words;
	unsigned int count[goods]; 
	unsigned int * bin_count = conf->bin_count;
	unsigned int * allocation = conf->curr_allocation = malloc(sizeof(int) *conf->words);
	//0 = bin, 1 index
	conf->allocation_id = malloc(sizeof(int) *conf->words);
	conf->allocation_id_index = 0;
	unsigned int * allocation_dummy = conf->allocation_dummy = malloc(sizeof(int) *conf->words);
	double max = 0;
	int x;	
	//initilize to 0	
	for(x = 0; x < goods; x++){
		count[x] = 0;
		conf->allocation_dummy[x] = 0;     
	}
	for(x = 0;x < words;x++) {
		allocation[x] = 0;
	}

	unsigned int bin_index[conf->goods];
	bin_index[0] = 0;

	for(x = 1; x< conf->goods;x++) {
		bin_index[x] = bin_count[x-1]+bin_index[x-1];
		printf("x %d bin_index %u\n",x,bin_index[x]);
	}
	printf("hello\n");
	int good = 0;//order[1][0]; //initiate at the lowest order good;
	int index_to_allocate =0;// = BID.id;
	conf->allocation_value =0;// BID.offer;
	

	while(1) {

	lbl_continue:
		;
//		getchar();
		// which word the good is located in
		int word_index = good/WORDSIZE; 
		// which bit position in the word the good represent
		int bit_index = good % WORDSIZE;
		int index = 0;
		while((allocation[word_index] & (1 << bit_index)) && (good < goods)) {
			printf("good %u not compat\n",good);
			good++;
			word_index = good/WORDSIZE;
			bit_index = good % WORDSIZE;			
		}

		if((good >= goods)) {
			printf("backtrack 2\n");
			goto backtrack;
		};
		
		do{	      
			//this means we have past the singleton
			if(count[good] >= bin_count[good]) {
				//hence reset this value right here ma'm
				printf("backtracking good %u id %u %u\n",
				       good,
				       conf->allocation_id[conf->allocation_id_index]);
				goto backtrack;
			}
			//printf("good %u count %u max %u\n",good,count[good], bin_count[good]);
			index = count[good];
			index_to_allocate = index+bin_index[good];			
			count[good] +=1;
			/* for(x=0;x<words;x++) { */
			/* 	is_conflict |= (allocation[x] & conf->allocation[x+index_to_allocate*words]); */
			/* 	printf("intersection %u, alloc %u, bid %u\n", */
			/* 	       (allocation[x] & conf->allocation[x+index_to_allocate*words]), */
			/* 	       allocation[x], */
			/* 	       conf->allocation[x+index_to_allocate*words]); */
			/* } */
			
		}while(is_conflict(conf,index_to_allocate));
		
		
		alloc_id(index_to_allocate,conf);
		printf("allocated id %u\n",conf->id[index_to_allocate]);
		print_binary(allocation,goods);
		if(conf->allocation_value > max) {
			max = conf->allocation_value;
			printf("New max of %.3lf \n",max);
			for(x=0;x<goods;x++) {
				printf("%u ",count[x]);
			}
			printf("\n");
		}			

		goto lbl_continue;
	backtrack:
		;
		int backtrack;
		do {

			backtrack = 0;			
			int good = dealloc_id(conf);
			
			/* printf("dealloc index %u bin %u aloc_index %u,bin count %u/%u\n", */
			/*        index_to_dealloc, */
			/*        good, */
			/*        conf->allocation_id_index, */
			/*        count[good], */
			/*        bin_count[good]); */

			if(count[good] >= conf->bin_count[good]) {
				printf("good %u\n",good);
				if(good == 0) {
					printf("New max of %.3lf \n",max);
					return;
				}

				count[good] =0;
				backtrack = 1;
			}
			good = 0;


		}while(backtrack);
	}
}


int main(int argc, char *argv[])   {
	int x;
	FILE * fp;
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	printf("hello\n");
	struct configuration * conf = get_configuration(fp);
	conf->singletons =0;
	unsigned int * have_singleton =  malloc(sizeof(int)*conf->goods);
	for(x =0; x < conf->goods;x++) {
		have_singleton[x] = 0;
	}
	unsigned int * bin_count = get_bincount(fp,conf,have_singleton);
	     
	for(x =0; x < conf->goods;x++) {
		if(!have_singleton[x]) {
			conf->singletons++;
			conf->bids++;
			bin_count[x]++;
		}
		printf("bin %d count %u\n",x,bin_count[x]);
	}
	
	fclose(fp);
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	allocate_all_bids(fp,conf,have_singleton,bin_count);
	fclose(fp);
	for(x=0;x < conf->bids;x++) {
		printf("x %d id %u, offer %.3f, bin %u, alloc %u\n",
		       x,
		       conf->id[x],
		       conf->offer[x],
		       conf->bin[x],
		       conf->allocation[x]);

	}	
	printf("words %u wordsize %lu\n",conf->words,WORDSIZE);

	calc_best(conf);
	printf("Bye\n");
exit(EXIT_SUCCESS);
}
