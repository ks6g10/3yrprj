#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#define WORDSIZE (sizeof(int_fast32_t)*8)

union asd {
	unsigned int as;
	char ds[3];
};

struct bid_ptr {
	unsigned int index;
	double value;
	unsigned int dummy;
};

struct bid2 {
	double offer;
	double average;
	unsigned int bin;
	unsigned int id;
	unsigned int dummy;
};
//typedef unsigned int int_fast32_t;

struct configuration {
	unsigned int * allocation;
	unsigned int * bin;
	unsigned int * id;
	unsigned int * dummies;
	unsigned int * bin_to_offset;
	unsigned int * bin_count;
	float * offer;
	float * score;
	float * average;
	unsigned int goods;
	unsigned int words;
	unsigned int bids;
	unsigned int dummy;
	unsigned int singletons;
};

struct bid_bin {
	unsigned int size;
	unsigned int good;
	double average;
	double score;
	struct bid2 * bids;
};

struct root_bid {
	unsigned int goods;
	struct bid_bin * bins;
};

/* struct configuration { */
/* 	unsigned int goods; */
/* 	unsigned int bids; */
/* 	unsigned int dummy; */
/* 	unsigned int words; */
/* 	unsigned int singletons; */
/* }; */

struct linked_bids {
	struct bid2 * this;
	struct linked_bids * next;
};

struct allocation {
	struct linked_bids * root;

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
					tmp_allocation[(good/WORDSIZE)] |= (1 << good);
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
			conf->allocation[index*2+x] = tmp_allocation[x];
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
				conf->allocation[index*2+y] = 0;
			}
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
		printf("id\n");
	}
	return ;	
}

struct bid_bin * allocate_bid_bin(FILE * fp,
				  unsigned int * bin_count,				  
				  struct configuration * conf,				
				  unsigned int (* bid_conf)[conf->words],
				  unsigned int * have_singelton){
				  
	struct bid_bin * bins = malloc(sizeof(struct bid_bin)*conf->goods);	
	char * head, * tail,*line; 
	unsigned long total_goods_count[conf->goods];
	unsigned int numbids_count[conf->goods];
	unsigned int goods[conf->goods];	
	unsigned int tmp_count[conf->goods];	
	size_t len =0;
	ssize_t read;
	int x;
	for(x = 0; x< conf->goods;x++) { 
		tmp_count[x] = 0;
		bins[x].bids = malloc((sizeof(struct bid2))*bin_count[x]);
		total_goods_count[x] = 0;
		bins[x].average = 0;
		numbids_count[x] = 0;
	}
	while ((read = getline(&line, &len, fp)) != -1) {		
		if(!isdigit(line[0])) {continue;}
		head = tail = line;
		while(*head != '\t') {head++;}


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
		double value = strtod(tail,&head);
		tail = head;
		head++;
		
		unsigned int dummy_good = 0; 
		unsigned int goods_count = 0;
		unsigned int bin_for_bid = 0;
		unsigned int good = 0;
		unsigned int tmp_conf[conf->words];
		for(x = 0;x<conf->words;x++) {
			tmp_conf[x] = 0;
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
					tmp_conf[good/32] |= (1 << good);
					goods[goods_count] = good;
				} else {
					dummy_good = good;
				}
				tail = head;
				goods_count++;
				if(goods_count == 1) {					
					bin_for_bid = good;
				}
			}
			head++;
		}
		
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].offer = value;
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].id = id;
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].dummy = dummy_good;
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].bin = bin_for_bid;
		double tmp_average = value/((double)goods_count);
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].average = tmp_average;
		if(tmp_average > bins[bin_for_bid].average) {
			bins[bin_for_bid].average = tmp_average;
		}

		for(x=0;x< conf->words;x++) {
			bid_conf[id][x] = 0;
		}

		for(x = 0;x < conf->words; x++) {
			bid_conf[id][x] = tmp_conf[x];
//			bins[bin_for_bid].bids[tmp_count[bin_for_bid]].conf[x] = tmp_conf[x];
		}
		for(x=0;x<goods_count;x++) {
			total_goods_count[goods[x]] += goods_count;
			numbids_count[goods[x]]++;
		}
		tmp_count[bin_for_bid] +=1;
		printf("id %d bin %u count %u value %.3lf\n",id,bin_for_bid,tmp_count[bin_for_bid],value);
	}
	int singleton_count = conf->bids-conf->singletons;
	int y;
	for(x =0;x< conf->goods;x++) {
		if(!have_singelton[x]) {
			bins[x].bids[tmp_count[x]].offer = 0;
			bins[x].bids[tmp_count[x]].id = singleton_count;
			bins[x].bids[tmp_count[x]].dummy = 0;
			bins[x].bids[tmp_count[x]].average = 0;
			bins[x].bids[tmp_count[x]].bin = x;
			total_goods_count[x] +=1; //add one more to the score stat
			numbids_count[x] +=1; // also add one more to the number of bids to the score stat
			for(y=0;y< conf->words;y++) {
				bid_conf[singleton_count][y] = 0;
			}
			int tmp_index = x/32;
			int tmp_bit = x%32;
			bid_conf[singleton_count][tmp_index] = (1<<tmp_bit);
			singleton_count++; // next singleton bid will have an consecutive bid id
		}
		double score =0;
		double avg;
		if(numbids_count[x]) {
		avg = ((double)total_goods_count[x])/((double)numbids_count[x]);
		score = ((double)numbids_count[x])/avg;
		}
		bins[x].score = score;
	}
	return bins;
}

unsigned int ** get_bin_order(struct bid_bin * bins,struct configuration * conf,unsigned int * bin_count) {
	unsigned int ** order = malloc(sizeof(int *)*2);	
	order[0] = malloc(sizeof(int)*conf->goods);
	order[1] = malloc(sizeof(int)*conf->goods);
	int x;
	for(x=0;x<conf->goods;x++) {
		//	score[x] = bins[x].score;
		order[0][x] = conf->goods;
	}
	
	for(x = 0;x<conf->goods;x++) {
		order[0][x] = x;
		order[1][x] = x;
		/* max_index = conf->goods; */
		/* max = 999999999999; */
		/* for(y=0;y<conf->goods;y++){ */
		/* 	if(score[y] < max && bin_count[y] >0) { */
		/* 		max_index = y;					 */
		/* 		max = score[y]; */
		/* 	} */
			
		/* } */
		/* if(max_index < conf->goods) { */
		/* 	score[max_index] = 99999999999999999;	 */
		/* 	order[0][max_index] = x; */
		/* 	order[1][x] = max_index; */
		/* } */
	}
	return order;
}

//wordsize in bits
#define WORD (32)
#define BIN (0)
#define INDEX (1)

int is_conflict(struct configuration * conf,
		unsigned int id,		  
		  unsigned int (* allocation)) {
	int x = 0;
	int ret = 0;
//	printf("hello id %u\n",id);
	

	for(x=0;(x<conf->words);x++)  {
		if((allocation[x] & conf->allocation[id*2+x]) !=0) {
			ret = 1;
			return ret;
		}
	}

	return (ret);
}

int is_dummy_conflict(unsigned int dummy,
		unsigned int (* dummies),
		unsigned int size) {
	int x;
	int ret = 0;
	//special case if there is no dummy good
	if(dummy ==0) {return ret;}

	for(x=0;x<size;x++) {
		if(dummy == dummies[x]) {
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

double get_estimate (struct configuration * conf,
		     struct bid_bin * bins,
		     unsigned int (* allocation)) {

	double estimate = 0;
	int x;
	for(x = 0;x < conf->goods;x++) {
		int bit = x % WORD;
		int index = x / WORD;
		if(!(allocation[index] & (1<<bit))) {
			estimate += bins[x].average;
		}
		
	}
	return estimate;
}

#define _BID(X) (bins[X].bids[count[X]])
#define BID (_BID(good))
void calc_best(struct configuration * conf) {
	const unsigned int goods = conf->goods;
	const unsigned int words = conf->words;
	unsigned int count[goods]; 
	unsigned int * bin_count = conf->bin_count;
	unsigned int allocation[conf->words];
	//0 = bin, 1 index
	unsigned int allocation_id[2][goods];
	unsigned int allocation_id_index = 0;
	unsigned int allocation_dummy[goods];
	double max = 0;
	int x;	
	//initilize to 0	
	for(x = 0; x < goods; x++){
		count[x] = 0;
		allocation_dummy[x] = 0;     
	}
	for(x = 0;x < words;x++) {
		allocation[x] = 0;
	}

	unsigned int bin_index[conf->goods];
	bin_index[0] = 0;

	for(x = 1; x< conf->goods;x++) {
		bin_index[x] = bin_count[x-1]+bin_index[x-1];
	}
	printf("hello\n");
	int good = 0;//order[1][0]; //initiate at the lowest order good;
	int index_to_allocate;// = BID.id;
	double allocation_value =0;// BID.offer;
	

	while(1) {

	lbl_continue:
		;

		// which word the good is located in
		int word_index = good/WORDSIZE; 
		// which bit position in the word the good represent
		int bit_index = good % WORDSIZE;
		int index = 0;
		while((allocation[word_index] & (1 << bit_index)) && (good < goods)) {
			good++;
			word_index = good/WORDSIZE;
			bit_index = good % WORDSIZE;
		}

		if(!(good < goods)) {
			goto backtrack;
		};

		do{
			//this means we have past the singleton
			if(count[good] >= bin_count[good]) {
				//hence reset this value right here ma'm
				count[good] = 0;
				goto backtrack;
			}
			//printf("good %u count %u max %u\n",good,count[good], bin_count[good]);
			index = count[good];
			index_to_allocate = index+bin_index[good];			
			count[good] +=1;
			
		}while(is_conflict(conf,
				   index_to_allocate,
				   allocation) ||
		       is_dummy_conflict(conf->dummies[index_to_allocate],
					 allocation_dummy,
					 allocation_id_index));
		
		
		//allocate bid
		for(x=0;x<words;x++) {
			allocation[x] |= conf->allocation[x+index_to_allocate*2];
		}
//		printf("allocating id %u\n",id_to_allocate);
 
		printf("alloc good %u,id %u offer %0.3lf\n",good,index_to_allocate,conf->offer[index_to_allocate]);
		allocation_value += conf->offer[index_to_allocate];// bid_value;

		if(allocation_value > max) {
			max = allocation_value;
			printf("New max of %.3lf \n",max);
						for(x=0;x<goods;x++) {
				printf("%u ",count[x]);
			}
			printf("\n");
		}
		//if(good == 5)
		if(good == 0) {


			printf("Good 0 count %u\n",index);

		}
//		printf("alloc %u max %u\n",count[5], bin_count[5]);
		allocation_dummy[allocation_id_index] = conf->dummies[index_to_allocate];
		//sets which allocation we added with reference to the bin and the index in the bin
		allocation_id[INDEX][allocation_id_index] = index;//count[good];//&(BID);
		allocation_id[BIN][allocation_id_index] = good;
		//how many allocation we have, decrement when we backtrack
		allocation_id_index++;
		goto lbl_continue;
	backtrack:
		;
		int backtrack;
		do {

			backtrack = 0;
			allocation_id_index--;
			allocation_dummy[allocation_id_index] = 0;
			int index = allocation_id[INDEX][allocation_id_index];			
			good = allocation_id[BIN][allocation_id_index];
			int index_to_dealloc = bin_index[good]+index;
			allocation_value -= conf->offer[index_to_dealloc];
//			id_to_allocate = bins[good].bids[index].id;
//allocation_id[allocation_id_index]->id;
//			good = allocation_id[allocation_id_index]->bin;
			//allocation_id[allocation_id_index] = NULL;
//			printf("dealloc good %u,id %u\n",good,id_to_allocate);
			/* for(x=good+1;x<goods;x++) { */
			/* 	count[x] =0; */
			/* } */
			
			printf("dealloc index %u bin %u aloc_index %u,bin count %u\n",
			       index_to_dealloc,
			       good,
			       allocation_id_index,
			       count[good]);
			for(x=0;x<words;x++) {
				allocation[x] &= ~(conf->allocation[index_to_dealloc*2+x]);
			}
			if(count[good] >= bin_count[good]) {			
				if(good == 0) {
					printf("New max of %.3lf \n",max);
					return;
				}

				count[good] =0;
				backtrack = 1;
			}
			good = 0;


		}while(backtrack);
//			printf("dealloc more\n");

//		goto backtrack;		
//		goto lbl_continue;
		;
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
		printf("id %u, offer %.3f, bin %u\n",conf->id[x],conf->offer[x],conf->bin[x]);

	}	
	
	/* unsigned int ** order;// = get_bin_order(bins,conf,bin_count); */
	/* printf("\nrow 0\t"); */
	/* for(x=0;x<conf->goods;x++) { */
	/* 	printf("%u\t",order[0][x]); */
	/* } */

	/* printf("\nrow 1\t"); */
	/* for(x=0;x<conf->goods;x++) { */
	/* 	printf("%u\t",order[1][x]); */
	/* } */
	/* printf("\n"); */

	/* for(x=0;x<conf->goods;x++) { */
	/* 	double score = conf->score[x]; */
	/* 	//printf("bin %d score %.3lf order %u count %u average %.3lf\n",x,score,order[0][x],bin_count[x],bins[x].average); */
	/* 	int y; */
	/* 	for(y =0; y < bin_count[x];y++) { */
	/* 		/\* unsigned int id = bins[x].bids[y].id; *\/ */
	/* 		/\* printf("Bin %d index %d id %u val %.03lf\n",x,y,id,bins[x].bids[y].offer); *\/ */
	/* 	} */
	/* } */
//	return;
	calc_best(conf);
	printf("Bye\n");
exit(EXIT_SUCCESS);
}
