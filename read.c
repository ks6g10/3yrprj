#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define WORDSIZE (32)

struct bid {
	double offer;
	unsigned int id;
	unsigned int bin;
	unsigned int goods;
	unsigned int dummy;
	unsigned int conf[20];
};

struct bid_ptr {
	unsigned int index;
	double value;
	unsigned int dummy;
};

struct bid2 {
	double offer;
	double average;
	unsigned int id;
//	unsigned int index;
	unsigned int dummy;
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

struct configuration {
	unsigned int goods;
	unsigned int bids;
	unsigned int dummy;
	unsigned int words;
	unsigned int singletons;
};

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
	struct configuration * ret = calloc(sizeof(struct configuration),1);
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
	ret->words = 1+(ret->goods-1)/WORDSIZE;
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
		float value = strtod(tail,&head);
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
		printf("id %d bin %u count %u\n",id,bin_for_bid,tmp_count[bin_for_bid]);
	}
	int singleton_count = conf->bids-conf->singletons;
	int y;
	for(x =0;x< conf->goods;x++) {
		if(!have_singelton[x]) {
			bins[x].bids[tmp_count[x]].offer = 0;
			bins[x].bids[tmp_count[x]].id = singleton_count;
			bins[x].bids[tmp_count[x]].dummy = 0;
			bins[x].bids[tmp_count[x]].average = 0;
			for(y=0;y< conf->words;y++) {
				bid_conf[singleton_count][y] = 0;
			}
			int tmp_index = x/32;
			int tmp_bit = x%32;
			bid_conf[singleton_count][tmp_index] = (1<<tmp_bit);
			singleton_count++;
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

unsigned int * get_bin_order(struct bid_bin * bins,struct configuration * conf,unsigned int * bin_count) {
	double score[conf->goods];
	unsigned int * order = malloc(sizeof(int)*conf->goods);	
	int x,y;
	for(x=0;x<conf->goods;x++) {
		score[x] = bins[x].score;
		order[x] = conf->goods;
	}
	double max = 0;
	unsigned int max_index =0;
	
	for(x = 0;x<conf->goods;x++) {
		max_index = conf->goods;
		max = 0;
		for(y=0;y<conf->goods;y++){
			if(score[y] > max && bin_count[y] >0) {
				max_index = y;					
				max = score[y];
			}
			
		}
		if(max_index < conf->goods) {
			score[max_index] = 0;		
			order[max_index] = x;
		}
	}
	return order;
}

//wordsize in bits
#define WORD (32)
#define BIN (0)
#define INDEX (1)

void calc_best(struct configuration * conf,
	       unsigned int * bin_count,
	       unsigned int * order,
	       struct bid_bin * bins) {
	const unsigned int goods = conf->goods;
	unsigned int allocation_count[goods];
	unsigned int allocation[conf->words];
	//0 = bin, 1 index
	unsigned int allocation_id[goods][2];
	unsigned int allocation_id_index = 0;
	unsigned int allocation_dummy[goods];
	double value = 0;;
	double max = 0;
	int x;	
	int bin_index = goods;
	int order_count = 0;
	for(x = 0; x < goods; x++){
		allocation_count[x] = 0;
		allocation_dummy[x] = 0;
		if(order_count == order[x]) {
		     bin_index = x;
		}	     
	}
	

	for(x=0;x<conf->words;x++) {
	     allocation[x] = 0;
	}
	printf("Order count index %u alloc %u\n",bin_index,allocation[0]);
	if(!(bin_index < goods)) {
	     printf("errorrr\n");
	     return;
	}
	while(1) {
	     	for(x = 0; x < goods; x++){
		     if(order_count == order[x]) {
			  bin_index = x;
		     }	     
		}
		int word_index = bin_index/WORDSIZE;
		int word_bit = bin_index%WORDSIZE;
		if((allocation[word_index] & (1<<word_bit)) == 0 && bin_count[bin_index]) {
			int status =0;			
			struct bid2 *bid;
		lbl_retry:
			status = 0;
			bid = &(bins[bin_index].bids[allocation_count[bin_index]]);
			printf("bin %u id %u order count %u bin count %u\n",bin_index,bid->id,order_count,allocation_count[bin_index]);
			for(x=0;x<conf->words && !status;x++) {
				//	status |= bid->conf[x] & allocation[x];
			}
			for(x = 0; x < allocation_id_index;x++) {
				if(allocation_dummy[x]) {
					status |= (bid->dummy == allocation_dummy[x]);
				}
			}
			if(status) {
				allocation_count[bin_index]++;
				if(allocation_count[bin_index] < bin_count[bin_index]) {
					goto lbl_retry;
				} else {
					order_count++;
				}
			} else {
				
				for(x=0; x < conf->words;x++) {
					//	allocation[x] |= bid->conf[x];
				}
				value += bid->offer;
				allocation_dummy[allocation_id_index] = bid->dummy;
				allocation_id[allocation_id_index][BIN] = bin_index;
				allocation_id[allocation_id_index][INDEX] = allocation_count[bin_index];
				allocation_id_index++;
				order_count++;
				if(value > max) {
					max = value;
					printf("new max %0.2lf",max);
				}
			}

		} else {
//		     printf("increment order count\n");
		     order_count++;
		}
		if(order_count > goods) {
			//backtrack
		     return;
		}

	}


}

double pi_max;

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

	unsigned int (*bid_conf)[conf->bids] = malloc(sizeof(int[conf->bids][conf->words]));	
	
	fclose(fp);
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	struct bid_bin * bins = allocate_bid_bin(fp,bin_count,conf,bid_conf,have_singleton);
	fclose(fp);
	unsigned int * order = get_bin_order(bins,conf,bin_count);
	for(x=0;x<conf->goods;x++) {
		double score = bins[x].score;
		printf("bin %d score %.3lf order %u count %u average %.3lf\n",x,score,order[x],bin_count[x],bins[x].average);
		int y;
		for(y =0; y < bin_count[x];y++) {
			unsigned int id = bins[x].bids[y].id;
			printf("Bin %d index %d id %u\n",x,y,id);
		}
	}

	calc_best(conf, bin_count,order,  bins);
	printf("Bye\n");
exit(EXIT_SUCCESS);
}
