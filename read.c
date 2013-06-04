#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
 	


#include <assert.h>
#define WORDSIZE (sizeof(unsigned int)*8)
#include <limits.h>
#include <float.h>

//typedef unsigned int int_fast32_t;

struct allocation {
	unsigned int a[2];
};

struct bid {
	struct allocation alloc;
	unsigned int id;
	unsigned int bin;
	unsigned int dummy;
	float offer;
	float average;
	struct bid * next;
	struct bid * prev;
};

struct configuration {
	struct allocation * allocation;
	unsigned int * bin;
	unsigned int * id;
	// each bids dummy, if any else is set to 0
	unsigned int * dummies;
	unsigned int * bin_count;
	unsigned int * max_offset;
	float * offer;
	float * average;
	// the allocation dummy	
	unsigned int goods;
	unsigned int words;
	unsigned int bids;
	unsigned int dummy;
	unsigned int singletons;
	unsigned int * allocation_id;
	unsigned int * allocation_dummy;
	unsigned int allocation_id_index; 
	float allocation_value;
};


 unsigned int ints =0;
/* unsigned int * MASK = void; */

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
	int got_dummy = 0;
	unsigned int all = 0;
	
	struct configuration * ret = malloc(sizeof(struct configuration));
	ret->goods = 0;
	ret->bids = 0;
	ret->dummy = 0;
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
			got_dummy = 1;
			printf("Number of dummy %u\n",ret->dummy);
//			ints = 1+(goods-1)/32;
		}
//			total_goods = goods + dummy;
		all = !!(ret->goods && ret->bids && got_dummy);
	}
	free(line);
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
	
	char * head =NULL;
	char * tail = NULL;
	char * line = NULL;
	size_t len =0;
	ssize_t read;
	printf("hello1\n");
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
		printf("abin %u count %u\n",which_bin,bin_count[which_bin]);
		int goods_count =1;
		while(*head != '#' && *head != '\0') {
			if(*head == '\t') {				
				goods_count++;
			}
			head++;
		}
		if(goods_count == 1) {			
			have_singelton[which_bin] = 1;
		}
//		printf("Bin %d count %u\n",which_bin,bin_count[which_bin]);
	}
	free(line);

	return bin_count;       
}


int compare_int(const void* p1, const void* p2)
{
//	struct bid b1 = 
	float i1 = ((struct bid*)p1)->average;
	float i2 = ((struct bid*)p2)->average;
	assert(((struct bid*)p1)->bin == ((struct bid*)p2)->bin);
return i1 > i2 ? -1 : i1 < i2 ? 1 : 0;

}


struct bid * remove_from_list(struct bid * curr,struct bid * root) {

	if(curr->prev) { // if current node is not the first in the list
		curr->prev->next = curr->next;// then point prev to the next in the list
	} else {
		root = curr->next;
	}
	if(curr->next) {//if current node is not last in the list
		curr->next->prev = curr->prev;// then point the next node to the prev
		assert(curr != curr->next->prev);
	}
	
	return root;
}

int get_next_best_good(struct configuration * conf, struct bid * curr) {
	int x;
	unsigned int total_goods_count[conf->goods];
	unsigned int numbids_count[conf->goods];
	for(x=0;x<conf->goods;x++) {
		total_goods_count[x] = numbids_count[x] = 0;
	}

	while(curr) {
		int goods_count = 0;
		for(x=0;x<conf->words;x++) {
			goods_count += __builtin_popcount(curr->alloc.a[x]);			
		}		
		for(x=0;x<conf->goods;x++) {
			int word_index = x/WORDSIZE;
			int bit_index = x%WORDSIZE;
			int result = !!(curr->alloc.a[word_index] & (1 << bit_index));
			total_goods_count[x] += result*goods_count;
			numbids_count[x] += result;
			
		}
		curr = curr->next;
	}
	int min_pos = 0;
	float min = FLT_MAX;
	for(x=0;x<conf->goods;x++) {
		float score = 0.0f;
		if(numbids_count[x]) {	
			double avg = ((double)total_goods_count[x])/((double)numbids_count[x]);
			score = ((double)numbids_count[x])/avg;
			//printf("x %d score %.3f\n",x,score);
			if(score < min) {
				min = score;
				min_pos = x;
			}
		}
	}

	return min_pos;
}


 

void allocate_all_bids(FILE * fp,
		       struct configuration * conf,
		       unsigned int * have_singelton,
		       unsigned int * bin_count) {	
	conf->allocation = malloc(sizeof(struct allocation)*conf->bids);

	conf->id = malloc(sizeof(unsigned int)*conf->bids);
	conf->dummies = malloc(sizeof(unsigned int)*conf->bids);
	conf->bin = malloc(sizeof(unsigned int)*conf->bids);
	conf->offer = malloc(sizeof(float)*conf->bids);
	conf->average = malloc(sizeof(float)*conf->bids);
	conf->max_offset = malloc(sizeof(unsigned int)*conf->goods);
	
	conf->bin_count = bin_count;
	
	char * head = NULL;
	char * tail = NULL;
	char *line = NULL;
	unsigned long total_goods_count[conf->goods];

	unsigned int numbids_count[conf->goods];
	unsigned int goods[conf->goods];	
	size_t len =0;
	ssize_t read;	
	int x;
	unsigned int bin_index[conf->goods];
	bin_index[0] = 0;
	total_goods_count[0] = numbids_count[0] = 0;
	for(x = 1; x< conf->goods;x++) {
		bin_index[x] = bin_count[x-1]+bin_index[x-1];
		total_goods_count[x] = numbids_count[x] = 0;
	}
	struct bid * tmp_bids = malloc(sizeof(struct bid)*conf->bids);
	struct bid * root = &tmp_bids[0];//malloc(sizeof(struct bid));
	struct bid * curr = root;
	curr->next = NULL;
	curr->prev = NULL;
	for(x= 1;x<conf->bids;x++) {
		curr->next = &tmp_bids[x];//malloc(sizeof(struct bid));
		curr->next->prev = curr;
		curr = curr->next;
		curr->next = NULL;
	}

	curr = root;
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
					goods[goods_count] = good;	
				} else {
					dummy = good;
				}
				tail = head;
				goods_count++;
			}
			head++;
			
		}
		if(dummy > 0) {
			goods_count--;
		}
		curr->average  = (float) offer/(goods_count);
		for(x=0;x<goods_count;x++) {
			total_goods_count[goods[x]] += goods_count;
			numbids_count[goods[x]]++;
		}
		curr->offer = offer;
		curr->bin = goods[0];
		curr->dummy = dummy;
		curr->id = id;
		for(x = 0;x<conf->words;x++) {
			curr->alloc.a[x] = tmp_allocation[x];
		}
		curr = curr->next;
//		printf("id %d bin %u count %u value %.3lf\n",bid_count,bin_for_bid,tmp_count[bin_for_bid],0);
		bin_index[goods[0]]++;
		//printf("hello\n");
	}
	free(line);
	float min = FLT_MAX;
	int min_pos = 0;
	int singleton_count = conf->bids-conf->singletons;	
	for(x =0;x< conf->goods;x++) {
		if(!have_singelton[x]) {
			int y;
			for(y = 1;y<conf->words;y++) {
				curr->alloc.a[x] = 0;
			}
			int word_index = x/WORDSIZE;
			int bit_index = x % WORDSIZE;			
			curr->alloc.a[word_index] = (1<< bit_index);
			curr->offer = curr->average = 0.0f;
			curr->dummy = 0;
			curr->bin = x;
			curr->id = singleton_count;

			total_goods_count[x] +=1; //add one more to the score stat
			numbids_count[x] +=1; // also add one more to the number of bids to the score stat
			singleton_count++; // next singleton bid will have an consecutive bid id
			curr = curr->next;
		}
		double score =0;
		double avg;
		if(numbids_count[x]) {	
			printf("x %d total good count %lu, numbids_count %d\n",x,
			       total_goods_count[x],numbids_count[x]);
		avg = ((double)total_goods_count[x])/((double)numbids_count[x]);
		score = ((double)numbids_count[x])/avg;
		}
		if(score < min) {
			min = score;
			min_pos = x;
		}
	}
	unsigned int bid_to_bit[conf->goods];

	for(x=0;x<conf->goods;x++) {
		bid_to_bit[x] = 0;	       
	}
	printf("min %.3f pos %d\n",min,min_pos);
	int bid_bit_count = -1;	
	struct bid * new_root = NULL;
	struct bid * new_curr = NULL;
	int bid_count = 0;
	while(root) {
		bid_count = 0;
		bid_bit_count++;
		bid_to_bit[min_pos] = bid_bit_count;
		curr = root;
	while(curr) {
		
		int word_index = min_pos/WORDSIZE;
		int bit_index = min_pos%WORDSIZE;
		struct bid * next = curr->next;
		if(curr->alloc.a[word_index] & (1 << bit_index)) {
			curr->bin = bid_bit_count;
			if(!new_root) {
				root = remove_from_list(curr,root);
				new_root = curr;
				//curr = curr->next;
				new_curr = new_root;
				new_curr->prev = NULL;				
				
			} else {
				root = remove_from_list(curr,root);
				new_curr->next = curr;
				new_curr->next->prev = new_curr;
				//curr = curr->next;
				
				new_curr = new_curr->next;
			}			
			new_curr->next = NULL;
			bid_count++;
		}		
		curr = next;

	}
	curr = root;
	conf->bin_count[bid_bit_count] = bid_count;
	conf->max_offset[bid_bit_count] = bid_count -1;
	min_pos = get_next_best_good(conf,curr);
	printf("min pos %u\n",min_pos);	
	}
	new_curr = new_root;
	while(new_curr) {
		struct allocation tmp;
		for(x=0;x<conf->words;x++) {
			tmp.a[x] = 0;
		}
		for(x=0;x<conf->goods;x++) {
			int bit_index = x % WORDSIZE;
			int word_index = x / WORDSIZE;
			if((new_curr->alloc.a[word_index] & (1 << bit_index))) {
				tmp.a[word_index] |= (1 << bid_to_bit[x]);
//				printf("good %d translation %d\n",x,bid_to_bit[x]);
			}
		}
		for(x=0;x<conf->words;x++) {
			new_curr->alloc.a[x] = tmp.a[x];
		}
		/* printf("id %u res %u\n",new_curr->id,new_curr->alloc.a[0]); */
		/* exit(0); */
		/* printf("%u bin %u\n",new_curr->id,bit_to_bid[new_curr->bin]); */
		new_curr = new_curr->next;
	}
	printf("total bids %u",conf->bids); 

	for(x=0;x<conf->goods;x++) {
		printf("x %d count %u\n",x,conf->bin_count[x]);
	}
	

//	exit(0);
	int y;

	int bin_index2[conf->goods];
	bin_index2[0] = 0;

	for(x = 1; x< conf->goods;x++) {
		bin_index2[x] = bin_count[x-1]+bin_index2[x-1];
	}

	x = 0;


	struct bid * lhead, * ltail;
	ltail = new_root;
	while(ltail) {
		int good = ltail->bin;
		lhead = ltail->next;
		while(lhead && lhead->bin == good) {
			if(lhead->average > ltail->average) {
				if(lhead->prev == ltail) {

					if(ltail->prev)
						ltail->prev->next = lhead;
					else
						new_root = lhead;
					if(lhead->next)
						lhead->next->prev = ltail;
					lhead->prev = ltail->prev;
					ltail->next = lhead->next;
					lhead->next = ltail;
					ltail->prev = lhead;
					struct bid * tmp = lhead;
					lhead= ltail;
					ltail = tmp;
				} else {
					
				struct bid *ltailprev,*ltailnext;
				assert(ltail->next != lhead);
				assert(lhead->prev != ltail);
				assert(lhead->next != ltail);
				assert(ltail->prev != lhead);
				ltailnext = ltail->next;
				ltailprev = ltail->prev;
				if(ltail->prev)
					ltail->prev->next = lhead;
				else
					new_root = lhead;
				if(ltail->next)
					ltail->next->prev = lhead;
				if(lhead->prev)
					lhead->prev->next = ltail;
				if(lhead->next)
					lhead->next->prev = ltail;
				ltail->next = lhead->next;
				ltail->prev = lhead->prev;
				lhead->next = ltailnext;
				lhead->prev = ltailprev;
				struct bid * tmp = lhead;
				lhead= ltail;
				ltail = tmp;
				}
			}
			lhead = lhead->next;
		}
		ltail = ltail->next;
	}
	new_curr = new_root;

	while(new_curr) {
	
		int index =x;
		for(y = 0; y< conf->words;y++) {
			conf->allocation[index].a[y] = new_curr->alloc.a[y];
			assert(conf->allocation[index].a[y] == new_curr->alloc.a[y]);
		}
		conf->bin[index] = new_curr->bin;
		assert(conf->allocation[index].a[conf->bin[index]/WORDSIZE] & (1<< conf->bin[index]));
		assert(x >= bin_index2[conf->bin[index]]);
		assert(x < bin_count[conf->bin[index]]+bin_index2[conf->bin[index]]);

		conf->offer[index] = new_curr->offer;
		conf->dummies[index] = new_curr->dummy;
		conf->id[index] = new_curr->id;
		conf->average[index] = new_curr->average;
		new_curr = new_curr->next;
		x++;		
	}
	free(tmp_bids);

	return ;	
}


//wordsize in bits
#define WORD (32)
#define BIN (0)
#define INDEX (1)

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

float h(struct configuration * conf,struct allocation  *  curr_allocation,int good,unsigned int * bin_index, unsigned int t) {
	int y,x;
	float value = 0.0f;
	assert(curr_allocation->a[0] == t);
	const int words = conf->words;
	for(y=0;y<conf->goods;y++){
		float partial_val = 0.0f;
		if((curr_allocation->a[y/WORDSIZE] & (1<< y%WORDSIZE)) == 0) {
			int count = 0;
			assert((curr_allocation->a[y/WORDSIZE] & (1<< y%WORDSIZE)) == 0);
			unsigned int status;
			int max_offset = conf->max_offset[y];
			for(status = 1;status != 0  && count <= max_offset ;count++) {

				status = 0;
				
				int index = count+bin_index[y]; 

				if(conf->dummies[index]) {
					for(x=0;x<conf->allocation_id_index;x++) {
						status |= (conf->allocation_dummy[x] == conf->dummies[index]);
					}
				}
				
				if(status) continue;
				for(x=0;x<words;x++) {				
					status |= (conf->allocation[index].a[x]) & (curr_allocation->a[x]); 
//					printf("status %u\n",status);
				}
				if(!status) {				
					partial_val = conf->average[index];
					break;
				}
				//			assert(count[good] <= conf->max_offset[good]);
			}


		}
		value += partial_val;
	}
	return value;
}


void calc_best2(struct configuration * conf) {
	struct allocation * curr_allocation = malloc(sizeof(struct allocation));	
	const unsigned int goods = conf->goods;
	const unsigned int words = conf->words;
	int x;
	float val = 0.0f;
	float max = 0.0f;
	unsigned int count[goods];	
	unsigned int * bin_count = conf->bin_count;	
	conf->allocation_id = malloc(sizeof(int) *conf->goods);
	conf->allocation_id_index = 0;
	conf->allocation_dummy = malloc(sizeof(int) *conf->goods);
	
	for(x=0;x<words;x++) {
		curr_allocation->a[x] = 0;
	}
	for(x= 0; x < goods; x++) {
		count[x] = 0;
	}
	int allocate = 1;
	int dealloc = 0;
	int good = 0;
	unsigned int bin_index[conf->goods];
	bin_index[0] = 0;
	for(x = 1; x< conf->goods;x++) {
		bin_index[x] = bin_count[x-1]+bin_index[x-1];
		assert(conf->bin[bin_index[x]] == x);
	}
	printf("count %d max %u a\n",count[good],bin_count[good]);
while(allocate || dealloc) {

	while(allocate) {
		// whilst the good have already been allocated
		while( (curr_allocation->a[good/WORDSIZE] & (1<< good%WORDSIZE) )  && 
		       good < goods) 
		{
			good++;
		}
		//if there are no more goods to allocate, exit allocation loop
		if(good >= goods) {
			allocate = 0;
			dealloc = 1;
			break;
		}
		int max_offset = conf->max_offset[good];
		int status = 1;
		//printf("%u\n",(conf->allocation[98].a[0]) & (curr_allocation->a[0]));
		// test if bid is allocatable, should not fail as we have have dummy bids
		for(;status != 0  && count[good] <= max_offset ;count[good]++) {

			//if(count[good] > conf->max_offset[good]){ break;}
			status = 0;
			assert((curr_allocation->a[good/WORDSIZE] & (1<< good%WORDSIZE)) == 0);
			int index = count[good]+bin_index[good]; 			
			assert(good == conf->bin[index]);
			if(conf->dummies[index]) {
				for(x=0;x<conf->allocation_id_index;x++) {
					status |= (conf->allocation_dummy[x] == conf->dummies[index]);
				}
			}

			if(status) continue;
			for(x=0;x<words;x++) {				
				status |= (conf->allocation[index].a[x]) & (curr_allocation->a[x]); 
			}
//			if(conf->id[index] == 75) {printf("hello %u\n",status);}
			assert(count[good] <= conf->max_offset[good]);
		}		
		if(count[good] > max_offset){
			assert(conf->bin[(count[good]+bin_index[good])] != good);
			count[good] = 0;
			allocate = 0;
			dealloc = 1;
			//	printf("dealloc max good %u\n",good);
			break;
		}
		//as we incremented it one more, reduce by one
		int which_bid_index = (count[good]-1)+bin_index[good]; 
		assert(which_bid_index < conf->bids);
		for(x=0;x<words;x++) {
			curr_allocation->a[x] |= conf->allocation[which_bid_index].a[x];
		}

		val += conf->offer[which_bid_index];
		//dummy bid for the bid we allocation

		conf->allocation_dummy[conf->allocation_id_index] = conf->dummies[which_bid_index];
		//which bid we allocated
		conf->allocation_id[conf->allocation_id_index] = which_bid_index;
		conf->allocation_id_index++;

		if(val > max) {
			max = val;
			printf("new max %.3f\n",max);
			printf("bid id separated by tab\n");
			for(x=0;x<conf->allocation_id_index;x++) {
				printf("%u\t",conf->id[conf->allocation_id[x]]);
			}
			printf("\n");
		}
//		printf("befor alloc %u\n",curr_allocation->a[0]);
/* 		float tmp = h(conf,curr_allocation,good,bin_index,curr_allocation->a[0]); */
/* 		if(val+tmp <= max) { */
/* 			printf("val %.3f tmp %.3f + %.3f max %.3f\n",val,tmp,val+tmp,max); */
/* //			count[good] = 0; */
/* 			allocate = 0; */
/* 			dealloc = 1; */
/* //				printf("dealloc max good %u\n",good); */
/* 			break; */
/* 		} */
	}

	while(dealloc) {
		dealloc = 0;
		allocate = 1;
		conf->allocation_id_index--;
		//printf("index %u\n",conf->allocation_id_index);
		int dealloc_index = conf->allocation_id[conf->allocation_id_index];
		val -= conf->offer[dealloc_index];
		for(x=0;x<words;x++) {
			curr_allocation->a[x] &= ~(conf->allocation[dealloc_index].a[x]);
		}
		int dealloc_good = conf->bin[dealloc_index];
		if(dealloc_good == 0) {
			printf("good 0 count %u\n",count[0]);
		}
		if(count[dealloc_good] >= conf->max_offset[dealloc_good]) {
			//	printf("re-de-alloc good %u\n",dealloc_good);
			count[dealloc_good] = 0;
			dealloc = 1;
			allocate = 0;
			if(dealloc_good == 0) {
				//incase it is the first good, meaning we have searched the tree
				dealloc = 0;
				
			}
		}
		if(conf->allocation_id_index) {
			int prev_dealloc_index = conf->allocation_id[conf->allocation_id_index-1];
			good = conf->bin[prev_dealloc_index];
		} else {
			good = 0;
		}
		/* printf("de-allocating good %u\n",dealloc_good); */
		/* int u; */
		/* for(u=0;u<conf->goods;u++) { */
		/* 	if(u == dealloc_good) { */
		/* 		printf("%u/%u<-\t",count[u],conf->bin_count[u]); */
		/* 	} else { */
		/* 		printf("%u/%u\t",count[u],conf->bin_count[u]); */
		/* 	} */
		/* } */
		/* printf("\n"); */
		//getchar();

	}

}
printf("max %.3f\n",max);
free(curr_allocation);
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
		printf("x %d id %u, offer %.3f, bin %u, alloc %u, bin_count %u\n",
		       x,
		       conf->id[x],
		       conf->offer[x],
		       conf->bin[x],
		       conf->allocation[x].a[0],
			conf->bin_count[conf->bin[x]]);

	}	
	printf("words %u wordsize %lu\n",conf->words,WORDSIZE);
	free(have_singleton);
	calc_best2(conf);
	free(conf->allocation);
	free(conf->bin);
	free(conf->id);
	free(conf->dummies);
	free(conf->bin_count);
	free(conf->max_offset);
	free(conf->offer);
	free(conf->average);
	free(conf->allocation_id);
	free(conf->allocation_dummy);
	free(conf);
	printf("Bye\n");
exit(EXIT_SUCCESS);
}
