/**
 em	 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>
#include <sys/types.h>
#include <assert.h>
#include <limits.h>
#include <float.h>


typedef uint32_t myint;
#define WORDSIZE (sizeof(myint)*8)
#define SIZE (1)
struct allocation {
	myint a[SIZE];
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

unsigned int ints = 0;
/* unsigned int * MASK = void; */

unsigned int next_index(unsigned int a_index) {
	return __builtin_ffs(a_index) - 1;
//	tmp  &= ~(1 << index);
//	upper_bound +=  v(bins[x],pi_conf,bin_counts[x]);

}

void print_binary(unsigned int * allocation, unsigned int goods) {
	int x;
	for (x = goods - 1; x >= 0; x--) {
		printf("%u", !!(*allocation & (1 << x)));
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
	unsigned int goods =0;
	unsigned int bids = 0;
	unsigned int dummy = 0;


	while ((read = getline(&line, &len, fp)) != -1 && !all) {
		if (line[0] == '%' || line[0] == '\n') {
			continue;
		}
		if (strncmp(line, s_goods, strlen(s_goods)) == 0) {
			goods = atoi(line + strlen(s_goods) + 1);
			printf("Number of goods %u\n", goods);
		} else if (strncmp(line, s_bids, strlen(s_bids)) == 0) {
			bids = atoi(line + strlen(s_bids) + 1);
			printf("Number of bids %u\n", bids);
		} else if (strncmp(line, s_dummy, strlen(s_dummy)) == 0) {
			dummy = atoi(line + strlen(s_dummy) + 1);
			got_dummy = 1;
			printf("Number of dummy %u\n", dummy);
//			ints = 1+(goods-1)/32;
		}
//			total_goods = goods + dummy;
		all = !!(goods && bids && got_dummy);
	}
	free(line);

	//if(goods <= 32) {
	struct configuration * ret = (struct configuration *) malloc(
					sizeof(struct configuration));


	//}

	ret->words = SIZE;
	return ret;
}

unsigned int * get_bincount(FILE * fp, struct configuration * conf,
		unsigned int * have_singelton) {
	unsigned int * bin_count = (unsigned int *) malloc(
			sizeof(int) * conf->goods);
	int x;
	for (x = 0; x < conf->goods; x++) {
		bin_count[x] = 0;
	}

	char * head = NULL;
	char * tail = NULL;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	printf("hello1\n");
	while ((read = getline(&line, &len, fp)) != -1) {
		if (!isdigit(line[0])) {
			continue;
		}
		head = tail = line;
		int tab_count = 0;
		head++;
		while (tab_count < 3) {
			if (*head == '\t') {
				tab_count++;
				if (tab_count <= 2) {
					tail = head;
					head++;
				}
			} else {
				head++;
			}
		}
		int which_bin = strtol(tail, &head, 10);
		head++;
		bin_count[which_bin]++;
		printf("abin %u count %u\n", which_bin, bin_count[which_bin]);
		int goods_count = 1;
		while (*head != '#' && *head != '\0') {
			if (*head == '\t') {
				goods_count++;
			}
			head++;
		}
		if (goods_count == 1) {
			have_singelton[which_bin] = 1;
		}
//		printf("Bin %d count %u\n",which_bin,bin_count[which_bin]);
	}
	free(line);

	return bin_count;
}

int compare_int(const void* p1, const void* p2) {
//	struct bid b1 =
	float i1 = ((struct bid*) p1)->average;
	float i2 = ((struct bid*) p2)->average;
	assert(((struct bid*)p1)->bin == ((struct bid*)p2)->bin);
	return i1 > i2 ? -1 : i1 < i2 ? 1 : 0;

}

struct bid * remove_from_list(struct bid * curr, struct bid * root) {

	if (curr->prev) { // if current node is not the first in the list
		curr->prev->next = curr->next; // then point prev to the next in the list
	} else {
		root = curr->next;
	}
	if (curr->next) { //if current node is not last in the list
		curr->next->prev = curr->prev; // then point the next node to the prev
		assert(curr != curr->next->prev);
	}

	return root;
}

int get_next_best_good(struct configuration * conf, struct bid * curr) {
	int x;
	unsigned int total_goods_count[conf->goods];
	unsigned int numbids_count[conf->goods];
	for (x = 0; x < conf->goods; x++) {
		total_goods_count[x] = numbids_count[x] = 0;
	}

	while (curr) {
		int goods_count = 0;
		for (x = 0; x < SIZE; x++) {
			goods_count += __builtin_popcount((unsigned int)curr->alloc.a[x]);
		}
		for (x = 0; x < conf->goods; x++) {
			int word_index = x / WORDSIZE;
			int bit_index = x % WORDSIZE;
			int result = !!(curr->alloc.a[word_index] & (1 << bit_index));
			total_goods_count[x] += result * goods_count;
			numbids_count[x] += result;

		}
		curr = curr->next;
	}
	int min_pos = 0;
	float min = FLT_MAX;
	for (x = 0; x < conf->goods; x++) {
		float score = 0.0f;
		if (numbids_count[x]) {
			double avg = ((double) total_goods_count[x])
					/ ((double) numbids_count[x]);
			score = ((double) numbids_count[x]) / avg;
			//printf("x %d score %.3f\n",x,score);
			if (score < min) {
				min = score;
				min_pos = x;
			}
		}
	}

	return min_pos;
}

void allocate_all_bids(FILE * fp, struct configuration * conf,
		unsigned int * have_singelton, unsigned int * bin_count) {
	conf->allocation = (struct allocation *) malloc(
			sizeof(struct allocation) * conf->bids);

	conf->id = (unsigned int *) malloc(sizeof(unsigned int) * conf->bids);
	conf->dummies = (unsigned int *) malloc(sizeof(unsigned int) * conf->bids);
	conf->bin = (unsigned int *) malloc(sizeof(unsigned int) * conf->bids);
	conf->offer = (float *) malloc(sizeof(float) * conf->bids);
	conf->average = (float *) malloc(sizeof(float) * conf->bids);
	conf->max_offset = (unsigned int *) malloc(
			sizeof(unsigned int) * conf->goods);

	conf->bin_count = bin_count;

	char * head = NULL;
	char * tail = NULL;
	char *line = NULL;
	unsigned long total_goods_count[conf->goods];

	unsigned int numbids_count[conf->goods];
	unsigned int goods[conf->goods];
	size_t len = 0;
	ssize_t read;
	int x;
	unsigned int bin_index[conf->goods];
	bin_index[0] = 0;
	total_goods_count[0] = numbids_count[0] = 0;
	for (x = 1; x < conf->goods; x++) {
		bin_index[x] = bin_count[x - 1] + bin_index[x - 1];
		total_goods_count[x] = numbids_count[x] = 0;
	}
	struct bid * tmp_bids = (struct bid *) malloc(
			sizeof(struct bid) * conf->bids);
	struct bid * root = &tmp_bids[0]; //malloc(sizeof(struct bid));
	struct bid * curr = root;
	curr->next = NULL;
	curr->prev = NULL;
	for (x = 1; x < conf->bids; x++) {
		curr->next = &tmp_bids[x]; //malloc(sizeof(struct bid));
		curr->next->prev = curr;
		curr = curr->next;
		curr->next = NULL;
	}

	curr = root;
	while ((read = getline(&line, &len, fp)) != -1) {
		if (!isdigit(line[0])) {
			continue;
		}
		head = tail = line;

		while (*head != '\t' && *head != '\0') {
			head++;
		}
		int id = strtol(tail, &head, 10);
		tail = head;
		head++;
		//get offer or value
		while (*head != '\t' && *head != '\0') {
			head++;
		}
		float offer = strtod(tail, &head);
		tail = head;
		head++;
		unsigned int goods_count = 0;
		unsigned int good = 0;
		unsigned int dummy = 0;
		unsigned int tmp_allocation[SIZE];
		for (x = 0; x < SIZE; x++) {
			tmp_allocation[x] = 0;
		}
		//reset the temporary goods array, used to determin the score
		goods_count = 0;
		for (x = 0; x < conf->goods; x++) {
			goods[x] = 0;
		}

		while (*head != '#' && *head != '\0') {
			if (*head == '\t') {
				good = strtol(tail, &head, 10);

				//sscanf(tail,"\t%u\t",tmp2);
				if (good < conf->goods) {
					tmp_allocation[(good / WORDSIZE)] += (1 << good);
					goods[goods_count] = good;
				} else {
					dummy = good;
				}
				tail = head;
				goods_count++;
			}
			head++;

		}
		if (dummy > 0) {
			goods_count--;
		}
		curr->average = (float) offer / (goods_count);
		for (x = 0; x < goods_count; x++) {
			total_goods_count[goods[x]] += goods_count;
			numbids_count[goods[x]]++;
		}
		curr->offer = offer;
		curr->bin = goods[0];
		curr->dummy = dummy;
		curr->id = id;
		for (x = 0; x < SIZE; x++) {
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
	int singleton_count = conf->bids - conf->singletons;
	for (x = 0; x < conf->goods; x++) {
		if (!have_singelton[x]) {
			int y;
			for (y = 1; y < SIZE; y++) {
				curr->alloc.a[x] = 0;
			}
			int word_index = x / WORDSIZE;
			int bit_index = x % WORDSIZE;
			curr->alloc.a[word_index] = (1 << bit_index);
			curr->offer = curr->average = 0.0f;
			curr->dummy = 0;
			curr->bin = x;
			curr->id = singleton_count;

			total_goods_count[x] += 1; //add one more to the score stat
			numbids_count[x] += 1; // also add one more to the number of bids to the score stat
			singleton_count++; // next singleton bid will have an consecutive bid id
			curr = curr->next;
		}
		double score = 0;
		double avg;
		if (numbids_count[x]) {
			printf("x %d total good count %lu, numbids_count %d\n", x,
					total_goods_count[x], numbids_count[x]);
			avg = ((double) total_goods_count[x]) / ((double) numbids_count[x]);
			score = ((double) numbids_count[x]) / avg;
		}
		if (score < min) {
			min = score;
			min_pos = x;
		}
	}
	unsigned int bid_to_bit[conf->goods];

	for (x = 0; x < conf->goods; x++) {
		bid_to_bit[x] = 0;
	}
	printf("min %.3f pos %d\n", min, min_pos);
	int bid_bit_count = -1;
	struct bid * new_root = NULL;
	struct bid * new_curr = NULL;
	int bid_count = 0;
	while (root) {
		bid_count = 0;
		bid_bit_count++;
		bid_to_bit[min_pos] = bid_bit_count;
		curr = root;
		while (curr) {

			int word_index = min_pos / WORDSIZE;
			int bit_index = min_pos % WORDSIZE;
			struct bid * next = curr->next;
			if (curr->alloc.a[word_index] & (1 << bit_index)) {
				curr->bin = bid_bit_count;
				if (!new_root) {
					root = remove_from_list(curr, root);
					new_root = curr;
					//curr = curr->next;
					new_curr = new_root;
					new_curr->prev = NULL;

				} else {
					root = remove_from_list(curr, root);
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
		conf->max_offset[bid_bit_count] = bid_count - 1;
		min_pos = get_next_best_good(conf, curr);
		printf("min pos %u\n", min_pos);
	}
	new_curr = new_root;
	while (new_curr) {
		struct allocation tmp;
		for (x = 0; x < SIZE; x++) {
			tmp.a[x] = 0;
		}
		for (x = 0; x < conf->goods; x++) {
			int bit_index = x % WORDSIZE;
			int word_index = x / WORDSIZE;
			if ((new_curr->alloc.a[word_index] & (1 << bit_index))) {
				tmp.a[word_index] |= (1 << bid_to_bit[x]);
//				printf("good %d translation %d\n",x,bid_to_bit[x]);
			}
		}
		for (x = 0; x < SIZE; x++) {
			new_curr->alloc.a[x] = tmp.a[x];
		}
		/* printf("id %u res %u\n",new_curr->id,new_curr->alloc.a[0]); */
		/* exit(0); */
		/* printf("%u bin %u\n",new_curr->id,bit_to_bid[new_curr->bin]); */
		new_curr = new_curr->next;
	}
	printf("total bids %u\n", conf->bids);

	for (x = 0; x < conf->goods; x++) {
		printf("%d %u\n", x, conf->bin_count[x]);
	}
	exit(0);

//	exit(0);
	int y;

	int bin_index2[conf->goods];
	bin_index2[0] = 0;

	for (x = 1; x < conf->goods; x++) {
		bin_index2[x] = bin_count[x - 1] + bin_index2[x - 1];
	}

	x = 0;

	struct bid * lhead, *ltail;
	ltail = new_root;
	while (ltail) {
		int good = ltail->bin;
		lhead = ltail->next;
		while (lhead && lhead->bin == good) {
			if (lhead->average > ltail->average) {
				if (lhead->prev == ltail) {

					if (ltail->prev)
						ltail->prev->next = lhead;
					else
						new_root = lhead;
					if (lhead->next)
						lhead->next->prev = ltail;
					lhead->prev = ltail->prev;
					ltail->next = lhead->next;
					lhead->next = ltail;
					ltail->prev = lhead;
					struct bid * tmp = lhead;
					lhead = ltail;
					ltail = tmp;
				} else {

					struct bid *ltailprev, *ltailnext;
					assert(ltail->next != lhead);
					assert(lhead->prev != ltail);
					assert(lhead->next != ltail);
					assert(ltail->prev != lhead);
					ltailnext = ltail->next;
					ltailprev = ltail->prev;
					if (ltail->prev)
						ltail->prev->next = lhead;
					else
						new_root = lhead;
					if (ltail->next)
						ltail->next->prev = lhead;
					if (lhead->prev)
						lhead->prev->next = ltail;
					if (lhead->next)
						lhead->next->prev = ltail;
					ltail->next = lhead->next;
					ltail->prev = lhead->prev;
					lhead->next = ltailnext;
					lhead->prev = ltailprev;
					struct bid * tmp = lhead;
					lhead = ltail;
					ltail = tmp;
				}
			}
			lhead = lhead->next;
		}
		ltail = ltail->next;
	}
	new_curr = new_root;

	while (new_curr) {

		int index = x;
		for (y = 0; y < SIZE; y++) {
			conf->allocation[index].a[y] = new_curr->alloc.a[y];
			assert(conf->allocation[index].a[y] == new_curr->alloc.a[y]);
		}
		conf->bin[index] = new_curr->bin;
		assert(
				conf->allocation[index].a[conf->bin[index]/WORDSIZE] & (1<< conf->bin[index]));
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

	return;
}

//wordsize in bits
#define WORD (32)
#define BIN (0)
#define INDEX (1)

void print_debug(struct configuration * conf, unsigned int (*allocation_count),
		unsigned int allocation_id_index, unsigned int low_order_good,
		unsigned int (*bin_count)) {
	return;
	printf("low order good %u allocation_id_index %u\n", low_order_good,
			allocation_id_index);
	int x;
	printf("allocation:\n");
	for (x = 0; x < conf->goods; x++) {
		printf("%u/%u\t", allocation_count[x], bin_count[x]);
	}
	printf("\n\n");

}

float h(struct configuration * conf, struct allocation * curr_allocation,
		int good, unsigned int * bin_index, unsigned int t) {
	int y, x;
	float value = 0.0f;
	assert(curr_allocation->a[0] == t);
	const int words = SIZE;
	for (y = 0; y < conf->goods; y++) {
		float partial_val = 0.0f;
		if ((curr_allocation->a[y / WORDSIZE] & (1 << y % WORDSIZE))== 0){
			int count = 0;
			assert((curr_allocation->a[y/WORDSIZE] & (1<< y%WORDSIZE)) == 0);
			unsigned int status;
			int max_offset = conf->max_offset[y];
			int index;
			for (; status != 0 && count <= max_offset; count++) {

				//if(count[good] > conf->max_offset[good]){ break;}
				status = 0;
				index = count + bin_index[y];

				assert(y == conf->bin[index]);
				if (conf->dummies[index]) {
					for (x = 0; x < conf->allocation_id_index; x++) {
						status |= (conf->allocation_dummy[x]
								== conf->dummies[index]);
					}
				}

				if (status)
					continue;
				for (x = 0; x < words; x++) {
					status |= (conf->allocation[index].a[x])
							& (curr_allocation->a[x]);
				}
				//if(conf->id[index] == 75) {printf("hello %u\n",status);}
			}
			if (!status) {
				partial_val = conf->average[index];
			}

		}
		value += partial_val;

	}
	printf("value %.3f\n", value);
	return value;
}

#define DEBUG (0)
#define H (1)
#define THREADS (1024)
#define WARPS (THREADS/32)
#define GOOD_I (0)
#define COUNT_I (1)
#define GOODS (32)
template<int uses_dummy,int single_word>
__global__ void calc_best2(unsigned int * max,
			   unsigned int * _max_index,
			   unsigned short * bin_index,
			   struct allocation * allocation,
			   unsigned int * offer,
			   unsigned int * dummies,
			   unsigned int * bin) {

	__shared__ struct allocation curr_allocation[WARPS];
	__shared__ unsigned short max_index[GOODS];
	__shared__ unsigned short count[WARPS][GOODS];
	__shared__ unsigned int value[WARPS];
	__shared__ short shared_vars[2];
	__shared__ unsigned int shared_max;
	__shared__ unsigned char allocation_id_index[WARPS];

	__shared__ unsigned short allocation_dummy[(WARPS)*uses_dummy][(GOODS-1)*uses_dummy];

	__shared__ unsigned short allocation_id[WARPS][GOODS-1];

	const char laneid = threadIdx.x  % 32;
	const char warpid = threadIdx.x / 32;
	//each warp reset the allocation
	if (laneid == 0) {
		curr_allocation[warpid] = allocation[blockIdx.x];
		value[warpid] = offer[blockIdx.x];
	}
	if (threadIdx.x < goods) {
		max_index[threadIdx.x] = _max_index[threadIdx.x];
	}
	if (laneid == 0) {
		allocation_id_index[warpid] = 0;
	}
	int x;
	for(x=laneid;x < (GOODS);x += 32) {
		
		count[warpid][x] = 0;
	}
	

	//start add the second bid ---------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------

	int good = 0;
	if (threadIdx.x == 0) {
		shared_max = *max;
		while ((curr_allocation[warpid].a[good / WORDSIZE]
			& (1 << good % WORDSIZE))&&
		       good < goods) {
			good++;
		}
		shared_vars[GOOD_I] = good;
		shared_vars[COUNT_I] = 32;
	}
	__syncthreads();
	good = shared_vars[GOOD_I];

	int binindex = bin_index[good];
	//if there is not enough bids in the start bin to continue
	//could later change such that it chooses the bin for which most bids exists
	if (warpid > max_index[good]) {
		return; //potentially fatal if i sync below
	}
	int index = binindex + warpid;
	int status = 1;

re_add_second_bid:
	;
	while (status) {
		status = 0;
		if (laneid < words) {
			status = curr_allocation[warpid].a[laneid]
				& allocation[index].a[laneid];
		}
		if (__any(status)) {
			if (laneid == 0) {
				index = atomicAdd((unsigned int *)&shared_vars[COUNT_I], 1);
			}
			index = __shfl(index, 0);
		}
		if (index > max_index[good]) {
			return;
		}
	}

	if (laneid < words) {
		curr_allocation[warpid].a[laneid] |= allocation[index].a[laneid];
	}

	if (laneid == 0) {
		value[warpid] += offer[index];
	}

	if (laneid == 0) {
		if(uses_dummy) {
			allocation_dummy[warpid][allocation_id_index[warpid]] = dummies[index];
		}
		//which bid we allocated
		allocation_id[warpid][allocation_id_index[warpid]] = index;
		allocation_id_index[warpid]++;

		if (value[warpid] > shared_max) {
			atomicMax(max, value[warpid]);
			atomicMax(&shared_max, value[warpid]);
		}
	}

	//end add the second bid -----------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------



	int allocate = 1; //could remove and just use goto
	int dealloc = 0;
	while (allocate || dealloc) {

		while (allocate) {
			// whilst the good have already been allocated
			status = 0;
			x = laneid;
			while (status == 0) {
				if (x < GOODS) { // if the good index is larger than the amount of goods
					if(single_word) {//template
						status = curr_allocation[warpid].a[0] & (1 << (x % WORDSIZE));
					} else {
						status = curr_allocation[warpid].a[x/WORDSIZE] & (1 << (x % WORDSIZE));
					}
				}
				status = __ballot((status == 0));

				if (status) {
					if (laneid == 0) {
						good = __ffs(status) - 1 + x; //which thread is the first to find empty good
					}
					good = __shfl(good, 0);
				} else {
					
					x += 32;
					if (x >= goods) {
						good = goods;
						break;
					}

				}
				
			}

			//if there are no more goods to allocate, exit allocation loop
			if (good >= goods) {
				//printf("dealloc full good\n");
				allocate = 0;
				dealloc = 1;
				break;
			}
			int max_offset = max_index[good];
			binindex = bin_index[good];

			index = 0;

			status = 0;
			while (status == 0) {

				index = count[warpid][good] + laneid;

				if(index >= max_offset) {
					status = 1;					
				}

				if (index < max_offset) {
					int y;
					if(uses_dummy) {
						unsigned int t_dummy = dummies[binindex + index];
						if (t_dummy) {
					
							for (y = 0; y < allocation_id_index[warpid]; x++) {
								status |= (allocation_dummy[warpid][y] == t_dummy);
							}
						}
					}
//					if(!__ballot(status == 0)) {
//						x += 32;
//						continue;
//					}
					for (y = 0; y < SIZE; y++) {
						status |= allocation[index + binindex].a[y]
							& (curr_allocation[warpid].a[y]);
					}
					status = __ballot((status == 0));
					
					if (status == 0) { //if no thread have found a compatible bid
						if (laneid == 0) {
							count[warpid][good] += 32;
						}						
					} else if (laneid == 0) {
						index = __ffs(status) - 1 + count[warpid][good];
						count[warpid][good] = index + 1;
					}

				}
				//if the count index is greater than the maximum offset
				if (count[warpid][good] > max_offset) {
					dealloc = 1;
					break;
				}

			}
			status = __shfl(status,0);
			
			//if count is greater than maximum offset and we did not find a suitable bid
			if (dealloc && !status) {
				if (laneid == 0) {
					count[warpid][good] = 0;
				}
				allocate = 0;
				break;
			}
			index = __shfl(index,0);
			
			//parallel --------------------------------------------------------------------------------------------------
			// add the goods
			if (laneid < SIZE) {
				curr_allocation[warpid].a[laneid] |=
					allocation[index].a[laneid];
			}

			//dummy bid for the bid we allocation
			if (laneid == 0) {
				value[warpid] += offer[index];
				if(uses_dummy) {
					allocation_dummy[warpid][allocation_id_index[warpid]] =
						dummies[index];
				}
				//which bid we allocated
				allocation_id[warpid][allocation_id_index[warpid]] = index;
				allocation_id_index[warpid]++;

				if (value[warpid] > shared_max) {										
					atomicMax(max, value[warpid]);
					atomicMax(&shared_max, value[warpid]);
					printf("new max %u\n", value[warpid]);
//				printf("bid id separated by tab\n");
//				for (x = 0; x < conf->allocation_id_index; x++) {
//					printf("%u\t", conf->id[conf->allocation_id[x]]);
//				}
//				printf("\n");
				}
			}
			//parallel --------------------------------------------------------------------------------------------------
		}

		while (dealloc) {
			dealloc = 0;
			allocate = 1;
			if (laneid == 0) {
				allocation_id_index[warpid] -= 1;
			}

			//printf("index %u\n",conf->allocation_id_index);
			int dealloc_index = allocation_id[warpid][allocation_id_index[warpid]];
			int dealloc_good = bin[dealloc_index];
			if (laneid == 0) {
				value[warpid] -= offer[dealloc_index];
			}

			if (laneid < SIZE) {
				curr_allocation[warpid].a[laneid] ^= allocation[dealloc_index].a[laneid];
			}



			if (count[warpid][dealloc_good] >= max_index[dealloc_good]) {
				//	printf("re-de-alloc good %u\n",dealloc_good);
				if (laneid == 0) {
					count[warpid][dealloc_good] = 0;
				}
				dealloc = 1;
				allocate = 0;
			}
			good = 0;
			if(allocation_id_index[warpid] == 0) {//if we deallocated the second bid
				dealloc = 0;
				allocate = 1;
				good = shared_vars[GOOD_I];
				if (laneid == 0) {
					index = atomicAdd((unsigned int *)&shared_vars[COUNT_I], 1);
				}
				index = __shfl(index, 0);				

				if(index > max_index[good]){
					return;
				}
				goto re_add_second_bid;								
			}			
		}

	}
}

/**
 * This macro checks return value of the CUDA runtime call and exits
 * the application if the call failed.
 */
#define CUDA_CHECK_RETURN(value) {											\
	cudaError_t _m_cudaStat = value;										\
	if (_m_cudaStat != cudaSuccess) {										\
		fprintf(stderr, "Error %s at line %d in file %s\n",					\
				cudaGetErrorString(_m_cudaStat), __LINE__, __FILE__);		\
		exit(1);															\
	} }

/**
 * Host function that prepares data array and passes it to the CUDA kernel.
 */
int setup_mem_and_run(struct configuration * conf) {
	void *d = NULL;
	int i;

	struct configuration transfer;
	unsigned short bin_index[conf->goods];
	bin_index[0] = 0;
	for(i = 1; i < conf->goods;i++) {
		bin_index[i] = bin_index[i-1] + conf->bin_count[i-1];
	}
	
	int bids = conf->bids;
	int goods = conf->goods;
	unsigned int * max;
	unsigned int _max = 0;

	unsigned short * d_bin_index;
	
	CUDA_CHECK_RETURN(cudaMalloc((void**) &d_bin_index, sizeof(unsigned short)));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &max, sizeof(unsigned int)));
	CUDA_CHECK_RETURN(
			cudaMalloc((void**) &transfer.allocation, sizeof(struct allocation) * bids));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &transfer.bin, sizeof(int) * bids));
	CUDA_CHECK_RETURN(cudaMalloc((void**) &transfer.id, sizeof(int) * bids));
	CUDA_CHECK_RETURN(
			cudaMalloc((void**) &transfer.dummies, sizeof(int) * bids));
	CUDA_CHECK_RETURN(
			cudaMalloc((void**) &transfer.offer, sizeof(float) * bids));
	CUDA_CHECK_RETURN(
			cudaMalloc((void**) &transfer.average, sizeof(float) * bids));
	CUDA_CHECK_RETURN(
			cudaMalloc((void**) &transfer.bin_count, sizeof(int) * goods));
	CUDA_CHECK_RETURN(
			cudaMalloc((void**) &transfer.max_offset, sizeof(int) * goods));
	//CUDA_CHECK_RETURN(cudaMalloc((void**) &d, sizeof(int) * WORK_SIZE));
	CUDA_CHECK_RETURN(
		cudaMemcpy(d_bin_index, bin_index, sizeof(unsigned short) * goods, cudaMemcpyHostToDevice));

	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.allocation, conf->allocation, sizeof(struct allocation) * bids, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.bin, conf->bin, sizeof(int) * bids, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.id, conf->id, sizeof(int) * bids, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.dummies, conf->dummies, sizeof(int) * bids, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.offer, conf->offer, sizeof(unsigned int) * bids, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.average, conf->average, sizeof(float) * bids, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.bin_count, conf->bin_count, sizeof(int) * goods, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(transfer.max_offset, conf->max_offset, sizeof(int) * goods, cudaMemcpyHostToDevice));
	CUDA_CHECK_RETURN(
			cudaMemcpy(max, &_max, sizeof(unsigned int)*1, cudaMemcpyHostToDevice));
	
	//CUDA_CHECK_RETURN(cudaMemcpy(d, idata, sizeof(int) * WORK_SIZE, cudaMemcpyHostToDevice));
	int blocks = conf->bin_count[0];
	calc_best2<0,1><<<blocks,THREADS,0>>>(max,transfer.max_offset,d_bin_index,transfer.allocation,transfer.offer,transfer.dummies,transfer.bin);
	
//	bitreverse<<<1, WORK_SIZE, WORK_SIZE * sizeof(int)>>>(d);

	CUDA_CHECK_RETURN(cudaThreadSynchronize());
	// Wait for the GPU launched work to complete
	CUDA_CHECK_RETURN(cudaGetLastError());
//	CUDA_CHECK_RETURN(
	//		cudaMemcpy(odata, d, sizeof(int) * WORK_SIZE, cudaMemcpyDeviceToHost));

//	for (i = 0; i < WORK_SIZE; i++)
//		printf("Input value: %u, device output: %u\n", idata[i], odata[i]);

	CUDA_CHECK_RETURN(
		cudaMemcpy(&_max , max, sizeof(unsigned int)*1, cudaMemcpyDeviceToHost));
	
	CUDA_CHECK_RETURN(cudaDeviceReset());

	return 0;
}

int main(int argc, char *argv[]) {
	int x;
	FILE * fp;
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	printf("hello\n");
	struct configuration * conf = get_configuration(fp);
	conf->singletons = 0;
	unsigned int * have_singleton = (unsigned int *) malloc(
			sizeof(int) * conf->goods);
	for (x = 0; x < conf->goods; x++) {
		have_singleton[x] = 0;
	}
	unsigned int * bin_count = get_bincount(fp, conf, have_singleton);

	for (x = 0; x < conf->goods; x++) {
		if (!have_singleton[x]) {
			conf->singletons++;
			conf->bids++;
			bin_count[x]++;
		}
		printf("bin %d count %u\n", x, bin_count[x]);
	}

	fclose(fp);
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	allocate_all_bids(fp, conf, have_singleton, bin_count);
	fclose(fp);
	for (x = 0; x < conf->bids; x++) {
		printf("x %d id %u, offer %.3f, bin %u, alloc %u, bin_count %u\n", x,
				conf->id[x], conf->offer[x], conf->bin[x],
				conf->allocation[x].a[0], conf->bin_count[conf->bin[x]]);

	}
	printf("words %u wordsize %lu\n", SIZE, WORDSIZE);
	free(have_singleton);
	setup_mem_and_run(conf);
	//calc_best2(conf);
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
