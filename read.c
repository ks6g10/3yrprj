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
	unsigned int conf[];
};

struct bid2 {
	double offer;
	double average;
	unsigned int id;
	unsigned int conf[20];
};

struct linked_bids {
	struct bid2 * this;
	struct linked_bids * next;
}

struct allocation {
	struct linked_bids * root;
	unsigned int price;
	unsigned int conf[INTS];

};


/* unsigned int ints =0; */
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

double * calc_score(unsigned int total_goods,
		    unsigned int bids,
		    unsigned int ints,
		    struct bid * src_ptr) {
	double * score =  malloc(sizeof(score)*total_goods);
	unsigned int * score_total = malloc(sizeof(score_total)*total_goods);
	unsigned int * score_c = malloc(sizeof(score_c)*total_goods);
	int x,y;
	for(x = 0;x < total_goods;x++) {
		score[x] = 0.0;
		score_total[x] = 0;
		score_c[x] = 0;		
	}
	for(x=0; x < bids; x++) {
		for(y=0;y<ints;y++) {
			unsigned int conf = src_ptr->conf[y];
			while(conf) {
				unsigned int index = next_index(conf);
				conf &= ~(1 << index);
				score_total[(y*WORDSIZE + index)] += src_ptr->goods;
				score_c[(y*WORDSIZE + index)]++;
			}
		}
		src_ptr++;
	}
	for(x = 0;x < total_goods;x++) {
		score[x] =((double) score_c[x])/((double) score_t[x]);
	}
	free(score_c);
	free(score_total);
	return score;
}

struct bid2 ** assign_bids_to_bins(unsigned int total_goods,
				   unsigned int bids,
				   unsigned int ints,
				   struct bid * src_ptr) {

	struct bid2 **bin =(struct bid2 **) malloc(sizeof(struct bid2 **)*(total_goods));
	for(x = 0; x < total_goods;x++) {
		(bin)[x] = malloc((sizeof(struct bid2)*bin_count[x]));
		if(!(bin[x])) {
			printf("Could not allocate memory as line %d\n",__LINE__);
			exit(1);
		}
		printf("count = %u, bin %u\n",bin_count[x],x);
		bin_count[x]=0;
	}
	int y;

	struct bid2 * dst_ptr;// = *bin;
	struct bid * src_ptr = bids_ptr;//NULL;// = *bin;
	for(x = 0; x < bids;x++) {
		//*src_ptr = bids_ptr[x];
		unsigned int index = bin_count[src_ptr->bin];
		bin_count[src_ptr->bin]++;
		dst_ptr = bin[src_ptr->bin];
		(dst_ptr[index]).offer = src_ptr->offer;
		(dst_ptr[index]).id = src_ptr->id;
		(dst_ptr[index]).average =(double) src_ptr->offer/((double) src_ptr->goods);
		for(y = 0; y < ints;y++) {
			(dst_ptr[index]).conf[y] = src_ptr->conf[y];
		}
		src_ptr++;
		printf("id %u val %lf\n",(dst_ptr[index]).id,(dst_ptr[index]).offer);		
	}
	return bin;
}

int main(int argc, char *argv[])   {
	const char * s_goods = "goods";
	const char * s_bids = "bids";
	const char * s_dummy = "dummy";
	FILE * fp;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	unsigned int goods = 0;
	unsigned int bids = 0;
	unsigned int dummy = 0;
	unsigned int all = 0;
	unsigned int total_goods = 0;
	
	unsigned int *tmp = NULL;
	struct bid * bids_ptr;
	unsigned int bids_count =0;
	unsigned int * bin_count = NULL;	
	int x;
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	printf("hello\n");
	while ((read = getline(&line, &len, fp)) != -1 && !all) {
		if(line[0] == '%' || line[0] == '\n') {
			continue;
		}
		if(!all) {
			if(strncmp(line,s_goods,strlen(s_goods)) == 0) {      
				goods = atoi(line+strlen(s_goods)+1);
				printf("Number of goods %u\n",goods);				
			} else
			if(strncmp(line,s_bids,strlen(s_bids)) == 0) {
				bids = atoi(line+strlen(s_bids)+1);
				
				printf("Number of bids %u\n",bids);
			} else
			if(strncmp(line,s_dummy,strlen(s_dummy)) == 0) {
				dummy = atoi(line+strlen(s_dummy)+1);
				printf("Number of dummy %u\n",dummy);
				ints = 1+(total_goods-1)/32;
			}
			total_goods = goods + dummy;
			all = !!(goods && bids && dummy);
		} 

	}

	bin_count = calloc(sizeof(bin_count),total_goods);
	tmp = malloc(sizeof(tmp)*ints);			
	if(!(tmp) || !(bin_count)) {
		printf("Could not allocate memory as line %d\n",__LINE__);
		exit(1);
	}

	for(x =0;x < total_goods;x++) {
		bin_count[x] = 0;
	}

	bids_ptr = calloc((sizeof(bids_ptr)+(1+(total_goods-1)/32)*sizeof(int)),bids);
	if(!(bids_ptr)) {
		printf("Could not allocate memory as line %d\n",__LINE__);
		exit(1);
	}
	
	while ((read = getline(&line, &len, fp)) != -1) {
			if(isdigit(line[0])) {
				int x;
				char * head;
				char * tail;
				double value;
				unsigned int goods_count = 0;
				unsigned int id;
				unsigned int tmp2;
				int bool = 0;
				for(x = 0; x < ints; x++) {
					tmp[x] = 0;
				}
				head = tail = line;
				while(*head != '\t' && *head != '\0') {
					head++;
				}  
				id = strtol(tail,&head,10);
				tail = head;
				head++;
				while(*head != '\t' && *head != '\0') {
					head++;
				}
				value = strtod(tail,&head);
				tail = head;
				head++;

				bool = 0;
				while(*head != '#' && *head != '\0') {
					if(*head == '\t') {
						goods_count++;
						tmp2 = strtol(tail,&head,10);
						//sscanf(tail,"\t%u\t",tmp2);
						tmp[tmp2/32] |= (1 << tmp2);
						tail = head;
						if(!bool) {
							bin_count[tmp2]++;
							printf("count %u at %u\n",bin_count[tmp2],tmp2);
							bool = 1;
							bids_ptr[bids_count].bin = tmp2;
							
						}
					}
					head++;
				}
				bids_ptr[bids_count].id = id;
				bids_ptr[bids_count].offer = value;
				bids_ptr[bids_count].goods = goods_count;
				for(x=0;x<ints;x++) {
					bids_ptr[bids_count].conf[x] = tmp[x];
				}
				printf("ID %u Value %f %lu :\n", bids_ptr[bids_count].id,bids_ptr[bids_count].offer,bids_ptr[bids_count].conf[0]+((unsigned long)bids_ptr[bids_count].conf[1] << 32));
				bids_count++;
				goods_count =0;

			}

		
	}
	
	struct bid2 **bin =(struct bid2 **) malloc(sizeof(struct bid2 **)*(total_goods));
	for(x = 0; x < total_goods;x++) {
		(bin)[x] = malloc((sizeof(struct bid2)*bin_count[x]));
		if(!(bin[x])) {
			printf("Could not allocate memory as line %d\n",__LINE__);
			exit(1);
		}
		printf("count = %u, bin %u\n",bin_count[x],x);
		bin_count[x]=0;
	}
	int y;

	struct bid2 * dst_ptr;// = *bin;
	struct bid * src_ptr = bids_ptr;//NULL;// = *bin;
	
	for(x = 0; x < bids;x++) {
		//*src_ptr = bids_ptr[x];
		unsigned int index = bin_count[src_ptr->bin];
		bin_count[src_ptr->bin]++;
		dst_ptr = bin[src_ptr->bin];
		(dst_ptr[index]).offer = src_ptr->offer;
		(dst_ptr[index]).id = src_ptr->id;
		(dst_ptr[index]).average =(double) src_ptr->offer/((double) src_ptr->goods);
		for(y = 0; y < ints;y++) {
			(dst_ptr[index]).conf[y] = src_ptr->conf[y];
		}
		src_ptr++;
		printf("id %u val %lf\n",(dst_ptr[index]).id,(dst_ptr[index]).offer);		
	}
	
	double * score = calc_score(total_goods,bids,ints,bids_ptr);
	free(bids_ptr);
	//bids_ptr++;
//	free(bids_ptr);
	for(x = 0; x < bin_count[0];x++) {
		dst_ptr = bin[0];
		printf("id %u val %lf\n",(dst_ptr[x]).id,(dst_ptr[x]).offer);
	}
 
	/* free(tmp); */
	/* int x; */
	/* for(x=0;x<bids;x++) { */
	/* 	free(&bids_ptr[x]); */
	/* } */

	/* if (line) */
	/* 	free(line); */
	exit(EXIT_SUCCESS);
}
