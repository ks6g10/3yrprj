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

struct bid2 {
	double offer;
	double average;
	unsigned int id;
	unsigned int dummy;
	unsigned int conf[20];
};

struct bid_bin {
	unsigned int size;
	unsigned int good;
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
};

struct linked_bids {
	struct bid2 * this;
	struct linked_bids * next;
};

struct allocation {
	struct linked_bids * root;
	unsigned int price;
	unsigned int conf[20];

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

/* double * calc_score(unsigned int total_goods, */
/* 		    unsigned int bids, */
/* 		    unsigned int ints, */
/* 		    struct bid * src_ptr) { */
/* 	double * score =  malloc(sizeof(score)*total_goods); */
/* 	unsigned int * score_total = malloc(sizeof(score_total)*total_goods); */
/* 	unsigned int * score_c = malloc(sizeof(score_c)*total_goods); */
/* 	int x,y; */
/* 	for(x = 0;x < total_goods;x++) { */
/* 		score[x] = 0.0; */
/* 		score_total[x] = 0; */
/* 		score_c[x] = 0;		 */
/* 	} */
/* 	for(x=0; x < bids; x++) { */
/* 		for(y=0;y<ints;y++) { */
/* 			unsigned int conf = src_ptr->conf[y]; */
/* 			while(conf) { */
/* 				unsigned int index = next_index(conf); */
/* 				conf &= ~(1 << index); */
/* 				score_total[(y*WORDSIZE + index)] += src_ptr->goods; */
/* 				score_c[(y*WORDSIZE + index)]++; */
/* 			} */
/* 		} */
/* 		src_ptr++; */
/* 	} */
/* 	for(x = 0;x < total_goods;x++) { */
/* 		score[x] =((double) score_c[x])/((double) score_total[x]); */
/* 	} */
/* //	free(score_c); */
/* //	free(score_total); */
/* 	return score; */
/* } */


/* struct bid2 ** assign_bids_to_bins(unsigned int total_goods, */
/* 				   unsigned int bids,				    */
/* 				   unsigned int* bin_count, */
/* 				   struct bid bids_ptr[]) { */
/* 	int x,y; */
/* 	printf("total_goods %u %d\n",total_goods,sizeof(struct bid2 *)); */
/* 	unsigned int tmp_count[total_goods]; */
/* 	struct bid2 **bin_ret =/\*NULL;//\*\/  (struct bid2**) malloc(sizeof(struct bid2 *)*(total_goods)); */
	
/* 	if(!bin_ret) { */
/* 			printf("Could not allocate memory at line %d\n",__LINE__); */
/* 			exit(1); */
/* 	} */
/* 	printf("total_goods %u\n",total_goods); */
	
/* 	printf("bincount %u\n",bin_count[1]); */
/* 	for(x = 0; x < total_goods;x++) { */
/* 	       struct bid2 * bin_tmp = calloc(sizeof(struct bid2),bin_count[x]); */
/* 	       bin_ret[x] = bin_tmp ; */
/* 		if(!bin_ret[x]) { */
/* 			printf("Could not allocate memory at line %d\n",__LINE__); */
/* 			exit(1); */
/* 		} */

	      
/* 		printf("count = %u, bin %u\n",bids_ptr[x].id,x); */
/* 		tmp_count[x]=0; */
/* 	} */

/* 	struct bid2 * dst_ptr;// = *bin; */
/* 	struct bid * src_ptr = bids_ptr;//NULL;// = *bin; */
/* 	for(x = 0; x < bids;x++) { */
/* 		//\*src_ptr = bids_ptr[x]; */
/* 		const unsigned int which_bin = bids_ptr->bin; */
/* 		const unsigned int index = tmp_count[which_bin]; */
/* 		printf("bin %u at %u\n",which_bin,index); */
/* 		dst_ptr = bin_ret[which_bin]; */
/* 		tmp_count[which_bin]++; */
/* 		(dst_ptr+index)->offer = bids_ptr->offer; */
/* 		dst_ptr[index].id = bids_ptr->id; */
/* //		(dst_ptr[index]).id = src_ptr->id; */
/* //		(dst_ptr[index]).average =(double) src_ptr->offer/((double) src_ptr->goods); */
/* 		bids_ptr++; */
/* 		for(y = 0; y < ints;y++) { */
/* //			(dst_ptr[index]).conf[y] = src_ptr->conf[y]; */
/* 		} */
/* //		src_ptr++; */
/* 		//	printf("id %u val %lf %u\n",(dst_ptr[index]).id,(bids_ptr[index]).offer,bids); */
/* 	} */
/* 	printf("bincount %u\n",bin_count[1]); */
/* 	return bin_ret; */
/* } */

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
	return ret;
}

unsigned int * get_bincount(FILE * fp,struct configuration * conf) {
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
		bin_count[which_bin]++;
//		printf("Bin %d count %u\n",which_bin,bin_count[which_bin]);
	}
	return bin_count;       
}

struct bid_bin * allocate_bid_bin(FILE * fp,
				  unsigned int * bin_count,				  
				  struct configuration * conf) {
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
		bins[x].bids = malloc(sizeof(struct bid2)*bin_count[x]);
		total_goods_count[x] = 0;
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
		unsigned int tmp_conf[20] = {0};
		for(x = 0;x<20;x++) {
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
				if(goods_count == 1) {					
					bin_for_bid = good;
				}
				goods_count++;
			}
			head++;
		}
		
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].offer = value;
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].id = id;
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].dummy = dummy_good;
		bins[bin_for_bid].bids[tmp_count[bin_for_bid]].average = value/((double)goods_count);
		for(x = 0;x < 20; x++) {
			bins[bin_for_bid].bids[tmp_count[bin_for_bid]].conf[x] = tmp_conf[x];
		}
		for(x=0;x<goods_count;x++) {
			total_goods_count[goods[x]] += goods_count;
			numbids_count[goods[x]]++;
		}
		tmp_count[bin_for_bid] +=1;
		printf("id %d bin %u count %u\n",id,bin_for_bid,tmp_count[bin_for_bid]);
	}
	for(x =0;x< conf->goods;x++) {
		double score =0;
		double avg;
		if(numbids_count[x]) {
		avg = ((double)total_goods_count[x])/((double)numbids_count[x]);
		score = ((double)numbids_count[x])/avg;
		}
		bins[x].score = score;
		printf("bin %d score %.3lf\n",x,score);
	}
	return bins;
}


int main(int argc, char *argv[])   {

	FILE * fp;
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	printf("hello\n");
	struct configuration * conf = get_configuration(fp);	
	unsigned int * bin_count = get_bincount(fp,conf);
	int x;
	for(x =0; x < conf->goods;x++) {
		printf("bin %d count %u\n",x,bin_count[x]);
	}
	fclose(fp);
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	struct bid_bin * bins = allocate_bid_bin(fp,bin_count,conf);
/* 	bin_count = malloc(sizeof(bin_count)*total_goods); */
/* 	tmp = malloc(sizeof(tmp)*ints);			 */
/* 	if(!(tmp) || !(bin_count)) { */
/* 		printf("Could not allocate memory at line %d\n",__LINE__); */
/* 		exit(1); */
/* 	} */

/* 	for(x =0;x < total_goods;x++) { */
/* 		bin_count[x] = 0; */
/* 	} */

/* 	bids_ptr = malloc((sizeof(bids_ptr))*bids); */
/* 	if(!(bids_ptr)) { */
/* 		printf("Could not allocate memory at line %d\n",__LINE__); */
/* 		exit(1); */
/* 	} */
/* 	struct bid * next_ptr = bids_ptr; */
/* 	while ((read = getline(&line, &len, fp)) != -1) { */
/* 			if(isdigit(line[0])) { */
/* 				int x; */
/* 				char * head; */
/* 				char * tail; */
/* 				double value; */
/* 				unsigned int goods_count = 0; */
/* 				unsigned int id; */
/* 				unsigned int tmp2; */
/* 				int bool = 0; */
/* 				for(x = 0; x < ints; x++) { */
/* 					tmp[x] = 0; */
/* 				} */
/* 				head = tail = line; */
/* 				//get id */
/* 				while(*head != '\t' && *head != '\0') { */
/* 					head++; */
/* 				}   */
/* 				id = strtol(tail,&head,10); */
/* 				tail = head; */
/* 				head++; */
/* 				//get offer or value */
/* 				while(*head != '\t' && *head != '\0') { */
/* 					head++; */
/* 				} */
/* 				value = strtod(tail,&head); */
/* 				tail = head; */
/* 				head++; */
/* 				unsigned int dummy_good = 0; */
/* 				bool = 0; */
/* 				while(*head != '#' && *head != '\0') { */
/* 					if(*head == '\t') { */
/* 						goods_count++; */
/* 						tmp2 = strtol(tail,&head,10); */
/* 						//sscanf(tail,"\t%u\t",tmp2); */
/* 						if(tmp2 < goods) { */
/* 							tmp[tmp2/32] |= (1 << tmp2); */
/* 						} else { */
/* 							dummy_good = tmp2; */
/* 						} */
/* 						tail = head; */
/* 						if(!bool) { */
/* 							bin_count[tmp2]++; */
/* 							printf("count %u at %u\n",bin_count[tmp2],tmp2); */
/* 							bool = 1; */
/* 							next_ptr->bin = tmp2; */
							
/* 						} */
/* 					} */
/* 					head++; */
/* 				} */
/* 				next_ptr->id = id; */
/* 				next_ptr->offer = value; */
/* 				next_ptr->goods = goods_count; */
/* 				next_ptr->dummy = dummy_good; */
/* 				for(x=0;x<ints;x++) { */
/* 					next_ptr->conf[x] = tmp[x]; */
/* 				} */
/* 				printf("ID %u Value %f %lu :\n", next_ptr->id,next_ptr->offer,next_ptr->conf[0]+((unsigned long)next_ptr->conf[1] << 32)); */
/* 				bids_count++; */
/* 				next_ptr++; */
/* 				goods_count =0; */

/* 			} */
/* 	} */
/* 	fclose(fp); */
/* 	struct root_bid * bin = malloc(sizeof(struct root_bid)); */
/* 	bin->goods = goods; */
/* 	bin->bins = malloc(sizeof(struct bid_bin)*bin->goods); */
/* 	for(x= 0;x< goods;x++){ */
/* 		bin->bins[x].size=bin_count[x]; */
/* 		bin->bins[x].good=x; */
/* 		bin->bins[x].bids = malloc(sizeof(struct bid2)*bin->bins[x].size); */
/* 		bin->bins[x].size = 0; */
/* 	 	printf("bin %d count %u\n",x,bin_count[x]); */
/* 	} */
/* 	next_ptr = bids_ptr; */
/* 	for(x=0;x<bids;x++) {		 */
/* 		unsigned int bidbin = next_ptr->bin; */
/* //		struct bid_bin tmp_bin = ; */
/* 		unsigned int index = bin->bins[bidbin].size; */
/* 		bin->bins[bidbin].size +=1; */
/* //		tmp_bin.size +=1; */
/* 		printf("bin %u count %u max %u\n",bidbin,index,bin_count[bidbin]); */
/* 		struct bid2 * tmp_bid = &(bin->bins[bidbin].bids[index]); */
/* 		tmp_bid->id = next_ptr->id; */
/* 		tmp_bid->offer = next_ptr->offer; */
/* 		tmp_bid->average =(double) next_ptr->offer/((double)next_ptr->goods); */
/* 		tmp_bid->dummy = next_ptr->dummy; */
/* 		//memcpy(&(tmp_bid->conf),next_ptr->conf,sizeof(unsigned int)*20); */
/* //		memcpy((*bin).bins[bidbin].bids[index],&tmp,sizeof(struct bid2)); */
		
/* 		next_ptr++; */
/* 	} */
/* 	printf("Hello\n"); */
/* 	struct bid_bin tmp_bin = bin->bins[x]; */
/* 	for(x =0;x < bin_count[0];x++){ */
/* 		printf("id %u\n",(tmp_bin.bids[x]).id); */
/* 	} */


//	struct bid2 **bin = assign_bids_to_bins(total_goods,bids_count,bin_count,bids_ptr);	
	printf("Bye\n");
	//double * score = calc_score(total_goods,bids,ints,bids_ptr);
//	free(bids_ptr);
	//bids_ptr++;
//	free(bids_ptr);



 
	/* free(tmp); */
	/* int x; */
	/* for(x=0;x<bids;x++) { */
	/* 	free(&bids_ptr[x]); */
	/* } */

	/* if (line) */
	/* 	free(line); */
	exit(EXIT_SUCCESS);
}
