#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


struct bid {
	double offer;
	unsigned int id;
	unsigned int bin;
	unsigned int conf[];
};

struct bid2 {
	double offer;
	unsigned int id;
	unsigned int conf[20];
};

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
	
	unsigned int ints =0;
	unsigned int *tmp;
	struct bid * bids_ptr;
	unsigned int bids_count =0;
	unsigned int * bin_count;
	int x;
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Could not open file\n");
		exit(EXIT_FAILURE);
	}
	printf("hello\n");
	while ((read = getline(&line, &len, fp)) != -1) {
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
				bids_ptr = calloc((sizeof(bids_ptr)+(1+(goods-1+dummy)/32)*sizeof(int)),bids);
				printf("Number of bids %u\n",bids);
			} else
			if(strncmp(line,s_dummy,strlen(s_dummy)) == 0) {
				dummy = atoi(line+strlen(s_dummy)+1);
				printf("Number of dummy %u\n",dummy);
				ints = 1+(goods+dummy-1)/32;
				tmp = malloc(sizeof(tmp)*ints);
				bin_count = malloc(sizeof(bin_count)*(goods+dummy));
				for(x =0;x < goods+dummy;x++) {
					bin_count[x] = 0;
				}
			}
			all = !!(goods && bids && dummy);
		} else {
			if(isdigit(line[0])) {
				int x;
				char * head;
				char * tail;
				double value;
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
//				sscanf(line,"%u\t%lf",id,value);	
				head++;
				/* while(count < 2 && *head != '\0') { */
				/* 	if(*head == '\t') { */
				/* 		count++; */
				/* 	} */
				/* 	head++; */
				/* } */
				/* count =0; */
				/* tail = head--; */
				bool = 0;
				while(*head != '#' && *head != '\0') {
					if(*head == '\t') {
						tmp2 = strtol(tail,&head,10);
						//sscanf(tail,"\t%u\t",tmp2);
						tmp[tmp2/32] |= (1 << tmp2);
						tail = head;
						if(!bool) {
							bin_count[tmp2]++;
							bool = 1;
							bids_ptr[bids_count].bin = tmp2;
							
						}
					}
					head++;
				}
				bids_ptr[bids_count].id = id;
				bids_ptr[bids_count].offer = value;
				for(x=0;x<ints;x++) {
					bids_ptr[bids_count].conf[x] = tmp[x];
				}
				printf("ID %u Value %f %lu :\n", bids_ptr[bids_count].id,bids_ptr[bids_count].offer,bids_ptr[bids_count].conf[0]+((unsigned long)bids_ptr[bids_count].conf[1] << 32));
				bids_count++;
			}

		}
	}
	
	struct bid2 **bin = malloc(sizeof(struct bid2 *)*(goods+dummy));
	for(x = 0; x < goods+dummy;x++) {
		(bin)[x] = malloc((sizeof(struct bid2)*bin_count[x]));
		bin_count[x]=0;
	}
	int y;
	for(x = 0; x < bids;x++) {
		unsigned int index = bids_ptr[x].bin;
		struct bid2 * tmp = bin[index];
		tmp[bin_count[index]].offer = 1;
//		bin[index][bin_count[index]].offer = 1;
			//	(bin[index]+bin_count[index])->offer = bids_ptr[x].offer;
//		bin[index][bin_count[index]].id = bids_ptr[x].id;
		for(y = 0; y < ints;y++) {
//			(bin[index])[bin_count[index]].conf[y] = bids_ptr[x].conf[y];
		}
		//free(&bin[index][bin_count[index]]);
		bin_count[index]++;
		
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
