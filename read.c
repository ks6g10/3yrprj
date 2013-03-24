#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


struct bid {
	double offer;
	unsigned int id;
	unsigned int conf[];
};

int main(void)   {
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
	unsigned int id;
	unsigned int ints =0;
	unsigned int *tmp;
	struct bid * bids_ptr;
	unsigned int bids_count =0;
	double value;
	fp = fopen("/home/ks6g10/Desktop/cats/f1.txt", "r");
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
				ints = 1+(goods-1)/32;
				tmp = malloc(sizeof(tmp)*ints);
			} else
			if(strncmp(line,s_bids,strlen(s_bids)) == 0) {
				bids = atoi(line+strlen(s_bids)+1);
				bids_ptr = calloc((sizeof(bids_ptr)+(1+(goods-1)/32)*sizeof(int)),bids);
				printf("Number of bids %u\n",bids);
			} else
			if(strncmp(line,s_dummy,strlen(s_dummy)) == 0) {
				dummy = atoi(line+strlen(s_dummy)+1);
				printf("Number of dummy %u\n",dummy);
			}
			all = !!(goods && bids && dummy);
		} else {
			if(isdigit(line[0])) {
				int x;
				char * head;
				char * tail;
				int count = 0;
				unsigned int tmp2;
				for(x = 0; x < ints; x++) {
					tmp[x] = 0;
				}
				sscanf(line,"%u\t%lf",&id,&value);	
				head = tail = line;
				while(count < 2 && *head != '\0') {
					if(*head == '\t') {
						count++;
					}
					head++;
				}
				count =0;
				tail = head--;
	
				while(*head != '#' && *head != '\0') {					
					if(*head == '\t') {
						sscanf(tail,"\t%u\t",&tmp2);
						if(tmp2 < goods) {
							tmp[tmp2/32] |= (1 << tmp2);
						}					
						tail = head;	
					}
					head++;
				}
				bids_ptr[bids_count].id = id;
				bids_ptr[bids_count].offer = value;
				for(x=0;x<ints;x++) {
					bids_ptr[bids_count].conf[x] = tmp[x];
				}
				printf("ID %u Value %f %u :\n", bids_ptr[bids_count].id,bids_ptr[bids_count].offer,bids_ptr[bids_count].conf[0]);
				bids_count++;
			}

		}
	}
	free(tmp);
	
	if (line)
		free(line);
	exit(EXIT_SUCCESS);
}
