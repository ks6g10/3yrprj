
#define NGOODS (32)
#define NBIDS (32)

struct bid_skel {
	unsigned int offer;
	unsigned int size;
	unsigned int id;
	double score;
};

struct bid {
	unsigned int conf[1+((NGOODS-1)/32)];	
	unsigned short quant[];
};

struct allocation {
	unsigned int price;
	unsigned short quant[NGOODS];
};

struct singleton {
	unsigned int bid;
	unsigned int good;
	unsigned short quant;
	float rev_p_quant;
};

struct sizebin {
	unsigned int size;
	bid bids[];
}

struct good {
	unsigned int quant;
	double score;
};

unsigned short q[NGOODS];

unsigned int * sizebin_ptr[NGOODS]
