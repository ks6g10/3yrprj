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

#define SINGLETON (int)
struct singleton {
	unsigned int good;
	unsigned int bid[];
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

struct sizebin * sizebin_ptr[NGOODS];
