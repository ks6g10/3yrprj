
#define NGOODS (32)

struct bid {
	unsigned int offer;
	unsigned int conf[1+((NGOODS-1)/32)];
	unsigned int count;
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

unsigned short q[NGOODS];
