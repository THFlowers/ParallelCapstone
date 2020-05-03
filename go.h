
#ifndef GO_HEADER
#define GO_HEADER

#define BOARD_SIZE 19
#define MAX_SQUARES (BOARD_SIZE*BOARD_SIZE)
#define TRUE 1
#define FALSE 0
// 524287 is (1<<BOARD_SIZE)-1 for a 19x19
#define SIZE_MASK 524287
//#define SIZE_MASK ((1<<BOARD_SIZE)-1)

typedef unsigned int binmatrix;
struct point { int row; int col; };
enum side {BLACK=1, WHITE=-1}; /* multiply by -1 to change turns */

unsigned int count_set_bits(unsigned int n);
void printBinary(int input);
void printBinmatrix(binmatrix* input);
void printBoard(binmatrix* black, binmatrix* white);

void prettyBoard(binmatrix* black, binmatrix* white);

struct point pick_random_move(binmatrix *black, binmatrix *white);
int play_move(
	binmatrix *black,
	binmatrix *white,
	binmatrix *old_black,
	binmatrix *old_white,
	unsigned int *black_score,
	unsigned int *white_score,
	struct point pos,
	enum side *turn);

binmatrix* bensons_pass_alive(binmatrix* black, binmatrix* white);

void estimate_score(
	binmatrix* black, // in
	binmatrix* white, // in
	unsigned int* blk_score,  // in/out,
	unsigned int* wht_score); // in/out

// unless noted, each function modifies self (left arg)
binmatrix* binmatrix_new();
binmatrix* binmatrix_set(binmatrix* input, struct point pos);
binmatrix* binmatrix_clear(binmatrix* input, struct point pos);
binmatrix* binmatrix_toggle(binmatrix* input, struct point pos);
binmatrix* binmatrix_all(binmatrix* input);
binmatrix* binmatrix_none(binmatrix* input);
binmatrix* binmatrix_copy(binmatrix* self, binmatrix* other);
binmatrix* binmatrix_clone(binmatrix* input); // copies!

binmatrix* binmatrix_and(binmatrix* self, binmatrix* other);
binmatrix* binmatrix_or(binmatrix* self, binmatrix* other);
binmatrix* binmatrix_xor(binmatrix* self, binmatrix* other);
binmatrix* binmatrix_minus(binmatrix* self, binmatrix* other);
binmatrix* binmatrix_expand_to(binmatrix* self, binmatrix* other);

binmatrix* binmatrix_inverse(binmatrix* input);
binmatrix* binmatrix_increment(binmatrix* input);
binmatrix* binmatrix_twos_complement(binmatrix* input);
binmatrix* binmatrix_first(binmatrix* input);
binmatrix* binmatrix_first_chain(binmatrix* input);
binmatrix* binmatrix_up(binmatrix* input);
binmatrix* binmatrix_down(binmatrix* input);
binmatrix* binmatrix_left(binmatrix* input);
binmatrix* binmatrix_right(binmatrix* input);
binmatrix* binmatrix_explode4(binmatrix* input);
binmatrix* binmatrix_grow4(binmatrix* input);
binmatrix* binmatrix_neighbors4(binmatrix* input);

int binmatrix_is(binmatrix* input, struct point pos);
int binmatrix_equals(binmatrix* self, binmatrix* other);
int binmatrix_is_empty(binmatrix* self);
int binmatrix_is_full(binmatrix* self);
unsigned int binmatrix_count(binmatrix* self);
#endif
