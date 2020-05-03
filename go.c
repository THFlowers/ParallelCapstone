#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "go.h"

#ifdef VECTOR
#include <immintrin.h>
#include "mmintrin.h"
#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"
#endif

#ifdef _OPENMP
extern unsigned int* random_store;
#endif

//// Binary Helper functions

// Variation of Brian Kernighan's algorithm
// Found at geeksforgeeks.org/count-set-bits-in-an-integer
unsigned int
count_set_bits(unsigned int n)
{
	unsigned int count = 0;
	while (n)
	{
		n &= (n-1);
		count++;
	}
	return count;
}

void
printBinary(const int input)
{
	unsigned char byte, mask;

	for (int i=sizeof(int)-1; i>=0; i--)
	{
		byte = (input >> i*8) & 0xFF;
		mask = 1 << 7;
		for (int j=0; j<8; j++)
		{
			if (byte & mask)
				putchar('1');
			else
				putchar('0');
			mask = mask >> 1;

			if (8*i+j == 20)
				putchar(' ');
		}
	}
	putchar('\n');
}

void
printBinmatrix(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
		printBinary(input[i]);
}

void
printBoard(binmatrix* black, binmatrix* white)
{
	for (int i=0; i<BOARD_SIZE; i++)
		printBinary(black[i]+white[i]);
}

void
prettyBoard(binmatrix* black, binmatrix* white)
{
	char blk_sym='X', wht_sym='O', emp_sym='*';
	unsigned char wht_byte, blk_byte, mask;
	int num_skip, max_skip=sizeof(int)*8-BOARD_SIZE;

	for (int line=0; line<BOARD_SIZE; line++)
	{
		num_skip = 0;
		for (int i=sizeof(int)-1; i>=0; i--)
		{
			wht_byte = (white[line] >> i*8) & 0xFF;
			blk_byte = (black[line] >> i*8) & 0xFF;
			mask = 1 << 7;
			for (int j=0; j<8; j++)
			{
				if (num_skip < max_skip)
					num_skip++;
				else
				{
					if (blk_byte & mask)
						putchar(blk_sym);
					else if (wht_byte & mask)
						putchar(wht_sym);
					else
						putchar(emp_sym);
				}
				mask = mask >> 1;
			}
		}
		putchar('\n');
	}
	putchar('\n');
}

//// Go play methods
struct point
pick_random_move(binmatrix *black, binmatrix *white)
{
	struct point pos;
	unsigned int row_sum;
	for (int i=0; i<100; i++)
	{
//#ifndef _OPENMP
		pos.row = rand() % BOARD_SIZE;
		pos.col = rand() % BOARD_SIZE;
/*
#else
		int my_thread = omp_get_thread_num();
		pos.row = rand_r(random_store+my_thread) % BOARD_SIZE;
		pos.col = rand_r(random_store+my_thread) % BOARD_SIZE;
#endif
 */
		row_sum = black[pos.row]+white[pos.row];
		if (!(row_sum & (1<<(BOARD_SIZE-pos.col-1))))
			return pos;
	}

	// Could not find random spot in a timely manner, find first empty spot and return that
	fprintf(stderr, "Random failed, brute force empty location\n");
	for (pos.row=0; pos.row<BOARD_SIZE; pos.row++)
	{
		row_sum = black[pos.row] + white[pos.row];
		int temp = 1;
		for (pos.col=0; pos.col<BOARD_SIZE; pos.col++)
		{
			if (row_sum & temp)
				return pos;
			temp = temp << 1;
		}
	}

	fprintf(stderr, "What the hell!?\n");
	exit(-999);
}

// algorithm by Kjeld Petersen, see binmatrix methods for more information
int
play_move(
	binmatrix *black,
	binmatrix *white,
	binmatrix *old_black,
	binmatrix *old_white,
	unsigned int *black_score,
	unsigned int *white_score,
	struct point pos,
	enum side *turn)
{
	int capture = 0;
	binmatrix *P, *S, *N = NULL;
	if (*turn == BLACK)
	{
		P = binmatrix_clone(black);
		binmatrix_set(P, pos);

		S = binmatrix_new();
		binmatrix_set(S, pos);

		// remember most binmatrix methods modify left arg
		binmatrix_neighbors4(S);
			// printf("S.neighbors4()\n"); printBinmatrix(S); printf("\n");
		binmatrix_and(S, white);
			// printf("S.and(W)\n"); printBinmatrix(S); printf("\n");
		if (binmatrix_is_empty(S))
			capture = 0;
		else
		{
			N = binmatrix_clone(binmatrix_expand_to(S, white));
				// printf("N.copy(S.expandTo(W)\n"); printBinmatrix(N); printf("\n");
			binmatrix_minus(binmatrix_neighbors4(N), P);
				// printf("N.neighbors4().minus(P)\n"); printBinmatrix(N); printf("\n");
			binmatrix_and(binmatrix_neighbors4(N), S);
				// printf("N.neighbors4().and(S)\n"); printBinmatrix(N); printf("\n");
			binmatrix_expand_to(N, S);
				// printf("N.expandTo(S)\n"); printBinmatrix(N); printf("\n");
			capture = binmatrix_count(binmatrix_xor(N, S));
				// printf("N.xor(S) -- %d captured\n", capture); printBinmatrix(N); printf("\n");
		}

		if (capture == 0)
		{
			binmatrix_set(binmatrix_none(S), pos);
			binmatrix_neighbors4(binmatrix_expand_to(S, P));
			binmatrix_minus(S, white);
			if (binmatrix_is_empty(S))
				return FALSE; // suicide
		}
		else
		{
			if (capture == 1 && binmatrix_equals(P, old_black))
				return FALSE; // ko
			*black_score += capture;
			binmatrix_copy(old_white, white);
			binmatrix_xor(white, N);
		}
		binmatrix_copy(old_black, black);
		binmatrix_set(black, pos);
	}
	else
	{
		P = binmatrix_clone(white);
		binmatrix_set(P, pos);

		S = binmatrix_new();
		binmatrix_set(S, pos);

		binmatrix_neighbors4(S);
		binmatrix_and(S, black);
		if (binmatrix_is_empty(S))
			capture = 0;
		else
		{
			N = binmatrix_clone(binmatrix_expand_to(S, black));
			binmatrix_minus(binmatrix_neighbors4(N), P);
			binmatrix_and(binmatrix_neighbors4(N), S);
			binmatrix_expand_to(N, S);
			capture = binmatrix_count(binmatrix_xor(N, S));
		}

		if (capture == 0)
		{
			binmatrix_set(binmatrix_none(S), pos);
			binmatrix_neighbors4(binmatrix_expand_to(S, P));
			binmatrix_minus(S, black);
			if (binmatrix_is_empty(S))
				return FALSE; // suicide
		}
		else
		{
			if (capture == 1 && binmatrix_equals(P, old_white))
				return FALSE; // ko
			*white_score += capture;
			binmatrix_copy(old_black, black);
			binmatrix_xor(black, N);
		}
		binmatrix_copy(old_white, white);
		binmatrix_set(white, pos);
	}

	free(P);
	free(S);
	if (N != NULL)
		free(N);

	*turn = -(*turn);
	return TRUE;
}

// Adaptated to the binmatrix technique by Kjeld Petersen,
// It is known to have a few incorrect edge cases (good enough for our uses though)
// returns alive black stones, swap order in call to get alive white stones)
binmatrix*
bensons_pass_alive(binmatrix* black, binmatrix* white)
{
	binmatrix* L = binmatrix_minus(binmatrix_neighbors4(binmatrix_clone(black)), white);
	binmatrix* t = binmatrix_or(binmatrix_clone(white), L); // original code forgot about dealloc
	binmatrix* P = binmatrix_expand_to(binmatrix_clone(L), t);
	binmatrix* E = binmatrix_minus(binmatrix_neighbors4(binmatrix_clone(P)), black);
	binmatrix* R = binmatrix_xor(binmatrix_expand_to(binmatrix_grow4(binmatrix_clone(E)), P), P);

	binmatrix* N = binmatrix_new();
	binmatrix* M = binmatrix_new();
	binmatrix* C = binmatrix_new();
	binmatrix* F = binmatrix_new();

	do
	{
		binmatrix_copy(C, R);
		while (FALSE == binmatrix_is_empty(C))
		{
			binmatrix_expand_to(binmatrix_first(binmatrix_copy(F, C)), R);
			binmatrix_xor(binmatrix_expand_to(binmatrix_neighbors4(binmatrix_xor(binmatrix_copy(t, R), F)), black), black);
			binmatrix_or(N, t);
			binmatrix_xor(C, F);
		}
		binmatrix_expand_to(binmatrix_neighbors4(binmatrix_copy(M, N)), R);
		binmatrix_minus(R, M);
	} while (FALSE == binmatrix_is_empty(M));

	binmatrix_xor(binmatrix_copy(t, black), N);

	free(C);
	free(E);
	free(F);
	free(L);
	free(M);
	free(N);
	free(P);
	free(R);
	
	return t;
}

// Basically Bouzy's 5/21 algorithm
//
// Algorithm originaly by Bruno Bouzy in his dissertation "Les Ensembles Flous Au jeu De Go"
// Found in Gnu Go documentation chapter 17 section 2 "Bouzy's 5/21 Algorithm" (hosted at www.delorie.com)
//
// Converts from binmatrix to short[][] (vals are 128 * BLACK/WHITE)
void
estimate_score(
	binmatrix* alive_black, // in
	binmatrix* alive_white, // in
	unsigned int* black_score, // in/out
	unsigned int* white_score) // in/out
{
	unsigned short i, j, count, adj_pos, adj_neg, adj_zero;
	short temp;
	short matrix[BOARD_SIZE][BOARD_SIZE]; // read in initial state, use for reading current state
	short new_matrix[BOARD_SIZE][BOARD_SIZE]; // write to new_matrix then copy values over to matrix, then wipe
	for (int i=0; i<BOARD_SIZE; i++)
	{
		for (j=0; j<BOARD_SIZE; j++)
		{
			if (alive_black[i] & (1<<(BOARD_SIZE-j-1)))
				matrix[i][j] = BLACK * 128;
			else if (alive_white[i] & (1<<(BOARD_SIZE-j-1)))
				matrix[i][j] = WHITE * 128;
			else
				matrix[i][j] = 0;
		}
	}

	/*
	printf("Initial State\n");
	for (int i=0; i<BOARD_SIZE; i++)
	{
		for (j=0; j<BOARD_SIZE; j++)
			printf("%4d ", matrix[i][j]);
		printf("\n");
	}
	printf("\n");
	*/

	// dilation
	for (count=0; count<5; count++)
	{
		memset(new_matrix, 0, BOARD_SIZE*BOARD_SIZE*sizeof(short));

		for (i=0; i<BOARD_SIZE; i++)
		{
			for (j=0; j<BOARD_SIZE; j++)
			{
				adj_pos = 0;
				adj_neg = 0;
				if (i-1 >=0)
				{
					temp = matrix[i-1][j];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
				}
				if (i+1 < BOARD_SIZE)
				{
					temp = matrix[i+1][j];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
				}
				if (j-1 >=0)
				{
					temp = matrix[i][j-1];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
				}
				if (j+1 < BOARD_SIZE)
				{
					temp = matrix[i][j+1];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
				}

				// don't mark empty space as owned by anyone
				if (matrix[i][j]==0 && adj_neg==0 && adj_pos==0)
					continue;
				
				if (matrix[i][j] >= 0 && adj_neg==0)
					new_matrix[i][j] = matrix[i][j] + adj_pos;
				else if (matrix[i][j] <= 0 && adj_pos==0)
					new_matrix[i][j] = matrix[i][j] - adj_neg;
			}
		}

		memcpy(matrix, new_matrix, BOARD_SIZE*BOARD_SIZE*sizeof(short));

	}

	/*
	printf("After %d dilations\n", count);
	for (int i2=0; i2<BOARD_SIZE; i2++)
	{
		for (int j2=0; j2<BOARD_SIZE; j2++)
			printf("%4d ", matrix[i2][j2]);
		printf("\n");
	}
	printf("\n");
	*/

	// erosion
	for (count=0; count<21; count++)
	{
		memset(new_matrix, 0, BOARD_SIZE*BOARD_SIZE*sizeof(short));

		for (i=0; i<BOARD_SIZE; i++)
		{
			for (j=0; j<BOARD_SIZE; j++)
			{
				adj_pos  = 0;
				adj_neg  = 0;
				adj_zero = 0;
				if (i-1 >=0)
				{
					temp = matrix[i-1][j];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
					else
						adj_zero++;
				}
				if (i+1 < BOARD_SIZE)
				{
					temp = matrix[i+1][j];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
					else
						adj_zero++;
				}
				if (j-1 >=0)
				{
					temp = matrix[i][j-1];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
					else
						adj_zero++;
				}
				if (j+1 < BOARD_SIZE)
				{
					temp = matrix[i][j+1];
					if (temp > 0)
						adj_pos++;
					else if (temp < 0)
						adj_neg++;
					else
						adj_zero++;
				}

				if (matrix[i][j] > 0)
				{
					new_matrix[i][j] = matrix[i][j] - (adj_neg+adj_zero);
					if (new_matrix[i][j] < 0)
						new_matrix[i][j] = 0;
				}
				else if (matrix[i][j] < 0)
				{
					new_matrix[i][j] = matrix[i][j] + (adj_pos+adj_zero);
					if (new_matrix[i][j] > 0)
						new_matrix[i][j] = 0;
				}
			}
		}

		memcpy(matrix, new_matrix, BOARD_SIZE*BOARD_SIZE*sizeof(short));
	}

	//printf("After %d Erosions\n", count);
	for (int i=0; i<BOARD_SIZE; i++)
	{
		for (j=0; j<BOARD_SIZE; j++)
		{
			// recall BLACK=1, WHITE=-1
			if (matrix[i][j] >= BLACK)
				*black_score++;
			else if (matrix[i][j] <= WHITE)
				*white_score++;
			//printf("%4d ", matrix[i][j]);
		}
		//printf("\n");
	}
	//printf("\n");
}


/* ----- binmatrix methods -----
 * binmatrix algorithms and general methodology for go rule enforcement by Kjeld Petersen
 * For more information refer to (https://senseis.xmp.net/?BinMatrix)
 * Available for use without permission, provided atrribution is given
 */

binmatrix*
binmatrix_new()
{
	binmatrix* temp = malloc(sizeof(binmatrix)*BOARD_SIZE);
	binmatrix_none(temp);
	return temp;
}

// Manipulators
binmatrix*
binmatrix_set(binmatrix* input, struct point pos)
{
	input[pos.row]|=1<<(BOARD_SIZE-pos.col-1);
	return input;
}

binmatrix*
binmatrix_clear(binmatrix* input, struct point pos)
{
	input[pos.row]|=1<<(BOARD_SIZE-pos.col-1);
	return input;
}

binmatrix*
binmatrix_toggle(binmatrix* input, struct point pos)
{
	input[pos.row]^=1<<(BOARD_SIZE-pos.col-1);
	return input;
}

binmatrix*
binmatrix_all(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
		input[i] = SIZE_MASK;
	return input;
}

binmatrix*
binmatrix_none(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
		input[i] = 0;
	return input;
}

binmatrix*
binmatrix_copy(binmatrix* self, binmatrix* other)
{
	for (int i=0; i<BOARD_SIZE; i++)
		self[i] = other[i];
	return self;
}

binmatrix* // copy, be sure to free() manually!
binmatrix_clone(binmatrix* input)
{
	binmatrix* copy = binmatrix_new();
	memcpy(copy, input, BOARD_SIZE*sizeof(binmatrix));
	return copy;
}

// Binary operators
binmatrix*
binmatrix_and(binmatrix* self, binmatrix* other)
{
	int i=0;
#ifdef VECTOR
	__m128i s,o, res;
	for (; i<BOARD_SIZE-4; i+=4)
	{
		s = _mm_loadu_si128(self+i);
		o = _mm_loadu_si128(other+i);
		res = _mm_and_si128(s,o);
		_mm_stream_si128(self+i, res);
	}
#endif
	for (; i<BOARD_SIZE; i++)
		self[i] &= other[i];

	return self;
}

binmatrix*
binmatrix_or(binmatrix* self, binmatrix* other)
{
	int i=0;
#ifdef VECTOR
	__m128i s,o, res;
	for (; i<BOARD_SIZE-4; i+=4)
	{
		s = _mm_loadu_si128(self+i);
		o = _mm_loadu_si128(other+i);
		res = _mm_or_si128(s,o);
		_mm_stream_si128(self+i, res);
	}
#endif
	for (; i<BOARD_SIZE; i++)
		self[i] |= other[i];

	return self;
}

binmatrix*
binmatrix_xor(binmatrix* self, binmatrix* other)
{
	int i=0;
#ifdef VECTOR
	__m128i s,o, res;
	for (; i<BOARD_SIZE-4; i+=4)
	{
		s = _mm_loadu_si128(self+i);
		o = _mm_loadu_si128(other+i);
		res = _mm_xor_si128(s,o);
		_mm_stream_si128(self+i, res);
	}
#endif
	for (; i<BOARD_SIZE; i++)
		self[i] ^= other[i];

	return self;
}

binmatrix*
binmatrix_minus(binmatrix* self, binmatrix* other)
{
	int i=0;
#ifdef VECTOR
	__m128i s,o, res;
	for (; i<BOARD_SIZE-4; i+=4)
	{
		s = _mm_loadu_si128(self+i);
		o = _mm_loadu_si128(other+i);
		res = _mm_andnot_si128(s,o);
		_mm_stream_si128(self+i, res);
	}
#endif
	for (; i<BOARD_SIZE; i++)
		self[i] &= ~other[i];

	return self;
}

binmatrix*
binmatrix_expand_to(binmatrix* self, binmatrix* other)
{
	unsigned int old_count = 0;
	do
	{
		old_count = binmatrix_count(self);
		binmatrix_grow4(self);
		binmatrix_and(self, other);
	}
	while (old_count != binmatrix_count(self));

	return self;
}

// Unary operators
binmatrix*
binmatrix_inverse(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
		input[i] = (~input[i]) & SIZE_MASK;
	return input;
}

binmatrix*
binmatrix_increment(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
	{
		input[i]++;
		if (input[i] < SIZE_MASK+1)
			break;
	}
	return input;
}

binmatrix*
binmatrix_twos_complement(binmatrix* input)
{
	binmatrix_inverse(input);
	binmatrix_increment(input);
	return input;
}

binmatrix*
binmatrix_first(binmatrix* input)
{
	binmatrix* complement = binmatrix_twos_complement(binmatrix_clone(input));
	binmatrix_and(input, complement);
	free(complement);
	return input;
}

binmatrix*
binmatrix_first_chain(binmatrix* input)
{
	// should be equivalent to the specification but supports chain calls
	binmatrix* copy = binmatrix_clone(input);
	binmatrix_expand_to(binmatrix_first(input), copy);
	free(copy);
	return input;
}

binmatrix*
binmatrix_up(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE-1; i++)
		input[i] = input[i+1];
	input[BOARD_SIZE] = 0;
	return input;
}

binmatrix*
binmatrix_down(binmatrix* input)
{
	for (int i=BOARD_SIZE-1; i>0; i--)
		input[i] = input[i-1];
	input[0] = 0;
	return input;
}

binmatrix*
binmatrix_left(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
		input[i] = (input[i] << 1) & SIZE_MASK;
	return input;
}

binmatrix*
binmatrix_right(binmatrix* input)
{
	for (int i=0; i<BOARD_SIZE; i++)
		input[i] = (input[i] >> 1) & SIZE_MASK;
	return input;
}

binmatrix*
binmatrix_explode4(binmatrix* input)
{
	binmatrix* up = binmatrix_up(binmatrix_clone(input));
	binmatrix* dn = binmatrix_down(binmatrix_clone(input));
	binmatrix* lf = binmatrix_left(binmatrix_clone(input));
	binmatrix* ri = binmatrix_right(binmatrix_clone(input));

	binmatrix_none(input);
	binmatrix_or(input, up);
	binmatrix_or(input, dn);
	binmatrix_or(input, lf);
	binmatrix_or(input, ri);

	free(up);
	free(dn);
	free(lf);
	free(ri);

	return input;
}

binmatrix*
binmatrix_grow4(binmatrix* input)
{
	binmatrix* copy = binmatrix_clone(input);
	binmatrix_or(binmatrix_explode4(input), copy);
	free(copy);
	return input;
}

binmatrix*
binmatrix_neighbors4(binmatrix* input)
{
	binmatrix* copy = binmatrix_clone(input);
	binmatrix_and(binmatrix_explode4(input), binmatrix_inverse(copy));
	free(copy);
	return input;
}

// Testing
int
binmatrix_is(binmatrix* input, struct point pos)
{
	return input[pos.row] & (1<<(BOARD_SIZE-pos.col-1));
}

int
binmatrix_equals(binmatrix* self, binmatrix* other)
{
	int res = 0;
	for (int i=0; i<BOARD_SIZE; i++)
	{
		res = self[i] ^ other[i];
		if (res != 0)
			return FALSE;
	}
	return TRUE;
}

int
binmatrix_is_empty(binmatrix* self)
{
	for (int i=0; i<BOARD_SIZE; i++)
		if (self[i] != 0)
			return FALSE;
	return TRUE;
}

int
binmatrix_is_full(binmatrix* self)
{
	for (int i=0; i<BOARD_SIZE; i++)
		if (self[i] != SIZE_MASK)
			return FALSE;
	return TRUE;
}

unsigned int
binmatrix_count(binmatrix* self)
{
	unsigned int count = 0;
	for (int i=0; i<BOARD_SIZE; i++)
		count += count_set_bits(self[i]);
	return count;
}
