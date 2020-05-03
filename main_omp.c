#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*
 * Moving to windows years later, below guards fails to detect/compile in CLion
 * However identical code detected in main_ai.c. Neither actually runs
 */

//#ifdef _OPENMP
#include <omp.h>
//#endif

#include "go.h"

// how long to run MCTS in seconds
#define INTERVAL 10
#define TREE_NODE_GROWTH_INCREMENT 80;

struct tree_node* tree_node_new();
void tree_node_destroy(struct tree_node* self);
unsigned int ucb1_search(struct tree_node* node);

int  read_in_board(FILE* file, binmatrix* black, binmatrix* white);
void run_mcts(binmatrix* black_in, binmatrix* white_in, binmatrix* old_black_in, binmatrix* old_white_in, enum side ai_side);
void random_simulation(binmatrix* black, binmatrix* white, binmatrix* old_black, binmatrix* old_white, unsigned int* blk_score, unsigned int* wht_score, unsigned int max_moves, enum side start_turn);

struct tree_node
{
	struct point location;
	int hits;
	int success;
	int num_children;
	int capacity;
	struct tree_node** children; // dyanamic array of tree_node*
	struct tree_node*  parent;
};
struct tree_node *root;

//# ifdef _OPENMP
// for ease of coding make rand_r intermediate values global
unsigned int* random_store;
//# endif

int
main(int argc, char* argv[])
{
	argc--; // not sure if safe, but makes lower calculations less dumb (off by one everywhere)
	if (argc < 2 || argc > 3)
	{
		fprintf(stderr, "usage: go_ai current_board.txt [previous_board.txt] BLACK | WHITE\n");
		exit(0);
	}

	binmatrix black[BOARD_SIZE], old_black[BOARD_SIZE];
	binmatrix white[BOARD_SIZE], old_white[BOARD_SIZE];

	FILE* current = fopen(argv[1], "r");
	if (current==NULL)
	{
		fprintf(stderr, "Error: opening fle %s (%s)\n", argv[1], strerror(errno));
		exit(-1);
	}
	read_in_board(current, black, white);
	fclose(current);

	if (argc == 3)
	{
		FILE* previous = fopen(argv[2], "r");
		if (previous==NULL)
		{
			fprintf(stderr, "Error: opening fle %s (%s)\n", argv[2], strerror(errno));
			exit(-1);
		}
		read_in_board(previous, old_black, old_white);
		fclose(previous);
	}

	enum side ai_side;
	if (strcasecmp(argv[argc], "black") == 0)
		ai_side = BLACK;
	else if (strcasecmp(argv[argc], "white") == 0)
		ai_side = WHITE;
	else
	{
		fprintf(stderr, "AI side must be one of BLACK or WHITE\n");
		exit(0);
	}

	int num_threads = 32;
//# ifdef _OPENMP
	omp_set_num_threads(omp_get_num_procs());
	random_store = malloc(num_threads * sizeof(unsigned int));
//# endif

	// add to command line later
	unsigned int black_score = 0;
	unsigned int white_score = 0;

	root = tree_node_new();
	root->hits = 1; // to prevent math errors, and we have one hit already (the move that got us to present state)

	printf("Working, check this for correctness in the mean time\n");
	prettyBoard(black, white);
	printf("\n");

	time_t start = time(NULL);
	time_t end = start + INTERVAL;

# pragma omp parallel num_threads(num_threads)
	{
//# ifdef _OPENMP
		int my_thread = omp_get_thread_num();
		random_store[my_thread]=(start^my_thread);
//# endif
		srand((unsigned) start);
		do
		{
			run_mcts(black, white, old_black, old_white, ai_side);
		}
		while (time(NULL) < end);
	}
	printf("MCTS done %d games played\n", root->hits);
	fflush(stdout);
	
	// evaluate top level of tree and pick best before updating and printing
	// Will segfault until tree update has been written, ucb1_search indexs into first position which is still a NULL pointer
	unsigned int index = ucb1_search(root);
	struct point chosen_move = root->children[index]->location;
	enum side turn = ai_side;
	play_move(black, white, old_black, old_white, &black_score, &white_score, chosen_move, &turn);

	if (ai_side==BLACK)
		printf("Black plays: ");
	else
		printf("White plays: ");
	printf("(%d,%d)\n", chosen_move.row, chosen_move.col);
	prettyBoard(black, white);
	
	tree_node_destroy(root);
	return 0;
}

int // empty_spots, here for convenience
read_in_board(FILE* file, binmatrix* black, binmatrix* white)
{
	char line[BOARD_SIZE];
	int empty_spots = 0;
	for (int i=0; i<BOARD_SIZE; i++)
	{
		fscanf(file, "%s\n", line);
		white[i]=0;
		black[i]=0;
		for (int j=0; j<BOARD_SIZE; j++)
		{
			char ch = line[j];
			if (ch == 'X' || ch == 'x')
				black[i] |= 1<<(19-j-1);
			else if (ch == 'O' || ch == 'o')
				white[i] |= 1<<(19-j-1);
			else
				empty_spots++;
		}
	}

	return empty_spots;
}

struct tree_node*
tree_node_new()
{
	struct tree_node *self = malloc(sizeof(struct tree_node));
	self->location.row = 0;
	self->location.col = 0;
	self->hits = 0;
	self->success = 0;
	self->num_children = 0;
	self->capacity = TREE_NODE_GROWTH_INCREMENT;
	self->children = malloc(sizeof(struct tree_node*)*self->capacity);
	self->parent = NULL;
	return self;
}

void
tree_node_destroy(struct tree_node* self)
{
	for (int i=0; i<self->num_children; i++)
		tree_node_destroy(self->children[i]);
	free(self->children);
	free(self);
}

unsigned int
ucb1_search(struct tree_node* node)
{
	unsigned int num_games = root->hits;
	unsigned int best_index = 0;
	double best_value = 0.0;
	double cur_value = 0.0;
	struct tree_node* cur;
	for (unsigned int i=0; i<node->num_children; i++)
	{
		cur = node->children[i];
		cur_value = sqrt((2.0*log(num_games)) / cur->hits);
		if (cur_value > best_value)
		{
			best_value = cur_value;
			best_index=i;
		}
	}
	
	return best_index;
}

int
tree_node_find_child(struct tree_node* node, struct point move)
{
	for (int i=0; i<node->num_children; i++)
	{
		//printf("stuck in find_child\n");
		struct point child_move = node->children[i]->location;
		if (child_move.row==move.row && child_move.col==move.col)
			return TRUE;
	}
	return FALSE;
}

void
tree_node_add_child(struct tree_node* parent, struct tree_node* new_child)
{
	if (parent->num_children+1 > parent->capacity)
	// kind of obvious this is a critical section, somehow missed it for awhile, never had any obvious errors
# pragma omp critical(growth)
	{
			parent->capacity += TREE_NODE_GROWTH_INCREMENT;
			struct tree_node** new_array = malloc(parent->capacity*sizeof(struct tree_node*));
			for (int i=0; i<parent->num_children; i++)
				new_array[i]=parent->children[i];
			free(parent->children);
			parent->children = new_array;
	}

	int index = parent->num_children;
	parent->children[index] = new_child;
	parent->num_children++;
}

// initial board state is feed in, updating by traversing root doesn't make sense if given board is 1) not empty and 2) data is already given to us
// copy board, progress through tree updating copy, when UCT says stop add new node or go straight to playouts, after playouts update tree with findings
void
run_mcts(binmatrix* black_in, binmatrix* white_in, binmatrix* old_black_in, binmatrix* old_white_in, enum side ai_side)
{
	binmatrix* black = binmatrix_clone(black_in);
	binmatrix* white = binmatrix_clone(white_in);

	binmatrix* old_black = binmatrix_clone(old_black_in);
	binmatrix* old_white = binmatrix_clone(old_white_in);

	unsigned int black_score = 0;
	unsigned int white_score = 0;

	enum side turn = ai_side; // main sets ai_side, tree descending gets this to proper value
	
	// descend through tree
	struct tree_node* cur = root;
	struct tree_node* new_leaf = NULL;
	enum side new_node_side;
	int free_spots = (19*19) - binmatrix_count(black)+binmatrix_count(white);
	// num_children==free_spots is a nice heuristic, but may be impossible to reach (suicide and ko) so 3/4 way there is a good compromise (plus get to rest of mcts faster)
	// free_spots >= 20 is for case where many remaining spots are suicide spots

# pragma omp critical
	{
		while (cur->num_children != 0 && cur->num_children >= (0.75*free_spots) && free_spots >= 20)
		{
			int index = ucb1_search(cur);
			cur = cur->children[index];
			play_move(black, white, old_black, old_white, &black_score, &white_score, cur->location, &turn);
			free_spots = (19*19) - binmatrix_count(black)+binmatrix_count(white);
		}

		// pick new move
		struct point move;
		new_node_side = turn; // save side of current turn for backprop
		int flag = 0;
		for (int attempts=0; attempts<100; attempts++)
		{
			//printf("stuck in pick\n");
			move = pick_random_move(white, black);
			// make sure move hasn't been played before
			if (tree_node_find_child(cur, move))
				continue;
			// make sure move is legal (also plays the move)
			if (play_move(black, white, old_black, old_white, &black_score, &white_score, move, &turn))
			{
				flag = 1;
				break;
			}
		}
		// failed to find an available move in a timely manner
		// Could not find random spot in a timely manner, find first empty spot and return that
		if (!flag)
		{
			for (move.row=0; move.row<BOARD_SIZE; move.row++)
			{
				int row_sum = black[move.row]+white[move.row];
				int temp = 1;
				for (move.col=0; move.col<BOARD_SIZE; move.col++)
				{
					if (row_sum & temp && TRUE == play_move(black, white, old_black, old_white, &black_score, &white_score, move, &turn))
					{
						//fprintf(stderr, "Found a sequential move!\n");
						flag = 1;
						break;
					}
					else
						temp = temp << 1;
				}

			if (flag)
				break;
			}

			// if we got here naturally and not via flag then panic
			if (!flag)
			{
				fprintf(stderr, "No moves available, dying now\n");
				exit(-999);
			}
		}

		// make new leaf node for new move
		new_leaf = tree_node_new();
		new_leaf->location = move;
		tree_node_add_child(cur, new_leaf);
		new_leaf->parent = cur;
	}

	// run simulation
	free_spots = (19*19) - (binmatrix_count(black)+(binmatrix_count(white)));
	free_spots /= 2;
	// if free_spots is very low, remaining spots might all mostly suicides so reduce number of moves in playout hopefully avoids infintie looping
	if (free_spots <= 20)
		free_spots /= 2;
	random_simulation(black, white, old_black, old_white, &black_score, &white_score, free_spots, turn);

	//printf("Results of random play\n");
	//prettyBoard(black, white);

	// run benson's
	//printf("Benson's pass alive stones\n");
	binmatrix* alive_black = bensons_pass_alive(black, white);
	binmatrix* alive_white = bensons_pass_alive(white, black);
	//prettyBoard(alive_black, alive_white);

	// estimate score
	estimate_score(alive_black, alive_white, &black_score, &white_score);
	//printf("finished scoring: B %d, W %d\n", black_score, white_score);
	
	// determine winner (tie counts as win, another code simplifier)
	enum side winner;
	if (black_score > white_score)
		winner = BLACK;
	else if (white_score > black_score)
		winner = WHITE;
	else
		winner = new_node_side;
	
	// back prop (winner turns get +1 win, all else just get +1 hits)
	
# pragma omp critical
	{
		turn = new_node_side;
		cur = new_leaf;
		while (cur != root)
		{
			//printf("Climbing tree\n");
			if (turn == winner)
				cur->success++;
			cur->hits++;

			cur = cur->parent;
			turn *= -1;
		}
		root->hits++;
	}

	free(black);
	free(white);
	free(old_black);
	free(old_white);
}

void
random_simulation(
	binmatrix* black, // in/out
	binmatrix* white, // in/out
	binmatrix* old_black, // in/out
	binmatrix* old_white, // in/out
	unsigned int* blk_score,
	unsigned int* wht_score,
	unsigned int max_moves,
	enum side start_turn)
{
	enum side cur_turn = start_turn;
	struct point move;

	for (unsigned int num_moves=0; num_moves < max_moves; num_moves++)
	{
		//printf("Stuck inside simulations\n");
		int flag = 0;
		for (int attempts=0; attempts<100; attempts++)
		{
			move = pick_random_move(black, white);
			if (TRUE == play_move(black, white, old_black, old_white, blk_score, wht_score, move, &cur_turn))
			{
				flag = 1;
				break;
			}
		}

		// failed to find an available move in a timely manner
		// Could not find random spot in a timely manner, find first empty spot and return that
		if (!flag)
		{
			fprintf(stderr, "Simulation failed, brute force first available move\n");
			for (move.row=0; move.row<BOARD_SIZE; move.row++)
			{
				int row_sum = black[move.row]+white[move.row];
				int temp = 1;
				for (move.col=0; move.col<BOARD_SIZE; move.col++)
				{
					if (row_sum & temp && TRUE == play_move(black, white, old_black, old_white, blk_score, wht_score, move, &cur_turn))
					{
						flag = 1;
						break;
					}
					else
						temp = temp << 1;
				}

				if (flag)
					break;
			}

			// if we got here naturally and not via flag then panic
			if (!flag)
			{
				fprintf(stderr, "No moves available, dying now\n");
				exit(-999);
			}
		}
	}
}
