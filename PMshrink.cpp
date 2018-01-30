#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "PMimplementation.h"


PerfectMatching::Node* PerfectMatching::FindBlossomRoot(Edge* a0, bool clear_is_outer_flag, int& size)
{
	Node* i;
	Node* j;
	Node* _i[2];
	Node* r;
	int branch;

	size = -1;

	_i[0] = ARC_HEAD(a0);
	_i[1] = ARC_TAIL(a0);
	branch = 0;
	while ( 1 )
	{
		if (!_i[branch]->is_outer) 
		{
			r = _i[branch]; 
			j = _i[1-branch];
			break; 
		}
		_i[branch]->is_outer = 0;
		if (_i[branch]->is_tree_root)
		{
			j = _i[branch];
			i = _i[1-branch];
			while (i->is_outer)
			{
				i->is_outer = 0;
				i = ARC_HEAD(i->match);
				i->is_outer = 0;
				size ++;
				i = ARC_HEAD(i->tree_parent);
			}
			r = i;
			break;
		}
		i = ARC_HEAD(_i[branch]->match);
		i->is_outer = 0;
		size ++;
		_i[branch] = ARC_HEAD(i->tree_parent);
		branch = 1 - branch;
	}
	i = r;
	while ( i != j )
	{
		i = ARC_HEAD(i->match);
		i->is_outer = 1;
		i = ARC_HEAD(i->tree_parent);
		i->is_outer = 1;
		size --;
	}
	if (!clear_is_outer_flag)
	{
		i = ARC_HEAD(a0);
		while (i != r)
		{
			i->is_outer = 1;
			i = ARC_HEAD(i->match);
			i->is_outer = 1;
			i = ARC_HEAD(i->tree_parent);
		}
		i = ARC_TAIL(a0);
		while (i != r)
		{
			i->is_outer = 1;
			i = ARC_HEAD(i->match);
			i->is_outer = 1;
			i = ARC_HEAD(i->tree_parent);
		}
		r->is_outer = 1;
	}
	return r;
}


void PerfectMatching::Shrink(Edge* a0)
{
	assert(a0->head[0]->is_outer && a0->head[1]->is_outer);
	assert(a0->head[0]->flag == 0 && a0->head[1]->flag == 0);

	double start_time = get_time();

	int branch, dir;
	Node* r;
	Node* i;
	Node* j;
	Node** i_ptr;
	Edge* a;
	Edge** a_inner_ptr;
	Arc* a_prev;
	Node* b = blossoms->New();
	int size;
	Edge* a_augment = NULL;
	Edge* b_match;

	b->first[0] = b->first[1] = NULL;

	// set is_outer=0 for all nodes in the blossom
	r = FindBlossomRoot(a0, true, size);
	Tree* t = r->tree;
	REAL eps = t->eps;

	i_ptr = &b->first_tree_child;
	i = ARC_HEAD(a0);
	branch = 0;
	while ( 1 )
	{
		if (i == r && branch) break;
		i->is_marked = 1;
		for (j=i->first_tree_child; j; j=j->tree_sibling)
		{
			if (j->is_outer)
			{
				*i_ptr = j; i_ptr = &j->tree_sibling; 
#ifdef LAZY_CONTRACTION
				Node* k = ARC_HEAD(j->match);
				a = ARC_TO_EDGE_PTR(k->tree_parent);
				dir = ARC_TO_EDGE_DIR(k->tree_parent);
				k = a->head[dir];
				dir ^= 1;
				MOVE_EDGE(k, b, a, dir);
#endif
			}
		}
		if (i == r)
		{
			if (branch) break;
			branch = 1;
			i = ARC_TAIL(a0);
		}
		else 
		{ 
			i = ARC_HEAD(i->match);
			i->is_marked = 1;
			if (i->is_blossom)
			{
				a = ARC_TO_EDGE_PTR(i->match);
				t->pq_blossoms.Remove(a, pq_buf);
				REAL tmp = a->slack; a->slack = i->y; i->y = tmp;
			}
			i = ARC_HEAD(i->tree_parent); 
		}
	}
	*i_ptr = NULL;

	// init b
	b->is_removed = 0;
	b->is_outer = 1;
	b->flag = 0;
	b->is_blossom = 1;
	b->is_tree_root = r->is_tree_root;
	b->is_processed = 1;
	b->tree = t;
	b->tree_sibling = r->tree_sibling;
	b->y = -eps;
	b->is_marked = 0;
	if (b->is_tree_root)
	{
		b->tree_root_prev = r->tree_root_prev;
		if (b->tree_sibling) b->tree_sibling->tree_root_prev = b;
		b->tree_root_prev->tree_sibling = b;
		b->tree->root = b;
		b_match = NULL;
	}
	else
	{
		b->match = r->match;
		i_ptr = & ARC_HEAD(b->match); i_ptr = & ARC_HEAD((*i_ptr)->tree_parent)->first_tree_child;
		while (*i_ptr != r) i_ptr = & ((*i_ptr)->tree_sibling);
		*i_ptr = b;
		b_match = ARC_TO_EDGE_PTR(b->match);
	}
	REAL b_match_slack;
	if (b_match && ARC_HEAD(b->match)->is_blossom)
	{
		b_match_slack = b_match->slack;
		b_match->slack = ARC_HEAD(b->match)->y;
	}

	// second pass over nodes in the blossom
	branch = 0;
	a_prev = EDGE_DIR_TO_ARC(a0, 0);
	i = ARC_HEAD(a_prev);
	while ( 1 )
	{
		// update Arc::next and Arc::head pointers
		if (i->flag == 0) i->y += eps;
		else              i->y -= eps;
		i->is_processed = 0;

#ifdef LAZY_CONTRACTION

		if (i->flag == 1)
		{
			Edge* a_prev;
			for (dir=0; dir<2; dir++)
			if (i->first[dir])
			{
				for (a_inner_ptr=&i->first[dir], a=*a_inner_ptr, a_prev=a->prev[dir], a_prev->next[dir]=NULL; a; a=*a_inner_ptr)
				{
					Node* j0 = a->head[dir];
					for (j=j0; !j->is_outer && !j->is_marked; j = j->blossom_parent) {}
					if (j != j0) { assert(j->flag == 0); int dir_rev = 1 - dir; MOVE_EDGE(j0, j, a, dir_rev); }
					if (j->is_marked) // "inner" arc
					{
						a_inner_ptr = &a->next[dir];
						a->prev[dir] = a_prev;
						a_prev = a;

						if (j->flag == 1) a->slack += eps;
					}
					else // "boundary" arc
					{
						*a_inner_ptr = a->next[dir];
						ADD_EDGE(b, a, dir);

						if (j->flag == 0 && j->tree != t) 
						{
							j->tree->pq_current->pq01[1-j->tree->dir_current].Remove(a, pq_buf);
							if (a->slack + eps <= j->tree->eps) a_augment = a;
						}
						a->slack += 2*eps;
						if (j->flag == 2) t->pq0.Add(a);
						else if (j->flag == 0)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq00.Add(a);
						}
						else if (j->tree != t)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
						}
					}
				}
				if (i->first[dir])
				{
					a_prev->next[dir] = i->first[dir];
					i->first[dir]->prev[dir] = a_prev;
				}
			}
		}

#else // !LAZY_CONTRACTION

		if (i->flag == 1)
		{
			for (dir=0; dir<2; dir++)
			if (i->first[dir])
			{
				for (a_inner_ptr=&i->first[dir], a=*a_inner_ptr; a; a=*a_inner_ptr)
				{
					j = a->head[dir];
					if (j->is_marked) // "inner" arc
					{
						a_inner_ptr = &a->next[dir];

						if (j->flag == 1) a->slack += eps;
					}
					else // "boundary" arc
					{
						*a_inner_ptr = a->next[dir];
						ADD_EDGE(b, a, dir);

						if (j->flag == 0 && j->tree != t) 
						{
							j->tree->pq_current->pq01[1-j->tree->dir_current].Remove(a, pq_buf);
							if (a->slack + eps <= j->tree->eps) a_augment = a;
						}
						a->slack += 2*eps;
						if (j->flag == 2) t->pq0.Add(a);
						else if (j->flag == 0)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq00.Add(a);
						}
						else if (j->tree != t)
						{
							if (!j->tree->pq_current) AddTreeEdge(t, j->tree);
							j->tree->pq_current->pq01[j->tree->dir_current].Add(a);
						}
					}
				}
				*a_inner_ptr = NULL;
			}
		}
		else  // i->flag == 0
		{
			for (dir=0; dir<2; dir++)
			if (i->first[dir])
			{
				for (a_inner_ptr=&i->first[dir], a=*a_inner_ptr; a; a=*a_inner_ptr)
				{
					j = a->head[dir];
					if (j->is_marked) // "inner" arc
					{
						a_inner_ptr = &a->next[dir];

						if (j->flag == 0 && i < j) { t->pq00.Remove(a, pq_buf); a->slack -= 2*eps; }
					}
					else // "boundary" arc
					{
						*a_inner_ptr = a->next[dir];
						ADD_EDGE(b, a, dir);
					}
				}
				*a_inner_ptr = NULL;
			}
		}

#endif // LAZY_CONTRACTION

		Arc* a_next = (i->flag == 0) ? i->match : i->tree_parent;
		i->blossom_parent = b;
		i->match = NULL;
#ifdef USE_GRANDPARENTS
		i->blossom_grandparent = b;
#endif
#ifdef LAZY_CONTRACTION
		i->blossom_selfloops = NULL;
#endif
		if (branch == 0)
		{
			i->blossom_sibling = a_next;
			if (i == r)
			{
				branch = 1;
				a_prev = ARC_REV(a0);
				i = ARC_HEAD(a_prev);
				if (i == r) break;
			}
			else
			{
				a_prev = i->blossom_sibling;
				i = ARC_HEAD(a_prev);
			}
		}
		else
		{
			i->blossom_sibling = ARC_REV(a_prev);
			a_prev = a_next;
			i = ARC_HEAD(a_prev);
			if (i == r) break;
		}
	}
	i->blossom_sibling = ARC_REV(a_prev);
	r->is_tree_root = 0;

	for (i=ARC_HEAD(r->blossom_sibling); ; i = ARC_HEAD(i->blossom_sibling))
	{
		i->is_marked = 0;
		ARC_TO_EDGE_PTR(i->blossom_sibling)->blossom_eps = eps;
		if (i == r) break;
	}

	if (b_match)
	{
		if (ARC_HEAD(b->match)->is_blossom)
		{
			b_match->slack = b_match_slack;
		}
#ifdef LAZY_CONTRACTION
		dir = ARC_TO_EDGE_DIR(b->match);
		assert(b_match->head[1-dir] == r);
		MOVE_EDGE(r, b, b_match, dir);
#endif
	}

	stat.shrink_count ++;
	blossom_num ++;
	stat.shrink_time += get_time() - start_time;

	if (a_augment) Augment(a_augment);
}

