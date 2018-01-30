/*
    PQ.h - implements Fibonacci priority queue

    Copyright 2008 Vladimir Kolmogorov (vnk@adastral.ucl.ac.uk)

    This software can be used for research purposes only. Commercial use is prohibited.
    Public redistribution of the code or its derivatives is prohibited.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef HFKSJHFKJHARBABDAKFAF
#define HFKSJHFKJHARBABDAKFAF

#include <assert.h>
#include <string.h>

template <typename REAL> class PriorityQueue
{
public:
	struct Item
	{
		REAL	slack;

		Item*	prevPQ;
		union
		{
			struct
			{
				Item*	nextPQ;
				Item*	parentPQ;
				Item*	childPQ;
				unsigned short degreePQ;
				unsigned short markPQ;
			};
			struct
			{
				REAL	blossom_eps;
				REAL	y_saved;
			};
		};
	};
	static void* AllocateBuf();
	static void DeallocateBuf(void* buf);

	static void ResetItem(Item* i);
	static bool isReset(Item* i);

	//////////////////////////////////////////////////////////

	void Reset();
	void Add(Item* i);
	void Remove(Item* i, void* buf);
	void Decrease(Item* i_old, Item* i_new, void* buf);
	Item* GetMin();

	//////////////////////////////////////////////////////////

	void Update(REAL delta);
	void Merge(PriorityQueue<REAL>& dest);

	// traversing items in the order they are stored (irrespective of slack).
	// The caller must go through all items, no other member functions can be called during the scan.
	Item* GetAndResetFirst();
	Item* GetAndResetNext();

	Item* GetFirst();
	Item* GetNext(Item* i);

	//////////////////////////////////////////////////////////

private:
	struct Buf
	{
		Item**	array;
		int		array_size;
	};
	Item*	rootPQ;
	void Merge(Buf* buf);
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline void* PriorityQueue<REAL>::AllocateBuf()
{
	Buf* buf = new Buf();
	buf->array_size = 16;
	buf->array = (Item**) malloc(buf->array_size*sizeof(Item*));
	memset(buf->array, 0, buf->array_size*sizeof(Item*));
	return buf;
}

template <typename REAL> inline void PriorityQueue<REAL>::DeallocateBuf(void* _buf)
{
	Buf* buf = (Buf*)_buf;
	free(buf->array);
	delete buf;
}

template <typename REAL> inline void PriorityQueue<REAL>::ResetItem(Item* i) 
{ 
	i->prevPQ = NULL;
}

template <typename REAL> inline bool PriorityQueue<REAL>::isReset(Item* i) 
{ 
	return (i->prevPQ == NULL);
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline void PriorityQueue<REAL>::Reset() 
{ 
	rootPQ = NULL; 
}

template <typename REAL> inline void PriorityQueue<REAL>::Add(Item* i)
{
	if (!rootPQ)
	{
		rootPQ = i;
		i->prevPQ = i->nextPQ = i;
	}
	else
	{
		i->nextPQ = rootPQ;
		i->prevPQ = rootPQ->prevPQ;
		rootPQ->prevPQ->nextPQ = i;
		rootPQ->prevPQ = i;
		if (rootPQ->slack > i->slack) rootPQ = i;
	}
	i->parentPQ = NULL;
	i->childPQ = NULL;
	i->degreePQ = 0;
	i->markPQ = 0;
}

template <typename REAL> inline void PriorityQueue<REAL>::Merge(Buf* buf)
{
	if (!rootPQ) return;
	rootPQ->parentPQ = NULL;
	if (rootPQ->nextPQ == rootPQ) return;

	Item* i = rootPQ;
	Item* next = i->nextPQ;
	Item* j;

	i->prevPQ->nextPQ = NULL;
	i->prevPQ = i->nextPQ = i;

	while ( 1 )
	{
		//int d = 0;
		//if (i->childPQ) { for (d++, j=i->childPQ->nextPQ; j!=i->childPQ; j=j->nextPQ) d++; }
		//assert(d == i->degreePQ);

		if (i->degreePQ >= buf->array_size)
		{
			int size_old = buf->array_size;
			buf->array_size = i->degreePQ + 1 + i->degreePQ/2;
			buf->array = (Item**)realloc(buf->array, buf->array_size*sizeof(Item*));
			memset(buf->array+size_old, 0, (buf->array_size-size_old)*sizeof(Item*));
		}
		if (buf->array[i->degreePQ] == NULL)
		{
			buf->array[i->degreePQ] = i;
			i->parentPQ = NULL;
			if (!next) break;
			i = next;
			next = i->nextPQ;
			i->nextPQ = rootPQ;
			i->prevPQ = rootPQ->prevPQ;
			rootPQ->prevPQ->nextPQ = i;
			rootPQ->prevPQ = i;
			if (rootPQ->slack > i->slack) rootPQ = i;
			continue;
		}
		if (buf->array[i->degreePQ]->slack >= i->slack) j = buf->array[i->degreePQ];
		else { j = i; i = buf->array[i->degreePQ]; }
		if (rootPQ == j) rootPQ = i;
		buf->array[i->degreePQ] = NULL;
		i->degreePQ ++;
		j->prevPQ->nextPQ = j->nextPQ;
		j->nextPQ->prevPQ = j->prevPQ;
		j->parentPQ = i;
		if (i->childPQ)
		{
			j->nextPQ = i->childPQ;
			j->prevPQ = i->childPQ->prevPQ;
			i->childPQ->prevPQ->nextPQ = j;
			i->childPQ->prevPQ = j;
		}
		else
		{
			i->childPQ = j;
			j->nextPQ = j->prevPQ = j;
		}
	}

	for (i=rootPQ->nextPQ; ; i=i->nextPQ) { buf->array[i->degreePQ] = NULL; if (i == rootPQ) break; }
}

template <typename REAL> inline void PriorityQueue<REAL>::Remove(Item* i, void* _buf)
{
	Item* i0 = i;
	if (i->parentPQ)
	while ( 1 )
	{
		// cut i
		Item* parent = i->parentPQ;
		if (i->prevPQ == i) parent->childPQ = NULL;
		else
		{
			i->prevPQ->nextPQ = i->nextPQ;
			i->nextPQ->prevPQ = i->prevPQ;
			parent->childPQ = i->nextPQ;
		}
		i->nextPQ = rootPQ;
		i->prevPQ = rootPQ->prevPQ;
		rootPQ->prevPQ->nextPQ = i;
		rootPQ->prevPQ = i;
		i->markPQ = 0;

		i = parent;
		i->degreePQ --;
		if (!i->parentPQ) break;
		if (i->markPQ == 0) { i->markPQ = 1; break; }
	}
	i = i0;

	Item* first_child = i->childPQ;
	Item* sibling;

	if (i->prevPQ == i) sibling = NULL;
	else
	{
		sibling = i->prevPQ;
		sibling->nextPQ = i->nextPQ;
		i->nextPQ->prevPQ = sibling;
	}

	if (i == rootPQ) rootPQ = sibling;

	if (first_child)
	{
		Item* last_child = first_child->prevPQ;
		if (rootPQ)
		{
			first_child->prevPQ = rootPQ;
			last_child->nextPQ = rootPQ->nextPQ;
			rootPQ->nextPQ->prevPQ = last_child;
			rootPQ->nextPQ = first_child;
		}
		else
		{
			first_child->prevPQ = last_child;
			last_child->nextPQ = first_child;
			rootPQ = first_child;
		}
	}
	Merge((Buf*)_buf);
	i->prevPQ = NULL;
}

template <typename REAL> inline void PriorityQueue<REAL>::Decrease(Item* i_old, Item* i_new, void* _buf)
{
	Item* i;

	i_new->parentPQ = i_old->parentPQ;
	i_new->childPQ = i_old->childPQ;
	if (i_old->childPQ) 
	{
		for (i=i_old->childPQ->nextPQ; ; i=i->nextPQ)
		{
			i->parentPQ = i_new;
			if (i == i_old->childPQ) break;
		}
	}
	if (i_old->prevPQ == i_old) i_new->prevPQ = i_new->nextPQ = i_new;
	else
	{
		i_old->prevPQ->nextPQ = i_new;
		i_old->nextPQ->prevPQ = i_new;
		i_new->prevPQ = i_old->prevPQ;
		i_new->nextPQ = i_old->nextPQ;
	}

	if (!i_old->parentPQ || i_old->parentPQ->slack <= i_new->slack)
	{
		if (i_old->parentPQ) i_old->parentPQ->childPQ = i_new;
		else if (rootPQ->slack > i_new->slack) rootPQ = i_new;
	}
	else
	{
		i = i_new;
		while ( 1 )
		{
			// cut i
			Item* parent = i->parentPQ;
			if (i->prevPQ == i) parent->childPQ = NULL;
			else
			{
				i->prevPQ->nextPQ = i->nextPQ;
				i->nextPQ->prevPQ = i->prevPQ;
				parent->childPQ = i->nextPQ;
			}
			i->nextPQ = rootPQ;
			i->prevPQ = rootPQ->prevPQ;
			rootPQ->prevPQ->nextPQ = i;
			rootPQ->prevPQ = i;
			i->markPQ = 0;

			i = parent;
			i->degreePQ --;
			if (!i->parentPQ) break;
			if (i->markPQ == 0) { i->markPQ = 1; break; }
		}
		Merge((Buf*)_buf);
	}
	i_old->prevPQ = NULL;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetMin()
{
	return rootPQ;
}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



template <typename REAL> inline void PriorityQueue<REAL>::Merge(PriorityQueue<REAL>& dest)
{
	if (!rootPQ) return;
	if (!dest.rootPQ) dest.rootPQ = rootPQ;
	else
	{
		Item* last = rootPQ->prevPQ;
		rootPQ->prevPQ = dest.rootPQ->prevPQ;
		last->nextPQ = dest.rootPQ;
		dest.rootPQ->prevPQ->nextPQ = rootPQ;
		dest.rootPQ->prevPQ = last;
		if (dest.rootPQ->slack > rootPQ->slack) dest.rootPQ = rootPQ;
	}
	rootPQ = NULL;
}



template <typename REAL> inline void PriorityQueue<REAL>::Update(REAL delta)
{
	if (!rootPQ) return;

	Item* i = rootPQ->nextPQ;
	while ( 1 )
	{
		i->slack += delta;

		if (i->childPQ) i = i->childPQ->nextPQ;
		else
		{
			while ( i->parentPQ && i == i->parentPQ->childPQ ) i = i->parentPQ;
			if (i == rootPQ) break;
			i = i->nextPQ;
		}
	}
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetAndResetFirst()
{
	if (rootPQ) rootPQ->prevPQ->nextPQ = NULL;
	return GetAndResetNext();
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetAndResetNext()
{
	if (!rootPQ) return NULL;
	Item* i = rootPQ;
	if (i->childPQ)
	{
		rootPQ = i->childPQ;
		rootPQ->prevPQ->nextPQ = i->nextPQ;
	}
	else rootPQ = i->nextPQ;
	i->prevPQ = NULL;
	return i;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetFirst()
{
	if (!rootPQ) return NULL;
	return rootPQ->nextPQ;
}

template <typename REAL> inline typename PriorityQueue<REAL>::Item* PriorityQueue<REAL>::GetNext(Item* i)
{
	if (i->childPQ) return i->childPQ->nextPQ;
	else
	{
		while ( i->parentPQ && i == i->parentPQ->childPQ ) i = i->parentPQ;
		if (i == rootPQ) return NULL;
		return i->nextPQ;
	}
}

#endif
