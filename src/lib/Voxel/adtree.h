#include <mystdlib.h>


#include <myadt.hpp>
// class DenseMatrix;
#include <gprim.hpp>

namespace netgen
{


  /* ******************************* ADTree ******************************* */


  void ADTree3 :: GetIntersecting (const float * bmin, 
				   const float * bmax,
				   Array<int> & pis) const
  {
    ArrayMem<ADTreeNode3*,1001> stack(1000);
    ArrayMem<int,1001> stackdir(1000);
    ADTreeNode3 * node;
    int dir, stacks;

    stack.SetSize (1000);
    stackdir.SetSize(1000);
    pis.SetSize(0);

    stack.Elem(1) = root;
    stackdir.Elem(1) = 0;
    stacks = 1;

    while (stacks)
      {
	node = stack.Get(stacks);
	dir = stackdir.Get(stacks); 
	stacks--;

	if (node->pi != -1)
	  {
	    if (node->data[0] >= bmin[0] && node->data[0] <= bmax[0] &&
		node->data[1] >= bmin[1] && node->data[1] <= bmax[1] &&
		node->data[2] >= bmin[2] && node->data[2] <= bmax[2])

	      pis.Append (node->pi);
	  }


	int ndir = dir+1;
	if (ndir == 3)
	  ndir = 0;

	if (node->left && bmin[dir] <= node->sep)
	  {
	    stacks++;
	    stack.Elem(stacks) = node->left;
	    stackdir.Elem(stacks) = ndir;
	  }
	if (node->right && bmax[dir] >= node->sep)
	  {
	    stacks++;
	    stack.Elem(stacks) = node->right;
	    stackdir.Elem(stacks) = ndir;
	  }
      }
  }

}
