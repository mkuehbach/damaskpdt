/*
	Copyright Max-Planck-Institut f\"uer Eisenforschung, GmbH, D\"usseldorf
	Programmer: Markus K\"uhbach, 2017

	damaskpdt --- the DAMASK parallelized property distribution tracker --- is an MPI/OpenMP/SIMD-
	parallelized tool for datamining the spatial distribution of state variable values at interfaces
	which evolve during straining in 3d full-field crystal plasticity simulations.

	It supplements the python post-processing routines of the DAMASK crystal plasticity solver
	package <https://damask.mpie.de/>. The novelty is the rigorous parallel and locality-aware
	spatial-distributed storing of snapshot data to make the processing of terabyte-scale
	crystal plasticity simulation snapshot data efficient, scalable.

	This enables to study and quantify distributional trends, which other techniques and
	analysis workflows have so far not practically been able to process.

	The source code was developed by Markus K\"uhbach during his postdoc with the Max-Planck-Institut f\"ur
	Eisenforschung GmbH in D\"usseldorf for the purpose of quantifying the potential accumulation of dislocation
	density, stresses, and strain in front of interfaces. The model implements a distance based approach
	m.kuehbach at mpie.de

	The authors gratefully acknowledge the financial support from the Deutsche Forschungsgemeinschaft
	(DFG) within the project "Consistent physically-based modeling of dynamic recrystallization
	under hot working conditions" project number 315419526.

	This file is part of damaskpdt.

	damaskpdt is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	damask_pdt is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with damask_pdt.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __PDT_H_AABBTREE_INTERFACE_H__
#define __PDT_H_AABBTREE_INTERFACE_H__

#include "PDT_Math.h"

/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2017 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/


/// Null node flag
#define NULL_NODE				(UINT32MX)
#define TREE_ERROR				(UINT32MX)

#ifdef EMPLOY_SINGLEPRECISION
	#define ASSERTATION_ACCURACY	(1.0e-6)
#else
	#define ASSERTATION_ACCURACY	(1.0e-12)
#endif

//MK::modifications considering that the voro_real world APT tip is 3d and does not show periodicity
//Axis-aligned bounding boxes (AABBs) store information for the minimum orthorhombic bounding-box for a triangle

struct trpl
{
	//a triple of voro_real numbers
	voro_real a;
	voro_real b;
	voro_real c;
	trpl() : a(0.f), b(0.f), c(0.f) {}
	trpl(const voro_real _a, const voro_real _b, const voro_real _c ) : a(_a), b(_b), c(_c) {}
	//~trpl(){}
};


class hedgesAABB
{
public:
	hedgesAABB() : surfaceArea(0.f) {
		//trpl objects get values assigned via their default constructor
	};
	hedgesAABB(const trpl lowerBound_, const trpl upperBound_);
	~hedgesAABB(){};

	voro_real computeSurfaceArea() const;
	voro_real getSurfaceArea() const;
	void merge(const hedgesAABB&, const hedgesAABB&);
	bool contains(const hedgesAABB&) const;
	bool overlaps(const hedgesAABB&) const;
	trpl computeCentre();
	
	//variables
	trpl lowerBound;
	trpl upperBound;
	trpl centre;
	voro_real surfaceArea;
};

/*! \brief A node of the hedgesAABB tree.

    Each node of the tree contains an hedgesAABB object which corresponds to a
    particle, or a group of particles, in the simulation box. The AABB
    objects of individual particles are "fattened" before they are stored
    to avoid having to continually update and rebalance the tree when
    displacements are small.

    Nodes are aware of their position within in the tree. The isLeaf member
    function allows the tree to query whether the node is a leaf, i.e. to
    determine whether it holds a single particle.
*/
struct Node
{
	Node() : parent(NULL_NODE), next(NULL_NODE), left(NULL_NODE), right(NULL_NODE), height(0), particle(0) {
		//##MK::add default constructor choices
	};
	hedgesAABB aabb;
	//node dependencies, parent and next
	unsigned int parent;
	unsigned int next;

	//indices of left and right child nodes
	unsigned int left;
	unsigned int right;

	/// Height of the node. This is 0 for a leaf and -1 for a free node.
	int height;

	/// The index of the particle that the node contains (leaf nodes only).
	unsigned int particle;

	//! Test whether the node is a leaf.

	bool isLeaf() const;
};

/*! \brief The dynamic hedgesAABB tree.

    The dynamic hedgesAABB tree is a hierarchical data structure that can be used
    to efficiently query overlaps between objects of arbitrary shape and
    size that lie inside of a simulation box. Support is provided for
    periodic and non-periodic boxes, as well as boxes with partial
    periodicity, e.g. periodic along specific axes.
	//MK::periodicty was eliminated
    */
class Tree
{
public:
	Tree(unsigned int nParticles = 16);
	//Tree(const trpl boxSize_, unsigned int nParticles = 16);
	//void setBoxSize(const trpl boxSize_);

	//! Insert a particle into the tree (point particle) with index, position and radius
	void insertParticle(unsigned int, const trpl, voro_real);

	//! Insert a particle into the tree (arbitrary shape with bounding box) with index, lower and upper
	void insertParticle(unsigned int, const trpl, const trpl);

	//! Remove a particle from the tree.
	void removeParticle(unsigned int);

	//MK::no updates of tree required

	//! Query the tree to find candidate interactions for a particle.
	std::vector<unsigned int> query(unsigned int);

	//! Query the tree to find candidate interactions for an hedgesAABB, in: particle index, out:: aabb, particles indices
	std::vector<unsigned int> query(unsigned int, const hedgesAABB&);

	//! Query the tree to find candidate interactions for an hedgesAABB, in an hedgesAABB, out particles indices
	std::vector<unsigned int> query(const hedgesAABB&);

	//! Get a particle hedgesAABB for particle index
	const hedgesAABB& getAABB(unsigned int);

	//! Get the height of the tree.
	unsigned int getHeight() const;

	//! Get the number of nodes in the tree
	unsigned int getNodeCount() const;

	//! Compute the maximum balance of the tree.
	unsigned int computeMaximumBalance() const;

	//! Compute the surface area ratio of the tree.
	voro_real computeSurfaceAreaRatio() const;

	//! Validate the tree.
	void validate() const;

	void report_tree( const string csv_fn );

private:
	/// The index of the root node.
	unsigned int root;

	/// The dynamic tree. //##MK::place better on heap...
	vector<Node> nodes;

	/// The current number of nodes in the tree.
	unsigned int nodeCount;

	/// The current node capacity.
	unsigned int nodeCapacity;

	/// The position of node at the top of the free list.
	unsigned int freeList;

	/// The size of the system in each dimension.
	//trpl boxSize;

	/// A map between particle and node indices.
	map<unsigned int, unsigned int> particleMap;

	//! Allocate a new node, returns index of the allocated node
	unsigned int allocateNode();

	//! Free an existing node, index of the node to be freed
	void freeNode(unsigned int);
	//unsigned int return

	//! Insert a leaf into the tree, index of the leaf node
	unsigned int insertLeaf(unsigned int);

	//! Remove a leaf from the tree, index of the leaf node
	void removeLeaf(unsigned int);

	//! Balance the tree, index of the node
	unsigned int balance(unsigned int);

	//! Compute the height of the tree, height of the entire tree.
	unsigned int computeHeight() const;

	//! Compute the height of a sub-tree, in: index of rootNode, out: height of sub-tree
	unsigned int computeHeight(unsigned int) const;

	//! Assert that the sub-tree has a valid structure.
	void validateStructure(unsigned int) const;

	//! Assert that the sub-tree has valid metrics.
	void validateMetrics(unsigned int) const;
};


#endif
