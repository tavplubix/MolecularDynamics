#pragma once
#include "Molecule.h"
#include <iterator>

class AbstractContainer;


template <typename T> 
class AbstractIterator :
	public std::iterator<std::forward_iterator_tag, T>
{
	friend class AbstractContainer;
	AbstractIterator(T* p);
	T* p;
public:
	virtual bool operator!=(AbstractIterator const& other) const { return p != other.p; };
	virtual bool operator==(AbstractIterator const& other) const { return p == other.p; };
	virtual typename AbstractIterator::reference operator*() const = 0;
	virtual AbstractIterator& operator++() = 0;
};




class AbstractContainer
{

public:
	typedef AbstractIterator<Molecule> MoleculeIterator;
	AbstractContainer();
	



};




