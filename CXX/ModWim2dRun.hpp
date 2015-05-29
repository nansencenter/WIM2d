/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t  -*- */


#ifndef __MODWIM2DRUN_H
#define __MODWIM2DRUN_H 1

/**
 *
 */
template<class T> class ModWim2dRun
{
public:

	/**
	 *
	 */
	typedef T value_type;

	/**
	 *
	 */
	typedef value_type &reference;

	/**
	 *
	 */
	typedef const reference const_reference;

	/**
	 *
	 */
	typedef unsigned int size_type;

	/**
	 * Creates a new circular buffer.
	 *
	 * @param capacity The number of items the buffer can contain (strictly
	 *                 greater than 1).
	 */
	circular_buffer(size_type capacity);

	/**
	 *
	 */
	~circular_buffer();

	/**
	 *
	 */
	inline reference at(size_type index)
	{
		// Reuse the implementation of at(size_type) const.
		return (reference) this->at(index);
	}

	/**
	 *
	 */
	inline const_reference at(size_type index) const
	{
		if (this->empty())
		{
			throw std::out_of_range("No such index");
		}

		// Reuse the implementation of operator[](size_type) const.
		return this[index];
	}

	/**
	 *
	 */
	inline size_type capacity() const
	{
		return this->_capacity;
	}

	/**
	 *
	 */
	void clear();

	/**
	 *
	 */
	inline size_type empty() const
	{
		return (this->size() == 0);
	}

	/**
	 *
	 */
	inline size_type full() const
	{
		return (this->size() == this->capacity());
	}

	/**
	 *
	 */
	inline reference operator[](size_type index)
	{
		// Reuse the implementation of operator[](size_type) const.
		return (reference) this[index];
	}

	/**
	 *
	 */
	inline const_reference operator[](size_type index) const
	{
		requires(index < this->size());

		index = (index + this->_start) % this->capacity();

		assert(index < this->capacity());

		return this->_buffer[index];
	}

	/**
	 *
	 */
	void pop_back();

	/**
	 *
	 */
	void pop_front();

	/**
	 *
	 */
	void push_back(const_reference item = value_type());

	/**
	 *
	 */
	void push_front(const_reference item = value_type());

	/**
	 *
	 */
	inline size_type size() const
	{
		return this->_size;
	}

private:

	/**
	 * This array contains the data.
	 */
	T *_buffer;

	/**
	 * The capacity of _buffer, i.e. the number of items which can be contained.
	 */
	size_type _capacity;

	/**
	 * The number of items currently in _buffer.
	 */
	size_type _size;

	/**
	 * The index in _buffer from which the items are placed.
	 */
	size_type _start;
};

#endif
