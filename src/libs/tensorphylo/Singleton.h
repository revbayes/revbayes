//    HOGAN is an implementation of a parallel Metropolis-Hastings algorithm
//    developped for evolutionnary biology model.
//    Copyright (C) 2016  Xavier Meyer
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
/*
 * Singleton.h
 *
 *  Created on: Nov 30, 2011
 *      Author: meyerx
 */

#ifndef SINGLETON_H_
#define SINGLETON_H_

#include <cstdio>

namespace Utils {
template <typename T>
class Singleton {

protected:
	// Constructor / Destructor
	Singleton() { }
	~Singleton() { }

public:

	// Function creating if needed the singleton and returning a pointer to it
	static T *getInstance() {
		if (_singleton == NULL){
			_singleton = new T;
		}
		return (static_cast<T*> (_singleton));
	}

	static void kill() {
		if (_singleton != NULL){
			delete _singleton;
			_singleton = NULL;
		}
	}

private:
	// Unique instance
	static T *_singleton;
};

template <typename T>
T *Singleton<T>::_singleton = NULL;

} // Utils

#endif /* SINGLETON_H_ */
