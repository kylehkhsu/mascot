/*! \file Helper.hh
 *  Contains helper functions.
 */

#ifndef HELPER_HH_
#define HELPER_HH_

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include "SymbolicSet.hh"

using std::cout;
using std::string;

namespace helper {

    /*! Ensures that a specified subdirectory exists.
     *  \param[in]  dirName     Name of the desired subdirectory.
     */
    template <class dir_type>
    void checkMakeDir(dir_type dirName) {
        DIR* dir = opendir(dirName);
        if (dir) {
            closedir(dir);
        }
        else if (ENOENT == errno) {
            int result = mkdir(dirName, 0777);
            (void) result;
        }
    }

    /*! Prints the contents of an array to the console.
     *  \param[in] array    Array.
     *  \param[in] size     Number of elements.
     */
    template<class array_type>
    void printArray (array_type array, int size) {
        for (int i = 0; i < size; i++) {
            cout << array[i] << ' ';
        }
        cout << '\n';
    }

    /*! Prints the contents of a 2D array to the console.
     *  \param[in]  array       2D array.
     *  \param[in]  size1       Number of rows.
     *  \param[in]  size2       Number of columns.
     */
    template<class array_type>
    void printArrayArray (array_type array, int size1, int size2) {
        for (int i = 0; i < size1; i++) {
            printArray(array[i], size2);
        }
        cout << '\n';
    }

    /*! Writes a vector of SymbolicSets into a subdirectory.
     *  \param[in] vec       Vector of SymbolicSets.
     *  \param[in] vecPrefix Subdirectory name + '/' + filename prefix.
     */
    template <class vec_type, class prefix_type>
    void saveVec(vec_type vec, prefix_type vecPrefix) {
        for (size_t i = 0; i < vec.size(); i++) {
            string Str = "";
            Str += vecPrefix;
            Str += std::to_string(i+1);
            Str += ".bdd";
            char Char[20];
            size_t Length = Str.copy(Char, Str.length() + 1);
            Char[Length] = '\0';
            //            cout << Char << '\n';
            vec[i]->writeToFile(Char);
        }
    }

    /*! Frees a vector's objects from heap memory.
     *  \param[in]  vec     Vector of objects.
     */
    template<class vec_type>
    void deleteVec(vec_type vec) {
        for (size_t i = 0; i < vec.size(); i++) {
            delete vec[i];
        }
    }

    /*! Prints a vector's SymbolicSets to the console.
     *  \param[in]  vec         Vector of SymbolicSets.
     *  \param[in]  vecPrefix   Name of a SymbolicSet.
     */
    template <class vec_type, class prefix_type>
    void printVec(vec_type vec, prefix_type vecPrefix) {
        for (size_t i = 0; i < vec.size(); i++) {
            cout << vecPrefix << ' ' << i << ":\n";
            vec[i]->printInfo(1);
        }
    }

    /*! Frees a vector's arrays from heap memory.
     *  \param[in]  vec     Vector of arrays.
     */
    template<class vec_type>
    void deleteVecArray(vec_type vec) {
        for (size_t i = 0; i < vec.size(); i++) {
            delete[] vec[i];
        }
    }
}

#endif /* HELPER_HH_ */
