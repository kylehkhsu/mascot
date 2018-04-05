/*! \file Helper.hh
 *  Contains helper functions.
 */

#ifndef HELPER_HH_
#define HELPER_HH_

#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
//#include "SymbolicSet.hh"

using std::clog;
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

    /*! Prints the contents of an array to a log file.
     *  \param[in] array    Array.
     *  \param[in] size     Number of elements.
     */
    template<class array_type>
    void printArray (array_type array, int size) {
        for (int i = 0; i < size; i++) {
            clog << array[i] << ' ';
        }
        clog << '\n';
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
            //vec[i]->writeToFile(Char);
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

    template<class vecVec_type>
    void deleteVecVec(vecVec_type vecVec) {
        for (size_t i = 0; i < vecVec.size(); i++) {
            deleteVec(*vecVec[i]);
            delete vecVec[i];
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

    template<class vec_t, class prefix_t>
    void printVecArray(vec_t vec, prefix_t prefix, int arrayElems) {
        for (size_t i = 0; i < vec.size(); i++) {
            clog << prefix << "[" << i << "]: ";
            printArray(vec[i], arrayElems);
        }
    }

    template <class vecVec_type, class prefix_type>
    void printVecVec(vecVec_type vecVec, prefix_type vecPrefix) {
        for (size_t i = 0; i < vecVec.size(); i++) {
            for (size_t j = 0; j < vecVec[i]->size(); j++) {
                cout << vecPrefix << '[' << i << "][" << j << "]:\n";
                (*vecVec[i])[j]->printInfo(1);
            }
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

    template<class vecVec_type>
    void deleteVecVecArray(vecVec_type vecVec) {
        for (size_t i = 0; i < vecVec.size(); i++) {
            deleteVecArray(*vecVec[i]);
            delete vecVec[i];
        }
    }

    /*! Throws a logic error along with specified message and closes the log file.
     *  \param[in]  msg     Error message to log to file.
     */
    template<class msg_type>
    void error(msg_type msg) {
        std::ostringstream os;
        os << msg;
        throw std::logic_error(os.str().c_str());
        fclose(stderr);
    }
}

#endif /* HELPER_HH_ */
