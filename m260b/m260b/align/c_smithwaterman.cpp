#define C_SMITHWATERMAN_MODULE
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// actually if we do this, we lose the ability to `release' the data
// (e.g. access to the flag attribute of the ndarray. Not good.)
#include <Python.h>
#include <arrayobject.h>
#include "npy_common.h"
#include "c_smithwaterman.h"
#include <iostream>
using namespace std;;

static PyMethodDef swmethods[] = {
    {"banded_sw_c", banded_sw_c, METH_VARARGS},
    {NULL, NULL} // has to be here
};


PyMODINIT_FUNC initc_smithwaterman(void) {
    (void) Py_InitModule("c_smithwaterman", swmethods);
    import_array();
}

int max3_(int a, int b, int c) {
    if (a > b) {
        if (a > c) {
            return a;
        }
        return c;
    } else {
        if (b > c) {
            return b;
        }
        return c;
    }
}

void set_(int* ar, int nr, int nc, int r, int c, int v) {
    ar[r * nc + c] = v;
}

int get_(int* ar, int nr, int nc, int r, int c) {
    return ar[r*nc + c];
}

static PyObject* banded_sw_c(PyObject *dummy, PyObject *args) {
    const char *ref, *read;
    PyArrayObject *np_score, *np_moves;
    int refsize, readsize, off_edge_penalty;
    if (!PyArg_ParseTuple(args, "ssiiiOO", &ref, &read, &refsize, &readsize,
                          &off_edge_penalty, &np_score, &np_moves)) {
        return NULL;
    }
    /*cout << "Parsed reference:" << endl;
    for(int i = 0; i < refsize; i++){
        cout << ref[i];
    }*/
    // create the score matrix and move matrix. 
    // note that because we're allocating these on the heap, they have to be
    // one-dimensional, but indexed via a stride
    //cout << "Set dims to " << dims[0] << " and " << dims[1] << endl;
    int *swscore = (int*) np_score->data;
    int *moves = (int*) np_moves->data;
    for (int r = 0; r < refsize + 1; r++) {
        for (int c = 0; c < readsize + 1; c ++) {
            set_(swscore, 1+refsize, 1+readsize, r, c, 0);
            set_(moves, 1+refsize, 1+readsize, r, c, 0);
        }
    }
    for (int c=1; c < readsize + 1; c++) {
        set_(swscore, 1+refsize, 1+readsize, 0, c, -off_edge_penalty);
    }
    /*cout << endl << "Initialized arrays" << endl;
    for (int r = 0; r < refsize + 1; r ++) {
        cout << endl;
        for (int c = 0; c < readsize + 1; c++) {
            cout << " ";
            cout << get_(swscore, 1+refsize, 1+readsize, r, c); 
        }
    }
    cout << endl;*/
    // unlike the implementation of smith waterman in python, here the matrix is
    // propagated by column rather than by row, and the band limits are placed
    // on the rows. This is because there may be multiple good starting locations
    // on the reference, so we need to wait until one is clearly better before
    // banding it away.
    int nrow = 1 + refsize;
    int ncol = 1 + readsize;
    int rowband_start = 1;
    int rowband_end = refsize + 1;
    int mvalue, dvalue, ivalue, bestvalue = -1;
    for (int col = 1; col < ncol; col++) {
	int colbest=-1;
        for(int row=rowband_start; row < rowband_end; row++) {
            // look in the header for these constants
            if (ref[row-1] == read[col-1]) {
                mvalue = get_(swscore, nrow, ncol, row-1, col-1) + MATCH_GAIN;
            } else {
                mvalue = get_(swscore, nrow, ncol, row-1, col-1) - MISMATCH_COST;
            }
            if (get_(moves, nrow, ncol, row-1, col) == 1) {
                dvalue = get_(swscore, nrow, ncol, row-1, col) - GAP_EXTENSION_COST;
            } else {
                dvalue = get_(swscore, nrow, ncol, row-1, col) - DEL_GAP_OPEN_COST;
            }
            if (get_(moves, nrow, ncol, row, col-1) == -1) {
                ivalue = get_(swscore, nrow, ncol, row, col-1) - GAP_EXTENSION_COST;
            } else {
                ivalue = get_(swscore, nrow, ncol, row, col-1) - INS_GAP_OPEN_COST;
            }
            bestvalue = max3_(mvalue, dvalue, ivalue);
            set_(swscore, nrow, ncol, row, col, bestvalue);
            if (ivalue == bestvalue) {
                set_(moves, nrow, ncol, row, col, -1);
            } else if (dvalue == bestvalue) {
                set_(moves, nrow, ncol, row, col, 1);
            }
            if (bestvalue > colbest) {
                colbest = bestvalue;
            }
            //cout << " " << bestvalue << endl;
        }
        if (colbest > BAND_START_BEST) {
            while(get_(swscore, nrow, ncol, rowband_start, col) < colbest - BAND_DIFF) {
                ++rowband_start;
            }
            while(get_(swscore, nrow, ncol, rowband_end-1, col) < colbest - BAND_DIFF) {
                --rowband_end;
            }
            if (rowband_end < refsize + 1) {
                ++rowband_end;  // allow a match on the diagonal
            }
        }
    }
    /*cout << endl;
    for (int r = 0; r < 1 + refsize; r++) {
        cout << endl; 
        for (int c = 0; c < 1 + readsize; c++) {
            cout << " ";
            cout << get_(swscore, nrow, ncol, r, c);
        }
    }
    cout << "Making numpy arrays" << endl;*/
    // construct numpy objects around the matrices
    //cout << "arrays reated. Trying to set flags." << endl;
    //cout << "building the return object" << endl;
    PyObject* ret = Py_BuildValue("OO", np_score, np_moves);
    //cout << "now freeing memory" << endl;
    //cout << "returning" << endl;
    return ret;
};
