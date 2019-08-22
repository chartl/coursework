#define RLM_MODULE
#include <Python.h>
#include <arrayobject.h>
#include <iostream>
#include <cmath>
#include "npy_common.h"
#include "rlm.h"

using namespace std;

const static double LARGE_VALUE = pow(2., 1020);
const static double LARGE_VALUE_LOG = log(LARGE_VALUE);

static PyMethodDef rlmmethods[] = {
    {"rlm_readlikelihood", rlm_readlikelihood, METH_VARARGS},
    {NULL, NULL}
};


PyMODINIT_FUNC initrlm(void) {
    (void) Py_InitModule("rlm", rlmmethods);
    import_array();
}


static PyObject* rlm_readlikelihood(PyObject* self, PyObject* args) {
    double *match, *ins, *del;
    PyArrayObject *m_ob, *i_ob, *d_ob;
    int n_row, n_col;
    char *read_bases, *haplo_bases, *read_quals;
    if (!PyArg_ParseTuple(args, "OOOsssii", &m_ob, &i_ob, &d_ob, &read_bases,
                          &haplo_bases, &read_quals, &n_row, &n_col)) {
        return NULL;
    }
    match = (double*) m_ob->data;
    ins = (double*) i_ob->data;
    del = (double*) d_ob->data;
    double probability = 0.0;  // total probability
    double lrm, lrd, lri; // last row, same column
    double lrcm, lrcd, lrci; // last row, last column
    double _init_del = LARGE_VALUE / n_row;
    for (int row = 0; row < n_row; row++) {
        for (int col = 0; col < n_col; col++) {
            double bprob = qual2prob(read_quals[col], haplo_bases[row] == read_bases[col]);
            if (row == 0) {
                lrm = 0.0;
                lrd = 0.0;
                lri = 0.0;
            } else {
                lrm = match[col];
                lrd = del[col];
                lri = ins[col];
            }
            if (col == 0) {
                match[col] = _init_del * bprob * DEL2MAT;
                del[col] = lrd * DEL_CONT;
                ins[col] = 0;
            } else {
                match[col] = bprob * ( lrcm * MAT2MAT + lrci * INS2MAT + lrcd * DEL2MAT );
                del[col] = lrm * MAT2DEL + lrd * DEL_CONT;
                ins[col] = match[col-1] * MAT2INS + ins[col-1] * INS_CONT;
            }
            lrcm = lrm;
            lrcd = lrd;
            lrci = lri;
        }
        probability = probability + match[n_col-1] + ins[n_col - 1];
    }
    return Py_BuildValue("d", log(probability) - LARGE_VALUE_LOG);
}


double qual2prob(char qual, bool match) {
    if (match) {
        return 1.0 - qual2prob(qual, false);
    }
    double res = pow(10.0, - ((int) qual - 33)/10.0);
    //cout << qual << " --> " << ((int) qual)-33 << " --> " << res << endl;
    return res;
}

