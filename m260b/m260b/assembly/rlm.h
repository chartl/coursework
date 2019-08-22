static PyObject* rlm_readlikelihood(PyObject* self, PyObject* args);
double qual2prob(char qual, bool match);
const static double MAT2MAT = 0.999;
const static double DEL2MAT = 0.95;
const static double INS2MAT = 0.95;
const static double MAT2DEL = 0.0005;
const static double MAT2INS = 0.0005;
const static double DEL_CONT = 0.05;
const static double INS_CONT = 0.05;
