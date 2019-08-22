static PyObject *banded_sw_c(PyObject *dummy, PyObject *args);
int max3_(int a, int b, int c);
void set_(int* ar, int nr, int nc, int r, int c, int v);
int get_(int* ar, int nr, int nc, int r, int c);

const static int MATCH_GAIN = 15;
const static int MISMATCH_COST = 9;
const static int DEL_GAP_OPEN_COST = 20;
const static int INS_GAP_OPEN_COST = 9;
const static int GAP_EXTENSION_COST = 3;
const static int BAND_START_BEST = 100;
const static int BAND_DIFF = 80;
