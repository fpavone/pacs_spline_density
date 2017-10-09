#include<R.h>
#include<Rinternals.h>

extern"C"{
  SEXP add ( SEXP aa, SEXP b) {
    SEXP result = PROTECT ( allocVector ( REALSXP , 1));
    REAL ( result )[0] = asReal (aa) + asReal (b);
    UNPROTECT (1);
    return result ;
    }
}
