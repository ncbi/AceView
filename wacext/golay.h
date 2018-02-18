int golay_decode(int correct_mode, int *errs, unsigned long *cw) ;
/* This function decodes codeword *cw in one of two modes. If correct_mode
   is nonzero, error correction is attempted, with *errs set to the number of
   bits corrected, and returning 0 if no errors exist, or 1 if parity errors
   exist. If correct_mode is zero, error detection is performed on *cw,
   returning 0 if no errors exist, 1 if an overall parity error exists, and
   2 if a codeword error exists. */
