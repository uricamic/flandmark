#ifndef lib_svmlight_format_h
#define lib_svmlight_format_h

#include <stdlib.h>
#include <stdint.h>

#define LIBSLF_MAXLINELEN 1000000

int32_t svmlight_format_parse_line(char *line, int32_t *label, uint32_t *feat_idx, double *feat_val);
int32_t svmlight_format_parse_line_doubley(char *line, double *label, uint32_t *feat_idx, double *feat_val);


#endif
