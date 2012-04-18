#include "lib_svmlight_format.h"

static int32_t next_occurence(char *line, int32_t start, char what)
{
  int32_t i;
  for(i=start; i < LIBSLF_MAXLINELEN && line[i] != '\0'; i++)
  {
    if(line[i] == what)
      return(i);
  }

  return(-1);
}

int32_t svmlight_format_parse_line(char *line, int32_t *label, uint32_t *feat_idx, double *feat_val)
{  
  int32_t beg, end, nnzf=0;

  end = next_occurence(line,0,' ');
  beg = end;
  if(end == -1) 
    return(-1);

  *label = (int32_t)atol(line);

  int go = 1;
  while(go) {
    end = next_occurence(line,beg,':');

    if(end == -1)
      return(nnzf);

    feat_idx[nnzf] = (uint32_t)atol(&line[beg]);

    beg = end + 1;

    end = next_occurence(line,beg,' ');
    if(end == -1) {
      end = next_occurence(line,beg,'\n');
      if(end == -1)
        return(-1);

      go = 0;
    }

    feat_val[nnzf] = atof(&line[beg]);    

    beg = end;

    nnzf++;
  }

  return(nnzf);
}


/* difference to svmlight_format_parse_line is that here the label is float */
int32_t svmlight_format_parse_line_doubley(char *line, double *label, uint32_t *feat_idx, double *feat_val)
{  
  int32_t beg, end, nnzf=0;

  end = next_occurence(line,0,' ');
  beg = end;
  if(end == -1) 
    return(-1);

  *label = (double)atof(line);

  int go = 1;
  while(go) {
    end = next_occurence(line,beg,':');

    if(end == -1)
      return(nnzf);

    feat_idx[nnzf] = (uint32_t)atol(&line[beg]);

    beg = end + 1;

    end = next_occurence(line,beg,' ');
    if(end == -1) {
      end = next_occurence(line,beg,'\n');
      if(end == -1)
        return(-1);

      go = 0;
    }

    feat_val[nnzf] = atof(&line[beg]);    

    beg = end;

    nnzf++;
  }

  return(nnzf);
}
