/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * Written (W) 2012 Vojtech Franc, Michal Uricar
 * Copyright (C) 2012 Vojtech Franc, Michal Uricar
 */

#include "liblbp.h"

/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------*/
void liblbp_pyr_features_sparse(t_index* vec, uint32_t vec_nDim, uint32_t* img, uint16_t img_nRows, uint16_t img_nCols)
{
    uint32_t offset, ww, hh, x, y, center, j, idx;
    uint8_t pattern;

    idx = 0;
    offset = 0;
    ww = img_nCols;
    hh = img_nRows;
    while(1)
    {
        for(x = 1; x < ww-1; x++)
        {
            for(y = 1; y< hh-1; y++)
            {
                pattern = 0;
                center = img[LIBLBP_INDEX(y,x,img_nRows)];
                if(img[LIBLBP_INDEX(y-1,x-1,img_nRows)] < center) pattern = pattern | 0x01;
                if(img[LIBLBP_INDEX(y-1,x,img_nRows)] < center)   pattern = pattern | 0x02;
                if(img[LIBLBP_INDEX(y-1,x+1,img_nRows)] < center) pattern = pattern | 0x04;
                if(img[LIBLBP_INDEX(y,x-1,img_nRows)] < center)   pattern = pattern | 0x08;
                if(img[LIBLBP_INDEX(y,x+1,img_nRows)] < center)   pattern = pattern | 0x10;
                if(img[LIBLBP_INDEX(y+1,x-1,img_nRows)] < center) pattern = pattern | 0x20;
                if(img[LIBLBP_INDEX(y+1,x,img_nRows)] < center)   pattern = pattern | 0x40;
                if(img[LIBLBP_INDEX(y+1,x+1,img_nRows)] < center) pattern = pattern | 0x80;

                //vec[offset+pattern]++;
                vec[idx++] = offset+pattern;
                offset += 256;
            }
        }
        if(vec_nDim <= idx)
          return;

        if(ww % 2 == 1) ww--;
        if(hh % 2 == 1) hh--;

        ww = ww/2;
        
        for(x=0; x < ww; x++)
          for(j=0; j < hh; j++)
            img[LIBLBP_INDEX(j,x,img_nRows)] = img[LIBLBP_INDEX(j,2*x,img_nRows)] + 
              img[LIBLBP_INDEX(j,2*x+1,img_nRows)];

        hh = hh/2;
        
        for(y=0; y < hh; y++)
          for(j=0; j < ww; j++)
            img[LIBLBP_INDEX(y,j,img_nRows)] = img[LIBLBP_INDEX(2*y,j,img_nRows)] + 
              img[LIBLBP_INDEX(2*y+1,j,img_nRows)];
    }
    return;
}

/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------*/
void liblbp_pyr_features(char *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols )
{
  uint32_t offset, ww, hh, x, y,center,j ;
  uint8_t pattern;

  offset=0;
/*  ww=win_W;*/
/*  hh=win_H;*/
  ww=img_nCols;
  hh=img_nRows;
  while(1)
  {
    for(x=1; x < ww-1; x++)
    {
      for(y=1; y< hh-1; y++)
      {
        pattern = 0;
        center = img[LIBLBP_INDEX(y,x,img_nRows)];
        if(img[LIBLBP_INDEX(y-1,x-1,img_nRows)] < center) pattern = pattern | 0x01;
        if(img[LIBLBP_INDEX(y-1,x,img_nRows)] < center)   pattern = pattern | 0x02;
        if(img[LIBLBP_INDEX(y-1,x+1,img_nRows)] < center) pattern = pattern | 0x04;
        if(img[LIBLBP_INDEX(y,x-1,img_nRows)] < center)   pattern = pattern | 0x08;
        if(img[LIBLBP_INDEX(y,x+1,img_nRows)] < center)   pattern = pattern | 0x10;
        if(img[LIBLBP_INDEX(y+1,x-1,img_nRows)] < center) pattern = pattern | 0x20;
        if(img[LIBLBP_INDEX(y+1,x,img_nRows)] < center)   pattern = pattern | 0x40;
        if(img[LIBLBP_INDEX(y+1,x+1,img_nRows)] < center) pattern = pattern | 0x80;

        vec[offset+pattern]++;
        offset += 256; 

      }
    }
    if(vec_nDim <= offset) 
      return;

    if(ww % 2 == 1) ww--;
    if(hh % 2 == 1) hh--;

    ww = ww/2;
    for(x=0; x < ww; x++)
      for(j=0; j < hh; j++)
        img[LIBLBP_INDEX(j,x,img_nRows)] = img[LIBLBP_INDEX(j,2*x,img_nRows)] + 
          img[LIBLBP_INDEX(j,2*x+1,img_nRows)];

    hh = hh/2;
    for(y=0; y < hh; y++)
      for(j=0; j < ww; j++)
        img[LIBLBP_INDEX(y,j,img_nRows)] = img[LIBLBP_INDEX(2*y,j,img_nRows)] + 
          img[LIBLBP_INDEX(2*y+1,j,img_nRows)];
    
  }

  return;
}


/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------*/
double liblbp_pyr_dotprod(double *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols)
{
  double dot_prod = 0;
  uint32_t offset=0;
  uint32_t ww, hh, center, x, y, j;
  uint8_t pattern;
  
/*  ww=win_W;*/
/*  hh=win_H;*/
  ww=img_nCols;
  hh=img_nRows;
  while(1)
  {
    for(x=1; x < ww-1; x++)
    {
      for(y=1; y< hh-1; y++)
      {
        pattern = 0;
        center = img[LIBLBP_INDEX(y,x,img_nRows)];
        if(img[LIBLBP_INDEX(y-1,x-1,img_nRows)] < center) pattern = pattern | 0x01;
        if(img[LIBLBP_INDEX(y-1,x,img_nRows)] < center)   pattern = pattern | 0x02;
        if(img[LIBLBP_INDEX(y-1,x+1,img_nRows)] < center) pattern = pattern | 0x04;
        if(img[LIBLBP_INDEX(y,x-1,img_nRows)] < center)   pattern = pattern | 0x08;
        if(img[LIBLBP_INDEX(y,x+1,img_nRows)] < center)   pattern = pattern | 0x10;
        if(img[LIBLBP_INDEX(y+1,x-1,img_nRows)] < center) pattern = pattern | 0x20;
        if(img[LIBLBP_INDEX(y+1,x,img_nRows)] < center)   pattern = pattern | 0x40;
        if(img[LIBLBP_INDEX(y+1,x+1,img_nRows)] < center) pattern = pattern | 0x80;

        dot_prod += vec[offset+pattern];
        offset += 256; 


      }
    }
    if(vec_nDim <= offset) 
      return(dot_prod);


    if(ww % 2 == 1) ww--;
    if(hh % 2 == 1) hh--;

    ww = ww/2;
    for(x=0; x < ww; x++)
      for(j=0; j < hh; j++)
        img[LIBLBP_INDEX(j,x,img_nRows)] = img[LIBLBP_INDEX(j,2*x,img_nRows)] + 
                                          img[LIBLBP_INDEX(j,2*x+1,img_nRows)];

    hh = hh/2;
    for(y=0; y < hh; y++)
      for(j=0; j < ww; j++)
        img[LIBLBP_INDEX(y,j,img_nRows)] = img[LIBLBP_INDEX(2*y,j,img_nRows)] + 
                                           img[LIBLBP_INDEX(2*y+1,j,img_nRows)];    
  }
 
  
}


/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------*/
void liblbp_pyr_addvec(int64_t *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols)
{
  uint32_t offset, ww, hh, x, y, center,j ;
  uint8_t pattern;

  offset=0;
/*  ww=win_W;*/
/*  hh=win_H;*/
  ww=img_nCols;
  hh=img_nRows;
  while(1)
  {
    for(x=1; x < ww-1; x++)
    {
      for(y=1; y< hh-1; y++)
      {
        pattern = 0;
        center = img[LIBLBP_INDEX(y,x,img_nRows)];
        if(img[LIBLBP_INDEX(y-1,x-1,img_nRows)] < center) pattern = pattern | 0x01;
        if(img[LIBLBP_INDEX(y-1,x,img_nRows)] < center)   pattern = pattern | 0x02;
        if(img[LIBLBP_INDEX(y-1,x+1,img_nRows)] < center) pattern = pattern | 0x04;
        if(img[LIBLBP_INDEX(y,x-1,img_nRows)] < center)   pattern = pattern | 0x08;
        if(img[LIBLBP_INDEX(y,x+1,img_nRows)] < center)   pattern = pattern | 0x10;
        if(img[LIBLBP_INDEX(y+1,x-1,img_nRows)] < center) pattern = pattern | 0x20;
        if(img[LIBLBP_INDEX(y+1,x,img_nRows)] < center)   pattern = pattern | 0x40;
        if(img[LIBLBP_INDEX(y+1,x+1,img_nRows)] < center) pattern = pattern | 0x80;

        vec[offset+pattern]++;
        offset += 256; 

      }
    }
    if(vec_nDim <= offset) 
      return;

    if(ww % 2 == 1) ww--;
    if(hh % 2 == 1) hh--;

    ww = ww/2;
    for(x=0; x < ww; x++)
      for(j=0; j < hh; j++)
        img[LIBLBP_INDEX(j,x,img_nRows)] = img[LIBLBP_INDEX(j,2*x,img_nRows)] + 
             img[LIBLBP_INDEX(j,2*x+1,img_nRows)];

    hh = hh/2;
    for(y=0; y < hh; y++)
      for(j=0; j < ww; j++)
        img[LIBLBP_INDEX(y,j,img_nRows)] = img[LIBLBP_INDEX(2*y,j,img_nRows)] + 
          img[LIBLBP_INDEX(2*y+1,j,img_nRows)];
    
  }

  return;
}



/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------*/
void liblbp_pyr_subvec(int64_t *vec, uint32_t vec_nDim, uint32_t *img, uint16_t img_nRows, uint16_t img_nCols)
{
  uint32_t offset, ww, hh, x, y,center,j ;
  uint8_t pattern;

  offset=0;
/*  ww=win_W;*/
/*  hh=win_H;*/
  ww=img_nCols;
  hh=img_nRows;
  while(1)
  {
    for(x=1; x < ww-1; x++)
    {
      for(y=1; y< hh-1; y++)
      {
        pattern = 0;
        center = img[LIBLBP_INDEX(y,x,img_nRows)];
        if(img[LIBLBP_INDEX(y-1,x-1,img_nRows)] < center) pattern = pattern | 0x01;
        if(img[LIBLBP_INDEX(y-1,x,img_nRows)] < center)   pattern = pattern | 0x02;
        if(img[LIBLBP_INDEX(y-1,x+1,img_nRows)] < center) pattern = pattern | 0x04;
        if(img[LIBLBP_INDEX(y,x-1,img_nRows)] < center)   pattern = pattern | 0x08;
        if(img[LIBLBP_INDEX(y,x+1,img_nRows)] < center)   pattern = pattern | 0x10;
        if(img[LIBLBP_INDEX(y+1,x-1,img_nRows)] < center) pattern = pattern | 0x20;
        if(img[LIBLBP_INDEX(y+1,x,img_nRows)] < center)   pattern = pattern | 0x40;
        if(img[LIBLBP_INDEX(y+1,x+1,img_nRows)] < center) pattern = pattern | 0x80;

        vec[offset+pattern]--;
        offset += 256; 

      }
    }
    if(vec_nDim <= offset) 
      return;

    if(ww % 2 == 1) ww--;
    if(hh % 2 == 1) hh--;

    ww = ww/2;
    for(x=0; x < ww; x++)
      for(j=0; j < hh; j++)
        img[LIBLBP_INDEX(j,x,img_nRows)] = img[LIBLBP_INDEX(j,2*x,img_nRows)] + 
          img[LIBLBP_INDEX(j,2*x+1,img_nRows)];

    hh = hh/2;
    for(y=0; y < hh; y++)
      for(j=0; j < ww; j++)
        img[LIBLBP_INDEX(y,j,img_nRows)] = img[LIBLBP_INDEX(2*y,j,img_nRows)] + 
          img[LIBLBP_INDEX(2*y+1,j,img_nRows)];
    
  }

  return;
}


/*-----------------------------------------------------------------------
  -----------------------------------------------------------------------*/

uint32_t liblbp_pyr_get_dim(uint16_t img_nRows, uint16_t img_nCols, uint16_t nPyramids)
{
  uint32_t w, h, N, i;

  for(w=img_nCols, h=img_nRows, N=0, i=0; i < nPyramids && LIBLBP_MIN(w,h) >= 3; i++)
  {
    N += (w-2)*(h-2);

    if(w % 2) w--;
    if(h % 2) h--;
    w = w/2;
    h = h/2;
  }
  return(256*N);
}

