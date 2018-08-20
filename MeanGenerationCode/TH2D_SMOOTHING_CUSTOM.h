// @(#)root/hist:$Id$
 // Author: Rene Brun   26/12/94
 
 /*************************************************************************
  * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
  * All rights reserved.                                                  *
  *                                                                       *
  * For the licensing terms see $ROOTSYS/LICENSE.                         *
  * For the list of contributors see $ROOTSYS/README/CREDITS.             *
  *************************************************************************/

 ////////////////////////////////////////////////////////////////////////////////
 /// Smooth bin contents of this 2-d histogram using kernel algorithms
 /// similar to the ones used in the raster graphics community.
 /// Bin contents in the active range are replaced by their smooth values.
 /// If Errors are defined via Sumw2, they are also scaled and computed.
 /// However, note the resulting errors will be correlated between different-bins, so
 /// the errors should not be used blindly to perform any calculation involving several bins,
 /// like fitting the histogram.  One would need to compute also the bin by bin correlation matrix.
 ///
 /// 3 kernels are proposed k5a, k5b and k3a.
 /// k5a and k5b act on 5x5 cells (i-2,i-1,i,i+1,i+2, and same for j)
 /// k5b is a bit more stronger in smoothing
 /// k3a acts only on 3x3 cells (i-1,i,i+1, and same for j).
 /// By default the kernel "k5a" is used. You can select the kernels "k5b" or "k3a"
 /// via the option argument.
 /// If TAxis::SetRange has been called on the x or/and y axis, only the bins
 /// in the specified range are smoothed.
 /// In the current implementation if the first argument is not used (default value=1).
 ///
 /// implementation by David McKee (dmckee@bama.ua.edu). Extended by Rene Brun
 


///CODE HAS BEEN UPDATED BY KHALID GAMEIL TO ADD FUNCTIONALITY
// Root
#ifndef TH2D_SMOOTHING_CUSTOM
#define TH2D_SMOOTHING_CUSTOM

#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TGFrame.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TPaveLabel.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TString.h"
#include "TStyle.h"
#include "TVirtualPad.h"

namespace TH2D_Custom{
 void Smooth(TH2D * h, Option_t *option)
 {
    Double_t k5a[5][5] =  { { 0, 0, 1, 0, 0 },
                            { 0, 2, 2, 2, 0 },
                            { 1, 2, 5, 2, 1 },
                            { 0, 2, 2, 2, 0 },
                            { 0, 0, 1, 0, 0 } };
    Double_t k5b[5][5] =  { { 0, 1, 2, 1, 0 },
                            { 1, 2, 4, 2, 1 },
                            { 2, 4, 8, 4, 2 },
                            { 1, 2, 4, 2, 1 },
                            { 0, 1, 2, 1, 0 } };
    Double_t k3a[3][3] =  { { 0, 0.1, 0 },
                            { 0.1, 1, 0.1 },
                            { 0, 0.1, 0 } };
 
    TString opt = option;
    opt.ToLower();
    Int_t ksize_x=5;
    Int_t ksize_y=5;
    Double_t *kernel = &k5a[0][0];
    if (opt.Contains("k5b")) kernel = &k5b[0][0];
    if (opt.Contains("k3a")) {
       kernel = &k3a[0][0];
       ksize_x=3;
       ksize_y=3;
    }
 
    // find i,j ranges
    Int_t ifirst = h->GetXaxis()->GetFirst();
    Int_t ilast  = h->GetXaxis()->GetLast();
    Int_t jfirst = h->GetYaxis()->GetFirst();
    Int_t jlast  = h->GetYaxis()->GetLast();
 
    // Determine the size of the bin buffer(s) needed
    Double_t nentries = h->GetEntries();
    Int_t nx = h->GetNbinsX();
    Int_t ny = h->GetNbinsY();
    Int_t bufSize  = (nx+2)*(ny+2);
    Double_t *buf  = new Double_t[bufSize];
    Double_t *ebuf = 0;
 
    // Copy all the data to the temporary buffers
    Int_t i,j,bin;
    for (i=ifirst; i<=ilast; i++){
       for (j=jfirst; j<=jlast; j++){
          bin = h->GetBin(i,j);
          buf[bin] = h->GetBinContent(bin);
          if (ebuf) ebuf[bin]=h->GetBinError(bin);
       }
    }
 
    // Kernel tail sizes (kernel sizes must be odd for this to work!)
    Int_t x_push = (ksize_x-1)/2;
    Int_t y_push = (ksize_y-1)/2;
 
    // main work loop
    for (i=ifirst; i<=ilast; i++){
       for (j=jfirst; j<=jlast; j++) {
          Double_t content = 0.0;
          Double_t error = 0.0;
          Double_t norm = 1;
 
          for (Int_t n=0; n<ksize_x; n++) {
             for (Int_t m=0; m<ksize_y; m++) {
                Int_t xb = i+(n-x_push);
                Int_t yb = j+(m-y_push);
                if ( (xb >= 1) && (xb <= nx) && (yb >= 1) && (yb <= ny) ) {
                   bin = h->GetBin(xb,yb);
                   Double_t k = kernel[n*ksize_y +m];
                   //if ( (k != 0.0 ) && (buf[bin] != 0.0) ) { // General version probably does not want the second condition
                   if ( k != 0.0 ) {
                      //norm    += k;
                      content += k*buf[bin];
                      //if (ebuf) error   += k*k*ebuf[bin]*ebuf[bin];
                   }
                }
             }
          }
 
          if ( norm != 0.0 ) {
             h->SetBinContent(i,j,content/norm);
             if (ebuf) {
                //error /= (norm*norm);
                //h->SetBinError(i,j,sqrt(error));
             }
          }
       }
    }
    h->SetEntries(nentries);
 
    delete [] buf;
    delete [] ebuf;
 }}
#endif 
