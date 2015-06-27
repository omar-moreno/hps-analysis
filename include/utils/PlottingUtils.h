/**
 *	@file	PlottingUtils.h
 *	@brief	A compilation of plotting utilities.
 *	@author Omar Moreno <omoreno1@ucsc.edu>
 *			Santa Cruz Institute for Paritcle Physics
 *			University of California, Santa Cruz
 *	@date	2013
 *
 */

#ifndef __PLOTTING_UTILS_H__
#define __PLOTTING_UTILS_H__

//------------------//
//--- C++ StdLib ---//
//------------------//
#include <string>
#include <sstream>

//------------//
//--- ROOT ---//
//------------//
#include <TColor.h>
#include <TStyle.h>

namespace PlottingUtils { 


	/**
	 *	Setup the color palette used by ROOT
	 */
    void setPalette(); 
   
	/**
	 *  Setup the style used by ROOT
	 */
	void setStyle(); 

	/**
	 *	Convert a numerical value to a string.  This method should only be 
	 *	used when the C++11 to_string function is not available.
	 *
	 *	@parameter value Numerical value.
	 *	@return A string object representation of value
	 */
    template <typename T>
    std::string toString(const T value) {
            std::stringstream ss; 
            ss << value; 
            return ss.str(); 
    }

    /**
     *  Instantiate a given number of pointers to 1D histograms of type T and
     *  return them as a vector.  The pointers must be deleted by the user.
     *  Note that the number of bins is calculated from the histogram range.
     *
     *  @param histogram_n - Number of histograms to instantiate
     *  @param histogram_name - The name of the histogram
     *  @param x_min - Minimum histogram value along the x-axis
     *  @param x_min - Maximum histogram value along the x-axis
     *  @return A vector of 1D histograms of type T
     */
    template <typename T>
    std::vector <T*> book1DHistograms(int histogram_n, std::string histogram_name,
                                      int x_min, int x_max) {
        int bins_x = x_max -x_min;
        return book1DHistograms <T> (histogram_n, histogram_name,
                                     bins_x, x_min, x_max);
    }


    /**
     *  Instantiate a given number of pointers to 1D histograms of type T and
     *  return them as a vector.  The pointers must be deleted by the user.
     *
     *  @param histogram_n - Number of histograms to instantiate
     *  @param histogram_name - The name of the histogram
     *  @param bins_x - Number of bins along the x-axis
     *  @param x_min - Minimum histogram value along the x-axis
     *  @param x_min - Maximum histogram value along the x-axis
     *  @return A vector of 1D histograms of type T
     */
    template <typename T>
    std::vector <T*> book1DHistograms(int histogram_n, std::string histogram_name,
                                      int bins_x, int x_min, int x_max) {
        std::vector <T*> histograms;
        for (int histogram = 0; histogram < histogram_n; ++histogram) {
            std::string name = histogram_name + ": " + toString <int> (histogram);
            histograms.push_back(new T(name.c_str(), name.c_str(),
                                       bins_x, x_min, x_max));
        }
        return histograms;
    }


    /**
     *  Instantiate a given number of pointers to 2D histograms of type T and
     *  return them as a vector.  The pointers must be deleted by the user.
     *  Note that the number of bins is calculated from the histogram range.
     *
     *  @param histogram_n - Number of histograms to instantiate
     *  @param histogram_name - The name of the histogram
     *  @param x_min - Minimum histogram value along the x-axis
     *  @param x_min - Maximum histogram value along the x-axis
     *  @param y_min - Minimum histogram value along the y-axis
     *  @param y_min - Maximum histogram value along the y-axis
     *  @return A vector of 2D histograms of type T
     */
    template <typename T>
    std::vector <T*> book2DHistograms(int histogram_n, std::string histogram_name,
    							     int x_min, int x_max,
    							     int y_min, int y_max) {
        int bins_x = x_max -x_min;
        int bins_y = y_max - y_min;
        return book2DHistograms <T> (histogram_n, histogram_name,
                              bins_x, x_min, x_max,
                              bins_y, y_min, y_max);
    }

    /**
     *  Instantiate a given number of pointers to 2D histograms of type T and
     *  return them as a vector.  The pointers must be deleted by the user.
     *
     *  @param histogram_n - Number of histograms to instantiate
     *  @param histogram_name - The name of the histogram
     *  @param bins_x - Number of bins along the x-axis
     *  @param x_min - Minimum histogram value along the x-axis
     *  @param x_min - Maximum histogram value along the x-axis
     *  @param bins_y - Number of bins along the y-axis
     *  @param y_min - Minimum histogram value along the y-axis
     *  @param y_min - Maximum histogram value along the y-axis
     *  @return A vector of 2D histograms of type T
     */
    template <typename T>
    std::vector <T*> book2DHistograms(int histogram_n, std::string histogram_name,
    							     int bins_x, int x_min, int x_max,
    							     int bins_y, int y_min, int y_max) {
        std::vector<T*> histograms;
        for (int histogram = 0; histogram < histogram_n; ++histogram) {
            std::string name = histogram_name + ": " + toString <int> (histogram);
            histograms.push_back(new T(name.c_str(), name.c_str(),
                                       bins_x, x_min, x_max,
                                       bins_y, y_min, y_max));
        }
        return histograms;
    }

	/**
	 *
	 */
    template <typename T>
        void adjust1DPlotRange(T* &plot, int bin_threshold)
        {

            int histo_min = plot->FindFirstBinAbove(bin_threshold, 1);
            int histo_max = plot->FindLastBinAbove(bin_threshold, 1);
            plot->GetXaxis()->SetRange(histo_min, histo_max);
        }

	/**
	 *
	 */
    template <typename T>
        void adjust2DPlotRange(T* &plot, int bin_threshold)
        {
            int histo_x_min = plot->FindFirstBinAbove(bin_threshold, 1);
            int histo_x_max = plot->FindLastBinAbove(bin_threshold, 1);
            int histo_y_min = plot->FindFirstBinAbove(bin_threshold, 2);
            int histo_y_max = plot->FindLastBinAbove(bin_threshold, 2);
            plot->GetXaxis()->SetRange(histo_x_min, histo_x_max);
            plot->GetYaxis()->SetRange(histo_y_min, histo_y_max);

        }

}; // PlottingUtils

#endif // __PLOTTING_UTILS_H__
