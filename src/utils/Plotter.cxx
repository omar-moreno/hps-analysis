
#include <Plotter.h>

Plotter::Plotter()
    : type("float"), 
      color(kAzure+2) {
}

Plotter::~Plotter() {
}

TH1* Plotter::build1DHistogram(std::string name, int n_bins, int x_min, int x_max) { 
    
    if (histogram1D_map[name] != NULL) { 
        throw std::runtime_error("A histogram with that name already exist!");
    }
    
    std::string root_name = name + "_" + std::to_string(rand()%10000);
    
    if (type.compare("float") == 0) {
        histogram1D_map[name] = new TH1F(root_name.c_str(), name.c_str(), n_bins, x_min, x_max);
    } else if (type.compare("double") == 0) { 
        histogram1D_map[name] = new TH1D(root_name.c_str(), name.c_str(), n_bins, x_min, x_max);
    }

    histogram1D_map[name]->SetLineColor(color);
}


TH2* Plotter::build2DHistogram(std::string name, int n_bins_x, int x_min, int x_max, 
        int n_bins_y, int y_min, int y_max) { 
    
    if (histogram2D_map[name] != NULL) { 
        throw std::runtime_error("A histogram with that name already exist!");
    }

    std::string root_name = name + "_" + std::to_string(rand()%10000);

    if (type.compare("float") == 0) {
        histogram2D_map[name] = new TH2F(root_name.c_str(), name.c_str(), n_bins_x, x_min, x_max,
                n_bins_y, y_min, y_max);
    } else if (type.compare("double") == 0) { 
        histogram2D_map[name] = new TH2D(root_name.c_str(), name.c_str(), n_bins_x, x_min, x_max,
                n_bins_y, y_min, y_max);
    }
    
    histogram2D_map[name]->SetLineColor(color);
}

TH1* Plotter::get1DHistogram(std::string name) { 
    
    if (histogram1D_map[name] == NULL) { 
        throw std::runtime_error("Histogram doesn't exist.");
    }

    return histogram1D_map[name];
}

TH2* Plotter::get2DHistogram(std::string name) { 

    if (histogram2D_map[name] == NULL) { 
        throw std::runtime_error("Histogram doesn't exist.");
    }

    return histogram2D_map[name]; 
}
