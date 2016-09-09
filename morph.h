/*
 * See Copyright Notice in atomorph.h
 */ 

#include "atomorph.h"
 
namespace am {

class morph {
    public:
    morph();
    ~morph();
    
    void        clear              ();
    bool        add_pixel          (size_t frame, pixel px);
    bool        add_frame          (size_t frame);
    size_t      get_pixel_count    (size_t frame);
    size_t      get_frame_key      (double t);
    size_t      get_blob_count     (size_t frame);
    size_t      get_blob_count     ();
    pixel       get_average_pixel  (size_t frame);
    pixel       get_average_pixel  (size_t frame, size_t blob);
    pixel       get_pixel          (size_t frame, size_t position);
    color       get_background     (uint16_t x, uint16_t y, double t);
    point       interpolate        (point pt1, point pt2, double pt1_weight);    
    pixel       interpolate        (pixel px1, pixel px2, double px1_weight);
    color       interpolate        (color c1, color c2, double c1_weight);
    color       interpolate        (color c1, color c2, double lag, double slope, double c1_weight);    
    const blob* get_pixels         (size_t blob, double t, std::vector<pixel> *to);
    void        get_pixels         (double t, std::vector<pixel> *to);
    void        set_seed           (unsigned seed);
    double      get_time           (size_t current_frame, size_t total_frames);
    double      normalize_time     (double t);

    uint16_t    get_width       ()                   {return width;}
    uint16_t    get_height      ()                   {return height;}
    size_t      get_frame_count ()                   {return frames.size();}
    unsigned    get_state       ()                   {return state;}
    void        next_state      ()                   {skip_state = true;}                       
    const blob* get_blob        (size_t f, size_t b) {return (has_blob(f, b) ? ((const blob*) frames[f].blobs[b]) : nullptr);}
    double      get_energy      ()                   {return energy;}
    pixel       blob2pixel      (const blob *bl);
    void        compute         ();                  // If not busy yet, sets worker thread to run indefinitely.
    void        iterate         (size_t iterations); // If not busy yet, sets worker thread to run the defined number of iterations.
    void        compute         (double seconds);    // If not busy yet, sets worker thread to run for the defined number of seconds.
    void        suspend         ();                  // Blocks the calling thread until morph is no longer busy. Returns when worker is paused.
    bool        suspend         (double timeout);    // Blocks the calling thread for timeout seconds at maximum. Returns false on timeout.
    bool        is_busy         () const {return (worker.is_running() && !worker.is_paused());}
    
    bool        synchronize     (); // Synchronizes data with the worker thread. Returns false when worker is busy at the moment.

    void        set_blob_delimiter  (unsigned char d) {blob_delimiter = d;} // Color distance formula.
    void        set_blob_threshold  (double        t) {blob_threshold = t;} // Color difference toleration.
    void        set_blob_max_size   (size_t        s) {blob_max_size  = s;} // Maximum size of a blob.
    void        set_blob_min_size   (size_t        s) {blob_min_size  = s;} // Minimum size of a blob.
    void        set_blob_box_grip   (uint16_t      g) {blob_box_grip  = g;} // Defines the bounding box.
    void        set_blob_box_samples(size_t        s) {blob_box_samples=s;} // Number of dust samples.
    void        set_blob_number     (size_t        n) {blob_number    = n;} // Preferred number of blobs.
    void        set_blob_rgba_weight(unsigned char w) {blob_rgba_weight=w;} // Blob RGBA weight.
    void        set_blob_size_weight(unsigned char w) {blob_size_weight=w;} // Blob size weight.
    void        set_blob_xy_weight  (unsigned char w) {blob_xy_weight  =w;} // Blob location weight.        
    void        set_degeneration    (size_t        d) {degeneration    =d;} // Solution degeneration period.
    void        set_motion          (unsigned char m) {motion          =m;} // How to interpolate movement.
    void        set_fading          (unsigned char f) {fading          =f;} // How to interpolate colors.
    void        set_threads         (size_t        t) {threads         =t;} // How many threads the worker thread can spawn.
    void        set_cycle_length    (size_t        c) {cycle_length    =c;} // How many times to morph per iteration.
    void        set_feather         (size_t        f) {feather         =f;} // Blob outer layers to fade out.
    void        set_keep_background (bool          k) {keep_background =k;} // Background is kept and cross-dissolved.
    void        set_finite          (bool          f) {finite          =f;} // Morph will not repeat itself.
    void        set_show_blobs      (unsigned      b) {show_blobs      =b;} // Determines how blobs are drawn.
    
    // These change identifier, causing the morph to restart:
    void        set_fluid           (unsigned      f);                      // Steps of fluid simulation per frame.    
    void        set_density         (uint16_t      d);                      // Key points per pixel.

    void        set_resolution      (uint16_t w, uint16_t h);

    private:
    inline bool   has_pixel     (size_t f,   size_t pos)       {return (has_frame(f) && frames[f].pixels.find(pos) != frames[f].pixels.end());}    
    inline bool   has_frame     (size_t f              ) const {return (frames.find(f) != frames.end());}    
    inline bool   has_blob      (size_t f,   size_t b  )       {return (has_frame(f) && frames[f].blobs.size() > b);}            
    void          refresh_frames();            

    blob * find_blob_by_group(size_t frame_key, size_t group);

    void          draw_atoms(double t, std::vector<pixel> *image);
    void          draw_fluid(double t, std::vector<pixel> *image);
    void          step_fluid(size_t frame_key, double t, double time);
    
    void          update_particle (Particle *p, const blob * blob_before, const blob * blob_after, point **points, 
                                   std::map<size_t, size_t>& src, std::map<size_t, size_t>& dest, double t, size_t y, 
                                   size_t y_next, size_t frame_key, size_t next_frame_key, size_t chain_key, double time);
            
    std::map<size_t, frame> frames; 
    std::map<size_t, chain> chains; // Key is the blob's group.
    
    unsigned char fading           =       PERLIN; // Color interpolation.
    unsigned char motion           =       SPLINE; // Trajectory interpolation.
    unsigned char blob_delimiter   =          HSP;
    double        blob_threshold   =          1.0;
    size_t        blob_number      =            1;
    size_t        blob_max_size    =     SIZE_MAX;
    size_t        blob_min_size    =            1;
    size_t        blob_box_samples =           10;
    uint16_t      blob_box_grip    =   UINT16_MAX;
    unsigned char blob_rgba_weight =            1;
    unsigned char blob_size_weight =            1;
    unsigned char blob_xy_weight   =            1;
    size_t        degeneration     =            0;
    uint16_t      density          =            1;
    size_t        threads          =            0;
    size_t        cycle_length     =         1000;
    size_t        feather          =            0;
    bool          keep_background  =        false;
    bool          blend_blobs      =        false;
    bool          finite           =        false;
    unsigned      fluidsteps       =            0;
    unsigned      show_blobs       =      TEXTURE;
    unsigned      seed;        
    
    // Minimum Bounding Box that fits all pixels of all frames:
    uint16_t bbox_x1 = UINT16_MAX;
    uint16_t bbox_y1 = UINT16_MAX;
    uint16_t bbox_x2 = 0;
    uint16_t bbox_y2 = 0;
    
    uint16_t width   = 0;
    uint16_t height  = 0;
    
    size_t        identifier=0; // Changed when anything changes in the key frames.
    unsigned      state     =STATE_BLOB_DETECTION;
    bool          skip_state=false;
    
    thread worker;    
    double energy=0.0;
    
    PerlinNoise lag_map;
    PerlinNoise slope_map;    
    
    FluidModel *fluid = nullptr;
    unsigned sparseness = 1;
    
    std::default_random_engine e1;
};

}

