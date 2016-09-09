/*
 * See Copyright Notice in atomorph.h
 */ 

#include "atomorph.h"
 
namespace am {

class thread {
    public:
    thread();
    ~thread();
    
    void         next_state()                      {skip_state = true;}
    unsigned     get_state()                       {return state;}        
    size_t       get_identifier()                  {return identifier;}
    void         set_identifier(size_t id)         {if (paused) identifier = id;}    
    frame       *get_frame(size_t nr)              {return (has_frame(nr) ? &(frames[nr]) : nullptr);}
    const std::map<size_t, chain> *get_chains()    {return (const std::map<size_t, chain> *) &chains;}
    unsigned     get_seed()                        {return seed;}
    double       get_energy();
    void         set_frame(size_t nr, frame* f);
    size_t       get_blob_count(size_t frame);
    size_t       get_blob_count();
    void         set_bbox(uint16_t x1, uint16_t y1, uint16_t x2, uint16_t y2);
    void         set_blob_weights(unsigned char rgba, unsigned char size, unsigned char xy);
        
    bool clear();        
    void set_seed(unsigned seed);
    void set_blob_delimiter(unsigned char d) {if (paused) blob_delimiter = d;} // Color distance formula.
    void set_blob_threshold(double        t) {if (paused) blob_threshold = t;} // Color difference toleration.
    void set_blob_max_size (size_t        s) {if (paused) blob_max_size  = s;} // Maximum size of a blob.
    void set_blob_min_size (size_t        s) {if (paused) blob_min_size  = s;} // Minimum size of a blob.
    void set_blob_box_grip (uint16_t      g) {if (paused) blob_box_grip  = g;} // Defines the bounding box.
    void set_blob_box_samples(size_t      s) {if (paused) blob_box_samples=s;} // Number of dust samples.
    void set_blob_number   (size_t        n) {if (paused) blob_number    = n;} // Preferred number of blobs.
    void set_degeneration  (size_t        d) {if (paused) degenerate     = d;} // Degeneration period.
    void set_density       (uint16_t      d) {if (paused) point_density  = d;} // Key points per pixel.
    void set_threads       (size_t        n) {if (paused) thread_count   = n;} // Number of threads to spawn.
    void set_cycle_length  (size_t        l) {if (paused) cycle_length   = l;} // Steps per iteration to repeat.
    
    bool is_running() const {return running;}
    bool is_paused()  const {return paused;}

    void stop()   { signal_stop  = true; while(running) std::this_thread::sleep_for(std::chrono::milliseconds(0)); step_thread.join();}
    void start()  { running      = true; step_thread = std::thread(&thread::run, this);}
    void pause()  { signal_pause = true; }
    void resume() { paused       = false;}       
    
    void start (size_t iterations);    
    void resume(size_t iterations);
    void start (double seconds);
    void resume(double seconds);     
    
    private:
    double blob_distance(const blob *b1, const blob *b2);
    pixel get_pixel(size_t frame, size_t position);

    bool          blobify       ();
    bool          unify         ();        
    bool          match         ();
    bool          morph         ();
    bool          blobify_frame (size_t f);
    bool          unify_frame   (size_t f);
    bool          init_morph    ();
    double        get_energy    (struct blob ***map);
    double        get_energy    (chain *ch);
    void          fix_volatiles (std::vector<blob *> *to_be_fixed);
    inline bool   has_frame     (size_t f              ) const {return (frames.find(f) != frames.end());}
    inline bool   has_pixel     (size_t f,   size_t pos)       {return (has_frame(f) && frames[f].pixels.find(pos) != frames[f].pixels.end());}
    inline bool   has_blob      (size_t f,   size_t b  )       {return (has_frame(f) && frames[f].blobs.size() > b);}
    inline bool   can_expand    (size_t frame, size_t from_pos, size_t to_pos);    
    
    void run();
    void step();    
    
    std::map<size_t, frame> frames;
    std::map<size_t, chain> chains; // Key is the blob's group.
    
    uint16_t point_density = 1; // How many key points to create per pixel at minimum.
    
    struct blob *** blob_map;
    size_t blob_map_w; // Number of blobs per key frame.
    size_t blob_map_h; // Number of key frames.
    double blob_map_e; // Energy of the current blob map.
    double best_blob_map_e;
    
    double chain_map_e;
    
    // Minimum Bounding Box that fits all pixels of all frames:
    uint16_t bbox_x1 = UINT16_MAX;
    uint16_t bbox_y1 = UINT16_MAX;
    uint16_t bbox_x2 = 0;
    uint16_t bbox_y2 = 0;
    uint32_t bbox_d  = 0; // Squared diagonal length.
    
    unsigned char blob_delimiter   = HSP;
    double        blob_threshold   = 1.0;
    size_t        blob_max_size    = SIZE_MAX;
    size_t        blob_min_size    = 1;
    uint16_t      blob_box_grip    = UINT16_MAX;        
    size_t        blob_box_samples = 10;    
    size_t        blob_number      = 1;   
    size_t        thread_count     = 0;
    size_t        cycle_length     = 1;
    
    double        blob_rgba_weight = 0.33;
    double        blob_size_weight = 0.33;
    double        blob_xy_weight   = 0.34;
    
    size_t   identifier = SIZE_MAX;
    unsigned state      = STATE_BLOB_DETECTION;
    unsigned seed       = 0;
    bool     skip_state = false;
    size_t   degenerate = 0;
    size_t   counter    = 0;
    double   best_e;
    bool     deviant; // True when a bad solution has been accepted lately.
    
    std::default_random_engine e1;
        
    std::atomic<size_t> iterations;
    std::atomic<double> seconds;
        
    std::atomic<bool> signal_stop;
    std::atomic<bool> signal_pause;
    std::atomic<bool> running;
    std::atomic<bool> paused;
    std::thread step_thread;
};

}


