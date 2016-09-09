/*
 * See Copyright Notice in atomorph.h
 */

#include <cstring>
#include "atomorph.h"

namespace am {

thread::thread() {
    signal_stop  = false;
    signal_pause = false;    
    running      = false;
    paused       = true;
    
    iterations   = 0;
    seconds      = 0.0;    
    
    blob_map     = nullptr;
    blob_map_w   = 0;
    blob_map_h   = 0;
    blob_map_e   = 0.0;
    chain_map_e  = 0.0;
    
    clear();
}

thread::~thread() {
    if (step_thread.joinable()) {
        pause();
        stop();
    }
    clear();
}

bool thread::clear() {
    if (running && !paused) return false;    
    
    if (blob_map) {
        for ( size_t i = 0; i < blob_map_w; ++i) {
            delete [] blob_map[i];
        }
        delete [] blob_map;
        blob_map   = nullptr;        
    }
    blob_map_w = 0;
    blob_map_h = 0;
    blob_map_e = 0.0;    
    best_blob_map_e = 0.0;    
    
    best_e      = std::numeric_limits<double>::infinity();
    chain_map_e = std::numeric_limits<double>::infinity();    
    deviant= false;
    
    while (!frames.empty()) {
        frame f = frames.begin()->second;
        while (!f.blobs.empty()) {
            delete f.blobs.back();
            f.blobs.pop_back();
        }
        frames.erase(frames.begin());
    }
    frames.clear();    
    
    while (!chains.empty()) {
        chain c = chains.begin()->second;
        clear_chain(&c);
        chains.erase(chains.begin());
    }
    chains.clear();        
    
    state      = STATE_BLOB_DETECTION;
    identifier = SIZE_MAX;    
    counter    = 1;
    
    bbox_x1 = UINT16_MAX;
    bbox_y1 = UINT16_MAX;
    bbox_x2 = 0;
    bbox_y2 = 0;
    bbox_d  = 0;
    
    return true;
}

void thread::set_seed(unsigned seed) {
    if (running && !paused) return;
    this->seed = seed;
    std::default_random_engine e(seed);
    e1 = e;
}

void thread::set_frame (size_t nr, frame * frame_pt) {
    if (running && !paused) return;
        
    frames[nr].pixels.insert(frame_pt->pixels.begin(), frame_pt->pixels.end());
    frames[nr].x = frame_pt->x;
    frames[nr].y = frame_pt->y;
    frames[nr].r = frame_pt->r;
    frames[nr].g = frame_pt->g;
    frames[nr].b = frame_pt->b;
    frames[nr].a = frame_pt->a;                    
    frames[nr].first_expansion=0; 
    frames[nr].first_dust     =0; 
    
    // Refresh frame indexes:
    std::map<size_t, frame>::iterator it;
    size_t index = 0;
    for (it=frames.begin(); it!=frames.end(); ++it) {
        it->second.index = index++;
    }                
}

void thread::resume(size_t iterations) {
    if (!is_paused()) return;
    this->iterations = iterations;
    resume();
}

void thread::resume(double seconds) {
    if (!is_paused()) return;
    this->seconds = seconds;
    resume();
}

void thread::start(size_t iterations) {
    if (is_running()) return;
    this->iterations = iterations;
    start();
    resume();
}

void thread::start(double seconds) {
    if (is_running()) return;
    this->seconds = seconds;
    start();
    resume();
}

void thread::step() {
    switch (state) {
        case STATE_BLOB_DETECTION:   if (blobify()) {
                                         state = STATE_BLOB_UNIFICATION;
                                         counter = 0;
                                     } 
                                     break;
        case STATE_BLOB_UNIFICATION: if (unify()) {
                                         state = STATE_BLOB_MATCHING;
                                         counter = 0;                                         
                                     }
                                     break;
        case STATE_BLOB_MATCHING:    if (skip_state || match()) {
                                         if (!init_morph()) {
                                             state = STATE_DONE;
                                         }
                                         else {
                                             state = STATE_ATOM_MORPHING;
                                         }
                                         counter = 0;
                                     }
                                     break;        
        case STATE_ATOM_MORPHING:    if (skip_state || morph()) {
                                         state = STATE_DONE;
                                         counter = 0;
                                     }
                                     break;
        default:                     break;
    }
    skip_state = false;
    counter++;
    return;
}

void thread::run() {
    std::chrono::steady_clock::time_point start,end;
    start = std::chrono::steady_clock::now();

    while (!signal_stop) {
        if (!signal_pause && !paused) {
            step();
            if (iterations > 0 && --iterations == 0) {
                // When the defined number of steps have been made, automatically pause the thread.
                signal_pause = true;
            }
            
            if (seconds > 0.0) {
                // When worked more than the defined number of seconds, automatically pause the thread.
                end = std::chrono::steady_clock::now();                 
                if (std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() > (seconds * 1000000000)) {
                    signal_pause = true;
                    seconds      = 0.0;
                }            
            }
        }
        else {
            paused       = true;
            signal_pause = false;
            std::this_thread::sleep_for(std::chrono::milliseconds(0));
            
            start = std::chrono::steady_clock::now();
        }                
    }
    running      = false;
    signal_stop  = false;
    signal_pause = false;
    paused       = true;    
}

bool thread::can_expand(size_t frame, size_t from_pos, size_t to_pos) {
    if (!has_pixel(frame, to_pos)) return false;    
    
    pixel p1 = get_pixel(frame, from_pos);
    pixel p2 = get_pixel(frame, to_pos);
    
    if (color_distance(p1.c, p2.c) <= blob_threshold) return true;
    
    return false;
}

pixel thread::get_pixel (size_t f, size_t pos) {
    if (!has_pixel(f,pos)) return create_pixel(0,0,0,0,0,0);
    return frames[f].pixels[pos];
}

// Returns false when not finished.
bool thread::blobify_frame(size_t frame_key) {
    frame* f = &frames[frame_key];
    size_t sz = f->blobs.size();
    size_t pos,to_pos,x,y,b;        
        
    if (sz == 0) {
        // First step is to interpret all pixels as atomic blobs.
        std::map<size_t, pixel>::iterator it;
        for (it=f->pixels.begin(); it!=f->pixels.end(); ++it) {
            pos = it->first;
            blob* new_blob = new blob;
            new_blob->surface.insert(pos); 
            new_blob->x = it->second.x;
            new_blob->y = it->second.y;
            new_blob->r = (it->second.c.r/255.0);
            new_blob->g = (it->second.c.g/255.0);
            new_blob->b = (it->second.c.b/255.0);
            new_blob->a = (it->second.c.a/255.0);
            new_blob->group = pos;
            x = pos % (UINT16_MAX+1);
            y = pos / (UINT16_MAX+1);

            // Add reasonable border.
            if (x < UINT16_MAX && has_pixel(frame_key, (to_pos = xy2pos(x+1, y  )))) new_blob->border.insert(to_pos);
            if (x >          0 && has_pixel(frame_key, (to_pos = xy2pos(x-1, y  )))) new_blob->border.insert(to_pos);
            if (y < UINT16_MAX && has_pixel(frame_key, (to_pos = xy2pos(x  , y+1)))) new_blob->border.insert(to_pos);
            if (y >          0 && has_pixel(frame_key, (to_pos = xy2pos(x  , y-1)))) new_blob->border.insert(to_pos);

            f->blobs.push_back(new_blob);
            f->owners[pos] = new_blob;
        }        
        
        std::shuffle(f->blobs.begin(), f->blobs.end(), e1);
        sz = f->blobs.size();

        for (size_t bb = 0; bb<sz; ++bb) {
            f->blobs[bb]->index = bb;
        }

        if (f->blobs.size() == 0 || blob_max_size == 1) return true;

        return false;
    }       
    
    Again: // For popping all sequential nullptrs from the end of the vector.
    sz = f->blobs.size();
    b  = 0;
    if (sz > 1) {
        // Always take the last blob from the initially shuffled vector
        // because popping the last element is fast.
        b = sz-1;
    }
    else return true;

    blob* bl = (f->blobs[b]);
    
    if (!bl || bl->surface.size() == 0) {
        // This blob is empty, remove it.
        f->blobs.pop_back();
        delete bl;
        goto Again;
    }    

    if (f->blobs.size() <= blob_number) return true;

    std::set<size_t>::iterator it;    
    std::vector<size_t> erased_border;
    std::set<blob *>    checked;
    bool   best_found = false;
    double best_dist  = 0.0;
    blob * best_match = nullptr;
    size_t best_equal = 0;
    
    if (bl->surface.size() >= blob_min_size) {
        // This blob has met its minimum size, prevent it from being unified later.
        bl->unified = true;
    }    
    
    bool check_threshold = true;    
    if (bl->surface.size() < blob_max_size) {
        IgnoreColourThreshold:
        
        for (it=bl->border.begin(); it!=bl->border.end(); ++it) {
            pos = *it;
            if (bl->surface.count(pos) != 0
            || f->owners.count(pos) == 0
            || checked.count(f->owners[pos]) != 0
            || !has_pixel(frame_key, pos)) {
                // Don't expand into itself and empty positions.
                if (check_threshold && bl->surface.size() >= blob_min_size) {
                    erased_border.push_back(pos);
                }
            }
            else {
                // Find out which blob currently owns this pixel.
                blob *owner = f->owners[pos];
                checked.insert(owner);
                if (owner != bl) {
                    double dist;
                    color c1 = create_color(owner->r,owner->g,owner->b,owner->a);
                    color c2 = create_color(bl->r, bl->g, bl->b, bl->a);
                    dist = color_distance(c1, c2);                                        
                    
                    if ((dist <= blob_threshold || !check_threshold) && (!best_found || dist <= best_dist)) {
                        if (best_found && dist == best_dist) {
                            std::uniform_int_distribution<size_t> uniform_dist_border (0, best_equal);    
                            best_equal++;
                            if (uniform_dist_border(e1) != 0) continue;                                                
                        }
                        else {
                            best_equal = 1;
                        } 
                        best_found = true;
                        best_match = owner;
                        best_dist  = dist;
                    }            
                }
            }
        }
    }    

    if (best_found) {
        // Merge with this blob.
        size_t surface_1 = bl->surface.size();
        size_t surface_2 = best_match->surface.size();
        double weight    = (surface_2) / double(surface_1 + surface_2);

        bl->x = ((1.0 - weight) * bl->x + weight * best_match->x);
        bl->y = ((1.0 - weight) * bl->y + weight * best_match->y);
        bl->r = ((1.0 - weight) * bl->r + weight * best_match->r);
        bl->g = ((1.0 - weight) * bl->g + weight * best_match->g);
        bl->b = ((1.0 - weight) * bl->b + weight * best_match->b);
        bl->a = ((1.0 - weight) * bl->a + weight * best_match->a);
        
        for (it=best_match->surface.begin(); it!=best_match->surface.end(); ++it) {
            f->owners[*it] = bl;
            bl->surface.insert(*it);
        }
        
        for (it=best_match->border.begin(); it!=best_match->border.end(); ++it) {
            bl->border.insert(*it);
        }                
 
        best_match->surface.clear();
        best_match->border .clear();                     

        size_t blob_index = best_match->index;
        f->blobs[blob_index] = bl;
        bl->index            = blob_index;
        bl->unified          = (bl->unified || best_match->unified);
        
        delete best_match;
        f->blobs.pop_back();     
    }
    else {        
        size_t bb = f->first_expansion;        
        
        if (bl->surface.size() < blob_min_size && bl->border.size() > 0 && check_threshold) {
            // This outlier has neighbours but it is undersized. Discard colour threshold constraint.
            check_threshold = false;
            checked.clear();
            goto IgnoreColourThreshold;
        }
        
        bl->border.clear(); // Don't attempt to expand any more.
        
        if (bb < sz) {
            // Swap it with last random blob.
            blob *other  = f->blobs[bb];            
            if (other) other->index = b;
            
            f->blobs[b]  = other;
            f->blobs[bb] = bl;
            bl->index    = bb;
            
            f->first_expansion++;
        }
        else return true; // None of the remaining blobs can expand any more.                     
    }

    // Remove erased borders.
    size_t esz = erased_border.size();
    for (size_t i=0; i<esz; ++i) {
        bl->border.erase(erased_border[i]);        
    }
    
    return false;
}

bool thread::blobify() {
    bool done = true;
    std::map<size_t, frame>::iterator it;
    
    for (it=frames.begin(); it!=frames.end(); ++it) {
        done &= blobify_frame(it->first);
    }
    
    return done;
}

// Returns false when not finished.
bool thread::unify_frame(size_t frame_key) {
    frame* f = &frames[frame_key];
    size_t fexp; // Only unify blobs that cannot expand any more.
    size_t sz;
    size_t bb;
    blob *bl, *ubl;
    blob *best;
    size_t score = 0;
    uint16_t ubl_x, ubl_y;
    size_t samples;
    bool done = true;
    size_t blobs_before = 0;
    std::set<size_t>::iterator it;

    // Remove all sequential nullptrs from the end of blobs.
    while (!f->blobs.empty() && f->blobs.back() == nullptr) {
        f->blobs.pop_back();
    }
    
    blobs_before = f->blobs.size();
    if (blob_box_samples == 0 || blobs_before <= blob_number) return true;
    
    Again:
    ubl     = nullptr;
    best    = nullptr;
    samples = 0;
    
    NotEnoughSamples:                
    sz      = f->blobs.size();        
    fexp    = f->first_expansion;        
    
    for (bb = f->first_dust; bb<fexp; ++bb) {        
        if ((bl = f->blobs[bb]) == nullptr) continue;

        if (bl == ubl) {
            // Cycle has ended.
            goto Samples;
        }

        if (bl->unified) {
            size_t fd = f->first_dust;
            if (bb != fd) {
                // Swap this already unified blob with current first dust
                // and increase the first dust index by 1.
                blob *fdb = f->blobs[fd];
                if (fdb) fdb->index = bb;
                bl->index = fd;               
                f->blobs[fd] = bl;
                f->blobs[bb] = fdb;
            }
            f->first_dust++;
            goto Again;            
        }
        
        if (samples >= blob_box_samples) {
            Samples:
            samples = 0;
            if (best == nullptr) {
                // Nothing suitable for merging.
                if (done) {
                    ubl->unified = true;
                }
                bb = ubl->index;
                ubl = nullptr;
                continue;
            }
            else {                    
                // Enough samples gathered, now merge with the best one.
                size_t bi = best->index;
                             
                size_t surface_1 = ubl->surface.size();
                size_t surface_2 = best->surface.size();
                double weight    = (surface_2) / double(surface_1 + surface_2);

                ubl->x = ((1.0 - weight) * ubl->x + weight * best->x);
                ubl->y = ((1.0 - weight) * ubl->y + weight * best->y);                
                ubl->r = ((1.0 - weight) * ubl->r + weight * best->r);
                ubl->g = ((1.0 - weight) * ubl->g + weight * best->g);
                ubl->b = ((1.0 - weight) * ubl->b + weight * best->b);
                ubl->a = ((1.0 - weight) * ubl->a + weight * best->a);                        
                
                for (it=best->surface.begin(); it!=best->surface.end(); ++it) {
                    f->owners[*it] = ubl;
                    ubl->surface.insert(*it);
                }                
         
                best->surface.clear();
                best->border .clear();
                ubl ->border .clear(); // Not needed because it's a pixel dust cloud anyway.
                
                f->blobs[bi] = nullptr; 
                
                delete best;
                best = nullptr;
                
                if (sz == fexp) {
                    // sz should always equal to fexp if everything 
                    // has been properly blobified before unifying.
                    // This optimization is only then possible.
                    if (f->blobs.back()) {
                        // Last element is not nullptr, swap it with blobs[bi]
                        // so that the latter can be safely popped later.
                        size_t last_index = f->blobs.back()->index;
                        f->blobs[bi] = f->blobs.back();
                        f->blobs[bi]->index = bi;
                        f->blobs[last_index] = nullptr;                        
                    }
                    
                    while (!f->blobs.empty() && f->blobs.back() == nullptr) {
                        f->blobs.pop_back();
                        sz--;
                        fexp--;
                        f->first_expansion--;
                    }
                }
                done = false;
                
                if (sz <= blob_number) {
                    ubl->unified = true;
                    return false;
                }
            }
            bb = ubl->index; // Next, take a blob right after ubl.            
            ubl = nullptr;
            continue;
        }
        
        if (ubl == nullptr) {        
            ubl = bl;
            ubl_x = round(ubl->x);
            ubl_y = round(ubl->y);
            continue;
        }
        
        samples++;
        
        std::pair<size_t, size_t> bounds_x = std::minmax(ubl_x, uint16_t(round(bl->x)));
        std::pair<size_t, size_t> bounds_y = std::minmax(ubl_y, uint16_t(round(bl->y)));
                        
        if (bounds_x.first + blob_box_grip >= bounds_x.second
        &&  bounds_y.first + blob_box_grip >= bounds_y.second) {
            size_t dx,dy,current_score;
            dx = (bounds_x.second - bounds_x.first);
            dy = (bounds_y.second - bounds_y.first);
            current_score = dx*dx+dy*dy;
            
            if (best == nullptr
            ||  current_score < score) {
                best = bl;
                score = current_score;
                continue;
            }            
        }
    }
    if (samples > 0 && ubl) {
        goto NotEnoughSamples;
    }
    
    return done;
}

bool thread::unify() {
    bool done = true;
    std::map<size_t, frame>::iterator it;
    
    for (it=frames.begin(); it!=frames.end(); ++it) {
        done = (done && unify_frame(it->first));
    }
    
    return done;
}

bool thread::match() {    
    if (blob_map == nullptr) {
        // Find the frame with the largest number of blobs.        
        std::map<size_t, frame>::iterator it;
        size_t blob_count = 0;
        size_t frame_count= 0;
        for (it=frames.begin(); it!=frames.end(); ++it) {
            frame *f = &(it->second);
            if (blob_count < f->blobs.size()) {
                blob_count = f->blobs.size();
            }
            frame_count++;
        }       
        
        blob_map_w = blob_count;
        blob_map_h = frame_count;
        
        size_t i;
        blob_map = new (std::nothrow) blob** [blob_map_w];
        if (blob_map) {
            for (i = 0 ; i < blob_map_w ; ++i ) {
                blob_map[i] = new (std::nothrow) blob* [blob_map_h];
                if (blob_map[i] == nullptr) {
                    for(size_t j = 0; j<=i; ++j) {
                        delete [] blob_map[j];
                    }
                    delete [] blob_map;
                    blob_map = nullptr;
                    return false;
                }
            }            
            
            // Fill blob map with pointers to real blobs.
            size_t f=0;
            for (it=frames.begin(); it!=frames.end(); ++it) {
                frame *fp = &(it->second);
                
                // Add empty blobs if needed.
                while (fp->blobs.size() < blob_count) {
                    blob* new_blob          = new blob;
                    new_blob->x             = fp->x;
                    new_blob->y             = fp->y;
                    new_blob->r             = fp->r;
                    new_blob->g             = fp->g;
                    new_blob->b             = fp->b;
                    new_blob->a             = fp->a;
                    new_blob->unified       = true;
                    new_blob->index         = fp->blobs.size(); 
                    fp->blobs.push_back(new_blob);
                }                
                
                size_t blobs = fp->blobs.size();
                   
                std::shuffle(fp->blobs.begin(), fp->blobs.end(), e1);
                   
                for (size_t b=0; b<blobs; ++b) {
                    blob_map[b][f] = fp->blobs[b];
                    fp->blobs[b]->group = b;
                    fp->blobs[b]->index = b;
                }   
                ++f;             
            }
            
            blob_map_e = get_energy(blob_map);            
        }
        
        return false;        
    }
    
    if (blob_map_h == 0 || blob_map_w <= 1) return true;
    
    // Pick 2 blobs randomly and swap them if it would decrease the energy.
    std::uniform_int_distribution<size_t> uniform_dist_x(0, blob_map_w - 1);
    std::uniform_int_distribution<size_t> uniform_dist_y(0, blob_map_h - 1);   
    size_t x1,x2;
    size_t y = uniform_dist_y(e1);
    size_t y_next = (y+1) % blob_map_h;
    size_t y_prev = (y>0 ? y-1 : blob_map_h - 1);

    x1 = uniform_dist_x(e1);    
    do {
        x2 = uniform_dist_x(e1);
    } while (x1 == x2);        
    
    blob* x1_y_prev = blob_map[x1][y_prev];
    blob* x2_y_prev = blob_map[x2][y_prev];
    blob* x1_y      = blob_map[x1][y     ];
    blob* x2_y      = blob_map[x2][y     ];
    blob* x1_y_next = blob_map[x1][y_next];
    blob* x2_y_next = blob_map[x2][y_next];

    bool x1_volatile = x1_y->surface.empty();
    bool x2_volatile = x2_y->surface.empty();    
    
    if (x1_volatile && x2_volatile) {
        return false; // No point to swap empty blobs.
    }    
    
    double x1_e_before = blob_distance(x1_y_prev, x1_y) + blob_distance(x1_y, x1_y_next);
    double x2_e_before = blob_distance(x2_y_prev, x2_y) + blob_distance(x2_y, x2_y_next);
    double x1_e_after  = blob_distance(x2_y_prev, x1_y) + blob_distance(x1_y, x2_y_next);
    double x2_e_after  = blob_distance(x1_y_prev, x2_y) + blob_distance(x2_y, x1_y_next);
    
    double c1 = x1_e_before + x2_e_before;
    double c2 = x1_e_after  + x2_e_after;

    if (c1 >= c2 || (degenerate && (counter % degenerate) == 0)) {
        blob *buf       = blob_map[x2][y];
        blob_map[x2][y] = blob_map[x1][y];
        blob_map[x1][y] = buf;                        
        double gain = c1 - c2;                            
        
        blob_map_e -= gain;
        if (blob_map_e < 0.0) blob_map_e  = 0.0;                        
        
        if (blob_map_e < best_e) {
            best_e          = blob_map_e;
            best_blob_map_e = best_e;
            counter         = 0;

            if (deviant) {
                // Refresh all groups.
                for (size_t j=0; j<blob_map_h; ++j) {
                    for (size_t i=0; i<blob_map_w; ++i) {
                        blob_map[i][j]->group = i;
                    }
                }
            }
            else {
                blob_map[x1][y]->group = x1;
                blob_map[x2][y]->group = x2;
            }
            deviant = false;
        }
        else deviant = true;
    }
    
    if (best_e == 0.0) return true;
    
    return false;
}

bool thread::init_morph() {
    size_t i,j;
    
    if (chains.empty() && blob_map) {
        chain_map_e = 0.0;        
        
        // Refresh blob map according to groups.
        // Needed in case deviant was left TRUE so 
        // that the blob map might be messed up.
        std::map<size_t, frame>::iterator it;
        size_t frame_index = 0;
        for (it=frames.begin(); it!=frames.end(); ++it) {
            std::vector<blob*> *blobs = &(it->second.blobs);
            size_t blob_count = blobs->size();
            for (i=0; i<blob_count; ++i) {
                blob_map[blobs->at(i)->group][frame_index] = blobs->at(i);
            }
            frame_index++;
        }            
        
        // Make sure that volatile blobs have their positions averaged.
        for (i=0; i<blob_map_w; ++i) {
            std::vector<blob *> to_be_fixed;
            for (j=0; j<blob_map_h; ++j) {
                to_be_fixed.push_back(blob_map[i][j]);
            }
            fix_volatiles(&to_be_fixed);
            to_be_fixed.clear();
        }
        
        for (i=0; i<blob_map_w; ++i) {
            size_t max_size=0;
            // Find out the largest blob.
            for (j=0; j<blob_map_h; ++j) {
                if (blob_map[i][j]->surface.size() > max_size) {
                    max_size = blob_map[i][j]->surface.size();
                }
            }
            size_t point_count = max_size * point_density;

            // Then create a blob chain that includes all the blobs that share the same group.
            for (j=0; j<blob_map_h; ++j) {               
                if (chains[i].points == nullptr) {                
                    if (!renew_chain( &(chains[i]), point_count, blob_map_h)) {
                        chains.erase(i);
                        return false;
                    }
                }
                
                // Remember the maximum surface area this blob can take:
                chains[i].max_surface = max_size;
                
                // Initialize the points in this chain at key frame j.
                pixel px;
                size_t p, pos;
                std::vector<point> points;
                blob  *bl = blob_map[i][j];
                point *pt;
                std::set<size_t>::iterator it;
                size_t frame_key;

                std::map<size_t, frame>::iterator fit = frames.begin();
                assert(frames.size() > j);
                std::advance(fit, j);
                frame_key = fit->first;

                std::uniform_int_distribution<uint8_t> uniform_dist_byte(0, UINT8_MAX); 

                p = 0;
                for (it=bl->surface.begin(); it!=bl->surface.end(); ++it) {
                    pos = *it;
                    px  = get_pixel(frame_key, pos);
                    
                    pt = &(chains[i].points[p][j]);
                    pt->s.x       = px.x;
                    pt->s.y       = px.y;
                    pt->s.flags   = HAS_PIXEL|HAS_FLUID;
                    pt->s.x_fract = 0;
                    pt->s.y_fract = 0;
                    points.push_back(*pt);                    
                    
                    ++p;
                }

                if (points.empty()) {
                    point volatile_point;
                    volatile_point.s.x      = std::round(bl->x);
                    volatile_point.s.y      = std::round(bl->y);
                    volatile_point.s.flags  = HAS_FLUID; // This is still needed for the trajectories of fluid particles.
                    volatile_point.s.x_fract= 0;
                    volatile_point.s.y_fract= 0;
                    
                    chains[i].points[p][j] = volatile_point;                    
                    points.push_back(volatile_point);                    
                    ++p;
                }

                std::shuffle(points.begin(), points.end(), e1);
                size_t psz = points.size();
                
                for (; p<chains[i].width; ++p) {
                    pt = &(chains[i].points[p][j]);
                    pt->s.x       =  points[p%psz].s.x;
                    pt->s.y       =  points[p%psz].s.y;
                    pt->s.flags   =  points[p%psz].s.flags;
                    pt->s.flags  &= ~HAS_FLUID; // Duplicates must not represent fluid particles.
                    pt->s.x_fract =  uniform_dist_byte(e1);
                    pt->s.y_fract =  uniform_dist_byte(e1);                    
                }     
                
#ifdef ATOMORPH_OPENCV
                // Initiate Locality Sensitive Hashing:
                //cv::flann::LinearIndexParams indexParams;
                cv::flann::KDTreeIndexParams indexParams(16);
                cv::Mat * features = new cv::Mat(chains[i].width, 2, CV_32F);
                if (features) {
                    for (size_t f=0; f<chains[i].width; ++f) {
                        float x,y;
                        point2xy(chains[i].points[f][j], &x, &y);
                        features->at<float>(f, 0) = x;
                        features->at<float>(f, 1) = y;
                        chains[i].places[f][j]    = f;
                    }

                    cv::flann::Index *kdtree = new cv::flann::Index(*features, indexParams, cvflann::FLANN_DIST_L1);
                    chains[i].kdtrees[j] = (void *) kdtree;
                    chains[i].feature[j] = (void *) features;
                }
#endif
                chains[i].energy = get_energy(&(chains[i]));
                chain_map_e     += chains[i].energy;
            }
        }
        if (chains.empty()) return false;
    }

    // These extra checks are done here in advance to have the morph
    // function as fast as possible, without wasting any time on the
    // sanity checks.
    if (blob_map_h == 0 || blob_map_w == 0) return false;
    
    size_t max_w = 0;
    for (i=0; i<blob_map_w; ++i) {
        chain *ch = &(chains[i]);
        if (ch->height == 0) return false;
        if (ch->width > max_w) max_w = ch->width;
    }
    
    if (max_w <= 1) return false;
    
    return true;
}

#ifdef ATOMORPH_OPENCV
bool thread::morph() {
    size_t x1, x2;
    chain *ch = &(chains[counter % blob_map_w]);
    
    if (ch->width <= 1) return false;
    
    std::uniform_int_distribution<size_t> uniform_dist_x(0, ch->width - 1);
    std::uniform_int_distribution<size_t> uniform_dist_y(0, ch->height- 1);   

    size_t y      = uniform_dist_y(e1);
    size_t y_next = (y+1) % ch->height;
    size_t y_prev = (y>0 ? y-1 : ch->height - 1);

    x1 = uniform_dist_x(e1);
    do {
        x2 = uniform_dist_x(e1);
    } while (x1 == x2);

    cv::flann::Index *kdtree = (cv::flann::Index *) ch->kdtrees[y_next];
    
    const size_t knn = 2;      // Number of nearest neighbors to search for
    std::vector<float> query1; // Search near the source
    std::vector<float> query2; // Search near the destination
    std::vector<int>   index1(knn);
    std::vector<int>   index2(knn);	
    std::vector<float> dist_1(knn);
    std::vector<float> dist_2(knn);	
    float qx1, qx2, qy1, qy2;

    point2xy(ch->points[x1][y] ,&qx1, &qy1);
    query1.push_back(qx1);
    query1.push_back(qy1);
   
    point2xy(ch->points[x1][y_next] ,&qx2, &qy2);    
    query2.push_back(qx2);
    query2.push_back(qy2);

    kdtree->knnSearch(query1, index1, dist_1, knn, cv::flann::SearchParams(1));
    kdtree->knnSearch(query2, index2, dist_2, knn, cv::flann::SearchParams(1));

    std::vector<size_t> x2s; // candidates
    while (!index1.empty()) {
	    size_t xn = ch->places[index1.back()][y_next];
	    if (xn != x1) {
	        x2s.push_back(xn);
        }
        index1.pop_back();
    }
    while (!index2.empty()) {
        size_t xn = ch->places[index2.back()][y_next];
        if (xn != x1) {
            x2s.push_back(xn);
        }
        index2.pop_back();
    }	
    x2s.push_back(x2);
    std::uniform_int_distribution<size_t> uniform_dist_x2s(0, x2s.size() - 1);
    x2 = x2s[uniform_dist_x2s(e1)];

    point **map = ch->points;
    
    point x1_y_prev = map[x1][y_prev];
    point x2_y_prev = map[x2][y_prev];
    point x1_y      = map[x1][y     ];
    point x2_y      = map[x2][y     ];
    point x1_y_next = map[x1][y_next];
    point x2_y_next = map[x2][y_next];    
    
    double x1_e_before = point_distance(x1_y_prev, x1_y) + point_distance(x1_y, x1_y_next);
    double x2_e_before = point_distance(x2_y_prev, x2_y) + point_distance(x2_y, x2_y_next);
    double x1_e_after  = point_distance(x2_y_prev, x1_y) + point_distance(x1_y, x2_y_next);
    double x2_e_after  = point_distance(x1_y_prev, x2_y) + point_distance(x2_y, x1_y_next);
    
    double c1 = x1_e_before + x2_e_before;
    double c2 = x1_e_after  + x2_e_after;    
 
    if (c1 >= c2) {
        // Make sure point map indexes are mapped correctly:
        size_t **plm = ch->places;
        plm[x2][y]   = x1;
        plm[x1][y]   = x2;

        point buf  = map[x2][y];
        map[x2][y] = map[x1][y];
        map[x1][y] = buf;                        
        double gain = c1 - c2;                            

        chain_map_e -= gain;
        if (chain_map_e < 0.0) chain_map_e  = 0.0;
    }    

    return false;
}
#else

std::mutex morph_mutex;
void morph_asynch(point **map, size_t width, size_t height, size_t repeat, unsigned rng_seed, double *gain) {
    std::mt19937 gen(rng_seed);

    std::uniform_int_distribution<size_t> uniform_dist_x(0, width -1);                
    std::uniform_int_distribution<size_t> uniform_dist_y(0, height-1);
                    
    size_t y      = uniform_dist_y(gen);
    size_t y_next = (y+1) % (height);
    size_t y_prev = (y>0 ? y-1 : height-1);     

    for (size_t i=0; i<repeat; ++i) {    
        size_t x1, x2;
        x1 = uniform_dist_x(gen);
        do {
            x2 = uniform_dist_x(gen);
        } while (x1 == x2);

        point x1_y_prev = map[x1][y_prev];
        point x2_y_prev = map[x2][y_prev];
        point x1_y      = map[x1][y     ];
        point x2_y      = map[x2][y     ];
        point x1_y_next = map[x1][y_next];
        point x2_y_next = map[x2][y_next];    
        
        double x1_e_before = point_distance(x1_y_prev, x1_y) + point_distance(x1_y, x1_y_next);
        double x2_e_before = point_distance(x2_y_prev, x2_y) + point_distance(x2_y, x2_y_next);
        double x1_e_after  = point_distance(x2_y_prev, x1_y) + point_distance(x1_y, x2_y_next);
        double x2_e_after  = point_distance(x1_y_prev, x2_y) + point_distance(x2_y, x1_y_next);
        
        double c1 = x1_e_before + x2_e_before;
        double c2 = x1_e_after  + x2_e_after;    
     
        if (c1 >= c2) {
            std::lock_guard<std::mutex> lock(morph_mutex);
            
            // Make sure no other thread has touched these points yet:
            if (sizeof(point) == 8) {
                if (map[x1][y].word != x1_y.word 
                ||  map[x2][y].word != x2_y.word) continue;
            }
            else {
                if (!std::memcmp(&map[x1][y], &x1_y, sizeof(point)) 
                ||  !std::memcmp(&map[x2][y], &x2_y, sizeof(point))) continue;
            }

            point buf  = map[x2][y];
            map[x2][y] = map[x1][y];
            map[x1][y] = buf;
            *gain += c1 - c2;
        }
    }    
}

bool thread::morph() {
    chain *ch = &(chains[counter % blob_map_w]);
    if (ch->width <= 1) return false;
    std::uniform_int_distribution<unsigned> uniform_dist_seed(0, std::numeric_limits<unsigned>::max());
    
    point **map = ch->points;    
    
    double gain = 0.0;
    if (thread_count == 0) {
        morph_asynch(map, ch->width, ch->height, cycle_length, uniform_dist_seed(e1), &gain);
    }
    else {
        std::vector<std::thread> workers;
        for (size_t i=0; i<thread_count; ++i) {
            workers.push_back(std::thread(morph_asynch, map, ch->width, ch->height, cycle_length, uniform_dist_seed(e1), &gain));
        }
        
        while (!workers.empty()) {
            workers.back().join();
            workers.pop_back();
        }
    }
    chain_map_e-=gain;
      
    return false;
}
#endif

size_t thread::get_blob_count (size_t f) {
    if (!has_frame(f)) return 0;
    return frames[f].blobs.size();
}

size_t thread::get_blob_count () {
    size_t count = 0;
    std::map<size_t, frame>::iterator it;
    
    for (it=frames.begin(); it!=frames.end(); ++it) {
        count += get_blob_count(it->first);
    }    
    
    return count;
}

double thread::get_energy(struct blob ***map) {
    if (blob_map_h == 0
    ||  blob_map_w == 0) return 0.0;
    
    blob *pframe_blob, *cframe_blob;
    double e=0.0;
    
    for (size_t i=0; i<blob_map_w; ++i) {
        pframe_blob = map[i][blob_map_h-1];
        
        for (size_t j=0; j<blob_map_h; ++j) {
            cframe_blob = map[i][j];
            
            e += blob_distance(pframe_blob, cframe_blob);
            
            pframe_blob = cframe_blob;
        }
    }

    return e;
}

double thread::get_energy(chain *ch){
    size_t w = ch->width;
    size_t h = ch->height;
    
    if (w  <= 1 || h == 0) return 0.0;
    
    point **map = ch->points;
    double e=0.0;
    for (size_t i=0; i<w; ++i) {
        for (size_t j=0; j<h; ++j) {
            size_t next = (j+1)%h;
            e += point_distance(map[i][j], map[i][next]);
        }
    }
    
    return e;
}

void thread::set_bbox(uint16_t x1, uint16_t y1, uint16_t x2, uint16_t y2) {
    if (running && !paused) return;    
    if (x1 > x2 || y1 > y2) return;
    
    bbox_x1 = x1;
    bbox_x2 = x2;
    bbox_y1 = y1;
    bbox_y2 = y2;
    
    int dx = bbox_x2 - bbox_x1;
    int dy = bbox_y2 - bbox_y1;
    bbox_d = dx*dx + dy*dy;
}

void thread::set_blob_weights(unsigned char rgba, unsigned char size, unsigned char xy) {
    if (running && !paused) return;    
    double sum = rgba*rgba + size*size + xy*xy;
    if (sum == 0.0) return;
    
    blob_rgba_weight = double(rgba*rgba)/sum; 
    blob_size_weight = double(size*size)/sum; 
    blob_xy_weight   = 1.0 - (blob_rgba_weight + blob_size_weight);
}

double thread::blob_distance(const blob *b1, const blob *b2) {
    size_t sz1 = b1->surface.size();
    size_t sz2 = b2->surface.size();
    size_t szs = sz1+sz2;
    double pix_dist = 0.0; 
    double col_dist = 0.0;
    double siz_dist = 0.0;
    
    if (szs > 0) {
        siz_dist = fabs(double(sz1)-sz2)/double(szs);
    }
    
    if (sz1 > 0 && sz2 > 0) {
        pixel p1,p2;
        p1 = create_pixel(b1->x, b1->y, create_color(b1->r, b1->g, b1->b, b1->a));
        p2 = create_pixel(b2->x, b2->y, create_color(b2->r, b2->g, b2->b, b2->a));
        pix_dist = sqrt(double(pixel_distance(p1, p2))/bbox_d);        
        col_dist = color_distance(p1.c, p2.c); 
    }

    return (blob_xy_weight   * pix_dist + 
            blob_rgba_weight * col_dist + 
            blob_size_weight * siz_dist);
}

double thread::get_energy() {
    if (!paused) return std::numeric_limits<double>::infinity();
    double e = -1.0;
    switch (state) {
        case STATE_BLOB_MATCHING: e = best_blob_map_e; break;        
        case STATE_ATOM_MORPHING: e = chain_map_e;     break;
        default:                                       break;
    }    
    return e;
}

void thread::fix_volatiles (std::vector<blob *> *to_be_fixed) {
    size_t sz = to_be_fixed->size();
    
    if (sz <= 1) return;
    
    size_t i  = 0;
    bool started = false;
    std::vector<blob *> volatiles;
    
    blob *first_static = nullptr;
    blob *previous_static = nullptr;

    while (1) {
        i = (i+1)%sz;
        blob *bl = to_be_fixed->at(i);
        bool empty = bl->surface.empty();
        
        if (!started) {
            if (empty) {
                if (i == 0) break; // All are empty.
                continue;
            }            
            started = true;
            first_static = bl;
            previous_static = bl;
            continue;
        }        
        
        if (empty) {
            volatiles.push_back(bl);
            continue;
        }
        
        if (!volatiles.empty()) {
            // Fix the range between previous_static and bl.
            size_t vsz = volatiles.size();
            for (size_t v=0; v<vsz; ++v) {
                double t = (v+1.0) / double(vsz+1.0);
                volatiles[v]->x = t*bl->x + (1.0 - t)*previous_static->x;
                volatiles[v]->y = t*bl->y + (1.0 - t)*previous_static->y;
            }
            volatiles.clear();
        }
        previous_static = bl;
        if (previous_static == first_static) break; // Cycle has ended.
    }
}

}

