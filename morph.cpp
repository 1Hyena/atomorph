/*
 * See Copyright Notice in atomorph.h
 */

#include "atomorph.h"

namespace am {

morph::morph() {   
    std::default_random_engine e(0);
    e1 = e;    
}

morph::~morph() {
    clear();
}

void morph::clear() {
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
    
    identifier= 0;
    energy    = 0.0;
    
    bbox_x1 = UINT16_MAX;
    bbox_y1 = UINT16_MAX;
    bbox_x2 = 0;
    bbox_y2 = 0;
    
    if (worker.is_running()) {
        worker.stop();
        worker.clear();
    }
    
    if (fluid) delete fluid;
}

void morph::compute() {
    if (is_busy()) return;    
    
    if (worker.is_running()) worker.resume();
    else                     worker.start();
}

void morph::iterate(size_t iterations) {
    if (is_busy()) return;
    if (worker.is_running()) worker.resume(iterations);
    else                     worker.start (iterations);   
}

void morph::compute(double seconds) {
    if (is_busy()) return;
    
    if (worker.is_running()) worker.resume(seconds);
    else                     worker.start(seconds);
}

void morph::suspend() {
    while (is_busy()) {
        // Attempt to suspend the worker.
        worker.pause();        
        std::this_thread::sleep_for(std::chrono::milliseconds(0));
    }    
}

bool morph::suspend(double timeout) {
    if (!is_busy()) return true;
    
    std::chrono::steady_clock::time_point t_start, t_end;    
    t_start = std::chrono::steady_clock::now();            
    
    while (1) {
        // Attempt to suspend the worker.
        worker.pause();
        if (worker.is_paused()) break;
        
        t_end = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start).count() > (timeout * 1000000000)) {
            return false;
        }
        
        std::this_thread::sleep_for(std::chrono::milliseconds(0));
    }    
    
    return true;
}

bool morph::synchronize() {
    if (!worker.is_paused()) return false;

    if (worker.get_seed() != seed) worker.set_seed(seed);
    
    worker.set_blob_delimiter  (blob_delimiter);
    worker.set_blob_threshold  (blob_threshold);
    worker.set_blob_max_size   (blob_max_size);
    worker.set_blob_min_size   (blob_min_size);
    worker.set_blob_box_grip   (blob_box_grip);
    worker.set_blob_box_samples(blob_box_samples);
    worker.set_blob_number     (blob_number);
    worker.set_degeneration    (degeneration);
    worker.set_density         (density);
    worker.set_threads         (threads);
    worker.set_cycle_length    (cycle_length);

    worker.set_bbox(bbox_x1, bbox_y1, bbox_x2, bbox_y2);
    worker.set_blob_weights(blob_rgba_weight, blob_size_weight, blob_xy_weight);
    
    if (identifier != worker.get_identifier()) {
        // Reset everything.
        worker.clear();
        
        std::map<size_t, frame>::iterator it;
        
        refresh_frames();
        for (it=frames.begin(); it!=frames.end(); ++it) {
            worker.set_frame(it->first, &(it->second));
        }
        worker.set_identifier(identifier);

        if (fluid) delete fluid;
        fluid = nullptr;
    }
    
    {
        // Load results.
        state  = worker.get_state();
        energy = worker.get_energy();
        
        std::map<size_t, frame>::iterator fit;        
        for (fit=frames.begin(); fit!=frames.end(); ++fit) {
            frame *to   = &(fit->second);
            frame *from = worker.get_frame(fit->first); 
            
            while (!to->blobs.empty()) {
                delete to->blobs.back();
                to->blobs.pop_back();
            }                   
            to->owners.clear();
                        
            for (size_t i=0; i<from->blobs.size(); ++i) {
                blob* b= from->blobs[i];
                if (!b) continue;
                
                blob* new_blob = new blob;                
                
                new_blob->index         = to->blobs.size();
                new_blob->r             = b->r;          
                new_blob->g             = b->g;
                new_blob->b             = b->b;
                new_blob->a             = b->a;
                new_blob->x             = b->x;
                new_blob->y             = b->y;                                                                                
                new_blob->surface       = b->surface;
                new_blob->group         = b->group;
                
                to->blobs.push_back(new_blob);
            }
        }
        
        const std::map<size_t, chain> *chains_from = worker.get_chains();        
        std::map<size_t, chain>::const_iterator cit;
        // Get chain updates from worker's chains.
        for (cit=chains_from->begin(); cit!=chains_from->end(); ++cit) {
            size_t chain_key = cit->first;
            
            // If this chain already exists, it might have different attributes.
            if ((chains[chain_key].width  != chains_from->at(chain_key).width
              || chains[chain_key].height != chains_from->at(chain_key).height)) {
                // Chains are different, free the target to allocate memory again later.
                if (!renew_chain(&(chains[chain_key]), cit->second.width, cit->second.height)) {
                    chains.erase(chain_key);
                    return false;
                }
            }
            
            // Normally the above memory allocation is not needed because the morph should 
            // not change its parameters very often, although the implementation allows it 
            // by essentially restarting the whole procedure whenever critical parameters 
            // change. For normal synchronization, only the below code is called to refresh
            // the key points and splines for example.
            for (size_t y=0; y<chains[chain_key].height; ++y) {
                for (size_t x=0; x<chains[chain_key].width; ++x) {
                    chains[chain_key].points[x][y] = chains_from->at(chain_key).points[x][y];
                    
                    double vx = chains[chain_key].points[x][y].s.x + chains[chain_key].points[x][y].s.x_fract / double(UINT8_MAX+1);
                    double vy = chains[chain_key].points[x][y].s.y + chains[chain_key].points[x][y].s.y_fract / double(UINT8_MAX+1);
                    
                    if (y == 0) chains[chain_key].splines[x].clearCPoints();
                    chains[chain_key].splines[x].AddSplinePoint(glnemo::Vec3D(vx, vy, 0.0));
                }
            }
            if (chains[chain_key].max_surface != chains_from->at(chain_key).max_surface) {
                // Number of particles in the fluid simulator has changed, reset the simulation.
                chains[chain_key].max_surface = chains_from->at(chain_key).max_surface;
                if (fluid) {
                    delete fluid;
                    fluid = nullptr;
                }
            }
        }
        // Delete chains that don't exist in worker's chains.
        std::vector<size_t> delete_keys;
        for (cit=chains.begin(); cit!=chains.end(); ++cit) {
            if (chains_from->find(cit->first) != chains_from->end()) continue;
            delete_keys.push_back(cit->first);
        }
        while (!delete_keys.empty()) {
            clear_chain(&(chains[delete_keys.back()]));
            chains.erase(delete_keys.back());
            delete_keys.pop_back();
        }
        
        worker.set_identifier(identifier);
        
        if (skip_state) {
            worker.next_state();
            skip_state = false;
        }
        
        // Initialize the fluid simulator if needed.
        if (!fluid && width > 0 && height > 0) {
            size_t particle_count = 0;
            
            for (cit=chains.begin(); cit!=chains.end(); ++cit) {
                particle_count += cit->second.max_surface;  
            }
            
            fluid = new (std::nothrow) FluidModel(width*sparseness+20, height*sparseness+20, particle_count);
            if (fluid) {                                
                Particle *particles = fluid->getParticles();
                
                for (size_t i=0; i<particle_count; ++i) {
                    Particle *p = &(particles[i]);
                    
                    p->clear();
                    p->x = width *sparseness/2.0 - 10;
                    p->y = height*sparseness/2.0 - 10;
                    p->u = 0.0;
                    p->v = 0.0;
                    p->gravity_x = p->x;
                    p->gravity_y = p->y;
                    p->freedom_r = 1.0;
                    p->R = 0.0; p->r = p->R;
                    p->G = 1.0; p->g = p->G; // Should never appear in morphs.
                    p->B = 0.0; p->b = p->B;
                    p->A = 1.0; p->a = p->A;
                    p->active   = false;
                    p->mature   = false;
                }
            }
        }
    }        
    
    return true;
}

void morph::refresh_frames() {
    std::map<size_t, frame>::iterator it;
    size_t index = 0;
    for (it=frames.begin(); it!=frames.end(); ++it) {
        it->second.index = index++;
        if (it->second.pixels.empty()) {
            it->second.x = (bbox_x1 <= bbox_x2 ? (bbox_x1 + bbox_x2)/2.0 : 0.0);
            it->second.y = (bbox_y1 <= bbox_y2 ? (bbox_y1 + bbox_y2)/2.0 : 0.0);
            it->second.r = 0.0;
            it->second.g = 0.0;
            it->second.b = 0.0;
            it->second.a = 0.0;
        }
    }
}

bool morph::add_frame(size_t frame_key) {
    if (frames.find(frame_key) != frames.end()) return false;
     
    frames[frame_key].x = (bbox_x1 <= bbox_x2 ? (bbox_x1 + bbox_x2)/2.0 : 0.0);
    frames[frame_key].y = (bbox_y1 <= bbox_y2 ? (bbox_y1 + bbox_y2)/2.0 : 0.0);
    
    refresh_frames();
    
    return true;
}

bool morph::add_pixel(size_t frame_key, pixel px) {
    if (blob_delimiter == HSP) {
        px.c = rgb_to_hsp(px.c);
    }

    add_frame(frame_key);

    am::frame *f = &(frames[frame_key]);
    
    auto pixelmap = &(f->pixels);
    size_t pos = px.y * (UINT16_MAX+1) + px.x;   
    
    (*pixelmap)[pos] = px;
    identifier++;
    
    if (px.x < bbox_x1) bbox_x1 = px.x;
    if (px.y < bbox_y1) bbox_y1 = px.y;
    if (px.x > bbox_x2) bbox_x2 = px.x;
    if (px.y > bbox_y2) bbox_y2 = px.y;
    
    if (pixelmap->size() == 1) {
        f->x = px.x;
        f->y = px.y;        
        f->r = px.c.r / 255.0;
        f->g = px.c.g / 255.0;
        f->b = px.c.b / 255.0;
        f->a = px.c.a / 255.0;
        return true;
    }
    
    double weight = 1.0 / double(pixelmap->size());
    f->x = (1.0 - weight)*f->x + weight*px.x;
    f->y = (1.0 - weight)*f->y + weight*px.y;
    f->r = (1.0 - weight)*f->r + weight*(px.c.r/255.0);
    f->g = (1.0 - weight)*f->g + weight*(px.c.g/255.0);
    f->b = (1.0 - weight)*f->b + weight*(px.c.b/255.0);
    f->a = (1.0 - weight)*f->a + weight*(px.c.a/255.0);    
    
    return true;
}

pixel morph::get_average_pixel (size_t f) {
    pixel px = create_pixel(0,0,0,0,0,0);
    if (!has_frame(f)) return px;
    px.x   = round(frames[f].x);
    px.y   = round(frames[f].y);
    px.c.r = round(255.0*frames[f].r);
    px.c.g = round(255.0*frames[f].g);
    px.c.b = round(255.0*frames[f].b);
    px.c.a = round(255.0*frames[f].a);
    
    if (blob_delimiter == HSP) {
        px.c = hsp_to_rgb(px.c);
    }
                        
    return px;
}

pixel morph::get_average_pixel (size_t f, size_t b) {
    pixel px = create_pixel(0,0,0,0,0,0);
    if (!has_frame(f) 
    ||  frames[f].blobs.size() <= b 
    ||  frames[f].blobs[b] == nullptr) {
        return px;
    }
    
    px.x   = round(frames[f].blobs[b]->x);
    px.y   = round(frames[f].blobs[b]->y);
    px.c.r = round(255.0*frames[f].blobs[b]->r);
    px.c.g = round(255.0*frames[f].blobs[b]->g);
    px.c.b = round(255.0*frames[f].blobs[b]->b);
    px.c.a = round(255.0*frames[f].blobs[b]->a);
    
    if (blob_delimiter == HSP) {
        px.c = hsp_to_rgb(px.c);
    }    
    
    return px;    
}

pixel morph::get_pixel (size_t f, size_t pos) {
    pixel px = create_pixel(0,0,0,0,0,0);
    if (!has_pixel(f,pos)) {
        // Asking a pixel from void gives
        // a fully transparent pixel.
        return px;
    }
    px = frames[f].pixels[pos];
        
    if (blob_delimiter == HSP) {
        px.c = hsp_to_rgb(px.c);
    }       
    
    return px;
}

size_t morph::get_pixel_count (size_t f) {
    if (!has_frame(f)) return 0;
    return frames[f].pixels.size();
}

size_t morph::get_blob_count (size_t f) {
    if (!has_frame(f)) return 0;
    return frames[f].blobs.size();
}

size_t morph::get_frame_key (double t) {
    double integ;
    t = std::modf(t, &integ);
    if (t < 0.0) t += 1.0;
    
    size_t f = t * frames.size();

    std::map<size_t, frame>::iterator it;
    
    size_t i=0;
    size_t frame = SIZE_MAX;
    for (it=frames.begin(); it!=frames.end(); ++it) {
        if (i++ == f) {
            frame = it->first;
            break;
        }
    }

    return frame;
}

size_t morph::get_blob_count () {
    size_t count = 0;
    std::map<size_t, frame>::iterator it;
    
    for (it=frames.begin(); it!=frames.end(); ++it) {
        count += get_blob_count(it->first);
    }    
    
    return count;
}

void morph::set_seed(unsigned seed) {
    std::default_random_engine e(seed);
    e1 = e;
    
    PerlinNoise l_map(seed);   lag_map   = l_map;
    PerlinNoise s_map(seed+1); slope_map = s_map;    
    
    this->seed = seed;    
}

double morph::normalize_time(double t) {
    double integ, time = std::modf(t, &integ);    
    if (time < 0.0) time += 1.0;
    return time;
}

const blob* morph::get_pixels(size_t blob_index, double time, std::vector<pixel> *to) {
    time = normalize_time(time);
    double t = time;

    if (frames.empty()) return nullptr;
    
    size_t frame_key  = get_frame_key(t);
    size_t blob_count = get_blob_count(frame_key);    
    
    if (blob_index >= blob_count) return nullptr;    

    double dt = 1.0 / double(frames.size());        
    t = std::max(0.0, (t - (frame_key * dt)) / dt);
    
    const blob * bl = get_blob(frame_key, blob_index);
    if (!bl) return nullptr;
    
    if (chains.empty()) {
        std::set<size_t>::iterator it;
        for (it=bl->surface.begin(); it!=bl->surface.end(); ++it) {
            pixel px = get_pixel(frame_key, *it);
            to->push_back(px);
        }
    }
    else {
        std::map<size_t, std::vector<color >> colors;
        std::map<size_t, std::vector<double>> weights;
        
        if (chains.find(bl->group) == chains.end()) return nullptr;
        
        size_t w = chains[bl->group].width;
        size_t h = chains[bl->group].height;
        point**p = chains[bl->group].points;
        
        size_t y = frames[frame_key].index;
        size_t y_next = (y+1)%h;
        assert(y < h);
        
        std::map<size_t, frame>::iterator it = frames.find(frame_key);
        ++it;
        if (it == frames.end()) it = frames.begin();
        size_t next_frame_key = it->first;
        
        for (size_t x=0; x<w; ++x) {
            pixel px1, px2;
            point pt1 = p[x][y];
            point pt2 = p[x][y_next];
            
            if (!point_has_pixel(pt1)
            &&  !point_has_pixel(pt2)) continue;
            
            if (point_has_pixel(pt1)
            && !point_has_pixel(pt2)) {
                px1 = get_pixel(frame_key, xy2pos(pt1.s.x, pt1.s.y));
                px2 = px1;
                px2.x = pt2.s.x;
                px2.y = pt2.s.y;
                px2.c.a = 0;
            }            
            else if (point_has_pixel(pt2)
                 && !point_has_pixel(pt1)) {
                px2 = get_pixel(next_frame_key, xy2pos(pt2.s.x, pt2.s.y));
                px1 = px2;
                px1.x = pt1.s.x;
                px1.y = pt1.s.y;
                px1.c.a = 0;
            }
            else {
                px1 = get_pixel(frame_key, xy2pos(pt1.s.x, pt1.s.y));
                px2 = get_pixel(next_frame_key, xy2pos(pt2.s.x, pt2.s.y));
            }
            point pt = pt1;
            
            if (motion == LINEAR) pt = interpolate(pt1, pt2, 1.0 - t);
            else if (motion == SPLINE && chains[bl->group].splines) {
                glnemo::Vec3D v = chains[bl->group].splines[x].GetInterpolatedSplinePoint(time);
                double fract, integ;
                fract = std::modf(v.x, &integ); pt.s.x = integ; pt.s.x_fract = std::round(fract * UINT8_MAX);
                fract = std::modf(v.y, &integ); pt.s.y = integ; pt.s.y_fract = std::round(fract * UINT8_MAX);
            }
                        
            pixel px;
            px.x = pt.s.x;
            px.y = pt.s.y;
            
            if (fading == PERLIN) {
                double        f = 8.0; // Frequency
                int     octaves = 8;   // Octaves
                double bbox_w   = bbox_x2 - bbox_x1 + 1.0;
                double bbox_h   = bbox_y2 - bbox_y1 + 1.0;
                double perlin_x = (((pt1.s.x-bbox_x1)*(UINT8_MAX+1)+pt1.s.x_fract) / double(bbox_w*(UINT8_MAX+1)))*f;
                double perlin_y = (((pt1.s.y-bbox_y1)*(UINT8_MAX+1)+pt1.s.y_fract) / double(bbox_h*(UINT8_MAX+1)))*f;
                double lag      = lag_map.  octaveNoise(perlin_x, perlin_y, octaves)*0.5 + 0.5;
                double slope    = slope_map.octaveNoise(perlin_x, perlin_y,       8)*0.5 + 0.5;

                px.c = interpolate(px1.c, px2.c, lag, slope, 1.0 - t);
            }
            else if (fading == COSINE) px.c = interpolate(px1.c, px2.c, 0.5, 0.5, 1.0 - t);
            else                       px.c = interpolate(px1.c, px2.c, 1.0 - t);

            if (px.x >= width || px.y >= height) {
                if (px.x > bbox_x2 || px.x < bbox_x1 
                ||  px.y > bbox_y2 || px.y < bbox_y1) continue;
            }
            
            // Bilinear interpolation:
            double total_weight= UINT8_MAX * UINT8_MAX;
            double weight_x1y1 = ((UINT8_MAX - pt.s.x_fract) * (UINT8_MAX - pt.s.y_fract)) / total_weight;
            double weight_x2y1 = (             pt.s.x_fract  * (UINT8_MAX - pt.s.y_fract)) / total_weight;
            double weight_x1y2 = ((UINT8_MAX - pt.s.x_fract) *              pt.s.y_fract ) / total_weight;
            double weight_x2y2 = (             pt.s.x_fract  *              pt.s.y_fract ) / total_weight;

            size_t pos = xy2pos(px.x, px.y);
            if (weight_x1y1 > 0.0) {
                colors [pos].push_back(px.c);
                weights[pos].push_back(weight_x1y1);
            }
            if ((px.x < bbox_x2 || px.x+1 < width) && weight_x2y1 > 0.0) {
                pos = xy2pos(px.x + 1, px.y);
                colors [pos].push_back(px.c);
                weights[pos].push_back(weight_x2y1);                
            }
            
            if ((px.y < bbox_y2 || px.y+1 < height) && weight_x1y2 > 0.0) {
                pos = xy2pos(px.x, px.y + 1);
                colors [pos].push_back(px.c);
                weights[pos].push_back(weight_x1y2);                
            }
            
            if (weight_x2y2 > 0.0) {
                if ((px.y   < bbox_y2 && px.x   < bbox_x2)
                ||  (px.y+1 < height  && px.x+1 < width)) {
                    pos = xy2pos(px.x + 1, px.y + 1);
                    colors [pos].push_back(px.c);
                    weights[pos].push_back(weight_x2y2);                
                }
            }
        }
        
        std::map<size_t, pixel> blob;
        std::map<size_t, std::vector<color>>::iterator cit;
        for (cit = colors.begin(); cit!=colors.end(); ++cit) {
            std::vector<color> *cs = &(cit->second);
            std::vector<double>*ws = &(weights[cit->first]);
            size_t vsz = cs->size();
            
            double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
            double weight_sum = 0.0;
            for (size_t c = 0; c<vsz; ++c) {
                weight_sum += ws->at(c);
                r += cs->at(c).r * ws->at(c);
                g += cs->at(c).g * ws->at(c);
                b += cs->at(c).b * ws->at(c);
                a += cs->at(c).a * ws->at(c);
            }
            double w = density > 0 ? double(vsz)/double(density) : 0.0;            
            if (w > 1.0) w = 1.0;

            r = std::round(r/weight_sum);
            g = std::round(g/weight_sum);
            b = std::round(b/weight_sum);
            a = std::round(w*(a/weight_sum));
            
            size_t pos = cit->first;
            uint16_t x = pos % (UINT16_MAX+1);
            uint16_t y = pos / (UINT16_MAX+1);
            
            pixel px = create_pixel(x,y,r,g,b,a);
            
            if (feather > 0) blob[pos] = px;
            else             to->push_back(px);
        }
        
        if (feather > 0) {
            std::vector<std::set<size_t>> layers;
            std::map<size_t, pixel> peeled_blob = blob;
            while (!peeled_blob.empty() && layers.size() < feather) {
                std::set<size_t> border;        
                std::map<size_t, pixel>::iterator pit;
                for (pit = peeled_blob.begin(); pit!=peeled_blob.end(); ++pit) {
                    size_t pos = pit->first;
                    uint16_t x = pos % (UINT16_MAX+1);
                    uint16_t y = pos / (UINT16_MAX+1);
                    
                    if (x == 0 || x == UINT16_MAX
                    ||  y == 0 || y == UINT16_MAX) {
                        border.insert(pos);
                        continue;
                    }
                    
                    if (peeled_blob.find(xy2pos(x+1, y  )) == peeled_blob.end()
                    ||  peeled_blob.find(xy2pos(x-1, y  )) == peeled_blob.end()
                    ||  peeled_blob.find(xy2pos(x,   y+1)) == peeled_blob.end()
                    ||  peeled_blob.find(xy2pos(x,   y-1)) == peeled_blob.end()) {
                        border.insert(pos);
                    }
                }
                std::set<size_t>::iterator bt;
                for (bt = border.begin(); bt != border.end(); ++bt) {
                    size_t pos = *bt;
                    peeled_blob.erase(pos);
                }
                
                layers.push_back(border);
            }
            
            size_t lsz = layers.size();
            for (size_t l=0; l<lsz; ++l) {
                std::set<size_t> *layer = &(layers[l]);
                std::set<size_t>::iterator lt;
                
                for (lt = layer->begin(); lt != layer->end(); ++lt) {
                    pixel px = blob[*lt];
                    px.c.a = std::round(double(px.c.a)*(double(l+1)/double(feather+1)));
                    to->push_back(px);
                }
            }
            
            std::map<size_t, pixel>::iterator pit;
            for (pit = peeled_blob.begin(); pit!=peeled_blob.end(); ++pit) {
                to->push_back(pit->second);
            }            
        }
    }
        
    return bl;
}

void morph::update_particle( Particle *p, const blob * blob_before, const blob * blob_after, point **points, 
                             std::map<size_t, size_t>& sources, std::map<size_t, size_t>& destinations, double t,
                             size_t y, size_t y_next, size_t frame_key, size_t next_frame_key, size_t chain_key, double time ) 
{
    double previous_surface = blob_before->surface.size();
    double next_surface     = blob_after->surface.size();

    size_t src_pos = p->source_pos;
    size_t dst_pos = p->destination_pos;
    size_t src_x   = sources[src_pos];
    size_t dst_x   = destinations[dst_pos];

    // If area is growing then gain the real source x from the destination.
    if (next_surface > previous_surface) {
        src_x = destinations[dst_pos];
    }

    {
        point pt1 = points[src_x][y];
        point pt2 = points[dst_x][y_next];
        
        pixel px1, px2;
        
        if (point_has_pixel(pt1)
        && !point_has_pixel(pt2)) {
            px1 = get_pixel(frame_key, xy2pos(pt1.s.x, pt1.s.y));
            px2 = px1;
            px2.x = pt2.s.x;
            px2.y = pt2.s.y;
            if (next_surface > 0.0) px2.c.a = 0;
        }            
        else if (point_has_pixel(pt2)
             && !point_has_pixel(pt1)) {
            px2 = get_pixel(next_frame_key, xy2pos(pt2.s.x, pt2.s.y));
            px1 = px2;
            px1.x = pt1.s.x;
            px1.y = pt1.s.y;
            if (previous_surface > 0.0) px1.c.a = 0;
        }
        else {
            px1 = get_pixel(frame_key, xy2pos(pt1.s.x, pt1.s.y));
            px2 = get_pixel(next_frame_key, xy2pos(pt2.s.x, pt2.s.y));
        }
        
        color c;
        
        if (fading == PERLIN) {
            double        f = 8.0; // Frequency
            int     octaves = 8;   // Octaves
            double bbox_w   = bbox_x2 - bbox_x1 + 1.0;
            double bbox_h   = bbox_y2 - bbox_y1 + 1.0;
            double perlin_x = (((pt1.s.x-bbox_x1)*(UINT8_MAX+1)+pt1.s.x_fract) / double(bbox_w*(UINT8_MAX+1)))*f;
            double perlin_y = (((pt1.s.y-bbox_y1)*(UINT8_MAX+1)+pt1.s.y_fract) / double(bbox_h*(UINT8_MAX+1)))*f;
            double lag      = lag_map.  octaveNoise(perlin_x, perlin_y, octaves)*0.5 + 0.5;
            double slope    = slope_map.octaveNoise(perlin_x, perlin_y,       8)*0.5 + 0.5;

            c = interpolate(px1.c, px2.c, lag, slope, 1.0 - t);
        }
        else if (fading == COSINE) c = interpolate(px1.c, px2.c, 0.5, 0.5, 1.0 - t);
        else                       c = interpolate(px1.c, px2.c, 1.0 - t);
        
        point pt = pt1;
        
        if (motion == LINEAR) pt = interpolate(pt1, pt2, 1.0 - t);
        else if (motion == SPLINE && chains[chain_key].splines) {
            glnemo::Vec3D v = chains[chain_key].splines[src_x].GetInterpolatedSplinePoint(time);
            double fract, integ;
            fract = std::modf(v.x, &integ); pt.s.x = integ; pt.s.x_fract = std::round(fract * UINT8_MAX);
            fract = std::modf(v.y, &integ); pt.s.y = integ; pt.s.y_fract = std::round(fract * UINT8_MAX);
        }
        
        float ptx, pty;
        point2xy(pt, &ptx, &pty);            

        p->gravity_x = std::min(std::max((ptx * sparseness) + 10.0, 1.0), double(width *sparseness + 10.0) );
        p->gravity_y = std::min(std::max((pty * sparseness) + 10.0, 1.0), double(height*sparseness + 10.0) );
        
        p->R = c.r/255.0;
        p->G = c.g/255.0;
        p->B = c.b/255.0;
        p->A = c.a/255.0;                        

        if (show_blobs == DISTINCT) {
            std::mt19937 gen(chain_key);
            std::uniform_real_distribution<double> uniform_dist_color(0.0, 1.0);
            p->R = uniform_dist_color(gen);
            p->G = uniform_dist_color(gen);
            p->B = uniform_dist_color(gen);
            p->A = 1.0;
        }
        else if (show_blobs == AVERAGE) {
            p->R = blob_before->r;
            p->G = blob_before->g;
            p->B = blob_before->b;
            p->A = blob_before->a;
        }

        p->r = t*p->R + (1.0-t)*p->r;
        p->g = t*p->G + (1.0-t)*p->g;
        p->b = t*p->B + (1.0-t)*p->b;
        p->a = t*p->A + (1.0-t)*p->a;

        if (!p->mature) {
            if (p->source_owner) {
                p->x = p->gravity_x;
                p->y = p->gravity_y;

                // These lines should only be called
                // for a newly created particle that was
                // first to occupy its source location:                                
                p->r = p->R;
                p->g = p->G;
                p->b = p->B;
                p->a = p->A;
                // Otherwise wrongly coloured single particles
                // start appearing during the morph when new
                // particles are created.
            }
            
            if (previous_surface == 0.0) p->strength = 0.1;
            else                         p->strength = 1.0;
            
            p->freedom_r = 1.0;
        }
    }
}

void morph::step_fluid(size_t frame_key, double t, double time) {
    // t is for current transition between 2 images, time is overall time across all key frames.
    size_t iptc = 0;
    size_t particle_count = fluid->get_particle_count();    
    Particle *particles   = fluid->getParticles();
    
    std::map<size_t, chain>::const_iterator cit;    
    for (cit=chains.begin(); cit!=chains.end(); ++cit) {
        size_t   chain_key  = cit->first;        
        size_t   w          = cit->second.width;
        size_t   h          = cit->second.height;
        point  **points     = cit->second.points;            
        
        size_t y = frames[frame_key].index;
        size_t y_next = (y+1)%h;
        
        std::map<size_t, frame>::iterator fit = frames.find(frame_key); ++fit;
        if (fit == frames.end()) fit = frames.begin();
        size_t next_frame_key = fit->first;
        
        const blob* blob_before = find_blob_by_group(frame_key,      chain_key);
        const blob* blob_after  = find_blob_by_group(next_frame_key, chain_key);

        double previous_surface = blob_before->surface.size();
        double next_surface     = blob_after->surface.size();        
        double current_surface  = std::round((1.0-t)*previous_surface + t*next_surface);
        
        size_t active_count = 0;
        size_t active_limit = current_surface;        
        
        {
            std::map<size_t, size_t> sources;      // pt1 has fluid, pt2 maybe not
            std::map<size_t, size_t> destinations; // pt2 has fluid, pt1 maybe not
            std::map<size_t, std::vector<Particle *>> particles_by_destination;
            std::map<size_t, std::vector<Particle *>> particles_by_source;

            std::vector<Particle *> free_particles;
            std::vector<Particle *> update_particles;
            
            for (size_t x=0; x<w && iptc < particle_count; ++x) {
                point pt1 = points[x][y];
                point pt2 = points[x][y_next];
                
                {
                    // Map this pixel as a possible destination.
                    // The vector it holds contains indexes for
                    // detailed lookups about this destination.
                    size_t destination_pos = xy2pos(pt2.s.x, pt2.s.y);
                    if (point_has_fluid(pt2)) {
                        destinations[destination_pos] = x;
                    }
                    
                    // Map this pixel as a possible source.
                    size_t source_pos = xy2pos(pt1.s.x, pt1.s.y);
                    if (point_has_fluid(pt1)) {
                        sources[source_pos] = x;
                    }
                }                 
                
                {            
                    Particle *p = &(particles[iptc++]);
                    if (p->active && p->frame_key == frame_key) {
                        // Map this particle as an existing particle.
                        active_count++;
                        particles_by_destination[p->destination_pos].push_back(p);
                        particles_by_source     [p->source_pos     ].push_back(p);
                        update_particles.push_back(p);
                    }
                    else {
                        p->active = false;
                        free_particles.push_back(p);
                    }
                }                
            }
            
            if (active_count < active_limit) {
                size_t to_create = active_limit - active_count;
                // Add new particles.
                
                while (to_create > 0) {
                    assert(!free_particles.empty());
                    Particle *p = free_particles.back();
                    free_particles.pop_back();
                    update_particles.push_back(p);
                    p->clear();
                    p->frame_key = frame_key;
                    
                    {
                        bool source_found = false;
                        std::map<size_t, size_t>::iterator it;
                        std::vector<size_t    > source_candidates;
                        
                        // First see if any of the source positions is still unused, 
                        // it is the highest priority to have these filled.                        
                        for (it=sources.begin(); it!=sources.end(); ++it) {
                            size_t pos = it->first;
                            
                            if (particles_by_source.find(pos) != particles_by_source.end()
                            && !particles_by_source[pos].empty()) {
                                // This source is already taken.
                                continue;
                            }
                            
                            // Free source was found, take it.
                            source_candidates.push_back(pos);
                        }
                        
                        if (!source_candidates.empty()) {
                            std::uniform_int_distribution<size_t> uniform_dist_candidate(0, source_candidates.size()-1);
                            size_t pos = source_candidates[uniform_dist_candidate(e1)];
                            assert(sources.find(pos) != sources.end());
                            size_t   x = sources[pos];
                            
                            p->source_pos = pos;
                            particles_by_source[pos].push_back(p);
                            source_found = true;
                         
                            {
                                // This source was unoccupied, take its default destination
                                // even if it leads to a place without HAS_FLUID flag and thus this
                                // particle must be destroyed sooner or later during this morph.
                                
                                point pt2 = points[x][y_next];
                                size_t destination_pos = xy2pos(pt2.s.x, pt2.s.y);
                                p->destination_pos = destination_pos;
                                particles_by_destination[destination_pos].push_back(p);
                                
                                assert(destinations.find(destination_pos) != destinations.end());
                            }
                        }
                        
                        if (!source_found) {
                            // All sources are already occupied, this new particle must now
                            // share the source with some other particle, but it will still
                            // must get an unused destination. Search for such destination.
                            bool destination_found = false;
                            for (it=destinations.begin(); it!=destinations.end(); ++it) {
                                size_t pos = it->first;
                                size_t   x = it->second;
                                
                                if (particles_by_destination.find(pos) != particles_by_destination.end()
                                && !particles_by_destination[pos].empty()) {
                                    // This destination is already taken.
                                    continue;
                                }
                                
                                // Free destination was found, take it.
                                p->destination_pos = pos;
                                particles_by_destination[pos].push_back(p);
                                destination_found = true;
                                
                                {
                                    // Destination was found but all the sources were already
                                    // occupied. Now see the pt1 of this destination even if
                                    // it does not have HAS_FLUID flag.
                                    
                                    point pt1 = points[x][y];
                                    size_t source_pos = xy2pos(pt1.s.x, pt1.s.y);
                                    p->source_pos = source_pos;
                                    particles_by_source[source_pos].push_back(p);
                                    source_found = true;
                                    
                                    assert(sources.find(source_pos) != sources.end());
                                }
                                
                                break;
                            }
                            
                            assert(destination_found);                            
                        }
                        
                        assert(source_found);
                    }
                    
                    {
                        size_t pos = p->source_pos;
                        
                        // Find a suitable starting position for this particle.
                        if (particles_by_source.find(pos) != particles_by_source.end()
                        &&  particles_by_source[pos].size() > 1) {
                            // Another particle exists with the same source,
                            // copy its location for seamless creation.

                            std::uniform_int_distribution<size_t> uniform_dist_p(0, particles_by_source[pos].size()-2);
                            Particle *parent = particles_by_source[pos].at(uniform_dist_p(e1));
                            assert(parent != p);
                            
                            p->copy_from(parent);
                            p->mature = parent->mature;
                            p->strength = parent->strength;
                            std::uniform_real_distribution<double> uniform_dist_fuzz(0.0, 1.0);
                            
                            p->x = std::min(p->x + uniform_dist_fuzz(e1), width * sparseness-10.0);
                            p->y = std::min(p->y + uniform_dist_fuzz(e1), height* sparseness-10.0);

                            p->source_owner = false;
                        }
                        else {
                            // No other particle currently exists with the
                            // same source location. Copy the attractor's
                            // current position for the initial setup of this
                            // new particle.
                            
                            p->mature = false;
                            p->source_owner = true;
                            
                            update_particle(p, blob_before, blob_after, points, sources, 
                                            destinations, t, y, y_next, frame_key, 
                                            next_frame_key, chain_key, time);
                        }
                    }
                    p->active = true;
                    to_create--;
                }
            }            
            else if (active_count > active_limit) {
                size_t to_delete = active_count - active_limit;
                // Delete particles.
                std::shuffle(update_particles.begin(), update_particles.end(), e1);
                while (to_delete > 0) {
                    assert(!update_particles.empty());
                    Particle *p = update_particles.back();
                    
                    p->active = false;
                    update_particles.pop_back();
                    
                    to_delete--;
                }
            }

            {
                // Particles either created or deleted by now, update their positions. 
                size_t usz = update_particles.size();
                for (size_t i=0; i < usz; ++i) {
                    Particle *p = update_particles[i];
                    if (!p->active) continue;
                    
                    update_particle(p, blob_before, blob_after, points, sources, 
                                    destinations, t, y, y_next, frame_key, 
                                    next_frame_key, chain_key, time);
                }
            }
        }
    }

    double freedom = std::sqrt(width*width + height*height) * double(sparseness) / 4.0;
    double freedom_radius = freedom*(1.0 - t);
    for (size_t step = 0; step < fluidsteps; ++step) fluid->step(fluidsteps - (step+1), freedom_radius, t);
}

void morph::draw_fluid(double time, std::vector<pixel> *image) {
    if (!fluid || chains.empty()) return;

    time = normalize_time(time);
    double t = time;

    if (frames.empty()) return;
    
    size_t frame_key = get_frame_key(t);
    double dt = 1.0 / double(frames.size());        
    t = std::max(0.0, (t - (frame_key * dt)) / dt);
    // t is 0.0 at the beginning of frame
    // t is 1.0 at the end of this frame      

    step_fluid(frame_key, t, time);

    std::map<size_t, color> image_buffer;

    Particle *particles   = fluid->getParticles();
    size_t particle_count = fluid->get_particle_count();
    
    std::map<size_t, std::vector<color     >> colors;
    std::map<size_t, std::vector<double    >> weights;    
    std::map<size_t, std::vector<Particle* >> startups;
    
    for (size_t i=0; i<particle_count; ++i) {
        Particle *p = &(particles[i]);
        if (p->active == false) continue;
        
        {
            // Correct the colours so they would not blur too much.
            // Blobs that appear from nothingness should be blended
            // with their surroundings a lot though to avoid sharp
            // gradient spots in the very beginnning of the morph.
            //double w = 1.0 / (100.0 * std::pow(t, 2.0)+1.0);
            double w = (1.0 / (100.0*std::pow((t-(1.0*(1.0-p->strength))), 2.0)+1.0)) * std::pow(p->strength, 2.0);

            p->r = (w * p->R)+(1.0 - w)*p->r;
            p->g = (w * p->G)+(1.0 - w)*p->g;
            p->b = (w * p->B)+(1.0 - w)*p->b;
            p->a = (w * p->A)+(1.0 - w)*p->a;
        }        
        
        double px, py;

             if (p->x <= 10.0)                    px = 0.0;
        else if (p->x >= width*sparseness + 10.0) px = (width-1);
        else                                      px = (p->x-10.0)/sparseness;
        
             if (p->y <= 10.0)                    py = 0.0;
        else if (p->y >= height*sparseness+ 10.0) py = (height-1);
        else                                      py = (p->y-10.0)/sparseness;

        uint16_t x = std::min(int(px), width - 1);
        uint16_t y = std::min(int(py), height- 1);        
        size_t pos = xy2pos(x, y);

        if (!p->mature) {
            p->mature = true;
        }

        double x_fract, y_fract;
        
        x_fract = std::modf(px, &px);
        y_fract = std::modf(py, &py);
        
        if (px >= width || py >= height) {
            if (px > bbox_x2 || px < bbox_x1 
            ||  py > bbox_y2 || py < bbox_y1) continue;
        }
        
        // Bilinear interpolation:
        double weight_x1y1 = ((1.0 - x_fract) * (1.0 - y_fract));
        double weight_x2y1 = (       x_fract  * (1.0 - y_fract));
        double weight_x1y2 = ((1.0 - x_fract) *        y_fract );
        double weight_x2y2 = (       x_fract  *        y_fract );

        color c = create_color(p->r, p->g, p->b, p->a);

        if (weight_x1y1 > 0.0) {
            colors [pos].push_back(c);
            weights[pos].push_back(weight_x1y1);
        }
        if ((x < bbox_x2 || x+1 < width) && weight_x2y1 > 0.0) {
            pos = xy2pos(x + 1, y);
            colors [pos].push_back(c);
            weights[pos].push_back(weight_x2y1);
        }
        
        if ((y < bbox_y2 || y+1 < height) && weight_x1y2 > 0.0) {
            pos = xy2pos(x, y + 1);
            colors [pos].push_back(c);
            weights[pos].push_back(weight_x1y2);
        }
        
        if (weight_x2y2 > 0.0) {
            if ((y   < bbox_y2 && x   < bbox_x2)
            ||  (y+1 < height  && x+1 < width)) {
                pos = xy2pos(x + 1, y + 1);
                colors [pos].push_back(c);
                weights[pos].push_back(weight_x2y2);
            }
        }
    }

    std::map<size_t, color> blob;
    std::map<size_t, std::vector<color>>::iterator cit;
    for (cit = colors.begin(); cit!=colors.end(); ++cit) {
        std::vector<color     > *cs = &(cit->second);
        std::vector<double    > *ws = &(weights[cit->first]);
        
        size_t vsz = cs->size();
        
        double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
        double weight_sum = 0.0;
        for (size_t c = 0; c<vsz; ++c) {        
            double weight = ws->at(c);

            weight_sum += weight;
            r += cs->at(c).r * weight;
            g += cs->at(c).g * weight;
            b += cs->at(c).b * weight;
            a += cs->at(c).a * weight;
        }

        r = std::round(r/weight_sum);
        g = std::round(g/weight_sum);
        b = std::round(b/weight_sum);
        a = std::round(a/weight_sum);

        if (keep_background) {
            double alpha = 1.0 - std::pow(t, 24.0);
            a *= alpha;
        }

        size_t pos = cit->first;
        uint16_t x = pos % (UINT16_MAX+1);
        uint16_t y = pos / (UINT16_MAX+1);
        
        pixel px = create_pixel(x,y,r,g,b,a);
        
        size_t imgpos = y*width + x;
        
        if (feather > 0) blob[pos] = px.c;
        else             image_buffer[imgpos] = px.c;
    }
    
    if (feather > 0) {        
        for (size_t layer = 0; layer < feather; ++layer) {
            double alpha = 1.0 - (double(layer + 1) / double(feather + 1));
            std::set<size_t> border;        
            std::map<size_t, color>::iterator it;
            for (it = blob.begin(); it!=blob.end(); ++it) {
                size_t pos = it->first;
                uint16_t x = pos % (UINT16_MAX+1);
                uint16_t y = pos / (UINT16_MAX+1);

                size_t ps;                
                if (x+1 <      width && blob.find( (ps = xy2pos(x+1, y  )) ) == blob.end()) border.insert(ps);
                if (y+1 <     height && blob.find( (ps = xy2pos(x,   y+1)) ) == blob.end()) border.insert(ps);
                if (x   >          0 && blob.find( (ps = xy2pos(x-1, y  )) ) == blob.end()) border.insert(ps);
                if (y   >          0 && blob.find( (ps = xy2pos(x,   y-1)) ) == blob.end()) border.insert(ps);
            }
            
            if (border.empty()) break;
            
            std::map<size_t, color> blob_buf;
            for (const auto& pos: border) {
                uint16_t x = pos % (UINT16_MAX+1);
                uint16_t y = pos / (UINT16_MAX+1);

                std::set<size_t> components;            
                size_t ps;                
                if (x+1 <      width && blob.find( (ps = xy2pos(x+1, y  )) ) != blob.end()) components.insert(ps);
                if (y+1 <     height && blob.find( (ps = xy2pos(x,   y+1)) ) != blob.end()) components.insert(ps);
                if (x   >          0 && blob.find( (ps = xy2pos(x-1, y  )) ) != blob.end()) components.insert(ps);
                if (y   >          0 && blob.find( (ps = xy2pos(x,   y-1)) ) != blob.end()) components.insert(ps);
                
                assert(!components.empty());
                
                double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
                double w = components.size();
                for (const auto& comp: components) {
                    color c = blob[comp];
                    r += (c.r/255.0) / w;
                    g += (c.g/255.0) / w;
                    b += (c.b/255.0) / w;
                    a += (c.a/255.0) / w;                                                
                }
                
                a *= alpha;
                
                blob_buf[pos] = create_color(r, g, b, a);
            }
            
            for (it = blob_buf.begin(); it!=blob_buf.end(); ++it) {
                blob[it->first] = it->second;
            }
        }
        
        std::map<size_t, color>::iterator it;
        for (it = blob.begin(); it!=blob.end(); ++it) {
            size_t pos = it->first;
            uint16_t x = pos % (UINT16_MAX+1);
            uint16_t y = pos / (UINT16_MAX+1);
            
            size_t imgpos = y*width + x;
            image_buffer[imgpos] = it->second;
        }
    }

    {
        std::map<size_t, color>::iterator cit;
        for (cit = image_buffer.begin(); cit!=image_buffer.end(); ++cit) {
            size_t pos = cit->first;
            color  c   = cit->second;

            if (keep_background) {
                color bgc = image->at(pos).c;
                double bgr = double(bgc.r)/255.0;
                double bgg = double(bgc.g)/255.0;
                double bgb = double(bgc.b)/255.0;
                double bga = double(bgc.a)/255.0;

                double r = c.r/255.0;
                double g = c.g/255.0;
                double b = c.b/255.0;
                double a = c.a/255.0;

                r  = a*r + (1.0-a)*bgr;
                g  = a*g + (1.0-a)*bgg;
                b  = a*b + (1.0-a)*bgb;
                a  = bga + (1.0-bga)*a;

                c.r = std::round(r*255.0);
                c.g = std::round(g*255.0);
                c.b = std::round(b*255.0);
                c.a = std::round(a*255.0);                                    
            }

            image->at(pos).c = c;
        }
    }
}

void morph::draw_atoms(double t, std::vector<pixel> *image) {
    std::map<size_t, std::vector<am::color >> colors;
    const am::blob *bl;    
    size_t b = 0;
    
    std::vector<pixel> pixels;

    while ( (bl = get_pixels(b++, t, &pixels)) != nullptr ) {
        size_t group  = bl->group;
        pixel  px_ave = blob2pixel(bl);
        
        while (!pixels.empty()) {
            pixel px   = pixels.back();            
            size_t pos = px.y*width + px.x;
            pixels.pop_back();            
            
            if (pos >= image->size()) continue;
            
            unsigned char rr,gg,bb,aa;
            if (show_blobs == DISTINCT) {
                std::mt19937 gen(group);
                std::uniform_int_distribution<unsigned char> uniform_dist_byte(0, 255);                
                rr = uniform_dist_byte(gen);
                gg = uniform_dist_byte(gen);
                bb = uniform_dist_byte(gen);
                aa = 255;
            }
            else if (show_blobs == AVERAGE) {
                rr = px_ave.c.r;
                gg = px_ave.c.g;
                bb = px_ave.c.b;
                aa = px_ave.c.a;
            }
            else {
                rr = px.c.r;
                gg = px.c.g;
                bb = px.c.b;
                aa = px.c.a;
            }
            
            if (aa == 0) continue;
            
            colors[pos].push_back(create_color(rr,gg,bb,aa));
        }
    }

    std::map<size_t, std::vector<color>>::iterator cit;
    for (cit = colors.begin(); cit!=colors.end(); ++cit) {
        size_t pos = cit->first;        
        std::vector<color> *cs = &(cit->second);
        size_t vsz = cs->size();
        
        double r    = 0.0, g    = 0.0, b    = 0.0, a    = 0.0;
        double rsum = 0.0, gsum = 0.0, bsum = 0.0, asum = 0.0;

        for (size_t c = 0; c<vsz; ++c) {
            double sr = cs->at(c).r/255.0;
            double sg = cs->at(c).g/255.0;
            double sb = cs->at(c).b/255.0;
            double sa = cs->at(c).a/255.0;

            asum += sa;
            rsum += sr*sa;
            gsum += sg*sa;
            bsum += sb*sa;
            
            if (c == 0) {
                r = sr;
                g = sg;
                b = sb;
                a = sa;
            }
            else {
                r  = sa*sr + (1.0-sa)*r;
                g  = sa*sg + (1.0-sa)*g;
                b  = sa*sb + (1.0-sa)*b;
                a  = a + (1.0-a)*(cs->at(c).a/255.0);
            }
        }
        
        if (blend_blobs) {
            r = rsum/asum;
            g = gsum/asum;
            b = bsum/asum;            
        }
        
        if (keep_background) {
            color bgc = image->at(pos).c;
            double bgr = double(bgc.r)/255.0;
            double bgg = double(bgc.g)/255.0;
            double bgb = double(bgc.b)/255.0;
            double bga = double(bgc.a)/255.0;
            
            r  = a*r + (1.0-a)*bgr;
            g  = a*g + (1.0-a)*bgg;
            b  = a*b + (1.0-a)*bgb;
            a  = bga + (1.0-bga)*a;
        }
        
        image->at(pos).c = create_color(r, g, b, a);
    }
}

void morph::get_pixels(double t, std::vector<pixel> *image) {
    image->clear();
    image->reserve(width*height);

    for (size_t y=0; y<height; ++y) {
        for (size_t x=0; x<width; ++x) {
            pixel px = create_pixel(x,y,0,0,0,0);
            if (keep_background) px.c = get_background(x, y, t);
            image->push_back(px);
        }
    }

    if (fluid && fluidsteps > 0) draw_fluid(t, image);
    else                         draw_atoms(t, image);
    
    return;
}

pixel morph::blob2pixel(const blob *bl) {
    pixel px = create_pixel(bl->x, bl->y, create_color(bl->r, bl->g, bl->b, bl->a));
    if (blob_delimiter == HSP) {
        px.c = hsp_to_rgb(px.c);
    }
    return px;
}

color morph::get_background(uint16_t x, uint16_t y, double time) {
    if (frames.empty()) return create_color(0.0, 0.0 ,0.0 ,0.0);
    size_t position = xy2pos(x, y);
    
    time = normalize_time(time);
    double t = time;
    
    size_t frame_key  = get_frame_key(t);
    std::map<size_t, frame>::iterator it = frames.find(frame_key);
    ++it;
    if (it == frames.end()) it = frames.begin();
    size_t next_frame_key = it->first;

    double dt = 1.0 / double(frames.size());        
    t = std::max(0.0, (t - (frame_key * dt)) / dt);

    color c1 = get_pixel(frame_key,      position).c;
    color c2 = get_pixel(next_frame_key, position).c;
    
    if (fading == PERLIN) {
        double f        = 8.0; // Frequency
        int    octaves  = 8;   // Octaves
        double bbox_w   = bbox_x2 - bbox_x1 + 1.0;
        double bbox_h   = bbox_y2 - bbox_y1 + 1.0;
        double perlin_x = ((x-bbox_x1) / double(bbox_w))*f;
        double perlin_y = ((y-bbox_y1) / double(bbox_h))*f;
        double lag      = lag_map.  octaveNoise(perlin_x, perlin_y, octaves)*0.5 + 0.5;
        double slope    = slope_map.octaveNoise(perlin_x, perlin_y,       8)*0.5 + 0.5;

        return interpolate(c1, c2, lag, slope, 1.0 - t);
    }
    else if (fading == COSINE) return interpolate(c1, c2, 0.5, 0.5, 1.0 - t);
    
    return interpolate(c1, c2, 1.0 - t);
}

color morph::interpolate(color c1, color c2, double lag, double slope, double str) {
    if (fading == COSINE || fading == PERLIN) {
        const double pi  = 3.14159265358;
        double s = (slope+0.1)/1.1;
        double l = (1.0 - s) * lag;
        
             if (str <=      l ) str = 0.0;
        else if (str >=   (l+s)) str = 1.0;
        else                     str = ((-cos((str-l)*(pi/s))+1.0)/2.0);    
    }

    return interpolate(c1, c2, str);
}

color morph::interpolate(color c1, color c2, double c1_weight) {
    color c;
    c.r = std::round(c1_weight * double(c1.r) + (1.0 - c1_weight) * double(c2.r));
    c.g = std::round(c1_weight * double(c1.g) + (1.0 - c1_weight) * double(c2.g));
    c.b = std::round(c1_weight * double(c1.b) + (1.0 - c1_weight) * double(c2.b));
    c.a = std::round(c1_weight * double(c1.a) + (1.0 - c1_weight) * double(c2.a));
    return c;
}

pixel morph::interpolate(pixel px1, pixel px2, double px1_weight) {
    pixel px;
    px.x   = std::round(px1_weight * double(px1.x)   + (1.0 - px1_weight) * double(px2.x));
    px.y   = std::round(px1_weight * double(px1.y)   + (1.0 - px1_weight) * double(px2.y));
    px.c.r = std::round(px1_weight * double(px1.c.r) + (1.0 - px1_weight) * double(px2.c.r));
    px.c.g = std::round(px1_weight * double(px1.c.g) + (1.0 - px1_weight) * double(px2.c.g));
    px.c.b = std::round(px1_weight * double(px1.c.b) + (1.0 - px1_weight) * double(px2.c.b));
    px.c.a = std::round(px1_weight * double(px1.c.a) + (1.0 - px1_weight) * double(px2.c.a));
    return px;
}

point morph::interpolate(point pt1, point pt2, double pt1_weight) {
    point pt;
    
    double x1 = (double(UINT8_MAX+1)*pt1.s.x + pt1.s.x_fract);
    double y1 = (double(UINT8_MAX+1)*pt1.s.y + pt1.s.y_fract);
    double x2 = (double(UINT8_MAX+1)*pt2.s.x + pt2.s.x_fract);
    double y2 = (double(UINT8_MAX+1)*pt2.s.y + pt2.s.y_fract);    
    double x  = (pt1_weight * x1 + (1.0 - pt1_weight) * x2);
    double y  = (pt1_weight * y1 + (1.0 - pt1_weight) * y2);
    
    pt.s.x = x / double(UINT8_MAX+1); pt.s.x_fract = x - (pt.s.x*(UINT8_MAX+1));
    pt.s.y = y / double(UINT8_MAX+1); pt.s.y_fract = y - (pt.s.y*(UINT8_MAX+1));
        
    return pt;
}

void morph::set_resolution(uint16_t w, uint16_t h) {
    width  = w;
    height = h;
}

double morph::get_time(size_t frame_index, size_t total_frames_out) {
    if (total_frames_out == 0) return 0.0;
    
    if (finite) {
        if (total_frames_out == 1) {
            if (get_frame_count() == 0) return 0.0;
            
            return (1.0 - (1.0/double(get_frame_count())))/2.0;
        }
    
        double t = double(frame_index) / double(total_frames_out-1);
        t*= 1.0 - (1.0/double(get_frame_count()));
        return t;
    }
    
    return double(frame_index) / double(total_frames_out);
}

blob * morph::find_blob_by_group(size_t frame_key, size_t group) {
    if (!has_frame(frame_key)) return nullptr;
    
    std::vector<blob*> *blobs = &(frames[frame_key].blobs);
    size_t bsz = blobs->size();
    
    for (size_t i=0; i< bsz; ++i) {
        if (blobs->at(i)->group == group) return blobs->at(i);
    }
    
    return nullptr;
}

void morph::set_fluid(unsigned f) {
    if (fluidsteps != f) {
        fluidsteps = f;
        if (fluidsteps > 0 && density > 1) {
            set_density(1);
        }
    }
}
    
void morph::set_density(uint16_t d) {
    if (density != d) {
        density = d;
        identifier++;
        if (fluidsteps != 0 && density > 1) {
            set_fluid(0);
        }
    }
}

}

