/*
 * See Copyright Notice in atomorph.h
 */
 
#include "atomorph.h"

#ifdef ATOMORPH_DEPRECATED
const char * am_get_version() {
    return "0.51";
}

size_t am_get_warning() {
    size_t flags = 0;
    if (sizeof(void *) < 8) flags = flags|AM_WARN_POINTER_SIZE;
    if (sizeof(AM_ATOM)> 8) flags = flags|AM_WARN_ATOM_SIZE;
    return flags;
}

AM_SCENE::AM_SCENE() {
    atoms      = 0;
    frames     = 0;
    map        = nullptr;
    candidates = nullptr;    
    splines    = nullptr;
}

AM_SCENE::~AM_SCENE() {
    clear();
}

void AM_SCENE::clear(void) {
    if (candidates) {
        delete[] candidates;
        candidates = nullptr;
    }
    if (splines) {
        delete[] splines;
        splines = nullptr;
    }
    if (map) {
        for ( size_t i = 0; i < atoms; ++i) {
            free (map[i]);
        }
        free(map);
        map    = NULL;        
        atoms  = 0;
        frames = 0;        
    }
}

bool AM_SCENE::init(size_t atoms, size_t frames) {
    if (map != nullptr) return false;

    this->atoms  = atoms;
    this->frames = frames;    

    size_t i;
    map = (AM_ATOM **) malloc( atoms * sizeof(AM_ATOM *) );
    if (!map) return false;
    
    for (i = 0 ; i < atoms ; ++i ) {
        map[i] = (AM_ATOM *) malloc( frames * sizeof(AM_ATOM) );
        if (map[i] == nullptr) {
            for(size_t j = 0; j<=i; ++j) {
                free(map[j]);
            }
            free(map);
            return false;
        } 
    }
    
    candidates = new (std::nothrow) std::vector<AM_ATOM>[frames];    
    if (!candidates) {
        clear();
        return false;
    }
    for (i = 0; i < frames; ++i) candidates[i].reserve(atoms);    
    
    splines = new (std::nothrow) glnemo::CRSpline[atoms];    
    if (!splines) {
        clear();
        return false;
    }
    for (i = 0; i < atoms; ++i) splines[i].clearCPoints();
    
    return true;
}

void AM_SCENE::renew_splines() {
    for (size_t i = 0; i < atoms; ++i) {
        splines[i].clearCPoints();
        for (size_t j = 0; j < frames; ++j) {    
            double x = map[i][j].x / 65536.0;
            double y = map[i][j].y / 65536.0;
            splines[i].AddSplinePoint(glnemo::Vec3D(x, y, 0.0));       
        }    
    }
}

void AM_SCENE::get_xy(size_t atom, double t, double *x, double *y, unsigned interpolation) const {
    if (interpolation == AM_SPLINE) {
        glnemo::Vec3D v = splines[atom].GetInterpolatedSplinePoint(t);
        *x = v.x;
        *y = v.y;    
        return;
    }
 
    size_t cur_frame  = get_current_frame(t);
    AM_ATOM a1 = map[atom][cur_frame];
    if (interpolation == AM_LINEAR) {
        size_t next_frame = (cur_frame+1) % frames;
        AM_ATOM a2 = map[atom][next_frame];
        
        double p = 1.0/frames;
        double str = 1.0 - (t - p*cur_frame)/p;
        if (str > 1.0) str = 1.0;
        if (str < 0.0) str = 0.0;    
        
        *x = (str*a1.x + (1.0-str)*a2.x)/65535.0;
        *y = (str*a1.y + (1.0-str)*a2.y)/65535.0;
        return;
    }
    
    *x = a1.x / 65535.0;
    *y = a1.y / 65535.0;
}

// Returns 1.0 when interpolation at t is close to key frames.
// Returns 0.0 when interpolation at t is the farthest from key frames.
double AM_SCENE::get_certainty(double t) const {
    const double pi  = 3.14159265358;
    size_t cur_frame = get_current_frame(t);
    
    double p = 1.0/frames;
    double str = 1.0 - (t - p*cur_frame)/p;
    if (str > 1.0) str = 1.0;
    if (str < 0.0) str = 0.0;

    return ((cos(str*(pi/0.5))+1.0)/2.0);        
}

size_t AM_SCENE::get_current_frame(double t) const {
    size_t cur_frame = t*frames;
    cur_frame %= frames;
    return cur_frame;
}

AM_COLOR AM_SCENE::get_rgba(size_t atom, double t, double lag, double slope, unsigned interpolation) const {
    size_t cur_frame  = get_current_frame(t);
    AM_COLOR color;   
    
    if (interpolation == AM_NONE) {     
        color.r = map[atom][cur_frame].r;
        color.g = map[atom][cur_frame].g;
        color.b = map[atom][cur_frame].b;
        color.a = map[atom][cur_frame].a;                                    
        return color;    
    }
    
    const double pi   = 3.14159265358;
    size_t next_frame = (cur_frame+1) % frames;
   
    // Linear interpolation:
    double p = 1.0/frames;
    double str = 1.0 - (t - p*cur_frame)/p;
    if (str > 1.0) str = 1.0;
    if (str < 0.0) str = 0.0;    

    if (interpolation == AM_COSINE || interpolation == AM_PERLIN) {
        double s = (slope+0.1)/1.1;
        double l = (1.0 - s) * lag;
        
             if (str <=      l ) str = 0.0;
        else if (str >=   (l+s)) str = 1.0;
        else                     str = ((-cos((str-l)*(pi/s))+1.0)/2.0);    
    }
    
    color.r = (str)*map[atom][cur_frame].r + (1.0-str)*map[atom][next_frame].r;
    color.g = (str)*map[atom][cur_frame].g + (1.0-str)*map[atom][next_frame].g;
    color.b = (str)*map[atom][cur_frame].b + (1.0-str)*map[atom][next_frame].b;
    color.a = (str)*map[atom][cur_frame].a + (1.0-str)*map[atom][next_frame].a;
    return color;
}

double AM_SCENE::get_current_path_length(size_t atom, double t) const {
    size_t cur_frame = get_current_frame(t);
    size_t next_frame = (cur_frame+1) % frames;
    
    return am_atom_distance(map[atom][cur_frame], map[atom][next_frame]);    
}

bool AM_SCENE::push_atom(size_t frame, AM_ATOM atom) {
    if (frame >= frames) return false;
    size_t sz_before = candidates[frame].size();
    candidates[frame].push_back(atom);
    if (candidates[frame].size() == sz_before) return false;
    return true;
}

bool AM_SCENE::elect_atoms() {
    for (size_t i=0; i<frames; ++i) {
        if (candidates[i].size() == 0) {
            // If some frame has no candidates at all, artificially add
            // candidates to the center of the frame to avoid crashing.
            candidates[i].push_back(am_create_atom(0.5,0.5,0,0,0,0));
        }
        std::random_shuffle( candidates[i].begin(), candidates[i].end() );

        size_t clone=0; // If not enough candidates, start cloning them.        
        for (size_t j=candidates[i].size(); (j != 0 && j < atoms); ++j) {
            candidates[i].push_back( candidates[i].at(clone++) );
        }
        if (candidates[i].size() < atoms) return false;
        
        for (size_t a=0; a<atoms; ++a) {
            map[a][i].x = candidates[i].at(a).x;
            map[a][i].y = candidates[i].at(a).y;
            map[a][i].r = candidates[i].at(a).r;
            map[a][i].g = candidates[i].at(a).g;
            map[a][i].b = candidates[i].at(a).b;
            map[a][i].a = candidates[i].at(a).a;                                                            
        }        
    }
    return true;
}

double AM_SCENE::get_cost() const {
    double cost = 0.0;

    for (size_t a=0; a<atoms;  ++a) {        
        for (size_t f=0; f<frames; ++f) {
            size_t nf = (f+1) % frames;             
            cost += am_atom_distance(map[a][f], map[a][nf]);
        }
    }
    return cost;
}

double AM_SCENE::get_path_length(size_t atom) const {
    if (atom >= atoms) return -1.0;
    
    double distance = 0.0;
    for (size_t f=0; f<frames; ++f) {
        size_t nf = (f+1) % frames;             
        distance += am_atom_distance(map[atom][f], map[atom][nf]);
    }    
    
    return distance;
}

double AM_SCENE::get_path_color(size_t atom) const {
    if (atom >= atoms || frames==0) return -1.0;
    
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
    double a = 0.0;
    
    for (size_t f=0; f<frames; ++f) {
        r+=map[atom][f].r/255.0;
        g+=map[atom][f].g/255.0;
        b+=map[atom][f].b/255.0;
        a+=map[atom][f].a/255.0;                                
    }
    r/=double(frames);
    g/=double(frames);
    b/=double(frames);
    a/=double(frames);                                
                
    return r*g*b*a;
}

void AM_SCENE::shuffle() {
    size_t f,a;
    std::vector<AM_ATOM> row;
    row.reserve(atoms);
    for (f=0; f<frames; ++f) {
        row.clear();
        for (a=0; a<atoms; ++a) {
            row.push_back(map[a][f]);
        }                
        std::random_shuffle( row.begin(), row.end() );
        for (a=0; a<atoms; ++a) {
            map[a][f] = row[a];
        }
    }
}

size_t AM_SCENE::get_sorted_atom_at(size_t position) {
    if (position >= sorted_atoms.size()) return 0;
    return sorted_atoms[position];
}

bool am_compare_pairs( const std::pair<size_t, double>& i, const std::pair<size_t, double>& j ) {
    return j.second < i.second;
}

void AM_SCENE::sort_atoms() {
    std::vector< std::pair<size_t, double> > paths;
    paths.reserve(atoms);
    
    for (size_t a=0; a<atoms; ++a) {
        paths.push_back(std::make_pair(a, get_path_color(a)));
    }
    
    std::sort(paths.begin(), paths.end(), am_compare_pairs); 
    
    sorted_atoms.clear();
    size_t sz = paths.size();
    for (size_t a=0; a<sz; ++a) {
        sorted_atoms.push_back(paths[a].first);
    }
    
    return;
}

AM_ATOM AM_SCENE::get_atom(size_t atom, size_t frame) const {
    return map[atom][frame];
}

const std::vector<AM_ATOM> *AM_SCENE::get_candidates(size_t frame) const {
    return (const std::vector<AM_ATOM> *) &(candidates[frame]);
}

bool AM_SCENE::copy_map_from(const AM_SCENE *scene) {
    size_t atoms  = scene->atom_count();
    size_t frames = scene->frame_count();    
    if (atoms == 0 || frames == 0) return false;
    
    if (atom_count() != atoms || frame_count() != frames) {
        clear();
        if (!init(atoms, frames)) return false;
    }            
    
    for (size_t j=0; j<frames; ++j) {
        for (size_t i=0; i<atoms; ++i) {
            map[i][j] = scene->get_atom(i,j);
        }
    }
    return true;
}

bool AM_SCENE::copy_candidates_from(const AM_SCENE *scene) {
    if (frame_count() != scene->frame_count()
    ||  atom_count()  != scene->atom_count()) return false;
    
    for (size_t i=0; i<frames; ++i) {
        candidates[i] = *(scene->get_candidates(i));
    }
    return true;    
}

bool AM_SCENE::swap_atoms(size_t frame, size_t atom1, size_t atom2) {
    if (frame >= frames || atom1 >= atoms || atom2 >= atoms) return false;
    if (atom1 == atom2) return true;
    
    AM_ATOM abuf = map[atom1][frame];
    map[atom1][frame] = map[atom2][frame];
    map[atom2][frame] = abuf;
    
    return true;
}

AM_ATOM am_create_atom(double x, double y, unsigned char r, unsigned char g, unsigned char b, unsigned char a) {
    uint16_t real_x, real_y;
         if (x >= 1.0) real_x = 65535;
    else if (x <= 0.0) real_x = 0;
    else               real_x = round(x*65536.0);

         if (y >= 1.0) real_y = 65535;
    else if (y <= 0.0) real_y = 0;
    else               real_y = round(y*65536.0);
    
    AM_ATOM atom;
    atom.x = real_x;
    atom.y = real_y;
    atom.r = r;
    atom.g = g;
    atom.b = b;
    atom.a = a;
    return atom;
}

double am_atom_gradient(AM_ATOM a1, AM_ATOM a2) {
    double gradient = sqrt(pow(abs(a1.r-a2.r)/8.0, 2.0) + pow(abs(a1.g-a2.g)/8.0, 2.0) + pow(abs(a1.b-a2.b)/8.0, 2.0)) / 55.209119491;
    return (gradient < 0.0 ? 0.0 : (gradient > 1.0 ? 1.0 : gradient));
}

double am_atom_distance(AM_ATOM a1, AM_ATOM a2) {
    double cost  = 0.0;
    cost  = sqrt(pow(fabs((a1.x - a2.x)/256.0),2.0) + pow(fabs((a1.y - a2.y)/256.0),2.0)) / 362.033147696;
    return (cost < 0.0 ? 0.0 : (cost > 1.0 ? 1.0 : cost));
}

AM_THREAD::AM_THREAD() {
    running      = false;
    signal_stop  = false;
    signal_pause = false;
    paused       = true;
    subcost      = nullptr;
    cost         = 0.0;    
}

AM_THREAD::~AM_THREAD() {
    if (step_thread.joinable()) {
        pause();
        stop();
    }
    clear();
}

bool AM_THREAD::init(const AM_SCENE *scene) {
    if (running && !paused) return false;        
    if (!this->scene.copy_map_from(scene)) return false;        
    
    // Original candidates are copied from the input scene.
    // Then individual atoms for this thread get elected so
    // that the map of this thread would be different than
    // of the input scene's.
    if (!this->scene.copy_candidates_from(scene)) return false;    
    if (!this->scene.elect_atoms()) return false;    
    
    size_t atoms  = scene->atom_count();
    subcost = new (std::nothrow) double[atoms];    
    if (!subcost) return false;
    
    cost = this->scene.get_cost();    
    for (size_t i = 0; i < atoms; ++i) subcost[i]=0.0;

    std::default_random_engine e(0);
    e1 = e;
    
    this->step_size = 1000;
    magic_exponent = 3.0;
    
    gradient_importance = 0.0;

    return true;
}

void AM_THREAD::set_seed(unsigned seed) {
    if (running && !paused) return;
    std::default_random_engine e(seed);
    e1 = e;
}

// Sane range: [100, 1000]
// Very high number will make the threads to respond to signals very slowly.
void AM_THREAD::set_step_size(int step_size) {
    if (running && !paused) return;
    this->step_size = step_size;    
}

// Sane range: [1.0, 3.0]
void AM_THREAD::set_magic_exponent(double exponent) {
    if (running && !paused) return;
    magic_exponent = exponent;
}

// Sane range: [0.0, 1.0]
void AM_THREAD::set_gradient_importance(double importance) {
    if (running && !paused) return;
    gradient_importance = importance;
}

bool AM_THREAD::clear() {
    if (running && !paused) return false;
    scene.clear();
    
    if (subcost) {
        delete[] subcost;    
        subcost = nullptr;
    }
    
    return true;
}

double AM_THREAD::chain_length(AM_ATOM a1, AM_ATOM a2, AM_ATOM a3) {
    return am_atom_distance(a1, a2) + am_atom_distance(a2, a3);
}

double AM_THREAD::chain_gradient(AM_ATOM a1, AM_ATOM a2, AM_ATOM a3) {
    return am_atom_gradient(a1, a2) + am_atom_gradient(a2, a3);
}

void AM_THREAD::step() {
    size_t frames = scene.frame_count();
    size_t atoms  = scene.atom_count();

    int max_tries=step_size;
    
    if (frames == 0 || atoms <= 1) return;

    std::uniform_int_distribution<size_t> uniform_dist_atoms (0, atoms  - 1);
    std::uniform_int_distribution<size_t> uniform_dist_frames(0, frames - 1);
    
    double cost_per_atom  = cost / atoms;
    double cost_per_frame = cost_per_atom / frames;

    size_t frame = uniform_dist_frames(e1);
    size_t frame_before = (frame == 0 ? frames-1 : frame-1);
    size_t frame_after  = (frame + 1) % frames;
    size_t atom1 = 0, atom2;
    AM_ATOM a1b,a1c,a1a;
    double chain1b;
    double gradient1b=0.0;

    EvalAtom1:
    if (max_tries-- <= 0) return;
    atom1 = uniform_dist_atoms(e1);
    a1b = scene.get_atom(atom1, frame_before);
    a1c = scene.get_atom(atom1, frame);
    a1a = scene.get_atom(atom1, frame_after);        
    chain1b    = chain_length  (a1b,a1c,a1a);
    
    if (gradient_importance > 0.0) {
        gradient1b = chain_gradient(a1b,a1c,a1a);
    }

    if (chain1b/2.0 < cost_per_frame) goto EvalAtom1;
    max_tries++;
    Again:
    if (max_tries-- <= 0) return;
    atom2 = uniform_dist_atoms(e1);                              
    if (atom1 == atom2) goto Again;
    
    {    
        AM_ATOM a2b = scene.get_atom(atom2, frame_before);
        AM_ATOM a2c = scene.get_atom(atom2, frame);
        AM_ATOM a2a = scene.get_atom(atom2, frame_after);        
        
        double chain2b    = chain_length(a2b,a2c,a2a);       
        double chain1a    = chain_length(a1b,a2c,a1a);
        double chain2a    = chain_length(a2b,a1c,a2a);
        
        double gain = (chain1b - chain1a) + (chain2b - chain2a);                
        
        double c1a3 = pow(chain1a, magic_exponent);
        double c2a3 = pow(chain2a, magic_exponent);
        double c1b3 = pow(chain1b, magic_exponent);
        double c2b3 = pow(chain2b, magic_exponent);
        
        if (gradient_importance > 0.0) {        
            double gradient2b = chain_gradient(a2b,a2c,a2a);
            double gradient1a = chain_gradient(a1b,a2c,a1a);
            double gradient2a = chain_gradient(a2b,a1c,a2a);
            
            double inv = 1.0 - gradient_importance;
            c1a3 = c1a3*inv + gradient_importance*gradient1a;
            c2a3 = c2a3*inv + gradient_importance*gradient2a;
            c1b3 = c1b3*inv + gradient_importance*gradient1b;
            c2b3 = c2b3*inv + gradient_importance*gradient2b;
        }
        
        if ((c1a3 + c2a3) < (c1b3 + c2b3)) {
            cost -= gain;
            scene.swap_atoms(frame, atom1, atom2);        
        }
        else goto Again;           
    }
}

void AM_THREAD::run() {
    while (!signal_stop) {
        if (!signal_pause && !paused) {
            step();
        }
        else {
            paused = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(0));
        }
    }
    running      = false;
    signal_stop  = false;
    signal_pause = false;
    paused       = true;    
}

bool AM_THREAD::fetch_scene(AM_SCENE *target) const {
    if (running && !paused) return false;    
    return target->copy_map_from(&scene);
}

AM_IMAGE::AM_IMAGE() {
    w = 0;
    h = 0;
    done = true;
    
    running      = false;
    signal_stop  = false;
    signal_pause = false;
    paused       = true;    
}

AM_IMAGE::~AM_IMAGE() {
    if (step_thread.joinable()) {
        pause();
        stop();
    }
}

void AM_IMAGE::run() {
    while (!signal_stop) {
        if (!signal_pause && !paused) {
            render();
            if (done) paused = true;            
        }
        else {
            paused = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(0));
        }
    }
    running      = false;
    signal_stop  = false;
    signal_pause = false;
    paused       = true;    
    color_interpolation = AM_NONE;
    path_interpolation  = AM_NONE;
}

void AM_IMAGE::render() {
    if (done) return;
    
    if (w == 0 || h == 0) {
        atoms.clear();
        done = true;
        return;
    }
    size_t atom_count = scene.atom_count();
    
    scene.renew_splines();
    atoms.clear(); 
    atoms.reserve(atom_count);    

    std::map<size_t, double> rs; // red components
    std::map<size_t, double> gs; // green components
    std::map<size_t, double> bs; // blue components
    std::map<size_t, double> as; // alpha components
    std::map<size_t, double> weight; // sum of weights

    for (size_t i=0; i<atom_count; ++i) {
        int x,y;
        double px,py;

        scene.get_xy(i, t, &px, &py, path_interpolation);
             if (px > 1.0) px = 1.0;
        else if (px < 0.0) px = 0.0;
             if (py > 1.0) py = 1.0;
        else if (py < 0.0) py = 0.0;
        x = std::min(size_t(round(px*w)),size_t(w-1));
        y = std::min(size_t(round(py*h)),size_t(h-1));

        AM_COLOR color;
        AM_ATOM atom = scene.get_atom(i, scene.get_current_frame(t));
        unsigned char r,g,b,a;
        
        if (color_interpolation == AM_PERLIN) {
            double        f = 8.0; // Frequency
            int     octaves = 8;   // Octaves
            double perlin_x = (atom.x / 65535.0)*f;
            double perlin_y = (atom.y / 65535.0)*f;
            double lag   = lag_map.  octaveNoise(perlin_x, perlin_y, octaves)*0.5 + 0.5;
            double slope = slope_map.octaveNoise(perlin_x,perlin_y,8)*0.5 + 0.5;
            color = scene.get_rgba(i, t, lag, slope, AM_PERLIN);
        }
        else color = scene.get_rgba(i, t, 0.5, 0.5, color_interpolation);
        
        r = color.r;
        g = color.g;
        b = color.b;
        a = color.a;
        
        size_t index = y*w + x;

        double d = 1.0 - scene.get_current_path_length(i, t);

        if (weight.find(index) == weight.end()) weight[index] = d;           else weight[index]+= d;        
        if (rs.find(index)     ==     rs.end()) rs    [index] = d*(r/255.0); else rs    [index]+= d*(r/255.0);   
        if (gs.find(index)     ==     gs.end()) gs    [index] = d*(g/255.0); else gs    [index]+= d*(g/255.0);
        if (bs.find(index)     ==     bs.end()) bs    [index] = d*(b/255.0); else bs    [index]+= d*(b/255.0);
        if (as.find(index)     ==     as.end()) as    [index] = d*(a/255.0); else as    [index]+= d*(a/255.0);                        
    }
    
    AM_ATOM atom;
    for (auto& kv : weight) {
        double p = kv.second;
        atom.r = 255*(rs[kv.first] / p);
        atom.g = 255*(gs[kv.first] / p);
        atom.b = 255*(bs[kv.first] / p);
        atom.a = 255*(as[kv.first] / p);
        atom.x = kv.first % w;
        atom.y = kv.first / w;
        atoms.push_back(atom);
    }    

    done = true;
}

bool AM_IMAGE::set_scene(const AM_SCENE *scene) {
    if (running && !paused) return false;        
    if (!this->scene.copy_map_from(scene)) return false;                

    done = false;

    return true;
}

bool AM_IMAGE::set_resolution(size_t width, size_t height) {
    if (running && !paused) return false;        
    w = width;
    h = height;               

    return true;
}

bool AM_IMAGE::set_time(double time) {
    if (running && !paused) return false;        
    t = time;             

    return true;
}

bool AM_IMAGE::set_seed(unsigned seed) {
    if (running && !paused) return false;
    this->seed = seed;

    PerlinNoise l_map(seed);   lag_map   = l_map;
    PerlinNoise s_map(seed+1); slope_map = s_map;

    return true;
}

bool AM_IMAGE::set_color_interpolation(unsigned method) {
    if (running && !paused) return false;        
    color_interpolation = method;

    return true;
}

bool AM_IMAGE::set_path_interpolation(unsigned method) {
    if (running && !paused) return false;        
    path_interpolation = method;

    return true;
}

bool AM_IMAGE::get_xy(size_t pixel, int *x, int *y) const {
    if (running && !paused) return false;
    if (pixel >= atoms.size()) return false;
    AM_ATOM a = atoms[pixel];
    
    *x = a.x;
    *y = a.y;
    return true;
}

bool AM_IMAGE::get_rgba(size_t pixel, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *a) const {
    if (running && !paused) return false;        
    if (pixel >= atoms.size()) return false;
    AM_ATOM atom = atoms[pixel];
    
    *r = atom.r;
    *g = atom.g;
    *b = atom.b;
    *a = atom.a;
    return true;
}

bool AM_IMAGE::fetch_pixels(std::vector<AM_ATOM> *to) const {
    if (running && !paused) return false;
    if (to == nullptr) return false;
    
    size_t sz = atoms.size();
    to->reserve(sz);
    
    for (size_t i=0; i<sz; ++i) {
        to->push_back(atoms[i]);
    }
    return true;
}

// AM_BLENDER:

AM_BLENDER::AM_BLENDER() {
    w = 0;
    h = 0;
    done = true;
    
    running          = false;
    signal_stop      = false;
    signal_pause     = false;
    paused           = true;    
    median_combining = false;
}

AM_BLENDER::~AM_BLENDER() {
    if (step_thread.joinable()) {
        pause();
        stop();
    }
    clear();
}

bool AM_BLENDER::clear() {
    if (running && !paused) return false;
    atoms.clear();
    done = true;    
    layers.clear();
    
    return true;
}

void AM_BLENDER::run() {
    while (!signal_stop) {
        if (!signal_pause && !paused) {
            render();
            if (done) paused = true;            
        }
        else {
            paused = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(0));
        }
    }
    running      = false;
    signal_stop  = false;
    signal_pause = false;
    paused       = true;    
}

void AM_BLENDER::render() {
    if (done) return;
    size_t layer_count = layers.size();
        
    if (w == 0 || h == 0 || layer_count == 0) {
        atoms.clear();
        done = true;
        return;
    }

    size_t atom_count = atoms.size();
    
    std::map<size_t, size_t> counts; // number of atoms    
    std::map<size_t, std::vector<unsigned char> > rs; // red components for median combining
    std::map<size_t, std::vector<unsigned char> > gs; // green components for median combining
    std::map<size_t, std::vector<unsigned char> > bs; // blue components for median combining
    std::map<size_t, std::vector<unsigned char> > as; // alpha components for median combining
    std::map<size_t, size_t > rsa; // red components for averaging
    std::map<size_t, size_t > gsa; // green components for averaging
    std::map<size_t, size_t > bsa; // blue components for averaging
    std::map<size_t, size_t > asa; // alpha components for averaging    
    std::vector<size_t> indices;

    indices.reserve(atom_count);

    double weight = 1.0 / layer_count;

    for (size_t i=0; i<atom_count; ++i) {
        AM_ATOM atom = atoms[i];        
        size_t index = atom.y*w + atom.x;
        
        indices.push_back(index);
        
        if (median_combining) {        
            rs[index].push_back(atom.r);
            gs[index].push_back(atom.g);
            bs[index].push_back(atom.b);
            as[index].push_back(atom.a);
        }
        else {        
            if (rsa.find(index) == rsa.end()) rsa[index]=(atom.r)*weight; else rsa[index] += atom.r*weight;   
            if (gsa.find(index) == gsa.end()) gsa[index]=(atom.g)*weight; else gsa[index] += atom.g*weight;
            if (bsa.find(index) == bsa.end()) bsa[index]=(atom.b)*weight; else bsa[index] += atom.b*weight;
            if (asa.find(index) == asa.end()) asa[index]=(atom.a)*weight; else asa[index] += atom.a*weight;
        }                     
        
        if (counts.find(index) == counts.end()) counts[index] = 1; else counts[index]++;
    }
    
    atoms.clear();

    AM_ATOM atom;    
    size_t sz = indices.size();
    for (size_t i=0; i<sz; ++i) {    
        size_t index = indices[i];
        size_t x = index % w;
        size_t y = index / w;
        
        if (counts[index] < layer_count) {
            // If less than 50% of layers have a pixel on this position,
            // calculate the number of neighbors and if not enough found
            // skip this pixel.
            size_t yw = y*w;
            size_t n  = 0;
            if (          counts.find(yw+x+1) != counts.end()) n++;
            if (x != 0 && counts.find(yw+x-1) != counts.end()) n++;            
            if (          counts.find(yw+w+x) != counts.end()) n++;
            if (y != 0 && counts.find(yw-w+x) != counts.end()) n++;            
                        
            if (n < 3) continue;
        }
        atom.x = x;
        atom.y = y;
        
        size_t sz = counts[index];
        unsigned char r,g,b,a;

        if (median_combining) {
            std::sort(rs[index].begin(), rs[index].end());
            std::sort(gs[index].begin(), gs[index].end());
            std::sort(bs[index].begin(), bs[index].end());
            std::sort(as[index].begin(), as[index].end());
                                                
            if (sz % 2 == 0) {
                r = (rs[index][sz/2 - 1] + rs[index][sz/2])/2;
                g = (gs[index][sz/2 - 1] + gs[index][sz/2])/2;
                b = (bs[index][sz/2 - 1] + bs[index][sz/2])/2;
                a = (as[index][sz/2 - 1] + as[index][sz/2])/2;                                    
            }
            else {
                r = rs[index][sz/2];
                g = gs[index][sz/2];
                b = bs[index][sz/2];
                a = as[index][sz/2];                                    
            }                                
        }
        else {
            double p = double(counts[index])/layer_count;
            
            r = rsa[index]/p; if (r > 255) r = 255;
            g = gsa[index]/p; if (g > 255) g = 255;
            b = bsa[index]/p; if (b > 255) b = 255;
            a = asa[index]/p; if (a > 255) a = 255;    
        }
        
        atom.r = r;
        atom.g = g;
        atom.b = b;
        atom.a = a;
        atoms.push_back(atom);
    }    

    done = true;
}

bool AM_BLENDER::add_image(const AM_IMAGE *img) {
    if (running && !paused) return false;        
    
    size_t pixels_before = atoms.size();
    img->fetch_pixels(&atoms);   

    done = false;
    layers.push_back(pixels_before);

    return true;
}

bool AM_BLENDER::set_resolution(size_t width, size_t height) {
    if (running && !paused) return false;        
    w = width;
    h = height;               

    return true;
}

bool AM_BLENDER::set_median_combining(bool value) {
    if (running && !paused) return false;        
    median_combining = value;

    return true;
}

bool AM_BLENDER::get_xy(size_t pixel, int *x, int *y) const {
    if (running && !paused) return false;
    if (pixel >= atoms.size()) return false;    
    AM_ATOM a = atoms[pixel];
    
    *x = a.x;
    *y = a.y;
    return true;
}

bool AM_BLENDER::get_rgba(size_t pixel, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *a) const {
    if (running && !paused) return false;        
    if (pixel >= atoms.size()) return false;
    AM_ATOM atom = atoms[pixel];
    
    *r = atom.r;
    *g = atom.g;
    *b = atom.b;
    *a = atom.a;
    return true;
}

#endif

namespace am {

pixel create_pixel(uint16_t x, uint16_t y, unsigned char r, unsigned char g, unsigned char b, unsigned char a) {   
    pixel px;
    px.x   = x;
    px.y   = y;
    px.c.r = r;
    px.c.g = g;
    px.c.b = b;
    px.c.a = a;
    return px;
}

pixel create_pixel(uint16_t x, uint16_t y, color c) {
    return create_pixel(x, y, c.r, c.g, c.b, c.a);
}

const char * get_version() {
    return "1.0";
}

size_t get_warning() {
    size_t flags = 0;
    if (sizeof(void *) < 8) flags = flags|WARN_POINTER_SIZE;
    if (sizeof(pixel)  > 8) flags = flags|WARN_PIXEL_SIZE;
    if (sizeof(point)  > 8) flags = flags|WARN_POINT_SIZE;    
    return flags;
}

bool uses_opencv() {
#ifdef ATOMORPH_OPENCV
	return true;
#endif
    return false;
}

// Clears the chain, freeing dynamically allocated memory.
// After calling this, chain can be safely deleted.
void clear_chain(chain *c) {
    if (c->points) {
        for (size_t i=0; i<c->width; ++i) {
            delete [] c->points[i];
        }
        delete [] c->points;
        c->points = nullptr;
    }
    
    if (c->places) {
        for (size_t i=0; i<c->width; ++i) {
            delete [] c->places[i];
        }
        delete [] c->places;
        c->places = nullptr;
    }    
    
    if (c->splines) {
        delete [] c->splines;
        c->splines = nullptr;
    }

#ifdef ATOMORPH_OPENCV    
    if (c->kdtrees) {
        for (size_t i=0; i<c->height; ++i) {
            if (c->kdtrees[i] == nullptr) continue;
            cv::flann::Index *kdtree = (cv::flann::Index *) c->kdtrees[i];
            delete kdtree;
        }
        delete [] c->kdtrees;
        c->kdtrees = nullptr;
    }
    
    if (c->feature) {
        for (size_t i=0; i<c->height; ++i) {
            if (c->feature[i] == nullptr) continue;
            cv::Mat *m = (cv::Mat *) c->feature[i];
            delete m;
        }
        delete [] c->feature;
        c->feature = nullptr;
    }    
#endif

    c->width =0;
    c->height=0;
    c->energy=0.0;
    c->max_surface=0;
}

// Returns false if memory allocation has failed, chain will be left in a cleared state.
bool renew_chain(chain *c, size_t width, size_t height) {
    clear_chain(c);
    
    c->points = new (std::nothrow) point* [width];
    if (c->points == nullptr) {
        return false;
    }
    
    for (size_t k = 0; k < width; ++k ) {
        c->points[k] = new (std::nothrow) point [height];
        if (c->points[k] == nullptr) {
            for (size_t t = 0; t<=k; ++t) {
                delete [] c->points[t];
            }
            delete [] c->points;
            // Unable to allocate enough memory.
            c->points = nullptr;
            return false;
        }
    }        

    c->places = new (std::nothrow) size_t* [width];
    if (c->places == nullptr) {
        clear_chain(c);
        return false;
    }
    
    for (size_t k = 0; k < width; ++k ) {
        c->places[k] = new (std::nothrow) size_t [height];
        if (c->places[k] == nullptr) {
            for (size_t t = 0; t<=k; ++t) {
                delete [] c->places[t];
            }
            delete [] c->places;
            // Unable to allocate enough memory.
            c->places = nullptr;
            return false;
        }
    }   

    c->splines = new (std::nothrow) glnemo::CRSpline[width];
    if (c->splines == nullptr) {
        clear_chain(c);
        return false;
    }
    
#ifdef ATOMORPH_OPENCV
    c->kdtrees = new (std::nothrow) void* [height];
    if (c->kdtrees == nullptr) {
        clear_chain(c);
        return false;
    }

    for (size_t i=0; i<height; ++i) {
        c->kdtrees[i] = nullptr;
    }

    c->feature = new (std::nothrow) void* [height];
    if (c->feature == nullptr) {
        clear_chain(c);
        return false;
    }

    for (size_t i=0; i<height; ++i) {
        c->feature[i] = nullptr;
    }
#endif

    c->width = width;
    c->height= height;
        
    return true;
}

}


