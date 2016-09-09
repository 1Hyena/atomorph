/*
 * See Copyright Notice in main.h
 */
#include <stdio.h>
#include <math.h>
#include <chrono>
#include <thread>

#include "main.h"
#include "lodepng.h"

MORPH_OPTIONS options;

int main(int argc, char **argv) {
    if (!init(argc, argv)) {
        fprintf(stderr, "Failed to initialize!\n");
        return -1;
    }    
    if (options.exit_flag) return 0;    

    std::chrono::steady_clock::time_point program_start;
    std::chrono::steady_clock::time_point program_end;
    
    program_start = std::chrono::steady_clock::now();
    {
        am::morph morph;
        
        // Options should be set before any operations with 
        // the morph instance. However, they can be changed
        // during the run time too.
        morph.set_blob_delimiter  (options.differ_blobs);
        morph.set_blob_max_size   (options.blob_max_size);
        morph.set_blob_min_size   (options.blob_min_size);    
        morph.set_blob_box_grip   (options.blob_box_grip);
        morph.set_blob_box_samples(options.blob_box_samples);    
        morph.set_blob_threshold  (options.blob_threshold);
        morph.set_blob_number     (options.blob_number);    
        morph.set_seed            (options.seed);
        morph.set_blob_rgba_weight(options.blob_rgba_weight);
        morph.set_blob_size_weight(options.blob_size_weight);
        morph.set_blob_xy_weight  (options.blob_xy_weight);        
        morph.set_degeneration    (options.degenerate);    
        morph.set_density         (options.density); // Higher than 1 sets fluid to 0.
        morph.set_motion          (options.motion);
        morph.set_fading          (options.fading);
        morph.set_threads         (options.threads);
        morph.set_cycle_length    (options.cycle_length);
        morph.set_feather         (options.feather);
        morph.set_keep_background (options.keep_background);
        morph.set_finite          (options.finite);
        morph.set_show_blobs      (options.show_blobs);
        morph.set_fluid           (options.fluid); // Higher than 0 sets density to 1.
        
        if (!load_files(&morph)) return -1;

        main_loop (&morph);
        save_files(&morph);
    }
    program_end = std::chrono::steady_clock::now();
    
    if (options.verbose) {
        size_t duration = std::chrono::duration_cast<std::chrono::microseconds>(program_end - program_start).count();
             if (duration <    1000) printf("Process took %lu microseconds to finish.\n", duration);
        else if (duration < 1000000) printf("Process took %lu milliseconds to finish.\n", duration/1000);
        else                         printf("Process took %lu seconds to finish.\n",      duration/1000000);
    }

    return 0;
}

bool init(int argc, char **argv) {
    if (true == (am::get_warning()&am::WARN_POINTER_SIZE)) {
        fprintf(stderr, "Pointer size is insufficiently small.\n");
    }
    if (true == (am::get_warning()&am::WARN_PIXEL_SIZE)) {
        fprintf(stderr, "Pixel size (%lu) is larger than optimal (%lu).\n",
            sizeof(am::pixel),
            sizeof(void *)
        );
    }
    if (true == (am::get_warning()&am::WARN_POINT_SIZE)) {
        fprintf(stderr, "Point size (%lu) is larger than optimal (%lu).\n",
            sizeof(am::point),
            sizeof(void *)
        );
    }
    if (am::uses_opencv()) {
        fprintf(stderr, "Experimental OpenCV optimizations are enabled.\n");
    }    
    if (argc == 1) {
        fprintf(stderr, "No arguments specified, try \"\e[1;33m%s --help\e[0m\".\n", argv[0]);
    }
    return options.parse(argc, argv);
}

bool fill_morph(am::morph *morph, size_t frame, std::vector<unsigned char> *image, unsigned width) {
    size_t sz = image->size();
    unsigned char r=0,g=0,b=0,a=0;
    size_t pixel = 0;
    bool empty = true;
    
    if (width > UINT16_MAX) return false;
    
    for (size_t j=0; j<sz; ++j) {
        switch (j%4) {
            case  0: r=image->at(j); break;
            case  1: g=image->at(j); break;
            case  2: b=image->at(j); break;
            default: {
                a=image->at(j);
                
                if (a == 0) {
                    pixel++; 
                    continue;
                }                        
                morph->add_pixel(frame, am::create_pixel(pixel%width, pixel/width, r, g, b, a));                        
                pixel++;
                empty = false;
                break;
            }
        }
    }
    
    if (empty) morph->add_frame(frame);
    
    return true;
}

void save_files(am::morph *morph) {
    for (unsigned f = 0; f<options.frames_out; ++f) {
        write_image(morph, f, morph->get_width(), morph->get_height());        
    }
}

bool load_files(am::morph *morph) {
    std::string buf;
    std::vector<unsigned char> image;
    unsigned width, height, error;    
    unsigned max_width=0, max_height=0;
    image.reserve(262144);
    size_t i;
    
    // Load input image files:
    for (i=0; i<options.files.size(); ++i) {
        buf=options.indir; buf.append("/"); buf.append(options.files[i]);
        
        if (options.verbose) printf("Loading %-30s ... ", buf.c_str());                
        error = lodepng::decode(image, width, height, buf.c_str());        
        if (error) {
            std::cerr << lodepng_error_text(error) << "." << std::endl;
            return false;
        }
        
        max_width  = std::max(max_width,  width);
        max_height = std::max(max_height, height);        
         
        if (options.verbose) printf("%ux%u image decoded.\n", width, height);

        if (!fill_morph(morph, i, &image, width) && options.verbose) {
            printf("Unable to fill %lu. frame.\n", i);
        }
        
        image.clear();
    }
    
    if (morph->get_frame_count()==0) {
        fprintf(stderr, "Error. Morph does not contain any key frames.\n");
        return false;
    }
    
    width = max_width;
    height= max_height;
    morph->set_resolution(width, height);

    return true;
}

void write_image(am::morph *morph, size_t frame_out, unsigned width, unsigned height) {
    double t = morph->get_time(frame_out, options.frames_out);

    char buf[1024];
    sprintf(buf, "%s/%s_%04lu.png", options.outdir.c_str(), options.file.c_str(), frame_out+1);
    if (options.verbose) printf("Rendering %4lu. frame ... ", frame_out+1);

    std::vector<unsigned char> image;
    image.reserve(4*width*height);
    image.resize (4*width*height, 0);

    std::vector<am::pixel> pixels;
    morph->get_pixels(t, &pixels);
    
    while (!pixels.empty()) {
        am::pixel px = pixels.back();
        pixels.pop_back();
                
        size_t pos = (px.y*width + px.x)*4;        
        if (pos >= image.size()) continue;
        
        image[pos + 0] = px.c.r;
        image[pos + 1] = px.c.g;
        image[pos + 2] = px.c.b;
        image[pos + 3] = px.c.a;
    }

    if (options.verbose) printf("writing %s ... ", buf);

    unsigned error = lodepng::encode(buf, image, width, height);
    if (error) {
        if (options.verbose) printf("[\e[1;31mFAIL\e[0m]\n");
        std::cerr << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
    }
    else if (options.verbose) printf("[\e[1;32mDONE\e[0m]\n");
}

void main_loop(am::morph *morph) {
    if (options.verbose) {
        printf("Blobifying %ux%u morph.\n", morph->get_width(), morph->get_height());
    }

    std::chrono::steady_clock::time_point start,end;    
    start = std::chrono::steady_clock::now();
    
    size_t blob_count    =   0;
    size_t largest_frame =   0;
    double last_energy   = 0.0;
    size_t match_blobs   =   0;
    size_t morph_atoms   =   0;
    size_t fs            = options.files.size();
    
    while (1) {
        morph->suspend();
        morph->synchronize();
        morph->compute();

        blob_count = 0;
        for (size_t i=0; i<fs; ++i) {
            size_t bc = morph->get_blob_count(i);
            if (bc >= blob_count) {
                blob_count = bc;
                largest_frame = i;
            }
        }
        
        unsigned morph_state = morph->get_state();
        
        end = std::chrono::steady_clock::now();     
        if (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() > 1000) {        
            if (morph_state == am::STATE_BLOB_DETECTION
            ||  morph_state == am::STATE_BLOB_UNIFICATION) {
                if (options.verbose) {
                    printf("%lu/%lu blob%s remaining on frame %lu.\n", 
                        blob_count, morph->get_blob_count(), 
                        blob_count == 1 ? "" : "s", largest_frame
                    );
                }
            }
            else if (morph_state == am::STATE_BLOB_MATCHING) {
                if (options.verbose) {
                    double e = morph->get_energy();
                    if (last_energy < e) printf("Matching blobs, absolute energy was %30.10f.\n", e);
                    else                 printf("Matching blobs, energy decreased by %30.10f.\n", last_energy - e);                    
                    last_energy = e;
                }
                if (++match_blobs >= options.match_time) {
                    morph->next_state(); 
                    last_energy = 0.0;
                }
            }
            else if (morph_state == am::STATE_ATOM_MORPHING) {
                if (options.verbose) {
                    double e = morph->get_energy();
                    if (last_energy < e) printf("Matching atoms, absolute energy was %30.2f.\n", e);
                    else                 printf("Matching atoms, energy decreased by %30.2f.\n", last_energy - e);                    
                    last_energy = e;
                }
                if (++morph_atoms >= options.morph_time) {
                    morph->next_state();
                }
            }
            else if (morph_state == am::STATE_DONE) {
                if (options.verbose) printf("All done!\n");
                break;
            }
            else {
                if (options.verbose) printf("Unknown state!\n");
                break;
            }
            start = std::chrono::steady_clock::now();
        }        
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
             
    if (options.verbose) {
        printf("Frame %lu had the most blobs (%lu).\n", largest_frame, blob_count);
    }
    
    return;
}

