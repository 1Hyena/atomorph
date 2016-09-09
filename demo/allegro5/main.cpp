/*
 * See Copyright Notice in main.h
 */
#include <stdio.h>
#include <math.h>

#include "main.h"

const float  FPS      =  60.0; // Maximum FPS.
const int    SCREEN_W =   800;
const int    SCREEN_H =   600;
const int    MORPH_W  =   128; // Width of the morph.  Should be at most the width of the input image.
const int    MORPH_H  =   128; // Height of the morph. Should be at most the height of the input image.
const int    ATOMS    = 10000; // Number of atoms used in one thread.
const size_t THREAD_N =     5; // Number of threads to use to find a perfect morph.
const size_t SLOWNESS =    50; // How many frames to render per animation cycle.

size_t morph_time  =     0;
int    view_frame  =     0;
bool pressed_keys[ALLEGRO_KEY_MAX];

size_t  active_thread    = 0;       // When render is ON, morph only one thread at a time.
int     color_fade       = AM_NONE; // Color interpolation method.
int     trajectory       = AM_NONE; // Atom trajectory interpolation method.
bool    median_combining = false;   // Noise reduction method. FALSE for averaging.
bool    stop_morphing    = false;   // To halt the morph time temporarily.
bool    no_render        = false;   // When TRUE no blending is done, just atom morphing.

ALLEGRO_DISPLAY     *display     = NULL;
ALLEGRO_EVENT_QUEUE *event_queue = NULL;
ALLEGRO_TIMER       *timer       = NULL;
ALLEGRO_FONT        *font        = NULL;
ALLEGRO_BITMAP      *morph_bmp   = NULL;   // Holds the final morph as a bitmap.
ALLEGRO_BITMAP      *thread_bmp[THREAD_N]; // Holds the results of the morphing threads.

// Helper function to initially populate the AM_SCENE object according to the provided
// image file.
bool fill_scene(AM_SCENE *scene, size_t frame, const char *png_file) {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    ALLEGRO_BITMAP * bmp = al_load_bitmap(png_file);
    int bmp_w = al_get_bitmap_width(bmp);
    int bmp_h = al_get_bitmap_height(bmp);

    al_lock_bitmap(bmp, ALLEGRO_PIXEL_FORMAT_ANY, ALLEGRO_LOCK_READONLY);

    for (int j=0; j<bmp_h; j++) {
        for (int i=0; i<bmp_w; i++) {
            ALLEGRO_COLOR c = al_get_pixel(bmp, i, j);

            unsigned char r,g,b,a;
            al_unmap_rgba(c, &r, &g, &b, &a);

            if (a == 0) continue;
            double px,py;
            px = double(i) / double(bmp_w);
            py = double(j) / double(bmp_h);
            scene->push_atom(frame, am_create_atom(px,py,r,g,b,a));
        }
    }

    al_unlock_bitmap(bmp);
    al_destroy_bitmap(bmp);
    return true;
}

int main(int argc, char **argv) {        
    if (!init(argc, argv)) {
        fprintf(stderr, "Failed to initialize!\n");
        return -1;
    }

    std::random_device rd;
    std::default_random_engine seed_engine(rd());
    std::uniform_int_distribution<unsigned> uniform_dist(1, std::numeric_limits<unsigned>::max());

    morph_bmp = al_create_bitmap(MORPH_W, MORPH_H);
    al_set_target_bitmap(morph_bmp);
    al_clear_to_color(al_map_rgba(0,0,0,0));

    AM_BLENDER blender; // Used to combine the thread results into the final morph.
    blender.set_resolution(MORPH_W, MORPH_H);
    blender.set_median_combining(median_combining);
    blender.start();

    AM_THREAD  scene_thread[THREAD_N]; // Each of these will morph its own version of the animation.
    AM_SCENE   scene_buf   [THREAD_N]; // Temporarily holds the last results of the morphing threads.
    AM_IMAGE   image_buf   [THREAD_N]; // Used to render the final image of the provided scene.

    {
        AM_SCENE scene; // Needed temporarily to store the raw input data.

        scene.init(ATOMS, 6); // Reserve 6 frames for this scene.
        fill_scene(&scene, 0, "../tests/data/battlelord_1.png");
        fill_scene(&scene, 1, "../tests/data/battlelord_2.png");
        fill_scene(&scene, 2, "../tests/data/battlelord_3.png");
        fill_scene(&scene, 3, "../tests/data/battlelord_4.png");
        fill_scene(&scene, 4, "../tests/data/battlelord_5.png");
        fill_scene(&scene, 5, "../tests/data/battlelord_6.png");

        for (size_t i=0; i<THREAD_N; ++i) {
            scene_buf   [i].init(scene.atom_count(), scene.frame_count());
            scene_thread[i].init(&scene); // Give the initial work to the worker thread.
            scene_thread[i].set_seed(uniform_dist(seed_engine)); // Seed for its RNG.
            scene_thread[i].set_step_size(100); // Number of iterations to make per step.
            scene_thread[i].set_magic_exponent(3.0); // Affects the generated trajectories.
            scene_thread[i].set_gradient_importance(0.0); // Importance of the color values.
            scene_thread[i].start();

            thread_bmp[i] = al_create_bitmap(MORPH_W, MORPH_H);
            al_set_target_bitmap(thread_bmp[i]);
            al_clear_to_color(al_map_rgba(0,0,0,0));

            image_buf[i].set_scene(&scene); // Prepares the image rendering object.
            image_buf[i].set_resolution(MORPH_W, MORPH_H);
            image_buf[i].set_seed(i); // Seed for its RNG (used in Perlin noise).
            image_buf[i].start();
        }
    }

    // Helper variables:
    bool   redraw         = true;
    bool   doexit         = false;
    bool   started        = false;
    bool   debug          = false;
    int    frame          = 0;
    double scene_cost     = 0.0; // Lower cost means better quality. Cost decreases over time.
    int    old_morph_time = 0;
    bool   refresh        = false;

    while(!doexit) {
        ALLEGRO_EVENT ev;
        al_wait_for_event(event_queue, &ev);

        if(ev.type == ALLEGRO_EVENT_TIMER) {
            redraw = true;
        }
        else if(ev.type == ALLEGRO_EVENT_DISPLAY_CLOSE) {
            break;
        }
        else if(ev.type == ALLEGRO_EVENT_KEY_DOWN) {
            pressed_keys[ev.keyboard.keycode] = true;
            switch(ev.keyboard.keycode) {
                case ALLEGRO_KEY_PAD_PLUS:  view_frame++; break;
                case ALLEGRO_KEY_PAD_MINUS: view_frame--; break;
                case ALLEGRO_KEY_D: debug = !debug; break;
                case ALLEGRO_KEY_M: median_combining = !median_combining; break;
                case ALLEGRO_KEY_S: stop_morphing = !stop_morphing; break;
                case ALLEGRO_KEY_R: no_render = !no_render;
                                    if (no_render) old_morph_time = morph_time;
                                    else {
                                        morph_time = old_morph_time;
                                        refresh = true;
                                    }
                                    break;
                case ALLEGRO_KEY_C:
                                         if (color_fade == AM_NONE)   color_fade = AM_LINEAR;
                                    else if (color_fade == AM_LINEAR) color_fade = AM_COSINE;
                                    else if (color_fade == AM_COSINE) color_fade = AM_PERLIN;
                                    else                              color_fade = AM_NONE;
                                    break;
                case ALLEGRO_KEY_T:
                                         if (trajectory == AM_NONE)   trajectory = AM_LINEAR;
                                    else if (trajectory == AM_LINEAR) trajectory = AM_SPLINE;
                                    else                              trajectory = AM_NONE;
                                    break;
                default: break;
            }
        }
        else if(ev.type == ALLEGRO_EVENT_KEY_UP) {
            pressed_keys[ev.keyboard.keycode] = false;
            switch(ev.keyboard.keycode) {
                case ALLEGRO_KEY_ESCAPE:
                    doexit = true;
                    break;
                case ALLEGRO_KEY_SPACE:
                    started = !started;

                    for (size_t i=0; i<THREAD_N; ++i) {
                        if (!scene_thread[i].is_running()) continue;
                        if (!started) scene_thread[i].pause();
                        else          scene_thread[i].resume();
                    }

                    break;
                default: break;
            }
        }
        else if(ev.type == ALLEGRO_EVENT_MOUSE_AXES ||
                ev.type == ALLEGRO_EVENT_MOUSE_ENTER_DISPLAY);
        else if(ev.type == ALLEGRO_EVENT_MOUSE_BUTTON_DOWN)  ;
        else if(ev.type == ALLEGRO_EVENT_MOUSE_BUTTON_UP)    ;

        if(redraw && al_is_event_queue_empty(event_queue)) {
            int fps;
            if ( (fps = calculate_fps()) == -1) continue;

            frame++;

            if (started) {
                size_t i;
                bool skip_render = false;

                if (morph_time%SLOWNESS == 0 || refresh) {
                    double cost       =0.0;

                    active_thread = (active_thread + 1) % THREAD_N;

                    for (i=0; i<THREAD_N; ++i) {
                        if (!scene_thread[i].is_paused()) {
                            // Before thread results can be read it must be paused.
                            scene_thread[i].pause();
                            cost += scene_thread[i].get_cost();
                            if (!no_render) {
                                scene_thread[i].fetch_scene(&(scene_buf[i]));
                            }
                        }

                        // When rendering is disabled all morphing threads will work,
                        // othwerwise only the active morphing thread works.
                        if (active_thread == i || no_render) scene_thread[i].resume();
                    }
                    scene_cost  = cost / THREAD_N;
                    skip_render = true;
                    refresh     = false;
                }
                else if (no_render) active_thread = (active_thread + 1) % THREAD_N;

                if (no_render) skip_render = true;

                bool slow_down = false; // Is set to TRUE when image rendering in not yet finished.
                                        // This is to slow down the animation rather than skip frames.

                if (!skip_render) {
                    bool all_images_done = true;
                    for (i=0; i<THREAD_N; ++i) {
                        if (!image_buf[i].is_paused()) {
                            all_images_done = false;
                            break;
                        }
                    }

                    if (all_images_done && blender.is_paused()) {
                        double t = (morph_time % SLOWNESS)/(double(SLOWNESS));

                        // Render blender image:
                        blend_morphs(&blender, morph_bmp);
                        blender.clear();
                        blender.set_median_combining(median_combining);

                        for (i=0; i<THREAD_N; ++i) {
                            // Render thread images:
                            render_morph(&image_buf[i], thread_bmp[i]);

                            // Give a new job to image blender thread:
                            blender.add_image(&image_buf[i]);

                            // Give new job to image thread:
                            image_buf[i].set_scene(&scene_buf[i]);
                            image_buf[i].set_time(t);
                            image_buf[i].set_color_interpolation(color_fade);
                            image_buf[i].set_path_interpolation(trajectory);
                            image_buf[i].resume();
                        }
                        blender.resume();
                    }
                    else slow_down = true; // Image rendering is lagging behind!
                }

                if (!stop_morphing) {
                    if (!slow_down && ++morph_time == SLOWNESS) morph_time = 0;
                }
            }

            redraw = false;

            double k = double(MORPH_H) /double(MORPH_W); // Aspect ratio.
            int merged_h = SCREEN_H - SCREEN_H/8;
            int merged_w = merged_h * k;

            al_set_target_bitmap(al_get_backbuffer(display));
            al_set_blender(ALLEGRO_ADD, ALLEGRO_ONE, ALLEGRO_ZERO);
            al_clear_to_color(al_map_rgb(128,128,128));
            al_set_blender(ALLEGRO_ADD, ALLEGRO_ONE, ALLEGRO_INVERSE_ALPHA);
            int mh = SCREEN_H / 8;
            int mw = mh * k;
            double s = double(SCREEN_W)/double(THREAD_N);

            // Draw the thread images to the upper side of the screen:
            for (size_t i=0; i<THREAD_N; ++i) {
                if (active_thread == i) {
                    al_draw_rectangle(s*i+s/2.0 - mw/2.0, 0.0,
                                      s*i+s/2.0 + mw/2.0, mh,
                                      al_map_rgb(192,0,0), 4.0);
                }
                al_draw_scaled_bitmap(thread_bmp[i], 0.0, 0.0, MORPH_W, MORPH_H,
                                      s*i+s/2.0 - mw/2.0, 0.0, mw, mh, 0);
            }

            // Draw the final morph to the center of the screen:
            al_draw_scaled_bitmap(morph_bmp, 0.0, 0.0, MORPH_W, MORPH_H,
                                  SCREEN_W/2 - merged_w/2, SCREEN_H - merged_h,
                                  merged_w, merged_h, 0);

            // Draw textual information:
            if (font!=NULL) {
                al_set_blender(ALLEGRO_ADD, ALLEGRO_ONE, ALLEGRO_INVERSE_ALPHA);

                al_draw_filled_rectangle(0, 0, SCREEN_W-SCREEN_W/4, 14, al_map_rgba(0,0,0,128));
                al_draw_textf(font, al_map_rgb(0,255,0), 0,  0, 0,
                    "FPS: %3d; Atoms: %d; Cost: %2.3f; Thread: %lu/%lu; Morph time: %lu;",
                    fps, ATOMS, scene_cost, active_thread+1, THREAD_N, morph_time
                );
                if (!started) {
                    al_draw_filled_rectangle(0, 0, SCREEN_W, SCREEN_H/4, al_map_rgba(0,0,0,128));

                    al_draw_textf(font, al_map_rgb(0,255,0),
                        SCREEN_W/2, SCREEN_H/8, ALLEGRO_ALIGN_CENTRE,
                        "AtoMorph v%s Demo by Erich Erstu, 2013",
                        am_get_version()
                    );

                    al_draw_filled_rectangle(0, SCREEN_H - SCREEN_H/4,
                        SCREEN_W, SCREEN_H,  al_map_rgba(0,0,0,128)
                    );

                    al_draw_textf(font, al_map_rgb(0,255,0), SCREEN_W/2,
                        SCREEN_H - SCREEN_H/8, ALLEGRO_ALIGN_CENTRE, "Press SPACE to start!"
                    );
                }
                else {
                    al_draw_filled_rectangle(0, 14, SCREEN_W/3, SCREEN_H/5, al_map_rgba(0,0,0,128));
                    al_draw_textf(font, al_map_rgb(0,255,0), 0, 12, 0,
                        "[M]edian combining: %s", median_combining ? "ON" : "OFF"
                    );

                    al_draw_textf(font, al_map_rgb(0,255,0), 0, 24, 0,
                        "[S]top morph time.%s", stop_morphing ? " (Stopped)" : ""
                    );

                    al_draw_textf(font, al_map_rgb(0,255,0), 0, 36, 0,
                        "[C]olor interpolation: %s",
                        color_fade == AM_NONE   ? "NONE"   :
                        color_fade == AM_LINEAR ? "LINEAR" :
                        color_fade == AM_COSINE ? "COSINE" :
                        color_fade == AM_PERLIN ? "PERLIN" : ""
                    );

                    al_draw_textf(font, al_map_rgb(0,255,0), 0, 48, 0,
                        "[T]rajectory interpolation: %s",
                        trajectory == AM_NONE   ? "NONE"   :
                        trajectory == AM_LINEAR ? "LINEAR" :
                        trajectory == AM_SPLINE ? "SPLINE" : ""
                    );

                    al_draw_textf(font, al_map_rgb(0,255,0), 0, 60, 0,
                        "[R]endering: %s", !no_render ? "ON" : "OFF");

                    al_draw_text (font, al_map_rgb(0,255,0), 0, 72, 0, "SPACE to pause.");
                    al_draw_text (font, al_map_rgb(0,255,0), 0, 84, 0, "ESC to exit.");
                }
            }
            al_flip_display();
        }
    }

    for (size_t t=0; t<THREAD_N; ++t) {
        scene_thread[t].stop();
        image_buf[t].stop();
        al_destroy_bitmap(thread_bmp[t]);
    }
    blender.stop();

    al_destroy_bitmap(morph_bmp);
    al_destroy_font(font);
    al_destroy_timer(timer);
    al_destroy_display(display);
    al_destroy_event_queue(event_queue);

    return 0;
}

void render_morph(AM_IMAGE *img, ALLEGRO_BITMAP *to) {
    // Clear old bitmap:
    al_set_target_bitmap(to);
    al_set_blender(ALLEGRO_DEST_MINUS_SRC, ALLEGRO_ALPHA, ALLEGRO_INVERSE_ALPHA);
    al_draw_filled_rectangle(0.0, 0.0, MORPH_W, MORPH_H, al_map_rgba(0,0,0,255));

    // Prepare to render:
    al_set_blender(ALLEGRO_ADD, ALLEGRO_ALPHA, ALLEGRO_INVERSE_ALPHA);
    ALLEGRO_LOCKED_REGION * lock = al_lock_bitmap(to,ALLEGRO_PIXEL_FORMAT_ANY,ALLEGRO_LOCK_READWRITE);
    if (lock == NULL) return;

    // Put the pixels:
    size_t pixels = img->pixel_count();
    for (size_t i=0; i<pixels; ++i) {
        int x,y;
        unsigned char r,g,b,a;

        img->get_xy(i, &x, &y);
        img->get_rgba(i, &r, &g, &b, &a);

        al_put_pixel(x, y, al_map_rgba(r,g,b,a));
    }

    // Finally unlock the bitmap:
    if (lock) al_unlock_bitmap(to);
}

void blend_morphs(AM_BLENDER *blender, ALLEGRO_BITMAP *to) {
    // Clear old bitmap:
    al_set_target_bitmap(to);
    al_set_blender(ALLEGRO_DEST_MINUS_SRC, ALLEGRO_ALPHA, ALLEGRO_INVERSE_ALPHA);
    al_draw_filled_rectangle(0.0, 0.0, MORPH_W, MORPH_H, al_map_rgba(0,0,0,255));

    // Prepare to render:
    al_set_blender(ALLEGRO_ADD, ALLEGRO_ALPHA, ALLEGRO_INVERSE_ALPHA);
    ALLEGRO_LOCKED_REGION * lock = al_lock_bitmap(to,ALLEGRO_PIXEL_FORMAT_ANY,ALLEGRO_LOCK_READWRITE);
    if (lock == NULL) return;

    // Put the pixels:
    size_t pixels = blender->pixel_count();
    for (size_t i=0; i<pixels; ++i) {
        int x,y;
        unsigned char r,g,b,a;

        blender->get_xy(i, &x, &y);
        blender->get_rgba(i, &r, &g, &b, &a);

        al_put_pixel(x, y, al_map_rgba(r,g,b,a));
    }

    // Finally unlock the bitmap:
    if (lock) al_unlock_bitmap(to);
}

bool init(int argc, char **argv) {
    if (true == (am_get_warning()&AM_WARN_POINTER_SIZE)) {
        fprintf(stderr, "Pointer size is insufficiently small.\n");
    }
    if (true == (am_get_warning()&AM_WARN_ATOM_SIZE)) {
        fprintf(stderr, "Atom size (%lu) is larger than optimal (%lu).\n",
            sizeof(AM_ATOM),
            sizeof(void *)
        );
    }

    if(!al_init()) {
        fprintf(stderr, "failed to initialize allegro!\n");
        return false;
    }

    if(!al_install_keyboard()) {
        fprintf(stderr, "failed to initialize the keyboard!\n");
        return false;
    }

    al_install_mouse();
    al_init_image_addon();
    al_init_font_addon();
    al_init_primitives_addon();

    timer = al_create_timer(1.0 / FPS);
    if(!timer) {
        fprintf(stderr, "failed to create timer!\n");
        return false;
    }

    display = al_create_display(SCREEN_W, SCREEN_H);
    if(!display) {
        fprintf(stderr, "failed to create display!\n");
        al_destroy_timer(timer);
        return false;
    }

    al_set_new_bitmap_flags(ALLEGRO_MAG_LINEAR|ALLEGRO_MIN_LINEAR);

    font = al_load_font("data/fixed_font.tga", 0, 0);
    if (font==NULL) {
        fprintf(stderr, "failed to load font!\n");
        al_destroy_display(display);
        al_destroy_timer(timer);
        return false;
    }

    al_set_target_bitmap(al_get_backbuffer(display));

    event_queue = al_create_event_queue();
    if(!event_queue) {
        fprintf(stderr, "failed to create event_queue!\n");
        al_destroy_display(display);
        al_destroy_timer(timer);
        al_destroy_font(font);
        return false;
    }

    al_register_event_source(event_queue, al_get_display_event_source(display));
    al_register_event_source(event_queue, al_get_timer_event_source(timer));
    al_register_event_source(event_queue, al_get_keyboard_event_source());
    al_register_event_source(event_queue, al_get_mouse_event_source());

    al_clear_to_color(al_map_rgb(0,0,0));
    al_draw_textf(font, al_map_rgb(0,255,0), SCREEN_W/2, SCREEN_H/2,
                  ALLEGRO_ALIGN_CENTRE,
                  "LOADING...");

    al_flip_display();

    al_start_timer(timer);
    calculate_fps();

    return true;
}

int round_int( double r ) {
    return (r > 0.0) ? (r + 0.5) : (r - 0.5);
}

int calculate_fps() {
    static int times = 0;
    static double old_time = 0.0;
    static double delta_sum = 0.0;
    static int old_fps = -1;

    static bool first = true;
    if (first) {
        first = false;
        old_time = al_get_time();
        return -1;
    }

    int rec_times    = 0;
    int max_times    = round_int(FPS);
    double new_time  = al_get_time();
    double delta     = new_time - old_time;
    delta_sum += delta;
    old_time   = new_time;
    double p  = delta_sum * max_times;
    rec_times = round_int(p);

    if (times > rec_times) {
        return -1;
    }
    times++;

    int fps = 0;
    if (delta_sum >= 1.0 || times>=max_times) {
        fps = times;
        old_fps = fps;
        times=0;
        delta_sum=0.0;
    }
    else {
        if (old_fps == -1) fps = times;
        else               fps = old_fps;
    }

    return fps;
}
