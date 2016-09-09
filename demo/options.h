/*
 * See Copyright Notice in main.h
 */
#include <getopt.h>

class MORPH_OPTIONS {
    public:
    MORPH_OPTIONS()  {}
    ~MORPH_OPTIONS() {}

    const int AS_TEXTURE = am::TEXTURE;
    const int AS_AVERAGE = am::AVERAGE;
    const int AS_DISTINCT= am::DISTINCT;

    const int BY_RGB     = am::RGB;
    const int BY_HSP     = am::HSP;
    
    const int MOTION_NONE   = am::NONE;
    const int MOTION_LINEAR = am::LINEAR;
    const int MOTION_SPLINE = am::SPLINE;
    
    const int FADING_NONE   = am::NONE;
    const int FADING_LINEAR = am::LINEAR;
    const int FADING_COSINE = am::COSINE;
    const int FADING_PERLIN = am::PERLIN;
    
    int verbose        =  0;
    int blend_blobs    =  0;
    int finite         =  0;
    int keep_background=  0;
    int show_blobs     =  AS_TEXTURE;
    int differ_blobs   =  BY_HSP;
    int motion         =  MOTION_SPLINE;
    int fading         =  FADING_PERLIN;
    int exit_flag      =  0;
    std::string name   = "";
    std::string indir  = ".";
    std::string outdir = ".";
    std::string file   = "untitled_morph";
    int height         =  0;
    int width          =  0;
    unsigned seed      =  0;
    unsigned frames_out=  0;
    unsigned match_time=  1;
    unsigned morph_time=  3;
    unsigned feather  =   0;
    unsigned fluid    =   0;
    
    unsigned cycle_length=100000;
    unsigned threads     =8;

    size_t        blob_number      = 1;
    size_t        blob_max_size    = SIZE_MAX;
    size_t        blob_min_size    = 1;
    uint16_t      blob_box_grip    = UINT16_MAX;
    size_t        blob_box_samples = 100;
    double        blob_threshold   = 1.0;
    unsigned char blob_rgba_weight = 1;
    unsigned char blob_size_weight = 1;
    unsigned char blob_xy_weight   = 1;
    size_t        degenerate       = 10000;
    uint16_t      density          = 2;

    bool frames_out_defined = false;

    std::vector<std::string> files;
    
    inline void print_usage (FILE* stream) {
        fprintf (stream, "Usage:   %s [options] <input files...>\n", name.c_str());
        fprintf (stream, "Example: %s data/otu_1.png data/otu_2.png data/otu_3.png --verbose\n", name.c_str());
        fprintf (stream, "Options:\n");        
        fprintf (stream,
            "      --blend-blobs             Blob intersections are blended together.\n"
            "  -E  --blob-feather      INT   Number of gradually transparent outer layers.\n"
            "  -c  --blob-rgba-weight[0-255] Importance of a blob's color when matching.\n"
            "  -z  --blob-size-weight[0-255] Importance of a blob's size when matching.\n"
            "  -p  --blob-xy-weight  [0-255] Importance of a blob's location when matching.\n"
            "  -b  --blobs             INT   Preferred number of blobs to keep.\n"
            "      --blobs-as-texture        A blob is textured by its pixels (default).\n"
            "      --blobs-as-average        A blob has the average color of its pixels.\n"
            "      --blobs-as-distinct       A blob has a random color.\n"
            "      --blobs-by-rgb            Differentiate blobs by RGB distance.\n"
            "      --blobs-by-hsp            Differentiate blobs by HSP distance (default).\n"
            "  -B  --blobs-max-size    INT   Maximum size of a single blob in pixels.\n"
            "  -m  --blobs-min-size    INT   Minimum size of a single blob in pixels.\n"
            "  -g  --blobs-box-grip    INT   Unifying grip bounds for undersized blobs.\n"
            "  -S  --blobs-box-samples INT   Number of samples to take for each dust box.\n"
            "  -t  --blobs-threshold [0-255] Maximum color difference for merging blobs.\n"
            "      --brief                   Print brief messages (default).\n"
            "  -C  --cycle-length      INT   Number of iterations per morph cycle.\n"
            "  -d  --degenerate        INT   Degeneration period for energy minimzation.\n"
            "  -D  --density           INT   Number of atoms per pixel at minimum.\n"            
            "      --fading-none             Colours are not interpolated.\n"
            "      --fading-linear           Colours are interpolated linearly.\n"
            "      --fading-cosine           Colours are cosine interpolated.\n"
            "      --fading-perlin           Use Perlin noise colour transition (default).\n"
            "  -f  --file              STR   Output files with prefix.\n"
            "      --finite                  Morph will not repeat itself seamlessly.\n"            
            "  -L  --fluid             INT   Number of fluid simulation steps per frame.\n"
            "  -F  --frames            INT   Number of frames to generate.\n"
            "  -h  --help                    Display this usage information.\n"
            "  -i  --indir             STR   Read input from this directory.\n"
            "      --keep-background         Morph background is cross-dissolved.\n"
            "  -M  --match-time        INT   Time in seconds given for blob matching.\n"
            "  -O  --morph-time        INT   Time in seconds given for atom morphing.\n"
            "      --motion-none             Positions are not interpolated.\n"
            "      --motion-linear           Positions are linearly interpolated.\n"
            "      --motion-spline           Uses Catmull-Rom splines (default).\n"            
            "  -o  --outdir            STR   Write output to this directory.\n"
            "  -s  --seed              INT   Seed for the random number generator.\n"
            "  -T  --threads           INT   Number of worker threads to spawn.\n"
            "      --verbose                 Print verbose messages.\n"
            "  -v  --version                 Show version information.\n"
        );
    }    
    
    inline bool parse(int argc, char **argv) {
        int c;
        name = argv[0];
        while (1) {
            static struct option long_options[] = {
                // These options set a flag:
                {"verbose",             no_argument,         &verbose,           1 },
                {"blend-blobs",         no_argument,     &blend_blobs,           1 },
                {"finite",              no_argument,          &finite,           1 },
                {"blobs-as-texture",    no_argument,      &show_blobs,  AS_TEXTURE },
                {"blobs-as-average",    no_argument,      &show_blobs,  AS_AVERAGE },
                {"blobs-as-distinct",   no_argument,      &show_blobs, AS_DISTINCT },
                {"blobs-by-rgb",        no_argument,    &differ_blobs,      BY_RGB },
                {"blobs-by-hsp",        no_argument,    &differ_blobs,      BY_HSP },
                {"brief",               no_argument,         &verbose,           0 },
                {"keep-background",     no_argument, &keep_background,           1 },
                // These options don't set a flag. We distinguish them by their indices:
                {"blob-feather",        required_argument,        0,            'E'},
                {"blob-rgba-weight",    required_argument,        0,            'c'},
                {"blob-size-weight",    required_argument,        0,            'z'},
                {"blob-xy-weight",      required_argument,        0,            'p'},                                
                {"blobs",               required_argument,        0,            'b'},
                {"blobs-max-size",      required_argument,        0,            'B'},
                {"blobs-min-size",      required_argument,        0,            'm'},
                {"blobs-box-grip",      required_argument,        0,            'g'},
                {"blobs-box-samples",   required_argument,        0,            'S'},
                {"blobs-threshold",     required_argument,        0,            't'},
                {"cycle-length",        required_argument,        0,            'C'},
                {"degenerate",          required_argument,        0,            'd'},
                {"density",             required_argument,        0,            'D'},
                {"fading-none",         no_argument,        &fading,   FADING_NONE },
                {"fading-linear",       no_argument,        &fading, FADING_LINEAR },
                {"fading-cosine",       no_argument,        &fading, FADING_COSINE },
                {"fading-perlin",       no_argument,        &fading, FADING_PERLIN },                
                {"file",                required_argument,        0,            'f'},
                {"fluid",               required_argument,        0,            'L'},
                {"frames",              required_argument,        0,            'F'},
                {"help",                no_argument,              0,            'h'},
                {"indir",               required_argument,        0,            'i'},
                {"match-time",          required_argument,        0,            'M'},
                {"morph-time",          required_argument,        0,            'O'},
                {"motion-none",         no_argument,        &motion,   MOTION_NONE },
                {"motion-linear",       no_argument,        &motion, MOTION_LINEAR },
                {"motion-spline",       no_argument,        &motion, MOTION_SPLINE },
                {"outdir",              required_argument,        0,            'o'},
                {"seed",                required_argument,        0,            's'},
                {"threads",             required_argument,        0,            'T'},
                {"version",             no_argument,              0,            'v'},
                {0,                     0,                        0,             0 }
            };
            
            // getopt_long stores the option index here.
            int option_index = 0;
            
            c = getopt_long(argc, argv, "b:B:c:C:d:D:E:m:M:O:g:p:S:t:T:f:L:F:hi:o:s:vz:", long_options, &option_index);
            
            /* Detect the end of the options. */
            if (c == -1) break;
            
            switch (c) {
                case 0: // If this option set a flag, do nothing else now.            
                    if (long_options[option_index].flag != 0) break;
                    printf ("option %s", long_options[option_index].name);
                    if (optarg) printf(" with arg %s", optarg); printf ("\n");
                    break;
                case 'E': feather         = atoi(optarg);                     break;    
                case 'c': blob_rgba_weight= std::min(atoi(optarg),255);       break;
                case 'z': blob_size_weight= std::min(atoi(optarg),255);       break;
                case 'p': blob_xy_weight  = std::min(atoi(optarg),255);       break;                                
                case 'b': blob_number     = (size_t) atoi(optarg);            break;
                case 'B': blob_max_size   = atoi(optarg);                     break;
                case 'C': cycle_length    = atoi(optarg);                     break;
                case 'd': degenerate      = (size_t) atoi(optarg);            break;
                case 'D': density         = (uint16_t) atoi(optarg);          break;
                case 'm': blob_min_size   = atoi(optarg);                     break;
                case 'M': match_time      = atoi(optarg);                     break;
                case 'O': morph_time      = atoi(optarg);                     break;
                case 'g': blob_box_grip   = atoi(optarg);                     break;
                case 'S': blob_box_samples= atoi(optarg);                     break;
                case 't': blob_threshold  = std::min(atoi(optarg),255)/255.0; break;
                case 's': seed            = atoi(optarg);                     break;
                case 'T': threads         = std::min(atoi(optarg),1024);      break;
                case 'L': fluid           = atoi(optarg);                     break;
                case 'F': frames_out      = atoi(optarg);                     break;
                case 'i': indir           = optarg;                           break;
                case 'o': outdir          = optarg;                           break;
                case 'f': file            = optarg;                           break;
                case 'h': print_usage(stdout); exit_flag = 1;                 break;
                case 'v':
                    printf("AtoMorph %s Copyright (C) 2013-2014 Erich Erstu\n", am::get_version());
                    exit_flag = 1;
                    break;
                case '?':
                    /* getopt_long already printed an error message. */
                    break;
                default: return false;
            }
        }

        while (optind < argc) {
            files.push_back(argv[optind++]);
        }
        if (frames_out == 0) {
            if (files.size() > 1) frames_out = files.size()*2;
            else                  frames_out = 1;
        }
        return true;
    }
};

