#include <cistem_config.h>


#include "../../core/core_headers.h"


#include "../../constants/constants.h"




class
        MatchTemplateApp : public MyApp {
  public:
    bool DoCalculation( );
    void DoInteractiveUserInput( );


    float GetMaxJobWaitTimeInSeconds( ) { return 360.0f; }

  private:
};

IMPLEMENT_APP(MatchTemplateApp)



// override the DoInteractiveUserInput

void MatchTemplateApp::DoInteractiveUserInput( ) {
    wxString input_search_images = "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/00040_3_0.mrc";
    wxString input_reconstruction =  "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/7ood.mrc";
    wxString cc_output_file; // cross correlation output file 
    wxString mip_output_file = "";
    wxString best_psi_output_file ="";
    wxString best_theta_output_file = "";
    wxString best_phi_output_file = "";
    wxString best_defocus_output_file ="";
    wxString best_pixel_size_output_file="" ;

    wxString output_histogram_file="";
    wxString correlation_std_output_file="";
    wxString correlation_avg_output_file="";
    wxString scaled_mip_output_file="";
    wxString starfile_file = "/home/useradmin/Match_PCA_template_repo/cisTEM/src/programs/match_pca_template/peak_angles.txt";


    float pixel_size              = 1.5f;
    float voltage_kV              = 300.0f;
    float spherical_aberration_mm = 2.7f;
    float amplitude_contrast      = 0.07f;
    float defocus1                = 13849.74f;
    float defocus2                = 13271.80f;
    ;
    float    defocus_angle  = -4.5;
    float    phase_shift = 0;
    float    low_resolution_limit      = 300.0;
    float    high_resolution_limit     = 3.0;
    float    angular_step              = 2.5;
    int      best_parameters_to_keep   = 20;
    float    defocus_search_range      = 500;
    float    defocus_step              = 50;
    float    pixel_size_search_range   = 0.1f;
    float    pixel_size_step           = 0.02f;
    float    padding                   = 1.0;
    bool     ctf_refinement            = false;
    float    particle_radius_angstroms = 0.0f;
    wxString my_symmetry               = "C1";
    float    in_plane_angular_step     = 0;
    bool     use_gpu_input             = false;
    int      max_threads               = 1; // Only used for the GPU code

    UserInput* my_input = new UserInput("MatchTemplate", 1.00);

    input_search_images         = my_input->GetFilenameFromUser("Input images to be searched", "The input image stack, containing the images that should be searched", "image_stack.mrc", true);
    input_reconstruction        = my_input->GetFilenameFromUser("Input template reconstruction", "The 3D reconstruction from which projections are calculated", "reconstruction.mrc", true);
    // mip_output_file             = my_input->GetFilenameFromUser("Output MIP file", "The file for saving the maximum intensity projection image", "mip.mrc", false);
    // scaled_mip_output_file      = my_input->GetFilenameFromUser("Output Scaled MIP file", "The file for saving the maximum intensity projection image divided by correlation variance", "mip_scaled.mrc", false);
    // best_psi_output_file        = my_input->GetFilenameFromUser("Output psi file", "The file for saving the best psi image", "psi.mrc", false);
    // best_theta_output_file      = my_input->GetFilenameFromUser("Output theta file", "The file for saving the best psi image", "theta.mrc", false);
    // best_phi_output_file        = my_input->GetFilenameFromUser("Output phi file", "The file for saving the best psi image", "phi.mrc", false);
    // best_defocus_output_file    = my_input->GetFilenameFromUser("Output defocus file", "The file for saving the best defocus image", "defocus.mrc", false);
    // best_pixel_size_output_file = my_input->GetFilenameFromUser("Output pixel size file", "The file for saving the best pixel size image", "pixel_size.mrc", false);
    // correlation_avg_output_file = my_input->GetFilenameFromUser("Correlation average value", "The file for saving the average value of all correlation images", "corr_average.mrc", false);
    // correlation_std_output_file = my_input->GetFilenameFromUser("Correlation variance output file", "The file for saving the variance of all correlation images", "corr_variance.mrc", false);
    // output_histogram_file       = my_input->GetFilenameFromUser("Output histogram of correlation values", "histogram of all correlation values", "histogram.txt", false);
    // pixel_size                  = my_input->GetFloatFromUser("Pixel size of images (A)", "Pixel size of input images in Angstroms", "1.0", 0.0);
    // voltage_kV                  = my_input->GetFloatFromUser("Beam energy (keV)", "The energy of the electron beam used to image the sample in kilo electron volts", "300.0", 0.0);
    // spherical_aberration_mm     = my_input->GetFloatFromUser("Spherical aberration (mm)", "Spherical aberration of the objective lens in millimeters", "2.7");
    // amplitude_contrast          = my_input->GetFloatFromUser("Amplitude contrast", "Assumed amplitude contrast", "0.07", 0.0, 1.0);
    // defocus1                    = my_input->GetFloatFromUser("Defocus1 (angstroms)", "Defocus1 for the input image", "10000", 0.0);
    // defocus2                    = my_input->GetFloatFromUser("Defocus2 (angstroms)", "Defocus2 for the input image", "10000", 0.0);
    // defocus_angle               = my_input->GetFloatFromUser("Defocus Angle (degrees)", "Defocus Angle for the input image", "0.0");
    // phase_shift                 = my_input->GetFloatFromUser("Phase Shift (degrees)", "Additional phase shift in degrees", "0.0");
    // //    low_resolution_limit = my_input->GetFloatFromUser("Low resolution limit (A)", "Low resolution limit of the data used for alignment in Angstroms", "300.0", 0.0);
    // high_resolution_limit = my_input->GetFloatFromUser("High resolution limit (A)", "High resolution limit of the data used for alignment in Angstroms", "8.0", 0.0);
    // angular_step          = my_input->GetFloatFromUser("Out of plane angular step (0.0 = set automatically)", "Angular step size for global grid search", "0.0", 0.0);
    // in_plane_angular_step = my_input->GetFloatFromUser("In plane angular step (0.0 = set automatically)", "Angular step size for in-plane rotations during the search", "0.0", 0.0);
    // //    best_parameters_to_keep = my_input->GetIntFromUser("Number of top hits to refine", "The number of best global search orientations to refine locally", "20", 1);
    // defocus_search_range    = my_input->GetFloatFromUser("Defocus search range (A)", "Search range (-value ... + value) around current defocus", "500.0", 0.0);
    // defocus_step            = my_input->GetFloatFromUser("Defocus step (A) (0.0 = no search)", "Step size used in the defocus search", "50.0", 0.0);
    // pixel_size_search_range = my_input->GetFloatFromUser("Pixel size search range (A)", "Search range (-value ... + value) around current pixel size", "0.1", 0.0);
    // pixel_size_step         = my_input->GetFloatFromUser("Pixel size step (A) (0.0 = no search)", "Step size used in the pixel size search", "0.01", 0.0);
    // padding                 = my_input->GetFloatFromUser("Padding factor", "Factor determining how much the input volume is padded to improve projections", "1.0", 1.0, 2.0);
    // //    ctf_refinement = my_input->GetYesNoFromUser("Refine defocus", "Should the particle defocus be refined?", "No");
    // particle_radius_angstroms = my_input->GetFloatFromUser("Mask radius for global search (A) (0.0 = max)", "Radius of a circular mask to be applied to the input images during global search", "0.0", 0.0);
    // my_symmetry               = my_input->GetSymmetryFromUser("Template symmetry", "The symmetry of the template reconstruction", "C1");


    // starfile_file         = my_input->GetFilenameFromUser("Starfile", "File containing the Phi and Theta and Psi values for search", "orientations.txt", false);
    cc_output_file            = my_input->GetFilenameFromUser("Output cross correlation filename","", "output_cc.mrc", false);


    int   first_search_position           = -1;
    int   last_search_position            = -1;
    int   image_number_for_gui            = 3;
    int   number_of_jobs_per_image_in_gui = 0;
    float min_peak_radius                 = 10.0f;

    wxString directory_for_results = "/dev/null"; // shouldn't be used in interactive
    wxString result_filename       = "/dev/null"; // shouldn't be used in interactive

    delete my_input;
    const char* jop_code_arg_string = "ttffffffffffifffffbfftttttttttftiiiitttfbitt";
    my_current_job.ManualSetArguments(jop_code_arg_string, input_search_images.ToUTF8( ).data( ),
                                      input_reconstruction.ToUTF8( ).data( ),
                                      pixel_size,
                                      voltage_kV,
                                      spherical_aberration_mm,
                                      amplitude_contrast,
                                      defocus1,
                                      defocus2,
                                      defocus_angle,
                                      low_resolution_limit,
                                      high_resolution_limit,
                                      angular_step,
                                      best_parameters_to_keep,
                                      defocus_search_range,
                                      defocus_step,
                                      pixel_size_search_range,
                                      pixel_size_step,
                                      padding,
                                      ctf_refinement,
                                      particle_radius_angstroms,
                                      phase_shift,
                                      mip_output_file.ToUTF8( ).data( ),
                                      best_psi_output_file.ToUTF8( ).data( ),
                                      best_theta_output_file.ToUTF8( ).data( ),
                                      best_phi_output_file.ToUTF8( ).data( ),
                                      best_defocus_output_file.ToUTF8( ).data( ),
                                      best_pixel_size_output_file.ToUTF8( ).data( ),
                                      scaled_mip_output_file.ToUTF8( ).data( ),
                                      correlation_avg_output_file.ToUTF8( ).data( ),
                                      my_symmetry.ToUTF8( ).data( ),
                                      in_plane_angular_step,
                                      output_histogram_file.ToUTF8( ).data( ),
                                      first_search_position,
                                      last_search_position,
                                      image_number_for_gui,
                                      number_of_jobs_per_image_in_gui,
                                      correlation_std_output_file.ToUTF8( ).data( ),
                                      directory_for_results.ToUTF8( ).data( ),
                                      result_filename.ToUTF8( ).data( ),
                                      min_peak_radius,
                                      use_gpu_input,
                                      max_threads
                                      ,
                                      starfile_file.ToUTF8( ).data( ),
                                      cc_output_file.ToUTF8( ).data( ));
}

// override the do calculation method which will be what is actually run..

bool MatchTemplateApp::DoCalculation( ) {

    bool is_rotated_by_90 = false;
    wxDateTime start_time = wxDateTime::Now( );

    wxString input_search_images_filename  = my_current_job.arguments[0].ReturnStringArgument( );
    wxString input_reconstruction_filename = my_current_job.arguments[1].ReturnStringArgument( );
    float    pixel_size                    = my_current_job.arguments[2].ReturnFloatArgument( );
    float    voltage_kV                    = my_current_job.arguments[3].ReturnFloatArgument( );
    float    spherical_aberration_mm       = my_current_job.arguments[4].ReturnFloatArgument( );
    float    amplitude_contrast            = my_current_job.arguments[5].ReturnFloatArgument( );
    float    defocus1                      = my_current_job.arguments[6].ReturnFloatArgument( );
    float    defocus2                      = my_current_job.arguments[7].ReturnFloatArgument( );
    float    defocus_angle                 = my_current_job.arguments[8].ReturnFloatArgument( );
    ;
    float    low_resolution_limit            = my_current_job.arguments[9].ReturnFloatArgument( );
    float    high_resolution_limit_search    = my_current_job.arguments[10].ReturnFloatArgument( );
    float    angular_step                    = my_current_job.arguments[11].ReturnFloatArgument( );
    int      best_parameters_to_keep         = my_current_job.arguments[12].ReturnIntegerArgument( );
    float    defocus_search_range            = my_current_job.arguments[13].ReturnFloatArgument( );
    float    defocus_step                    = my_current_job.arguments[14].ReturnFloatArgument( );
    float    pixel_size_search_range         = my_current_job.arguments[15].ReturnFloatArgument( );
    float    pixel_size_step                 = my_current_job.arguments[16].ReturnFloatArgument( );
    float    padding                         = my_current_job.arguments[17].ReturnFloatArgument( );
    bool     ctf_refinement                  = my_current_job.arguments[18].ReturnBoolArgument( );
    float    particle_radius_angstroms       = my_current_job.arguments[19].ReturnFloatArgument( );
    float    phase_shift                     = my_current_job.arguments[20].ReturnFloatArgument( );
    wxString mip_output_file                 = my_current_job.arguments[21].ReturnStringArgument( );
    wxString best_psi_output_file            = my_current_job.arguments[22].ReturnStringArgument( );
    wxString best_theta_output_file          = my_current_job.arguments[23].ReturnStringArgument( );
    wxString best_phi_output_file            = my_current_job.arguments[24].ReturnStringArgument( );
    wxString best_defocus_output_file        = my_current_job.arguments[25].ReturnStringArgument( );
    wxString best_pixel_size_output_file     = my_current_job.arguments[26].ReturnStringArgument( );
    wxString scaled_mip_output_file          = my_current_job.arguments[27].ReturnStringArgument( );
    wxString correlation_avg_output_file     = my_current_job.arguments[28].ReturnStringArgument( );
    wxString my_symmetry                     = my_current_job.arguments[29].ReturnStringArgument( );
    float    in_plane_angular_step           = my_current_job.arguments[30].ReturnFloatArgument( );
    wxString output_histogram_file           = my_current_job.arguments[31].ReturnStringArgument( );
    int      first_search_position           = my_current_job.arguments[32].ReturnIntegerArgument( );
    int      last_search_position            = my_current_job.arguments[33].ReturnIntegerArgument( );
    int      image_number_for_gui            = my_current_job.arguments[34].ReturnIntegerArgument( );
    int      number_of_jobs_per_image_in_gui = my_current_job.arguments[35].ReturnIntegerArgument( );
    wxString correlation_std_output_file     = my_current_job.arguments[36].ReturnStringArgument( );
    wxString directory_for_results           = my_current_job.arguments[37].ReturnStringArgument( );
    wxString result_output_filename          = my_current_job.arguments[38].ReturnStringArgument( );
    float    min_peak_radius                 = my_current_job.arguments[39].ReturnFloatArgument( );
    bool     use_gpu                         = my_current_job.arguments[40].ReturnBoolArgument( );
    int      max_threads                     = my_current_job.arguments[41].ReturnIntegerArgument( );
    wxString starfile_file = my_current_job.arguments[42].ReturnStringArgument( );
    wxString  cc_output_file = my_current_job.arguments[43].ReturnStringArgument( );



    ParameterMap parameter_map; // needed for euler search init
    //for (int i = 0; i < 5; i++) {parameter_map[i] = true;}
    parameter_map.SetAllTrue( );

    float outer_mask_radius;
    float current_psi;
    float psi_step;
    float psi_max;
    float psi_start;
    float histogram_step;

    float expected_threshold;
    float actual_number_of_ccs_calculated;

    int current_bin;

    float  temp_float;
    float  variance;
    double temp_double;
    double temp_double_array[5];
    float  factor_score;

    int  number_of_rotations;
    long total_correlation_positions;
    long current_correlation_position;
    long total_correlation_positions_per_thread;
    long pixel_counter;

    int current_search_position;
    int current_x;
    int current_y;

    int factorizable_x;
    int factorizable_y;
    int factor_result_pos;
    int factor_result_neg;

    int defocus_i;
    int size_i;

    int i;

    long   original_input_image_x;
    long   original_input_image_y;
    int    remove_npix_from_edge = 0;
    double sqrt_input_pixels;

    EulerSearch     global_euler_search;
    AnglesAndShifts angles;


    int number_of_search_positions = 0;
    //RD : To open the text file containing the orientations from Starfile
    NumericTextFile starfile_binning(starfile_file, OPEN_TO_READ, 0);
    NumericTextFile pixel_file;
    ImageFile input_search_image_file;
    ImageFile input_reconstruction_file;

    Curve whitening_filter;
    Curve number_of_terms;

    input_search_image_file.OpenFile(input_search_images_filename.ToStdString( ), false);
    input_reconstruction_file.OpenFile(input_reconstruction_filename.ToStdString( ), false);

    //
    remove_npix_from_edge = myroundint(particle_radius_angstroms / pixel_size);
    //    wxPrintf("Removing %d pixels around the edge.\n", remove_npix_from_edge);

    Image input_image;
    Image padded_reference;
    Image input_reconstruction;
    Image template_reconstruction;
    Image current_projection;
    Image padded_projection;

    Image projection_filter;

    Image max_intensity_projection;

    Image best_psi;
    Image best_theta;
    Image best_phi;
    Image best_defocus;
    Image best_pixel_size;

    Image correlation_pixel_sum_image;
    Image correlation_pixel_sum_of_squares_image;

    Image temp_image;

    MRCFile cc_output;

    input_image.ReadSlice(&input_search_image_file, 1);

    // Resize input image to be factorizable by small numbers
    original_input_image_x = input_image.logical_x_dimension;
    original_input_image_y = input_image.logical_y_dimension;
    factorizable_x         = input_image.logical_x_dimension;
    factorizable_y         = input_image.logical_y_dimension;

    bool      DO_FACTORIZATION                       = true;
    bool      MUST_BE_POWER_OF_TWO                   = false; // Required for half-precision xforms
    bool      MUST_BE_FACTOR_OF_FOUR                 = true; // May be faster
    const int max_number_primes                      = 6;
    int       primes[max_number_primes]              = {2, 3, 5, 7, 9, 13};
    float     max_reduction_by_fraction_of_reference = 0.000001f; // FIXME the cpu version is crashing when the image is reduced, but not the GPU
    float     max_increas_by_fraction_of_image       = 0.1f;
    int       max_padding                            = 0; // To restrict histogram calculation
    float     histogram_padding_trim_rescale; // scale the counts to

    // for 5760 this will return
    // 5832 2     2     2     3     3     3     3     3     3 - this is ~ 10% faster than the previous solution BUT
    if ( DO_FACTORIZATION ) {
        for ( i = 0; i < max_number_primes; i++ ) {

            factor_result_neg = ReturnClosestFactorizedLower(original_input_image_x, primes[i], true, MUST_BE_FACTOR_OF_FOUR);
            factor_result_pos = ReturnClosestFactorizedUpper(original_input_image_x, primes[i], true, MUST_BE_FACTOR_OF_FOUR);

            if ( (float)(original_input_image_x - factor_result_neg) < (float)input_reconstruction_file.ReturnXSize( ) * max_reduction_by_fraction_of_reference ) {
                factorizable_x = factor_result_neg;
                break;
            }
            if ( (float)(-original_input_image_x + factor_result_pos) < (float)input_image.logical_x_dimension * max_increas_by_fraction_of_image ) {
                factorizable_x = factor_result_pos;
                break;
            }
        }
        factor_score = FLT_MAX;
        for ( i = 0; i < max_number_primes; i++ ) {

            factor_result_neg = ReturnClosestFactorizedLower(original_input_image_y, primes[i], true, MUST_BE_FACTOR_OF_FOUR);
            factor_result_pos = ReturnClosestFactorizedUpper(original_input_image_y, primes[i], true, MUST_BE_FACTOR_OF_FOUR);

            if ( (float)(original_input_image_y - factor_result_neg) < (float)input_reconstruction_file.ReturnYSize( ) * max_reduction_by_fraction_of_reference ) {
                factorizable_y = factor_result_neg;
                break;
            }
            if ( (float)(-original_input_image_y + factor_result_pos) < (float)input_image.logical_y_dimension * max_increas_by_fraction_of_image ) {
                factorizable_y = factor_result_pos;
                break;
            }
        }
        if ( factorizable_x - original_input_image_x > max_padding )
            max_padding = factorizable_x - original_input_image_x;
        if ( factorizable_y - original_input_image_y > max_padding )
            max_padding = factorizable_y - original_input_image_y;

        if ( ReturnThreadNumberOfCurrentThread( ) == 0 ) {
            wxPrintf("old x, y = %i %i\n  new x, y = %i %i\n", input_image.logical_x_dimension, input_image.logical_y_dimension, factorizable_x, factorizable_y);
        }

        input_image.Resize(factorizable_x, factorizable_y, 1, input_image.ReturnAverageOfRealValuesOnEdges( ));

        input_reconstruction.ReadSlices(&input_reconstruction_file, 1, input_reconstruction_file.ReturnNumberOfSlices( ));
        if ( padding != 1.0f ) {
            input_reconstruction.Resize(input_reconstruction.logical_x_dimension * padding, input_reconstruction.logical_y_dimension * padding, input_reconstruction.logical_z_dimension * padding, input_reconstruction.ReturnAverageOfRealValuesOnEdges( ));
        }

    }

    padded_reference.Allocate(input_image.logical_x_dimension, input_image.logical_y_dimension, 1);
  
    
    padded_reference.SetToConstant(0.0f);






    sqrt_input_pixels = sqrt((double)(input_image.logical_x_dimension * input_image.logical_y_dimension));

   

    CTF input_ctf;
    input_ctf.Init(voltage_kV, spherical_aberration_mm, amplitude_contrast, defocus1, defocus2, defocus_angle, 0.0, 0.0, 0.0, pixel_size, deg_2_rad(phase_shift));

    // assume cube

    current_projection.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
    projection_filter.Allocate(input_reconstruction_file.ReturnXSize( ), input_reconstruction_file.ReturnXSize( ), false);
    template_reconstruction.Allocate(input_reconstruction.logical_x_dimension, input_reconstruction.logical_y_dimension, input_reconstruction.logical_z_dimension, true);
    if ( padding != 1.0f )
        padded_projection.Allocate(input_reconstruction_file.ReturnXSize( ) * padding, input_reconstruction_file.ReturnXSize( ) * padding, false);

    // angular step

    float mask_radius_search;
    if ( particle_radius_angstroms < 1.0f ) {
        mask_radius_search = 200.0f;
    } // This was the original default value.
    else
        mask_radius_search = particle_radius_angstroms;

    if ( angular_step <= 0 ) {
        angular_step = CalculateAngularStep(high_resolution_limit_search, mask_radius_search);
    }

    if ( in_plane_angular_step <= 0 ) {
        psi_step = rad_2_deg(pixel_size / mask_radius_search);
        psi_step = 360.0 / int(360.0 / psi_step + 0.5);
    }
    else {
        psi_step = in_plane_angular_step;
    }

    //psi_start = psi_step / 2.0 * global_random_number_generator.GetUniformRandom();
    psi_start = 0.0f;
    psi_max   = 360.0f;

    //psi_step = 5;

    //wxPrintf("psi_start = %f, psi_max = %f, psi_step = %f\n", psi_start, psi_max, psi_step);

    // search grid

    global_euler_search.InitGrid(my_symmetry, angular_step, 0.0f, 0.0f, psi_max, psi_step, psi_start, pixel_size / high_resolution_limit_search, parameter_map, best_parameters_to_keep);


    // Changing the following 2024-4-5 Bah
    // float orientations[starfile_binning.number_of_lines];

    // This is no good for two reasons:
    // 1) you are allocating an array on the stack which does work for some compilers, but not all. It is generally considered that the
    // size of a stack allocated array should be known at compile time.
    // 2) even if you knew the size, you are allocating an array that is the size of the file, but what you want is an array
    // that has the same number of records per line. These sort of details make or break c++ code.
    std::vector<float> orientations(starfile_binning.records_per_line);

    number_of_search_positions                     = starfile_binning.number_of_lines;
    wxPrintf("Number of searches: %i", number_of_search_positions  );
    global_euler_search.number_of_search_positions = number_of_search_positions;
    Allocate2DFloatArray(global_euler_search.list_of_search_parameters, number_of_search_positions, 3);

    for ( int counter = 0; counter < starfile_binning.number_of_lines; counter++ ) {
        starfile_binning.ReadLine(orientations.data( ));
        global_euler_search.list_of_search_parameters[counter][0] = orientations.at(0);
        global_euler_search.list_of_search_parameters[counter][1] = orientations.at(1);
        global_euler_search.list_of_search_parameters[counter][2] = orientations.at(2);
    }
    starfile_binning.Close( );

    



    // global_euler_search.CalculateGridSearchPositions(false);
    // for now, I am assuming the MTF has been applied already.
    // work out the filter to just whiten the image..

    whitening_filter.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));
    number_of_terms.SetupXAxis(0.0, 0.5 * sqrtf(2.0), int((input_image.logical_x_dimension / 2.0 + 1.0) * sqrtf(2.0) + 1.0));

    wxDateTime my_time_out;
    wxDateTime my_time_in;

    // remove outliers
    // This won't work for movie frames (13.0 is used in unblur) TODO use poisson stats
    input_image.ReplaceOutliersWithMean(5.0f);
    input_image.ForwardFFT( );
    input_image.SwapRealSpaceQuadrants( );

    input_image.ZeroCentralPixel( );
    input_image.Compute1DPowerSpectrumCurve(&whitening_filter, &number_of_terms);
    whitening_filter.SquareRoot( );
    whitening_filter.Reciprocal( );
    whitening_filter.MultiplyByConstant(1.0f / whitening_filter.ReturnMaximumValue( ));

    input_image.ApplyCurveFilter(&whitening_filter);
    input_image.ZeroCentralPixel( );
    // Note: we are dividing by the sqrt of the sum of squares, so the variance in the images 1/N, not 1. This is where the need to multiply the mips by sqrt(N) comes from.
    // Dividing by sqrt(input_image.ReturnSumOfSquares() / N) would result in a properly normalized CCC value.
    input_image.DivideByConstant(sqrtf(input_image.ReturnSumOfSquares( )));
    //input_image.QuickAndDirtyWriteSlice("/tmp/white.mrc", 1);
    //exit(-1);

    // count total searches (lazy)

    total_correlation_positions  = number_of_search_positions ;
    current_correlation_position = 0;

    // if running locally, search over all of them

  
        first_search_position = 0;
        last_search_position  = number_of_search_positions   -1;




    total_correlation_positions = number_of_search_positions ;
    total_correlation_positions_per_thread = total_correlation_positions;

    number_of_rotations = 0;



    ProgressBar* my_progress;

    //Loop over ever search position

    wxPrintf("\n\tFor image id %i\n", image_number_for_gui);
    wxPrintf("Searching %i positions on the Euler sphere (first-last: %i-%i)\n", last_search_position - first_search_position, first_search_position, last_search_position);
    // wxPrintf("Searching %i rotations per position.\n", number_of_rotations);
    wxPrintf("There are %li correlation positions total.\n\n", total_correlation_positions);

    wxPrintf("Performing Search...\n\n");

    actual_number_of_ccs_calculated = 0.0;

    wxDateTime overall_start;
    wxDateTime overall_finish;
    overall_start = wxDateTime::Now( );

cc_output.OpenFile(cc_output_file.ToStdString( ), true);


        input_reconstruction.ChangePixelSize(&template_reconstruction, 1, 0.001f, true);
        //    template_reconstruction.ForwardFFT();
        template_reconstruction.ZeroCentralPixel( );
        template_reconstruction.SwapRealSpaceQuadrants( );




            // make the projection filter, which will be CTF * whitening filter
            input_ctf.SetDefocus(defocus1  / pixel_size, defocus2 / pixel_size, deg_2_rad(defocus_angle));
            //            input_ctf.SetDefocus((defocus1 + 200) / pixel_size, (defocus2 + 200) / pixel_size, deg_2_rad(defocus_angle));
            projection_filter.CalculateCTFImage(input_ctf);
            projection_filter.ApplyCurveFilter(&whitening_filter);

            //            projection_filter.QuickAndDirtyWriteSlices("/tmp/projection_filter.mrc",1,projection_filter.logical_z_dimension,true,1.5);

// MRCFile projection_output;
// wxString projection_output_file = "incorrect__1_projection.mrc";
// projection_output.OpenFile(projection_output_file.ToStdString( ), true);
// template_reconstruction.WriteSlice(&projection_output, 1);


// MRCFile template_output;
// wxString template_output_file = "incorrect__1_template.mrc";
// template_output.OpenFile(template_output_file.ToStdString( ), true);

// MRCFile current_output;
// wxString current_output_file = "incorrect__1_current.mrc";
// current_output.OpenFile(current_output_file.ToStdString( ), true);
             std::string input_mrc_filename  = "/home/useradmin/Project_cisTEM/final_peak_template.mrc";
             MRCFile mrc_file(input_mrc_filename);

            for ( current_search_position = first_search_position; current_search_position <= last_search_position; current_search_position++ ) {
                //loop over each rotation angle

       // wxPrintf("1");

                    //angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], global_euler_search.list_of_search_parameters[current_search_position][2], 0.0, 0.0);
                    //                    angles.Init(130.0, 30.0, 199.5, 0.0, 0.0);
                    angles.Init(global_euler_search.list_of_search_parameters[current_search_position][0], global_euler_search.list_of_search_parameters[current_search_position][1], global_euler_search.list_of_search_parameters[current_search_position][2], 0.0 ,0.0);
                    //angles.Init(297.75, 50.0, 325.0, 0.0 ,0.0);

                    if ( padding != 1.0f ) {
                        template_reconstruction.ExtractSlice(padded_projection, angles, 1.0f, false);
                        padded_projection.SwapRealSpaceQuadrants( );
                        padded_projection.BackwardFFT( );
                        padded_projection.ClipInto(&current_projection);
                        current_projection.ForwardFFT( );
                    }
                   
                   
                    else {
                        //wxPrintf("2");
                        // template_reconstruction.ExtractSlice(current_projection, angles, 1.0f, false);
                        // current_projection.SwapRealSpaceQuadrants( );
                        // current_projection.BackwardFFT( );
                        // current_projection.WriteSlice(&cc_output, current_search_position+1);
                        // return true;
                        // // continue;
                        current_projection.ReadSlice(&mrc_file, current_search_position+1);
                        current_projection.ForwardFFT(false );
                        //current_projection.SwapRealSpaceQuadrants( );
                        // current_projection.QuickAndDirtyWriteSlice("t4.mrc", current_search_position+1, pixel_size);
                        //return true;
                        
                    }

//template_reconstruction.WriteSlice(&template_output, current_search_position +1);
                    current_projection.MultiplyPixelWise(projection_filter);

                    current_projection.BackwardFFT( );
                    //current_projection.ReplaceOutliersWithMean(6.0f);

                    current_projection.AddConstant(-current_projection.ReturnAverageOfRealValuesOnEdges( ));

                    // We want a variance of 1 in the padded FFT. Scale the small SumOfSquares (which is already divided by n) and then re-divide by N.
                    variance = current_projection.ReturnSumOfSquares( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels - powf(current_projection.ReturnAverageOfRealValues( ) * current_projection.number_of_real_space_pixels / padded_reference.number_of_real_space_pixels, 2);
                    current_projection.DivideByConstant(sqrtf(variance));
                    current_projection.ClipIntoLargerRealSpace2D(&padded_reference);

//current_projection.WriteSlice(&current_output, current_search_position +1);
                    // Note: The real space variance is set to 1.0 (for the padded size image) and that results in a variance of N in the FFT do to the scaling of the FFT,
                    // but the FFT values are divided by 1/N so the variance becomes N / (N^2) = is 1/N
                    padded_reference.ForwardFFT( );
                    // Zeroing the central pixel is probably not doing anything useful...
                    padded_reference.ZeroCentralPixel( );

#ifdef MKL
                    // Use the MKL
                    vmcMulByConj(padded_reference.real_memory_allocated / 2, reinterpret_cast<MKL_Complex8*>(input_image.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), reinterpret_cast<MKL_Complex8*>(padded_reference.complex_values), VML_EP | VML_FTZDAZ_ON | VML_ERRMODE_IGNORE);
#else
                    for ( pixel_counter = 0; pixel_counter < padded_reference.real_memory_allocated / 2; pixel_counter++ ) {
                        padded_reference.complex_values[pixel_counter] = conj(padded_reference.complex_values[pixel_counter]) * input_image.complex_values[pixel_counter];
                    }
#endif

                    // Note: the cross correlation will have variance 1/N (the product of variance of the two FFTs assuming the means are both zero and the distributions independent.)
                    // Taking the inverse FFT scales this variance by N resulting in a MIP with variance 1
                    padded_reference.BackwardFFT( );

        // temp_image.CopyFrom(&padded_reference);
        // temp_image.Resize(original_input_image_x, original_input_image_y, current_search_position+1, 1);
        // temp_image.QuickAndDirtyWriteSlice(cc_output_file.ToStdString( ), current_search_position+1, pixel_size);
                    // padded_reference.Resize(original_input_image_x, original_input_image_y, 1);
                    // padded_reference.MultiplyByConstant((float)sqrt_input_pixels);
                    // padded_reference.WriteSlice(&cc_output, current_search_position +1);


                   padded_reference.WriteSlice(&cc_output, current_search_position+1);
                   return true;
                  


                   
                    current_projection.is_in_real_space = false;
                    padded_reference.is_in_real_space   = true;


                
                wxPrintf("\n\n\tIteration: %i\n", current_search_position);
 
        

    }
        wxPrintf("\n\n\tTimings: Overall: %s\n", (wxDateTime::Now( ) - overall_start).Format( ));
 
                // cc_output.SetPixelSize(pixel_size);
                // cc_output.WriteHeader( );




        // max_intensity_projection.MultiplyByConstant((float)sqrt_input_pixels);


        // temp_image.CopyFrom(&max_intensity_projection);
        // temp_image.Resize(original_input_image_x, original_input_image_y, 1, temp_image.ReturnAverageOfRealValuesOnEdges( ));
        // temp_image.QuickAndDirtyWriteSlice(mip_output_file.ToStdString( ), 1, pixel_size);

    return true;
}

